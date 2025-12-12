from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Tuple, Union, TYPE_CHECKING, Any, cast, IO
import contextlib
import os

if TYPE_CHECKING:
    from ase import Atoms


def _max_workers(task_count: int) -> int:
    if task_count <= 1:
        return 1
    return max(1, min(task_count, 24))


def _chunk_size(max_workers: int) -> int:
    return max(1, max_workers * 4)


def _parse_tess_header(mm: Any, file_size: int) -> Dict[str, Any]:
    import msgpack
    import msgpack_numpy as m
    
    # Enable numpy array deserialization for header contents that may include numpy arrays
    m.patch()

    header_start = int.from_bytes(mm[file_size - 8:file_size], 'little')
    header_bytes = mm[header_start:file_size - 8]
    header = msgpack.unpackb(header_bytes, raw=False, strict_map_key=False)

    format_version = header.get('format_version', 1)
    if format_version not in {2, 3, 4, 5}:
        raise ValueError("Invalid format version")

    # Normalize frame_offsets to a list of tuples for convenience
    if 'frame_offsets' in header and isinstance(header['frame_offsets'], list):
        header['frame_offsets'] = [tuple(offset) for offset in header['frame_offsets']]

    return header


def _read_tess_header(filename: str) -> Dict[str, Any]:

    import mmap

    with open(filename, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            return _parse_tess_header(mm, mm.size())
        finally:
            mm.close()


def _decode_frame(frame_bytes: bytes, header: Dict[str, Any], template_atoms: Optional['Atoms'] = None) -> 'Atoms':

    import msgpack
    import msgpack_numpy as m
    from atomict.io.atoms import dict_to_atoms
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    import zlib
    import lz4.block

    # Enable numpy array deserialization for any numpy contents
    m.patch()

    compression_info = header.get('compression')
    if isinstance(compression_info, dict):
        compression_type = compression_info.get('type', 'none')
    elif isinstance(compression_info, str):
        compression_type = compression_info
    else:
        compression_type = 'none'
    compression_type = (compression_type or 'none').lower()

    # Decompress if needed
    if compression_type == 'zlib':
        frame_bytes = zlib.decompress(frame_bytes)
    elif compression_type == 'lz4':
        frame_bytes = lz4.block.decompress(frame_bytes)

    # Unpack frame
    frame_dict = msgpack.unpackb(frame_bytes, raw=False, strict_map_key=False)

    format_version = header.get('format_version', 1)

    # Fast path for v4/v5 using static header data
    if format_version in {4, 5}:
        static_data = header.get('static_data', {})

        # Either copy provided template or construct one from static data
        if template_atoms is not None:
            atoms = template_atoms.copy()
        else:
            numbers = static_data['numbers']
            masses = static_data['masses']
            pbc = static_data['pbc']
            atoms = Atoms(numbers=numbers, pbc=pbc)
            atoms.set_masses(masses)
            if 'cell' in static_data and static_data['cell'] is not None:
                atoms.set_cell(static_data['cell'])
            if 'tags' in static_data and static_data['tags'] is not None:
                atoms.set_tags(static_data['tags'])
            if 'initial_charges' in static_data and static_data['initial_charges'] is not None:
                atoms.set_initial_charges(static_data['initial_charges'])
            if 'initial_magmoms' in static_data and static_data['initial_magmoms'] is not None:
                atoms.set_initial_magnetic_moments(static_data['initial_magmoms'])

        # Per-frame updates
        atoms.set_positions(frame_dict['positions'])
        if static_data.get('cell_changes', True) and 'cell' in frame_dict:
            atoms.set_cell(frame_dict['cell'])

        if format_version == 5:
            if 'tags' in frame_dict:
                atoms.set_tags(frame_dict['tags'])
            if 'initial_charges' in frame_dict:
                atoms.set_initial_charges(frame_dict['initial_charges'])
            if 'initial_magmoms' in frame_dict:
                atoms.set_initial_magnetic_moments(frame_dict['initial_magmoms'])
            if 'momenta' in frame_dict:
                atoms.set_momenta(frame_dict['momenta'])
            if 'forces' in frame_dict:
                atoms.arrays['forces'] = frame_dict['forces']
            if 'info' in frame_dict and isinstance(frame_dict['info'], dict):
                atoms.info.update(frame_dict['info'])
            if 'custom_arrays' in frame_dict and isinstance(frame_dict['custom_arrays'], dict):
                for k, v in frame_dict['custom_arrays'].items():
                    atoms.arrays[k] = v

            # Calculator results (merge stress if present)
            calc_data = {}
            if 'calc_results' in frame_dict and isinstance(frame_dict['calc_results'], dict):
                calc_data = dict(frame_dict['calc_results'])
            if 'stress' in frame_dict:
                calc_data['stress'] = frame_dict['stress']
            if len(calc_data) > 0:
                calc = SinglePointCalculator(atoms)
                for key, value in calc_data.items():
                    if key != 'name':
                        calc.results[key] = value
                if calc.results:
                    atoms.calc = calc

        return atoms

    # Legacy path: v1-v3 store full atoms dict per frame
    atoms_list, _ = dict_to_atoms(frame_dict)
    return atoms_list[0]

def write_tess(
    atoms: Union['Atoms', List['Atoms']],
    filename: Union[str, IO[bytes], os.PathLike],
    metadata: Optional[Dict] = None,
    compression: Optional[str] = 'none',
    compression_level: int = 0,
) -> None:

    import msgpack
    import msgpack_numpy as m
    from ase import Atoms
    import zlib
    import lz4.block
    import numpy as np

    # Enable numpy array serialization
    m.patch()

    if isinstance(atoms, Atoms):
        frames_list = [atoms]
    else:
        frames_list = list(atoms)

    compression_mode = (compression or 'none').lower()
    if compression_mode not in {'none', 'zlib', 'lz4'}:
        raise ValueError(
            f"Unsupported compression mode '{compression_mode}' for tess format"
        )
    if compression_mode == 'zlib' and not (0 <= compression_level <= 9):
        raise ValueError("compression_level must be between 0 and 9 for zlib compression")

    # Collect global unique_symbols
    unique_symbols_set = set()
    for a in frames_list:
        unique_symbols_set.update(a.get_chemical_symbols())
    unique_symbols = sorted(list(unique_symbols_set))
    unique_symbols_lookup = {s: i for i, s in enumerate(unique_symbols)}

    # Check if cell changes across frames
    first_cell = np.array(frames_list[0].get_cell(), dtype=np.float32)
    cell_changes = any(
        not np.allclose(np.array(a.get_cell(), dtype=np.float32), first_cell, rtol=1e-6)
        for a in frames_list[1:20]  # Sample first 20 frames
    )
    
    # Pre-extract positions as float32 arrays for all frames
    positions_list = [a.get_positions().astype(np.float32) for a in frames_list]
    if cell_changes:
        cells_list = [np.array(a.get_cell(), dtype=np.float32) for a in frames_list]

    # Ensure constant atom count and species across frames for tess format
    first_numbers = np.asarray(frames_list[0].get_atomic_numbers())
    for a in frames_list[1:]:
        if len(a) != len(frames_list[0]) or not np.array_equal(np.asarray(a.get_atomic_numbers()), first_numbers):
            raise ValueError("tess format requires constant atom count and species across frames")

    # Determine static vs per-frame ASE arrays
    first_tags = np.asarray(frames_list[0].get_tags())
    tags_static = all(np.array_equal(np.asarray(a.get_tags()), first_tags) for a in frames_list)

    first_init_charges = np.asarray(frames_list[0].get_initial_charges())
    init_charges_static = all(np.array_equal(np.asarray(a.get_initial_charges()), first_init_charges) for a in frames_list)

    first_init_magmoms = np.asarray(frames_list[0].get_initial_magnetic_moments())
    init_magmoms_static = all(np.array_equal(np.asarray(a.get_initial_magnetic_moments()), first_init_magmoms) for a in frames_list)

    # Per-frame dynamic arrays/props presence checks
    has_any_momenta = any(np.any(np.abs(a.get_momenta()) > 1e-10) for a in frames_list)
    has_any_forces = any('forces' in a.arrays for a in frames_list)
    has_any_stress = any(hasattr(a, 'stress') and getattr(a, 'stress') is not None for a in frames_list)

    # Helper to extract calculator results per frame (aligned with atoms_to_dict)
    def _calc_results_for(a):
        calc_data = {}
        calc_found = False
        if hasattr(a, 'calc') and a.calc is not None:
            calc_found = True
            calc_data['name'] = a.calc.__class__.__name__
            for prop in ['energy', 'free_energy', 'forces', 'stress', 'dipole', 'charges', 'magmom', 'magmoms']:
                try:
                    if hasattr(a.calc, 'results') and prop in a.calc.results:
                        calc_data[prop] = a.calc.results[prop]
                    else:
                        value = a.calc.get_property(prop, a)
                        if value is not None:
                            calc_data[prop] = value
                except Exception:
                    pass
        if not calc_found and hasattr(a, 'info'):
            calc_name = a.info.get('_calc_name')
            if calc_name:
                calc_data['name'] = calc_name
                for key, value in a.info.items():
                    if key.startswith('_calc_') and key != '_calc_name':
                        calc_data[key[6:]] = value
        return calc_data if len(calc_data) > 0 else None

    has_any_calc = any(_calc_results_for(a) is not None for a in frames_list)

    # Helper to serialize info dict values (to_dict/todict if available)
    def _serialize_info(info: Dict) -> Dict:
        out = {}
        for k, v in (info or {}).items():
            if hasattr(v, 'to_dict') and callable(v.to_dict):
                out[k] = v.to_dict()
            elif hasattr(v, 'todict') and callable(v.todict):
                out[k] = v.todict()
            else:
                out[k] = v
        return out
    
    # Optimized frame serialization - pack only essential per-frame data  
    def _pack_frame(i: int) -> bytes:
        # Only pack data that changes between frames
        a = frames_list[i]
        frame_dict = {'positions': positions_list[i]}
        if cell_changes:
            frame_dict['cell'] = cells_list[i]
        # Per-frame arrays
        if has_any_momenta:
            frame_dict['momenta'] = a.get_momenta()
        if not tags_static:
            frame_dict['tags'] = a.get_tags()
        if not init_charges_static:
            frame_dict['initial_charges'] = a.get_initial_charges()
        if not init_magmoms_static:
            frame_dict['initial_magmoms'] = a.get_initial_magnetic_moments()
        # Forces/stress
        if has_any_forces and 'forces' in a.arrays:
            frame_dict['forces'] = a.arrays['forces']
        if has_any_stress and hasattr(a, 'stress') and getattr(a, 'stress') is not None:
            frame_dict['stress'] = np.asarray(getattr(a, 'stress'), dtype=np.float64)
        # Calculator results
        if has_any_calc:
            calc_data = _calc_results_for(a)
            if calc_data is not None:
                frame_dict['calc_results'] = calc_data
        # Info dict
        if hasattr(a, 'info') and a.info:
            frame_dict['info'] = _serialize_info(a.info)
        # Custom arrays
        standard_arrays = {'numbers', 'positions', 'momenta', 'masses', 'tags', 'charges', 'forces'}
        custom_arrays = {}
        for key, value in a.arrays.items():
            if key not in standard_arrays:
                custom_arrays[key] = value
        if custom_arrays:
            frame_dict['custom_arrays'] = custom_arrays
        return cast(bytes, msgpack.packb(frame_dict, use_bin_type=True))

    if isinstance(filename, (str, os.PathLike)):
        f_ctx = open(filename, 'wb', buffering=1024*1024)
        f_to_close = True
    else:
        f_ctx = filename
        f_to_close = False

    with contextlib.ExitStack() as stack:
        if f_to_close:
            f = stack.enter_context(f_ctx)
        else:
            f = f_ctx

        frame_offsets: List[Tuple[int, int]] = []
        uncompressed_lengths: List[int] = []
        if compression_mode == 'zlib' or compression_mode == 'lz4':
            # Pre-pack all frames (fast operation)
            packed_frames = [_pack_frame(i) for i in range(len(frames_list))]
            
            worker_count = _max_workers(len(packed_frames))
            
            def _compress_frame(packed: bytes) -> bytes:
                if compression_mode == 'zlib':
                    return zlib.compress(packed, compression_level)
                else:  # lz4
                    return lz4.block.compress(packed, store_size=True)

            initial = f.tell()
            current = initial

            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                compressed_frames = list(executor.map(_compress_frame, packed_frames))
            
            # Write all compressed frames
            for packed, compressed in zip(packed_frames, compressed_frames):
                f.write(compressed)
                frame_offsets.append((current, len(compressed)))
                uncompressed_lengths.append(len(packed))
                current += len(compressed)

        else:
            # Write frames directly without intermediate buffer
            initial = f.tell()
            current = initial
            for i in range(len(frames_list)):
                packed = _pack_frame(i)
                l = len(packed)
                f.write(packed)
                frame_offsets.append((current, l))
                current += l

        # Store static data in header (data that's the same for all frames)
        first_frame = frames_list[0]
        symbols = first_frame.get_chemical_symbols()
        symbols_idx = np.array([unique_symbols_lookup[s] for s in symbols], dtype=np.uint16)
        
        # Write header
        header_dict: Dict[str, Any] = {
            'format_version': 5,
            'metadata': dict(metadata or {}),
            'unique_symbols': unique_symbols,
            'num_frames': len(frames_list),
            'frame_offsets': frame_offsets,
            # Static frame data stored once in header
            'static_data': {
                'n_atoms': len(first_frame),
                'symbols': symbols_idx,
                'numbers': first_frame.get_atomic_numbers(),
                'masses': first_frame.get_masses(),
                'pbc': first_frame.get_pbc(),
                'cell': first_cell,
                'cell_changes': cell_changes,
            }
        }
        # Add static arrays when invariant across frames
        if tags_static:
            header_dict['static_data']['tags'] = first_tags
        if init_charges_static:
            header_dict['static_data']['initial_charges'] = first_init_charges
        if init_magmoms_static:
            header_dict['static_data']['initial_magmoms'] = first_init_magmoms
        if compression_mode != 'none':
            comp_dict: Dict[str, Any] = {'type': compression_mode}
            if compression_mode == 'zlib':
                comp_dict['level'] = int(compression_level)
            header_dict['compression'] = comp_dict
            header_dict['frame_uncompressed_lengths'] = uncompressed_lengths
        header_bytes = cast(bytes, msgpack.packb(header_dict, use_bin_type=True))
        header_start = f.tell()
        f.write(header_bytes)

        # Write header offset
        f.write(header_start.to_bytes(8, 'little'))


def read_tess(filename: Union[str, IO[bytes]], frames_indices: Optional[List[int]] = None) -> Tuple[List['Atoms'], Dict]:

    import msgpack
    import msgpack_numpy as m
    from atomict.io.atoms import dict_to_atoms
    from ase import Atoms
    import mmap
    import zlib
    import lz4.block
    import numpy as np
    from ase.calculators.singlepoint import SinglePointCalculator

    # Enable numpy array deserialization
    m.patch()

    # Manage file context
    if isinstance(filename, (str, os.PathLike)):
        f_ctx = open(filename, 'rb')
        f_to_close = True
    else:
        f_ctx = filename
        f_to_close = False

    with contextlib.ExitStack() as stack:
        if f_to_close:
            f = stack.enter_context(f_ctx)
        else:
            f = f_ctx

        # Determine view and size
        mm: Any = None
        view: Any = None
        size: int = 0
        
        # Try mmap first (efficient for large files)
        if hasattr(f, 'fileno'):
            try:
                mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                stack.callback(mm.close)
                
                # Advise the kernel we'll read sequentially
                if hasattr(mm, 'madvise'):
                    import mmap as mmap_module
                    if hasattr(mmap_module, 'MADV_SEQUENTIAL'):
                        mm.madvise(mmap_module.MADV_SEQUENTIAL)
                
                view = mm
                size = mm.size()
            except Exception:
                # Fallback if mmap fails (e.g. pipe)
                pass
        
        if view is None:
            # Try getting buffer (BytesIO)
            if hasattr(f, 'getbuffer'):
                view = f.getbuffer()
                size = view.nbytes
            else:
                # Fallback: read all into memory
                if hasattr(f, 'seekable') and f.seekable():
                    f.seek(0)
                content = f.read()
                view = content
                size = len(content)

        # Read header
        header = _parse_tess_header(view, size)

        format_version = header.get('format_version', 1)
        compression_info = header.get('compression')
        if isinstance(compression_info, dict):
            compression_type = compression_info.get('type', 'none')
        elif isinstance(compression_info, str):
            compression_type = compression_info
        else:
            compression_type = 'none'
        compression_type = (compression_type or 'none').lower()
        if compression_type not in {'none', 'zlib', 'lz4'}:
            raise ValueError(f"Unsupported compression type '{compression_type}' in tess file")

        metadata = header.get('metadata', {})
        frame_offsets = [tuple(offset) for offset in header['frame_offsets']]
        num_frames = len(frame_offsets)

        # Determine which frames to load
        if frames_indices is not None:
            # Validate and normalize indices
            frames_to_load = []
            for idx in frames_indices:
                if idx < 0 or idx >= num_frames:
                    raise IndexError(f"Frame index {idx} out of range [0, {num_frames})")
                frames_to_load.append(idx)
            # Filter offsets to only requested frames
            filtered_offsets = [(i, frame_offsets[i]) for i in frames_to_load]
        else:
            frames_to_load = list(range(num_frames))
            filtered_offsets = [(i, offset) for i, offset in enumerate(frame_offsets)]

        # Pre-allocate the result list based on requested frames
        atoms_list: List[Optional[Atoms]] = [None] * len(frames_to_load)

        # Handle format version 4/5 with static data
        if format_version in {4, 5}:
            static_data = header.get('static_data', {})
            unique_symbols = header['unique_symbols']
            symbols_idx = static_data['symbols']
            symbols = [unique_symbols[i] for i in symbols_idx]
            numbers = static_data['numbers']
            masses = static_data['masses']
            pbc = static_data['pbc']
            static_cell = static_data.get('cell')
            cell_changes = static_data.get('cell_changes', True)
            static_tags = static_data.get('tags')
            static_init_charges = static_data.get('initial_charges')
            static_init_magmoms = static_data.get('initial_magmoms')

            # Create a template Atoms object to clone for efficiency
            template_atoms = Atoms(numbers=numbers, pbc=pbc)
            template_atoms.set_masses(masses)
            if static_cell is not None:
                template_atoms.set_cell(static_cell)
            if static_tags is not None:
                template_atoms.set_tags(static_tags)
            if static_init_charges is not None:
                template_atoms.set_initial_charges(static_init_charges)
            if static_init_magmoms is not None:
                template_atoms.set_initial_magnetic_moments(static_init_magmoms)
        else:
            convert = dict_to_atoms
            static_data = None
            cell_changes = True

        # Process frames with minimal overhead
        if len(filtered_offsets) > 0:
            # Unified decode path using helper; still parallelized
            def _decode(idx_offset: Tuple[int, Tuple[int, int]]) -> Tuple[int, 'Atoms']:
                orig_idx, (start, length) = idx_offset
                raw_slice = view[start:start + length]
                return orig_idx, _decode_frame(raw_slice, header, template_atoms if format_version in {4, 5} else None)

            worker_count = _max_workers(len(filtered_offsets))
            chunk_size = _chunk_size(worker_count)

            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                for chunk_start in range(0, len(filtered_offsets), chunk_size):
                    offset_chunk = filtered_offsets[chunk_start:chunk_start + chunk_size]
                    for result_idx, (orig_idx, atoms) in enumerate(
                        executor.map(_decode, offset_chunk), start=0
                    ):
                        atoms_list[chunk_start + result_idx] = atoms
    
    return cast(List[Atoms], atoms_list), metadata
