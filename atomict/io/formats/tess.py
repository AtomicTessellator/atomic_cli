from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Tuple, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


def _max_workers(task_count: int) -> int:
    if task_count <= 1:
        return 1
    return max(1, min(task_count, 24))


def _chunk_size(max_workers: int) -> int:
    return max(1, max_workers * 4)


def write_tess(
    atoms: Union['Atoms', List['Atoms']],
    filename: str,
    metadata: Optional[Dict] = None,
    compression: Optional[str] = 'none',
    compression_level: int = 0,
) -> None:

    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
        import zlib
        import lz4.block
        import numpy as np
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[utils]` to use msgpack I/O")

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
    
    # Optimized frame serialization - pack only essential per-frame data  
    def _pack_frame(i: int) -> bytes:
        # Only pack data that changes between frames
        frame_dict = {'positions': positions_list[i]}
        if cell_changes:
            frame_dict['cell'] = cells_list[i]
        return msgpack.packb(frame_dict, use_bin_type=True)

    with open(filename, 'wb', buffering=1024*1024) as f:
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
        header_dict: Dict[str, object] = {
            'format_version': 4,
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
        if compression_mode != 'none':
            header_dict['compression'] = {
                'type': compression_mode,
            }
            if compression_mode == 'zlib':
                header_dict['compression']['level'] = compression_level
            header_dict['frame_uncompressed_lengths'] = uncompressed_lengths
        header_bytes = msgpack.packb(header_dict, use_bin_type=True)
        header_start = f.tell()
        f.write(header_bytes)

        # Write header offset
        f.write(header_start.to_bytes(8, 'little'))


def read_tess(filename: str, frames_indices: Optional[List[int]] = None) -> Tuple[List['Atoms'], Dict]:

    try:
        import msgpack
        import msgpack_numpy as m
        from atomict.io.atoms import dict_to_atoms
        from ase import Atoms
        import mmap
        import zlib
        import lz4.block
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[utils]` to use msgpack I/O")

    # Enable numpy array deserialization
    m.patch()

    with open(filename, 'rb') as f:
        # Use mmap with larger read-ahead hint
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            # Advise the kernel we'll read sequentially
            if hasattr(mm, 'madvise'):
                import mmap as mmap_module
                if hasattr(mmap_module, 'MADV_SEQUENTIAL'):
                    mm.madvise(mmap_module.MADV_SEQUENTIAL)
            
            file_size = mm.size()
            header_start = int.from_bytes(mm[file_size - 8:file_size], 'little')
            header_bytes = mm[header_start:file_size - 8]
            header = msgpack.unpackb(header_bytes, raw=False, strict_map_key=False)

            format_version = header.get('format_version', 1)
            if format_version not in {2, 3, 4}:
                raise ValueError("Invalid format version")

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
            atoms_list = [None] * len(frames_to_load)
            
            # Handle format version 4 with static data
            if format_version == 4:
                static_data = header.get('static_data', {})
                unique_symbols = header['unique_symbols']
                symbols_idx = static_data['symbols']
                symbols = [unique_symbols[i] for i in symbols_idx]
                numbers = static_data['numbers']
                masses = static_data['masses']
                pbc = static_data['pbc']
                static_cell = static_data.get('cell')
                cell_changes = static_data.get('cell_changes', True)
                
                # Create a template Atoms object to clone for efficiency
                import numpy as np
                template_atoms = Atoms(numbers=numbers, pbc=pbc)
                template_atoms.set_masses(masses)
                if static_cell is not None:
                    template_atoms.set_cell(static_cell)
            else:
                convert = dict_to_atoms
                static_data = None
                cell_changes = True
            
            # Process frames with minimal overhead
            if compression_type == 'zlib' and len(filtered_offsets) > 0:

                def _decompress(idx_offset: Tuple[int, Tuple[int, int]]) -> Tuple[int, bytes]:
                    orig_idx, (start, length) = idx_offset
                    compressed = mm[start:start + length]
                    return orig_idx, zlib.decompress(compressed)

                worker_count = _max_workers(len(filtered_offsets))
                chunk_size = _chunk_size(worker_count)

                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    for chunk_start in range(0, len(filtered_offsets), chunk_size):
                        offset_chunk = filtered_offsets[chunk_start:chunk_start + chunk_size]
                        for result_idx, (orig_idx, frame_bytes) in enumerate(
                            executor.map(_decompress, offset_chunk), start=0
                        ):
                            frame_dict = msgpack.unpackb(
                                frame_bytes, raw=False, strict_map_key=False
                            )
                            if format_version == 4:
                                # Fast path for v4: clone template and update only changing data
                                atoms = template_atoms.copy()
                                atoms.set_positions(frame_dict['positions'])
                                if cell_changes and 'cell' in frame_dict:
                                    atoms.set_cell(frame_dict['cell'])
                                atoms_list[chunk_start + result_idx] = atoms
                            else:
                                atoms_list[chunk_start + result_idx] = convert(frame_dict)[0]
            elif compression_type == 'lz4' and len(filtered_offsets) > 0:

                def _decompress(idx_offset: Tuple[int, Tuple[int, int]]) -> Tuple[int, bytes]:
                    orig_idx, (start, length) = idx_offset
                    compressed = mm[start:start + length]
                    return orig_idx, lz4.block.decompress(compressed)

                worker_count = _max_workers(len(filtered_offsets))
                chunk_size = _chunk_size(worker_count)

                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    for chunk_start in range(0, len(filtered_offsets), chunk_size):
                        offset_chunk = filtered_offsets[chunk_start:chunk_start + chunk_size]
                        for result_idx, (orig_idx, frame_bytes) in enumerate(
                            executor.map(_decompress, offset_chunk), start=0
                        ):
                            frame_dict = msgpack.unpackb(
                                frame_bytes, raw=False, strict_map_key=False
                            )
                            if format_version == 4:
                                atoms = template_atoms.copy()
                                atoms.set_positions(frame_dict['positions'])
                                if cell_changes and 'cell' in frame_dict:
                                    atoms.set_cell(frame_dict['cell'])
                                atoms_list[chunk_start + result_idx] = atoms
                            else:
                                atoms_list[chunk_start + result_idx] = convert(frame_dict)[0]
            else:
                for result_idx, (orig_idx, (start, length)) in enumerate(filtered_offsets):
                    frame_slice = mm[start:start + length]
                    frame_dict = msgpack.unpackb(
                        frame_slice, raw=False, strict_map_key=False
                    )
                    if format_version == 4:
                        atoms = template_atoms.copy()
                        atoms.set_positions(frame_dict['positions'])
                        if cell_changes and 'cell' in frame_dict:
                            atoms.set_cell(frame_dict['cell'])
                        atoms_list[result_idx] = atoms
                    else:
                        atoms_list[result_idx] = convert(frame_dict)[0]
        finally:
            mm.close()
    
    return atoms_list, metadata
