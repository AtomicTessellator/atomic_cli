import os
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Tuple, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


def _max_workers(task_count: int) -> int:
    if task_count <= 1:
        return 1
    cpu_total = os.cpu_count() or 1
    return max(1, min(task_count, cpu_total))


def write_tess(
    atoms: Union['Atoms', List['Atoms']],
    filename: str,
    metadata: Optional[Dict] = None,
    compression: Optional[str] = 'zlib',
    compression_level: int = 1,
) -> None:

    from atomict.io.msgpack import atoms_to_dict

    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
        import zlib
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[utils]` to use msgpack I/O")

    # Enable numpy array serialization
    m.patch()

    if isinstance(atoms, Atoms):
        frames_list = [atoms]
    else:
        frames_list = list(atoms)

    compression_mode = (compression or 'none').lower()
    if compression_mode not in {'none', 'zlib'}:
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

    # Write frames
    frame_offsets: List[Tuple[int, int]] = []
    uncompressed_lengths: List[int] = []
    with open(filename, 'wb') as f:
        if compression_mode == 'zlib':
            frame_dicts = [atoms_to_dict([atom_frame]) for atom_frame in frames_list]

            def _pack_and_compress(frame_dict: Dict) -> Tuple[bytes, int]:
                packed = msgpack.packb(frame_dict, use_bin_type=True)
                return zlib.compress(packed, compression_level), len(packed)

            with ThreadPoolExecutor(max_workers=_max_workers(len(frame_dicts))) as executor:
                results = list(executor.map(_pack_and_compress, frame_dicts))

            for compressed_bytes, original_length in results:
                uncompressed_lengths.append(original_length)
                start = f.tell()
                f.write(compressed_bytes)
                frame_offsets.append((start, len(compressed_bytes)))
        else:
            for atom_frame in frames_list:
                frame_dict = atoms_to_dict([atom_frame])
                frame_bytes = msgpack.packb(frame_dict, use_bin_type=True)
                start = f.tell()
                f.write(frame_bytes)
                frame_offsets.append((start, len(frame_bytes)))
        
        # Write header
        header_dict: Dict[str, object] = {
            'format_version': 3 if compression_mode != 'none' else 2,
            'metadata': dict(metadata or {}),
            'unique_symbols': unique_symbols,
            'num_frames': len(frames_list),
            'frame_offsets': frame_offsets,
        }
        if compression_mode != 'none':
            header_dict['compression'] = {
                'type': compression_mode,
                'level': compression_level,
            }
            header_dict['frame_uncompressed_lengths'] = uncompressed_lengths
        header_bytes = msgpack.packb(header_dict, use_bin_type=True)
        header_start = f.tell()
        f.write(header_bytes)
        
        # Write header offset
        f.write(header_start.to_bytes(8, 'little'))


def read_tess(filename: str) -> Tuple[List['Atoms'], Dict]:

    try:
        import msgpack
        import msgpack_numpy as m
        from atomict.io.atoms import dict_to_atoms
        import mmap
        import zlib
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
            if format_version not in {2, 3}:
                raise ValueError("Invalid format version")

            compression_info = header.get('compression')
            if isinstance(compression_info, dict):
                compression_type = compression_info.get('type', 'none')
            elif isinstance(compression_info, str):
                compression_type = compression_info
            else:
                compression_type = 'none'
            compression_type = (compression_type or 'none').lower()
            if compression_type not in {'none', 'zlib'}:
                raise ValueError(f"Unsupported compression type '{compression_type}' in tess file")

            metadata = header.get('metadata', {})
            frame_offsets = [tuple(offset) for offset in header['frame_offsets']]
            num_frames = len(frame_offsets)
            uncompressed_lengths = header.get('frame_uncompressed_lengths', [])

            # Pre-allocate the result list
            atoms_list = [None] * num_frames
            convert = dict_to_atoms
            
            # Process frames with minimal overhead
            if compression_type == 'zlib' and num_frames:
                def _decompress(offset: Tuple[int, int]) -> bytes:
                    start, length = offset
                    compressed = mm[start:start + length]
                    return zlib.decompress(compressed)

                with ThreadPoolExecutor(max_workers=_max_workers(num_frames)) as executor:
                    decompressed_frames = list(executor.map(_decompress, frame_offsets))

                for i, frame_bytes in enumerate(decompressed_frames):
                    frame_dict = msgpack.unpackb(
                        frame_bytes, raw=False, strict_map_key=False
                    )
                    atoms_list[i] = convert(frame_dict)[0]
            else:
                for i, (start, length) in enumerate(frame_offsets):
                    frame_slice = mm[start:start + length]
                    frame_dict = msgpack.unpackb(
                        frame_slice, raw=False, strict_map_key=False
                    )
                    atoms_list[i] = convert(frame_dict)[0]
        finally:
            mm.close()
    
    return atoms_list, metadata
