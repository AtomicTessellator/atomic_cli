import os
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
    compression: Optional[str] = 'zlib',
    compression_level: int = 1,
) -> None:

    from atomict.io.msgpack import atoms_to_dict

    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
        import zlib
        import lz4.block
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

    with open(filename, 'wb') as f:
        frame_offsets: List[Tuple[int, int]] = []
        uncompressed_lengths: List[int] = []
        if compression_mode == 'zlib' or compression_mode == 'lz4':
            worker_count = _max_workers(len(frames_list))
            num_tasks = worker_count * 4
            group_size = max(1, (len(frames_list) + num_tasks - 1) // num_tasks)
            groups = [frames_list[i:i + group_size] for i in range(0, len(frames_list), group_size)]

            def _process_group(group: List['Atoms']) -> Tuple[bytes, List[int], List[int]]:
                group_compressed = []
                group_ulens = []
                for atom_frame in group:
                    frame_dict = atoms_to_dict([atom_frame])
                    packed = msgpack.packb(frame_dict, use_bin_type=True)
                    if compression_mode == 'zlib':
                        compressed = zlib.compress(packed, compression_level)
                    else:  # lz4
                        compressed = lz4.block.compress(packed, store_size=True)
                    group_compressed.append(compressed)
                    group_ulens.append(len(packed))
                group_data = b''.join(group_compressed)
                group_lengths = [len(c) for c in group_compressed]
                return group_data, group_ulens, group_lengths

            initial = f.tell()
            current = initial

            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                for group_data, group_ulens, group_lengths in executor.map(_process_group, groups):
                    f.write(group_data)
                    for l, u in zip(group_lengths, group_ulens):
                        frame_offsets.append((current, l))
                        uncompressed_lengths.append(u)
                        current += l

        else:
            compressed_list: List[bytes] = []
            for atom_frame in frames_list:
                frame_dict = atoms_to_dict([atom_frame])
                packed = msgpack.packb(frame_dict, use_bin_type=True)
                compressed_list.append(packed)

            # Compute offsets
            initial = f.tell()
            current = initial
            for cb in compressed_list:
                l = len(cb)
                frame_offsets.append((current, l))
                current += l

            # Write all frames
            f.write(b''.join(compressed_list))

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
            }
            if compression_mode == 'zlib':
                header_dict['compression']['level'] = compression_level
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
            if compression_type not in {'none', 'zlib', 'lz4'}:
                raise ValueError(f"Unsupported compression type '{compression_type}' in tess file")

            metadata = header.get('metadata', {})
            frame_offsets = [tuple(offset) for offset in header['frame_offsets']]
            num_frames = len(frame_offsets)

            # Pre-allocate the result list
            atoms_list = [None] * num_frames
            convert = dict_to_atoms
            
            # Process frames with minimal overhead
            if compression_type == 'zlib' and num_frames:

                def _decompress(offset: Tuple[int, int]) -> bytes:
                    start, length = offset
                    compressed = mm[start:start + length]
                    return zlib.decompress(compressed)

                worker_count = _max_workers(num_frames)
                chunk_size = _chunk_size(worker_count)

                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    for chunk_start in range(0, num_frames, chunk_size):
                        offset_chunk = frame_offsets[chunk_start:chunk_start + chunk_size]
                        for i, frame_bytes in enumerate(
                            executor.map(_decompress, offset_chunk), start=chunk_start
                        ):
                            frame_dict = msgpack.unpackb(
                                frame_bytes, raw=False, strict_map_key=False
                            )
                            atoms_list[i] = convert(frame_dict)[0]
            elif compression_type == 'lz4' and num_frames:

                def _decompress(offset: Tuple[int, int]) -> bytes:
                    start, length = offset
                    compressed = mm[start:start + length]
                    return lz4.block.decompress(compressed)

                worker_count = _max_workers(num_frames)
                chunk_size = _chunk_size(worker_count)

                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    for chunk_start in range(0, num_frames, chunk_size):
                        offset_chunk = frame_offsets[chunk_start:chunk_start + chunk_size]
                        for i, frame_bytes in enumerate(
                            executor.map(_decompress, offset_chunk), start=chunk_start
                        ):
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
