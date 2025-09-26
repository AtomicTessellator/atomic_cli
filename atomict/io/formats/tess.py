from typing import Union, List, Dict, TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


def write_tess(atoms: Union['Atoms', List['Atoms']], filename: str, metadata: Dict = None):

    from atomict.io.msgpack import atoms_to_dict

    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[utils]` to use msgpack I/O")

    # Enable numpy array serialization
    m.patch()

    if isinstance(atoms, Atoms):
        frames_list = [atoms]
    else:
        frames_list = atoms
    
    # Collect global unique_symbols
    unique_symbols_set = set()
    for a in frames_list:
        unique_symbols_set.update(a.get_chemical_symbols())
    unique_symbols = sorted(list(unique_symbols_set))
    
    # Write frames
    frame_offsets = []
    with open(filename, 'wb') as f:
        for atom_frame in frames_list:
            frame_dict = atoms_to_dict([atom_frame])
            frame_bytes = msgpack.packb(frame_dict, use_bin_type=True)
            start = f.tell()
            f.write(frame_bytes)
            frame_offsets.append((start, len(frame_bytes)))
        
        # Write header
        header_dict = {
            'format_version': 2,
            'metadata': metadata or {},
            'unique_symbols': unique_symbols,
            'num_frames': len(frames_list),
            'frame_offsets': frame_offsets
        }
        header_bytes = msgpack.packb(header_dict, use_bin_type=True)
        header_start = f.tell()
        f.write(header_bytes)
        
        # Write header offset
        f.write(header_start.to_bytes(8, 'little'))


def read_tess(filename: str) -> tuple[List['Atoms'], Dict]:

    try:
        import msgpack
        import msgpack_numpy as m
        from atomict.io.atoms import dict_to_atoms
        import mmap
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

            if header['format_version'] != 2:
                raise ValueError("Invalid format version")

            metadata = header['metadata']
            frame_offsets = header['frame_offsets']
            num_frames = len(frame_offsets)

            # Pre-allocate the result list
            atoms_list = [None] * num_frames
            convert = dict_to_atoms
            
            # Create a reusable Unpacker for better performance
            unpacker = msgpack.Unpacker(raw=False, strict_map_key=False)
            
            # Process frames with minimal overhead
            for i, (start, length) in enumerate(frame_offsets):
                # Direct slice without creating intermediate bytes object
                unpacker.feed(mm[start:start + length])
                frame_dict = next(unpacker)
                atoms_list[i] = convert(frame_dict)[0]
        finally:
            mm.close()
    
    return atoms_list, metadata
