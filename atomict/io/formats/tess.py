from typing import Union, List, Dict


def write_tess(atoms: Union['ase.Atoms', List['ase.Atoms']], filename: str, metadata: Dict = None):

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
            frame_dict = atoms_to_dict([atom_frame], selective=False)
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


def read_tess(filename: str) -> tuple[List['ase.Atoms'], Dict]:

    try:
        import msgpack
        import msgpack_numpy as m
        from atomict.io.atoms import dict_to_atoms
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[utils]` to use msgpack I/O")

    # Enable numpy array deserialization
    m.patch()
    
    with open(filename, 'rb') as f:
        f.seek(0, 2)
        file_size = f.tell()
        f.seek(file_size - 8)
        header_start = int.from_bytes(f.read(8), 'little')
        f.seek(header_start)
        header_bytes = f.read(file_size - header_start - 8)
        header = msgpack.unpackb(header_bytes, raw=False)
        
        if header['format_version'] != 2:
            raise ValueError("Invalid format version")
        
        metadata = header['metadata']
        frame_offsets = header['frame_offsets']
        
        atoms_list = []
        for start, length in frame_offsets:
            f.seek(start)
            frame_bytes = f.read(length)
            frame_dict = msgpack.unpackb(frame_bytes, raw=False)
            atoms = dict_to_atoms(frame_dict)
            atoms_list.append(atoms[0])
    
    return atoms_list, metadata
