from ase.io.aims import check_convergence, get_aims_out_chunks, get_header_chunk
from ase.utils import reader


# We have to write our own function because he read_aims_output function inside ASE only supports
# reading a single image from the output file. We want to support reading all images.
@reader
def read_aims_output(fd, index=None, non_convergence_ok=False):
    """Import FHI-aims output files with all data available, i.e.
    relaxations, MD information, force information etc etc etc."""
    header_chunk = get_header_chunk(fd)
    chunks = list(get_aims_out_chunks(fd, header_chunk))
    check_convergence(chunks, non_convergence_ok)

    # Relaxations have an additional footer chunk due to how it is split
    if header_chunk.is_relaxation:
        images = [chunk.atoms for chunk in chunks[:-1]]
    else:
        images = [chunk.atoms for chunk in chunks]

    if index is None:
        return images
    else:
        return images[index]
