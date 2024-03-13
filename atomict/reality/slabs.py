from atomict.api import get


def get_slab(slab_id: str):
    return get(f"simulation/slabs/download/{slab_id}/")
