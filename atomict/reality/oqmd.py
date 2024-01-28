from atomict.api import get


def get_bulk(bulk_id: str):
    return get(f"oqmd/structure/{bulk_id}/")
