import json

from atomict.api import get


def get_ea_exploration(exploration_id: str):
    return get(f"api/ea-exploration/{exploration_id}/")
