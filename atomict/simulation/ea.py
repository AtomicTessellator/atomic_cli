import json

from atomict.api import get


def get_ea_exploration(exploration_id: str):
    return get(f"api/ea-exploration/{exploration_id}/")


def get_ea_exploration_sample(sample_id: str):
    return get(f"api/ea-exploration-sample/{sample_id}/")


def get_ea_exploration_samples(exploration_id: str):
    return get(f"api/ea-exploration-sample/?exploration={exploration_id}")


def get_ea_exploration_analysis(analysis_id: str):
    return get(f"api/ea-exploration-analysis/{analysis_id}/")


def get_ea_exploration_analysis_file(analysis_file_id: str):
    return get(f"api/ea-exploration-analysis-file/{analysis_file_id}/")
