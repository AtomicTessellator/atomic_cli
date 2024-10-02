from atomict.api import get, post


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


def create_exploration_sample(
    exploration_id: str, simulation_id: str, strain: float = None, matrix: int = None
):
    """
    Create an exploration sample

    exploration_id: str - The ID of the exploration to associate the sample with
    simulation_id: str - The ID of the simulation to associate with the exploration
    strain: float - The strain to associate with the sample
    matrix: int - The matrix to associate with the sample
    """

    result = post(
        "api/ea-exploration-sample/",
        payload={
            "exploration": exploration_id,
            "simulation": simulation_id,
            "strain": strain,
            "matrix": matrix,
        },
    )

    return result
