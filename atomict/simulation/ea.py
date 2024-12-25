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



def associate_user_upload_with_ea_exploration(user_upload_id: str, analysis_id: str):
    return post(
        "api/ea-exploration-analysis-file/",
        payload={"user_upload": user_upload_id, "analysis": analysis_id},
    )


def create_exploration_sample(
    exploration_id: str,
    simulation_id: str = None,
    mlrelax_id: str = None,
    strain: float = None,
    matrix: int = None,
):
    """
    Create an exploration sample

    exploration_id: str - The ID of the exploration to associate the sample with
    simulation_id: str - The ID of the simulation to associate with the exploration
    strain: float - The strain to associate with the sample
    matrix: int - The matrix to associate with the sample
    """

    if simulation_id is None and mlrelax_id is None:
        raise ValueError("Either simulation_id or mlrelax_id must be provided")

    payload = {
        "exploration": exploration_id,
        "strain": strain,
        "matrix": matrix,
    }

    if simulation_id:
        payload["simulation"] = simulation_id
    elif mlrelax_id:
        payload["mlrelax"] = mlrelax_id

    return post(
        "api/ea-exploration-sample/",
        payload=payload,
    )
