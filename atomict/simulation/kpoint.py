from atomict.api import get, post, patch


def get_kpoint_exploration(simulation_id: str, *, api_root: str = None, token: str = None):
    """
    Get kpoints for a simulation
    """
    result = get(f"api/kpoint-exploration/{simulation_id}/", api_root=api_root, token=token)
    return result


def update_kpoint_exploration(exploration_id: str, fields: dict, *, api_root: str = None, token: str = None):
    """
    Update KPoint exploration
    """
    result = patch(f"api/kpoint-exploration/{exploration_id}/", payload=fields, api_root=api_root, token=token)
    return result


def get_kpoint_simulation_list(exploration_id: str, *, api_root: str = None, token: str = None):
    """
    Get kpoints for a simulation
    """
    result = get(f"api/kpoint-simulation/?exploration__id={exploration_id}", api_root=api_root, token=token)
    return result


def get_kpoint_simulation(simulation_id: str, *, api_root: str = None, token: str = None):
    """
    Get KPoint simulation
    """
    result = get(f"api/kpoint-simulation/{simulation_id}/", api_root=api_root, token=token)
    return result


def create_kpoint_simulation(
    exploration_id: str,
    simulation_id: str,
    k_points: list[float],
    *,
    api_root: str = None,
    token: str = None,
):
    """
    Create KPoint simulation
    """
    result = post(
        "api/kpoint-simulation/",
        payload={
            "exploration_id": exploration_id,
            "simulation_id": simulation_id,
            "k_points": k_points,
        },
        api_root=api_root,
        token=token,
    )
    return result


def update_kpoint_simulation(simulation_id: str, fields: dict, *, api_root: str = None, token: str = None):
    """
    Update KPoint simulation
    """
    result = patch(f"api/kpoint-simulation/{simulation_id}/", payload=fields, api_root=api_root, token=token)
    return result


def get_kpoint_analysis(analysis_id: str, *, api_root: str = None, token: str = None):
    """
    Get KPoint analysis
    """
    result = get(f"api/kpoint-analysis/{analysis_id}/", api_root=api_root, token=token)
    return result


def update_kpoint_analysis(analysis_id: str, fields: dict, *, api_root: str = None, token: str = None):
    """
    Update KPoint analysis
    """
    result = patch(f"api/kpoint-analysis/{analysis_id}/", payload=fields, api_root=api_root, token=token)
    return result


def delete_kpoint_simulations(exploration_id: str, *, api_root: str = None, token: str = None):
    """Delete all FHIAims simulations linked to the exploration via KPointSimulation.

    Calls the server-side reset_simulations endpoint which cascade-deletes
    KPointSimulation records by deleting their FHIAimsSimulation parents.
    Safe to call when no simulations exist (returns {"deleted": 0}).
    Used as a saga compensation and as a clean-slate step before retrying Stage 1.
    """
    result = post(
        f"api/kpoint-exploration/{exploration_id}/reset_simulations/",
        payload={},
        api_root=api_root,
        token=token,
    )
    return result
