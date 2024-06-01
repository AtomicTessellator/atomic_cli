from atomict.api import get, post


def get_qe_simulation(simulation_id: str):
    """
    Get a Quantum Espresso simulation
    """
    result = get(f"api/qe-simulation/{simulation_id}/")
    return result


def associate_user_upload_with_qe_simulation(user_upload_id: str, qe_simulation_id: str):
    """
    Associate a user upload with a Quantum Espresso simulation
    """
    result = post(
        "api/qe-simulation-file/",
        payload={
            "user_upload": user_upload_id,
            "simulation": qe_simulation_id
        }    
    )
    return result
