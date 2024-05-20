from atomict.api import post


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
