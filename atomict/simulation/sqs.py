from atomict.api import get


def get_simulation(simulation_id: str, full: bool = False):
    """
    Get a SQS simulation
    """
    if full:
        result = get(f"api/sqs-exploration/{simulation_id}/?full=true")
    else:
        result = get(f"api/sqs-exploration/{simulation_id}/")

    return result
