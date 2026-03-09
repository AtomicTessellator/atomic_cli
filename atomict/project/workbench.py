from atomict.api import get


def get_workbench(project_id: str) -> dict:
    """Get project workbench details.
    
    Args:
        project_id (str): The project identifier.
    
    Returns:
        dict: The API response payload.
    """
    response = get(f"simulation/workbench/?project_id={project_id}")
    return response
