from atomict.api import get


def get_workbench(project_id: str) -> dict:
    response = get(f"simulation/workbench/?project_id={project_id}")
    return response
