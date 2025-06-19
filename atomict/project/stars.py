from atomict.api import delete, get, post


def create_project_star(project_id: str) -> dict:
    """Star a project for the current user.
    
    Args:
        project_id: ID of the project to star
        
    Returns:
        dict: Created star record
    """
    payload = {
        "project": project_id,
    }

    response = post(
        "api/project-star/", payload, extra_headers={"Content-Type": "application/json"}
    )
    return response


def delete_project_star(star_id: str) -> dict:
    """Unstar a project by removing the star record.
    
    Args:
        star_id: ID of the star record to delete
        
    Returns:
        dict: Response from deletion
    """
    response = delete(f"api/project-star/{star_id}/")
    return response


def list_project_stars() -> dict:
    """List all starred projects for the current user.
    
    Returns:
        dict: List of star records
    """
    response = get("api/project-star/")
    return response


def get_project_star(star_id: str) -> dict:
    """Get a single star record by ID.
    
    Args:
        star_id: ID of the star record to retrieve
        
    Returns:
        dict: Star record details
    """
    response = get(f"api/project-star/{star_id}/")
    return response
