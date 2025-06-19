from atomict.api import delete, get, patch, post


def create_project(name: str, description: str = None, thumbnail_smiles: str = None) -> dict:

    payload = {
        "name": name,
        "description_html": description,
        "thumbnail_smiles": thumbnail_smiles,
    }

    response = post(
        "api/project/", payload, extra_headers={"Content-Type": "application/json"})
    return response


def delete_project(project_id: str) -> dict:
    response = delete(f"api/project/{project_id}/")

    return response


def project_exists(name: str) -> bool:
    response = get("api/project/")
    
    # Filter results by exact name match since backend doesn't support name filtering
    for project in response.get('results', []):
        if project.get('name') == name:
            return True
    return False


def get_project_by_name(name: str) -> dict:
    response = get("api/project/")
    
    # Filter results by exact name match since backend doesn't support name filtering
    for project in response.get('results', []):
        if project.get('name') == name:
            return project
    
    raise ValueError(f"Project with name '{name}' not found")


def get_project(project_id: str) -> dict:
    response = get(f"api/project/{project_id}/")
    return response


def list_projects(search: str = None, ordering: str = None, **filters) -> dict:
    params = {}
    if search is not None:
        params['search'] = search
    if ordering is not None:
        params['ordering'] = ordering
    params.update(filters)
    
    query_string = '&'.join([f"{k}={v}" for k, v in params.items()])
    url = "api/project/"
    if query_string:
        url += f"?{query_string}"
    
    response = get(url)
    return response


def update_project(project_id: str, name: str = None, description: str = None, thumbnail_smiles: str = None) -> dict:
    payload = {}
    if name is not None:
        payload["name"] = name
    if description is not None:
        payload["description_html"] = description
    if thumbnail_smiles is not None:
        payload["thumbnail_smiles"] = thumbnail_smiles
    
    response = patch(f"api/project/{project_id}/", payload)
    return response
