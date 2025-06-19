from atomict.api import delete, get, post


def create_project_molecule(project_id: str, smiles: str, name: str = None) -> dict:
    """Create a new project molecule.
    
    Args:
        project_id: ID of the project
        smiles: SMILES string for the molecule
        name: Optional name/heading for the molecule
        
    Returns:
        API response dict containing the created molecule
    """
    payload = {
        "project": project_id,
        "smiles": smiles,
    }
    
    if name is not None:
        payload["heading"] = name  # Backend uses 'heading' field

    response = post("api/project-molecule/", payload)
    return response


def delete_project_molecule(molecule_id: str) -> dict:
    """Delete a project molecule.
    
    Args:
        molecule_id: ID of the molecule to delete
        
    Returns:
        API response dict
    """
    response = delete(f"api/project-molecule/{molecule_id}/")
    return response


def list_project_molecules(project_id: str = None) -> dict:
    """List project molecules, optionally filtered by project.
    
    Args:
        project_id: Optional project ID to filter by
        
    Returns:
        API response dict containing list of molecules
    """
    url = "api/project-molecule/"
    if project_id:
        url += f"?project__id={project_id}"
    response = get(url)
    return response


def get_project_molecule(molecule_id: str) -> dict:
    """Get a single project molecule.
    
    Args:
        molecule_id: ID of the molecule to retrieve
        
    Returns:
        API response dict containing the molecule
    """
    response = get(f"api/project-molecule/{molecule_id}/")
    return response
