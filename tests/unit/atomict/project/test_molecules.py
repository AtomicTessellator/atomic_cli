from unittest.mock import patch

import pytest

from atomict.project.molecules import (
    create_project_molecule,
    delete_project_molecule,
    get_project_molecule,
    list_project_molecules,
)


class TestCreateProjectMolecule:
    """Test project molecule creation function"""

    @patch("atomict.project.molecules.post")
    def test_create_project_molecule_basic(self, mock_post):
        """Test basic project molecule creation without name"""
        mock_post.return_value = {
            "id": "mol-123",
            "smiles": "CCO",
            "project": "proj-456",
        }

        result = create_project_molecule(project_id="proj-456", smiles="CCO")

        # Verify API call
        mock_post.assert_called_once_with(
            "api/project-molecule/",
            {
                "project": "proj-456",
                "smiles": "CCO",
            },
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {"id": "mol-123", "smiles": "CCO", "project": "proj-456"}

    @patch("atomict.project.molecules.post")
    def test_create_project_molecule_with_name(self, mock_post):
        """Test project molecule creation with name parameter"""
        mock_post.return_value = {
            "id": "mol-123",
            "smiles": "CCO",
            "name": "Ethanol",
            "project": "proj-456",
        }

        result = create_project_molecule(
            project_id="proj-456", smiles="CCO", name="Ethanol"
        )

        # Verify API call
        mock_post.assert_called_once_with(
            "api/project-molecule/",
            {
                "project": "proj-456",
                "smiles": "CCO",
                "name": "Ethanol",
            },
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {
            "id": "mol-123",
            "smiles": "CCO",
            "name": "Ethanol",
            "project": "proj-456",
        }

    @patch("atomict.project.molecules.post")
    def test_create_project_molecule_complex_smiles(self, mock_post):
        """Test project molecule creation with complex SMILES string"""
        complex_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        mock_post.return_value = {"id": "mol-789", "smiles": complex_smiles}

        result = create_project_molecule(
            project_id="proj-456", smiles=complex_smiles, name="Caffeine"
        )

        # Verify complex SMILES is passed correctly
        mock_post.assert_called_once_with(
            "api/project-molecule/",
            {
                "project": "proj-456",
                "smiles": complex_smiles,
                "name": "Caffeine",
            },
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {"id": "mol-789", "smiles": complex_smiles}


class TestDeleteProjectMolecule:
    """Test project molecule deletion function"""

    @patch("atomict.project.molecules.delete")
    def test_delete_project_molecule(self, mock_delete):
        """Test project molecule deletion"""
        mock_delete.return_value = {"detail": "Successfully deleted"}

        result = delete_project_molecule(molecule_id="mol-123")

        # Verify API call
        mock_delete.assert_called_once_with("api/project-molecule/mol-123/")

        assert result == {"detail": "Successfully deleted"}

    @patch("atomict.project.molecules.delete")
    def test_delete_project_molecule_with_uuid(self, mock_delete):
        """Test project molecule deletion with UUID format"""
        uuid_id = "550e8400-e29b-41d4-a716-446655440000"
        mock_delete.return_value = {}

        result = delete_project_molecule(molecule_id=uuid_id)

        # Verify API call with UUID
        mock_delete.assert_called_once_with(f"api/project-molecule/{uuid_id}/")

        assert result == {}


class TestListProjectMolecules:
    """Test project molecules listing function"""

    @patch("atomict.project.molecules.get")
    def test_list_project_molecules_all(self, mock_get):
        """Test listing all project molecules without filter"""
        mock_response = {
            "results": [
                {"id": "mol-1", "smiles": "CCO", "name": "Ethanol"},
                {"id": "mol-2", "smiles": "CC", "name": "Ethane"},
            ]
        }
        mock_get.return_value = mock_response

        result = list_project_molecules()

        # Verify API call without filter
        mock_get.assert_called_once_with("api/project-molecule/")

        assert result == mock_response

    @patch("atomict.project.molecules.get")
    def test_list_project_molecules_with_project_filter(self, mock_get):
        """Test listing project molecules filtered by project"""
        mock_response = {
            "results": [{"id": "mol-1", "smiles": "CCO", "project": "proj-123"}]
        }
        mock_get.return_value = mock_response

        result = list_project_molecules(project_id="proj-123")

        # Verify API call with project filter
        mock_get.assert_called_once_with("api/project-molecule/?project=proj-123")

        assert result == mock_response

    @patch("atomict.project.molecules.get")
    def test_list_project_molecules_with_none_project_filter(self, mock_get):
        """Test listing project molecules with None project_id (should not filter)"""
        mock_response = {"results": []}
        mock_get.return_value = mock_response

        result = list_project_molecules(project_id=None)

        # Verify API call without filter when project_id is None
        mock_get.assert_called_once_with("api/project-molecule/")

        assert result == mock_response

    @patch("atomict.project.molecules.get")
    def test_list_project_molecules_with_empty_string_project_filter(self, mock_get):
        """Test listing project molecules with empty string project_id (should not filter)"""
        mock_response = {"results": []}
        mock_get.return_value = mock_response

        result = list_project_molecules(project_id="")

        # Verify API call without filter when project_id is empty string
        mock_get.assert_called_once_with("api/project-molecule/")

        assert result == mock_response


class TestGetProjectMolecule:
    """Test project molecule retrieval function"""

    @patch("atomict.project.molecules.get")
    def test_get_project_molecule(self, mock_get):
        """Test retrieving a single project molecule"""
        mock_response = {
            "id": "mol-123",
            "smiles": "CCO",
            "name": "Ethanol",
            "project": "proj-456",
        }
        mock_get.return_value = mock_response

        result = get_project_molecule(molecule_id="mol-123")

        # Verify API call
        mock_get.assert_called_once_with("api/project-molecule/mol-123/")

        assert result == mock_response

    @patch("atomict.project.molecules.get")
    def test_get_project_molecule_with_uuid(self, mock_get):
        """Test retrieving a project molecule with UUID format ID"""
        uuid_id = "550e8400-e29b-41d4-a716-446655440000"
        mock_response = {
            "id": uuid_id,
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "name": "Caffeine",
        }
        mock_get.return_value = mock_response

        result = get_project_molecule(molecule_id=uuid_id)

        # Verify API call with UUID
        mock_get.assert_called_once_with(f"api/project-molecule/{uuid_id}/")

        assert result == mock_response

    @patch("atomict.project.molecules.get")
    def test_get_project_molecule_without_name(self, mock_get):
        """Test retrieving a project molecule that has no name"""
        mock_response = {
            "id": "mol-789",
            "smiles": "CO",
            "name": None,
            "project": "proj-456",
        }
        mock_get.return_value = mock_response

        result = get_project_molecule(molecule_id="mol-789")

        # Verify API call
        mock_get.assert_called_once_with("api/project-molecule/mol-789/")

        assert result == mock_response
