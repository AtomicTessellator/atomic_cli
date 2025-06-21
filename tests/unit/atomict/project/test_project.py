from unittest.mock import patch

import pytest

from atomict.project.project import (
    create_project,
    delete_project,
    get_project,
    get_project_by_name,
    list_projects,
    project_exists,
    update_project,
)


class TestCreateProject:
    """Test project creation function"""

    @patch("atomict.project.project.post")
    def test_create_project_basic(self, mock_post):
        """Test basic project creation with required parameters"""
        mock_post.return_value = {"id": "project-123", "name": "Test Project"}

        result = create_project(name="Test Project")

        mock_post.assert_called_once_with(
            "api/project/",
            {
                "name": "Test Project",
                "description_html": None,
                "thumbnail_smiles": None,
            },
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {"id": "project-123", "name": "Test Project"}

    @patch("atomict.project.project.post")
    def test_create_project_with_all_params(self, mock_post):
        """Test project creation with all parameters"""
        mock_post.return_value = {"id": "project-456", "name": "Full Project"}

        result = create_project(
            name="Full Project",
            description="A detailed description",
            thumbnail_smiles="CCO",
        )

        mock_post.assert_called_once_with(
            "api/project/",
            {
                "name": "Full Project",
                "description_html": "A detailed description",
                "thumbnail_smiles": "CCO",
            },
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {"id": "project-456", "name": "Full Project"}


class TestDeleteProject:
    """Test project deletion function"""

    @patch("atomict.project.project.delete")
    def test_delete_project(self, mock_delete):
        """Test project deletion"""
        mock_delete.return_value = {"success": True}

        result = delete_project("project-123")

        mock_delete.assert_called_once_with("api/project/project-123/")
        assert result == {"success": True}


class TestProjectExists:
    """Test project existence check function"""

    @patch("atomict.project.project.get")
    def test_project_exists_true(self, mock_get):
        """Test project exists returns True when project found"""
        mock_get.return_value = {"count": 1, "results": [{"id": "project-123"}]}

        result = project_exists("Test Project")

        mock_get.assert_called_once_with("api/project/?name=Test Project")
        assert result is True

    @patch("atomict.project.project.get")
    def test_project_exists_false(self, mock_get):
        """Test project exists returns False when no project found"""
        mock_get.return_value = {"count": 0, "results": []}

        result = project_exists("Nonexistent Project")

        mock_get.assert_called_once_with("api/project/?name=Nonexistent Project")
        assert result is False


class TestGetProjectByName:
    """Test get project by name function"""

    @patch("atomict.project.project.get")
    def test_get_project_by_name(self, mock_get):
        """Test get project by name returns first result"""
        expected_project = {"id": "project-123", "name": "Test Project"}
        mock_get.return_value = {"count": 1, "results": [expected_project]}

        result = get_project_by_name("Test Project")

        mock_get.assert_called_once_with("api/project/?name=Test Project")
        assert result == expected_project


class TestGetProject:
    """Test get project by ID function"""

    @patch("atomict.project.project.get")
    def test_get_project(self, mock_get):
        """Test get project by ID"""
        expected_project = {"id": "project-123", "name": "Test Project"}
        mock_get.return_value = expected_project

        result = get_project("project-123")

        mock_get.assert_called_once_with("api/project/project-123/")
        assert result == expected_project


class TestListProjects:
    """Test list projects function"""

    @patch("atomict.project.project.get")
    def test_list_projects_no_params(self, mock_get):
        """Test list projects with no parameters"""
        expected_response = {"count": 2, "results": [{"id": "1"}, {"id": "2"}]}
        mock_get.return_value = expected_response

        result = list_projects()

        mock_get.assert_called_once_with("api/project/")
        assert result == expected_response

    @patch("atomict.project.project.get")
    def test_list_projects_with_search(self, mock_get):
        """Test list projects with search parameter"""
        expected_response = {"count": 1, "results": [{"id": "1", "name": "Test"}]}
        mock_get.return_value = expected_response

        result = list_projects(search="Test")

        mock_get.assert_called_once_with("api/project/?search=Test")
        assert result == expected_response

    @patch("atomict.project.project.get")
    def test_list_projects_with_ordering(self, mock_get):
        """Test list projects with ordering parameter"""
        expected_response = {"count": 2, "results": [{"id": "1"}, {"id": "2"}]}
        mock_get.return_value = expected_response

        result = list_projects(ordering="-created_at")

        mock_get.assert_called_once_with("api/project/?ordering=-created_at")
        assert result == expected_response

    @patch("atomict.project.project.get")
    def test_list_projects_with_filters(self, mock_get):
        """Test list projects with additional filters"""
        expected_response = {"count": 1, "results": [{"id": "1"}]}
        mock_get.return_value = expected_response

        result = list_projects(owner="user123", active=True)

        # URL should contain both filter parameters
        call_args = mock_get.call_args[0][0]
        assert call_args.startswith("api/project/?")
        assert "owner=user123" in call_args
        assert "active=True" in call_args
        assert result == expected_response

    @patch("atomict.project.project.get")
    def test_list_projects_with_all_params(self, mock_get):
        """Test list projects with all parameters combined"""
        expected_response = {"count": 1, "results": [{"id": "1"}]}
        mock_get.return_value = expected_response

        result = list_projects(
            search="Test", ordering="-name", owner="user123", active=True
        )

        call_args = mock_get.call_args[0][0]
        assert call_args.startswith("api/project/?")
        assert "search=Test" in call_args
        assert "ordering=-name" in call_args
        assert "owner=user123" in call_args
        assert "active=True" in call_args
        assert result == expected_response


class TestUpdateProject:
    """Test project update function"""

    @patch("atomict.project.project.patch")
    def test_update_project_name_only(self, mock_patch):
        """Test updating project with name only"""
        mock_patch.return_value = {"id": "project-123", "name": "Updated Name"}

        result = update_project("project-123", name="Updated Name")

        mock_patch.assert_called_once_with(
            "api/project/project-123/",
            {"name": "Updated Name"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": "project-123", "name": "Updated Name"}

    @patch("atomict.project.project.patch")
    def test_update_project_description_only(self, mock_patch):
        """Test updating project with description only"""
        mock_patch.return_value = {"id": "project-123", "description_html": "New desc"}

        result = update_project("project-123", description="New desc")

        mock_patch.assert_called_once_with(
            "api/project/project-123/",
            {"description_html": "New desc"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": "project-123", "description_html": "New desc"}

    @patch("atomict.project.project.patch")
    def test_update_project_thumbnail_only(self, mock_patch):
        """Test updating project with thumbnail_smiles only"""
        mock_patch.return_value = {"id": "project-123", "thumbnail_smiles": "CCO"}

        result = update_project("project-123", thumbnail_smiles="CCO")

        mock_patch.assert_called_once_with(
            "api/project/project-123/",
            {"thumbnail_smiles": "CCO"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": "project-123", "thumbnail_smiles": "CCO"}

    @patch("atomict.project.project.patch")
    def test_update_project_all_params(self, mock_patch):
        """Test updating project with all parameters"""
        mock_patch.return_value = {
            "id": "project-123",
            "name": "New Name",
            "description_html": "New Description",
            "thumbnail_smiles": "CCO",
        }

        result = update_project(
            "project-123",
            name="New Name",
            description="New Description",
            thumbnail_smiles="CCO",
        )

        mock_patch.assert_called_once_with(
            "api/project/project-123/",
            {
                "name": "New Name",
                "description_html": "New Description",
                "thumbnail_smiles": "CCO",
            },
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {
            "id": "project-123",
            "name": "New Name",
            "description_html": "New Description",
            "thumbnail_smiles": "CCO",
        }

    @patch("atomict.project.project.patch")
    def test_update_project_no_params(self, mock_patch):
        """Test updating project with no optional parameters (empty payload)"""
        mock_patch.return_value = {"id": "project-123"}

        result = update_project("project-123")

        mock_patch.assert_called_once_with(
            "api/project/project-123/",
            {},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {"id": "project-123"}
