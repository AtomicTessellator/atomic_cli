"""Unit tests for project stars operations"""

import unittest
from unittest.mock import patch, MagicMock
import pytest

from atomict.project.stars import (
    create_project_star,
    delete_project_star,
    list_project_stars,
    get_project_star,
)


class TestProjectStars(unittest.TestCase):
    """Test project stars functions"""

    @patch("atomict.project.stars.post")
    def test_create_project_star_success(self, mock_post):
        """Test successful project star creation"""
        mock_post.return_value = {
            "id": "star-123",
            "project": "proj-456",
            "user": "user-789",
            "created_at": "2024-01-01T00:00:00Z",
        }

        result = create_project_star(project_id="proj-456")

        mock_post.assert_called_once_with(
            "api/project-star/",
            {"project": "proj-456"},
            extra_headers={"Content-Type": "application/json"},
        )
        assert result == {
            "id": "star-123",
            "project": "proj-456",
            "user": "user-789",
            "created_at": "2024-01-01T00:00:00Z",
        }

    @patch("atomict.project.stars.delete")
    def test_delete_project_star_success(self, mock_delete):
        """Test successful project star deletion"""
        mock_delete.return_value = {"status": "deleted"}

        result = delete_project_star(star_id="star-123")

        mock_delete.assert_called_once_with("api/project-star/star-123/")
        assert result == {"status": "deleted"}

    @patch("atomict.project.stars.get")
    def test_list_project_stars_success(self, mock_get):
        """Test successful project stars listing"""
        mock_get.return_value = {
            "count": 2,
            "results": [
                {
                    "id": "star-123",
                    "project": "proj-456",
                    "user": "user-789",
                    "created_at": "2024-01-01T00:00:00Z",
                },
                {
                    "id": "star-124",
                    "project": "proj-457",
                    "user": "user-789",
                    "created_at": "2024-01-02T00:00:00Z",
                },
            ],
        }

        result = list_project_stars()

        mock_get.assert_called_once_with("api/project-star/")
        assert result["count"] == 2
        assert len(result["results"]) == 2
        assert result["results"][0]["id"] == "star-123"

    @patch("atomict.project.stars.get")
    def test_get_project_star_success(self, mock_get):
        """Test successful project star retrieval"""
        mock_get.return_value = {
            "id": "star-123",
            "project": "proj-456",
            "user": "user-789",
            "created_at": "2024-01-01T00:00:00Z",
        }

        result = get_project_star(star_id="star-123")

        mock_get.assert_called_once_with("api/project-star/star-123/")
        assert result == {
            "id": "star-123",
            "project": "proj-456",
            "user": "user-789",
            "created_at": "2024-01-01T00:00:00Z",
        }

    @patch("atomict.project.stars.post")
    def test_create_project_star_payload_mapping(self, mock_post):
        """Test that project_id is correctly mapped to 'project' field"""
        mock_post.return_value = {"id": "star-123"}

        create_project_star(project_id="test-project-123")

        # Verify the payload uses 'project' field, not 'project_id'
        expected_payload = {"project": "test-project-123"}
        mock_post.assert_called_once_with(
            "api/project-star/",
            expected_payload,
            extra_headers={"Content-Type": "application/json"},
        )

    @patch("atomict.project.stars.get")
    def test_list_project_stars_empty_response(self, mock_get):
        """Test listing project stars with empty response"""
        mock_get.return_value = {"count": 0, "results": []}

        result = list_project_stars()

        mock_get.assert_called_once_with("api/project-star/")
        assert result["count"] == 0
        assert result["results"] == []

    @patch("atomict.project.stars.post")
    def test_create_project_star_api_error(self, mock_post):
        """Test create project star with API error"""
        mock_post.side_effect = Exception("API Error")

        with pytest.raises(Exception, match="API Error"):
            create_project_star(project_id="proj-456")

    @patch("atomict.project.stars.delete")
    def test_delete_project_star_api_error(self, mock_delete):
        """Test delete project star with API error"""
        mock_delete.side_effect = Exception("Not Found")

        with pytest.raises(Exception, match="Not Found"):
            delete_project_star(star_id="nonexistent-star")

    @patch("atomict.project.stars.get")
    def test_get_project_star_api_error(self, mock_get):
        """Test get project star with API error"""
        mock_get.side_effect = Exception("Not Found")

        with pytest.raises(Exception, match="Not Found"):
            get_project_star(star_id="nonexistent-star")


if __name__ == "__main__":
    unittest.main()
