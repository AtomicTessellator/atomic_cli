"""Unit tests for project note operations"""

import unittest
from unittest.mock import patch, MagicMock
import pytest

from atomict.project.note import (
    create_project_note,
    delete_project_note,
    get_project_note,
    list_project_notes,
    update_project_note,
)


class TestCreateProjectNote(unittest.TestCase):
    """Test create_project_note function"""

    @patch("atomict.project.note.post")
    def test_create_project_note_minimal(self, mock_post):
        """Test project note creation with minimal parameters"""
        mock_post.return_value = {"id": "note-123", "title": "Test Note"}

        result = create_project_note(
            project_id="proj-456",
            title="Test Note",
            content="<p>Test content</p>"
        )

        mock_post.assert_called_once_with(
            "api/project-note/",
            {
                "project": "proj-456",
                "title": "Test Note", 
                "content_html": "<p>Test content</p>",
                "show_description": True,
            },
            extra_headers={"Content-Type": "application/json"}
        )
        assert result == {"id": "note-123", "title": "Test Note"}

    @patch("atomict.project.note.post")
    def test_create_project_note_with_show_description_false(self, mock_post):
        """Test project note creation with show_description disabled"""
        mock_post.return_value = {"id": "note-123", "title": "Test Note"}

        result = create_project_note(
            project_id="proj-456",
            title="Test Note", 
            content="<p>Test content</p>",
            show_description=False
        )

        mock_post.assert_called_once_with(
            "api/project-note/",
            {
                "project": "proj-456",
                "title": "Test Note",
                "content_html": "<p>Test content</p>",
                "show_description": False,
            },
            extra_headers={"Content-Type": "application/json"}
        )
        assert result == {"id": "note-123", "title": "Test Note"}


class TestDeleteProjectNote(unittest.TestCase):
    """Test delete_project_note function"""

    @patch("atomict.project.note.delete")
    def test_delete_project_note(self, mock_delete):
        """Test project note deletion"""
        mock_delete.return_value = {"deleted": True}

        result = delete_project_note("note-123")

        mock_delete.assert_called_once_with("api/project-note/note-123/")
        assert result == {"deleted": True}


class TestGetProjectNote(unittest.TestCase):
    """Test get_project_note function"""

    @patch("atomict.project.note.get")
    def test_get_project_note(self, mock_get):
        """Test project note retrieval"""
        expected_note = {
            "id": "note-123",
            "title": "Test Note",
            "content_html": "<p>Test content</p>",
            "show_description": True
        }
        mock_get.return_value = expected_note

        result = get_project_note("note-123")

        mock_get.assert_called_once_with("api/project-note/note-123/")
        assert result == expected_note


class TestListProjectNotes(unittest.TestCase):
    """Test list_project_notes function"""

    @patch("atomict.project.note.get")
    def test_list_project_notes_no_filter(self, mock_get):
        """Test listing all project notes without project filter"""
        expected_notes = {
            "results": [
                {"id": "note-1", "title": "Note 1"},
                {"id": "note-2", "title": "Note 2"}
            ]
        }
        mock_get.return_value = expected_notes

        result = list_project_notes()

        mock_get.assert_called_once_with("api/project-note/")
        assert result == expected_notes

    @patch("atomict.project.note.get")
    def test_list_project_notes_with_project_filter(self, mock_get):
        """Test listing project notes filtered by project"""
        expected_notes = {
            "results": [
                {"id": "note-1", "title": "Project Note 1", "project": "proj-456"}
            ]
        }
        mock_get.return_value = expected_notes

        result = list_project_notes(project_id="proj-456")

        mock_get.assert_called_once_with("api/project-note/?project=proj-456")
        assert result == expected_notes

    @patch("atomict.project.note.get")
    def test_list_project_notes_with_none_project_id(self, mock_get):
        """Test listing project notes with None project_id (should not add filter)"""
        expected_notes = {"results": []}
        mock_get.return_value = expected_notes

        result = list_project_notes(project_id=None)

        mock_get.assert_called_once_with("api/project-note/")
        assert result == expected_notes


class TestUpdateProjectNote(unittest.TestCase):
    """Test update_project_note function"""

    @patch("atomict.project.note.patch")
    def test_update_project_note_all_fields(self, mock_patch):
        """Test project note update with all fields"""
        mock_patch.return_value = {"id": "note-123", "title": "Updated Note"}

        result = update_project_note(
            note_id="note-123",
            title="Updated Note",
            content="<p>Updated content</p>",
            show_description=False
        )

        mock_patch.assert_called_once_with(
            "api/project-note/note-123/",
            {
                "title": "Updated Note",
                "content_html": "<p>Updated content</p>",
                "show_description": False,
            }
        )
        assert result == {"id": "note-123", "title": "Updated Note"}

    @patch("atomict.project.note.patch")
    def test_update_project_note_title_only(self, mock_patch):
        """Test project note update with title only"""
        mock_patch.return_value = {"id": "note-123", "title": "New Title"}

        result = update_project_note(note_id="note-123", title="New Title")

        mock_patch.assert_called_once_with(
            "api/project-note/note-123/",
            {"title": "New Title"}
        )
        assert result == {"id": "note-123", "title": "New Title"}

    @patch("atomict.project.note.patch")
    def test_update_project_note_content_only(self, mock_patch):
        """Test project note update with content only"""
        mock_patch.return_value = {"id": "note-123", "content_html": "<p>New content</p>"}

        result = update_project_note(note_id="note-123", content="<p>New content</p>")

        mock_patch.assert_called_once_with(
            "api/project-note/note-123/",
            {"content_html": "<p>New content</p>"}
        )
        assert result == {"id": "note-123", "content_html": "<p>New content</p>"}

    @patch("atomict.project.note.patch")
    def test_update_project_note_show_description_only(self, mock_patch):
        """Test project note update with show_description only"""
        mock_patch.return_value = {"id": "note-123", "show_description": True}

        result = update_project_note(note_id="note-123", show_description=True)

        mock_patch.assert_called_once_with(
            "api/project-note/note-123/",
            {"show_description": True}
        )
        assert result == {"id": "note-123", "show_description": True}

    @patch("atomict.project.note.patch")
    def test_update_project_note_no_fields(self, mock_patch):
        """Test project note update with no fields (empty payload)"""
        mock_patch.return_value = {"id": "note-123"}

        result = update_project_note(note_id="note-123")

        mock_patch.assert_called_once_with(
            "api/project-note/note-123/",
            {}
        )
        assert result == {"id": "note-123"}

    @patch("atomict.project.note.patch")
    def test_update_project_note_content_html_field_mapping(self, mock_patch):
        """Test that content parameter maps to content_html field"""
        mock_patch.return_value = {"id": "note-123"}

        update_project_note(note_id="note-123", content="<h1>Test</h1>")

        # Verify the field mapping: content -> content_html
        mock_patch.assert_called_once_with(
            "api/project-note/note-123/",
            {"content_html": "<h1>Test</h1>"}
        )


if __name__ == "__main__":
    unittest.main()
