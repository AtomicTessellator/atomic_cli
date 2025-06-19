"""
Integration tests for project notes functionality.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)

To run: PYTHONPATH=/path/to/atomic_cli uv run pytest tests/integration/test_project_notes.py -v -m integration
"""

import os
import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.project.note import (
    create_project_note,
    delete_project_note,
    get_project_note,
    list_project_notes,
    update_project_note,
)


@pytest.fixture(scope="session", autouse=True)
def setup_authentication():
    """Setup authentication for integration tests"""
    load_dotenv()
    username = os.getenv("ATOMICT_USERNAME")
    password = os.getenv("ATOMICT_PASSWORD")

    if not username or not password:
        pytest.skip("ATOMICT_USERNAME and ATOMICT_PASSWORD must be set in .env file")

    try:
        token = authenticate(username, password)
        os.environ["AT_TOKEN"] = token
        return token
    except Exception as e:
        pytest.skip(f"Authentication failed: {e}")


@pytest.fixture(scope="session")
def test_project_id():
    """Get test project ID from environment"""
    project_id = os.getenv("TEST_PROJECT_ID")
    if not project_id:
        pytest.skip("TEST_PROJECT_ID environment variable must be set")
    return project_id


@pytest.fixture
def cleanup_project_notes():
    """Track created project notes for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for note_id in created_ids:
        try:
            delete_project_note(note_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestProjectNoteIntegration:
    """Integration tests for project note CRUD operations"""

    def test_create_project_note_basic(
        self, test_project_id, cleanup_project_notes
    ):
        """Test creating a basic project note"""
        result = create_project_note(
            project_id=test_project_id,
            title="Integration Test Note",
            content="<p>This is a test note created by integration tests.</p>",
        )

        # Track for cleanup
        cleanup_project_notes.append(result["id"])

        # Verify response structure
        assert "id" in result
        assert result["title"] == "Integration Test Note"
        assert result["content_html"] == "<p>This is a test note created by integration tests.</p>"
        assert result["show_description"] is True  # default value
        assert result["project"] == test_project_id

        # Verify timestamps
        assert "created_at" in result
        assert "updated_at" in result

    def test_create_project_note_with_show_description_false(
        self, test_project_id, cleanup_project_notes
    ):
        """Test creating a project note with show_description=False (backend overrides to True when content exists)"""
        result = create_project_note(
            project_id=test_project_id,
            title="Note with Content",
            content="<p>This note has content, so backend sets show_description=True.</p>",
            show_description=False,
        )

        # Track for cleanup
        cleanup_project_notes.append(result["id"])

        # Backend automatically sets show_description=True when content is not empty
        assert result["show_description"] is True
        assert result["title"] == "Note with Content"

    def test_get_project_note(self, test_project_id, cleanup_project_notes):
        """Test retrieving a project note"""
        # Create a note first
        created_note = create_project_note(
            project_id=test_project_id,
            title="Get Test Note",
            content="<p>Content for get test.</p>",
        )

        note_id = created_note["id"]
        cleanup_project_notes.append(note_id)

        # Get the note
        retrieved_note = get_project_note(note_id)

        assert retrieved_note["id"] == note_id
        assert retrieved_note["title"] == "Get Test Note"
        assert retrieved_note["content_html"] == "<p>Content for get test.</p>"
        assert retrieved_note["project"] == test_project_id
        assert "created_at" in retrieved_note
        assert "updated_at" in retrieved_note

    def test_list_project_notes_with_filter(
        self, test_project_id, cleanup_project_notes
    ):
        """Test listing project notes filtered by project"""
        # Create multiple notes
        note1 = create_project_note(
            project_id=test_project_id,
            title="List Test Note 1",
            content="<p>First note for list test.</p>",
        )
        note2 = create_project_note(
            project_id=test_project_id,
            title="List Test Note 2", 
            content="<p>Second note for list test.</p>",
        )

        # Track for cleanup
        cleanup_project_notes.extend([note1["id"], note2["id"]])

        # List notes for this project
        notes_response = list_project_notes(project_id=test_project_id)

        # Handle paginated response format
        if isinstance(notes_response, dict) and "results" in notes_response:
            notes = notes_response["results"]
        else:
            notes = notes_response

        # Find our created notes
        created_note_ids = {note1["id"], note2["id"]}
        found_notes = [note for note in notes if note["id"] in created_note_ids]

        assert len(found_notes) == 2
        assert {note["title"] for note in found_notes} == {"List Test Note 1", "List Test Note 2"}

    def test_list_project_notes_all(self):
        """Test listing all project notes without filter"""
        notes_response = list_project_notes()

        # Should return some response (could be empty or have notes)
        assert notes_response is not None
        
        # Handle paginated response format
        if isinstance(notes_response, dict):
            if "results" in notes_response:
                assert isinstance(notes_response["results"], list)
                assert "count" in notes_response or "next" in notes_response
            else:
                assert isinstance(notes_response, dict)
        else:
            assert isinstance(notes_response, list)

    def test_update_project_note_title_only(
        self, test_project_id, cleanup_project_notes
    ):
        """Test updating only the title of a project note"""
        # Create a note first
        created_note = create_project_note(
            project_id=test_project_id,
            title="Original Title",
            content="<p>Original content.</p>",
        )

        note_id = created_note["id"]
        cleanup_project_notes.append(note_id)

        # Update only the title
        update_result = update_project_note(
            note_id=note_id,
            title="Updated Title",
        )

        # Verify the update
        updated_note = get_project_note(note_id)
        assert updated_note["title"] == "Updated Title"
        assert updated_note["content_html"] == "<p>Original content.</p>"  # Should remain unchanged
        assert updated_note["show_description"] == created_note["show_description"]  # Should remain unchanged

    def test_update_project_note_content_only(
        self, test_project_id, cleanup_project_notes
    ):
        """Test updating only the content of a project note"""
        # Create a note first
        created_note = create_project_note(
            project_id=test_project_id,
            title="Original Title",
            content="<p>Original content.</p>",
        )

        note_id = created_note["id"]
        cleanup_project_notes.append(note_id)

        # Update only the content
        update_result = update_project_note(
            note_id=note_id,
            content="<p>Updated content with <strong>formatting</strong>.</p>",
        )

        # Verify the update
        updated_note = get_project_note(note_id)
        assert updated_note["title"] == "Original Title"  # Should remain unchanged
        assert updated_note["content_html"] == "<p>Updated content with <strong>formatting</strong>.</p>"
        assert updated_note["show_description"] == created_note["show_description"]  # Should remain unchanged

    def test_update_project_note_show_description_only(
        self, test_project_id, cleanup_project_notes
    ):
        """Test updating only the show_description flag of a project note"""
        # Create a note first
        created_note = create_project_note(
            project_id=test_project_id,
            title="Original Title",
            content="<p>Original content.</p>",
            show_description=True,
        )

        note_id = created_note["id"]
        cleanup_project_notes.append(note_id)

        # Update only the show_description flag
        update_result = update_project_note(
            note_id=note_id,
            show_description=False,
        )

        # Verify the update
        updated_note = get_project_note(note_id)
        assert updated_note["title"] == "Original Title"  # Should remain unchanged
        assert updated_note["content_html"] == "<p>Original content.</p>"  # Should remain unchanged
        # Backend automatically sets show_description=True when content is not empty
        assert updated_note["show_description"] is True

    def test_update_project_note_all_fields(
        self, test_project_id, cleanup_project_notes
    ):
        """Test updating all fields of a project note"""
        # Create a note first
        created_note = create_project_note(
            project_id=test_project_id,
            title="Original Title",
            content="<p>Original content.</p>",
            show_description=True,
        )

        note_id = created_note["id"]
        cleanup_project_notes.append(note_id)

        # Update all fields
        update_result = update_project_note(
            note_id=note_id,
            title="Completely Updated Title",
            content="<p>Completely updated content with <em>emphasis</em>.</p>",
            show_description=False,
        )

        # Verify the updates
        updated_note = get_project_note(note_id)
        assert updated_note["title"] == "Completely Updated Title"
        assert updated_note["content_html"] == "<p>Completely updated content with <em>emphasis</em>.</p>"
        # Backend automatically sets show_description=True when content is not empty
        assert updated_note["show_description"] is True

    def test_delete_project_note(self, test_project_id):
        """Test deleting a project note"""
        # Create a note first
        created_note = create_project_note(
            project_id=test_project_id,
            title="Note to Delete",
            content="<p>This note will be deleted.</p>",
        )

        note_id = created_note["id"]

        # Verify it exists
        retrieved_note = get_project_note(note_id)
        assert retrieved_note["id"] == note_id

        # Delete the note
        delete_result = delete_project_note(note_id)

        # Verify it's deleted (should raise exception or return 404)
        with pytest.raises(Exception):
            get_project_note(note_id)


@pytest.mark.integration
class TestProjectNoteWorkflow:
    """End-to-end workflow tests for project notes"""

    def test_complete_project_note_workflow(self, test_project_id):
        """Test complete project note workflow: create → get → update → list → delete"""
        
        # Step 1: Create note
        created_note = create_project_note(
            project_id=test_project_id,
            title="Workflow Test Note",
            content="<p>Initial content for workflow test.</p>",
            show_description=True,
        )

        note_id = created_note["id"]
        assert created_note["title"] == "Workflow Test Note"
        assert created_note["show_description"] is True

        # Step 2: Get note details
        retrieved_note = get_project_note(note_id)
        assert retrieved_note["id"] == note_id
        assert retrieved_note["title"] == "Workflow Test Note"
        assert retrieved_note["content_html"] == "<p>Initial content for workflow test.</p>"

        # Step 3: Update note with new content and settings
        update_project_note(
            note_id=note_id,
            title="Updated Workflow Test Note",
            content="<p>Updated content with <strong>bold text</strong>.</p>",
            show_description=False,
        )

        # Step 4: Verify updates
        updated_note = get_project_note(note_id)
        assert updated_note["title"] == "Updated Workflow Test Note"
        assert updated_note["content_html"] == "<p>Updated content with <strong>bold text</strong>.</p>"
        # Backend automatically sets show_description=True when content is not empty
        assert updated_note["show_description"] is True

        # Step 5: Verify note appears in project listing
        notes_list = list_project_notes(project_id=test_project_id)
        if isinstance(notes_list, dict) and "results" in notes_list:
            notes = notes_list["results"]
        else:
            notes = notes_list

        found_note = next((note for note in notes if note["id"] == note_id), None)
        assert found_note is not None
        assert found_note["title"] == "Updated Workflow Test Note"

        # Step 6: Delete note
        delete_project_note(note_id)

        # Step 7: Verify deletion
        with pytest.raises(Exception):
            get_project_note(note_id)

    def test_multiple_notes_same_project(self, test_project_id):
        """Test creating and managing multiple notes for the same project"""
        note_ids = []

        try:
            # Create multiple notes
            for i in range(3):
                note = create_project_note(
                    project_id=test_project_id,
                    title=f"Multi-Note Test {i+1}",
                    content=f"<p>Content for note {i+1}.</p>",
                    show_description=(i % 2 == 0),  # Alternate show_description
                )
                note_ids.append(note["id"])

            # List notes and verify all are present
            notes_list = list_project_notes(project_id=test_project_id)
            if isinstance(notes_list, dict) and "results" in notes_list:
                notes = notes_list["results"]
            else:
                notes = notes_list

            created_notes = [note for note in notes if note["id"] in note_ids]
            assert len(created_notes) == 3

            # Verify each note has correct properties
            for i, note in enumerate(sorted(created_notes, key=lambda x: x["title"])):
                expected_title = f"Multi-Note Test {i+1}"
                assert note["title"] == expected_title
                assert note["content_html"] == f"<p>Content for note {i+1}.</p>"
                # Backend automatically sets show_description=True when content is not empty
                assert note["show_description"] is True

            # Update one note and verify others are unchanged
            update_project_note(
                note_id=note_ids[1],
                title="Updated Middle Note",
                content="<p>This note was updated.</p>",
            )

            # Verify the update and that others remain unchanged
            updated_note = get_project_note(note_ids[1])
            assert updated_note["title"] == "Updated Middle Note"
            assert updated_note["content_html"] == "<p>This note was updated.</p>"

            other_note = get_project_note(note_ids[0])
            assert other_note["title"] == "Multi-Note Test 1"  # Should be unchanged

        finally:
            # Cleanup all created notes
            for note_id in note_ids:
                try:
                    delete_project_note(note_id)
                except Exception:
                    pass  # Best effort cleanup


@pytest.mark.integration
class TestProjectNoteValidation:
    """Integration tests for project note validation and error handling"""

    def test_get_nonexistent_note(self):
        """Test getting a non-existent note raises appropriate error"""
        with pytest.raises(Exception):
            get_project_note("non-existent-note-id")

    def test_update_nonexistent_note(self):
        """Test updating a non-existent note raises appropriate error"""
        with pytest.raises(Exception):
            update_project_note(
                note_id="non-existent-note-id",
                title="Updated Title",
            )

    def test_delete_nonexistent_note(self):
        """Test deleting a non-existent note raises appropriate error"""
        with pytest.raises(Exception):
            delete_project_note("non-existent-note-id")

    def test_create_note_invalid_project(self):
        """Test creating a note with invalid project ID raises appropriate error"""
        with pytest.raises(Exception):
            create_project_note(
                project_id="invalid-project-id",
                title="Test Note",
                content="<p>Test content.</p>",
            )

    def test_list_notes_invalid_project(self):
        """Test listing notes with invalid project ID - should raise validation error"""
        # Backend validates UUID format and raises APIValidationError for invalid UUID
        from atomict.exceptions import APIValidationError
        
        with pytest.raises(APIValidationError) as exc_info:
            list_project_notes(project_id="invalid-project-id")
        
        # Should contain UUID validation error
        error_data = exc_info.value.args[0]
        assert "project__id" in error_data
        assert "valid UUID" in error_data["project__id"][0]
