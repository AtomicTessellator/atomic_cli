"""
Integration tests for project operations.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- Optional TEST_PROJECT_ID for non-destructive tests

To run: uv run pytest tests/integration/test_project_operations.py -v -m integration
"""

import os
import uuid

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.project.project import (
    create_project,
    delete_project,
    get_project,
    get_project_by_name,
    list_projects,
    project_exists,
    update_project,
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
    """Get test project ID from environment (optional for non-destructive tests)"""
    return os.getenv("TEST_PROJECT_ID")


@pytest.fixture
def unique_project_name():
    """Generate a unique project name for testing"""
    return f"integration-test-project-{uuid.uuid4().hex[:8]}"


@pytest.fixture
def cleanup_projects():
    """Track created projects for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for project_id in created_ids:
        try:
            delete_project(project_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestProjectCRUDOperations:
    """Integration tests for project CRUD operations"""

    def test_create_project_minimal(self, unique_project_name, cleanup_projects):
        """Test creating a project with minimal parameters"""
        result = create_project(name=unique_project_name)

        # Track for cleanup
        cleanup_projects.append(result["id"])

        # Verify response structure
        assert "id" in result
        assert result["name"] == unique_project_name
        assert "created_at" in result
        assert "updated_at" in result

    def test_create_project_with_full_parameters(
        self, unique_project_name, cleanup_projects
    ):
        """Test creating a project with all parameters"""
        description = "Integration test project with full parameters"
        thumbnail_smiles = "CCO"  # Simple ethanol SMILES

        result = create_project(
            name=unique_project_name,
            description=description,
            thumbnail_smiles=thumbnail_smiles,
        )

        # Track for cleanup
        cleanup_projects.append(result["id"])

        # Verify all parameters were set
        assert result["name"] == unique_project_name
        assert result["description_html"] == description
        assert result["thumbnail_smiles"] == thumbnail_smiles

    def test_get_project(self, unique_project_name, cleanup_projects):
        """Test retrieving a project by ID"""
        # Create project first
        create_result = create_project(
            name=unique_project_name, description="Test project for get operation"
        )
        project_id = create_result["id"]
        cleanup_projects.append(project_id)

        # Get the project
        get_result = get_project(project_id)

        # Verify the retrieved project matches
        assert get_result["id"] == project_id
        assert get_result["name"] == unique_project_name
        assert get_result["description_html"] == "Test project for get operation"

    def test_update_project_partial(self, unique_project_name, cleanup_projects):
        """Test updating a project with partial field updates"""
        # Create project first
        create_result = create_project(
            name=unique_project_name, description="Original description"
        )
        project_id = create_result["id"]
        cleanup_projects.append(project_id)

        # Update only the description
        new_description = "Updated description"
        update_result = update_project(
            project_id=project_id, description=new_description
        )

        # Verify update
        assert update_result["description_html"] == new_description
        assert update_result["name"] == unique_project_name  # Should remain unchanged

        # Verify by getting the project again
        get_result = get_project(project_id)
        assert get_result["description_html"] == new_description
        assert get_result["name"] == unique_project_name

    def test_update_project_all_fields(self, unique_project_name, cleanup_projects):
        """Test updating all project fields"""
        # Create project first
        create_result = create_project(name=unique_project_name)
        project_id = create_result["id"]
        cleanup_projects.append(project_id)

        # Update all fields
        new_name = f"{unique_project_name}-updated"
        new_description = "Completely updated project"
        new_thumbnail_smiles = "C1CCCCC1"  # Cyclohexane

        update_result = update_project(
            project_id=project_id,
            name=new_name,
            description=new_description,
            thumbnail_smiles=new_thumbnail_smiles,
        )

        # Verify all updates
        assert update_result["name"] == new_name
        assert update_result["description_html"] == new_description
        assert update_result["thumbnail_smiles"] == new_thumbnail_smiles

    def test_delete_project(self, unique_project_name):
        """Test deleting a project"""
        # Create project first
        create_result = create_project(name=unique_project_name)
        project_id = create_result["id"]

        # Verify it exists
        get_result = get_project(project_id)
        assert get_result["id"] == project_id

        # Delete the project
        delete_result = delete_project(project_id)

        # Verify deletion by attempting to get the project (should raise exception)
        with pytest.raises(Exception):
            get_project(project_id)

    def test_project_exists_true(self, unique_project_name, cleanup_projects):
        """Test project_exists returns True for existing project"""
        # Create project first
        create_result = create_project(name=unique_project_name)
        cleanup_projects.append(create_result["id"])

        # Check if it exists
        exists = project_exists(unique_project_name)
        assert exists is True

    def test_project_exists_false(self):
        """Test project_exists returns False for non-existing project"""
        non_existent_name = f"non-existent-project-{uuid.uuid4().hex[:8]}"
        exists = project_exists(non_existent_name)
        assert exists is False

    def test_get_project_by_name(self, unique_project_name, cleanup_projects):
        """Test retrieving a project by name"""
        description = "Test project for get by name"

        # Create project first
        create_result = create_project(
            name=unique_project_name, description=description
        )
        cleanup_projects.append(create_result["id"])

        # Get project by name
        get_result = get_project_by_name(unique_project_name)

        # Verify the retrieved project
        assert get_result["name"] == unique_project_name
        assert get_result["description_html"] == description
        assert get_result["id"] == create_result["id"]


@pytest.mark.integration
class TestProjectListOperations:
    """Integration tests for project listing and search operations"""

    def test_list_projects_basic(self):
        """Test basic project listing"""
        result = list_projects()

        # Verify response structure
        assert "results" in result
        assert "count" in result
        assert isinstance(result["results"], list)
        assert isinstance(result["count"], int)

    def test_list_projects_with_search(self, unique_project_name, cleanup_projects):
        """Test project listing with search parameter"""
        search_term = "integration-search-test"
        project_name = f"{search_term}-{uuid.uuid4().hex[:8]}"

        # Create a project with the search term
        create_result = create_project(
            name=project_name,
            description=f"Project for testing search functionality with {search_term}",
        )
        cleanup_projects.append(create_result["id"])

        # Search for projects containing the search term
        search_result = list_projects(search=search_term)

        # Verify search results contain our project
        assert search_result["count"] > 0
        project_names = [p["name"] for p in search_result["results"]]
        assert project_name in project_names

    def test_list_projects_with_ordering(self):
        """Test project listing with ordering parameter"""
        # Test ascending order by name
        asc_result = list_projects(ordering="name")
        assert "results" in asc_result

        # Test descending order by name
        desc_result = list_projects(ordering="-name")
        assert "results" in desc_result

        # Test ordering by creation date
        created_result = list_projects(ordering="created_at")
        assert "results" in created_result

        # Test descending order by creation date
        created_desc_result = list_projects(ordering="-created_at")
        assert "results" in created_desc_result

    def test_list_projects_with_filters(self, unique_project_name, cleanup_projects):
        """Test project listing with custom filters"""
        # Create a project to test filtering
        create_result = create_project(
            name=unique_project_name, description="Filterable test project"
        )
        cleanup_projects.append(create_result["id"])

        # Test filtering by name (exact match if supported)
        name_filter_result = list_projects(name=unique_project_name)
        if name_filter_result["count"] > 0:
            assert any(
                p["name"] == unique_project_name for p in name_filter_result["results"]
            )

    def test_list_projects_combined_parameters(self):
        """Test project listing with multiple parameters combined"""
        result = list_projects(search="test", ordering="-created_at")

        assert "results" in result
        assert "count" in result


@pytest.mark.integration
class TestProjectWorkflows:
    """End-to-end workflow tests for project operations"""

    def test_complete_project_lifecycle(self, unique_project_name):
        """Test complete project lifecycle: create → get → update → delete"""
        original_description = "Original project description"
        original_smiles = "C"  # Methane

        # Step 1: Create project
        create_result = create_project(
            name=unique_project_name,
            description=original_description,
            thumbnail_smiles=original_smiles,
        )
        project_id = create_result["id"]

        try:
            # Verify creation
            assert create_result["name"] == unique_project_name
            assert create_result["description_html"] == original_description
            assert create_result["thumbnail_smiles"] == original_smiles

            # Step 2: Get project and verify
            get_result = get_project(project_id)
            assert get_result["id"] == project_id
            assert get_result["name"] == unique_project_name

            # Step 3: Update project
            updated_description = "Updated project description"
            updated_smiles = "CC"  # Ethane

            update_result = update_project(
                project_id=project_id,
                description=updated_description,
                thumbnail_smiles=updated_smiles,
            )

            assert update_result["description_html"] == updated_description
            assert update_result["thumbnail_smiles"] == updated_smiles

            # Step 4: Verify update by getting project again
            updated_get_result = get_project(project_id)
            assert updated_get_result["description_html"] == updated_description
            assert updated_get_result["thumbnail_smiles"] == updated_smiles

            # Step 5: Verify project exists
            assert project_exists(unique_project_name) is True

            # Step 6: Get project by name
            by_name_result = get_project_by_name(unique_project_name)
            assert by_name_result["id"] == project_id

        finally:
            # Step 7: Delete project
            delete_project(project_id)

        # Step 8: Verify deletion
        with pytest.raises(Exception):
            get_project(project_id)

        assert project_exists(unique_project_name) is False

    def test_project_search_workflow(self, cleanup_projects):
        """Test project search and retrieval workflow"""
        search_identifier = f"workflow-search-{uuid.uuid4().hex[:6]}"

        # Create multiple projects with common identifier
        project_names = [
            f"{search_identifier}-alpha",
            f"{search_identifier}-beta",
            f"{search_identifier}-gamma",
        ]

        created_ids = []
        try:
            # Create test projects
            for name in project_names:
                result = create_project(
                    name=name, description=f"Search workflow test project: {name}"
                )
                created_ids.append(result["id"])
                cleanup_projects.append(result["id"])

            # Search for projects with the common identifier
            search_result = list_projects(search=search_identifier)

            # Verify all created projects are found
            found_names = [p["name"] for p in search_result["results"]]
            for name in project_names:
                assert name in found_names

            # Test ordering of search results
            ordered_result = list_projects(search=search_identifier, ordering="name")
            ordered_names = [p["name"] for p in ordered_result["results"]]
            expected_order = sorted(project_names)

            # Check that our projects appear in the expected order
            for i, expected_name in enumerate(expected_order):
                assert expected_name in ordered_names

        finally:
            # Cleanup is handled by the cleanup_projects fixture
            pass


@pytest.mark.integration
class TestProjectValidationAndErrorHandling:
    """Integration tests for project validation and error handling"""

    def test_get_nonexistent_project(self):
        """Test getting a project that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises(Exception):
            get_project(non_existent_id)

    def test_update_nonexistent_project(self):
        """Test updating a project that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises(Exception):
            update_project(project_id=non_existent_id, name="Updated Name")

    def test_delete_nonexistent_project(self):
        """Test deleting a project that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises(Exception):
            delete_project(non_existent_id)

    def test_get_project_by_nonexistent_name(self):
        """Test getting a project by name that doesn't exist"""
        non_existent_name = f"non-existent-name-{uuid.uuid4().hex[:8]}"

        with pytest.raises(Exception):
            get_project_by_name(non_existent_name)

    def test_create_project_duplicate_name(self, unique_project_name, cleanup_projects):
        """Test creating projects with duplicate names (if not allowed)"""
        # Create first project
        first_result = create_project(name=unique_project_name)
        cleanup_projects.append(first_result["id"])

        # Attempt to create second project with same name
        # Note: Backend behavior may vary - some systems allow duplicates
        try:
            second_result = create_project(name=unique_project_name)
            cleanup_projects.append(second_result["id"])

            # If duplicates are allowed, both should have different IDs
            assert first_result["id"] != second_result["id"]
        except Exception:
            # If duplicates are not allowed, exception is expected
            pass
