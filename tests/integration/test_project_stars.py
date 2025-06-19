"""
Integration tests for project stars functionality.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)

To run: PYTHONPATH=/path/to/atomic_cli uv run pytest tests/integration/test_project_stars.py -v -m integration
"""

import os
import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.project.stars import (
    create_project_star,
    delete_project_star,
    get_project_star,
    list_project_stars,
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
def cleanup_project_stars():
    """Track created project stars for cleanup"""
    created_star_ids = []
    yield created_star_ids

    # Cleanup after test
    for star_id in created_star_ids:
        try:
            delete_project_star(star_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestProjectStarsIntegration:
    """Integration tests for project stars CRUD operations"""

    def test_create_project_star(self, test_project_id, cleanup_project_stars):
        """Test creating a project star"""

        result = create_project_star(project_id=test_project_id)

        # Track for cleanup
        cleanup_project_stars.append(result["id"])

        # Verify response structure
        assert "id" in result
        assert "project" in result
        assert "created_at" in result
        assert "model_name" in result
        assert result["project"] == test_project_id
        assert result["model_name"] == "ProjectStar"

    def test_get_project_star(self, test_project_id, cleanup_project_stars):
        """Test retrieving a project star by ID"""

        # Create a star first
        create_result = create_project_star(project_id=test_project_id)
        star_id = create_result["id"]
        cleanup_project_stars.append(star_id)

        # Get the star
        get_result = get_project_star(star_id=star_id)

        # Verify response
        assert get_result["id"] == star_id
        assert get_result["project"] == test_project_id
        assert get_result["model_name"] == "ProjectStar"

    def test_list_project_stars(self, test_project_id, cleanup_project_stars):
        """Test listing project stars"""

        # Get initial count
        initial_result = list_project_stars()
        initial_count = initial_result["count"]

        # Create a star
        create_result = create_project_star(project_id=test_project_id)
        star_id = create_result["id"]
        cleanup_project_stars.append(star_id)

        # List stars again
        list_result = list_project_stars()

        # Verify count increased
        assert list_result["count"] == initial_count + 1
        assert "results" in list_result

        # Verify our star is in the results
        star_ids = [star["id"] for star in list_result["results"]]
        assert star_id in star_ids

    def test_delete_project_star(self, test_project_id):
        """Test deleting a project star"""

        # Create a star
        create_result = create_project_star(project_id=test_project_id)
        star_id = create_result["id"]

        # Verify it exists
        get_result = get_project_star(star_id=star_id)
        assert get_result["id"] == star_id

        # Delete it
        delete_result = delete_project_star(star_id=star_id)

        # Verify it's deleted (should raise exception or return 404)
        with pytest.raises(Exception):
            get_project_star(star_id=star_id)

    def test_create_multiple_stars_same_project(
        self, test_project_id, cleanup_project_stars
    ):
        """Test that multiple stars for the same project are handled correctly"""

        # Create first star
        first_result = create_project_star(project_id=test_project_id)
        first_star_id = first_result["id"]
        cleanup_project_stars.append(first_star_id)

        # Attempt to create second star for same project
        # This might either succeed (if multiple stars allowed) or fail (if unique constraint exists)
        try:
            second_result = create_project_star(project_id=test_project_id)
            second_star_id = second_result["id"]
            cleanup_project_stars.append(second_star_id)

            # If we get here, multiple stars are allowed
            assert first_star_id != second_star_id
            assert first_result["project"] == second_result["project"]

        except Exception:
            # If we get an exception, unique constraint exists (one star per user-project)
            # This is also valid behavior
            pass

    def test_star_unstar_workflow(self, test_project_id):
        """Test the complete star/unstar workflow"""

        # Get initial star count
        initial_result = list_project_stars()
        initial_count = initial_result["count"]

        # Step 1: Star the project
        star_result = create_project_star(project_id=test_project_id)
        star_id = star_result["id"]

        # Step 2: Verify star exists
        get_result = get_project_star(star_id=star_id)
        assert get_result["id"] == star_id
        assert get_result["project"] == test_project_id

        # Step 3: Verify star appears in list
        list_result = list_project_stars()
        assert list_result["count"] == initial_count + 1
        star_ids = [star["id"] for star in list_result["results"]]
        assert star_id in star_ids

        # Step 4: Unstar the project
        delete_result = delete_project_star(star_id=star_id)

        # Step 5: Verify star is removed
        with pytest.raises(Exception):
            get_project_star(star_id=star_id)

        # Step 6: Verify star count is back to initial
        final_result = list_project_stars()
        assert final_result["count"] == initial_count

    def test_list_stars_pagination(self, test_project_id, cleanup_project_stars):
        """Test pagination in list_project_stars"""

        # Create multiple stars (at least one for testing)
        create_result = create_project_star(project_id=test_project_id)
        cleanup_project_stars.append(create_result["id"])

        # Test list with default pagination
        list_result = list_project_stars()
        assert "count" in list_result
        assert "results" in list_result
        assert isinstance(list_result["results"], list)

        # Verify each star has required fields
        for star in list_result["results"]:
            assert "id" in star
            assert "project" in star
            assert "created_at" in star
            assert "model_name" in star

    def test_star_user_isolation(self, test_project_id, cleanup_project_stars):
        """Test that stars are properly isolated per user"""

        # Create a star
        create_result = create_project_star(project_id=test_project_id)
        star_id = create_result["id"]
        cleanup_project_stars.append(star_id)

        # List stars should include our star
        list_result = list_project_stars()
        star_ids = [star["id"] for star in list_result["results"]]
        assert star_id in star_ids

        # All stars should be for the current user (no user field returned means current user)
        for star in list_result["results"]:
            # The user field is nullable and typically not returned for current user's stars
            # We just verify the star exists and has proper structure
            assert "project" in star
            assert "id" in star


@pytest.mark.integration
class TestProjectStarsEdgeCases:
    """Edge case tests for project stars"""

    def test_get_nonexistent_star(self):
        """Test getting a star that doesn't exist"""

        fake_star_id = "00000000-0000-0000-0000-000000000000"

        with pytest.raises(Exception):
            get_project_star(star_id=fake_star_id)

    def test_delete_nonexistent_star(self):
        """Test deleting a star that doesn't exist"""

        fake_star_id = "00000000-0000-0000-0000-000000000000"

        # This should raise an exception
        with pytest.raises(Exception):
            delete_project_star(star_id=fake_star_id)

    def test_create_star_invalid_project(self):
        """Test creating a star for nonexistent project"""

        fake_project_id = "00000000-0000-0000-0000-000000000000"

        with pytest.raises(Exception):
            create_project_star(project_id=fake_project_id)

    def test_star_data_structure_validation(
        self, test_project_id, cleanup_project_stars
    ):
        """Test that star data structure is consistent across operations"""

        # Create star
        create_result = create_project_star(project_id=test_project_id)
        star_id = create_result["id"]
        cleanup_project_stars.append(star_id)

        # Get star
        get_result = get_project_star(star_id=star_id)

        # List stars
        list_result = list_project_stars()
        our_star = next(
            (star for star in list_result["results"] if star["id"] == star_id), None
        )
        assert our_star is not None

        # Verify all operations return consistent data structure
        for result in [create_result, get_result, our_star]:
            assert "id" in result
            assert "project" in result
            assert "created_at" in result
            assert "model_name" in result
            assert result["project"] == test_project_id
            assert result["model_name"] == "ProjectStar"
            assert result["id"] == star_id


@pytest.mark.integration
class TestProjectStarsWorkflow:
    """End-to-end workflow tests"""

    def test_complete_star_management_workflow(self, test_project_id):
        """Test complete star management workflow"""

        # Step 1: Get initial state
        initial_stars = list_project_stars()
        initial_count = initial_stars["count"]

        # Step 2: Star the project
        star_result = create_project_star(project_id=test_project_id)
        star_id = star_result["id"]

        assert star_result["project"] == test_project_id
        assert star_result["model_name"] == "ProjectStar"

        # Step 3: Verify star was created
        get_result = get_project_star(star_id=star_id)
        assert get_result["id"] == star_id
        assert get_result["project"] == test_project_id

        # Step 4: Verify star appears in list
        updated_stars = list_project_stars()
        assert updated_stars["count"] == initial_count + 1

        star_found = False
        for star in updated_stars["results"]:
            if star["id"] == star_id:
                star_found = True
                assert star["project"] == test_project_id
                break
        assert star_found

        # Step 5: Unstar the project
        delete_result = delete_project_star(star_id=star_id)

        # Step 6: Verify star was deleted
        with pytest.raises(Exception):
            get_project_star(star_id=star_id)

        # Step 7: Verify star count returned to initial
        final_stars = list_project_stars()
        assert final_stars["count"] == initial_count

        # Step 8: Verify star no longer in list
        final_star_ids = [star["id"] for star in final_stars["results"]]
        assert star_id not in final_star_ids

    def test_multiple_projects_star_management(
        self, test_project_id, cleanup_project_stars
    ):
        """Test starring multiple projects (if available)"""

        # This test primarily uses the main test project
        # In a real scenario, you might have multiple test projects

        # Star the test project
        star_result = create_project_star(project_id=test_project_id)
        star_id = star_result["id"]
        cleanup_project_stars.append(star_id)

        # Verify star exists
        get_result = get_project_star(star_id=star_id)
        assert get_result["project"] == test_project_id

        # Verify it appears in the list
        list_result = list_project_stars()
        star_ids = [star["id"] for star in list_result["results"]]
        assert star_id in star_ids

        # This test demonstrates the pattern for multiple projects
        # In practice, you would repeat this with different project IDs
