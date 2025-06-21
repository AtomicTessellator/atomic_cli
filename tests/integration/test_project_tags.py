"""
Integration tests for project tags functionality.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)

To run: PYTHONPATH=/path/to/atomic_cli uv run pytest tests/integration/test_project_tags.py -v -m integration
"""

import os
import uuid

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.project.tags import (
    VALID_TAG_COLOURS,
    create_project_tag,
    create_tag,
    delete_project_tag,
    delete_project_tag_association,
    get_project_tag,
    get_tag_by_name,
    list_project_tag_associations,
    list_project_tags,
    project_tag_exists,
    tag_exists,
    update_project_tag,
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
def cleanup_tags():
    """Track created tags for cleanup"""
    created_tag_ids = []
    yield created_tag_ids

    # Cleanup after test
    for tag_id in created_tag_ids:
        try:
            delete_project_tag(tag_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def cleanup_associations():
    """Track created tag-project associations for cleanup"""
    created_association_ids = []
    yield created_association_ids

    # Cleanup after test
    for association_id in created_association_ids:
        try:
            delete_project_tag_association(association_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def test_tag_name():
    """Generate unique tag name for testing"""
    return f"test-tag-{uuid.uuid4().hex[:8]}"


@pytest.fixture
def test_tag_color():
    """Get a valid tag color for testing"""
    return "bg-primary"


@pytest.mark.integration
class TestProjectTagCRUD:
    """Test CRUD operations for project tags"""

    def test_create_tag_success(self, test_tag_name, test_tag_color, cleanup_tags):
        """Test successful tag creation"""
        response = create_tag(test_tag_name, test_tag_color)

        assert response is not None
        tag_id = response.get("id")
        assert tag_id is not None
        cleanup_tags.append(tag_id)

        assert response["tag"] == test_tag_name
        assert response["color"] == test_tag_color

    def test_create_tag_invalid_color(self, test_tag_name):
        """Test tag creation with invalid color"""
        with pytest.raises(ValueError, match="Invalid tag color"):
            create_tag(test_tag_name, "invalid-color")

    def test_create_tag_valid_colors(self, cleanup_tags):
        """Test tag creation with all valid colors"""
        for i, color in enumerate(VALID_TAG_COLOURS[:3]):  # Test first 3 colors
            tag_name = f"color-test-{i}-{uuid.uuid4().hex[:8]}"
            response = create_tag(tag_name, color)

            assert response["color"] == color
            cleanup_tags.append(response["id"])

    def test_get_tag_by_name(self, test_tag_name, test_tag_color, cleanup_tags):
        """Test retrieving tag by name"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        cleanup_tags.append(created_tag["id"])

        # Get tag by name
        retrieved_tag = get_tag_by_name(test_tag_name)

        assert retrieved_tag["id"] == created_tag["id"]
        assert retrieved_tag["tag"] == test_tag_name
        assert retrieved_tag["color"] == test_tag_color

    def test_tag_exists_true(self, test_tag_name, test_tag_color, cleanup_tags):
        """Test tag_exists returns True for existing tag"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        cleanup_tags.append(created_tag["id"])

        # Check existence
        assert tag_exists(test_tag_name) is True

    def test_tag_exists_false(self):
        """Test tag_exists returns False for non-existing tag"""
        non_existent_name = f"non-existent-{uuid.uuid4().hex}"
        assert tag_exists(non_existent_name) is False

    def test_get_project_tag(self, test_tag_name, test_tag_color, cleanup_tags):
        """Test retrieving tag by ID"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Get tag by ID
        retrieved_tag = get_project_tag(tag_id)

        assert retrieved_tag["id"] == tag_id
        assert retrieved_tag["tag"] == test_tag_name
        assert retrieved_tag["color"] == test_tag_color

    def test_update_project_tag_name(self, test_tag_name, test_tag_color, cleanup_tags):
        """Test updating tag name"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Update tag name
        new_name = f"updated-{test_tag_name}"
        updated_tag = update_project_tag(tag_id, tag=new_name)

        assert updated_tag["tag"] == new_name
        assert updated_tag["color"] == test_tag_color

    @pytest.mark.xfail(
        reason="Backend bug: ProjectTagViewSet.perform_update() fails when updating only color field (KeyError on serializer.validated_data['tag'])"
    )
    def test_update_project_tag_color(
        self, test_tag_name, test_tag_color, cleanup_tags
    ):
        """Test updating tag color"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Update tag color
        new_color = "bg-danger"
        updated_tag = update_project_tag(tag_id, color=new_color)

        assert updated_tag["tag"] == test_tag_name
        assert updated_tag["color"] == new_color

    def test_update_project_tag_both(self, test_tag_name, test_tag_color, cleanup_tags):
        """Test updating both tag name and color"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Update both name and color
        new_name = f"updated-{test_tag_name}"
        new_color = "bg-success"
        updated_tag = update_project_tag(tag_id, tag=new_name, color=new_color)

        assert updated_tag["tag"] == new_name
        assert updated_tag["color"] == new_color

    def test_update_project_tag_invalid_color(
        self, test_tag_name, test_tag_color, cleanup_tags
    ):
        """Test updating tag with invalid color"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Try to update with invalid color
        with pytest.raises(ValueError, match="Invalid tag color"):
            update_project_tag(tag_id, color="invalid-color")

    def test_update_project_tag_no_params(
        self, test_tag_name, test_tag_color, cleanup_tags
    ):
        """Test updating tag with no parameters"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Try to update with no parameters
        with pytest.raises(ValueError, match="At least one parameter"):
            update_project_tag(tag_id)

    def test_list_project_tags(self, test_tag_name, test_tag_color, cleanup_tags):
        """Test listing all project tags"""
        # Create a test tag
        created_tag = create_tag(test_tag_name, test_tag_color)
        cleanup_tags.append(created_tag["id"])

        # List all tags
        response = list_project_tags()

        assert "results" in response
        assert "count" in response
        assert response["count"] > 0

        # Check that our created tag is in the list
        tag_ids = [tag["id"] for tag in response["results"]]
        assert created_tag["id"] in tag_ids

    def test_delete_project_tag(self, test_tag_name, test_tag_color):
        """Test deleting a project tag"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]

        # Delete the tag
        delete_response = delete_project_tag(tag_id)
        # Note: delete might return empty response or 204 status

        # Verify tag no longer exists
        assert tag_exists(test_tag_name) is False


@pytest.mark.integration
class TestProjectTagAssociations:
    """Test project-tag association operations"""

    def test_create_project_tag_association(
        self,
        test_project_id,
        test_tag_name,
        test_tag_color,
        cleanup_tags,
        cleanup_associations,
    ):
        """Test creating project-tag association"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Create association
        association = create_project_tag(test_project_id, tag_id)

        assert association is not None
        association_id = association.get("id")
        assert association_id is not None
        cleanup_associations.append(association_id)

        assert association["project"] == test_project_id
        assert association["project_tag"] == tag_id

    def test_project_tag_exists_true(
        self,
        test_project_id,
        test_tag_name,
        test_tag_color,
        cleanup_tags,
        cleanup_associations,
    ):
        """Test project_tag_exists returns True for existing association"""
        # Create tag first
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Create association
        association = create_project_tag(test_project_id, tag_id)
        cleanup_associations.append(association["id"])

        # Check existence
        assert project_tag_exists(test_project_id, tag_id) is True

    def test_project_tag_exists_false(
        self, test_project_id, test_tag_name, test_tag_color, cleanup_tags
    ):
        """Test project_tag_exists returns False for non-existing association"""
        # Create tag but don't associate it
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        # Check non-existent association
        assert project_tag_exists(test_project_id, tag_id) is False

    def test_list_project_tag_associations_by_project(
        self,
        test_project_id,
        test_tag_name,
        test_tag_color,
        cleanup_tags,
        cleanup_associations,
    ):
        """Test listing associations filtered by project"""
        # Create tag and association
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        association = create_project_tag(test_project_id, tag_id)
        cleanup_associations.append(association["id"])

        # List associations for this project
        response = list_project_tag_associations(project_id=test_project_id)

        assert "results" in response
        assert "count" in response
        assert response["count"] > 0

        # Check that our association is in the list
        association_ids = [assoc["id"] for assoc in response["results"]]
        assert association["id"] in association_ids

    def test_list_project_tag_associations_by_tag(
        self,
        test_project_id,
        test_tag_name,
        test_tag_color,
        cleanup_tags,
        cleanup_associations,
    ):
        """Test listing associations filtered by tag"""
        # Create tag and association
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        association = create_project_tag(test_project_id, tag_id)
        cleanup_associations.append(association["id"])

        # List associations for this tag
        response = list_project_tag_associations(tag_id=tag_id)

        assert "results" in response
        assert "count" in response
        assert response["count"] > 0

        # Check that our association is in the list
        association_ids = [assoc["id"] for assoc in response["results"]]
        assert association["id"] in association_ids

    def test_list_project_tag_associations_both_filters(
        self,
        test_project_id,
        test_tag_name,
        test_tag_color,
        cleanup_tags,
        cleanup_associations,
    ):
        """Test listing associations filtered by both project and tag"""
        # Create tag and association
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        association = create_project_tag(test_project_id, tag_id)
        cleanup_associations.append(association["id"])

        # List associations with both filters
        response = list_project_tag_associations(
            project_id=test_project_id, tag_id=tag_id
        )

        assert "results" in response
        assert "count" in response
        assert response["count"] > 0

        # Should find exactly our association
        found_association = None
        for assoc in response["results"]:
            if assoc["id"] == association["id"]:
                found_association = assoc
                break

        assert found_association is not None
        assert found_association["project"] == test_project_id
        assert found_association["project_tag"] == tag_id

    def test_list_project_tag_associations_no_filters(self):
        """Test listing all associations without filters"""
        response = list_project_tag_associations()

        assert "results" in response
        assert "count" in response
        # Should return all associations (count could be 0 if no associations exist)

    def test_delete_project_tag_association(
        self, test_project_id, test_tag_name, test_tag_color, cleanup_tags
    ):
        """Test deleting a project-tag association"""
        # Create tag and association
        created_tag = create_tag(test_tag_name, test_tag_color)
        tag_id = created_tag["id"]
        cleanup_tags.append(tag_id)

        association = create_project_tag(test_project_id, tag_id)
        association_id = association["id"]

        # Delete the association
        delete_response = delete_project_tag_association(association_id)
        # Note: delete might return empty response or 204 status

        # Verify association no longer exists
        assert project_tag_exists(test_project_id, tag_id) is False


@pytest.mark.integration
class TestProjectTagWorkflows:
    """Test complete workflows combining multiple operations"""

    @pytest.mark.xfail(
        reason="Backend bug: ProjectTagViewSet.perform_update() fails when updating only color field (KeyError on serializer.validated_data['tag'])"
    )
    def test_full_tag_workflow(self, test_project_id):
        """Test complete workflow: create tag → associate → update → delete"""
        tag_name = f"workflow-tag-{uuid.uuid4().hex[:8]}"
        initial_color = "bg-primary"

        # Step 1: Create tag
        created_tag = create_tag(tag_name, initial_color)
        tag_id = created_tag["id"]

        try:
            # Step 2: Verify tag exists
            assert tag_exists(tag_name) is True

            # Step 3: Create project association
            association = create_project_tag(test_project_id, tag_id)
            association_id = association["id"]

            try:
                # Step 4: Verify association exists
                assert project_tag_exists(test_project_id, tag_id) is True

                # Step 5: Update tag
                new_color = "bg-success"
                updated_tag = update_project_tag(tag_id, color=new_color)
                assert updated_tag["color"] == new_color

                # Step 6: Verify tag is still associated
                assert project_tag_exists(test_project_id, tag_id) is True

                # Step 7: Delete association
                delete_project_tag_association(association_id)
                assert project_tag_exists(test_project_id, tag_id) is False

            except Exception:
                # Cleanup association if it exists
                try:
                    delete_project_tag_association(association_id)
                except Exception:
                    pass
                raise

            # Step 8: Delete tag
            delete_project_tag(tag_id)
            assert tag_exists(tag_name) is False

        except Exception:
            # Cleanup tag if it exists
            try:
                delete_project_tag(tag_id)
            except Exception:
                pass
            raise

    def test_multiple_associations_same_tag(self, test_project_id):
        """Test associating the same tag with multiple projects (if possible)"""
        tag_name = f"multi-assoc-{uuid.uuid4().hex[:8]}"

        # Create tag
        created_tag = create_tag(tag_name, "bg-info")
        tag_id = created_tag["id"]

        try:
            # Create association with test project
            association1 = create_project_tag(test_project_id, tag_id)

            try:
                # Verify association exists
                assert project_tag_exists(test_project_id, tag_id) is True

                # Test listing associations for this tag
                tag_associations = list_project_tag_associations(tag_id=tag_id)
                assert tag_associations["count"] >= 1

                # Find our association in the results
                found = False
                for assoc in tag_associations["results"]:
                    if (
                        assoc["project"] == test_project_id
                        and assoc["project_tag"] == tag_id
                    ):
                        found = True
                        break
                assert found is True

                # Cleanup association
                delete_project_tag_association(association1["id"])

            except Exception:
                # Cleanup association
                try:
                    delete_project_tag_association(association1["id"])
                except Exception:
                    pass
                raise

            # Cleanup tag
            delete_project_tag(tag_id)

        except Exception:
            # Cleanup tag
            try:
                delete_project_tag(tag_id)
            except Exception:
                pass
            raise
