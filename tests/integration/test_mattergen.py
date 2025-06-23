import os
import uuid
import pytest
from dotenv import load_dotenv
from atomict.simulation.mattergen import (
    create_mattergen,
    get_mattergen,
    associate_user_upload_with_mattergen,
    DEFAULT_BATCH_SIZE,
    DEFAULT_NUM_BATCHES,
)
from atomict.exceptions import APIValidationError, PermissionDenied
from atomict.auth import authenticate


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


@pytest.mark.integration
class TestMattergenIntegration:

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test environment"""
        # Check required environment variables
        self.project_id = os.environ.get("TEST_PROJECT_ID")
        if not self.project_id:
            pytest.skip("TEST_PROJECT_ID not set")

    def test_create_mattergen_draft(self):
        """Test creating a MatterGen exploration in DRAFT mode"""
        exploration_name = f"Test MatterGen Draft {uuid.uuid4().hex[:8]}"

        result = create_mattergen(
            project_id=self.project_id,
            name=exploration_name,
            description="Integration test exploration",
            batch_size=8,
            num_batches=2,
        )

        # Verify response structure
        assert "id" in result
        assert result["name"] == exploration_name
        assert result["description"] == "Integration test exploration"
        assert result["batch_size"] == 8
        assert result["num_batches"] == 2
        assert "task" in result
        assert result["task"]["status"] == 0  # DRAFT status

        # Store ID for cleanup and further tests
        self.exploration_id = result["id"]

    def test_create_mattergen_with_diffusion_guidance(self):
        """Test creating exploration with diffusion guidance factor"""
        exploration_name = f"Test MatterGen Guidance {uuid.uuid4().hex[:8]}"

        result = create_mattergen(
            project_id=self.project_id,
            name=exploration_name,
            diffusion_guidance_factor=25,
        )

        assert "id" in result
        assert result["diffusion_guidance_factor"] == 25

    def test_create_mattergen_defaults(self):
        """Test creating exploration with default parameters"""
        exploration_name = f"Test MatterGen Defaults {uuid.uuid4().hex[:8]}"

        result = create_mattergen(project_id=self.project_id, name=exploration_name)

        assert result["batch_size"] == DEFAULT_BATCH_SIZE
        assert result["num_batches"] == DEFAULT_NUM_BATCHES
        assert result["description"] == ""
        assert result["task"]["status"] == 0  # DRAFT

    def test_get_mattergen(self):
        """Test retrieving a MatterGen exploration"""
        # First create an exploration
        exploration_name = f"Test MatterGen Get {uuid.uuid4().hex[:8]}"
        create_result = create_mattergen(
            project_id=self.project_id,
            name=exploration_name,
            description="Test get functionality",
        )

        exploration_id = create_result["id"]

        # Then retrieve it
        get_result = get_mattergen(exploration_id)

        assert get_result["id"] == exploration_id
        assert get_result["name"] == exploration_name
        assert get_result["description"] == "Test get functionality"
        assert "task" in get_result

    def test_create_mattergen_invalid_project(self):
        """Test error handling for invalid project ID"""
        with pytest.raises((APIValidationError, PermissionDenied)):
            create_mattergen(
                project_id="non-existent-project-id", name="Test Exploration"
            )

    def test_field_mappings_verification(self):
        """Test that field mappings work correctly with backend"""
        exploration_name = f"Test Field Mapping {uuid.uuid4().hex[:8]}"

        # This test verifies that our field mapping (project_id -> "project") works
        # If the mapping is wrong, this will fail with an API error
        result = create_mattergen(
            project_id=self.project_id,
            name=exploration_name,
            batch_size=4,
            num_batches=3,
            diffusion_guidance_factor=15,
        )

        # Verify all fields were set correctly
        assert result["name"] == exploration_name
        assert result["batch_size"] == 4
        assert result["num_batches"] == 3
        assert result["diffusion_guidance_factor"] == 15

        # Verify task has correct project reference
        assert result["task"]["project"]["id"] == self.project_id

    def test_associate_user_upload_with_mattergen_fields(self):
        """Test that associate function uses correct field names"""
        # This test verifies the field mapping for the association function
        # We don't actually perform the association since we'd need valid file uploads
        # but we can test that the function calls the API with correct field names

        # Create a mock exploration first
        exploration_name = f"Test Association Fields {uuid.uuid4().hex[:8]}"
        exploration = create_mattergen(
            project_id=self.project_id, name=exploration_name
        )

        exploration_id = exploration["id"]

        # Test that the function fails gracefully with invalid IDs
        # This verifies the field names are correct (user_upload_id, exploration_id)
        with pytest.raises((APIValidationError, PermissionDenied)):
            associate_user_upload_with_mattergen(
                user_upload_id="non-existent-upload", exploration_id=exploration_id
            )


@pytest.mark.integration
class TestMattergenLaunchIntegration:
    """Tests for LAUNCH functionality - requires cluster configuration"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test environment"""
        self.project_id = os.environ.get("TEST_PROJECT_ID")
        if not self.project_id:
            pytest.skip("TEST_PROJECT_ID not set")

    @pytest.mark.skipif(
        not os.environ.get("TEST_CLUSTER_ID"),
        reason="TEST_CLUSTER_ID not set - skipping LAUNCH tests",
    )
    def test_create_mattergen_launch(self):
        """Test creating and launching a MatterGen exploration"""
        cluster_id = os.environ.get("TEST_CLUSTER_ID")
        exploration_name = f"Test MatterGen Launch {uuid.uuid4().hex[:8]}"

        result = create_mattergen(
            project_id=self.project_id,
            name=exploration_name,
            action="LAUNCH",
            extra_kwargs={"selected_cluster": cluster_id},
        )

        # Verify exploration was created and launched
        assert result["name"] == exploration_name
        assert result["task"]["status"] in [1, 2]  # READY or RUNNING status

        # Note: In real scenarios, you might want to clean up launched jobs
        # but for integration tests, we let them run or fail naturally


@pytest.mark.integration
class TestMattergenEdgeCases:
    """Test edge cases and error conditions"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test environment"""
        self.project_id = os.environ.get("TEST_PROJECT_ID")
        if not self.project_id:
            pytest.skip("TEST_PROJECT_ID not set")

    def test_very_long_name(self):
        """Test exploration with very long name"""
        long_name = "Very Long Exploration Name " + "X" * 200

        # This might succeed or fail depending on backend field limits
        # The test verifies our SDK handles the response appropriately
        try:
            result = create_mattergen(project_id=self.project_id, name=long_name)
            # If successful, verify the name was stored
            assert result["name"] == long_name
        except APIValidationError:
            # If validation fails, that's expected behavior
            pass

    def test_unicode_characters(self):
        """Test exploration with unicode characters"""
        unicode_name = f"Test MatterGen 🧪⚛️ {uuid.uuid4().hex[:8]}"

        result = create_mattergen(
            project_id=self.project_id,
            name=unicode_name,
            description="Testing unicode: 数学 αβγ 🔬",
        )

        assert result["name"] == unicode_name
        assert "数学" in result["description"]

    def test_extreme_batch_values(self):
        """Test with large batch values"""
        exploration_name = f"Test Large Batches {uuid.uuid4().hex[:8]}"

        result = create_mattergen(
            project_id=self.project_id,
            name=exploration_name,
            batch_size=128,
            num_batches=10,
        )

        assert result["batch_size"] == 128
        assert result["num_batches"] == 10
