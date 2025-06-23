import pytest
from unittest.mock import patch
from atomict.simulation.mattergen import (
    create_mattergen,
    DEFAULT_BATCH_SIZE,
    DEFAULT_NUM_BATCHES,
)
from atomict.exceptions import APIValidationError


class TestCreateMattergen:

    def test_valid_draft_creation(self):
        """Test creating a MatterGen exploration with valid parameters"""
        with patch("atomict.simulation.mattergen.post") as mock_post:
            mock_post.return_value = {"id": "test-id", "name": "test-exploration"}

            create_mattergen(
                project_id="test-project",
                name="Test Exploration",
                description="Test description",
            )

            mock_post.assert_called_once_with(
                "api/mattergen-exploration/",
                payload={
                    "project": "test-project",
                    "name": "Test Exploration",
                    "description": "Test description",
                    "batch_size": DEFAULT_BATCH_SIZE,
                    "num_batches": DEFAULT_NUM_BATCHES,
                    "action": "DRAFT",
                },
            )

    def test_valid_launch_creation(self):
        """Test creating a MatterGen exploration with LAUNCH action"""
        with patch("atomict.simulation.mattergen.post") as mock_post:
            mock_post.return_value = {"id": "test-id", "name": "test-exploration"}

            create_mattergen(
                project_id="test-project", name="Test Exploration", action="LAUNCH"
            )

            expected_payload = {
                "project": "test-project",
                "name": "Test Exploration",
                "description": "",
                "batch_size": DEFAULT_BATCH_SIZE,
                "num_batches": DEFAULT_NUM_BATCHES,
                "action": "LAUNCH",
            }
            mock_post.assert_called_once_with(
                "api/mattergen-exploration/", payload=expected_payload
            )

    def test_custom_batch_parameters(self):
        """Test creating exploration with custom batch parameters"""
        with patch("atomict.simulation.mattergen.post") as mock_post:
            mock_post.return_value = {"id": "test-id"}

            create_mattergen(
                project_id="test-project",
                name="Test Exploration",
                batch_size=32,
                num_batches=5,
            )

            call_payload = mock_post.call_args[1]["payload"]
            assert call_payload["batch_size"] == 32
            assert call_payload["num_batches"] == 5

    def test_diffusion_guidance_factor(self):
        """Test creating exploration with diffusion guidance factor"""
        with patch("atomict.simulation.mattergen.post") as mock_post:
            mock_post.return_value = {"id": "test-id"}

            create_mattergen(
                project_id="test-project",
                name="Test Exploration",
                diffusion_guidance_factor=50,
            )

            call_payload = mock_post.call_args[1]["payload"]
            assert call_payload["diffusion_guidance_factor"] == 50

    def test_diffusion_guidance_factor_not_included_when_none(self):
        """Test that diffusion_guidance_factor is not included when None"""
        with patch("atomict.simulation.mattergen.post") as mock_post:
            mock_post.return_value = {"id": "test-id"}

            create_mattergen(
                project_id="test-project",
                name="Test Exploration",
                diffusion_guidance_factor=None,
            )

            call_payload = mock_post.call_args[1]["payload"]
            assert "diffusion_guidance_factor" not in call_payload

    def test_extra_kwargs(self):
        """Test passing extra parameters via extra_kwargs"""
        with patch("atomict.simulation.mattergen.post") as mock_post:
            mock_post.return_value = {"id": "test-id"}

            create_mattergen(
                project_id="test-project",
                name="Test Exploration",
                extra_kwargs={"selected_cluster": "cluster-1"},
            )

            call_payload = mock_post.call_args[1]["payload"]
            assert call_payload["selected_cluster"] == "cluster-1"

    def test_project_field_mapping(self):
        """Test that project_id maps to 'project' field for LaunchableViewSet"""
        with patch("atomict.simulation.mattergen.post") as mock_post:
            mock_post.return_value = {"id": "test-id"}

            create_mattergen(project_id="test-project-123", name="Test Exploration")

            call_payload = mock_post.call_args[1]["payload"]
            assert call_payload["project"] == "test-project-123"
            assert "project_id" not in call_payload


class TestCreateMattergenValidation:

    def test_invalid_action(self):
        """Test validation error for invalid action"""
        with pytest.raises(
            APIValidationError, match="Action must be 'DRAFT' or 'LAUNCH'"
        ):
            create_mattergen(
                project_id="test-project", name="Test Exploration", action="INVALID"
            )

    def test_zero_batch_size(self):
        """Test validation error for zero batch_size"""
        with pytest.raises(APIValidationError, match="batch_size must be positive"):
            create_mattergen(
                project_id="test-project", name="Test Exploration", batch_size=0
            )

    def test_negative_batch_size(self):
        """Test validation error for negative batch_size"""
        with pytest.raises(APIValidationError, match="batch_size must be positive"):
            create_mattergen(
                project_id="test-project", name="Test Exploration", batch_size=-1
            )

    def test_zero_num_batches(self):
        """Test validation error for zero num_batches"""
        with pytest.raises(APIValidationError, match="num_batches must be positive"):
            create_mattergen(
                project_id="test-project", name="Test Exploration", num_batches=0
            )

    def test_negative_num_batches(self):
        """Test validation error for negative num_batches"""
        with pytest.raises(APIValidationError, match="num_batches must be positive"):
            create_mattergen(
                project_id="test-project", name="Test Exploration", num_batches=-5
            )

    def test_zero_diffusion_guidance_factor(self):
        """Test validation error for zero diffusion_guidance_factor"""
        with pytest.raises(
            APIValidationError, match="diffusion_guidance_factor must be positive"
        ):
            create_mattergen(
                project_id="test-project",
                name="Test Exploration",
                diffusion_guidance_factor=0,
            )

    def test_negative_diffusion_guidance_factor(self):
        """Test validation error for negative diffusion_guidance_factor"""
        with pytest.raises(
            APIValidationError, match="diffusion_guidance_factor must be positive"
        ):
            create_mattergen(
                project_id="test-project",
                name="Test Exploration",
                diffusion_guidance_factor=-10,
            )


class TestConstants:

    def test_default_constants(self):
        """Test that default constants are set correctly"""
        assert DEFAULT_BATCH_SIZE == 16
        assert DEFAULT_NUM_BATCHES == 1
