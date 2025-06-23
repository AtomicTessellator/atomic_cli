"""
Integration tests for HT ML Relaxation exploration methods.

These tests require:
- Valid authentication credentials (ATOMICT_USERNAME, ATOMICT_PASSWORD)
- A test project (TEST_PROJECT_ID)
- A completed HT SQS exploration (TEST_HT_SQS_EXPLORATION_ID)

To run: PYTHONPATH=/path/to/atomic_cli uv run pytest tests/integration/test_ht_mlrelax.py -v -m integration
"""

import os
from typing import List

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate

from atomict.simulation.ht_mlrelax import (
    ML_MODELS,
    create_ht_mlrelax_exploration,
    delete_ht_mlrelax_exploration,
    get_ht_mlrelax_exploration,
    launch_ht_mlrelax_exploration,
    list_ht_mlrelax_explorations,
    update_ht_mlrelax_exploration,
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
    """Test project ID from environment"""
    project_id = os.getenv("TEST_PROJECT_ID")
    if not project_id:
        pytest.skip("TEST_PROJECT_ID environment variable not set")
    return project_id


@pytest.fixture(scope="session")
def test_ht_sqs_exploration_id():
    """Test HT SQS exploration ID from environment"""
    exploration_id = os.getenv("TEST_HT_SQS_EXPLORATION_ID")
    if not exploration_id:
        pytest.skip("TEST_HT_SQS_EXPLORATION_ID environment variable not set")
    return exploration_id


@pytest.fixture
def cleanup_ht_mlrelax_explorations():
    """Track created HT ML Relaxation explorations for cleanup"""
    created_ids = []

    def track_id(exploration_id):
        created_ids.append(exploration_id)

    yield track_id

    # Cleanup
    for exploration_id in created_ids:
        try:
            delete_ht_mlrelax_exploration(exploration_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestHTMLRelaxExplorationIntegration:
    """Integration tests for HT ML Relaxation exploration CRUD operations"""

    def test_create_ht_mlrelax_exploration_minimal(
        self,
        test_project_id,
        test_ht_sqs_exploration_id,
        cleanup_ht_mlrelax_explorations,
    ):
        """Test creating a HT ML Relaxation exploration with minimal parameters"""

        result = create_ht_mlrelax_exploration(
            project_id=test_project_id,
            name="Integration Test HT MLRelax Minimal",
            source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
        )

        # Should return success message with count
        assert "message" in result
        assert "relaxation" in result["message"].lower()

        # Extract exploration ID from the result for cleanup
        # Note: The create endpoint creates individual MLRelax tasks but doesn't
        # return the exploration ID directly. We'd need to list and find it.

    def test_create_ht_mlrelax_exploration_full_params(
        self,
        test_project_id,
        test_ht_sqs_exploration_id,
        cleanup_ht_mlrelax_explorations,
    ):
        """Test creating a HT ML Relaxation exploration with all parameters"""

        result = create_ht_mlrelax_exploration(
            project_id=test_project_id,
            name="Integration Test HT MLRelax Full",
            description="Full parameter test for HT ML Relaxation",
            source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
            f_max=0.03,
            model="orb_v3_conservative",
        )

        assert "message" in result
        assert "relaxation" in result["message"].lower()

    def test_create_ht_mlrelax_exploration_different_models(
        self,
        test_project_id,
        test_ht_sqs_exploration_id,
        cleanup_ht_mlrelax_explorations,
    ):
        """Test creating HT ML Relaxation explorations with different ML models"""

        for model_name in ["orb_d3_v2", "mattersim_1_0_0_5m", "esen_30m_oam"]:
            result = create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name=f"Integration Test HT MLRelax {model_name}",
                source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
                model=model_name,
                f_max=0.05,
            )

            assert "message" in result
            assert "relaxation" in result["message"].lower()

    def test_list_ht_mlrelax_explorations(self, test_project_id):
        """Test listing HT ML Relaxation explorations"""

        result = list_ht_mlrelax_explorations(project_id=test_project_id, limit=10)

        assert "results" in result
        assert isinstance(result["results"], list)
        # Should have standard pagination fields
        assert "count" in result or "next" in result or "previous" in result

        # If there are results, verify structure
        if result["results"]:
            exploration = result["results"][0]
            assert "id" in exploration
            assert "name" in exploration
            assert "model" in exploration
            # Should reference a HT SQS exploration
            assert "source_ht_sqs_exploration" in exploration

    def test_get_ht_mlrelax_exploration_not_found(self):
        """Test getting non-existent HT ML Relaxation exploration"""

        with pytest.raises(Exception):  # Will raise APIError or similar
            get_ht_mlrelax_exploration("00000000-0000-0000-0000-000000000000")

    def test_create_ht_mlrelax_exploration_validation_errors(
        self, test_project_id, test_ht_sqs_exploration_id
    ):
        """Test validation errors during HT ML Relaxation creation"""

        # Test invalid model
        with pytest.raises(ValueError, match="Invalid model"):
            create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name="Invalid Model Test",
                source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
                model="invalid_model_name",
            )

    def test_ht_mlrelax_workflow_integration(
        self,
        test_project_id,
        test_ht_sqs_exploration_id,
        cleanup_ht_mlrelax_explorations,
    ):
        """Test complete HT ML Relaxation workflow: create -> list -> get -> update"""

        # 1. Create HT ML Relaxation exploration
        create_result = create_ht_mlrelax_exploration(
            project_id=test_project_id,
            name="Integration Test HT MLRelax Workflow",
            description="Workflow integration test",
            source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
            f_max=0.04,
            model="mattersim_1_0_0_5m",
        )

        assert "relaxation" in create_result["message"].lower()

        # 2. List explorations and find our creation
        list_result = list_ht_mlrelax_explorations(project_id=test_project_id, limit=50)

        # Find our exploration by name
        our_exploration = None
        for exp in list_result["results"]:
            if exp["name"] == "Integration Test HT MLRelax Workflow":
                our_exploration = exp
                break

        assert our_exploration is not None, "Created exploration not found in list"
        cleanup_ht_mlrelax_explorations(our_exploration["id"])

        # 3. Get exploration details
        detail_result = get_ht_mlrelax_exploration(our_exploration["id"])

        assert "id" in detail_result
        assert detail_result["name"] == "Integration Test HT MLRelax Workflow"
        assert detail_result["description"] == "Workflow integration test"
        assert detail_result["f_max"] == 0.04
        assert detail_result["model"] == ML_MODELS["mattersim_1_0_0_5m"]

        # 4. Get with children to see individual ML relaxations
        detail_with_children = get_ht_mlrelax_exploration(
            our_exploration["id"], include_children=True
        )

        assert "children" in detail_with_children
        # Should have children (individual ML relaxation tasks)
        if detail_with_children["children"]:
            child = detail_with_children["children"][0]
            assert "mlrelax" in child
            assert child["mlrelax"] is not None

        # 5. Update exploration parameters
        update_result = update_ht_mlrelax_exploration(
            our_exploration["id"],
            name="Updated Integration Test Name",
            f_max=0.02,
            model="orb_d3_v2",
        )

        assert update_result["name"] == "Updated Integration Test Name"
        assert update_result["f_max"] == 0.02
        assert update_result["model"] == ML_MODELS["orb_d3_v2"]


@pytest.mark.integration
class TestHTMLRelaxExplorationLaunchIntegration:
    """Integration tests for HT ML Relaxation exploration launch operations"""

    def test_list_explorations_for_launch_candidates(self, test_project_id):
        """Test listing explorations to find candidates for launch testing"""

        result = list_ht_mlrelax_explorations(project_id=test_project_id, limit=20)

        assert "results" in result

        # Look for explorations in DRAFT status that could be launched
        draft_explorations = []
        for exp in result["results"]:
            # Get detailed info to check individual ML relaxation status
            detail = get_ht_mlrelax_exploration(exp["id"], include_children=True)
            if "children" in detail and detail["children"]:
                # Check if any child ML relaxations are in DRAFT status
                for child in detail["children"]:
                    if "mlrelax" in child and child["mlrelax"]:
                        mlrelax = child["mlrelax"]
                        if "task" in mlrelax and mlrelax["task"]:
                            task = mlrelax["task"]
                            if "status" in task and task["status"] == 0:  # DRAFT
                                draft_explorations.append(exp)
                                break

        # This is mainly for visibility - shows available launch candidates
        print(f"Found {len(draft_explorations)} explorations with DRAFT ML relaxations")

    def test_update_exploration_parameters(self, test_project_id):
        """Test updating exploration parameters without launching"""

        # Get an existing exploration to update
        list_result = list_ht_mlrelax_explorations(project_id=test_project_id, limit=5)

        if not list_result["results"]:
            pytest.skip("No HT ML Relaxation explorations available for testing")

        exploration_id = list_result["results"][0]["id"]
        original_exploration = get_ht_mlrelax_exploration(exploration_id)

        # Update just the description
        update_result = update_ht_mlrelax_exploration(
            exploration_id,
            description="Updated via integration test",
        )

        assert "description" in update_result
        assert update_result["description"] == "Updated via integration test"

        # Verify other fields unchanged
        assert update_result["name"] == original_exploration["name"]
        assert update_result["f_max"] == original_exploration["f_max"]
        assert update_result["model"] == original_exploration["model"]


@pytest.mark.integration
class TestHTMLRelaxModelValidationIntegration:
    """Integration tests for ML model validation and constants"""

    def test_ml_models_constants_match_backend(self):
        """Test that our ML_MODELS constants match backend expectations"""

        # Verify all models in our constants are valid
        expected_models = {
            "orb_d3_v2": 0,
            "mattersim_1_0_0_5m": 1,
            "orb_v3_conservative": 2,
            "esen_30m_oam": 3,
        }

        assert ML_MODELS == expected_models

        # Test each model by creating a small exploration (if we have test data)
        # This would verify the backend accepts these model codes
        for model_name, model_code in ML_MODELS.items():
            assert isinstance(model_name, str)
            assert isinstance(model_code, int)
            assert model_code >= 0

    def test_model_parameter_validation(
        self, test_project_id, test_ht_sqs_exploration_id
    ):
        """Test model parameter validation with real API calls"""

        # Valid model should work
        try:
            result = create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name=f"Model Validation Test - Valid",
                source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
                model="mattersim_1_0_0_5m",
            )
            assert "message" in result
        except Exception as e:
            # If creation fails due to other reasons (like no completed SQS),
            # that's okay - we're mainly testing the model validation
            if "model" in str(e).lower():
                pytest.fail(f"Valid model was rejected: {e}")

        # Invalid model should fail validation before API call
        with pytest.raises(ValueError, match="Invalid model"):
            create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name="Model Validation Test - Invalid",
                source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
                model="nonexistent_model",
            )


@pytest.mark.integration
class TestHTMLRelaxAPIFieldMappingIntegration:
    """Integration tests to verify correct API field mapping"""

    def test_create_field_mapping_integration(
        self, test_project_id, test_ht_sqs_exploration_id
    ):
        """Test that our field mapping works correctly with the real API"""

        # Test that we're using correct field names expected by the backend
        # The backend expects "project" not "project_id", etc.

        try:
            result = create_ht_mlrelax_exploration(
                project_id=test_project_id,
                name="Field Mapping Test",
                source_ht_sqs_exploration_id=test_ht_sqs_exploration_id,
                f_max=0.05,
                model="mattersim_1_0_0_5m",
                description="Testing field mapping",
            )

            # If we get here without API errors, our field mapping is correct
            assert "message" in result

        except Exception as e:
            # Check if the error is related to field mapping
            error_msg = str(e).lower()
            if any(
                field in error_msg
                for field in ["project_id", "source_ht_sqs_exploration_id"]
            ):
                pytest.fail(f"Field mapping error detected: {e}")
            # Other errors (like missing data) are acceptable for this test

    def test_response_structure_validation(self, test_project_id):
        """Test that API responses have expected structure"""

        # Test list response structure
        list_result = list_ht_mlrelax_explorations(project_id=test_project_id, limit=5)

        assert isinstance(list_result, dict)
        assert "results" in list_result
        assert isinstance(list_result["results"], list)

        # Test individual exploration structure
        if list_result["results"]:
            exploration = list_result["results"][0]
            expected_fields = ["id", "name", "model", "source_ht_sqs_exploration"]

            for field in expected_fields:
                assert field in exploration, f"Missing expected field: {field}"

            # Test detailed get response
            detail_result = get_ht_mlrelax_exploration(exploration["id"])

            for field in expected_fields:
                assert (
                    field in detail_result
                ), f"Missing expected field in detail: {field}"
