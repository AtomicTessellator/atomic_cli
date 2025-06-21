"""
Integration tests for SQS exploration methods.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)

To run: uv run pytest tests/integration/test_sqs.py -v -m integration
"""

import os
from pathlib import Path

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.simulation.sqs import (
    create_sqs_exploration,
    delete_sqs_exploration,
    get_simulation,
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


@pytest.fixture(scope="session")
def test_structure_id():
    """Get test structure ID from environment (FHI-aims, MLRelax, or UserUpload ID)"""
    structure_id = os.getenv("TEST_STRUCTURE_ID")
    if not structure_id:
        pytest.skip("TEST_STRUCTURE_ID environment variable must be set")
    return structure_id


@pytest.fixture(scope="session")
def test_structure_type():
    """Get structure type from environment (defaults to fhiaims)"""
    return os.getenv("TEST_STRUCTURE_TYPE", "fhiaims")


@pytest.fixture
def cleanup_sqs_explorations():
    """Track created SQS explorations for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for sqs_id in created_ids:
        try:
            delete_sqs_exploration(sqs_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestSQSExplorationIntegration:
    """Integration tests for SQS exploration CRUD operations"""

    def test_create_sqs_exploration_draft(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_sqs_explorations,
    ):
        """Test creating a SQS exploration in DRAFT mode"""
        target_concentrations = [
            {"element": "Al", "weight": 0.75},
            {"element": "Ni", "weight": 0.25},
        ]

        result = create_sqs_exploration(
            project_id=test_project_id,
            name="Integration Test SQS Al-Ni",
            target_concentrations=target_concentrations,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
            description="Integration test SQS exploration for Al-Ni alloy",
        )

        # Track for cleanup
        cleanup_sqs_explorations.append(result["id"])

        # Verify response structure
        assert "id" in result
        assert result["name"] == "Integration Test SQS Al-Ni"
        assert (
            result["description"] == "Integration test SQS exploration for Al-Ni alloy"
        )
        assert result["max_size"] == 8
        assert result["cluster_cutoffs"] == [4.0, 4.0]
        assert result["task"]["status"] == 0  # DRAFT status

        # Verify target concentrations
        assert len(result["target_concentrations"]) == 2
        concentration_elements = {c["element"] for c in result["target_concentrations"]}
        assert concentration_elements == {"Al", "Ni"}

        # Verify weights
        weights = {c["element"]: c["weight"] for c in result["target_concentrations"]}
        assert weights["Al"] == 0.75
        assert weights["Ni"] == 0.25

        # Test getting the created exploration
        retrieved = get_simulation(result["id"])
        assert retrieved["id"] == result["id"]
        assert retrieved["name"] == "Integration Test SQS Al-Ni"

    def test_create_sqs_exploration_with_custom_params(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_sqs_explorations,
    ):
        """Test creating SQS exploration with custom parameters"""
        target_concentrations = [
            {"element": "Fe", "weight": 0.5},
            {"element": "Cr", "weight": 0.5},
        ]

        result = create_sqs_exploration(
            project_id=test_project_id,
            name="Custom Params SQS",
            target_concentrations=target_concentrations,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
            description="Custom parameters test",
            max_size=12,
            cluster_cutoffs=[3.5, 4.5],
        )

        # Track for cleanup
        cleanup_sqs_explorations.append(result["id"])

        # Verify custom parameters
        assert result["max_size"] == 12
        assert result["cluster_cutoffs"] == [3.5, 4.5]
        assert result["name"] == "Custom Params SQS"

    def test_create_and_delete_sqs_exploration(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test creating and then deleting a SQS exploration"""
        target_concentrations = [
            {"element": "Ti", "weight": 0.6},
            {"element": "Al", "weight": 0.4},
        ]

        # Create exploration
        create_result = create_sqs_exploration(
            project_id=test_project_id,
            name="Delete Test SQS",
            target_concentrations=target_concentrations,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
        )

        exploration_id = create_result["id"]

        # Verify it exists
        get_result = get_simulation(exploration_id)
        assert get_result["id"] == exploration_id

        # Delete it
        delete_result = delete_sqs_exploration(exploration_id)
        # Delete endpoint may return empty response or success message

        # Verify it's deleted (should raise exception or return 404)
        with pytest.raises(Exception):
            get_simulation(exploration_id)

    def test_get_simulation_with_params(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_sqs_explorations,
    ):
        """Test getting SQS simulation with additional parameters"""
        target_concentrations = [{"element": "Cu", "weight": 1.0}]

        # Create exploration first
        create_result = create_sqs_exploration(
            project_id=test_project_id,
            name="Params Test SQS",
            target_concentrations=target_concentrations,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
        )

        cleanup_sqs_explorations.append(create_result["id"])
        exploration_id = create_result["id"]

        # Test getting with full=True
        full_result = get_simulation(exploration_id, full=True)
        assert full_result["id"] == exploration_id

        # Test getting with additional parameters
        param_result = get_simulation(exploration_id, full=True, expand="task")
        assert param_result["id"] == exploration_id

    def test_create_sqs_exploration_validation_errors(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test that validation errors are properly raised"""

        # Test invalid action
        with pytest.raises(ValueError, match="Action must be 'DRAFT' or 'LAUNCH'"):
            create_sqs_exploration(
                project_id=test_project_id,
                name="Invalid Action",
                target_concentrations=[{"element": "Cu", "weight": 1.0}],
                structure_id=test_structure_id,
                structure_type=test_structure_type,
                action="INVALID",
            )

        # Test invalid structure type
        with pytest.raises(ValueError, match="structure_type must be one of"):
            create_sqs_exploration(
                project_id=test_project_id,
                name="Invalid Structure Type",
                target_concentrations=[{"element": "Cu", "weight": 1.0}],
                structure_id=test_structure_id,
                structure_type="invalid_type",
            )

        # Test invalid concentration sum
        with pytest.raises(
            ValueError, match="Sum of target concentrations must equal 1.0"
        ):
            create_sqs_exploration(
                project_id=test_project_id,
                name="Invalid Sum",
                target_concentrations=[
                    {"element": "Cu", "weight": 0.3},
                    {"element": "Zn", "weight": 0.3},  # Sum = 0.6
                ],
                structure_id=test_structure_id,
                structure_type=test_structure_type,
            )

    def test_create_sqs_exploration_different_structure_types(
        self, test_project_id, test_structure_id, cleanup_sqs_explorations
    ):
        """Test creating SQS explorations with different structure types"""
        target_concentrations = [{"element": "Fe", "weight": 1.0}]

        # Test userupload type (default)
        result_userupload = create_sqs_exploration(
            project_id=test_project_id,
            name="UseruploadType Test SQS",
            target_concentrations=target_concentrations,
            structure_id=test_structure_id,
            structure_type="userupload",
            action="DRAFT",
        )

        cleanup_sqs_explorations.append(result_userupload["id"])
        assert result_userupload["name"] == "UseruploadType Test SQS"


@pytest.mark.integration
class TestSQSWorkflow:
    """End-to-end workflow tests"""

    def test_complete_sqs_workflow(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test complete SQS exploration workflow"""
        target_concentrations = [
            {"element": "Mg", "weight": 0.3},
            {"element": "Al", "weight": 0.7},
        ]

        # Step 1: Create exploration
        create_result = create_sqs_exploration(
            project_id=test_project_id,
            name="Workflow Test SQS",
            target_concentrations=target_concentrations,
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
            description="Complete workflow test",
        )

        exploration_id = create_result["id"]
        assert create_result["task"]["status"] == 0  # DRAFT status

        # Step 2: Get exploration details
        get_result = get_simulation(exploration_id, full=True)
        assert get_result["id"] == exploration_id
        assert get_result["name"] == "Workflow Test SQS"

        # Step 3: Verify target concentrations
        concentrations = get_result["target_concentrations"]
        assert len(concentrations) == 2
        elements = {c["element"]: c["weight"] for c in concentrations}
        assert elements["Mg"] == 0.3
        assert elements["Al"] == 0.7

        # Step 4: Clean up
        delete_result = delete_sqs_exploration(exploration_id)

        # Step 5: Verify deletion (should raise exception or return 404)
        with pytest.raises(Exception):
            get_simulation(exploration_id)
