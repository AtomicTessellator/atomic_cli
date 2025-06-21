"""
Integration tests for EA exploration methods.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)
- A valid test structure_id (set TEST_STRUCTURE_ID env var)

To run: uv run pytest tests/integration/test_ea.py -v -m integration
"""

import os

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.simulation.ea import (
    create_ea_exploration,
    create_ea_exploration_analysis,
    delete_ea_exploration,
    delete_ea_exploration_analysis,
    get_ea_exploration,
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
    """Get test structure ID from environment"""
    structure_id = os.getenv("TEST_RELAXED_STRUCTURE_ID")
    if not structure_id:
        pytest.skip("TEST_RELAXED_STRUCTURE_ID environment variable must be set")
    return structure_id


@pytest.fixture(scope="session")
def test_structure_type():
    """Get structure type from environment (defaults to mlrelax)"""
    return os.getenv("TEST_RELAXED_STRUCTURE_TYPE", "mlrelax")


@pytest.fixture(scope="session")
def test_cluster_id():
    """Get first available cluster for LAUNCH tests"""
    try:
        from atomict.cli.core.client import get_client

        client = get_client()
        clusters = client.get_all("/api/k8s-cluster/")

        if not clusters:
            pytest.skip("No clusters available for LAUNCH tests")

        return clusters[0]["id"]
    except Exception as e:
        pytest.skip(f"Failed to fetch clusters: {e}")


@pytest.fixture
def cleanup_ea_explorations():
    """Track created EA explorations for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for exploration_id in created_ids:
        try:
            delete_ea_exploration(exploration_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def cleanup_ea_analyses():
    """Track created EA analyses for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for analysis_id in created_ids:
        try:
            delete_ea_exploration_analysis(analysis_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestEAExplorationIntegration:
    """Integration tests for EA exploration CRUD operations"""

    def test_create_ea_exploration_draft(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_ea_explorations,
    ):
        """Test creating an EA exploration in DRAFT mode"""

        result = create_ea_exploration(
            project_id=test_project_id,
            name="Integration Test EA Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
            description="Created by integration test",
        )

        # Track for cleanup
        cleanup_ea_explorations.append(result["id"])

        # Verify response structure
        assert "id" in result
        assert result["name"] == "Integration Test EA Exploration"
        assert result["description"] == "Created by integration test"
        assert result["task"]["status"] == 0  # DRAFT status

        # Verify default parameters
        assert result["strains_list"] == [-0.06, -0.03, 0.03, 0.06]
        assert result["stress_algorithm"] == 2  # ASESS
        assert result["stress_method"] == 1  # Dynamic
        assert result["calculator"] == 0  # FHI-AIMS
        assert result["num_last_samples"] == 1000
        assert result["make_conventional_cell"] is False
        assert result["remove_spurious_distortions"] is True

        # Test getting the created exploration
        retrieved = get_ea_exploration(result["id"])
        assert retrieved["id"] == result["id"]
        assert retrieved["name"] == "Integration Test EA Exploration"

    def test_create_ea_exploration_with_custom_parameters(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_ea_explorations,
    ):
        """Test EA exploration creation with custom parameters"""

        result = create_ea_exploration(
            project_id=test_project_id,
            name="Custom Params EA Test",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            strains_list=[-0.08, -0.04, 0.04, 0.08, 0.12],
            stress_algorithm=2,  # ASESS
            stress_method=1,  # Dynamic
            calculator=0,  # FHI-AIMS
            num_last_samples=500,
            make_conventional_cell=True,
            remove_spurious_distortions=False,
            add_vacuum=15,
            is_ht=False,
        )

        # Track for cleanup
        cleanup_ea_explorations.append(result["id"])

        # Verify custom parameters were set correctly
        retrieved = get_ea_exploration(result["id"])
        assert retrieved["strains_list"] == [-0.08, -0.04, 0.04, 0.08, 0.12]
        assert retrieved["stress_algorithm"] == 2
        assert retrieved["stress_method"] == 1
        assert retrieved["num_last_samples"] == 500
        assert retrieved["make_conventional_cell"] is True
        assert retrieved["remove_spurious_distortions"] is False
        assert retrieved["add_vacuum"] == 15
        assert retrieved["is_ht"] is False

    def test_create_ea_exploration_with_launch_action(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        test_cluster_id,
        cleanup_ea_explorations,
    ):
        """Test EA exploration creation with LAUNCH action"""

        result = create_ea_exploration(
            project_id=test_project_id,
            name="Launch Test EA Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="LAUNCH",
            strains_list=[-0.05, -0.025, 0.025, 0.05],
            stress_algorithm=1,  # OHESS
            stress_method=0,  # Static
            calculator=0,  # FHI-AIMS
            extra_kwargs={"selected_cluster": test_cluster_id},
        )

        # Track for cleanup
        cleanup_ea_explorations.append(result["id"])

        # Check that task was created (should have task field)
        retrieved = get_ea_exploration(result["id"])
        assert "task" in retrieved

    def test_create_and_delete_ea_exploration(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test creating and then deleting an EA exploration"""

        # Create exploration
        create_result = create_ea_exploration(
            project_id=test_project_id,
            name="Delete Test EA Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
        )

        exploration_id = create_result["id"]

        # Verify it exists
        get_result = get_ea_exploration(exploration_id)
        assert get_result["id"] == exploration_id

        # Delete it
        delete_result = delete_ea_exploration(exploration_id)

        # Verify it's deleted (should raise exception or return 404)
        with pytest.raises(Exception):
            get_ea_exploration(exploration_id)

    def test_ea_exploration_different_structure_types(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_ea_explorations,
    ):
        """Test EA exploration creation with the configured structure type"""

        # Test with configured structure type (mlrelax)
        result = create_ea_exploration(
            project_id=test_project_id,
            name="MLrelax Structure Test",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
        )

        # Track for cleanup
        cleanup_ea_explorations.append(result["id"])

        # Verify structure reference is set correctly
        retrieved = get_ea_exploration(result["id"])
        if test_structure_type == "mlrelax":
            assert "starting_structure_mlrelax" in retrieved
            assert retrieved["starting_structure_mlrelax"] is not None
        elif test_structure_type == "userupload":
            assert "starting_structure_userupload" in retrieved
            assert retrieved["starting_structure_userupload"] is not None
        elif test_structure_type == "fhiaims":
            assert "starting_structure" in retrieved
            assert retrieved["starting_structure"] is not None

    def test_ea_exploration_validation_errors(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test that validation errors are properly raised"""

        # Test invalid action
        with pytest.raises(ValueError, match="Action must be 'DRAFT' or 'LAUNCH'"):
            create_ea_exploration(
                project_id=test_project_id,
                name="Invalid Action",
                structure_id=test_structure_id,
                structure_type=test_structure_type,
                action="INVALID",
            )

        # Test invalid structure type
        with pytest.raises(ValueError, match="structure_type must be one of"):
            create_ea_exploration(
                project_id=test_project_id,
                name="Invalid Structure Type",
                structure_id=test_structure_id,
                structure_type="invalid_type",
            )

        # Test invalid strains list (too few elements)
        with pytest.raises(
            ValueError, match="strains_list must contain at least 4 elements"
        ):
            create_ea_exploration(
                project_id=test_project_id,
                name="Invalid Strains",
                structure_id=test_structure_id,
                structure_type=test_structure_type,
                strains_list=[0.01, 0.02],  # Only 2 elements
            )

    def test_ea_exploration_analysis_workflow(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_ea_explorations,
        cleanup_ea_analyses,
    ):
        """Test complete EA exploration analysis workflow"""

        # Create exploration first
        exploration_result = create_ea_exploration(
            project_id=test_project_id,
            name="Analysis Test EA Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
        )

        exploration_id = exploration_result["id"]
        cleanup_ea_explorations.append(exploration_id)

        # Create analysis
        analysis_result = create_ea_exploration_analysis(
            exploration_id=exploration_id,
            compute_directional_properties=True,
            action="DRAFT",
        )

        assert "id" in analysis_result
        analysis_id = analysis_result["id"]
        cleanup_ea_analyses.append(analysis_id)

        # Delete analysis
        delete_ea_exploration_analysis(analysis_id)
        cleanup_ea_analyses.remove(analysis_id)


@pytest.mark.integration
class TestEAWorkflow:
    """End-to-end workflow tests"""

    def test_complete_ea_workflow(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test complete EA exploration workflow"""

        # Step 1: Create exploration
        create_result = create_ea_exploration(
            project_id=test_project_id,
            name="Workflow Test EA",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
            description="Complete workflow test",
            strains_list=[-0.04, -0.02, 0.02, 0.04],
            stress_algorithm=2,  # ASESS
            stress_method=1,  # Dynamic
        )

        exploration_id = create_result["id"]
        assert create_result["task"]["status"] == 0  # DRAFT status

        # Step 2: Get exploration details
        get_result = get_ea_exploration(exploration_id)
        assert get_result["id"] == exploration_id
        assert get_result["name"] == "Workflow Test EA"

        # Step 3: Verify parameters
        assert get_result["strains_list"] == [-0.04, -0.02, 0.02, 0.04]
        assert get_result["stress_algorithm"] == 2
        assert get_result["stress_method"] == 1

        # Step 4: Clean up
        delete_result = delete_ea_exploration(exploration_id)

        # Step 5: Verify deletion (should raise exception or return 404)
        with pytest.raises(Exception):
            get_ea_exploration(exploration_id)
