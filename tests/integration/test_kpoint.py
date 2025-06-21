"""
Integration tests for K-point convergence methods.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)
- A valid test structure_id (set TEST_RELAXED_STRUCTURE_ID env var)

To run: PYTHONPATH=/path/to/atomic_cli uv run pytest tests/integration/test_kpoint.py -v -m integration
"""

import os

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.simulation.kpoint import (
    STRUCTURE_FIELD_MAP,
    create_kpoint_analysis,
    create_kpoint_exploration,
    delete_kpoint_analysis,
    delete_kpoint_exploration,
    delete_kpoint_simulation,
    get_kpoint_analysis,
    get_kpoint_exploration,
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
def cleanup_kpoint_explorations():
    """Track created K-point explorations for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for exploration_id in created_ids:
        try:
            delete_kpoint_exploration(exploration_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def cleanup_kpoint_analyses():
    """Track created K-point analyses for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for analysis_id in created_ids:
        try:
            delete_kpoint_analysis(analysis_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestKPointExplorationIntegration:
    """Integration tests for K-point exploration CRUD operations"""

    def test_create_kpoint_exploration_draft(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_kpoint_explorations,
    ):
        """Test creating a K-point exploration in DRAFT mode"""

        result = create_kpoint_exploration(
            project_id=test_project_id,
            name="Integration Test K-point Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
        )

        # Track for cleanup
        cleanup_kpoint_explorations.append(result["id"])

        # Verify response structure
        assert "id" in result
        assert result["name"] == "Integration Test K-point Exploration"
        assert result["task"]["status"] == 0  # DRAFT status

        # Verify default parameters
        assert result["k_point_range_lower"] == 3
        assert result["k_point_range_upper"] == 6
        assert result["evenly_spaced_kpoints"] is False

        # Verify structure reference is set correctly
        if test_structure_type == "mlrelax":
            assert "starting_structure_mlrelax" in result
            assert result["starting_structure_mlrelax"] is not None
        elif test_structure_type == "userupload":
            assert "starting_structure_userupload" in result
            assert result["starting_structure_userupload"] is not None
        elif test_structure_type == "fhiaims":
            assert "starting_structure" in result
            assert result["starting_structure"] is not None

        # Test getting the created exploration
        retrieved = get_kpoint_exploration(result["id"])
        assert retrieved["id"] == result["id"]
        assert retrieved["name"] == "Integration Test K-point Exploration"

    def test_create_kpoint_exploration_with_custom_parameters(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_kpoint_explorations,
    ):
        """Test K-point exploration creation with custom parameters"""

        result = create_kpoint_exploration(
            project_id=test_project_id,
            name="Custom K-point Parameters Test",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            k_point_range_lower=2,
            k_point_range_upper=8,
            evenly_spaced_kpoints=True,
            action="DRAFT",
        )

        # Track for cleanup
        cleanup_kpoint_explorations.append(result["id"])

        # Verify custom parameters were set correctly
        retrieved = get_kpoint_exploration(result["id"])
        assert retrieved["k_point_range_lower"] == 2
        assert retrieved["k_point_range_upper"] == 8
        assert retrieved["evenly_spaced_kpoints"] is True

    def test_create_kpoint_exploration_with_launch_action(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        test_cluster_id,
        cleanup_kpoint_explorations,
    ):
        """Test K-point exploration creation with LAUNCH action"""

        result = create_kpoint_exploration(
            project_id=test_project_id,
            name="Launch Test K-point Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="LAUNCH",
            k_point_range_lower=4,
            k_point_range_upper=7,
            extra_kwargs={"selected_cluster": test_cluster_id},
        )

        # Track for cleanup
        cleanup_kpoint_explorations.append(result["id"])

        # Check that task was created (should have task field)
        retrieved = get_kpoint_exploration(result["id"])
        assert "task" in retrieved

    def test_create_and_delete_kpoint_exploration(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test creating and then deleting a K-point exploration"""

        # Create exploration
        create_result = create_kpoint_exploration(
            project_id=test_project_id,
            name="Delete Test K-point Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
        )

        exploration_id = create_result["id"]

        # Verify it exists
        get_result = get_kpoint_exploration(exploration_id)
        assert get_result["id"] == exploration_id

        # Delete it
        delete_result = delete_kpoint_exploration(exploration_id)

        # Verify it's deleted (should raise exception or return 404)
        with pytest.raises(Exception):
            get_kpoint_exploration(exploration_id)

    def test_kpoint_exploration_different_structure_types(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_kpoint_explorations,
    ):
        """Test K-point exploration creation with the configured structure type"""

        # Test with configured structure type
        result = create_kpoint_exploration(
            project_id=test_project_id,
            name=f"{test_structure_type.title()} Structure K-point Test",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
        )

        # Track for cleanup
        cleanup_kpoint_explorations.append(result["id"])

        # Verify structure reference is set correctly
        retrieved = get_kpoint_exploration(result["id"])
        if test_structure_type == "mlrelax":
            assert "starting_structure_mlrelax" in retrieved
            assert retrieved["starting_structure_mlrelax"] is not None
        elif test_structure_type == "userupload":
            assert "starting_structure_userupload" in retrieved
            assert retrieved["starting_structure_userupload"] is not None
        elif test_structure_type == "fhiaims":
            assert "starting_structure" in retrieved
            assert retrieved["starting_structure"] is not None

    def test_kpoint_exploration_validation_errors(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test that validation errors are properly raised"""

        # Test invalid action
        with pytest.raises(Exception):
            create_kpoint_exploration(
                project_id=test_project_id,
                name="Invalid Action",
                structure_id=test_structure_id,
                structure_type=test_structure_type,
                action="INVALID",
            )

        # Test invalid structure type
        with pytest.raises(Exception):
            create_kpoint_exploration(
                project_id=test_project_id,
                name="Invalid Structure Type",
                structure_id=test_structure_id,
                structure_type="invalid_type",
            )

        # Test invalid k-point range
        with pytest.raises(Exception):
            create_kpoint_exploration(
                project_id=test_project_id,
                name="Invalid K-point Range",
                structure_id=test_structure_id,
                structure_type=test_structure_type,
                k_point_range_lower=8,
                k_point_range_upper=4,  # Lower > upper
            )

    def test_kpoint_analysis_workflow(
        self,
        test_project_id,
        test_structure_id,
        test_structure_type,
        cleanup_kpoint_explorations,
        cleanup_kpoint_analyses,
    ):
        """Test complete K-point exploration analysis workflow"""

        # Create exploration first
        exploration_result = create_kpoint_exploration(
            project_id=test_project_id,
            name="Analysis Test K-point Exploration",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
        )

        exploration_id = exploration_result["id"]
        cleanup_kpoint_explorations.append(exploration_id)

        # Create analysis
        analysis_result = create_kpoint_analysis(
            exploration_id=exploration_id,
        )

        assert "id" in analysis_result
        analysis_id = analysis_result["id"]
        cleanup_kpoint_analyses.append(analysis_id)

        # Verify analysis was created
        retrieved_analysis = get_kpoint_analysis(analysis_id)
        assert retrieved_analysis["id"] == analysis_id

        # Delete analysis
        delete_kpoint_analysis(analysis_id)
        cleanup_kpoint_analyses.remove(analysis_id)

        # Verify analysis is deleted
        with pytest.raises(Exception):
            get_kpoint_analysis(analysis_id)


@pytest.mark.integration
class TestKPointWorkflow:
    """End-to-end workflow tests"""

    def test_complete_kpoint_workflow(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test complete K-point exploration workflow"""

        # Step 1: Create exploration
        create_result = create_kpoint_exploration(
            project_id=test_project_id,
            name="Workflow Test K-point",
            structure_id=test_structure_id,
            structure_type=test_structure_type,
            action="DRAFT",
            k_point_range_lower=2,
            k_point_range_upper=5,
            evenly_spaced_kpoints=True,
        )

        exploration_id = create_result["id"]
        assert create_result["task"]["status"] == 0  # DRAFT status

        # Step 2: Get exploration details
        get_result = get_kpoint_exploration(exploration_id)
        assert get_result["id"] == exploration_id
        assert get_result["name"] == "Workflow Test K-point"

        # Step 3: Verify parameters
        assert get_result["k_point_range_lower"] == 2
        assert get_result["k_point_range_upper"] == 5
        assert get_result["evenly_spaced_kpoints"] is True

        # Step 4: Create analysis
        analysis_result = create_kpoint_analysis(exploration_id=exploration_id)
        analysis_id = analysis_result["id"]

        # Step 5: Verify analysis
        analysis_get = get_kpoint_analysis(analysis_id)
        assert analysis_get["id"] == analysis_id

        # Step 6: Clean up analysis first
        delete_kpoint_analysis(analysis_id)

        # Step 7: Clean up exploration
        delete_result = delete_kpoint_exploration(exploration_id)

        # Step 8: Verify deletion (should raise exception or return 404)
        with pytest.raises(Exception):
            get_kpoint_exploration(exploration_id)

    def test_structure_field_mapping_validation(
        self, test_project_id, test_structure_id, test_structure_type
    ):
        """Test that structure field mapping works correctly"""

        # Verify all structure types in mapping are valid
        for structure_type in STRUCTURE_FIELD_MAP.keys():
            if structure_type == test_structure_type:
                # Test the configured structure type
                result = create_kpoint_exploration(
                    project_id=test_project_id,
                    name=f"Field Mapping Test {structure_type}",
                    structure_id=test_structure_id,
                    structure_type=structure_type,
                )

                # Clean up
                delete_kpoint_exploration(result["id"])

                # If we get here, the mapping worked
                assert True
