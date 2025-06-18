"""
Integration tests for EA exploration methods.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid test project_id (set TEST_PROJECT_ID env var)
- A valid test structure_id (set TEST_STRUCTURE_ID env var)

To run: uv run pytest tests/integration/test_ea.py -v -m integration
"""

import os
import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.simulation.ea import (
    create_ea_exploration,
    delete_ea_exploration,
    get_ea_exploration,
    create_ea_exploration_analysis,
    delete_ea_exploration_analysis,
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
    structure_id = os.getenv("TEST_STRUCTURE_ID")
    if not structure_id:
        pytest.skip("TEST_STRUCTURE_ID environment variable must be set")
    return structure_id


@pytest.fixture
def cleanup_explorations():
    """Fixture to track created explorations for cleanup"""
    created_explorations = []

    yield created_explorations

    # Cleanup any explorations that were created during tests
    for exploration_id in created_explorations:
        try:
            delete_ea_exploration(exploration_id)
        except Exception as e:
            print(f"Failed to cleanup exploration {exploration_id}: {e}")


@pytest.fixture
def cleanup_analyses():
    """Fixture to track created analyses for cleanup"""
    created_analyses = []

    yield created_analyses

    # Cleanup any analyses that were created during tests
    for analysis_id in created_analyses:
        try:
            delete_ea_exploration_analysis(analysis_id)
        except Exception as e:
            print(f"Failed to cleanup analysis {analysis_id}: {e}")


@pytest.mark.integration
class TestEAExplorationIntegration:
    """Integration tests for EA exploration operations"""

    def test_create_and_delete_ea_exploration_basic(
        self, test_project_id, test_structure_id, cleanup_explorations
    ):
        """Test basic EA exploration creation and deletion"""

        # Create EA exploration
        exploration_data = create_ea_exploration(
            project_id=test_project_id,
            name="Integration Test EA Exploration",
            structure_id=test_structure_id,
            structure_type="userupload",
            description="Created by integration test",
        )

        assert "id" in exploration_data
        assert exploration_data["name"] == "Integration Test EA Exploration"
        assert exploration_data["description"] == "Created by integration test"

        exploration_id = exploration_data["id"]
        cleanup_explorations.append(exploration_id)

        # Verify we can retrieve the exploration
        retrieved = get_ea_exploration(exploration_id)
        assert retrieved["id"] == exploration_id
        assert retrieved["name"] == "Integration Test EA Exploration"

        # Delete the exploration
        delete_ea_exploration(exploration_id)
        # Remove from cleanup list since we deleted it
        cleanup_explorations.remove(exploration_id)

        # Verify deletion worked (should raise exception or return error)
        with pytest.raises(Exception):
            get_ea_exploration(exploration_id)

    def test_create_ea_exploration_with_launch_action(
        self, test_project_id, test_structure_id, cleanup_explorations
    ):
        """Test EA exploration creation with LAUNCH action"""

        exploration_data = create_ea_exploration(
            project_id=test_project_id,
            name="Launch Test EA Exploration",
            structure_id=test_structure_id,
            action="LAUNCH",
            strains_list=[-0.05, -0.025, 0.025, 0.05],
            stress_algorithm=1,  # OHESS
            stress_method=0,  # Static
            calculator=0,  # FHI-AIMS
        )

        assert "id" in exploration_data
        exploration_id = exploration_data["id"]
        cleanup_explorations.append(exploration_id)

        # Check that task was created (should have task field)
        retrieved = get_ea_exploration(exploration_id)
        assert "task" in retrieved

    def test_create_ea_exploration_with_custom_parameters(
        self, test_project_id, test_structure_id, cleanup_explorations
    ):
        """Test EA exploration creation with various parameter combinations"""

        exploration_data = create_ea_exploration(
            project_id=test_project_id,
            name="Custom Params EA Test",
            structure_id=test_structure_id,
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

        assert "id" in exploration_data
        exploration_id = exploration_data["id"]
        cleanup_explorations.append(exploration_id)

        # Verify parameters were set correctly
        retrieved = get_ea_exploration(exploration_id)
        assert retrieved["strains_list"] == [-0.08, -0.04, 0.04, 0.08, 0.12]
        assert retrieved["stress_algorithm"] == 2
        assert retrieved["stress_method"] == 1
        assert retrieved["num_last_samples"] == 500
        assert retrieved["make_conventional_cell"] is True
        assert retrieved["remove_spurious_distortions"] is False
        assert retrieved["add_vacuum"] == 15
        assert retrieved["is_ht"] is False

    def test_ea_exploration_analysis_workflow(
        self, test_project_id, test_structure_id, cleanup_explorations, cleanup_analyses
    ):
        """Test complete EA exploration analysis workflow"""

        # Create exploration first
        exploration_data = create_ea_exploration(
            project_id=test_project_id,
            name="Analysis Test EA Exploration",
            structure_id=test_structure_id,
        )

        exploration_id = exploration_data["id"]
        cleanup_explorations.append(exploration_id)

        # Create analysis
        analysis_data = create_ea_exploration_analysis(
            exploration_id=exploration_id,
            compute_directional_properties=True,
            action="DRAFT",
        )

        assert "id" in analysis_data
        analysis_id = analysis_data["id"]
        cleanup_analyses.append(analysis_id)

        # Delete analysis
        delete_ea_exploration_analysis(analysis_id)
        cleanup_analyses.remove(analysis_id)

    def test_multiple_structure_types(
        self, test_project_id, test_structure_id, cleanup_explorations
    ):
        """Test EA exploration creation with different structure types"""

        # Test with userupload structure (default)
        exploration_data = create_ea_exploration(
            project_id=test_project_id,
            name="UserUpload Structure Test",
            structure_id=test_structure_id,
            structure_type="userupload",
        )

        assert "id" in exploration_data
        cleanup_explorations.append(exploration_data["id"])

        # Verify structure reference is set correctly
        retrieved = get_ea_exploration(exploration_data["id"])
        assert "starting_structure_userupload" in retrieved
        assert retrieved["starting_structure_userupload"] is not None

    def test_ea_exploration_validation_errors(self, test_project_id, test_structure_id):
        """Test that API validation errors are properly handled"""

        # Test with invalid project_id
        with pytest.raises(Exception):
            create_ea_exploration(
                project_id="invalid-project-id",
                name="Should Fail",
                structure_id=test_structure_id,
            )

        # Test with invalid structure_id
        with pytest.raises(Exception):
            create_ea_exploration(
                project_id=test_project_id,
                name="Should Fail",
                structure_id="invalid-structure-id",
            )
