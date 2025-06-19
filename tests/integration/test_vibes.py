"""
Integration tests for VIBES phonon simulation methods.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- A valid project_id for testing (set TEST_PROJECT_ID env var)
- A valid FHI-aims structure_id (set TEST_AIMS_STRUCTURE_ID env var)

To run: uv run pytest tests/integration/test_vibes.py -v -m integration
"""

import os
import pytest
import time
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.simulation.vibes import (
    create_vibes_simulation,
    get_vibes_simulation,
    delete_vibes_simulation,
    associate_file_with_vibes_simulation,
    get_vibes_simulation_files,
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
def test_fhiaims_structure_id():
    """Get test FHI-aims structure ID from environment"""
    structure_id = os.getenv("TEST_AIMS_STRUCTURE_ID")
    if not structure_id:
        pytest.skip("TEST_AIMS_STRUCTURE_ID environment variable must be set")
    return structure_id


@pytest.mark.integration
class TestVibesSimulationIntegration:
    """Integration tests for VIBES simulation operations"""

    def test_create_vibes_simulation_draft(
        self, test_project_id, test_fhiaims_structure_id
    ):
        """Test creating a VIBES simulation in DRAFT mode"""
        simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
            action="DRAFT",
            name="Integration Test VIBES DRAFT",
            description="Test VIBES simulation created by integration tests",
            calculator_parameters={"xc": "pw-lda"},
            calculator_kpoints={"density": 3.5},
            calculator_basis_set={"default": "light"},
        )

        assert "id" in simulation
        assert simulation["name"] == "Integration Test VIBES DRAFT"
        assert (
            simulation["description"]
            == "Test VIBES simulation created by integration tests"
        )
        assert simulation.get("calculator_parameters") == {"xc": "pw-lda"}
        assert simulation.get("calculator_kpoints") == {"density": 3.5}
        assert simulation.get("calculator_basis_set") == {"default": "light"}

        # Clean up
        delete_vibes_simulation(simulation["id"])

    def test_create_vibes_simulation_launch(
        self, test_project_id, test_fhiaims_structure_id
    ):
        """Test creating a VIBES simulation with LAUNCH action"""
        simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
            action="LAUNCH",
            name="Integration Test VIBES LAUNCH",
            description="Test VIBES simulation with launch",
            calculator_parameters={"xc": "pbe"},
            calculator_kpoints={"density": 2.0},
            calculator_basis_set={"default": "tight"},
        )

        assert "id" in simulation
        assert simulation["name"] == "Integration Test VIBES LAUNCH"

        # Verify the simulation was created with a task
        if "task" in simulation and simulation["task"]:
            assert "status" in simulation["task"]

        # Clean up
        delete_vibes_simulation(simulation["id"])

    def test_get_vibes_simulation(self, test_project_id, test_fhiaims_structure_id):
        """Test retrieving a VIBES simulation"""
        # Create a simulation first
        created_simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
            action="DRAFT",
            name="Test Get VIBES Simulation",
        )

        simulation_id = created_simulation["id"]

        # Get the simulation
        retrieved_simulation = get_vibes_simulation(simulation_id)

        assert retrieved_simulation["id"] == simulation_id
        assert retrieved_simulation["name"] == "Test Get VIBES Simulation"
        assert "created_at" in retrieved_simulation
        assert "updated_at" in retrieved_simulation

        # Clean up
        delete_vibes_simulation(simulation_id)

    def test_delete_vibes_simulation(self, test_project_id, test_fhiaims_structure_id):
        """Test deleting a VIBES simulation"""
        # Create a simulation first
        simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
            action="DRAFT",
            name="Test Delete VIBES Simulation",
        )

        simulation_id = simulation["id"]

        # Delete the simulation
        delete_result = delete_vibes_simulation(simulation_id)

        # Verify deletion (API may return different responses)
        # The exact response depends on the backend implementation
        assert delete_result is not None

    def test_vibes_simulation_minimal_parameters(
        self, test_project_id, test_fhiaims_structure_id
    ):
        """Test creating VIBES simulation with only required parameters"""
        simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
        )

        assert "id" in simulation
        # Default action should be DRAFT
        assert simulation.get("name") is None or simulation.get("name") == ""
        assert (
            simulation.get("description") is None or simulation.get("description") == ""
        )

        # Clean up
        delete_vibes_simulation(simulation["id"])

    def test_vibes_simulation_calculator_parameters(
        self, test_project_id, test_fhiaims_structure_id
    ):
        """Test VIBES simulation with various calculator parameters"""
        calc_params = {
            "xc": "pbe",
            "relativistic": "atomic_zora scalar",
            "spin": "collinear",
        }
        calc_kpoints = {"density": 4.0, "offset": [0.0, 0.0, 0.0]}
        calc_basis = {"default": "tier1", "H": "tier2"}

        simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
            action="DRAFT",
            name="Test Calculator Parameters",
            calculator_parameters=calc_params,
            calculator_kpoints=calc_kpoints,
            calculator_basis_set=calc_basis,
        )

        assert "id" in simulation
        assert simulation.get("calculator_parameters") == calc_params
        assert simulation.get("calculator_kpoints") == calc_kpoints
        assert simulation.get("calculator_basis_set") == calc_basis

        # Clean up
        delete_vibes_simulation(simulation["id"])

    def test_vibes_simulation_workflow(
        self, test_project_id, test_fhiaims_structure_id
    ):
        """Test complete VIBES simulation workflow"""
        # Step 1: Create simulation
        simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
            action="DRAFT",
            name="Workflow Test Simulation",
            description="Testing complete workflow",
        )

        simulation_id = simulation["id"]
        assert "id" in simulation

        # Step 2: Retrieve simulation details
        retrieved_simulation = get_vibes_simulation(simulation_id)
        assert retrieved_simulation["id"] == simulation_id
        assert retrieved_simulation["name"] == "Workflow Test Simulation"

        # Step 3: Get simulation files (should be empty initially)
        files = get_vibes_simulation_files(simulation_id)
        assert isinstance(files, (list, dict))
        if isinstance(files, dict) and "results" in files:
            assert len(files["results"]) == 0
        elif isinstance(files, list):
            assert len(files) == 0

        # Step 4: Clean up
        delete_vibes_simulation(simulation_id)


@pytest.mark.integration
class TestVibesSimulationValidation:
    """Integration tests for VIBES simulation validation"""

    def test_invalid_action_parameter(self, test_project_id, test_fhiaims_structure_id):
        """Test that invalid action parameter raises ValueError"""
        with pytest.raises(ValueError, match="Action must be 'DRAFT' or 'LAUNCH'"):
            create_vibes_simulation(
                project_id=test_project_id,
                starting_structure_id=test_fhiaims_structure_id,
                action="INVALID_ACTION",
            )

    def test_invalid_project_id(self, test_fhiaims_structure_id):
        """Test behavior with invalid project ID"""
        # This should raise an API error from the backend
        with pytest.raises(
            Exception
        ):  # Specific exception depends on backend implementation
            create_vibes_simulation(
                project_id="invalid-project-id",
                starting_structure_id=test_fhiaims_structure_id,
                action="DRAFT",
            )

    def test_invalid_structure_id(self, test_project_id):
        """Test behavior with invalid structure ID"""
        # This should raise an API error from the backend
        with pytest.raises(
            Exception
        ):  # Specific exception depends on backend implementation
            create_vibes_simulation(
                project_id=test_project_id,
                starting_structure_id="invalid-structure-id",
                action="DRAFT",
            )


@pytest.mark.integration
class TestVibesSimulationFileOperations:
    """Integration tests for VIBES simulation file operations"""

    def test_get_vibes_simulation_files_empty(
        self, test_project_id, test_fhiaims_structure_id
    ):
        """Test getting files from a simulation with no associated files"""
        simulation = create_vibes_simulation(
            project_id=test_project_id,
            starting_structure_id=test_fhiaims_structure_id,
            action="DRAFT",
            name="File Test Simulation",
        )

        simulation_id = simulation["id"]

        # Get files (should be empty)
        files = get_vibes_simulation_files(simulation_id)

        # Handle different possible response formats
        if isinstance(files, dict):
            if "results" in files:
                assert len(files["results"]) == 0
            elif "count" in files:
                assert files["count"] == 0
        elif isinstance(files, list):
            assert len(files) == 0

        # Clean up
        delete_vibes_simulation(simulation_id)
