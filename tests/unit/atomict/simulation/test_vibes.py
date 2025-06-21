from unittest.mock import Mock, patch

import pytest

from atomict.simulation.vibes import (
    associate_file_with_vibes_simulation,
    create_vibes_simulation,
    delete_vibes_simulation,
    get_vibes_simulation,
    get_vibes_simulation_files,
)


@pytest.fixture
def mock_get():
    with patch("atomict.simulation.vibes.get") as mock:
        yield mock


@pytest.fixture
def mock_post():
    with patch("atomict.simulation.vibes.post") as mock:
        yield mock


@pytest.fixture
def mock_delete():
    with patch("atomict.simulation.vibes.delete") as mock:
        yield mock


class TestCreateVibesSimulation:
    def test_create_vibes_simulation_minimal_params(self, mock_post):
        """Test creating a VIBES simulation with minimal required parameters"""
        mock_post.return_value = {"id": "vibes_123"}

        result = create_vibes_simulation(
            project_id="proj_456", starting_structure_id="struct_789", action="DRAFT"
        )

        assert result == {"id": "vibes_123"}
        mock_post.assert_called_once_with(
            "api/vibes-simulation/",
            {
                "project_id": "proj_456",
                "starting_structure_id": "struct_789",
                "action": "DRAFT",
            },
            extra_headers={"Content-Type": "application/json"},
        )

    def test_create_vibes_simulation_all_params(self, mock_post):
        """Test creating a VIBES simulation with all parameters"""
        mock_post.return_value = {"id": "vibes_456"}

        calculator_params = {"xc": "pw-lda"}
        calculator_kpoints = {"density": 3.5}
        calculator_basis = {"default": "light"}
        extra_kwargs = {"selected_cluster": "cluster_123"}

        result = create_vibes_simulation(
            project_id="proj_123",
            starting_structure_id="struct_456",
            action="LAUNCH",
            name="Test VIBES",
            description="Test description",
            calculator_parameters=calculator_params,
            calculator_kpoints=calculator_kpoints,
            calculator_basis_set=calculator_basis,
            extra_kwargs=extra_kwargs,
        )

        assert result == {"id": "vibes_456"}
        mock_post.assert_called_once_with(
            "api/vibes-simulation/",
            {
                "project_id": "proj_123",
                "starting_structure_id": "struct_456",
                "action": "LAUNCH",
                "name": "Test VIBES",
                "description": "Test description",
                "calculator_parameters": calculator_params,
                "calculator_kpoints": calculator_kpoints,
                "calculator_basis_set": calculator_basis,
                "selected_cluster": "cluster_123",
            },
            extra_headers={"Content-Type": "application/json"},
        )

    def test_create_vibes_simulation_valid_actions(self, mock_post):
        """Test creating VIBES simulations with valid actions"""
        valid_actions = ["DRAFT", "LAUNCH"]

        for action in valid_actions:
            mock_post.return_value = {"id": f"vibes_{action.lower()}"}

            result = create_vibes_simulation(
                project_id="proj_123", starting_structure_id="struct_456", action=action
            )

            assert result == {"id": f"vibes_{action.lower()}"}

    def test_create_vibes_simulation_invalid_action(self, mock_post):
        """Test error handling for invalid action"""
        with pytest.raises(ValueError, match="Action must be 'DRAFT' or 'LAUNCH'"):
            create_vibes_simulation(
                project_id="proj_123",
                starting_structure_id="struct_456",
                action="INVALID",
            )

        mock_post.assert_not_called()

    def test_create_vibes_simulation_default_action(self, mock_post):
        """Test creating VIBES simulation with default action"""
        mock_post.return_value = {"id": "vibes_default"}

        result = create_vibes_simulation(
            project_id="proj_123", starting_structure_id="struct_456"
        )

        assert result == {"id": "vibes_default"}
        expected_payload = mock_post.call_args[0][1]
        assert expected_payload["action"] == "DRAFT"


class TestGetVibesSimulation:
    def test_get_vibes_simulation_without_params(self, mock_get):
        """Test getting a VIBES simulation without additional parameters"""
        mock_get.return_value = {"id": "vibes_123", "name": "test_vibes"}

        result = get_vibes_simulation("vibes_123")

        assert result == {"id": "vibes_123", "name": "test_vibes"}
        mock_get.assert_called_once_with("api/vibes-simulation/vibes_123/")

    def test_get_vibes_simulation_with_params(self, mock_get):
        """Test getting a VIBES simulation with additional parameters"""
        mock_get.return_value = {"id": "vibes_123", "name": "test_vibes"}

        result = get_vibes_simulation("vibes_123", expand=True, include_files=True)

        assert result == {"id": "vibes_123", "name": "test_vibes"}
        mock_get.assert_called_once_with(
            "api/vibes-simulation/vibes_123/?expand=True&include_files=True"
        )

    def test_get_vibes_simulation_with_single_param(self, mock_get):
        """Test getting a VIBES simulation with a single parameter"""
        mock_get.return_value = {"id": "vibes_456"}

        result = get_vibes_simulation("vibes_456", status="completed")

        assert result == {"id": "vibes_456"}
        mock_get.assert_called_once_with(
            "api/vibes-simulation/vibes_456/?status=completed"
        )


class TestDeleteVibesSimulation:
    def test_delete_vibes_simulation(self, mock_delete):
        """Test deleting a VIBES simulation"""
        mock_delete.return_value = {"status": "deleted"}

        result = delete_vibes_simulation("vibes_123")

        assert result == {"status": "deleted"}
        mock_delete.assert_called_once_with("api/vibes-simulation/vibes_123/")


class TestAssociateFileWithVibesSimulation:
    def test_associate_file_success(self, mock_post):
        """Test successful association of file with VIBES simulation"""
        mock_post.return_value = {"id": "association_123"}

        result = associate_file_with_vibes_simulation("upload_456", "vibes_789")

        assert result == {"id": "association_123"}
        mock_post.assert_called_once_with(
            "api/vibes-simulation-file/",
            payload={"user_upload": "upload_456", "simulation": "vibes_789"},
        )


class TestGetVibesSimulationFiles:
    def test_get_vibes_simulation_files(self, mock_get):
        """Test getting files associated with a VIBES simulation"""
        expected_files = [
            {"id": "file_1", "filename": "vibes_input.yaml"},
            {"id": "file_2", "filename": "phonopy.conf"},
        ]
        mock_get.return_value = expected_files

        result = get_vibes_simulation_files("vibes_123")

        assert result == expected_files
        mock_get.assert_called_once_with(
            "api/vibes-simulation-file/?simulation=vibes_123"
        )

    def test_get_vibes_simulation_files_empty(self, mock_get):
        """Test getting files when no files are associated"""
        mock_get.return_value = []

        result = get_vibes_simulation_files("vibes_456")

        assert result == []
        mock_get.assert_called_once_with(
            "api/vibes-simulation-file/?simulation=vibes_456"
        )


class TestVibesModule:
    """Integration-style tests for the VIBES module"""

    @patch("atomict.simulation.vibes.get")
    @patch("atomict.simulation.vibes.post")
    @patch("atomict.simulation.vibes.delete")
    def test_vibes_workflow_integration(self, mock_delete, mock_post, mock_get):
        """Test a typical workflow using multiple VIBES functions"""
        # Mock responses for workflow
        mock_post.side_effect = [
            {"id": "vibes_123"},  # create_vibes_simulation
            {"id": "file_assoc_456"},  # associate_file_with_vibes_simulation
        ]
        mock_get.side_effect = [
            {"id": "vibes_123", "status": "completed"},  # get_vibes_simulation
            [
                {"id": "file_1", "filename": "vibes_output.yaml"}
            ],  # get_vibes_simulation_files
        ]
        mock_delete.return_value = {"status": "deleted"}

        # Create a VIBES simulation
        sim_result = create_vibes_simulation(
            project_id="proj_123",
            starting_structure_id="struct_456",
            action="LAUNCH",
            name="Integration Test VIBES",
        )
        assert sim_result["id"] == "vibes_123"

        # Get the created simulation
        sim_data = get_vibes_simulation("vibes_123")
        assert sim_data["status"] == "completed"

        # Associate a file with the simulation
        file_assoc = associate_file_with_vibes_simulation("upload_789", "vibes_123")
        assert file_assoc["id"] == "file_assoc_456"

        # Get associated files
        files = get_vibes_simulation_files("vibes_123")
        assert len(files) == 1
        assert files[0]["filename"] == "vibes_output.yaml"

        # Delete the simulation
        delete_result = delete_vibes_simulation("vibes_123")
        assert delete_result["status"] == "deleted"

    def test_calculator_parameters_handling(self, mock_post):
        """Test proper handling of calculator parameter dictionaries"""
        mock_post.return_value = {"id": "vibes_calc_test"}

        # Test with comprehensive calculator settings
        calc_params = {
            "xc": "pw-lda",
            "relativistic": "atomic_zora scalar",
            "spin": "collinear",
        }
        calc_kpoints = {"density": 3.5, "offset": [0.0, 0.0, 0.0]}
        calc_basis = {"default": "light", "H": "tier1"}

        result = create_vibes_simulation(
            project_id="proj_calc_test",
            starting_structure_id="struct_calc_test",
            calculator_parameters=calc_params,
            calculator_kpoints=calc_kpoints,
            calculator_basis_set=calc_basis,
        )

        assert result == {"id": "vibes_calc_test"}

        # Verify that calculator parameters are passed correctly
        call_args = mock_post.call_args[0][1]
        assert call_args["calculator_parameters"] == calc_params
        assert call_args["calculator_kpoints"] == calc_kpoints
        assert call_args["calculator_basis_set"] == calc_basis
