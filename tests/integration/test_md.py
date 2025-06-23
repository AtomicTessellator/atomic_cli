"""
Integration tests for ab-initio MD simulation functionality.
Tests all core MD functions against the live API.
"""

import os
import pytest
from atomict.auth import resolve_token
from atomict.simulation.md import (
    create_md_simulation,
    get_md_simulation,
    delete_md_simulation,
    associate_user_upload_with_md_simulation,
)
from atomict.project import list_projects
from atomict.user.files import get_user_uploads


@pytest.fixture(scope="module")
def auth_setup():
    """Set up authentication for all tests."""
    token = resolve_token()
    os.environ["AT_TOKEN"] = token
    return token


@pytest.fixture(scope="module")
def test_project_id(auth_setup):
    """Get a test project ID."""
    projects = list_projects()
    assert projects and projects.get("results"), "No projects found"
    return projects["results"][0]["id"]


@pytest.fixture(scope="module")
def test_structure_id(auth_setup):
    """Get a test structure ID."""
    structures = get_user_uploads()
    structure_list = (
        structures if isinstance(structures, list) else structures.get("results", [])
    )
    assert structure_list, "No structures found for testing"
    return structure_list[0]["id"]


class TestMDSimulationBasics:
    """Test basic MD simulation operations."""

    def test_create_basic_md_simulation(self, test_project_id, test_structure_id):
        """Test creating a basic MD simulation."""
        md_data = {
            "project": test_project_id,
            "source_geometry_id": test_structure_id,
            "name": "Test Basic MD",
            "description": "Basic integration test",
            "mode": 0,  # Single temperature
            "temperature": 300.0,
            "num_steps": 1000,
            "timestep_fs": 1.0,
            "calculator": 2,  # ORB v3 Conservative
            "ensemble_type": 0,  # NPT
            "thermostat_type": 0,  # Langevin
            "barostat_type": 0,  # Berendsen
            "barostat_time_fs": 100.0,
            "thermostat_time_fs": 100.0,
        }

        result = create_md_simulation(md_data)

        # Verify creation
        assert result["id"]
        assert result["name"] == "Test Basic MD"
        assert result["mode"] == 0
        assert result["temperature"] == 300.0
        assert result["num_steps"] == 1000
        assert result["calculator"] == 2
        assert result["ensemble_type"] == 0
        assert result["task"]["status"] == 0  # DRAFT

        # Clean up
        delete_md_simulation(result["id"])

    def test_create_annealing_md_simulation(self, test_project_id, test_structure_id):
        """Test creating an annealing MD simulation."""
        md_data = {
            "project": test_project_id,
            "source_geometry_id": test_structure_id,
            "name": "Test Annealing MD",
            "description": "Annealing integration test",
            "mode": 1,  # Annealing
            "high_temperature": 800.0,
            "low_temperature": 300.0,
            "num_intervals": 3,
            "high_temperature_steps": 500,
            "annealing_steps": 1500,
            "timestep_fs": 1.0,
            "calculator": 2,
            "ensemble_type": 0,  # NPT
            "thermostat_type": 1,  # Berendsen
            "barostat_type": 0,
            "barostat_time_fs": 100.0,
            "thermostat_time_fs": 50.0,
        }

        result = create_md_simulation(md_data)

        # Verify annealing parameters
        assert result["id"]
        assert result["mode"] == 1
        assert result["high_temperature"] == 800.0
        assert result["low_temperature"] == 300.0
        assert result["num_intervals"] == 3
        assert result["high_temperature_steps"] == 500
        assert result["annealing_steps"] == 1500

        # Clean up
        delete_md_simulation(result["id"])

    def test_nvt_ensemble_simulation(self, test_project_id, test_structure_id):
        """Test creating NVT ensemble simulation."""
        md_data = {
            "project": test_project_id,
            "source_geometry_id": test_structure_id,
            "name": "Test NVT MD",
            "description": "NVT ensemble test",
            "mode": 0,
            "temperature": 500.0,
            "num_steps": 2000,
            "timestep_fs": 0.5,
            "calculator": 2,
            "ensemble_type": 1,  # NVT
            "thermostat_type": 3,  # Nosé-Hoover
            "thermostat_time_fs": 25.0,
            "velocity_init": 0,  # Maxwell-Boltzmann
            "random_seed": 42,
        }

        result = create_md_simulation(md_data)

        # Verify NVT parameters
        assert result["ensemble_type"] == 1
        assert result["thermostat_type"] == 3
        assert result["velocity_init"] == 0
        assert result["random_seed"] == 42
        assert result["timestep_fs"] == 0.5

        # Clean up
        delete_md_simulation(result["id"])


class TestMDSimulationManagement:
    """Test MD simulation management operations."""

    def test_retrieve_md_simulation(self, test_project_id, test_structure_id):
        """Test retrieving MD simulation details."""
        # Create simulation
        md_data = {
            "project": test_project_id,
            "source_geometry_id": test_structure_id,
            "name": "Test Retrieve MD",
            "description": "Retrieval test",
            "mode": 0,
            "temperature": 350.0,
            "num_steps": 1500,
            "timestep_fs": 1.0,
            "calculator": 2,
            "ensemble_type": 0,
            "thermostat_type": 0,
            "barostat_type": 0,
            "barostat_time_fs": 100.0,
            "thermostat_time_fs": 100.0,
        }

        created = create_md_simulation(md_data)
        simulation_id = created["id"]

        # Retrieve simulation
        retrieved = get_md_simulation(simulation_id)

        # Verify retrieval
        assert retrieved["id"] == simulation_id
        assert retrieved["name"] == "Test Retrieve MD"
        assert retrieved["temperature"] == 350.0
        assert retrieved["num_steps"] == 1500
        assert "task" in retrieved
        assert "source_geometry" in retrieved

        # Clean up
        delete_md_simulation(simulation_id)

    def test_delete_md_simulation(self, test_project_id, test_structure_id):
        """Test deleting MD simulation."""
        # Create simulation
        md_data = {
            "project": test_project_id,
            "source_geometry_id": test_structure_id,
            "name": "Test Delete MD",
            "mode": 0,
            "temperature": 300.0,
            "num_steps": 1000,
            "timestep_fs": 1.0,
            "calculator": 2,
            "ensemble_type": 0,
            "thermostat_type": 0,
            "barostat_type": 0,
            "barostat_time_fs": 100.0,
            "thermostat_time_fs": 100.0,
        }

        created = create_md_simulation(md_data)
        simulation_id = created["id"]

        # Delete simulation
        delete_result = delete_md_simulation(simulation_id)

        # Verify deletion (function should not raise exception)
        assert delete_result is not None  # May return empty response


class TestMDParameterVariations:
    """Test various MD parameter combinations."""

    def test_advanced_thermostat_options(self, test_project_id, test_structure_id):
        """Test different thermostat configurations."""
        thermostat_configs = [
            {"thermostat_type": 0, "friction": 0.005},  # Langevin
            {"thermostat_type": 1, "thermostat_time_fs": 75.0},  # Berendsen
            {"thermostat_type": 2, "thermostat_time_fs": 50.0},  # Bussi
            {"thermostat_type": 3, "thermostat_time_fs": 25.0},  # Nosé-Hoover
        ]

        simulation_ids = []

        for i, config in enumerate(thermostat_configs):
            md_data = {
                "project": test_project_id,
                "source_geometry_id": test_structure_id,
                "name": f"Test Thermostat {i}",
                "mode": 0,
                "temperature": 300.0,
                "num_steps": 500,
                "timestep_fs": 1.0,
                "calculator": 2,
                "ensemble_type": 0,
                "barostat_type": 0,
                "barostat_time_fs": 100.0,
                **config,
            }

            result = create_md_simulation(md_data)
            simulation_ids.append(result["id"])

            # Verify thermostat configuration
            assert result["thermostat_type"] == config["thermostat_type"]
            if "friction" in config:
                assert result["friction"] == config["friction"]
            if "thermostat_time_fs" in config:
                assert result["thermostat_time_fs"] == config["thermostat_time_fs"]

        # Clean up all simulations
        for sim_id in simulation_ids:
            delete_md_simulation(sim_id)

    def test_output_format_options(self, test_project_id, test_structure_id):
        """Test different output format configurations."""
        output_configs = [
            {"output_format": 0, "dump_interval": 50},  # XYZ
            {
                "output_format": 1,
                "dump_interval": 25,
                "save_stress": True,
            },  # Extended XYZ
            {"output_format": 2, "dump_interval": 100},  # NetCDF
        ]

        simulation_ids = []

        for i, config in enumerate(output_configs):
            md_data = {
                "project": test_project_id,
                "source_geometry_id": test_structure_id,
                "name": f"Test Output {i}",
                "mode": 0,
                "temperature": 300.0,
                "num_steps": 1000,
                "timestep_fs": 1.0,
                "calculator": 2,
                "ensemble_type": 0,
                "thermostat_type": 0,
                "barostat_type": 0,
                "barostat_time_fs": 100.0,
                "thermostat_time_fs": 100.0,
                **config,
            }

            result = create_md_simulation(md_data)
            simulation_ids.append(result["id"])

            # Verify output configuration
            assert result["output_format"] == config["output_format"]
            assert result["dump_interval"] == config["dump_interval"]
            if "save_stress" in config:
                assert result["save_stress"] == config["save_stress"]

        # Clean up all simulations
        for sim_id in simulation_ids:
            delete_md_simulation(sim_id)

    def test_velocity_initialization_options(self, test_project_id, test_structure_id):
        """Test velocity initialization configurations."""
        md_data = {
            "project": test_project_id,
            "source_geometry_id": test_structure_id,
            "name": "Test Velocity Init",
            "mode": 0,
            "temperature": 400.0,
            "num_steps": 1000,
            "timestep_fs": 1.0,
            "calculator": 2,
            "ensemble_type": 1,  # NVT
            "thermostat_type": 0,
            "thermostat_time_fs": 100.0,
            "velocity_init": 0,  # Maxwell-Boltzmann
            "random_seed": 12345,
            "velocity_scaling_factor": 1.2,
        }

        result = create_md_simulation(md_data)

        # Verify velocity parameters
        assert result["velocity_init"] == 0
        assert result["random_seed"] == 12345
        assert result["velocity_scaling_factor"] == 1.2

        # Clean up
        delete_md_simulation(result["id"])


class TestMDErrorHandling:
    """Test MD simulation error handling."""

    def test_missing_required_fields(self, test_project_id):
        """Test error handling for missing required fields."""
        with pytest.raises(Exception):  # Should raise APIValidationError or similar
            create_md_simulation(
                {
                    "project": test_project_id,
                    # Missing source_geometry_id - should fail
                    "name": "Invalid MD",
                    "mode": 0,
                    "temperature": 300.0,
                }
            )

    def test_invalid_parameter_values(self, test_project_id, test_structure_id):
        """Test error handling for invalid parameter values."""
        with pytest.raises(Exception):  # Should raise validation error
            create_md_simulation(
                {
                    "project": test_project_id,
                    "source_geometry_id": test_structure_id,
                    "name": "Invalid Parameters",
                    "mode": 0,
                    "temperature": -100.0,  # Invalid negative temperature
                    "num_steps": 1000,
                    "calculator": 999,  # Invalid calculator
                }
            )
