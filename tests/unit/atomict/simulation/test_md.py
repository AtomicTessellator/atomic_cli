from unittest.mock import Mock, patch

import pytest

from atomict.simulation.md import (
    associate_user_upload_with_md_simulation,
    create_md_simulation,
    create_md_simulation_file,
    delete_md_simulation,
    delete_md_simulation_file,
    get_md_simulation,
    get_md_simulation_file,
)


@pytest.fixture
def mock_get():
    with patch("atomict.simulation.md.get") as mock:
        yield mock


@pytest.fixture
def mock_post():
    with patch("atomict.simulation.md.post") as mock:
        yield mock


@pytest.fixture
def mock_delete():
    with patch("atomict.simulation.md.delete") as mock:
        yield mock


class TestGetMDSimulation:
    def test_get_md_simulation(self, mock_get):
        """Test getting a MD simulation by ID"""
        mock_get.return_value = {
            "id": "md_123",
            "name": "test_md",
            "mode": 0,
            "temperature": 300.0,
            "task": {"status": 0}
        }

        result = get_md_simulation("md_123")

        assert result == {
            "id": "md_123",
            "name": "test_md",
            "mode": 0,
            "temperature": 300.0,
            "task": {"status": 0}
        }
        mock_get.assert_called_once_with("api/md/md_123/")


class TestGetMDSimulationFile:
    def test_get_md_simulation_file(self, mock_get):
        """Test getting a MD simulation file by ID"""
        mock_get.return_value = {
            "id": "file_456",
            "md_id": "md_123",
            "user_upload_id": "upload_789"
        }

        result = get_md_simulation_file("file_456")

        assert result == {
            "id": "file_456",
            "md_id": "md_123",
            "user_upload_id": "upload_789"
        }
        mock_get.assert_called_once_with("api/md-file/file_456/")


class TestCreateMDSimulation:
    def test_create_basic_md_simulation(self, mock_post):
        """Test creating a basic MD simulation with minimal parameters"""
        mock_post.return_value = {"id": "new_md_123", "mode": 0, "temperature": 300.0}
        
        md_data = {
            "project": "proj_456",
            "source_geometry_id": "struct_789",
            "name": "Basic MD Test",
            "mode": 0,
            "temperature": 300.0,
            "num_steps": 10000,
            "timestep_fs": 1.0,
            "calculator": 2,
            "ensemble_type": 0,
            "thermostat_type": 0,
            "barostat_type": 0,
            "barostat_time_fs": 100.0,
            "thermostat_time_fs": 100.0,
        }

        result = create_md_simulation(md_data)

        assert result == {"id": "new_md_123", "mode": 0, "temperature": 300.0}
        mock_post.assert_called_once_with("api/md/", md_data)

    def test_create_annealing_md_simulation(self, mock_post):
        """Test creating an annealing MD simulation"""
        mock_post.return_value = {
            "id": "new_md_456", 
            "mode": 1, 
            "high_temperature": 800.0,
            "low_temperature": 300.0
        }
        
        md_data = {
            "project": "proj_456",
            "source_geometry_id": "struct_789",
            "name": "Annealing MD Test",
            "mode": 1,
            "high_temperature": 800.0,
            "low_temperature": 300.0,
            "num_intervals": 5,
            "high_temperature_steps": 2000,
            "annealing_steps": 8000,
            "timestep_fs": 1.0,
            "calculator": 2,
            "ensemble_type": 0,
            "thermostat_type": 1,
            "barostat_type": 0,
        }

        result = create_md_simulation(md_data)

        assert result == {
            "id": "new_md_456", 
            "mode": 1, 
            "high_temperature": 800.0,
            "low_temperature": 300.0
        }
        mock_post.assert_called_once_with("api/md/", md_data)

    def test_create_nvt_md_simulation(self, mock_post):
        """Test creating an NVT ensemble MD simulation"""
        mock_post.return_value = {
            "id": "new_md_789", 
            "ensemble_type": 1,
            "thermostat_type": 3
        }
        
        md_data = {
            "project": "proj_456",
            "source_geometry_id": "struct_789",
            "name": "NVT MD Test",
            "mode": 0,
            "temperature": 500.0,
            "num_steps": 15000,
            "timestep_fs": 0.5,
            "calculator": 2,
            "ensemble_type": 1,  # NVT
            "thermostat_type": 3,  # Nosé-Hoover
            "thermostat_time_fs": 25.0,
            "velocity_init": 0,
            "random_seed": 42,
        }

        result = create_md_simulation(md_data)

        assert result == {
            "id": "new_md_789", 
            "ensemble_type": 1,
            "thermostat_type": 3
        }
        mock_post.assert_called_once_with("api/md/", md_data)

    def test_create_md_simulation_with_advanced_parameters(self, mock_post):
        """Test creating MD simulation with advanced parameter options"""
        mock_post.return_value = {"id": "new_md_advanced"}
        
        md_data = {
            "project": "proj_456",
            "source_geometry_id": "struct_789",
            "name": "Advanced MD Test",
            "description": "MD with advanced parameters",
            "mode": 0,
            "temperature": 400.0,
            "num_steps": 20000,
            "timestep_fs": 0.5,
            "calculator": 3,  # eSEN 30M OAM
            "ensemble_type": 0,
            "thermostat_type": 2,  # Bussi-Donadio-Parrinello
            "barostat_type": 1,    # Parrinello-Rahman
            "barostat_time_fs": 150.0,
            "thermostat_time_fs": 75.0,
            "friction": 0.005,
            "cell_handling": 1,
            "auto_convert_cell": True,
            "velocity_init": 0,
            "random_seed": 12345,
            "velocity_scaling_factor": 1.1,
            "dump_interval": 25,
            "output_format": 1,  # Extended XYZ
            "save_stress": True,
            "fixed_atoms": "1,2,5-10",
            "region_based_thermostat": True,
            "vacuum_padding": 5.0,
            "slab_dipole_correction": True,
            "charge": 0,
            "spin_multiplicity": 1,
        }

        result = create_md_simulation(md_data)

        assert result == {"id": "new_md_advanced"}
        mock_post.assert_called_once_with("api/md/", md_data)


class TestCreateMDSimulationFile:
    def test_create_md_simulation_file(self, mock_post):
        """Test creating a MD simulation file"""
        mock_post.return_value = {"id": "new_file_123"}
        
        file_data = {
            "md_id": "md_456",
            "user_upload_id": "upload_789",
            "description": "Output trajectory file"
        }

        result = create_md_simulation_file(file_data)

        assert result == {"id": "new_file_123"}
        mock_post.assert_called_once_with("api/md-file/", file_data)


class TestAssociateUserUploadWithMDSimulation:
    def test_associate_user_upload_with_md_simulation(self, mock_post):
        """Test associating a user upload with an MD simulation"""
        mock_post.return_value = {"id": "association_123"}

        result = associate_user_upload_with_md_simulation(
            user_upload_id="upload_456",
            simulation_id="md_789"
        )

        assert result == {"id": "association_123"}
        mock_post.assert_called_once_with(
            "api/md-file/",
            payload={"user_upload_id": "upload_456", "md_id": "md_789"}
        )


class TestDeleteMDSimulation:
    def test_delete_md_simulation(self, mock_delete):
        """Test deleting a MD simulation"""
        mock_delete.return_value = None

        result = delete_md_simulation("md_123")

        assert result is None
        mock_delete.assert_called_once_with("api/md/md_123/")


class TestDeleteMDSimulationFile:
    def test_delete_md_simulation_file(self, mock_delete):
        """Test deleting a MD simulation file"""
        mock_delete.return_value = None

        result = delete_md_simulation_file("file_456")

        assert result is None
        mock_delete.assert_called_once_with("api/md-file/file_456/")


class TestMDParameterValidation:
    """Test parameter validation scenarios (mocked behavior)"""

    def test_thermostat_type_options(self, mock_post):
        """Test different thermostat type configurations"""
        thermostat_configs = [
            {"thermostat_type": 0, "name": "Langevin"},
            {"thermostat_type": 1, "name": "Berendsen"},
            {"thermostat_type": 2, "name": "Bussi-Donadio-Parrinello"},
            {"thermostat_type": 3, "name": "Nosé-Hoover"},
        ]

        for i, config in enumerate(thermostat_configs):
            mock_post.return_value = {"id": f"md_{i}", "thermostat_type": config["thermostat_type"]}
            
            md_data = {
                "project": "proj_123",
                "source_geometry_id": "struct_456",
                "name": f"Thermostat Test {config['name']}",
                "mode": 0,
                "temperature": 300.0,
                "num_steps": 1000,
                "thermostat_type": config["thermostat_type"],
                "ensemble_type": 0,
                "barostat_type": 0,
            }

            result = create_md_simulation(md_data)
            
            assert result["thermostat_type"] == config["thermostat_type"]
            mock_post.assert_called_with("api/md/", md_data)

    def test_ensemble_type_options(self, mock_post):
        """Test NPT vs NVT ensemble configurations"""
        ensemble_configs = [
            {"ensemble_type": 0, "name": "NPT"},
            {"ensemble_type": 1, "name": "NVT"},
        ]

        for i, config in enumerate(ensemble_configs):
            mock_post.return_value = {"id": f"md_{i}", "ensemble_type": config["ensemble_type"]}
            
            md_data = {
                "project": "proj_123",
                "source_geometry_id": "struct_456",
                "name": f"Ensemble Test {config['name']}",
                "mode": 0,
                "temperature": 300.0,
                "num_steps": 1000,
                "ensemble_type": config["ensemble_type"],
                "thermostat_type": 0,
            }
            
            # Only add barostat for NPT
            if config["ensemble_type"] == 0:
                md_data["barostat_type"] = 0
                md_data["barostat_time_fs"] = 100.0

            result = create_md_simulation(md_data)
            
            assert result["ensemble_type"] == config["ensemble_type"]
            mock_post.assert_called_with("api/md/", md_data)

    def test_output_format_options(self, mock_post):
        """Test different output format configurations"""
        output_configs = [
            {"output_format": 0, "name": "XYZ"},
            {"output_format": 1, "name": "Extended XYZ"},
            {"output_format": 2, "name": "NetCDF"},
        ]

        for i, config in enumerate(output_configs):
            mock_post.return_value = {"id": f"md_{i}", "output_format": config["output_format"]}
            
            md_data = {
                "project": "proj_123",
                "source_geometry_id": "struct_456",
                "name": f"Output Test {config['name']}",
                "mode": 0,
                "temperature": 300.0,
                "num_steps": 1000,
                "output_format": config["output_format"],
                "dump_interval": 50,
                "ensemble_type": 0,
                "thermostat_type": 0,
                "barostat_type": 0,
            }

            result = create_md_simulation(md_data)
            
            assert result["output_format"] == config["output_format"]
            mock_post.assert_called_with("api/md/", md_data)

    def test_calculator_options(self, mock_post):
        """Test different calculator configurations"""
        calculator_configs = [
            {"calculator": 2, "name": "ORB v3 Conservative OMAT"},
            {"calculator": 3, "name": "eSEN 30M OAM"},
        ]

        for i, config in enumerate(calculator_configs):
            mock_post.return_value = {"id": f"md_{i}", "calculator": config["calculator"]}
            
            md_data = {
                "project": "proj_123",
                "source_geometry_id": "struct_456",
                "name": f"Calculator Test {config['name']}",
                "mode": 0,
                "temperature": 300.0,
                "num_steps": 1000,
                "calculator": config["calculator"],
                "ensemble_type": 0,
                "thermostat_type": 0,
                "barostat_type": 0,
            }

            result = create_md_simulation(md_data)
            
            assert result["calculator"] == config["calculator"]
            mock_post.assert_called_with("api/md/", md_data)
