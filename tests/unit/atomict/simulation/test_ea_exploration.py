import pytest
from unittest.mock import patch
from atomict.simulation.ea import (
    create_ea_exploration,
    delete_ea_exploration,
    delete_ea_exploration_sample,
    create_ea_exploration_analysis,
    delete_ea_exploration_analysis,
    create_ea_exploration_analysis_file,
    delete_ea_exploration_analysis_file,
    STRESS_ALGORITHMS,
    STRESS_METHODS,
    CALCULATORS,
)


class TestEAExplorationCreation:
    """Test EA exploration creation function"""

    @patch("atomict.simulation.ea.post")
    def test_create_ea_exploration_basic(self, mock_post):
        """Test basic EA exploration creation"""
        mock_post.return_value = {"id": "test-id", "status": "DRAFT"}

        result = create_ea_exploration(
            project_id="proj-123", name="Test EA Exploration", structure_id="struct-456"
        )

        # Verify API call
        mock_post.assert_called_once_with(
            "api/ea-exploration/",
            {
                "project": "proj-123",
                "name": "Test EA Exploration",
                "description": None,
                "strains_list": [-0.06, -0.03, 0.03, 0.06],
                "stress_algorithm": 2,  # ASESS default
                "stress_method": 1,  # Dynamic default
                "calculator": 0,  # FHI-AIMS default
                "num_last_samples": 1000,
                "make_conventional_cell": False,
                "remove_spurious_distortions": True,
                "add_vacuum": None,
                "is_ht": False,
                "action": "DRAFT",
                "starting_structure_userupload_id": "struct-456",
            },
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {"id": "test-id", "status": "DRAFT"}

    @patch("atomict.simulation.ea.post")
    def test_create_ea_exploration_with_all_params(self, mock_post):
        """Test EA exploration creation with all parameters"""
        mock_post.return_value = {"id": "test-id", "status": "READY"}

        create_ea_exploration(
            project_id="proj-123",
            name="Full EA Test",
            structure_id="struct-456",
            structure_type="mlrelax",
            action="LAUNCH",
            strains_list=[-0.08, -0.04, 0.04, 0.08],
            stress_algorithm=1,  # OHESS
            stress_method=0,  # Static
            calculator=3,  # ORB_V3
            num_last_samples=2000,
            description="Full test exploration",
            make_conventional_cell=True,
            remove_spurious_distortions=False,
            add_vacuum=10,
            is_ht=True,
        )

        # Verify payload structure
        call_args = mock_post.call_args
        payload = call_args[0][1]

        assert payload["starting_structure_mlrelax_id"] == "struct-456"
        assert payload["action"] == "LAUNCH"
        assert payload["strains_list"] == [-0.08, -0.04, 0.04, 0.08]
        assert payload["stress_algorithm"] == 1
        assert payload["stress_method"] == 0
        assert payload["calculator"] == 3
        assert payload["is_ht"] is True

    def test_create_ea_exploration_validation_errors(self):
        """Test validation errors in EA exploration creation"""

        # Test invalid action
        with pytest.raises(ValueError, match="Action must be 'DRAFT' or 'LAUNCH'"):
            create_ea_exploration("proj", "name", "struct", action="INVALID")

        # Test invalid structure_type
        with pytest.raises(ValueError, match="structure_type must be one of"):
            create_ea_exploration("proj", "name", "struct", structure_type="invalid")

        # Test invalid add_vacuum negative
        with pytest.raises(ValueError, match="add_vacuum must be between 0 and 100"):
            create_ea_exploration("proj", "name", "struct", add_vacuum=-5)

        # Test invalid add_vacuum
        with pytest.raises(ValueError, match="add_vacuum must be between 0 and 100"):
            create_ea_exploration("proj", "name", "struct", add_vacuum=150)

        # Test invalid strains_list length
        with pytest.raises(
            ValueError, match="strains_list must contain at least 4 elements"
        ):
            create_ea_exploration("proj", "name", "struct", strains_list=[0.1, 0.2])

    @patch("atomict.simulation.ea.post")
    def test_create_ea_exploration_structure_types(self, mock_post):
        """Test different structure types map to correct API fields"""
        mock_post.return_value = {"id": "test"}

        # Test FHI-AIMS structure
        create_ea_exploration("proj", "name", "struct", structure_type="fhiaims")
        payload = mock_post.call_args[0][1]
        assert "starting_structure_id" in payload
        assert payload["starting_structure_id"] == "struct"

        mock_post.reset_mock()

        # Test MLRelax structure
        create_ea_exploration("proj", "name", "struct", structure_type="mlrelax")
        payload = mock_post.call_args[0][1]
        assert "starting_structure_mlrelax_id" in payload
        assert payload["starting_structure_mlrelax_id"] == "struct"

        mock_post.reset_mock()

        # Test UserUpload structure (default)
        create_ea_exploration("proj", "name", "struct")
        payload = mock_post.call_args[0][1]
        assert "starting_structure_userupload_id" in payload
        assert payload["starting_structure_userupload_id"] == "struct"


class TestEAExplorationDeletion:
    """Test EA exploration deletion functions"""

    @patch("atomict.simulation.ea.delete")
    def test_delete_ea_exploration(self, mock_delete):
        """Test EA exploration deletion"""
        mock_delete.return_value = {"success": True}

        result = delete_ea_exploration("exploration-123")

        mock_delete.assert_called_once_with("api/ea-exploration/exploration-123/")
        assert result == {"success": True}

    @patch("atomict.simulation.ea.delete")
    def test_delete_ea_exploration_sample(self, mock_delete):
        """Test EA exploration sample deletion"""
        mock_delete.return_value = {"success": True}

        result = delete_ea_exploration_sample("sample-456")

        mock_delete.assert_called_once_with("api/ea-exploration-sample/sample-456/")
        assert result == {"success": True}


class TestEAExplorationAnalysis:
    """Test EA exploration analysis functions"""

    @patch("atomict.simulation.ea.post")
    def test_create_ea_exploration_analysis_basic(self, mock_post):
        """Test basic EA exploration analysis creation"""
        mock_post.return_value = {"id": "analysis-123", "status": "DRAFT"}

        result = create_ea_exploration_analysis("exploration-456")

        mock_post.assert_called_once_with(
            "api/ea-exploration-analysis/",
            {
                "exploration_id": "exploration-456",
                "compute_directional_properties": True,
                "action": "DRAFT",
            },
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {"id": "analysis-123", "status": "DRAFT"}

    @patch("atomict.simulation.ea.post")
    def test_create_ea_exploration_analysis_with_params(self, mock_post):
        """Test EA exploration analysis creation with custom parameters"""
        mock_post.return_value = {"id": "analysis-123", "status": "READY"}

        create_ea_exploration_analysis(
            exploration_id="exploration-456",
            compute_directional_properties=False,
            action="LAUNCH",
        )

        payload = mock_post.call_args[0][1]
        assert payload["compute_directional_properties"] is False
        assert payload["action"] == "LAUNCH"

    def test_create_ea_exploration_analysis_validation(self):
        """Test validation in analysis creation"""
        with pytest.raises(ValueError, match="Action must be 'DRAFT' or 'LAUNCH'"):
            create_ea_exploration_analysis("exploration-123", action="INVALID")

    @patch("atomict.simulation.ea.delete")
    def test_delete_ea_exploration_analysis(self, mock_delete):
        """Test EA exploration analysis deletion"""
        mock_delete.return_value = {"success": True}

        result = delete_ea_exploration_analysis("analysis-789")

        mock_delete.assert_called_once_with("api/ea-exploration-analysis/analysis-789/")
        assert result == {"success": True}


class TestEAExplorationAnalysisFiles:
    """Test EA exploration analysis file functions"""

    @patch("atomict.simulation.ea.post")
    def test_create_ea_exploration_analysis_file(self, mock_post):
        """Test EA exploration analysis file creation"""
        mock_post.return_value = {"id": "file-123", "status": "created"}

        result = create_ea_exploration_analysis_file(
            analysis_id="analysis-456", user_upload_id="upload-789"
        )

        mock_post.assert_called_once_with(
            "api/ea-exploration-analysis-file/",
            {"analysis_id": "analysis-456", "user_upload_id": "upload-789"},
            extra_headers={"Content-Type": "application/json"},
        )

        assert result == {"id": "file-123", "status": "created"}

    @patch("atomict.simulation.ea.delete")
    def test_delete_ea_exploration_analysis_file(self, mock_delete):
        """Test EA exploration analysis file deletion"""
        mock_delete.return_value = {"success": True}

        result = delete_ea_exploration_analysis_file("file-123")

        mock_delete.assert_called_once_with(
            "api/ea-exploration-analysis-file/file-123/"
        )
        assert result == {"success": True}


class TestEAExplorationConstants:
    """Test EA exploration constants"""

    def test_stress_algorithms_constant(self):
        """Test stress algorithms constant values"""
        assert STRESS_ALGORITHMS["ULICS"] == 0
        assert STRESS_ALGORITHMS["OHESS"] == 1
        assert STRESS_ALGORITHMS["ASESS"] == 2

    def test_stress_methods_constant(self):
        """Test stress methods constant values"""
        assert STRESS_METHODS["Static"] == 0
        assert STRESS_METHODS["Dynamic"] == 1

    def test_calculators_constant(self):
        """Test calculators constant values"""
        assert CALCULATORS["FHI-AIMS"] == 0
        assert CALCULATORS["ORB_V3_CONSERVATIVE"] == 3
        assert CALCULATORS["ESEN_30M_OAM"] == 4
