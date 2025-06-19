"""Unit tests for K-point convergence operations"""

import unittest
from unittest.mock import patch, MagicMock
import pytest

from atomict.simulation.kpoint import (
    create_kpoint_exploration,
    delete_kpoint_exploration,
    delete_kpoint_simulation,
    create_kpoint_analysis,
    delete_kpoint_analysis,
    STRUCTURE_FIELD_MAP,
)
from atomict.exceptions import APIValidationError


class TestKPointExploration(unittest.TestCase):
    """Test K-point exploration functions"""

    @patch("atomict.simulation.kpoint.post")
    def test_create_kpoint_exploration_success(self, mock_post):
        """Test successful k-point exploration creation"""
        mock_post.return_value = {"id": "test-id", "status": "created"}

        result = create_kpoint_exploration(
            project_id="proj-123",
            name="Test K-point Study",
            structure_id="struct-456",
            structure_type="userupload",
        )

        mock_post.assert_called_once_with(
            "api/kpoint-exploration/",
            payload={
                "project": "proj-123",
                "name": "Test K-point Study",
                "k_point_range_lower": 3,
                "k_point_range_upper": 6,
                "evenly_spaced_kpoints": False,
                "action": "DRAFT",
                "starting_structure_userupload_id": "struct-456",
            },
        )
        assert result == {"id": "test-id", "status": "created"}

    @patch("atomict.simulation.kpoint.post")
    def test_create_kpoint_exploration_with_launch(self, mock_post):
        """Test k-point exploration creation with launch action"""
        mock_post.return_value = {"id": "test-id", "status": "launched"}

        result = create_kpoint_exploration(
            project_id="proj-123",
            name="Test Study",
            structure_id="struct-456",
            structure_type="mlrelax",
            k_point_range_lower=2,
            k_point_range_upper=8,
            evenly_spaced_kpoints=True,
            action="LAUNCH",
            extra_kwargs={"selected_cluster": "cluster-1"},
        )

        expected_payload = {
            "project": "proj-123",
            "name": "Test Study",
            "k_point_range_lower": 2,
            "k_point_range_upper": 8,
            "evenly_spaced_kpoints": True,
            "action": "LAUNCH",
            "starting_structure_mlrelax_id": "struct-456",
            "selected_cluster": "cluster-1",
        }
        mock_post.assert_called_once_with("api/kpoint-exploration/", payload=expected_payload)

    def test_create_kpoint_exploration_invalid_action(self):
        """Test validation for invalid action"""
        with pytest.raises(APIValidationError) as exc_info:
            create_kpoint_exploration(
                project_id="proj-123",
                name="Test Study",
                structure_id="struct-456",
                action="INVALID",
            )
        assert "Invalid action 'INVALID'" in str(exc_info.value)

    def test_create_kpoint_exploration_invalid_structure_type(self):
        """Test validation for invalid structure type"""
        with pytest.raises(APIValidationError) as exc_info:
            create_kpoint_exploration(
                project_id="proj-123",
                name="Test Study",
                structure_id="struct-456",
                structure_type="invalid_type",
            )
        assert "Invalid structure_type 'invalid_type'" in str(exc_info.value)

    def test_create_kpoint_exploration_invalid_k_point_range(self):
        """Test validation for invalid k-point range"""
        with pytest.raises(APIValidationError) as exc_info:
            create_kpoint_exploration(
                project_id="proj-123",
                name="Test Study",
                structure_id="struct-456",
                k_point_range_lower=6,
                k_point_range_upper=3,
            )
        assert "k_point_range_lower must be less than k_point_range_upper" in str(
            exc_info.value
        )

    @patch("atomict.simulation.kpoint.delete")
    def test_delete_kpoint_exploration(self, mock_delete):
        """Test k-point exploration deletion"""
        mock_delete.return_value = {"deleted": True}

        result = delete_kpoint_exploration("exploration-123")

        mock_delete.assert_called_once_with("api/kpoint-exploration/exploration-123/")
        assert result == {"deleted": True}


class TestKPointSimulation(unittest.TestCase):
    """Test K-point simulation functions"""

    @patch("atomict.simulation.kpoint.delete")
    def test_delete_kpoint_simulation(self, mock_delete):
        """Test k-point simulation deletion"""
        mock_delete.return_value = {"deleted": True}

        result = delete_kpoint_simulation("simulation-123")

        mock_delete.assert_called_once_with("api/kpoint-simulation/simulation-123/")
        assert result == {"deleted": True}


class TestKPointAnalysis(unittest.TestCase):
    """Test K-point analysis functions"""

    @patch("atomict.simulation.kpoint.post")
    def test_create_kpoint_analysis_basic(self, mock_post):
        """Test basic k-point analysis creation"""
        mock_post.return_value = {"id": "analysis-123", "status": "created"}

        result = create_kpoint_analysis("exploration-456")

        mock_post.assert_called_once_with(
            "api/kpoint-analysis/", payload={"exploration": "exploration-456"}
        )
        assert result == {"id": "analysis-123", "status": "created"}

    @patch("atomict.simulation.kpoint.post")
    def test_create_kpoint_analysis_with_extras(self, mock_post):
        """Test k-point analysis creation with extra parameters"""
        mock_post.return_value = {"id": "analysis-123", "status": "created"}

        result = create_kpoint_analysis(
            "exploration-456", extra_kwargs={"custom_param": "value"}
        )

        expected_payload = {"exploration": "exploration-456", "custom_param": "value"}
        mock_post.assert_called_once_with("api/kpoint-analysis/", payload=expected_payload)

    @patch("atomict.simulation.kpoint.delete")
    def test_delete_kpoint_analysis(self, mock_delete):
        """Test k-point analysis deletion"""
        mock_delete.return_value = {"deleted": True}

        result = delete_kpoint_analysis("analysis-123")

        mock_delete.assert_called_once_with("api/kpoint-analysis/analysis-123/")
        assert result == {"deleted": True}


class TestStructureFieldMapping(unittest.TestCase):
    """Test structure field mapping constants"""

    def test_structure_field_map_completeness(self):
        """Test that all expected structure types are mapped"""
        expected_types = ["fhiaims", "mlrelax", "userupload"]
        assert set(STRUCTURE_FIELD_MAP.keys()) == set(expected_types)

    def test_structure_field_map_values(self):
        """Test that field mappings are correct"""
        assert STRUCTURE_FIELD_MAP["fhiaims"] == "starting_structure_id"
        assert STRUCTURE_FIELD_MAP["mlrelax"] == "starting_structure_mlrelax_id"
        assert STRUCTURE_FIELD_MAP["userupload"] == "starting_structure_userupload_id"


if __name__ == "__main__":
    unittest.main()
