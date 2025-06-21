from unittest.mock import patch

import pytest

from atomict.simulation.ht_sqs import (
    create_ht_sqs_exploration,
    delete_ht_sqs_exploration,
    get_ht_sqs_child,
    get_ht_sqs_exploration,
    list_ht_sqs_children,
    list_ht_sqs_explorations,
)


class TestHTSQSExploration:
    """Unit tests for HT SQS exploration functions"""

    def test_create_ht_sqs_exploration_success(self):
        """Test successful HT SQS exploration creation"""
        with patch("atomict.simulation.ht_sqs.post") as mock_post:
            mock_post.return_value = {
                "message": "HTSQSExploration created successfully"
            }

            target_concentrations = [
                {"element": "Fe", "weight": 0.75},
                {"element": "Ni", "weight": 0.25},
            ]
            generated_permutations = [{"Fe": 0.8, "Ni": 0.2}, {"Fe": 0.7, "Ni": 0.3}]

            result = create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=target_concentrations,
                generated_permutations=generated_permutations,
                structure_id="test-structure",
            )

            assert result == {"message": "HTSQSExploration created successfully"}
            mock_post.assert_called_once()

            # Verify payload structure
            args, kwargs = mock_post.call_args
            payload = args[1]
            assert payload["project"] == "test-project"
            assert payload["name"] == "Test HT SQS"
            assert payload["target_concentrations"] == target_concentrations
            assert payload["generated_permutations"] == generated_permutations
            assert payload["starting_structure_userupload_id"] == "test-structure"
            assert payload["num_structures_per_permutation"] == 1  # default

    def test_create_ht_sqs_exploration_with_mlrelax_structure(self):
        """Test HT SQS exploration creation with MLRelax structure"""
        with patch("atomict.simulation.ht_sqs.post") as mock_post:
            mock_post.return_value = {
                "message": "HTSQSExploration created successfully"
            }

            result = create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[{"element": "Fe", "weight": 1.0}],
                generated_permutations=[{"Fe": 1.0}],
                structure_id="test-structure",
                structure_type="mlrelax",
            )

            args, kwargs = mock_post.call_args
            payload = args[1]
            assert payload["starting_structure_mlrelax_id"] == "test-structure"
            assert "starting_structure_userupload_id" not in payload

    def test_create_ht_sqs_exploration_invalid_structure_type(self):
        """Test validation of invalid structure type"""
        with pytest.raises(ValueError, match="structure_type must be one of"):
            create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[{"element": "Fe", "weight": 1.0}],
                generated_permutations=[{"Fe": 1.0}],
                structure_id="test-structure",
                structure_type="invalid",
            )

    def test_create_ht_sqs_exploration_invalid_target_concentrations(self):
        """Test validation of target concentrations"""
        # Empty concentrations
        with pytest.raises(ValueError, match="target_concentrations cannot be empty"):
            create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[],
                generated_permutations=[{"Fe": 1.0}],
                structure_id="test-structure",
            )

        # Concentrations don't sum to 1.0
        with pytest.raises(
            ValueError, match="Sum of target concentrations must equal 1.0"
        ):
            create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[
                    {"element": "Fe", "weight": 0.6},
                    {"element": "Ni", "weight": 0.6},
                ],
                generated_permutations=[{"Fe": 0.6, "Ni": 0.4}],
                structure_id="test-structure",
            )

        # Invalid concentration format
        with pytest.raises(
            ValueError, match="Each target concentration must be a dict"
        ):
            create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[{"element": "Fe"}],  # missing weight
                generated_permutations=[{"Fe": 1.0}],
                structure_id="test-structure",
            )

    def test_create_ht_sqs_exploration_invalid_permutations(self):
        """Test validation of generated permutations"""
        # Empty permutations
        with pytest.raises(ValueError, match="generated_permutations cannot be empty"):
            create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[{"element": "Fe", "weight": 1.0}],
                generated_permutations=[],
                structure_id="test-structure",
            )

        # Permutation doesn't sum to 1.0
        with pytest.raises(
            ValueError, match="Permutation 0 concentrations must sum to 1.0"
        ):
            create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[{"element": "Fe", "weight": 1.0}],
                generated_permutations=[{"Fe": 0.8}],  # doesn't sum to 1.0
                structure_id="test-structure",
            )

    def test_create_ht_sqs_exploration_invalid_num_structures(self):
        """Test validation of num_structures_per_permutation"""
        with pytest.raises(
            ValueError, match="num_structures_per_permutation must be >= 1"
        ):
            create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[{"element": "Fe", "weight": 1.0}],
                generated_permutations=[{"Fe": 1.0}],
                structure_id="test-structure",
                num_structures_per_permutation=0,
            )

    def test_get_ht_sqs_exploration_simple(self):
        """Test getting HT SQS exploration without children"""
        with patch("atomict.simulation.ht_sqs.get") as mock_get:
            mock_get.return_value = {"id": "123", "name": "test_ht_sqs"}

            result = get_ht_sqs_exploration("123")

            assert result == {"id": "123", "name": "test_ht_sqs"}
            mock_get.assert_called_once_with("api/ht-sqs-exploration/123/")

    def test_get_ht_sqs_exploration_with_children(self):
        """Test getting HT SQS exploration with children"""
        with patch("atomict.simulation.ht_sqs.get") as mock_get:
            mock_get.return_value = {
                "id": "123",
                "name": "test_ht_sqs",
                "children": [{"id": "child1"}, {"id": "child2"}],
            }

            result = get_ht_sqs_exploration("123", children=True)

            assert result["children"] is not None
            mock_get.assert_called_once_with(
                "api/ht-sqs-exploration/123/?children=true"
            )

    def test_list_ht_sqs_explorations_simple(self):
        """Test listing HT SQS explorations without filters"""
        with patch("atomict.simulation.ht_sqs.get") as mock_get:
            mock_get.return_value = {"results": [{"id": "1"}, {"id": "2"}]}

            result = list_ht_sqs_explorations()

            assert len(result["results"]) == 2
            mock_get.assert_called_once_with("api/ht-sqs-exploration/")

    def test_list_ht_sqs_explorations_with_filters(self):
        """Test listing HT SQS explorations with project filter"""
        with patch("atomict.simulation.ht_sqs.get") as mock_get:
            mock_get.return_value = {"results": [{"id": "1"}]}

            result = list_ht_sqs_explorations(
                project_id="test-project", limit=10, offset=0
            )

            mock_get.assert_called_once_with(
                "api/ht-sqs-exploration/?project__id=test-project&limit=10&offset=0"
            )

    def test_delete_ht_sqs_exploration(self):
        """Test deleting HT SQS exploration"""
        with patch("atomict.simulation.ht_sqs.delete") as mock_delete:
            mock_delete.return_value = {"message": "Deleted successfully"}

            result = delete_ht_sqs_exploration("123")

            assert result == {"message": "Deleted successfully"}
            mock_delete.assert_called_once_with("api/ht-sqs-exploration/123/")

    def test_get_ht_sqs_child(self):
        """Test getting HT SQS child"""
        with patch("atomict.simulation.ht_sqs.get") as mock_get:
            mock_get.return_value = {"id": "123", "sqs_exploration": {"id": "sqs-123"}}

            result = get_ht_sqs_child("123")

            assert result["sqs_exploration"]["id"] == "sqs-123"
            mock_get.assert_called_once_with("api/ht-sqs/123/")

    def test_list_ht_sqs_children(self):
        """Test listing HT SQS children with parent filter"""
        with patch("atomict.simulation.ht_sqs.get") as mock_get:
            mock_get.return_value = {"results": [{"id": "child1"}, {"id": "child2"}]}

            result = list_ht_sqs_children(ht_sqs_exploration_id="parent-123")

            assert len(result["results"]) == 2
            mock_get.assert_called_once_with(
                "api/ht-sqs/?ht_sqs_exploration=parent-123"
            )

    def test_create_with_extra_kwargs(self):
        """Test HT SQS creation with extra parameters"""
        with patch("atomict.simulation.ht_sqs.post") as mock_post:
            mock_post.return_value = {
                "message": "HTSQSExploration created successfully"
            }

            result = create_ht_sqs_exploration(
                project_id="test-project",
                name="Test HT SQS",
                target_concentrations=[{"element": "Fe", "weight": 1.0}],
                generated_permutations=[{"Fe": 1.0}],
                structure_id="test-structure",
                extra_kwargs={"custom_field": "custom_value"},
            )

            args, kwargs = mock_post.call_args
            payload = args[1]
            assert payload["custom_field"] == "custom_value"
