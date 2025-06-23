from unittest.mock import patch

import pytest

from atomict.simulation.ht_mlrelax import (
    ML_MODELS,
    create_ht_mlrelax_exploration,
    delete_ht_mlrelax_exploration,
    get_ht_mlrelax_exploration,
    launch_ht_mlrelax_exploration,
    list_ht_mlrelax_explorations,
    update_ht_mlrelax_exploration,
)


class TestHTMLRelaxExploration:
    """Unit tests for HT ML Relaxation exploration functions"""

    def test_create_ht_mlrelax_exploration_success(self):
        """Test successful HT ML Relaxation exploration creation"""
        with patch("atomict.simulation.ht_mlrelax.post") as mock_post:
            mock_post.return_value = {"message": "Created 5 relaxations successfully"}

            result = create_ht_mlrelax_exploration(
                project_id="test-project",
                name="Test HT MLRelax",
                source_ht_sqs_exploration_id="test-sqs-exploration",
                f_max=0.05,
            )

            assert result == {"message": "Created 5 relaxations successfully"}
            mock_post.assert_called_once()

            # Verify payload structure
            args, kwargs = mock_post.call_args
            payload = args[1]
            assert payload["project"] == "test-project"
            assert payload["name"] == "Test HT MLRelax"
            assert payload["source_ht_sqs_exploration"] == "test-sqs-exploration"
            assert payload["f_max"] == 0.05
            assert payload["model"] == ML_MODELS["mattersim_1_0_0_5m"]  # default

    def test_create_ht_mlrelax_exploration_with_custom_model(self):
        """Test HT ML Relaxation exploration creation with custom model"""
        with patch("atomict.simulation.ht_mlrelax.post") as mock_post:
            mock_post.return_value = {"message": "Created 3 relaxations successfully"}

            result = create_ht_mlrelax_exploration(
                project_id="test-project",
                name="Test HT MLRelax",
                source_ht_sqs_exploration_id="test-sqs-exploration",
                model="orb_d3_v2",
            )

            args, kwargs = mock_post.call_args
            payload = args[1]
            assert payload["model"] == ML_MODELS["orb_d3_v2"]

    def test_create_ht_mlrelax_exploration_invalid_model(self):
        """Test validation of invalid model"""
        with pytest.raises(ValueError, match="Invalid model 'invalid_model'"):
            create_ht_mlrelax_exploration(
                project_id="test-project",
                name="Test HT MLRelax",
                source_ht_sqs_exploration_id="test-sqs-exploration",
                model="invalid_model",
            )

    def test_create_ht_mlrelax_exploration_with_all_params(self):
        """Test HT ML Relaxation exploration creation with all parameters"""
        with patch("atomict.simulation.ht_mlrelax.post") as mock_post:
            mock_post.return_value = {"message": "Created 2 relaxations successfully"}

            result = create_ht_mlrelax_exploration(
                project_id="test-project",
                name="Full Test HT MLRelax",
                source_ht_sqs_exploration_id="test-sqs-exploration",
                f_max=0.03,
                model="esen_30m_oam",
                description="Test description",
                extra_kwargs={"custom_field": "custom_value"},
            )

            args, kwargs = mock_post.call_args
            payload = args[1]
            assert payload["project"] == "test-project"
            assert payload["name"] == "Full Test HT MLRelax"
            assert payload["source_ht_sqs_exploration"] == "test-sqs-exploration"
            assert payload["f_max"] == 0.03
            assert payload["model"] == ML_MODELS["esen_30m_oam"]
            assert payload["description"] == "Test description"
            assert payload["custom_field"] == "custom_value"

    def test_get_ht_mlrelax_exploration_simple(self):
        """Test getting HT ML Relaxation exploration without children"""
        with patch("atomict.simulation.ht_mlrelax.get") as mock_get:
            mock_get.return_value = {
                "id": "123",
                "name": "test_ht_mlrelax",
                "f_max": 0.05,
                "model": 1,
            }

            result = get_ht_mlrelax_exploration("123")

            assert result == {
                "id": "123",
                "name": "test_ht_mlrelax",
                "f_max": 0.05,
                "model": 1,
            }
            mock_get.assert_called_once_with("api/ht-mlrelax/123/")

    def test_get_ht_mlrelax_exploration_with_children(self):
        """Test getting HT ML Relaxation exploration with children"""
        with patch("atomict.simulation.ht_mlrelax.get") as mock_get:
            mock_get.return_value = {
                "id": "123",
                "name": "test_ht_mlrelax",
                "children": [
                    {"id": "child1", "mlrelax": {"id": "mlrelax1"}},
                    {"id": "child2", "mlrelax": {"id": "mlrelax2"}},
                ],
            }

            result = get_ht_mlrelax_exploration("123", include_children=True)

            assert result["children"] is not None
            assert len(result["children"]) == 2
            mock_get.assert_called_once_with("api/ht-mlrelax/123/?children=true")

    def test_get_ht_mlrelax_exploration_with_extra_params(self):
        """Test getting HT ML Relaxation exploration with additional parameters"""
        with patch("atomict.simulation.ht_mlrelax.get") as mock_get:
            mock_get.return_value = {"id": "123", "name": "test_ht_mlrelax"}

            result = get_ht_mlrelax_exploration(
                "123", include_children=True, extra_param="value"
            )

            # Check that the URL contains both parameters (order may vary)
            args, kwargs = mock_get.call_args
            url = args[0]
            assert url.startswith("api/ht-mlrelax/123/?")
            assert "children=true" in url
            assert "extra_param=value" in url

    def test_list_ht_mlrelax_explorations_simple(self):
        """Test listing HT ML Relaxation explorations without filters"""
        with patch("atomict.simulation.ht_mlrelax.get") as mock_get:
            mock_get.return_value = {
                "results": [{"id": "1", "name": "test1"}, {"id": "2", "name": "test2"}]
            }

            result = list_ht_mlrelax_explorations()

            assert len(result["results"]) == 2
            mock_get.assert_called_once_with("api/ht-mlrelax/")

    def test_list_ht_mlrelax_explorations_with_filters(self):
        """Test listing HT ML Relaxation explorations with filters"""
        with patch("atomict.simulation.ht_mlrelax.get") as mock_get:
            mock_get.return_value = {"results": [{"id": "1"}]}

            result = list_ht_mlrelax_explorations(
                project_id="test-project", limit=10, offset=0
            )

            mock_get.assert_called_once_with(
                "api/ht-mlrelax/?project__id=test-project&limit=10&offset=0"
            )

    def test_list_ht_mlrelax_explorations_with_extra_params(self):
        """Test listing HT ML Relaxation explorations with extra parameters"""
        with patch("atomict.simulation.ht_mlrelax.get") as mock_get:
            mock_get.return_value = {"results": []}

            result = list_ht_mlrelax_explorations(
                project_id="test-project", custom_filter="value"
            )

            # Check that the URL contains both parameters (order may vary)
            args, kwargs = mock_get.call_args
            url = args[0]
            assert url.startswith("api/ht-mlrelax/?")
            assert "project__id=test-project" in url
            assert "custom_filter=value" in url

    def test_update_ht_mlrelax_exploration_partial(self):
        """Test partial update of HT ML Relaxation exploration"""
        with patch("atomict.simulation.ht_mlrelax.patch") as mock_patch:
            mock_patch.return_value = {
                "id": "123",
                "name": "Updated Name",
                "f_max": 0.03,
            }

            result = update_ht_mlrelax_exploration(
                exploration_id="123", name="Updated Name", f_max=0.03
            )

            assert result["name"] == "Updated Name"
            assert result["f_max"] == 0.03
            mock_patch.assert_called_once_with(
                "api/ht-mlrelax/123/", {"name": "Updated Name", "f_max": 0.03}
            )

    def test_update_ht_mlrelax_exploration_with_model(self):
        """Test updating HT ML Relaxation exploration with model change"""
        with patch("atomict.simulation.ht_mlrelax.patch") as mock_patch:
            mock_patch.return_value = {"id": "123", "model": 2}

            result = update_ht_mlrelax_exploration(
                exploration_id="123", model="orb_v3_conservative"
            )

            args, kwargs = mock_patch.call_args
            payload = args[1]
            assert payload["model"] == ML_MODELS["orb_v3_conservative"]

    def test_update_ht_mlrelax_exploration_invalid_model(self):
        """Test update with invalid model"""
        with pytest.raises(ValueError, match="Invalid model 'invalid_model'"):
            update_ht_mlrelax_exploration(exploration_id="123", model="invalid_model")

    def test_update_ht_mlrelax_exploration_with_extra_kwargs(self):
        """Test update with extra parameters"""
        with patch("atomict.simulation.ht_mlrelax.patch") as mock_patch:
            mock_patch.return_value = {"id": "123"}

            result = update_ht_mlrelax_exploration(
                exploration_id="123",
                name="Test",
                extra_kwargs={"custom_field": "custom_value"},
            )

            args, kwargs = mock_patch.call_args
            payload = args[1]
            assert payload["name"] == "Test"
            assert payload["custom_field"] == "custom_value"

    def test_delete_ht_mlrelax_exploration(self):
        """Test deleting HT ML Relaxation exploration"""
        with patch("atomict.simulation.ht_mlrelax.delete") as mock_delete:
            mock_delete.return_value = {"message": "Deleted successfully"}

            result = delete_ht_mlrelax_exploration("123")

            assert result == {"message": "Deleted successfully"}
            mock_delete.assert_called_once_with("api/ht-mlrelax/123/")

    def test_launch_ht_mlrelax_exploration_simple(self):
        """Test launching HT ML Relaxation exploration without cluster"""
        with patch("atomict.simulation.ht_mlrelax.patch") as mock_patch:
            mock_patch.return_value = {"message": "Launched successfully"}

            result = launch_ht_mlrelax_exploration("123")

            assert result == {"message": "Launched successfully"}
            mock_patch.assert_called_once_with(
                "api/ht-mlrelax/123/", {"action": "LAUNCH"}
            )

    def test_launch_ht_mlrelax_exploration_with_cluster(self):
        """Test launching HT ML Relaxation exploration with cluster"""
        with patch("atomict.simulation.ht_mlrelax.patch") as mock_patch:
            mock_patch.return_value = {"message": "Launched successfully"}

            result = launch_ht_mlrelax_exploration("123", cluster_id="cluster-456")

            args, kwargs = mock_patch.call_args
            payload = args[1]
            assert payload["action"] == "LAUNCH"
            assert payload["selected_cluster"] == "cluster-456"

    def test_launch_ht_mlrelax_exploration_with_extra_kwargs(self):
        """Test launching HT ML Relaxation exploration with extra parameters"""
        with patch("atomict.simulation.ht_mlrelax.patch") as mock_patch:
            mock_patch.return_value = {"message": "Launched successfully"}

            result = launch_ht_mlrelax_exploration(
                "123", extra_kwargs={"billing_org": "test-org"}
            )

            args, kwargs = mock_patch.call_args
            payload = args[1]
            assert payload["action"] == "LAUNCH"
            assert payload["billing_org"] == "test-org"

    def test_ml_models_constants(self):
        """Test ML_MODELS constants are properly defined"""
        expected_models = {
            "orb_d3_v2": 0,
            "mattersim_1_0_0_5m": 1,
            "orb_v3_conservative": 2,
            "esen_30m_oam": 3,
        }

        assert ML_MODELS == expected_models

        # Test all models are valid integers
        for model_name, model_code in ML_MODELS.items():
            assert isinstance(model_name, str)
            assert isinstance(model_code, int)
            assert model_code >= 0
