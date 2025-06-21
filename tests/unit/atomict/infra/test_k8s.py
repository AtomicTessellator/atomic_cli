from unittest.mock import patch

import pytest

from atomict.infra.k8s import get_cluster, list_clusters


class TestK8sClusterManagement:
    """Test Kubernetes cluster management functions."""

    @patch("atomict.infra.k8s.get")
    def test_list_clusters_success(self, mock_get):
        """Test successful cluster listing."""
        mock_get.return_value = {
            "results": [
                {"id": "cluster-1", "name": "test-cluster", "public": True},
                {"id": "cluster-2", "name": "private-cluster", "public": False},
            ]
        }

        result = list_clusters()

        assert len(result) == 2
        assert result[0]["name"] == "test-cluster"
        mock_get.assert_called_once_with("api/k8s-cluster/")

    @patch("atomict.infra.k8s.get")
    def test_list_clusters_with_filters(self, mock_get):
        """Test cluster listing with search and filters."""
        mock_get.return_value = {"results": []}

        list_clusters(search="test", ordering="-created_at", public=True)

        expected_path = "api/k8s-cluster/?search=test&ordering=-created_at&public=True"
        mock_get.assert_called_once_with(expected_path)

    @patch("atomict.infra.k8s.get")
    def test_list_clusters_empty_results(self, mock_get):
        """Test cluster listing with empty results."""
        mock_get.return_value = {"results": []}

        result = list_clusters()

        assert result == []

    @patch("atomict.infra.k8s.get")
    def test_list_clusters_direct_response(self, mock_get):
        """Test cluster listing when API returns direct response (not paginated)."""
        mock_get.return_value = [{"id": "cluster-1", "name": "test-cluster"}]

        result = list_clusters()

        assert isinstance(result, list)
        assert result[0]["name"] == "test-cluster"

    @patch("atomict.infra.k8s.get")
    def test_get_cluster_success(self, mock_get):
        """Test successful cluster retrieval."""
        mock_get.return_value = {
            "id": "cluster-1",
            "name": "test-cluster",
            "public": True,
            "has_gpu": False,
            "loki_url": "http://loki:3100",
            "created_at": "2024-01-01T00:00:00Z",
        }

        result = get_cluster("cluster-1")

        assert result["name"] == "test-cluster"
        assert result["public"] is True
        assert result["has_gpu"] is False
        mock_get.assert_called_once_with("api/k8s-cluster/cluster-1/")

    @patch("atomict.infra.k8s.get")
    def test_get_cluster_minimal_response(self, mock_get):
        """Test cluster retrieval with minimal response."""
        mock_get.return_value = {"id": "cluster-1", "name": "minimal-cluster"}

        result = get_cluster("cluster-1")

        assert result["id"] == "cluster-1"
        assert result["name"] == "minimal-cluster"
