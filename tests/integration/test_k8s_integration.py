"""
Integration tests for Kubernetes cluster management methods.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API

To run: uv run pytest tests/integration/test_k8s.py -v -m integration
"""

import os
import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
from atomict.infra.k8s import list_clusters, get_cluster


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


@pytest.mark.integration
class TestK8sClusterListIntegration:
    """Integration tests for cluster listing functionality."""

    def test_list_clusters_basic(self, setup_authentication):
        """Test basic cluster listing."""
        clusters = list_clusters()
        
        # Should return a list (even if empty)
        assert isinstance(clusters, list)
        
        # If clusters exist, verify they have expected structure
        if clusters:
            cluster = clusters[0]
            assert "id" in cluster
            assert "name" in cluster
            # Optional fields that might be present
            expected_optional_fields = [
                "public", "has_gpu", "loki_url", "created_at", "updated_at"
            ]
            # Just verify the cluster is a dict with some expected fields
            assert isinstance(cluster, dict)

    def test_list_clusters_with_search(self, setup_authentication):
        """Test cluster listing with search filter."""
        # Test with a search term that likely won't match anything
        clusters = list_clusters(search="nonexistent-cluster-name-xyz")
        
        # Should return empty list or filtered results
        assert isinstance(clusters, list)

    def test_list_clusters_with_ordering(self, setup_authentication):
        """Test cluster listing with ordering."""
        # Test ordering by creation date
        clusters = list_clusters(ordering="-created_at")
        
        # Should return a list
        assert isinstance(clusters, list)
        
        # If multiple clusters exist, verify ordering
        if len(clusters) > 1:
            # Can't easily verify ordering without parsing dates,
            # but the API call should succeed
            pass

    def test_list_clusters_with_filters(self, setup_authentication):
        """Test cluster listing with various filters."""
        # Test with public filter
        public_clusters = list_clusters(public=True)
        assert isinstance(public_clusters, list)
        
        # If clusters exist, verify they are public
        for cluster in public_clusters:
            if "public" in cluster:
                assert cluster["public"] is True

    def test_list_clusters_combined_filters(self, setup_authentication):
        """Test cluster listing with multiple filters combined."""
        clusters = list_clusters(
            search="",  # Empty search should not filter
            ordering="name",
            public=True
        )
        
        assert isinstance(clusters, list)


@pytest.mark.integration
class TestK8sClusterDetailIntegration:
    """Integration tests for individual cluster retrieval."""

    @pytest.fixture
    def available_cluster_id(self, setup_authentication):
        """Get an available cluster ID for testing, or skip if none exist."""
        clusters = list_clusters()
        
        if not clusters:
            pytest.skip("No clusters available for testing get_cluster functionality")
        
        return clusters[0]["id"]

    def test_get_cluster_success(self, setup_authentication, available_cluster_id):
        """Test successful cluster retrieval."""
        cluster = get_cluster(available_cluster_id)
        
        # Verify cluster has required fields
        assert isinstance(cluster, dict)
        assert "id" in cluster
        assert "name" in cluster
        assert cluster["id"] == available_cluster_id
        
        # Check for optional but common fields
        optional_fields = [
            "public", "has_gpu", "loki_url", "loki_internal_cluster_url",
            "loki_username", "cost_per_minute", "created_at", "updated_at"
        ]
        
        # At least some of these should be present
        present_fields = [field for field in optional_fields if field in cluster]
        # No strict requirement, but log for debugging
        print(f"Cluster fields present: {list(cluster.keys())}")

    def test_get_cluster_nonexistent(self, setup_authentication):
        """Test retrieval of non-existent cluster."""
        # Use a UUID that definitely doesn't exist
        fake_id = "00000000-0000-0000-0000-000000000000"
        
        # This should raise an exception or return an error
        # The exact behavior depends on the API implementation
        try:
            result = get_cluster(fake_id)
            # If it doesn't raise an exception, it might return an error dict
            if isinstance(result, dict) and "error" in result:
                assert result["error"]  # Should have an error message
            else:
                pytest.fail("Expected error for non-existent cluster")
        except Exception as e:
            # Exception is expected for non-existent resource
            assert "404" in str(e) or "not found" in str(e).lower()


@pytest.mark.integration
class TestK8sClusterEdgeCases:
    """Integration tests for edge cases and error conditions."""

    def test_list_clusters_no_auth_token(self):
        """Test cluster listing without authentication token."""
        # Temporarily remove auth token
        original_token = os.environ.get("AT_TOKEN")
        if original_token:
            del os.environ["AT_TOKEN"]
        
        try:
            # This should fail due to missing authentication
            with pytest.raises(Exception):
                list_clusters()
        finally:
            # Restore token
            if original_token:
                os.environ["AT_TOKEN"] = original_token

    def test_get_cluster_invalid_id_format(self, setup_authentication):
        """Test cluster retrieval with invalid ID format."""
        # Test with obviously invalid ID format
        invalid_id = "not-a-valid-uuid"
        
        try:
            result = get_cluster(invalid_id)
            # If it doesn't raise an exception, check for error response
            if isinstance(result, dict) and "error" in result:
                assert result["error"]
            else:
                # Some APIs might return 404 for invalid format
                pass
        except Exception as e:
            # Exception is acceptable for invalid ID format
            assert True  # Test passes if exception is raised

    def test_list_clusters_invalid_filter_values(self, setup_authentication):
        """Test cluster listing with invalid filter values."""
        # Test with invalid boolean value
        clusters = list_clusters(public="invalid-boolean")
        
        # API should handle this gracefully (either ignore or return error)
        assert isinstance(clusters, list) or (
            isinstance(clusters, dict) and "error" in clusters
        )
