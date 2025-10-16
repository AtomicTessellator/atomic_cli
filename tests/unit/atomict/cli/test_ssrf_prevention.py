"""Unit tests for SSRF prevention in pagination."""
import pytest
from unittest.mock import MagicMock, patch

from atomict.cli.core.client import APIClient


@pytest.fixture
def client():
    """Create an API client for testing."""
    return APIClient(base_url="https://api.atomictessellator.com")


class TestSSRFPrevention:
    """Test SSRF prevention in pagination URL handling."""
    
    def test_sanitize_relative_url(self, client):
        """Test that relative URLs pass through unchanged."""
        relative_url = "/api/projects/?page=2"
        result = client._sanitize_pagination_url(relative_url)
        assert result == relative_url
    
    def test_sanitize_valid_absolute_url(self, client):
        """Test that valid absolute URLs matching base_url are converted to relative."""
        absolute_url = "https://api.atomictessellator.com/api/projects/?page=2"
        result = client._sanitize_pagination_url(absolute_url)
        assert result == "/api/projects/?page=2"
    
    def test_sanitize_rejects_external_http_url(self, client):
        """Test that external HTTP URLs are rejected."""
        malicious_url = "http://httpbin.org/ip"
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client._sanitize_pagination_url(malicious_url)
    
    def test_sanitize_rejects_external_https_url(self, client):
        """Test that external HTTPS URLs are rejected."""
        malicious_url = "https://evil.com/api/data"
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client._sanitize_pagination_url(malicious_url)
    
    def test_sanitize_rejects_localhost(self, client):
        """Test that localhost URLs are rejected (SSRF to local services)."""
        localhost_url = "http://localhost:8080/admin"
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client._sanitize_pagination_url(localhost_url)
    
    def test_sanitize_rejects_metadata_endpoint(self, client):
        """Test that cloud metadata endpoints are rejected (AWS IMDS attack)."""
        metadata_url = "http://169.254.169.254/latest/meta-data/"
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client._sanitize_pagination_url(metadata_url)
    
    def test_sanitize_rejects_file_scheme(self, client):
        """Test that file:// URLs are rejected."""
        file_url = "file:///etc/passwd"
        with pytest.raises(ValueError, match="Invalid URL scheme"):
            client._sanitize_pagination_url(file_url)
    
    def test_sanitize_rejects_ftp_scheme(self, client):
        """Test that ftp:// URLs are rejected."""
        ftp_url = "ftp://ftp.example.com/file.txt"
        with pytest.raises(ValueError, match="Invalid URL scheme"):
            client._sanitize_pagination_url(ftp_url)
    
    def test_sanitize_handles_base_url_with_trailing_slash(self, client):
        """Test handling when base_url comparison involves trailing slashes."""
        # Ensure our client's base_url doesn't have trailing slash
        assert not client.base_url.endswith('/')
        
        # Valid URL should still work
        valid_url = "https://api.atomictessellator.com/api/test"
        result = client._sanitize_pagination_url(valid_url)
        assert result == "/api/test"
    
    def test_sanitize_rejects_similar_domain(self, client):
        """Test that similar but different domains are rejected."""
        # Subdomain that's not the base_url
        similar_url = "https://evil.atomictessellator.com/api/data"
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client._sanitize_pagination_url(similar_url)
        
        # Domain with base_url as suffix
        similar_url2 = "https://fakapi.atomictessellator.com/api/data"
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client._sanitize_pagination_url(similar_url2)


class TestPaginationWithSSRFPrevention:
    """Test pagination method with SSRF prevention."""
    
    @patch.object(APIClient, 'get')
    def test_paginate_with_valid_absolute_next(self, mock_get, client):
        """Test pagination works with valid absolute next URLs."""
        # First response with next
        mock_get.return_value = {
            "results": [{"id": 1}, {"id": 2}],
            "next": "https://api.atomictessellator.com/api/projects/?page=2"
        }
        
        # Second response (last page)
        mock_get.side_effect = [
            {
                "results": [{"id": 1}, {"id": 2}],
                "next": "https://api.atomictessellator.com/api/projects/?page=2"
            },
            {
                "results": [{"id": 3}],
                "next": None
            }
        ]
        
        results = list(client.paginate("/api/projects/"))
        
        assert len(results) == 3
        assert results[0]["id"] == 1
        assert results[1]["id"] == 2
        assert results[2]["id"] == 3
    
    @patch.object(APIClient, 'get')
    def test_paginate_with_relative_next(self, mock_get, client):
        """Test pagination works with relative next URLs."""
        mock_get.side_effect = [
            {
                "results": [{"id": 1}],
                "next": "/api/projects/?page=2"
            },
            {
                "results": [{"id": 2}],
                "next": None
            }
        ]
        
        results = list(client.paginate("/api/projects/"))
        
        assert len(results) == 2
        assert mock_get.call_count == 2
    
    @patch.object(APIClient, 'get')
    def test_paginate_rejects_malicious_next(self, mock_get, client):
        """Test that pagination rejects malicious next URLs."""
        # Response with malicious next URL
        mock_get.return_value = {
            "results": [{"id": 1}],
            "next": "http://169.254.169.254/latest/meta-data/"
        }
        
        # Should raise ValueError when trying to process malicious URL
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            list(client.paginate("/api/projects/"))
    
    @patch.object(APIClient, 'get')
    def test_paginate_rejects_external_next(self, mock_get, client):
        """Test that pagination rejects external next URLs."""
        mock_get.return_value = {
            "results": [{"id": 1}],
            "next": "https://evil.com/api/steal-data"
        }
        
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            list(client.paginate("/api/projects/"))
    
    @patch.object(APIClient, 'get')
    def test_get_all_with_ssrf_protection(self, mock_get, client):
        """Test get_all method respects SSRF protection."""
        mock_get.side_effect = [
            {
                "results": [{"id": 1}],
                "next": "http://localhost:8080/admin"  # Malicious
            }
        ]
        
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client.get_all("/api/projects/")


class TestEdgeCases:
    """Test edge cases in SSRF prevention."""
    
    def test_empty_url(self, client):
        """Test handling of empty URL."""
        result = client._sanitize_pagination_url("")
        assert result == ""
    
    def test_url_with_port_in_base_url(self):
        """Test client with port in base_url."""
        client = APIClient(base_url="https://api.example.com:8443")
        
        # Valid URL with port
        valid_url = "https://api.example.com:8443/api/data"
        result = client._sanitize_pagination_url(valid_url)
        assert result == "/api/data"
        
        # Invalid URL (different port)
        invalid_url = "https://api.example.com:9443/api/data"
        with pytest.raises(ValueError, match="Invalid pagination URL"):
            client._sanitize_pagination_url(invalid_url)
    
    def test_url_with_query_params(self, client):
        """Test URLs with complex query parameters."""
        url = "https://api.atomictessellator.com/api/projects/?page=2&limit=10&sort=name"
        result = client._sanitize_pagination_url(url)
        assert result == "/api/projects/?page=2&limit=10&sort=name"
    
    def test_url_with_fragment(self, client):
        """Test URLs with fragments."""
        url = "https://api.atomictessellator.com/api/projects/#section"
        result = client._sanitize_pagination_url(url)
        assert result == "/api/projects/#section"

