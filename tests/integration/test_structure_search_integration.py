import os

import pytest

from atomict.structure.search import (
    create_discovery_query,
    delete_discovery_query,
    get_discovery_query,
    list_discovery_queries,
    search_structures,
    search_structures_by_element,
    update_discovery_query,
)


@pytest.mark.integration
class TestStructureSearchIntegration:
    """Integration tests for structure search functionality"""

    def test_search_structures_basic(self):
        """Test basic ASE database structure search"""
        # Test simple formula search
        result = search_structures("formula=H2O", limit=5)

        assert "status" in result
        assert result["status"] == "OK"
        assert "results" in result
        assert isinstance(result["results"], list)
        assert len(result["results"]) <= 5

    def test_search_structures_by_element_basic(self):
        """Test basic element-based structure search"""
        # Test simple element search
        result = search_structures_by_element("Cu")

        assert "results" in result
        assert isinstance(result["results"], list)

        # Check structure of results if any returned
        if result["results"]:
            first_result = result["results"][0]
            assert "mp_id" in first_result
            assert "cif_structure" in first_result
            assert "ase_atoms" in first_result
            assert "data" in first_result

    def test_search_structures_by_element_not_observed(self):
        """Test element search including non-observed structures"""
        result = search_structures_by_element("Li", exp_observed=False)

        assert "results" in result
        assert isinstance(result["results"], list)


@pytest.mark.integration
class TestDiscoveryQueryIntegration:
    """Integration tests for discovery query CRUD operations"""

    @pytest.fixture
    def test_project_id(self):
        """Get test project ID from environment"""
        project_id = os.getenv("TEST_PROJECT_ID")
        if not project_id:
            pytest.skip("TEST_PROJECT_ID environment variable not set")
        return project_id

    @pytest.fixture
    def sample_discovery_query(self, test_project_id):
        """Create a sample discovery query for testing"""
        query_data = create_discovery_query(
            project_id=test_project_id,
            name="Test Integration Query",
            query="SELECT 1 as test_column",
            description="Created by integration test",
        )

        yield query_data

        # Cleanup: delete the query after test
        try:
            delete_discovery_query(query_data["id"])
        except Exception:
            pass  # Ignore cleanup errors

    def test_create_discovery_query(self, test_project_id):
        """Test creating a discovery query"""
        result = create_discovery_query(
            project_id=test_project_id,
            name="Test Query Creation",
            query="SELECT 1 as test",
            description="Integration test query",
        )

        assert "id" in result
        assert result["name"] == "Test Query Creation"
        assert result["query"] == "SELECT 1 as test"
        assert result["description"] == "Integration test query"

        # Cleanup
        delete_discovery_query(result["id"])

    def test_get_discovery_query(self, sample_discovery_query):
        """Test getting a discovery query by ID"""
        query_id = sample_discovery_query["id"]

        result = get_discovery_query(query_id)

        assert result["id"] == query_id
        assert result["name"] == "Test Integration Query"
        assert result["query"] == "SELECT 1 as test_column"

    def test_list_discovery_queries(self, test_project_id, sample_discovery_query):
        """Test listing discovery queries"""
        # List all queries (user has access to)
        all_result = list_discovery_queries()
        assert "results" in all_result or isinstance(all_result, list)

        # List queries for specific project
        project_result = list_discovery_queries(project_id=test_project_id)
        assert "results" in project_result or isinstance(project_result, list)

        # Should find our sample query in project results
        if "results" in project_result:
            queries = project_result["results"]
        else:
            queries = project_result

        query_ids = [q["id"] for q in queries if isinstance(q, dict) and "id" in q]
        assert sample_discovery_query["id"] in query_ids

    def test_update_discovery_query(self, sample_discovery_query):
        """Test updating a discovery query"""
        query_id = sample_discovery_query["id"]

        result = update_discovery_query(
            query_id=query_id,
            name="Updated Test Query",
            description="Updated by integration test",
        )

        assert result["id"] == query_id
        assert result["name"] == "Updated Test Query"
        assert result["description"] == "Updated by integration test"
        # Original query should remain unchanged
        assert result["query"] == "SELECT 1 as test_column"

    def test_delete_discovery_query(self, test_project_id):
        """Test deleting a discovery query"""
        # Create a query to delete
        created = create_discovery_query(
            project_id=test_project_id,
            name="Query to Delete",
            query="SELECT 1",
            description="Will be deleted",
        )

        query_id = created["id"]

        # Delete the query
        result = delete_discovery_query(query_id)

        # Should return success response
        assert result == {"success": True}

        # Verify it's deleted by trying to get it (should fail)
        with pytest.raises(Exception):  # Could be various error types
            get_discovery_query(query_id)


@pytest.mark.integration
class TestErrorHandling:
    """Test error handling in integration scenarios"""

    def test_search_structures_invalid_query(self):
        """Test structure search with invalid query"""
        # ASE database gracefully handles invalid queries by returning empty results
        result = search_structures("invalid_syntax_here")
        assert result["status"] == "OK"
        assert isinstance(result["results"], list)

    def test_search_by_element_invalid(self):
        """Test element search with invalid element"""
        with pytest.raises(ValueError, match="Invalid element symbol"):
            search_structures_by_element("InvalidElement")

    def test_get_nonexistent_discovery_query(self):
        """Test getting non-existent discovery query"""
        with pytest.raises(Exception):  # Could be 404 error or similar
            get_discovery_query("nonexistent-id-12345")

    def test_create_discovery_query_invalid_project(self):
        """Test creating discovery query with invalid project ID"""
        with pytest.raises(Exception):  # Backend validation error
            create_discovery_query(
                project_id="invalid-project-id", name="Test Query", query="SELECT 1"
            )


if __name__ == "__main__":
    # Allow running integration tests directly
    pytest.main([__file__, "-v", "-m", "integration"])
