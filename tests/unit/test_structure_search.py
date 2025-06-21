from unittest.mock import MagicMock, patch

import pytest

from atomict.structure.search import (
    EXAMPLE_QUERIES,
    create_discovery_query,
    delete_discovery_query,
    get_discovery_query,
    list_discovery_queries,
    search_structures,
    search_structures_by_element,
    update_discovery_query,
    validate_ase_query,
    validate_element_symbol,
)


class TestValidation:
    """Test validation functions"""

    def test_validate_element_symbol_valid(self):
        """Test valid element symbols"""
        assert validate_element_symbol("H") is True
        assert validate_element_symbol("He") is True
        assert validate_element_symbol("Cu") is True
        assert validate_element_symbol("Fe") is True

    def test_validate_element_symbol_invalid(self):
        """Test invalid element symbols"""
        assert validate_element_symbol("") is False
        assert validate_element_symbol("h") is False  # lowercase first
        assert validate_element_symbol("HE") is False  # uppercase second
        assert validate_element_symbol("Abc") is False  # too long
        assert validate_element_symbol("123") is False
        assert validate_element_symbol(None) is False

    def test_validate_ase_query_basic(self):
        """Test basic ASE query validation"""
        assert validate_ase_query("formula=Cu2O") is True
        assert validate_ase_query("spacegroup=225") is True
        assert validate_ase_query("") is False
        assert validate_ase_query(None) is False


class TestSearchStructures:
    """Test ASE database search function"""

    @patch("atomict.structure.search.get")
    def test_search_structures_basic(self, mock_get):
        """Test basic structure search"""
        mock_get.return_value = {"status": "OK", "results": []}

        result = search_structures("formula=Cu2O")

        mock_get.assert_called_once_with("structure-search/?query=formula=Cu2O")
        assert result == {"status": "OK", "results": []}

    @patch("atomict.structure.search.get")
    def test_search_structures_with_limit(self, mock_get):
        """Test structure search with custom limit"""
        mock_get.return_value = {"status": "OK", "results": []}

        search_structures("formula=Cu2O", limit=10)

        mock_get.assert_called_once_with("structure-search/?query=formula=Cu2O")

    def test_search_structures_invalid_limit(self):
        """Test structure search with invalid limit"""
        with pytest.raises(
            ValueError, match="Limit must be an integer between 1 and 30"
        ):
            search_structures("formula=Cu2O", limit=0)

        with pytest.raises(
            ValueError, match="Limit must be an integer between 1 and 30"
        ):
            search_structures("formula=Cu2O", limit=31)

        with pytest.raises(
            ValueError, match="Limit must be an integer between 1 and 30"
        ):
            search_structures("formula=Cu2O", limit="10")

    @patch("atomict.structure.search.get")
    def test_search_structures_with_extra_kwargs(self, mock_get):
        """Test structure search with extra parameters"""
        mock_get.return_value = {"status": "OK", "results": []}

        search_structures("formula=Cu2O", extra_kwargs={"custom_param": "value"})

        expected_url = "structure-search/?query=formula=Cu2O&custom_param=value"
        mock_get.assert_called_once_with(expected_url)


class TestSearchStructuresByElement:
    """Test element-based search function"""

    @patch("atomict.structure.search.get")
    def test_search_by_element_basic(self, mock_get):
        """Test basic element search"""
        mock_get.return_value = {"results": []}

        result = search_structures_by_element("Cu")

        mock_get.assert_called_once_with(
            "reality/matprj/element-structures/?element=Cu"
        )
        assert result == {"results": []}

    @patch("atomict.structure.search.get")
    def test_search_by_element_exp_observed_false(self, mock_get):
        """Test element search with exp_observed=False"""
        mock_get.return_value = {"results": []}

        search_structures_by_element("Cu", exp_observed=False)

        expected_url = (
            "reality/matprj/element-structures/?element=Cu&exp_observed=false"
        )
        mock_get.assert_called_once_with(expected_url)

    def test_search_by_element_invalid_element(self):
        """Test element search with invalid element"""
        with pytest.raises(ValueError, match="Invalid element symbol: 'xyz'"):
            search_structures_by_element("xyz")

    @patch("atomict.structure.search.get")
    def test_search_by_element_with_extra_kwargs(self, mock_get):
        """Test element search with extra parameters"""
        mock_get.return_value = {"results": []}

        search_structures_by_element("Cu", extra_kwargs={"custom_param": "value"})

        expected_url = (
            "reality/matprj/element-structures/?element=Cu&custom_param=value"
        )
        mock_get.assert_called_once_with(expected_url)


class TestDiscoveryQuery:
    """Test discovery query CRUD operations"""

    @patch("atomict.structure.search.post")
    def test_create_discovery_query_basic(self, mock_post):
        """Test basic discovery query creation"""
        mock_post.return_value = {"id": "123", "name": "test"}

        result = create_discovery_query("proj123", "Test Query", "SELECT * FROM test")

        expected_payload = {
            "project": "proj123",
            "name": "Test Query",
            "query": "SELECT * FROM test",
            "description": None,
        }
        mock_post.assert_called_once_with(
            "api/discovery-query/", payload=expected_payload
        )
        assert result == {"id": "123", "name": "test"}

    @patch("atomict.structure.search.post")
    def test_create_discovery_query_with_description(self, mock_post):
        """Test discovery query creation with description"""
        mock_post.return_value = {"id": "123"}

        create_discovery_query(
            "proj123",
            "Test Query",
            "SELECT * FROM test",
            description="Test description",
        )

        expected_payload = {
            "project": "proj123",
            "name": "Test Query",
            "query": "SELECT * FROM test",
            "description": "Test description",
        }
        mock_post.assert_called_once_with(
            "api/discovery-query/", payload=expected_payload
        )

    def test_create_discovery_query_missing_required(self):
        """Test discovery query creation with missing required fields"""
        with pytest.raises(ValueError, match="project_id is required"):
            create_discovery_query("", "Test Query", "SELECT * FROM test")

        with pytest.raises(ValueError, match="name is required"):
            create_discovery_query("proj123", "", "SELECT * FROM test")

        with pytest.raises(ValueError, match="query is required"):
            create_discovery_query("proj123", "Test Query", "")

    @patch("atomict.structure.search.get")
    def test_get_discovery_query(self, mock_get):
        """Test getting discovery query by ID"""
        mock_get.return_value = {"id": "123", "name": "test"}

        result = get_discovery_query("123")

        mock_get.assert_called_once_with("api/discovery-query/123/")
        assert result == {"id": "123", "name": "test"}

    @patch("atomict.structure.search.get")
    def test_get_discovery_query_with_params(self, mock_get):
        """Test getting discovery query with parameters"""
        mock_get.return_value = {"id": "123"}

        get_discovery_query("123", include_results=True)

        mock_get.assert_called_once_with(
            "api/discovery-query/123/?include_results=True"
        )

    @patch("atomict.structure.search.get")
    def test_list_discovery_queries_all(self, mock_get):
        """Test listing all discovery queries"""
        mock_get.return_value = {"results": []}

        result = list_discovery_queries()

        mock_get.assert_called_once_with("api/discovery-query/")
        assert result == {"results": []}

    @patch("atomict.structure.search.get")
    def test_list_discovery_queries_by_project(self, mock_get):
        """Test listing discovery queries filtered by project"""
        mock_get.return_value = {"results": []}

        list_discovery_queries(project_id="proj123")

        mock_get.assert_called_once_with("api/discovery-query/?project__id=proj123")

    @patch("atomict.structure.search.patch")
    def test_update_discovery_query(self, mock_put):
        """Test updating discovery query"""
        mock_put.return_value = {"id": "123", "name": "updated"}

        result = update_discovery_query("123", name="Updated Name")

        expected_payload = {"name": "Updated Name"}
        mock_put.assert_called_once_with(
            "api/discovery-query/123/", payload=expected_payload
        )
        assert result == {"id": "123", "name": "updated"}

    def test_update_discovery_query_no_fields(self):
        """Test updating discovery query with no fields"""
        with pytest.raises(
            ValueError, match="At least one field must be provided for update"
        ):
            update_discovery_query("123")

    @patch("atomict.structure.search.delete")
    def test_delete_discovery_query(self, mock_delete):
        """Test deleting discovery query"""
        mock_delete.return_value = {"success": True}

        result = delete_discovery_query("123")

        mock_delete.assert_called_once_with("api/discovery-query/123/")
        assert result == {"success": True}


class TestConstants:
    """Test constants and examples"""

    def test_example_queries_exist(self):
        """Test that example queries are defined"""
        assert isinstance(EXAMPLE_QUERIES, dict)
        assert "formula" in EXAMPLE_QUERIES
        assert "elements" in EXAMPLE_QUERIES
        assert "space_group" in EXAMPLE_QUERIES
        assert len(EXAMPLE_QUERIES) > 0
