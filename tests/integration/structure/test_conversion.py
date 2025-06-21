"""Integration tests for structure conversion functions."""

import pytest
import os

from atomict.auth import resolve_token
from atomict.structure.conversion import create_supercell, create_conventional_cell
from atomict.exceptions import APIValidationError, PermissionDenied


@pytest.fixture(scope="session", autouse=True)
def setup_authentication():
    """Set up authentication for integration tests."""
    token = resolve_token()
    os.environ["AT_TOKEN"] = token


@pytest.mark.integration
class TestCreateSupercellIntegration:
    """Integration tests for create_supercell function."""

    @pytest.fixture
    def test_structure_id(self):
        """Get a test structure ID from environment."""
        structure_id = os.getenv("TEST_STRUCTURE_ID")
        if not structure_id:
            pytest.skip("TEST_STRUCTURE_ID environment variable not set")
        return structure_id

    def test_create_supercell_success(self, test_structure_id):
        """Test successful supercell creation with real API."""
        result = create_supercell(test_structure_id, [2, 2, 2])

        assert isinstance(result, dict)
        assert result.get("response") == "ok"

    def test_create_supercell_string_dimensions(self, test_structure_id):
        """Test supercell creation with string dimensions."""
        result = create_supercell(test_structure_id, "2 2 2")

        assert isinstance(result, dict)
        assert result.get("response") == "ok"

    def test_create_supercell_different_dimensions(self, test_structure_id):
        """Test supercell creation with non-uniform dimensions."""
        result = create_supercell(test_structure_id, [1, 2, 3])

        assert isinstance(result, dict)
        assert result.get("response") == "ok"

    def test_create_supercell_invalid_structure_id(self):
        """Test error handling for non-existent structure ID."""
        fake_id = "00000000-0000-0000-0000-000000000000"

        with pytest.raises((APIValidationError, PermissionDenied)):
            create_supercell(fake_id, [2, 2, 2])

    def test_create_supercell_malformed_id(self):
        """Test error handling for malformed structure ID."""
        with pytest.raises(APIValidationError):
            create_supercell("not-a-uuid", [2, 2, 2])


@pytest.mark.integration
class TestCreateConventionalCellIntegration:
    """Integration tests for create_conventional_cell function."""

    @pytest.fixture
    def test_structure_id(self):
        """Get a test structure ID from environment."""
        structure_id = os.getenv("TEST_STRUCTURE_ID")
        if not structure_id:
            pytest.skip("TEST_STRUCTURE_ID environment variable not set")
        return structure_id

    def test_create_conventional_cell_success(self, test_structure_id):
        """Test successful conventional cell creation with real API."""
        result = create_conventional_cell(test_structure_id)

        assert isinstance(result, dict)
        assert result.get("response") == "ok"

    def test_create_conventional_cell_invalid_structure_id(self):
        """Test error handling for non-existent structure ID."""
        fake_id = "00000000-0000-0000-0000-000000000000"

        with pytest.raises((APIValidationError, PermissionDenied)):
            create_conventional_cell(fake_id)
