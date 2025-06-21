"""Unit tests for structure conversion functions."""

import pytest
from unittest.mock import patch, MagicMock

from atomict.structure.conversion import create_supercell, create_conventional_cell
from atomict.exceptions import APIValidationError


class TestCreateSupercell:
    """Test cases for create_supercell function."""

    @patch("atomict.structure.conversion.post")
    def test_create_supercell_with_list_dimensions(self, mock_post):
        """Test creating supercell with list dimensions."""
        mock_post.return_value = {"response": "ok"}

        result = create_supercell("test-uuid-123", [2, 3, 4])

        mock_post.assert_called_once_with(
            "conversion/supercell/",
            {"user_upload": "test-uuid-123", "supercell": "2 3 4"},
        )
        assert result == {"response": "ok"}

    @patch("atomict.structure.conversion.post")
    def test_create_supercell_with_string_dimensions(self, mock_post):
        """Test creating supercell with string dimensions."""
        mock_post.return_value = {"response": "ok"}

        result = create_supercell("test-uuid-123", "2 2 2")

        mock_post.assert_called_once_with(
            "conversion/supercell/",
            {"user_upload": "test-uuid-123", "supercell": "2 2 2"},
        )
        assert result == {"response": "ok"}

    @patch("atomict.structure.conversion.post")
    def test_create_supercell_with_extra_kwargs(self, mock_post):
        """Test creating supercell with additional parameters."""
        mock_post.return_value = {"response": "ok"}

        result = create_supercell("test-uuid-123", [1, 2, 3], extra_param="value")

        mock_post.assert_called_once_with(
            "conversion/supercell/",
            {
                "user_upload": "test-uuid-123",
                "supercell": "1 2 3",
                "extra_param": "value",
            },
        )
        assert result == {"response": "ok"}

    def test_create_supercell_invalid_list_length(self):
        """Test error handling for invalid list length."""
        with pytest.raises(APIValidationError, match="exactly 3 values"):
            create_supercell("test-uuid", [2, 3])

        with pytest.raises(APIValidationError, match="exactly 3 values"):
            create_supercell("test-uuid", [1, 2, 3, 4])

    def test_create_supercell_invalid_list_values(self):
        """Test error handling for invalid list values."""
        with pytest.raises(APIValidationError, match="positive integers"):
            create_supercell("test-uuid", [0, 2, 3])

        with pytest.raises(APIValidationError, match="positive integers"):
            create_supercell("test-uuid", [-1, 2, 3])

        with pytest.raises(APIValidationError, match="positive integers"):
            create_supercell("test-uuid", [1.5, 2, 3])

    def test_create_supercell_invalid_string_format(self):
        """Test error handling for invalid string format."""
        with pytest.raises(APIValidationError, match="exactly 3 values"):
            create_supercell("test-uuid", "2 3")

        with pytest.raises(APIValidationError, match="exactly 3 values"):
            create_supercell("test-uuid", "1 2 3 4")

        with pytest.raises(APIValidationError, match="integers"):
            create_supercell("test-uuid", "a b c")

        with pytest.raises(APIValidationError, match="positive integers"):
            create_supercell("test-uuid", "0 2 3")

    def test_create_supercell_invalid_dimension_type(self):
        """Test error handling for invalid dimension types."""
        with pytest.raises(
            APIValidationError, match="list of integers or space-separated string"
        ):
            create_supercell("test-uuid", 123)

        with pytest.raises(
            APIValidationError, match="list of integers or space-separated string"
        ):
            create_supercell("test-uuid", {"a": 1, "b": 2, "c": 3})

    def test_create_supercell_invalid_user_upload_id(self):
        """Test error handling for invalid user_upload_id."""
        with pytest.raises(APIValidationError, match="non-empty string"):
            create_supercell("", [2, 2, 2])

        with pytest.raises(APIValidationError, match="non-empty string"):
            create_supercell(None, [2, 2, 2])

        with pytest.raises(APIValidationError, match="non-empty string"):
            create_supercell(123, [2, 2, 2])


class TestCreateConventionalCell:
    """Test cases for create_conventional_cell function."""

    @patch("atomict.structure.conversion.post")
    def test_create_conventional_cell_basic(self, mock_post):
        """Test basic conventional cell creation."""
        mock_post.return_value = {"response": "ok"}

        result = create_conventional_cell("test-uuid-123")

        mock_post.assert_called_once_with(
            "conversion/conventional/", {"user_upload": "test-uuid-123"}
        )
        assert result == {"response": "ok"}

    @patch("atomict.structure.conversion.post")
    def test_create_conventional_cell_with_kwargs(self, mock_post):
        """Test conventional cell creation with extra parameters."""
        mock_post.return_value = {"response": "ok"}

        result = create_conventional_cell("test-uuid-123", extra_param="value")

        mock_post.assert_called_once_with(
            "conversion/conventional/",
            {"user_upload": "test-uuid-123", "extra_param": "value"},
        )
        assert result == {"response": "ok"}

    def test_create_conventional_cell_invalid_user_upload_id(self):
        """Test error handling for invalid user_upload_id."""
        with pytest.raises(APIValidationError, match="non-empty string"):
            create_conventional_cell("")

        with pytest.raises(APIValidationError, match="non-empty string"):
            create_conventional_cell(None)

        with pytest.raises(APIValidationError, match="non-empty string"):
            create_conventional_cell(123)
