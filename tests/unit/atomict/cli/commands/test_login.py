from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from atomict.cli.commands.login import _login


@pytest.fixture
def mock_client():
    with patch("atomict.cli.core.client.APIClient") as mock:
        client = mock.return_value
        client.post = MagicMock()
        client.auth = MagicMock()
        yield client


def test_login_command_fails_with_error_message(mock_client):
    # Setup
    runner = CliRunner()
    # Mock the auth method to return an empty dict instead of None
    mock_client.auth.return_value = {}

    # Use the context manager to ensure APIClient() returns our mock
    with patch("atomict.cli.commands.login.APIClient", return_value=mock_client):
        # Execute
        result = runner.invoke(_login, input="testuser\ntestpass")

    # Assert
    assert result.exit_code == 1
    assert "failed to authenticate" in result.stderr.lower()
