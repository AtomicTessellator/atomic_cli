from unittest.mock import MagicMock, patch
import pytest
from click.testing import CliRunner
from atomict.cli.commands.login import _login
import httpx


@pytest.fixture
def mock_client():
    with patch('atomict.cli.core.client.APIClient') as mock:
        client = mock.return_value
        client.post = MagicMock()
        client.auth = MagicMock()
        yield client


def test_login_command_fails_with_error_message(mock_client):
    # Setup
    runner = CliRunner(mix_stderr=False)
    mock_client.auth.side_effect = httpx.HTTPStatusError(
        message="401 Unauthorized",
        request=httpx.Request("POST", "api-auth/"),
        response=httpx.Response(
            status_code=401,
            request=httpx.Request("POST", "api-auth/"),
            json={"detail": "Invalid credentials"}
        )
    )
    
    # Execute
    result = runner.invoke(_login, input='testuser\ntestpass')
    
    # Assert
    assert result.exit_code == 1
    assert "auth" in result.stderr.lower()
    assert "failed" in result.stderr.lower()
