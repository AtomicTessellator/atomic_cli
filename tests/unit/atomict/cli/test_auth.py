from unittest.mock import MagicMock, patch

import pytest

from atomict.cli.core.client import get_client


@pytest.fixture
def mock_config():
    with patch("atomict.cli.core.client.Config") as MockConfig:
        config = MockConfig.return_value
        # Set initial instance attributes
        config.username = "testuser"
        config.password = "testpass"
        config.token = None
        # Mock the methods
        config.ensure_auth = MagicMock()
        config.save_token = MagicMock()
        yield config


@pytest.fixture
def mock_client():
    with patch("atomict.cli.core.client.APIClient") as mock:
        client = mock.return_value
        client.post = MagicMock()
        client.set_token = MagicMock()
        client.set_auth = MagicMock()
        client.auth = MagicMock()
        yield client


def test_get_client_with_existing_token(mock_config, mock_client):
    # Setup
    mock_config.token = "existing_token"

    # Execute
    client = get_client()

    # Assert
    mock_client.post.assert_not_called()
    mock_config.ensure_auth.assert_not_called()


def test_get_client_auth_flow(mock_config, mock_client):
    # Setup
    mock_config.token = None
    mock_config.username = "testuser"
    mock_config.password = "testpass"
    mock_client.auth.return_value = {"token": "new_token"}
    mock_client._token = None  # Add this! We need to ensure the token is None

    # Create a new mock for the APIClient class itself
    with patch("atomict.cli.core.client.APIClient") as MockAPIClient:
        # Configure the mock to return our mock_client when instantiated with our config
        MockAPIClient.return_value = mock_client

        # Execute
        client = get_client()

        # Assert
        mock_client.set_token.assert_called_once_with("new_token")
        mock_client.auth.assert_called_once_with()
        mock_config.save_token.assert_called_once_with("new_token")
        mock_client.set_token.assert_called_once_with("new_token")

        # Verify APIClient was constructed with our config
        MockAPIClient.assert_called_once_with(config=mock_config)
