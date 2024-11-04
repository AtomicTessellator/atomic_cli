from unittest.mock import patch, MagicMock
import pytest
from atomict.cli.core.client import get_client


@pytest.fixture
def mock_config():
    with patch('atomict.cli.core.client.Config') as MockConfig:
        config = MockConfig.return_value
        # Set initial instance attributes
        config.username = 'testuser'
        config.password = 'testpass'
        config.token = None
        # Mock the methods
        config.ensure_auth = MagicMock()
        config.save_token = MagicMock()
        yield config


@pytest.fixture
def mock_client():
    with patch('atomict.cli.core.client.APIClient') as mock:
        client = mock.return_value
        client.post = MagicMock()
        client.set_token = MagicMock()
        yield client


def test_get_client_with_existing_token(mock_config, mock_client):
    # Setup
    mock_config.token = 'existing_token'
    
    # Execute
    client = get_client()
    
    # Assert
    mock_client.set_token.assert_called_once_with('existing_token')
    mock_client.post.assert_not_called()
    mock_config.ensure_auth.assert_not_called()


def test_get_client_auth_flow(mock_config, mock_client):
    # Setup
    mock_config.token = None
    mock_client.post.return_value = {'token': 'new_token'}
    
    # Execute
    client = get_client()
    
    # Assert
    mock_config.ensure_auth.assert_called_once()
    mock_client.post.assert_called_once_with(
        'api-auth/',
        {'username': 'testuser', 'password': 'testpass'}
    )
    mock_config.save_token.assert_called_once_with('new_token')
    mock_client.set_token.assert_called_once_with('new_token')
