"""Unit tests for secure_storage module."""
import os
from unittest.mock import MagicMock, patch

import pytest

from atomict.secure_storage import (
    delete_token,
    get_token,
    migrate_plaintext_token,
    store_token,
)


@pytest.fixture
def mock_keyring():
    """Mock keyring for testing."""
    with patch("atomict.secure_storage._get_keyring") as mock:
        mock_kr = MagicMock()
        mock.return_value = mock_kr
        yield mock_kr


@pytest.fixture
def no_keyring():
    """Mock no keyring available (fallback to encrypted storage)."""
    with patch("atomict.secure_storage._get_keyring", return_value=None):
        yield


def test_store_and_get_token_with_keyring(mock_keyring):
    """Test storing and retrieving token using keyring."""
    test_token = "test_token_123"
    
    # Store token
    store_token(test_token)
    
    # Verify keyring.set_password was called
    mock_keyring.set_password.assert_called_once_with("atomict", "token", test_token)
    
    # Mock keyring.get_password to return the token
    mock_keyring.get_password.return_value = test_token
    
    # Retrieve token
    retrieved_token = get_token()
    
    # Verify keyring.get_password was called
    mock_keyring.get_password.assert_called_once_with("atomict", "token")
    assert retrieved_token == test_token


def test_delete_token_with_keyring(mock_keyring):
    """Test deleting token using keyring."""
    delete_token()
    
    # Verify keyring.delete_password was called
    mock_keyring.delete_password.assert_called_once_with("atomict", "token")


def test_store_and_get_token_with_encryption(no_keyring, tmp_path):
    """Test storing and retrieving token using encrypted file storage."""
    with patch("atomict.secure_storage._get_encrypted_storage_path", return_value=str(tmp_path / ".secure_tokens")):
        test_token = "test_token_encrypted_456"
        
        # Store token
        store_token(test_token)
        
        # Verify file was created
        storage_path = tmp_path / ".secure_tokens"
        assert storage_path.exists()
        
        # Verify file permissions (600)
        file_mode = os.stat(storage_path).st_mode & 0o777
        assert file_mode == 0o600
        
        # Retrieve token
        retrieved_token = get_token()
        assert retrieved_token == test_token


def test_delete_token_with_encryption(no_keyring, tmp_path):
    """Test deleting token from encrypted file storage."""
    with patch("atomict.secure_storage._get_encrypted_storage_path", return_value=str(tmp_path / ".secure_tokens")):
        test_token = "test_token_to_delete"
        
        # Store token
        store_token(test_token)
        storage_path = tmp_path / ".secure_tokens"
        assert storage_path.exists()
        
        # Delete token
        delete_token()
        
        # Verify file was deleted
        assert not storage_path.exists()


def test_get_token_not_found(no_keyring, tmp_path):
    """Test retrieving non-existent token."""
    with patch("atomict.secure_storage._get_encrypted_storage_path", return_value=str(tmp_path / ".secure_tokens")):
        # Try to get token that doesn't exist
        retrieved_token = get_token()
        assert retrieved_token is None


def test_migrate_plaintext_token(no_keyring, tmp_path):
    """Test migrating a plaintext token to secure storage."""
    with patch("atomict.secure_storage._get_encrypted_storage_path", return_value=str(tmp_path / ".secure_tokens")):
        plaintext_token = "plaintext_token_789"
        
        # Migrate token
        migrate_plaintext_token(plaintext_token)
        
        # Verify token was stored
        retrieved_token = get_token()
        assert retrieved_token == plaintext_token


def test_keyring_fallback_on_error(mock_keyring, tmp_path):
    """Test fallback to encrypted storage when keyring fails."""
    with patch("atomict.secure_storage._get_encrypted_storage_path", return_value=str(tmp_path / ".secure_tokens")):
        # Make keyring.set_password raise an exception
        mock_keyring.set_password.side_effect = Exception("Keyring error")
        
        test_token = "test_token_fallback"
        
        # Store token (should fallback to encrypted storage)
        store_token(test_token)
        
        # Verify file was created (fallback)
        storage_path = tmp_path / ".secure_tokens"
        assert storage_path.exists()
        
        # Make keyring.get_password also fail
        mock_keyring.get_password.side_effect = Exception("Keyring error")
        
        # Retrieve token (should use encrypted storage)
        retrieved_token = get_token()
        assert retrieved_token == test_token


def test_encryption_key_deterministic():
    """Test that encryption key generation is deterministic."""
    from atomict.secure_storage import _get_encryption_key
    
    # Get key twice
    key1 = _get_encryption_key()
    key2 = _get_encryption_key()
    
    # Verify they're the same (deterministic)
    assert key1._signing_key == key2._signing_key

