"""
Secure storage module for sensitive data like authentication tokens.

This module provides secure storage mechanisms using OS-native credential managers
via the keyring library, with fallback to encrypted file storage when keyring is
unavailable.
"""
import json
import logging
import os
import platform
from typing import Optional

logger = logging.getLogger(__name__)

# Service name for keyring
SERVICE_NAME = "atomict"


def _get_keyring():
    """
    Get the keyring backend, returning None if unavailable or if it's a null backend.
    """
    try:
        import keyring
        from keyring.backends.fail import Keyring as FailKeyring
        
        backend = keyring.get_keyring()
        # Check if backend is a null/fail backend (means no secure storage available)
        if isinstance(backend, FailKeyring):
            logger.debug("Keyring is using a fail/null backend, falling back to encrypted storage")
            return None
        return keyring
    except (ImportError, RuntimeError) as e:
        logger.debug(f"Keyring unavailable: {e}, falling back to encrypted storage")
        return None


def _get_encryption_key():
    """
    Generate or retrieve an encryption key for fallback storage.
    
    The key is derived from machine-specific information to provide some
    level of obfuscation. This is not as secure as OS credential managers
    but better than plaintext.
    """
    from cryptography.fernet import Fernet
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
    import base64
    
    # Use machine-specific info as key material
    machine_id = platform.node() + platform.machine()
    
    # Derive a key using PBKDF2HMAC
    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=b'atomict_salt_v1',  # Static salt for deterministic key
        iterations=100000,
    )
    key = base64.urlsafe_b64encode(kdf.derive(machine_id.encode()))
    return Fernet(key)


def _get_encrypted_storage_path():
    """Get the path for encrypted token storage."""
    system = platform.system()
    if system == "Darwin":  # macOS
        base_path = os.path.join(
            os.path.expanduser("~/Library/Application Support"), "atomict"
        )
    elif system == "Windows":
        base_path = os.path.join(os.environ["APPDATA"], "atomict")
    else:  # Linux and other Unix-like systems
        config_path = os.path.join(os.path.expanduser("~"), ".config")
        base_path = os.path.join(config_path, "atomict")
    
    if not os.path.exists(base_path):
        os.makedirs(base_path, mode=0o700)  # Restrict to user only
    
    return os.path.join(base_path, ".secure_tokens")


def store_token(token: str) -> None:
    """
    Store a token securely.
    
    First attempts to use OS keyring, falls back to encrypted file storage.
    
    Args:
        token: The authentication token to store
    """
    keyring = _get_keyring()
    
    if keyring:
        try:
            keyring.set_password(SERVICE_NAME, "token", token)
            logger.debug("Token stored in system keyring")
            return
        except Exception as e:
            logger.warning(f"Failed to store token in keyring: {e}, falling back to encrypted storage")
    
    # Fallback to encrypted file storage
    try:
        fernet = _get_encryption_key()
        encrypted_token = fernet.encrypt(token.encode())
        
        storage_path = _get_encrypted_storage_path()
        with open(storage_path, 'wb') as f:
            f.write(encrypted_token)
        
        # Set strict permissions (user read/write only)
        os.chmod(storage_path, 0o600)
        logger.debug("Token stored in encrypted file storage")
    except Exception as e:
        logger.error(f"Failed to store token securely: {e}")
        raise RuntimeError(f"Unable to store token securely: {e}")


def get_token() -> Optional[str]:
    """
    Retrieve a stored token securely.
    
    First attempts to retrieve from OS keyring, falls back to encrypted file storage.
    
    Returns:
        The stored token, or None if not found
    """
    keyring = _get_keyring()
    
    if keyring:
        try:
            token = keyring.get_password(SERVICE_NAME, "token")
            if token:
                logger.debug("Token retrieved from system keyring")
                return token
        except Exception as e:
            logger.warning(f"Failed to retrieve token from keyring: {e}, trying encrypted storage")
    
    # Fallback to encrypted file storage
    try:
        storage_path = _get_encrypted_storage_path()
        if not os.path.exists(storage_path):
            return None
        
        with open(storage_path, 'rb') as f:
            encrypted_token = f.read()
        
        fernet = _get_encryption_key()
        token = fernet.decrypt(encrypted_token).decode()
        logger.debug("Token retrieved from encrypted file storage")
        return token
    except FileNotFoundError:
        return None
    except Exception as e:
        logger.warning(f"Failed to retrieve token from encrypted storage: {e}")
        return None


def delete_token() -> None:
    """
    Delete a stored token from all storage locations.
    """
    keyring = _get_keyring()
    
    # Try to delete from keyring
    if keyring:
        try:
            keyring.delete_password(SERVICE_NAME, "token")
            logger.debug("Token deleted from system keyring")
        except Exception as e:
            logger.debug(f"No token in keyring to delete or error: {e}")
    
    # Try to delete from encrypted file storage
    try:
        storage_path = _get_encrypted_storage_path()
        if os.path.exists(storage_path):
            os.remove(storage_path)
            logger.debug("Token deleted from encrypted file storage")
    except Exception as e:
        logger.debug(f"Failed to delete token from encrypted storage: {e}")


def migrate_plaintext_token(plaintext_token: str) -> None:
    """
    Migrate a plaintext token to secure storage.
    
    Args:
        plaintext_token: The plaintext token to migrate
    """
    if plaintext_token:
        store_token(plaintext_token)
        logger.info("Migrated plaintext token to secure storage")

