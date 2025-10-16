import json
import logging
import os
import platform

from atomict.secure_storage import (
    delete_token as secure_delete_token,
    get_token as secure_get_token,
    migrate_plaintext_token,
    store_token as secure_store_token,
)

DEFAULT_ATOMICT_API_ROOT = "https://api.atomictessellator.com"

logger = logging.getLogger(__name__)


def store(key, value):
    """Store a configuration value. Tokens are stored securely."""
    # Special handling for tokens - use secure storage
    if key == "token":
        secure_store_token(value)
        return
    
    # For other keys, use the config file
    cfg_path = get_config_path("config.json")

    if os.path.exists(cfg_path):
        with open(get_config_path("config.json"), 'r') as f:
            config = json.load(f)
            config[key] = value
    else:
        config = {key: value}

    with open(get_config_path("config.json"), 'w') as f:
        json.dump(config, f)
    
    # Set strict permissions on config file
    os.chmod(cfg_path, 0o600)


def get(key):
    """Retrieve a configuration value. Tokens are retrieved from secure storage."""
    # Special handling for tokens - use secure storage
    if key == "token":
        token = secure_get_token()
        if token:
            return token
        
        # Migration path: check if token exists in old plaintext storage
        try:
            cfg_path = get_config_path("config.json")
            if os.path.exists(cfg_path):
                with open(cfg_path) as f:
                    config = json.load(f)
                    plaintext_token = config.get("token")
                    if plaintext_token:
                        logger.info("Migrating plaintext token to secure storage")
                        migrate_plaintext_token(plaintext_token)
                        # Remove from plaintext storage
                        del config["token"]
                        with open(cfg_path, 'w') as fw:
                            json.dump(config, fw)
                        os.chmod(cfg_path, 0o600)
                        return plaintext_token
        except Exception as e:
            logger.warning(f"Failed to migrate plaintext token: {e}")
        
        return None
    
    # For other keys, use the config file
    with open(get_config_path("config.json")) as f:
        config = json.load(f)
        return config.get(key)


def delete(key):
    """Delete a configuration value. Tokens are deleted from secure storage."""
    # Special handling for tokens - use secure storage
    if key == "token":
        secure_delete_token()
        # Also try to remove from plaintext storage if it exists
        try:
            cfg_path = get_config_path("config.json")
            if os.path.exists(cfg_path):
                with open(cfg_path) as f:
                    config = json.load(f)
                if "token" in config:
                    del config["token"]
                    with open(cfg_path, "w") as f:
                        json.dump(config, f)
                    os.chmod(cfg_path, 0o600)
        except Exception as e:
            logger.debug(f"Error cleaning up plaintext token: {e}")
        return
    
    # For other keys, use the config file
    with open(get_config_path("config.json")) as f:
        config = json.load(f)
        del config[key]

    with open(get_config_path("config.json"), "w") as f:
        json.dump(config, f)
    
    os.chmod(get_config_path("config.json"), 0o600)


def clear():
    """Clear all configuration, including secure token storage."""
    # Delete token from secure storage
    secure_delete_token()
    
    # Delete config file
    path = get_config_path("config.json")
    if os.path.exists(path):
        os.remove(path)


def get_config_path(filename):
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

    # Create directory if it doesn't exist with secure permissions
    if not os.path.exists(base_path):
        os.makedirs(base_path, mode=0o700)  # Only user can read/write/execute

    return os.path.join(base_path, filename)


def check_api_root():
    api_root = os.environ.get("ATOMICT_API_ROOT")

    if api_root is not None:
        # Make sure it's a valid URL without a trailing slash
        if not api_root.startswith("http"):
            raise ValueError(
                "ATOMICT_API_ROOT must be a valid URL, currently: {api_root}"
            )
        if api_root.endswith("/"):
            raise ValueError(
                "ATOMICT_API_ROOT must not end with a trailing slash, currently: {api_root}"
            )


def check_environment():
    check_api_root()
