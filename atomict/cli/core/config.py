# cli/core/config.py
import json
import logging
import os
from pathlib import Path
from typing import Optional

import click
from rich.prompt import Prompt

from atomict.secure_storage import (
    get_token as secure_get_token,
    migrate_plaintext_token,
    store_token as secure_store_token,
)

CONFIG_DIR = Path.home() / ".config" / "atomict"
CONFIG_FILE = CONFIG_DIR / "config.json"

logger = logging.getLogger(__name__)


class Config:
    def __init__(self):
        self.username: Optional[str] = os.getenv("AT_USER")
        self.password: Optional[str] = os.getenv("AT_PASS")
        self.token: Optional[str] = None
        self.load_config()

    def load_config(self):
        """Load config from secure storage and config file."""
        # First, try to get token from secure storage
        self.token = secure_get_token()
        
        # If no token in secure storage, check for migration from plaintext
        if not self.token and CONFIG_FILE.exists():
            try:
                with open(CONFIG_FILE) as f:
                    config_data = json.load(f)
                    plaintext_token = config_data.get("token")
                    if plaintext_token:
                        logger.info("Migrating plaintext token to secure storage")
                        migrate_plaintext_token(plaintext_token)
                        self.token = plaintext_token
                        
                        # Remove token from plaintext storage
                        del config_data["token"]
                        with open(CONFIG_FILE, "w") as fw:
                            json.dump(config_data, fw)
                        # Set strict permissions
                        os.chmod(CONFIG_FILE, 0o600)
            except json.JSONDecodeError:
                click.secho("Warning: Invalid config file", fg="yellow", err=True)
            except Exception as e:
                logger.warning(f"Failed to migrate plaintext token: {e}")

    def ensure_auth(self):
        """Ensure we have authentication credentials"""
        if not self.username:
            self.username = Prompt.ask("Enter your username")
        if not self.password:
            self.password = Prompt.ask("Enter your password", password=True)

    def save_token(self, token: str):
        """Save token to secure storage"""
        CONFIG_DIR.mkdir(parents=True, exist_ok=True, mode=0o700)
        secure_store_token(token)
        self.token = token
