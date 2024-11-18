# cli/core/config.py
import json
import os
from pathlib import Path
from typing import Optional

import click
from rich.prompt import Prompt

CONFIG_DIR = Path.home() / ".config" / "atomict"
CONFIG_FILE = CONFIG_DIR / "config.json"


class Config:
    def __init__(self):
        self.username: Optional[str] = os.getenv("AT_USER")
        self.password: Optional[str] = os.getenv("AT_PASS")
        self.token: Optional[str] = None
        self.load_config()

    def load_config(self):
        """Load config from file. """
        if CONFIG_FILE.exists():
            try:
                with open(CONFIG_FILE) as f:
                    config_data = json.load(f)
                    self.token = config_data.get("token")
            except json.JSONDecodeError:
                click.secho("Warning: Invalid config file", fg="yellow", err=True)

    def ensure_auth(self):
        """Ensure we have authentication credentials"""
        if not self.username:
            self.username = Prompt.ask("Enter your username")
        if not self.password:
            self.password = Prompt.ask("Enter your password", password=True)

    def save_token(self, token: str):
        """Save token to config file"""
        CONFIG_DIR.mkdir(parents=True, exist_ok=True)
        with open(CONFIG_FILE, "w") as f:
            json.dump({"token": token}, f)
        self.token = token
