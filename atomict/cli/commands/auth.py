import sys

import click
from rich.console import Console

from ..core.client import APIClient
from ..core.config import Config


@click.command()
@click.option("--username", prompt=True)
@click.option("--password", prompt=True, hide_input=True)
def login(username: str, password: str):
    """Log in to your account and store authentication token."""
    client = APIClient()
    config = Config()
    console = Console()
    client.set_auth(username, password)

    # Attempt authentication
    response = client.post("api-auth/", {"username": username, "password": password})
    token = response.get("token")

    if token:
        config.save_token(token)
        client.set_token(token)
        console.print("[green]✓ Successfully logged in![/green]")
    else:
        Console(stderr=True).print(
            "[red]✗ Failed to authenticate. Please check your credentials.[/red]"
        )
        sys.exit(1)
