import sys

import click
from rich.console import Console

from ..core.client import APIClient
from ..core.config import Config, CONFIG_FILE


@click.command(name="login",
               help=f"""Log in to your account and store an authentication token to {CONFIG_FILE}.
                       Previous tokens will be overwritten.""")
@click.option("--username", prompt=True)
@click.option("--password", prompt=True, hide_input=True)
def _login(username: str, password: str):
    client = APIClient()
    response = client.auth(username, password)
    token = response.get("token")
    redirect = response.get("redirect")

    if redirect in ['terms', 'verify_email']:
        Console().print(
            f"[yellow]Terms and account verification need to be accepted. Please use a browser to login: {client.base_url}/login[/yellow]"
        )
        sys.exit(1)

    if token:
        Config().save_token(token)
        Console().print("[green]✓ Successfully logged in![/green]")
    else:
        Console(stderr=True).print(
            "[red]✗ Failed to authenticate. Please check your credentials and/or environment variables.[/red]"
        )
        sys.exit(1)