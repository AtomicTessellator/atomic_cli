import click
from rich.console import Console

from ..core.config import Config, CONFIG_FILE


@click.command(name="token",
               help="Display the current authentication token if it exists.")
def _token():
    config = Config()
    console = Console()
    
    if config.token:
        console.print(f"[green]Current token:[/green] {config.token}")
    else:
        console.print("[yellow]No authentication token found.[/yellow]")
        console.print(f"[yellow]Please run 'tess login' to authenticate and save a token to {CONFIG_FILE}[/yellow]")

