# cli/commands/organization.py
import json
import logging
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()
logger = logging.getLogger(__name__)


@click.group(name="org")
def org_group():
    """Manage organizations"""
    pass


@org_group.command()
@click.argument("id", required=False)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def get(id: Optional[str] = None, json_output: bool = False):
    """Get organization details or list all organizations"""
    # WIP
    from atomict.cli.core.client import get_client

    client = get_client()

    if id:
        result = client.get(f"/api/organisation/{id}/")
        if json_output:
            click.echo(json.dumps(result, indent=2))
            return

        console.print(Panel(f"[bold]Organization Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")

        # Show users if present
        if "users" in result:
            console.print("\n[bold]Users[/bold]")
            user_table = Table(show_header=True)
            user_table.add_column("ID")
            user_table.add_column("Username")
            user_table.add_column("Role")

            for user in result["users"]:
                user_table.add_row(
                    str(user["id"]),
                    user.get("username", "N/A"),
                    user.get("role", "N/A"),
                )
            console.print(user_table)
    else:
        results = client.get_all("/api/organisation/")
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        table = Table(show_header=True)
        table.add_column("ID")
        table.add_column("Name")
        table.add_column("Users")

        for org in results:
            table.add_row(
                str(org["id"]), org.get("name", "N/A"), str(len(org.get("users", [])))
            )
        console.print(table)


# Add other org commands...

org = org_group
