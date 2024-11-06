# cli/commands/ai_assist.py
import json
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from atomict.cli.core.client import get_client

console = Console()


@click.group(name="ai-assist")
def ai_assist_group():
    """Manage AI assistance features (WIP)"""
    pass


@ai_assist_group.command()
@click.argument("id", required=False)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def get_accepted(id: Optional[str] = None, json_output: bool = False):
    """Get accepted AI suggestions or list all"""
    client = get_client()

    if id:
        result = client.get(f"/api/ai-assist-accepted-code/{id}/")
        if json_output:
            click.echo(json.dumps(result, indent=2))
            return

        console.print(Panel(f"[bold]Accepted AI Suggestion[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Created: {result.get('created_at', 'N/A')}")
        console.print(f"\nCode:\n{result.get('code', 'N/A')}")
    else:
        results = client.get_all("/api/ai-assist-accepted-code/")
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        table = Table(show_header=True)
        table.add_column("ID")
        table.add_column("Created")
        table.add_column("User")

        for item in results:
            table.add_row(
                str(item["id"]), item.get("created_at", "N/A"), item.get("user", "N/A")
            )
        console.print(table)


@ai_assist_group.command()
@click.argument("id", required=False)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def get_votes(id: Optional[str] = None, json_output: bool = False):
    """Get AI suggestion votes or list all"""
    client = get_client()

    if id:
        result = client.get(f"/api/ai-assist-vote/{id}/")
        if json_output:
            click.echo(json.dumps(result, indent=2))
            return

        console.print(Panel(f"[bold]AI Suggestion Vote[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Vote: {result.get('vote', 'N/A')}")
        console.print(f"Created: {result.get('created_at', 'N/A')}")
    else:
        results = client.get_all("/api/ai-assist-vote/")
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        table = Table(show_header=True)
        table.add_column("ID")
        table.add_column("Vote")
        table.add_column("Created")
        table.add_column("User")

        for item in results:
            table.add_row(
                str(item["id"]),
                str(item.get("vote", "N/A")),
                item.get("created_at", "N/A"),
                item.get("user", "N/A"),
            )
        console.print(table)


@ai_assist_group.command()
@click.option(
    "--vote", type=click.Choice(["up", "down"]), required=True, help="Vote up or down"
)
@click.option("--suggestion-id", required=True, help="ID of the AI suggestion")
def vote(vote: str, suggestion_id: str):
    """Vote on an AI suggestion"""
    # TODO: iron out the server side of this
    client = get_client()

    data = {"vote": 1 if vote == "up" else -1, "suggestion": suggestion_id}

    result = client.post("/api/ai-assist-vote/", data)
    console.print(f"[green]Vote recorded: {result['id']}[/green]")


ai_assist = ai_assist_group
