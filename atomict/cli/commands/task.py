# cli/commands/task.py
from typing import Optional

import click
from rich.table import Table
from rich.console import Console
import json
from rich.panel import Panel

from atomict.cli.commands.helpers import get_status_string, format_datetime


console = Console()


@click.group(name='task')
def task():
    """Manage tasks and their status"""
    pass


@task.command()
@click.argument('id')
def cancel(id: str):
    """Cancel a running task"""
    from atomict.cli.core.client import get_client

    client = get_client()
    client.post(f'/api/tasks/{id}/cancel/', {})
    console.print(f"[green]Task {id} has been cancelled[/green]")

@task.command()
@click.argument('id')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def status(id: str, json_output: bool):
    """Get detailed task status"""
    from atomict.cli.core.client import get_client

    client = get_client()
    result = client.get(f'/api/tasks/{id}/status/')

    if json_output:
        click.echo(json.dumps(result, indent=2))
        return

    console.print(f"[bold]Task Status Details[/bold]")
    console.print(f"ID: {result['id']}")
    console.print(f"Status: {result.get('status', 'N/A')}")
    console.print(f"Progress: {result.get('progress', 'N/A')}")
    if result.get('error'):
        console.print(f"[red]Error: {result['error']}[/red]")


@task.command()
@click.argument('id', required=False)
# @click.option('--limit', type=int, help='Number of results to return')
@click.option('--status', type=click.Choice(['pending', 'running', 'completed', 'failed']),
              help='Filter by status')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get(id: Optional[str] = None, limit: Optional[int] = None, status: Optional[str] = None, json_output: bool = False):
    """
    Get task details. If no ID is provided, lists all tasks.
    """
    from atomict.cli.core.client import get_client
    client = get_client()
    params = {}
    # as it stands if you pass limit, you will just get page sizes of 20 instead of 20 total objects
    # if limit:
    #     params['limit'] = limit
    if status:
        params['status'] = status

    if id:
        result = client.get(f'/api/tasks/{id}/')
        if json_output:
            click.echo(json.dumps(result, indent=2))
            return

        console.print(Panel("[bold]Task Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Type: {result.get('task_type', 'N/A')}")
        console.print(f"Status: {result.get('status', 'N/A')}")
        # ... detailed view with nested data ...

    else:
        results = client.get_all('/api/tasks/', params=params)
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        table = Table(show_header=True)
        table.add_column("ID")
        table.add_column("Type")
        table.add_column("Status")
        table.add_column("Created")

        for task in results:
            table.add_row(
                task.get('id', 'N/A'),
                task.get('task_type', 'N/A'),
                get_status_string(task.get('status')),
                format_datetime(task.get('created_at', 'N/A'))
            )
        console.print(table)
