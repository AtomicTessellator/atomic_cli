# cli/commands/task.py
from typing import Optional

import click
from rich.console import Console
import json

from atomict.cli.core.client import get_client
from atomict.cli.commands.helpers import get_status_string, format_datetime
from atomict.cli.commands.common import table_0
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name='task')
def task():
    """Manage tasks and their status"""
    pass


@task.command()
@click.argument('id')
def cancel(id: str):
    """Cancel a running task"""
    client = get_client()
    client.post(f'/api/tasks/{id}/cancel/', {})
    console.print(f"[green]Task {id} has been cancelled[/green]")


@task.command()
@click.argument('id')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def status(id: str, json_output: bool):
    """Get detailed task status"""
    client = get_client()
    result = client.get(f'/api/tasks/{id}/status/')

    if json_output:
        click.echo(json.dumps(result, indent=2))
        return

    console.print("[bold]Task Status Details[/bold]")
    console.print(f"ID: {result['id']}")
    console.print(f"Status: {result.get('status', 'N/A')}")
    console.print(f"Progress: {result.get('progress', 'N/A')}")
    if result.get('error'):
        console.print(f"[red]Error: {result['error']}[/red]")


@task.command()
@click.argument('id', required=False)
@click.option('--search', help='Search tasks by ID, type, status, or error message')
@click.option('--ordering', type=click.Choice(['created_at', '-created_at', 'status', '-status', 
                                             'task_type', '-task_type']), 
              default='-created_at',
              help='Order results by field (prefix with - for descending)')
@click.option('--status', type=click.Choice(['pending', 'running', 'completed', 'failed']),
              help='Filter by status')
@click.option('--all', 'fetch_all', is_flag=True, help='Fetch all results')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get(id: Optional[str] = None, search: Optional[str] = None, 
        ordering: Optional[str] = '-created_at', 
        status: Optional[str] = None,
        fetch_all: bool = False,
        json_output: bool = False):
    """
    Get task details. If no ID is provided, lists all tasks.

    Search:
        Search across multiple fields including task ID, type, status, and error messages.
        Example: --search="training" will find tasks with "training" in any searchable field.

    Ordering options:
        created_at: Order by creation date (ascending)
        -created_at: Order by creation date (descending) [default]
        status: Order by status (ascending)
        -status: Order by status (descending)
        task_type: Order by task type (ascending)
        -task_type: Order by task type (descending)
    """
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

        # console.print(Panel("[bold]Task Details[/bold]"))
        console.print(f"ID:     {result['id']}")
        console.print(f"Type:   {result.get('task_type', 'N/A')}")
        console.print(f"Model Name: {result.get('model_name', 'N/A')}")
        console.print(f"Status: {result.get('status', 'N/A')}")
        # ... detailed view with nested data ...

    else:
        params = {}
        if search:
            params['search'] = search
        if ordering:
            params['ordering'] = ordering
        if status:
            params['status'] = status

        if fetch_all:
            results = client.get_all('/api/tasks/', params=params)
        else:
            results = client.get('/api/tasks/', params=params)
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        results, footer_string = get_pagination_info(results)
        table = table_0
        table.title = "Tasks"
        table.caption = footer_string

        table.add_column("ID")
        table.add_column("Type")
        # table.add_column("Input Params")  # usually blank ATM
        table.add_column("Status")
        table.add_column("Errors")
        table.add_column("Created")

        for task in results:
            table.add_row(
                task.get('id', 'N/A'),
                task.get('task_type', 'N/A'),
                # task.get('input_params', 'N/A'),
                get_status_string(task.get('status')),
                task.get('error', 'N/A'),
                format_datetime(task.get('created_at', 'N/A'))
            )
        console.print(table)
