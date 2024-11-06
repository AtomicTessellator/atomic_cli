# cli/commands/task.py
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="task")
def task():
    """Manage tasks and their status"""
    pass


@task.command()
@click.argument("id")
def cancel(id: str):
    """Cancel a running task"""
    client = get_client()

    data = {"status": 6}

    client.patch(f"/api/tasks/{id}/", data=data)
    console.print(f"[green]Task {id} has been cancelled[/green]")


@task.command()
@click.argument("id", required=False)
@click.option("--search", help="Search tasks by ID, type, status, or error message")
@click.option(
    "--ordering",
    type=click.Choice(
        ["created_at", "-created_at", "status", "-status", "task_type", "-task_type"]
    ),
    default="-created_at",
    help="Order results by field (prefix with - for descending)",
)
@click.option(
    "--status",
    type=click.Choice(["pending", "running", "completed", "failed", "aborted"]),
    help="Filter by status",
)
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def get(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = "-created_at",
    status: Optional[str] = None,
    fetch_all: bool = False,
    json_output: bool = False,
):
    """Get task details or list all tasks

    When no ID is provided, lists all tasks with optional filtering and ordering.
    """
    client = get_client()

    if id:
        result = client.get(f"/api/tasks/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]Task Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Type: {result.get('task_type', 'N/A')}")
        console.print(f"Model Name: {result.get('model_name', 'N/A')}")
        console.print(f"Status: {get_status_string(result.get('status', 'N/A'))}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at'))}")
        if result.get("error"):
            console.print(f"[red]Error: {result['error']}[/red]")
        if result.get("input_params"):
            console.print("\n[bold]Input Parameters:[/bold]")
            console.print_json(data=result["input_params"])
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering
        if status:
            params["status"] = status

        if fetch_all:
            results = client.get_all("/api/tasks/", params=params)
        else:
            results = client.get("/api/tasks/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Type", "task_type", None),
            ("Status", "status", get_status_string),
            ("Error", "error", None),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No tasks found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="Tasks", caption=footer_string
        )

        console.print(table)


@task.command()
@click.argument("task_id")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def get_status_history(task_id: str, json_output: bool = False):
    """Get status history for a task"""
    client = get_client()
    results = client.get("/api/task-status-history/", params={"id": task_id})

    if json_output:
        console.print_json(data=results)
        return

    columns = [
        ("NewStatus", "new_status", get_status_string),
        ("Timestamp", "timestamp", format_datetime),
    ]

    items, footer_string = get_pagination_info(results)

    if not items:
        console.print(f"[white]No status history found for task {task_id}[/white]")
        return

    table = create_table(
        columns=columns,
        items=items,
        title=f"Status History: {task_id}",
        caption=footer_string,
    )

    console.print(table)
