# cli/commands/k8s.py
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="k8s")
def k8s_group():
    """Manage Kubernetes jobs and clusters"""
    pass


@k8s_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get cluster details or list all clusters"""
    client = get_client()
    console = Console()

    if id:
        cluster = client.get(f"/api/k8s-cluster/{id}/")
        if json_output:
            console.print_json(data=cluster)
            return

        # Format single cluster output
        console.print(Panel(f"[bold]Kubernetes Cluster Details[/bold]"))
        console.print(f"ID: {cluster['id']}")
        console.print(f"Name: {cluster.get('name', 'N/A')}")
        console.print(f"Loki URL: {cluster.get('loki_url', 'N/A')}")
        console.print(f"Active: {cluster.get('active', False)}")
        if cluster.get("description"):
            console.print(f"Description: {cluster['description']}")
        console.print(f"Created: {format_datetime(cluster.get('created_at'))}")
        console.print(f"Updated: {format_datetime(cluster.get('updated_at'))}")

    else:
        params = {}
        if search:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        # Add filter parameters
        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/k8s-cluster/", params=params)
        else:
            results = client.get("/api/k8s-cluster/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Loki URL", "loki_url", None),
            ("Active", "active", str),
            ("Created", "created_at", format_datetime),
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No clusters found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="Kubernetes Clusters",
            caption=footer_string,
        )

        console.print(table)


@k8s_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_job(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get K8s job details or list all jobs"""
    client = get_client()
    console = Console()

    if id:
        job = client.get(f"/api/k8s-job/{id}/")
        if json_output:
            console.print_json(data=job)
            return

        # Format single job output
        console.print(Panel(f"[bold]K8s Job Details[/bold]"))
        console.print(f"ID: {job['id']}")
        console.print(f"Name: {job.get('name', 'N/A')}")
        console.print(f"Status: {job.get('status', 'N/A')}")
        console.print(f"CPU Request: {job.get('cpu_request', 'N/A')}")
        console.print(f"CPU Limit: {job.get('cpu_limit', 'N/A')}")
        console.print(f"Memory Request: {job.get('memory_request', 'N/A')}")
        console.print(f"Memory Limit: {job.get('memory_limit', 'N/A')}")

    else:
        params = {}
        if search:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        # Add filter parameters
        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/k8s-job/", params=params)
        else:
            results = client.get("/api/k8s-job/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Status", "status", None),
            ("CPU Req", "cpu_request", None),
            ("CPU Limit", "cpu_limit", None),
            ("Mem Req", "memory_request", None),
            ("Mem Limit", "memory_limit", None),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No jobs found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="Kubernetes Jobs", caption=footer_string
        )

        console.print(table)


# TBD: not supported
# @k8s.command()
@click.argument("id")
def logs(id: str):
    """Stream logs from a K8s job"""
    client = get_client()

    try:
        # TODO:
        for log_line in client.stream(f"/api/k8s-job/{id}/logs/"):
            console.print(log_line.get("message", ""))
    except KeyboardInterrupt:
        console.print("\n[yellow]Stopped log streaming[/yellow]")
