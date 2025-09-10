from typing import List, Optional

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="kpoint")
def kpoint_group():
    """Manage K-point simulations"""
    pass


@kpoint_group.command()
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
    """Get K-point simulation details or list all simulations"""
    client = get_client()
    console = Console()
    params = {}

    if id:
        simulation = client.get(f"/api/kpoint-simulation/{id}/")
        if json_output:
            console.print_json(data=simulation, params=params)
            return

        # Format single simulation output
        console.print(Panel(f"[bold]K-point Simulation Details[/bold]"))
        console.print(f"ID: {simulation['id']}")
        if simulation.get("exploration"):
            console.print(f"Exploration: {simulation['exploration'].get('name', 'N/A')}")
            status = get_status_string(simulation["exploration"].get("status"))
            console.print(f"Status: {status}")
        console.print(f"K-points: {simulation.get('k_points', [])}")
        if simulation.get("simulation"):
            console.print(
                f"Simulation Name: {simulation['simulation'].get('name', 'N/A')}"
            )
            status = get_status_string(simulation["simulation"].get("status"))
            console.print(f"Status: {status}")

    else:
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
                click.echo(f"Invalid filter format: {f}. Use field=value", err=True)
                return

        if fetch_all:
            results = client.get_all("/api/kpoint-simulation/", params=params)
        else:
            results = client.get("/api/kpoint-simulation/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            (
                "Exploration Name",
                "exploration",
                None
            ),
            (
                "Simulation Name",
                "simulation",
                lambda x: x.get("name") if isinstance(x, dict) else None,
            ),
            ("K-points", "k_points", str),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No simulations found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="K-point Simulations",
            caption=footer_string,
        )

        console.print(table)


@kpoint_group.command()
@click.option("--exploration", required=True, help="Exploration ID")
@click.option(
    "--k-points", required=True, multiple=True, type=float, help="K-point values"
)
def create(exploration: str, k_points: List[float]):
    """Create a new K-point simulation"""
    client = get_client()
    data = {"exploration": exploration, "k_points": list(k_points)}

    simulation = client.post("/api/kpoint-simulation/", data=data)
    console.print(f"[green]Created simulation with ID: {simulation['id']}[/green]")


@kpoint_group.command()
@click.argument("id")
def delete(id: str):
    """Delete a K-point simulation"""
    client = get_client()
    client.delete(f"/api/kpoint-simulation/{id}/")
    console.print(f"[green]Deleted simulation {id}[/green]")


@kpoint_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_exploration(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get K-point exploration details or list all explorations"""
    client = get_client()
    console = Console()

    if id:
        exploration = client.get(f"/api/kpoint-exploration/{id}/")
        if json_output:
            console.print_json(data=exploration)
            return

        # Format single exploration output
        console.print(Panel(f"[bold]K-point Exploration Details[/bold]"))
        console.print(f"ID: {exploration['id']}")
        console.print(f"Name: {exploration.get('name', 'N/A')}")
        console.print(f"Description: {exploration.get('description', 'N/A')}")
        console.print(f"Created: {format_datetime(exploration.get('created_at'))}")
        console.print(f"Updated: {format_datetime(exploration.get('updated_at'))}")
        if exploration.get("task"):
            console.print(
                f"Status: {get_status_string(exploration['task'].get('status'))}"
            )

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
                click.echo(f"Invalid filter format: {f}. Use field=value", err=True)
                return

        if fetch_all:
            results = client.get_all("/api/kpoint-exploration/", params=params)
        else:
            results = client.get("/api/kpoint-exploration/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            (
                "Status",
                "task",
                lambda x: (
                    get_status_string(x.get("status")) if isinstance(x, dict) else None
                ),
            ),
            ("Created", "created_at", format_datetime),
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No explorations found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="K-point Explorations",
            caption=footer_string,
        )

        console.print(table)


@kpoint_group.command()
@click.option("--name", required=True, help="Exploration name")
@click.option("--description", help="Exploration description")
def create_exploration(name: str, description: Optional[str] = None):
    """Create a new K-point exploration"""
    client = get_client()
    data = {
        "name": name,
    }
    if description:
        data["description"] = description

    exploration = client.post("/api/kpoint-exploration/", data=data)
    console.print(f"[green]Created exploration with ID: {exploration['id']}[/green]")


@kpoint_group.command()
@click.argument("id")
def delete_exploration(id: str):
    """Delete a K-point exploration"""
    client = get_client()
    client.delete(f"/api/kpoint-exploration/{id}/")
    console.print(f"[green]Deleted exploration {id}[/green]")


@kpoint_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_analysis(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get K-point analysis details or list all analyses"""
    client = get_client()

    if id:
        analysis = client.get(f"/api/kpoint-analysis/{id}/")
        if json_output:
            console.print_json(data=analysis)
            return

        # Format single analysis output
        console.print(Panel(f"[bold]K-point Analysis Details[/bold]"))
        console.print(f"ID: {analysis['id']}")
        console.print(f"Created: {format_datetime(analysis['created_at'])}")
        if analysis.get("exploration"):
            console.print(
                f"Exploration: {analysis['exploration'].get('name', 'N/A')} ({analysis['exploration']['id']})"
            )
        if analysis.get("task"):
            console.print(
                f"Status: {get_status_string(analysis['task'].get('status'))}"
            )
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
                click.echo(f"Invalid filter format: {f}. Use field=value", err=True)
                return

        if fetch_all:
            results = client.get_all("/api/kpoint-analysis/", params=params)
        else:
            results = client.get("/api/kpoint-analysis/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            (
                "Exploration",
                "exploration",
                lambda x: x.get("name") if isinstance(x, dict) else None,
            ),
            (
                "Status",
                "task",
                lambda x: (
                    get_status_string(x.get("status")) if isinstance(x, dict) else None
                ),
            ),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No analyses found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="K-point Analyses",
            caption=footer_string,
        )

        console.print(table)


@kpoint_group.command()
@click.argument("analysis")
@click.option("--exploration", required=True, help="Exploration ID")
def create_analysis(exploration: str, analysis: str):
    """Create a new K-point analysis. Takes JSON analysis data as main argument."""
    client = get_client()
    # TODO: what data would a user provide?
    data = {
        "analysis": analysis,
    }
    if exploration:
        data["exploration"] = exploration

    analysis = client.post("/api/kpoint-analysis/", data=data)
    console.print(f"[green]Created analysis with ID: {analysis['id']}[/green]")


@kpoint_group.command()
@click.argument("id")
def delete_analysis(id: str):
    """Delete a K-point analysis"""
    client = get_client()
    client.delete(f"/api/kpoint-analysis/{id}/")
    console.print(f"[green]Deleted analysis {id}[/green]")
