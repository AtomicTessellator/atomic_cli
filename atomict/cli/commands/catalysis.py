# cli/commands/exploration.py
from typing import Optional, Tuple

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info


@click.group(name="catalysis")
def catalysis_group():
    """Manage catalysis explorations"""
    pass


@catalysis_group.command()
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
    """Get catalysis exploration details or list all explorations"""
    # WIP
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/catalysis-exploration/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]Catalysis Exploration Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
        console.print(f"Status: {get_status_string(result.get('status', 'N/A'))}")

        if result.get("target_concentrations"):
            console.print("\n[bold]Target Concentrations[/bold]")
            columns = [
                ("Element", "element", None),
                ("Weight", "weight", None),
            ]

            table = create_table(
                columns=columns,
                items=result["target_concentrations"],
                title="Target Concentrations",
            )
            console.print(table)

        if result.get("results"):
            console.print("\n[bold]Exploration Results[/bold]")
            columns = [
                ("Structure ID", "id", None),
                ("Energy", "energy", lambda x: str(x) if x is not None else "N/A"),
                ("Status", "status", get_status_string),
            ]

            table = create_table(
                columns=columns, items=result["results"], title="Results"
            )
            console.print(table)
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

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
            results = client.get_all("/api/catalysis-exploration/", params=params)
        else:
            results = client.get("/api/catalysis-exploration/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Status", "status", get_status_string),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No catalysis explorations found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="Catalysis Explorations",
            caption=footer_string,
        )

        console.print(table)


@catalysis_group.command()
@click.option("--name", required=True, help="Name of the exploration")
@click.option("--structure", required=True, help="Starting structure ID")
@click.option("--project", required=True, help="Project ID")
@click.option("--element", "elements", multiple=True, help="Target elements")
@click.option("--weight", "weights", multiple=True, type=float, help="Target weights")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create(
    name: str,
    structure: str,
    project: str,
    elements: Tuple[str],
    weights: Tuple[float],
    json_output: bool = False,
):
    """Create a new catalysis exploration"""
    client = get_client()
    console = Console()

    if len(elements) != len(weights):
        console.print(
            "[red]Error: Number of elements must match number of weights[/red]"
        )
        return

    data = {
        "name": name,
        "starting_structure": structure,
        "project": project,
        "target_concentrations": [
            {"element": elem, "weight": weight}
            for elem, weight in zip(elements, weights)
        ],
    }

    result = client.post("/api/catalysis-exploration/", data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created catalysis exploration with ID: {result['id']}[/green]"
        )
