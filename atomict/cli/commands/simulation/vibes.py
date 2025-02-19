from typing import Optional

import click
from rich.console import Console

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="vibes")
def vibes_group():
    """Manage VIBES simulations"""
    pass


@vibes_group.command()
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
    """Get simulation details or list all simulations"""
    client = get_client()

    if id:
        simulation = client.get(f"/api/vibes-simulation/{id}/")
        if json_output:
            console.print_json(data=simulation)
            return
        # Format single simulation output
        console.print(f"ID: {simulation['id']}")
        console.print(f"Name: {simulation.get('name', 'N/A')}")
        console.print(f"Created: {format_datetime(simulation['created_at'])}")
        console.print(f"Updated: {format_datetime(simulation['updated_at'])}")
        console.print(f"Calculator Parameters: {simulation.get('calculator_parameters', {})}")
        console.print(f"Calculator K-points: {simulation.get('calculator_kpoints', {})}")
        console.print(f"Calculator Basis Set: {simulation.get('calculator_basis_set', {})}")
        if simulation.get("task"):
            status = get_status_string(simulation["task"].get("status"))
            console.print(f"Status: {status}")
        if simulation.get("starting_structure"):
            console.print(f"Starting Structure ID: {simulation['starting_structure']['id']}")
    else:
        params = {}
        if search:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(f"Invalid filter format: {f}. Use field=value", err=True)
                return

        if fetch_all:
            results = client.get_all("/api/vibes-simulation/", params=params)
        else:
            results = client.get("/api/vibes-simulation/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Status", "task", lambda x: get_status_string(x.get("status") if isinstance(x, dict) else None)),
            ("Created", "created_at", format_datetime),
            ("Updated", "updated_at", format_datetime),
            ("Starting Structure", "starting_structure", lambda x: x.get("id") if isinstance(x, dict) else None),
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
            title="VIBES Simulations",
            caption=footer_string,
        )

        console.print(table)


@vibes_group.command()
@click.option("--name", help="Simulation name")
@click.option("--description", help="Simulation description")
@click.option("--starting-structure-id", help="ID of the starting FHI-aims structure")
@click.option(
    "--xc",
    type=str,
    help="Phonopy.in setting: Exchange-correlation functional (e.g. 'pw-lda')",
)
@click.option(
    "--kpoint-density",
    type=float,
    help="Phonopy.in setting: K-point density (e.g. 3.5)",
)
@click.option(
    "--basis-set",
    type=str,
    help="Phonopy.in setting: Basis set (e.g. 'light')",
)
def create(
    name: Optional[str],
    description: Optional[str],
    starting_structure_id: Optional[str],
    xc: Optional[str],
    kpoint_density: Optional[float],
    basis_set: Optional[str],
):
    """Create a new VIBES simulation"""
    client = get_client()
    data = {}
    
    if name:
        data["name"] = name
    if description:
        data["description"] = description
    if starting_structure_id:
        data["starting_structure"] = starting_structure_id
    if xc:
        data["calculator_parameters"] = {"xc": xc}
    if kpoint_density:
        data["calculator_kpoints"] = {"density": kpoint_density}
    if basis_set:
        data["calculator_basis_set"] = {"default": basis_set}

    simulation = client.post("/api/vibes-simulation/", data=data)
    console.print(f"[green]Created simulation with ID: {simulation['id']}[/green]")


@vibes_group.command()
@click.argument("id")
def delete(id: str):
    """Delete a VIBES simulation"""
    client = get_client()
    client.delete(f"/api/vibes-simulation/{id}/")
    console.print(f"[green]Deleted simulation {id}[/green]")
