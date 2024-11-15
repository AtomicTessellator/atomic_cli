from typing import Optional

import click
from rich.console import Console

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="fhiaims")
def fhiaims_group():
    """Manage FHI-aims simulations"""
    pass


@fhiaims_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
# filter should be enum to match DRF config
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
        simulation = client.get(f"/api/fhiaims-simulation/{id}/")
        if json_output:
            console.print_json(data=simulation)
            return
        # Format single simulation output
        console.print(f"ID: {simulation['id']}")
        console.print(f"Name: {simulation.get('name', 'N/A')}")
        console.print(f"Finite diff displacement: {simulation.get('finite_diff_displacement', 'N/A')}")
        console.print(f"Created: {format_datetime(simulation['created_at'])}")
        if simulation.get("task"):
            status = get_status_string(simulation["task"].get("status"))
            console.print(f"Status: {status}")
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

        # Updated this section
        if fetch_all:
            # list
            results = client.get_all("/api/fhiaims-simulation/", params=params)
        else:
            # dict
            results = client.get("/api/fhiaims-simulation/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        # Define columns based on what's already in the table_0 setup
        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Finite diff displacement", "finite_diff_displacement", None),
            ("Task status", "task", lambda x: get_status_string(x.get("status"))),
            ("Created", "created_at", format_datetime),
            (
                "Status",
                "task",
                lambda x: get_status_string(
                    x.get("status") if isinstance(x, dict) else None
                ),
            ),
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
            title="FHI-aims Simulations",
            caption=footer_string,
        )

        console.print(table)


@fhiaims_group.command()
@click.option("--name", help="Simulation name")
@click.option("--description", help="Simulation description")
@click.option("--control-file", required=True, help="Control file content")
@click.option("--geometry-file", required=True, help="Geometry file content")
@click.option(
    "--generate-finite-diff",
    is_flag=True,
    help="Generate finite difference displacements",
)
@click.option(
    "--finite-diff-displacement",
    type=float,
    help="Finite difference displacement value",
)
def create(
    name: Optional[str],
    description: Optional[str],
    control_file: str,
    geometry_file: str,
    generate_finite_diff: bool = False,
    finite_diff_displacement: Optional[float] = None,
):
    """Create a new FHI-aims simulation"""
    client = get_client()
    data = {
        "control_file": control_file,
        "geometry_file": geometry_file,
        "generate_finite_diff_displacements": generate_finite_diff,
    }
    if name:
        data["name"] = name
    if description:
        data["description"] = description
    if finite_diff_displacement is not None:
        data["finite_diff_displacement"] = finite_diff_displacement

    simulation = client.post("/api/fhiaims-simulation/", data=data)
    console.print(f"[green]Created simulation with ID: {simulation['id']}[/green]")


@fhiaims_group.command()
@click.argument("id")
def delete(id: str):
    """Delete a FHI-aims simulation"""
    client = get_client()
    client.delete(f"/api/fhiaims-simulation/{id}/")
    console.print(f"[green]Deleted simulation {id}[/green]")


@fhiaims_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_files(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get simulation file details or list all simulation files"""
    client = get_client()

    if id:
        file = client.get(f"/api/fhiaims-simulation-file/{id}/")
        if json_output:
            console.print_json(data=file)
            return

        # Format single file output
        console.print(f"ID: {file['id']}")
        console.print(f"Simulation ID: {file['simulation']['id']}")
        if file.get("user_upload"):
            console.print(f"File Name: {file['user_upload'].get('users_name', 'N/A')}")
            console.print(
                f"Original Name: {file['user_upload'].get('orig_name', 'N/A')}"
            )
            console.print(f"Size: {file['user_upload'].get('size', 'N/A')} bytes")
            console.print(
                f"Uploaded: {format_datetime(file['user_upload'].get('uploaded'))}"
            )
            if file["user_upload"].get("users_description"):
                console.print(
                    f"Description: {file['user_upload']['users_description']}"
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
            results = client.get_all("/api/fhiaims-simulation-file/", params=params)
        else:
            results = client.get("/api/fhiaims-simulation-file/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            (
                "Simulation",
                "simulation",
                lambda x: x.get("id") if isinstance(x, dict) else None,
            ),
            (
                "File Name",
                "user_upload",
                lambda x: x.get("users_name") if isinstance(x, dict) else None,
            ),
            (
                "Original Name",
                "user_upload",
                lambda x: x.get("orig_name") if isinstance(x, dict) else None,
            ),
            (
                "Size (bytes)",
                "user_upload",
                lambda x: (
                    str(x.get("size"))
                    if isinstance(x, dict) and x.get("size") is not None
                    else None
                ),
            ),
            (
                "Uploaded",
                "user_upload",
                lambda x: (
                    format_datetime(x.get("uploaded")) if isinstance(x, dict) else None
                ),
            ),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No simulation files found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="FHI-aims Simulation Files",
            caption=footer_string,
        )

        console.print(table)
