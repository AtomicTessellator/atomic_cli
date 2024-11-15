from typing import List, Optional

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info


@click.group(name="sqs")
def sqs_group():
    """Manage Special Quasirandom Structure (SQS) explorations"""
    pass


@sqs_group.command()
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
    """Get SQS exploration details or list all explorations"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/sqs-exploration/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]SQS Exploration Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")
        console.print(f"Description: {result.get('description', 'N/A')}")
        console.print(
            f"Task status: {get_status_string(result['task'].get('status', 'N/A'))}"
        )
        console.print(f"Created: {format_datetime(result['created_at'])}")
        console.print(f"Updated: {format_datetime(result.get('updated_at'))}")
        console.print(
            f"Max Size: {'Auto' if result.get('auto_max_size') else result.get('max_size', 'N/A')}"
        )
        console.print(
            f"Atom Count Limit: {result.get('atom_count_upper_limit', 'N/A')}"
        )
        console.print(f"Cluster Cutoffs: {result.get('cluster_cutoffs', [])}")

        if result.get("task"):
            console.print(f"Status: {get_status_string(result['task'].get('status'))}")
            console.print(f"Running on: {result['task'].get('running_on', 'N/A')}")
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
                console.print(f"[red]Invalid filter format: {f}. Use field=value[/red]")
                return

        if fetch_all:
            results = client.get_all("/api/sqs-exploration/", params=params)
        else:
            results = client.get("/api/sqs-exploration/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            (
                "Project",
                "task",
                lambda x: x.get("project", "N/A") if isinstance(x, dict) else None,
            ),
            ("Max Size", "auto_max_size", lambda x: "Auto" if x else None),
            ("Cluster Cutoffs", "cluster_cutoffs", None),
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
                f"[white]No SQS explorations found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="SQS Explorations",
            caption=footer_string,
        )

        console.print(table)


@sqs_group.command()
@click.option("--name", help="Exploration name")
@click.option("--description", help="Exploration description")
@click.option("--project", required=True, help="Project ID")
@click.option("--auto-max-size", is_flag=True, help="Automatically determine max size")
@click.option("--max-size", type=int, help="Maximum structure size")
@click.option("--atom-count-limit", type=int, help="Upper limit for atom count")
@click.option("--cutoffs", type=float, multiple=True, help="Cluster cutoff values")
@click.option("--starting-structure", required=True, help="FHIAims simulation ID")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create(
    name: Optional[str],
    description: Optional[str],
    project: str,
    auto_max_size: bool = False,
    max_size: Optional[int] = None,
    atom_count_limit: Optional[int] = None,
    cutoffs: Optional[List[float]] = None,
    starting_structure: str = None,
    json_output: bool = False,
):
    """Create a new SQS exploration"""
    client = get_client()
    console = Console()

    data = {"project": project, "auto_max_size": auto_max_size}
    if name:
        data["name"] = name
    if description:
        data["description"] = description
    if max_size is not None:
        data["max_size"] = max_size
    if atom_count_limit is not None:
        data["atom_count_upper_limit"] = atom_count_limit
    if cutoffs:
        data["cluster_cutoffs"] = list(cutoffs)
    if starting_structure:
        data["starting_structure"] = starting_structure

    result = client.post("/api/sqs-exploration/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"[green]Created SQS exploration with ID: {result['id']}[/green]")


@sqs_group.command()
@click.argument("id")
def delete(id: str):
    """Delete an SQS exploration"""
    client = get_client()
    client.delete(f"/api/sqs-exploration/{id}/")
    console = Console()
    console.print(f"[green]Deleted SQS exploration {id}[/green]")


@sqs_group.command()
@click.argument("id", required=False)
@click.option("--exploration", help="SQS Exploration ID to filter by")
@click.option("--ordering", help="Field to order results by")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_target_concentration(
    id: Optional[str] = None,
    exploration: Optional[str] = None,
    ordering: Optional[str] = None,
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get SQS target concentration details or list all concentrations

    When listing concentrations, --exploration must be provided.
    """
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/sqs-target-concentration/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]SQS Target Concentration Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Exploration: {result['exploration']}")
        console.print(f"Element: {result.get('element', 'N/A')}")
        console.print(f"Concentration: {result.get('concentration', 'N/A')}")
        console.print(f"Weight: {result.get('weight', 'N/A')}")
    else:
        if not exploration:
            console.print(
                "[red]--exploration is required when listing concentrations[/red]"
            )
            return

        params = {"exploration": exploration}
        if ordering:
            params["ordering"] = ordering

        if fetch_all:
            results = client.get_all("/api/sqs-target-concentration/", params=params)
        else:
            results = client.get("/api/sqs-target-concentration/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Element", "element", None),
            ("Concentration", "concentration", None),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No target concentrations found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="SQS Target Concentrations",
            caption=footer_string,
        )

        console.print(table)


@sqs_group.command()
@click.option("--exploration", required=True, help="SQS Exploration ID")
@click.option("--element", required=True, help="Element symbol")
@click.option(
    "--concentration", required=True, type=float, help="Target concentration (0-1)"
)
@click.option("--weight", type=float, help="Weight for this target concentration")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_target_concentration(
    exploration: str,
    element: str,
    concentration: float,
    weight: float,
    json_output: bool = False,
):
    """Create a new SQS target concentration entry"""
    client = get_client()
    console = Console()

    if not 0 <= concentration <= 1:
        console.print("[red]Concentration must be between 0 and 1[/red]")
        return

    data = {
        "exploration": exploration,
        "element": element,
        "concentration": concentration,
        "weight": weight,
    }

    result = client.post("/api/sqs-target-concentration/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created target concentration for {element} with ID: {result['id']}[/green]"
        )


@sqs_group.command()
@click.argument("id")
def delete_target_concentration(id: str):
    """Delete a SQS target concentration entry"""
    client = get_client()
    client.delete(f"/api/sqs-target-concentration/{id}/")
    console = Console()
    console.print(f"[green]Deleted target concentration {id}[/green]")


@sqs_group.command()
@click.option("--exploration", required=True, help="SQS Exploration ID")
@click.option("--user-upload", required=True, help="Upload ID of the file")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_simulation_file(
    exploration: str, user_upload: str, json_output: bool = False
):
    """Create a new SQS simulation file entry"""
    client = get_client()
    console = Console()

    data = {
        "exploration": exploration,
        "user_upload": user_upload,
    }

    result = client.post("/api/sqs-simulation-file/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created SQS simulation file with ID: {result['id']}[/green]"
        )


@sqs_group.command()
@click.argument("id", required=False)
@click.option("--exploration", help="SQS Exploration ID to filter by")
@click.option("--ordering", help="Field to order results by")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_simulation_file(
    id: Optional[str] = None,
    exploration: Optional[str] = None,
    ordering: Optional[str] = None,
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get SQS simulation file details or list all files

    When listing files, --exploration must be provided.
    """
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/sqs-simulation-file/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]SQS Simulation File Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Exploration: {result['exploration']}")
        console.print(f"Upload ID: {result['user_upload'].get('id')}")
        console.print(f"File Name: {result['user_upload'].get('orig_name')}")
        console.print(f"File Size: {result['user_upload'].get('size')} bytes")
        console.print(
            f"Uploaded: {format_datetime(result['user_upload'].get('uploaded'))}"
        )
    else:
        if not exploration:
            console.print("[red]--exploration is required when listing files[/red]")
            return

        params = {"exploration": exploration}
        if ordering:
            params["ordering"] = ordering

        if fetch_all:
            results = client.get_all("/api/sqs-simulation-file/", params=params)
        else:
            results = client.get("/api/sqs-simulation-file/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            (
                "Upload ID",
                "user_upload",
                lambda x: x["id"] if isinstance(x, dict) else None,
            ),
            (
                "File Name",
                "user_upload",
                lambda x: x["orig_name"] if isinstance(x, dict) else None,
            ),
            (
                "Size",
                "user_upload",
                lambda x: f"{x['size']} bytes" if isinstance(x, dict) else None,
            ),
            (
                "Uploaded",
                "user_upload",
                lambda x: (
                    format_datetime(x["uploaded"]) if isinstance(x, dict) else None
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
            title="SQS Simulation Files",
            caption=footer_string,
        )

        console.print(table)


@sqs_group.command()
@click.argument("id")
def delete_simulation_file(id: str):
    """Delete a SQS simulation file"""
    client = get_client()
    client.delete(f"/api/sqs-simulation-file/{id}/")
    console = Console()
    console.print(f"[green]Deleted SQS simulation file {id}[/green]")
