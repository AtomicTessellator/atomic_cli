from pathlib import Path
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


PSEUDO_LIBRARY_CHOICES = ["efficiency", "precision", "sg15_oncv"]
EXECUTABLE_CHOICES = ["pw.x", "dos.x", "projwfc.x", "bands.x"]


def _read_input_file(input_file: Optional[str], input_path: Optional[str]) -> str:
    """Resolve the QE input file content from either an inline string or a path."""
    if input_file and input_path:
        raise click.UsageError(
            "Provide either --input-file or --input-path, not both"
        )
    if input_path:
        path = Path(input_path)
        if not path.is_file():
            raise click.UsageError(f"Input file not found: {input_path}")
        return path.read_text()
    if input_file:
        return input_file
    raise click.UsageError("One of --input-file or --input-path is required")


@click.group(name="qe")
def qe_group():
    """Manage Quantum ESPRESSO simulations"""
    pass


@qe_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option(
    "--filter", "filters", multiple=True, help="Filter in format field=value"
)
@click.option("--project", help="Filter by project ID")
@click.option(
    "--include-results",
    is_flag=True,
    help="Include parsed results (atoms, pwscf, dos, ...) when fetching a single sim",
)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    project: Optional[str] = None,
    include_results: bool = False,
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get QE simulation details or list all simulations"""
    client = get_client()

    if id:
        params = {}
        if include_results:
            params["include_results"] = "true"
        simulation = client.get(f"/api/qe-simulation/{id}/", params=params)
        if json_output:
            console.print_json(data=simulation)
            return

        console.print(Panel("[bold]Quantum ESPRESSO Simulation[/bold]"))
        console.print(f"ID: {simulation['id']}")
        console.print(f"Name: {simulation.get('name', 'N/A')}")
        console.print(f"Description: {simulation.get('description', 'N/A')}")
        console.print(f"Executable: {simulation.get('executable', 'N/A')}")
        console.print(
            f"Pseudo library: {simulation.get('pseudopotential_library', 'N/A')}"
        )
        console.print(f"Calc type: {simulation.get('calc_type', 'N/A')}")
        console.print(
            f"Convergence achieved: {simulation.get('convergence_achieved', 'N/A')}"
        )
        console.print(f"Created: {format_datetime(simulation.get('created_at'))}")
        if simulation.get("task"):
            status = get_status_string(simulation["task"].get("status"))
            console.print(f"Status: {status}")
        return

    params = {}
    if search:
        params["search"] = search
    if ordering:
        params["ordering"] = ordering
    if project:
        params["task__project"] = project

    for f in filters:
        try:
            field, value = f.split("=", 1)
            params[field] = value
        except ValueError:
            click.echo(f"Invalid filter format: {f}. Use field=value", err=True)
            return

    if fetch_all:
        results = client.get_all("/api/qe-simulation/", params=params)
    else:
        results = client.get("/api/qe-simulation/", params=params)

    if json_output:
        console.print_json(data=results)
        return

    columns = [
        ("ID", "id", None),
        ("Name", "name", None),
        ("Executable", "executable", None),
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
        title="Quantum ESPRESSO Simulations",
        caption=footer_string,
    )
    console.print(table)


@qe_group.command()
@click.option("--project", required=True, help="Project ID to create the simulation in")
@click.option("--name", help="Simulation name")
@click.option("--description", help="Simulation description")
@click.option(
    "--input-file",
    help="Inline QE input file content (use --input-path for a file on disk)",
)
@click.option(
    "--input-path",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to a QE input file (.in / .pwi) to read",
)
@click.option(
    "--executable",
    type=click.Choice(EXECUTABLE_CHOICES),
    help="Override the QE executable (auto-derived from input if omitted)",
)
@click.option(
    "--pseudo-library",
    type=click.Choice(PSEUDO_LIBRARY_CHOICES),
    help=(
        "Public pseudopotential library to auto-fetch (efficiency, precision, "
        "sg15_oncv). Omit to attach UPFs yourself via qe-simulation-file."
    ),
)
@click.option(
    "--parent-simulation",
    help="Parent QE simulation ID (e.g. SCF parent for an NSCF/bands run)",
)
@click.option(
    "--cluster",
    "k8s_cluster",
    help="K8s cluster ID to deploy to (required for --launch unless --batch-config is set)",
)
@click.option(
    "--batch-config",
    help="GCP batch config ID to deploy to (required for --launch unless --cluster is set)",
)
@click.option(
    "--launch",
    is_flag=True,
    help="Immediately launch the simulation after creating it (default: DRAFT)",
)
@click.option("--json-output", is_flag=True, help="Output created simulation as JSON")
def create(
    project: str,
    name: Optional[str],
    description: Optional[str],
    input_file: Optional[str],
    input_path: Optional[str],
    executable: Optional[str],
    pseudo_library: Optional[str],
    parent_simulation: Optional[str],
    k8s_cluster: Optional[str],
    batch_config: Optional[str],
    launch: bool,
    json_output: bool,
):
    """Create (and optionally launch) a Quantum ESPRESSO simulation"""
    content = _read_input_file(input_file, input_path)

    if launch and not (k8s_cluster or batch_config):
        raise click.UsageError(
            "--launch requires --cluster (k8s) or --batch-config (GCP)"
        )

    data = {
        "project": project,
        "input_file": content,
        "action": "LAUNCH" if launch else "DRAFT",
    }
    if name:
        data["name"] = name
    if description:
        data["description"] = description
    if executable:
        data["executable"] = executable
    if pseudo_library:
        data["pseudopotential_library"] = pseudo_library
    if parent_simulation:
        data["parent_simulation"] = parent_simulation
    if k8s_cluster:
        data["task__k8s_cluster"] = k8s_cluster
    if batch_config:
        data["task__batch_config"] = batch_config

    client = get_client()
    simulation = client.post("/api/qe-simulation/", data=data)

    if json_output:
        console.print_json(data=simulation)
        return

    verb = "Launched" if launch else "Created"
    console.print(
        f"[green]{verb} QE simulation with ID: {simulation['id']}[/green]"
    )


@qe_group.command()
@click.argument("id")
@click.option(
    "--cluster",
    "k8s_cluster",
    help="K8s cluster ID to deploy to (required unless --batch-config is set or the sim already has one)",
)
@click.option(
    "--batch-config",
    help="GCP batch config ID to deploy to (required unless --cluster is set or the sim already has one)",
)
@click.option("--json-output", is_flag=True, help="Output updated simulation as JSON")
def launch(
    id: str,
    k8s_cluster: Optional[str],
    batch_config: Optional[str],
    json_output: bool,
):
    """Launch an existing draft QE simulation"""
    data = {"action": "LAUNCH"}
    if k8s_cluster:
        data["task__k8s_cluster"] = k8s_cluster
    if batch_config:
        data["task__batch_config"] = batch_config

    client = get_client()
    simulation = client.patch(f"/api/qe-simulation/{id}/", data=data)

    if json_output:
        console.print_json(data=simulation)
        return

    console.print(f"[green]Launched QE simulation {id}[/green]")


@qe_group.command()
@click.argument("id")
def delete(id: str):
    """Delete a QE simulation"""
    client = get_client()
    client.delete(f"/api/qe-simulation/{id}/")
    console.print(f"[green]Deleted QE simulation {id}[/green]")


@qe_group.command(name="get-files")
@click.argument("id", required=False)
@click.option("--simulation", help="Filter files by parent simulation ID")
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option(
    "--filter", "filters", multiple=True, help="Filter in format field=value"
)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_files(
    id: Optional[str] = None,
    simulation: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get QE simulation file details or list all simulation files"""
    client = get_client()

    if id:
        file = client.get(f"/api/qe-simulation-file/{id}/")
        if json_output:
            console.print_json(data=file)
            return

        console.print(f"ID: {file['id']}")
        sim = file.get("simulation") or {}
        console.print(f"Simulation ID: {sim.get('id', 'N/A')}")
        if file.get("user_upload"):
            uu = file["user_upload"]
            console.print(f"File Name: {uu.get('users_name', 'N/A')}")
            console.print(f"Original Name: {uu.get('orig_name', 'N/A')}")
            console.print(f"Size: {uu.get('size', 'N/A')} bytes")
            console.print(f"Uploaded: {format_datetime(uu.get('uploaded'))}")
            if uu.get("users_description"):
                console.print(f"Description: {uu['users_description']}")
        return

    params = {}
    if simulation:
        params["simulation_uuid"] = simulation
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
        results = client.get_all("/api/qe-simulation-file/", params=params)
    else:
        results = client.get("/api/qe-simulation-file/", params=params)

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
        title="Quantum ESPRESSO Simulation Files",
        caption=footer_string,
    )
    console.print(table)
