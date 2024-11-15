from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info


@click.group(name="ea")
def soecexploration_group():
    """Manage EA / SOEC explorations and related resources"""
    pass


@soecexploration_group.command()
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
    """Get EA exploration details or list all explorations"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/ea-exploration/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]EA Exploration Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")
        console.print(f"Status: {get_status_string(result.get('status'))}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at'))}")

        if result.get("parameters"):
            console.print("\n[bold]Parameters[/bold]")
            console.print_json(data=result["parameters"])
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
            results = client.get_all("/api/ea-exploration/", params=params)
        else:
            results = client.get("/api/ea-exploration/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Strains list", "strains_list", None),
            ("Stress algorithm", "stress_algorithm", None),
            ("Stress method", "stress_method", None),
            ("Num last samples", "num_last_samples", None),
            ("Status", "status", get_status_string),
            ("Created", "created_at", format_datetime),
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No EA explorations found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="EA Explorations", caption=footer_string
        )

        console.print(table)


@soecexploration_group.command()
@click.option("--name", required=True, help="Exploration name")
@click.option("--project", required=True, help="Project ID")
@click.option(
    "--relaxed-structure-simulation", help="ID of the relaxed structure simulation"
)
@click.option("--relaxed-structure", help="ID of the relaxed structure upload")
@click.option("--k8s-cluster", help="K8S cluster ID (required for launch)")
@click.option("--launch", is_flag=True, help="Launch the exploration immediately")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create(
    name: str,
    project: str,
    relaxed_structure_simulation: Optional[str] = None,
    relaxed_structure: Optional[str] = None,
    k8s_cluster: Optional[str] = None,
    launch: bool = False,
    json_output: bool = False,
):
    """Create a new SOEC exploration

    Either --relaxed-structure-simulation or --relaxed-structure must be provided.
    If --launch is specified, --k8s-cluster must also be provided.
    """
    client = get_client()
    console = Console()

    # Validate required fields for launch
    if launch and not k8s_cluster:
        console.print("[red]K8S cluster ID is required when launching[/red]")
        return

    # Validate relaxed structure source
    if not (relaxed_structure_simulation or relaxed_structure):
        console.print(
            "[red]Either relaxed-structure-simulation or relaxed-structure must be provided[/red]"
        )
        return

    if relaxed_structure_simulation and relaxed_structure:
        console.print(
            "[red]Cannot specify both relaxed-structure-simulation and relaxed-structure[/red]"
        )
        return

    # Prepare request data
    data = {
        "name": name,
        "project": project,
    }

    if relaxed_structure_simulation:
        data["relaxed_structure_simulation"] = relaxed_structure_simulation
    if relaxed_structure:
        data["relaxed_structure"] = relaxed_structure
    if k8s_cluster:
        data["k8s_cluster"] = k8s_cluster
    if launch:
        data["action"] = "LAUNCH"

    result = client.post("/api/ea-exploration/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created SOEC exploration '{name}' with ID: {result['id']}[/green]"
        )
        if launch:
            console.print("[green]Exploration has been queued for launch[/green]")


@soecexploration_group.command()
@click.argument("id")
def delete(id: str):
    """Delete an EA exploration"""
    client = get_client()
    client.delete(f"/api/ea-exploration/{id}/")
    console = Console()
    console.print(f"[green]Deleted EA exploration {id}[/green]")


@soecexploration_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_sample(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get SOEC exploration sample details or list all samples"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/ea-exploration-sample/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]SOEC Sample Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(
            f"Exploration: {result['exploration'].get('name', 'N/A')} ({result['exploration']['id']})"
        )
        if result.get("task"):
            console.print(f"Status: {get_status_string(result['task'].get('status'))}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at'))}")

        if result.get("parameters"):
            console.print("\n[bold]Parameters[/bold]")
            console.print_json(data=result["parameters"])
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
            results = client.get_all("/api/ea-exploration-sample/", params=params)
        else:
            results = client.get("/api/ea-exploration-sample/", params=params)

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
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No samples found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="SOEC Samples", caption=footer_string
        )

        console.print(table)


@soecexploration_group.command()
@click.option("--exploration", required=True, help="SOEC Exploration ID")
@click.option("--simulation", required=True, help="Simulation ID")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_sample(exploration: str, simulation: str, json_output: bool = False):
    """Create a new SOEC exploration sample"""
    client = get_client()
    console = Console()

    data = {
        "exploration": exploration,
        "simulation": simulation,
    }

    result = client.post("/api/ea-exploration-sample/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created and launched SOEC sample with ID: {result['id']}[/green]"
        )


@soecexploration_group.command()
@click.argument("id")
def delete_sample(id: str):
    """Delete a SOEC exploration sample"""
    client = get_client()
    client.delete(f"/api/ea-exploration-sample/{id}/")
    console = Console()
    console.print(f"[green]Deleted SOEC sample {id}[/green]")


@soecexploration_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term (e.g. name)")
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
    """Get SOEC exploration analysis details or list all analyses"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/ea-exploration-analysis/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]SOEC Analysis Details[/bold]"))
        console.print(f"ID: {result['id']}")
        # console.print(f"Exploration: {result['exploration'].get('name', 'N/A')} ({result['exploration']['id']})")
        console.print(f"Exploration: {result['exploration']}")
        if result.get("task"):
            console.print(f"Status: {get_status_string(result['task'].get('status'))}")
            console.print(f"Progress: {result['task'].get('progress')}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at'))}")
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
            results = client.get_all("/api/ea-exploration-analysis/", params=params)
        else:
            results = client.get("/api/ea-exploration-analysis/", params=params)

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
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No analyses found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="SOEC Analyses", caption=footer_string
        )

        console.print(table)


@soecexploration_group.command()
@click.option("--exploration", required=True, help="SOEC Exploration ID")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_analysis(exploration: str, json_output: bool = False):
    """Create a new SOEC exploration analysis"""
    client = get_client()
    console = Console()

    data = {
        "exploration": exploration,
    }

    result = client.post("/api/ea-exploration-analysis/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created and launched SOEC analysis with ID: {result['id']}[/green]"
        )


@soecexploration_group.command()
@click.argument("id")
def delete_analysis(id: str):
    """Delete a SOEC exploration analysis"""
    client = get_client()
    client.delete(f"/api/ea-exploration-analysis/{id}/")
    console = Console()
    console.print(f"[green]Deleted SOEC analysis {id}[/green]")


@soecexploration_group.command()
@click.option("--analysis", required=True, help="SOEC Analysis ID")
@click.option("--file-id", required=True, help="Upload ID of the file")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_analysis_file(analysis: str, file_id: str, json_output: bool = False):
    """Create a new SOEC exploration analysis file entry"""
    client = get_client()
    console = Console()

    data = {
        "analysis": analysis,
        "user_upload": file_id,
    }

    result = client.post("/api/ea-exploration-analysis-file/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created SOEC analysis file with ID: {result['id']}[/green]"
        )


@soecexploration_group.command()
@click.argument("id")
def delete_analysis_file(id: str):
    """Delete a SOEC exploration analysis file"""
    client = get_client()
    client.delete(f"/api/ea-exploration-analysis-file/{id}/")
    console = Console()
    console.print(f"[green]Deleted SOEC analysis file {id}[/green]")


@soecexploration_group.command()
@click.argument("id", required=False)
@click.option("--analysis", help="Analysis ID to filter by")
@click.option("--file-id", help="Exploration ID to filter by")
@click.option("--ordering", help="Field to order results by")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_analysis_file(
    id: Optional[str] = None,
    analysis: Optional[str] = None,
    file_id: Optional[str] = None,
    ordering: Optional[str] = None,
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get SOEC exploration analysis file details or list all files

    Either --analysis or --file-id must be provided when listing files.
    """
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/ea-exploration-analysis-file/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]SOEC Analysis File Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Analysis: {result['analysis']}")
        console.print(f"File: {result.get('user_upload', {}).get('id')}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
    else:
        if not (analysis or file_id):
            console.print(
                "[red]Either --analysis or --exploration must be provided[/red]"
            )
            return

        if analysis and file_id:
            console.print("[red]Cannot specify both --analysis and --exploration[/red]")
            return

        params = {}
        if analysis:
            params["analysis"] = analysis
        if file_id:
            params["user_upload"] = file_id
        if ordering:
            params["ordering"] = ordering

        if fetch_all:
            results = client.get_all(
                "/api/soec-exploration-analysis-file/", params=params
            )
        else:
            results = client.get("/api/soec-exploration-analysis-file/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Analysis", "analysis", None),
            ("File", "file", None),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No analysis files found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="SOEC Analysis Files",
            caption=footer_string,
        )

        console.print(table)
