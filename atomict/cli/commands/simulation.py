# cli/commands/simulation.py
import logging
from typing import Any, Dict, Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress
from rich.table import Table

console = Console()
logger = logging.getLogger(__name__)


@click.group(name="simulation")
def simulation():
    """Manage simulations"""
    pass


def format_simulation_table(simulations: list[Dict[str, Any]]) -> Table:
    """Create a table with the most relevant simulation fields"""
    table = Table(show_header=True)

    # Key fields we want to show in the list view
    columns = {
        "ID": lambda x: str(x["id"]),
        "Name": lambda x: x.get("name", "N/A"),
        "Status": lambda x: x.get("task", {}).get("status", "N/A"),
        "Progress": lambda x: f"{x.get('task', {}).get('progress', 0)}%",
        "Results": lambda x: x.get("num_results", "0"),
        "Created": lambda x: x.get("created_at", "N/A"),
    }

    # Add columns
    for col_name in columns:
        table.add_column(col_name)

    # Add rows
    for sim in simulations:
        row = [columns[col_name](sim) for col_name in columns]
        table.add_row(*row)

    return table


def format_simulation_detail(sim: Dict[str, Any]) -> None:
    """Display detailed simulation information"""
    # Main info panel
    console.print(Panel(f"[bold]{sim.get('name', 'Unnamed Simulation')}[/bold]"))

    # Basic info table
    basic_info = Table(show_header=True, title="Basic Information")
    basic_info.add_column("Property")
    basic_info.add_column("Value")

    basic_fields = {
        "ID": "id",
        "Name": "name",
        "Created": "created_at",
        "Edited": "edited_at",
        "Bulk ID": "bulk_id",
        "Slab ID": "slab_id",
        "Calculator": "calculator_class",
        "Max Force": "fmax",
        "Steps": "steps",
    }

    for label, field in basic_fields.items():
        value = sim.get(field, "N/A")
        basic_info.add_row(label, str(value))

    console.print(basic_info)
    console.print()

    # Task status if present
    if task := sim.get("task"):
        task_info = Table(show_header=True, title="Task Status")
        task_info.add_column("Property")
        task_info.add_column("Value")

        task_fields = {
            "Status": "status",
            "Progress": "progress",
            "Running On": "running_on",
            "Error": "error",
        }

        for label, field in task_fields.items():
            value = task.get(field, "N/A")
            task_info.add_row(label, str(value))

        console.print(task_info)
        console.print()

    # Results summary if present
    results = Table(show_header=True, title="Results Summary")
    results.add_column("Metric")
    results.add_column("Count")

    result_fields = {
        "Total Results": "num_results",
        "Dissociated": "num_results_dissociated",
        "Desorbed": "num_results_desorbed",
        "Surface Changed": "num_results_surface_changed",
        "Intercalated": "num_results_intercalated",
    }

    for label, field in result_fields.items():
        value = sim.get(field, "0")
        results.add_row(label, str(value))

    console.print(results)


@simulation.command()
@click.argument("id", required=False)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def get(id: Optional[str] = None, json_output: bool = False):
    """Get simulation details or list all simulations"""
    from atomict.cli.core.client import get_client

    client = get_client()

    if id:
        result = client.get(f"/api/catalysis-simulation/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        format_simulation_detail(result)
    else:
        results = client.get_all("/api/catalysis-simulation/")
        if json_output:
            console.print_json(data=results)
            return

        console.print(format_simulation_table(results))


@simulation.command()
@click.argument("id")
@click.option("--output", "-o", help="Output file path")
def download_results(id: str, output: Optional[str]):
    """Download simulation results file"""
    from pathlib import Path

    from atomict.cli.core.client import get_client

    client = get_client()

    with console.status("[bold green]Downloading simulation results..."):
        response = client.get(f"/api/catalysis-simulation-results/{id}/download/")

        # Get filename from response headers or use provided output path
        if output:
            filename = Path(output)
        else:
            content_disposition = response.headers.get("content-disposition", "")
            filename = (
                Path(content_disposition.split("filename=")[-1].strip('"'))
                if content_disposition
                else Path(f"simulation_{id}_results.dat")
            )

        total_size = int(response.headers.get("content-length", 0))

        with open(filename, "wb") as f:
            with Progress() as progress:
                task = progress.add_task("[cyan]Downloading...", total=total_size)

                for chunk in response.iter_bytes(chunk_size=8192):
                    f.write(chunk)
                    progress.update(task, advance=len(chunk))

        console.print(f"[green]Results downloaded to: {filename}[/green]")
