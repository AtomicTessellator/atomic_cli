# cli/commands/simulation/catalysis.py
from typing import Any, Dict, Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from atomict.cli.commands.helpers import (
    create_detail_table,
    format_datetime,
    format_task_status,
)
from atomict.cli.core.client import get_client

console = Console()


@click.group(name="catalysis")
def catalysis_group():
    """Manage catalysis simulations"""
    pass


def format_catalysis_table(simulations: list[Dict[str, Any]]) -> Table:
    """Create a table with catalysis simulation fields"""
    table = Table(show_header=True)

    columns = {
        "ID": lambda x: str(x["id"]),
        "Name": lambda x: x.get("name", "N/A"),
        "Status": lambda x: x.get("task", {}).get("status", "N/A"),
        "Progress": lambda x: f"{x.get('task', {}).get('progress', 0)}%",
        "Results": lambda x: x.get("num_results", "0"),
        "Created": lambda x: format_datetime(x.get("created_at", "")),
    }

    for col_name in columns:
        table.add_column(col_name)

    for sim in simulations:
        row = [columns[col_name](sim) for col_name in columns]
        table.add_row(*row)

    return table


def format_catalysis_detail(sim: Dict[str, Any]) -> None:
    """Display detailed catalysis simulation information"""
    # Main info panel
    console.print(Panel(f"[bold]{sim.get('name', 'Unnamed Simulation')}[/bold]"))

    # Basic info
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
    console.print(create_detail_table("Basic Information", sim, basic_fields))

    # Task status
    if task_table := format_task_status(sim.get("task")):
        console.print(task_table)

    # Results summary
    results_fields = {
        "Total Results": "num_results",
        "Dissociated": "num_results_dissociated",
        "Desorbed": "num_results_desorbed",
        "Surface Changed": "num_results_surface_changed",
        "Intercalated": "num_results_intercalated",
    }
    console.print(create_detail_table("Results Summary", sim, results_fields))


def format_results_detail(results: Dict[str, Any]) -> None:
    """Display detailed simulation results"""
    # Energy and convergence info
    energy_table = Table(show_header=True, title="Energy Results")
    energy_table.add_column("Property")
    energy_table.add_column("Value")

    energy_fields = {
        "Total Energy": "energy",
        "ML Energy": "ml_energy",
        "Min ML Energy": "min_ml_energy",
        "Converged": lambda r: "✓" if r.get("converged") else "✗",
        "Force Calls": "force_calls",
    }

    for label, field in energy_fields.items():
        if callable(field):
            value = field(results)
        else:
            value = results.get(field, "N/A")
        energy_table.add_row(label, str(value))

    console.print(energy_table)

    # Structure information
    if structure := results.get("final_structure"):
        structure_table = Table(show_header=True, title="Final Structure")
        structure_table.add_column("Property")
        structure_table.add_column("Value")

        structure_fields = {
            "Atoms": "num_atoms",
            "Cell": "cell",
            "Positions": "positions",
            "Forces": "forces",
        }

        for label, field in structure_fields.items():
            value = structure.get(field, "N/A")
            structure_table.add_row(label, str(value))

        console.print(structure_table)


@catalysis_group.command()
@click.argument("id", required=False)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def get(id: Optional[str] = None, json_output: bool = False):
    """Get catalysis simulation details or list all"""
    # TODO: support more options, filtering, and search
    client = get_client()

    if id:
        result = client.get(f"/api/catalysis-simulation/{id}/")
        if json_output:
            console.print_json(data=result)
            return
        format_catalysis_detail(result)
    else:
        results = client.get_all("/api/catalysis-simulation/")
        if json_output:
            console.print_json(data=results)
            return

        # Add pagination info
        if isinstance(results, dict):
            items = results.get("results", [])
            count = results.get("count")
            page_size = len(items)
            if count and page_size:
                console.print(f"\nShowing page 1 of {(count - 1) // page_size + 1}")
                console.print(f"Total items: {count}")
                if results.get("next"):
                    console.print("Use --all to fetch all results")
        else:
            items = results

        console.print(format_catalysis_table(items))  # Make sure we print the table


@catalysis_group.group(name="results")
def results_group():
    """Manage catalysis simulation results"""
    pass


@results_group.command(name="get")
@click.argument("id", required=False)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--simulation-id", help="Filter by simulation ID")
@click.option("--converged", type=bool, help="Filter by convergence status")
def results_get(
    id: Optional[str] = None,
    json_output: bool = False,
    simulation_id: Optional[str] = None,
    converged: Optional[bool] = None,
):
    """Get simulation results details or list all results"""
    client = get_client()

    params = {}
    if simulation_id:
        params["simulation"] = simulation_id
    if converged is not None:
        params["converged"] = converged

    if id:
        result = client.get(f"/api/catalysis-simulation-results/{id}/")
        if json_output:
            console.print_json(data=result)
            return
        format_results_detail(result)
    else:
        results = client.get_all("/api/catalysis-simulation-results/", params)
        if json_output:
            console.print_json(data=results)
            return

        table = Table(show_header=True)
        table.add_column("ID")
        table.add_column("Simulation")
        table.add_column("Energy")
        table.add_column("ML Energy")
        table.add_column("Converged")
        table.add_column("Created")

        for res in results:
            table.add_row(
                str(res["id"]),
                str(res.get("simulation", "N/A")),
                str(res.get("energy", "N/A")),
                str(res.get("ml_energy", "N/A")),
                "✓" if res.get("converged") else "✗",
                format_datetime(res.get("created_at", "")),
            )
        console.print(table)


# Add results group to catalysis group
catalysis_group.add_command(results_group)

catalysis = catalysis_group
