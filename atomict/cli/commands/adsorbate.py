# cli/commands/adsorbate.py
from typing import Optional, Tuple

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="adsorbate")
def adsorbate_group():
    """Manage adsorbates"""
    pass


@adsorbate_group.command()
@click.argument("id", required=False)
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--search", help="Search term")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--ordering", help="Field to order by (prefix with - for descending)")
@click.option(
    "--all",
    "fetch_all",
    is_flag=True,
    help="Fetch all results (otherwise returns first page)",
)
def get(
    id: Optional[str],
    search: Optional[str] = None,
    filters: Optional[Tuple] = None,  # Fixed syntax error and moved default value
    ordering: Optional[str] = None,
    fetch_all: bool = False,
    json_output: bool = False,
):
    """List adsorbates with optional filtering and search"""
    client = get_client()

    if id:
        result = client.get(f"/api/adsorbate/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel("[bold]Adsorbate Details[/bold]"))
        console.print(f"ID: {result['id']}")
        # Add other relevant fields
    else:
        # Build query parameters
        params = {}
        if search:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering
        if filters:
            for f in filters:
                try:
                    key, value = f.split("=")
                    params[key] = value
                except ValueError:
                    click.echo(
                        f"[red]Invalid filter format: {f}. Use field=value[/red]",
                        err=True,
                    )
                    return

        if fetch_all:
            # Get all results using the paginate method
            results = client.get_all("/api/adsorbate/", params=params)
        else:
            # Get just the first page
            results = client.get("/api/adsorbate/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        # Define columns with optional formatters
        columns = [
            ("ID", "id", None),  # Using default str formatter
            ("SMILES", "smiles", None),
            ("Binding Indices", "binding_indices", None),
            ("Reaction String", "reaction_string", None),
        ]
        items, footer_string = get_pagination_info(results)
        # Create and display the table
        table = create_table(
            columns=columns, items=items, title="Adsorbates", caption=footer_string
        )

        console.print(table)


@adsorbate_group.command()
@click.option("--ase-atoms", required=True, help="ASE atoms string representation")
@click.option("--smiles", help="SMILES string representation of the molecule")
@click.option(
    "--binding-indices",
    multiple=True,
    type=int,
    help="Binding atom indices (can be specified multiple times)",
)
@click.option("--reaction-string", help="Reaction string representation")
def create(
    ase_atoms: Optional[str] = None,
    smiles: Optional[str] = None,
    binding_indices: Optional[tuple[int, ...]] = None,
    reaction_string: Optional[str] = None,
):
    """
    Create a new adsorbate.

    Examples:
        tess adsorbate create --smiles "CC(=O)O" --binding-indices 1 --binding-indices 2
        tess adsorbate create --ase-atoms "Atoms(...)" --reaction-string "A + B -> C"
    """
    client = get_client()

    data = {
        "ase_atoms": ase_atoms,  # this needs work on the server side
        "smiles": smiles,
        "binding_indices": list(binding_indices) if binding_indices else None,
        "reaction_string": reaction_string,
    }
    data = {k: v for k, v in data.items() if v is not None}

    result = client.post("/api/adsorbate/", data)
    console.print(f"[green]Created adsorbate {result['id']}[/green]")


@adsorbate_group.command()
@click.argument("id")
@click.option("--ase-atoms", help="ASE atoms string representation")
@click.option("--smiles", help="SMILES string representation of the molecule")
@click.option(
    "--binding-indices",
    multiple=True,
    type=int,
    help="Binding atom indices (can be specified multiple times)",
)
@click.option("--reaction-string", help="Reaction string representation")
def update(
    id: str,
    ase_atoms: Optional[str] = None,
    smiles: Optional[str] = None,
    binding_indices: Optional[tuple[int, ...]] = None,
    reaction_string: Optional[str] = None,
):
    """
    Update an existing adsorbate.

    Examples:
        tess adsorbate update 123 --smiles "CC(=O)O" --binding-indices 1 --binding-indices 2
        tess adsorbate update 456 --ase-atoms "Atoms(...)" --reaction-string "A + B -> C"
    """
    client = get_client()

    # Build update data excluding None values
    data = {
        "ase_atoms": ase_atoms,
        "smiles": smiles,
        "binding_indices": list(binding_indices) if binding_indices else None,
        "reaction_string": reaction_string,
    }
    data = {k: v for k, v in data.items() if v is not None}

    result = client.put(f"/api/adsorbate/{id}/", data)
    console.print(f"[green]Updated adsorbate {result['id']}[/green]")


@adsorbate_group.command()
@click.argument("id")
def delete(id: str):
    """Delete an adsorbate"""
    client = get_client()

    client.delete(f"/api/adsorbate/{id}/")
    console.print(f"[green]Deleted adsorbate {id}[/green]")
