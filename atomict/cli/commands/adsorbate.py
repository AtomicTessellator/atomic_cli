# cli/commands/adsorbate.py
import click
from rich.table import Table
from rich.console import Console
from rich.panel import Panel
import json
from typing import Optional
from atomict.cli.core.client import get_client

console = Console()


@click.group(name='adsorbate')
def adsorbate():
    """Manage adsorbates"""
    pass


@adsorbate.command()
@click.argument('id', required=False)
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get(id: Optional[str] = None, json_output: bool = False):
    """Get adsorbate details or list all adsorbates"""
    client = get_client()

    if id:
        result = client.get(f'/api/adsorbate/{id}/')
        if json_output:
            click.echo(json.dumps(result, indent=2))
            return

        console.print(Panel(f"[bold]Adsorbate Details[/bold]"))
        console.print(f"ID: {result['id']}")
        # Add other relevant fields
    else:
        results = client.get_all('/api/adsorbate/')
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        table = Table(show_header=True)
        table.add_column("ID")
        # Add other relevant columns
        
        for item in results:
            table.add_row(str(item['id']))
            # Add other relevant fields


@adsorbate.command()
@click.option('--ase-atoms', required=True, help='ASE atoms string representation')
@click.option('--smiles', help='SMILES string representation of the molecule')
@click.option('--binding-indices', multiple=True, type=int, 
              help='Binding atom indices (can be specified multiple times)')
@click.option('--reaction-string', help='Reaction string representation')
def create(ase_atoms: Optional[str] = None,
          smiles: Optional[str] = None,
          binding_indices: Optional[tuple[int, ...]] = None,
          reaction_string: Optional[str] = None):
    """
    Create a new adsorbate.

    Examples:
        at adsorbate create --smiles "CC(=O)O" --binding-indices 1 --binding-indices 2
        at adsorbate create --ase-atoms "Atoms(...)" --reaction-string "A + B -> C"
    """
    client = get_client()
    
    data = {
        'ase_atoms': ase_atoms,  # this needs work on the server side
        'smiles': smiles,
        'binding_indices': list(binding_indices) if binding_indices else None,
        'reaction_string': reaction_string
    }
    data = {k: v for k, v in data.items() if v is not None}
    
    result = client.post('/api/adsorbate/', data)
    console.print(f"[green]Created adsorbate {result['id']}[/green]")


@adsorbate.command()
@click.argument('id')
@click.option('--ase-atoms', help='ASE atoms string representation')
@click.option('--smiles', help='SMILES string representation of the molecule')
@click.option('--binding-indices', multiple=True, type=int, 
              help='Binding atom indices (can be specified multiple times)')
@click.option('--reaction-string', help='Reaction string representation')
def update(id: str, 
          ase_atoms: Optional[str] = None,
          smiles: Optional[str] = None,
          binding_indices: Optional[tuple[int, ...]] = None,
          reaction_string: Optional[str] = None):
    """
    Update an existing adsorbate.

    Examples:
        at adsorbate update 123 --smiles "CC(=O)O" --binding-indices 1 --binding-indices 2
        at adsorbate update 456 --ase-atoms "Atoms(...)" --reaction-string "A + B -> C"
    """
    client = get_client()
    
    # Build update data excluding None values
    data = {
        'ase_atoms': ase_atoms,
        'smiles': smiles,
        'binding_indices': list(binding_indices) if binding_indices else None,
        'reaction_string': reaction_string
    }
    data = {k: v for k, v in data.items() if v is not None}
    
    result = client.put(f'/api/adsorbate/{id}/', data)
    console.print(f"[green]Updated adsorbate {result['id']}[/green]")


@adsorbate.command()
@click.argument('id')
def delete(id: str):
    """Delete an adsorbate"""
    client = get_client()
    
    client.delete(f'/api/adsorbate/{id}/')
    console.print(f"[green]Deleted adsorbate {id}[/green]")
