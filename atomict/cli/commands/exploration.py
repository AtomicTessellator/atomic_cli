# cli/commands/exploration.py
import click
from rich.table import Table
from rich.console import Console
from rich.panel import Panel
import json
from typing import Optional

console = Console()

@click.group(name='exploration')
def exploration():
    """Manage SQS explorations and analyses"""
    pass


@exploration.command()
@click.argument('id', required=False)
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get(id: Optional[str] = None, json_output: bool = False):
    """
    Get exploration details. If no ID is provided, lists all explorations.
    """
    from atomict.cli.core.client import get_client
    client = get_client()

    if id:
        result = client.get(f'/api/catalysis-exploration/{id}/')
        if json_output:
            click.echo(json.dumps(result, indent=2))
            return

        console.print(Panel(f"[bold]Exploration Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")
        console.print(f"Created: {result.get('created_at', 'N/A')}")

        # Target concentrations table
        if 'target_concentrations' in result:
            console.print("\n[bold]Target Concentrations[/bold]")
            conc_table = Table(show_header=True)
            conc_table.add_column("Element")
            conc_table.add_column("Weight")
            
            for conc in result['target_concentrations']:
                conc_table.add_row(
                    conc['element'],
                    str(conc['weight'])
                )
            console.print(conc_table)
            console.print()

        # Results table
        if 'results' in result:
            console.print("[bold]Exploration Results[/bold]")
            results_table = Table(show_header=True)
            results_table.add_column("Structure ID")
            results_table.add_column("Energy")
            results_table.add_column("Status")
            
            for res in result['results']:
                results_table.add_row(
                    str(res['id']),
                    str(res.get('energy', 'N/A')),
                    res.get('status', 'N/A')
                )
            console.print(results_table)

    else:
        results = client.get_all('/api/catalysis-exploration/')
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        table = Table(show_header=True)
        table.add_column("ID")
        table.add_column("Name")
        table.add_column("Created")
        table.add_column("Status")

        for exp in results:
            table.add_row(
                str(exp['id']),
                exp.get('name', 'N/A'),
                exp.get('created_at', 'N/A'),
                exp.get('status', 'N/A')
            )
        console.print(table)


@exploration.command()
@click.option('--name', required=True, help='Name of the exploration')
@click.option('--structure', required=True, help='Starting structure ID')
@click.option('--element', 'elements', multiple=True, help='Target elements')
@click.option('--weight', 'weights', multiple=True, type=float, help='Target weights')
def create(name: str, structure: str, elements: tuple, weights: tuple):
    """Create a new SQS exploration"""
    if len(elements) != len(weights):
        console.print("[red]Error: Number of elements must match number of weights[/red]")
        return

    data = {
        'name': name,
        'starting_structure': structure,
        'target_concentrations': [
            {'element': elem, 'weight': weight}
            for elem, weight in zip(elements, weights)
        ]
    }

    from atomict.cli.core.client import get_client
    client = get_client()
    result = client.post('/api/catalysis-exploration/', data)    
    console.print(f"[green]Created exploration {result['id']}[/green]")
