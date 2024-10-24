import click
from rich.console import Console
from rich.table import Table
from typing import Optional, List

from atomict.cli.core.client import get_client
from atomict.cli.commands.helpers import get_status_string


@click.group(name='kpoint')
def kpoint_group():
    """Manage K-point simulations"""
    pass


@kpoint_group.command()
@click.argument('id', required=False)
@click.option('--search', help='Search term')
@click.option('--ordering', help='Field to order results by')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get(id: Optional[str] = None, search: Optional[str] = None,
        ordering: Optional[str] = None, json_output: bool = False):
    """Get K-point simulation details or list all simulations"""
    client = get_client()
    
    if id:
        simulation = client.get(f'/api/kpoint-simulation/{id}/')
        if json_output:
            click.echo(simulation)
        else:
            console = Console()
            console.print(f"ID: {simulation['id']}")
            console.print(f"Exploration: {simulation['exploration']}")
            console.print(f"K-points: {simulation.get('k_points', [])}")
            if simulation.get('simulation'):
                status = get_status_string(simulation['simulation'].get('status'))
                console.print(f"Status: {status}")
    else:
        params = {}
        if search:
            params['search'] = search
        if ordering:
            params['ordering'] = ordering
            
        simulations = client.get('/api/kpoint-simulation/', params=params)
        
        if json_output:
            click.echo(simulations)
        else:
            console = Console()
            table = Table(show_header=True)
            table.add_column("ID")
            table.add_column("Exploration")
            table.add_column("K-points")
            table.add_column("Status")
            
            for sim in simulations:
                status_code = sim.get('simulation', {}).get('status')
                status = get_status_string(status_code)
                table.add_row(
                    str(sim['id']),
                    str(sim['exploration']),
                    str(sim.get('k_points', [])),
                    status
                )
            console.print(table)


@kpoint_group.command()
@click.option('--exploration', required=True, help='Exploration ID')
@click.option('--k-points', required=True, multiple=True, type=float, help='K-point values')
def create(exploration: str, k_points: List[float]):
    """Create a new K-point simulation"""
    client = get_client()
    data = {
        'exploration': exploration,
        'k_points': list(k_points)
    }
        
    simulation = client.post('/api/kpoint-simulation/', json=data)
    click.echo(f"Created simulation with ID: {simulation['id']}")


@kpoint_group.command()
@click.argument('id')
def delete(id: str):
    """Delete a K-point simulation"""
    client = get_client()
    client.delete(f'/api/kpoint-simulation/{id}/')
    click.echo(f"Deleted simulation {id}")
