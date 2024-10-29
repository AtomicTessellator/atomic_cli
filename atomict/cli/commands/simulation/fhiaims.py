import click
from rich.console import Console
from rich.table import Table
from rich.box import ASCII_DOUBLE_HEAD
from typing import Optional
import json

from atomict.cli.commands.helpers import format_datetime
from atomict.cli.core.client import get_client
from atomict.cli.commands.common import table_0
from atomict.cli.core.utils import get_pagination_info


def get_status_string(status_code: Optional[int]) -> str:
    """Convert status code to human readable string"""
    status_map = {
        0: "Draft",
        1: "Ready",
        2: "Running",
        3: "Completed",
        4: "Error",
        5: "Paused",
        6: "User Aborted"
    }
    if status_code is None:
        return "N/A"
    return status_map.get(status_code, "Unknown")


@click.group(name='fhiaims')
def fhiaims_group():
    """Manage FHI-aims simulations"""
    pass


@fhiaims_group.command()
@click.argument('id', required=False)
@click.option('--search', help='Search term')
@click.option('--ordering', help='Field to order results by')
# filter should be enum to match DRF config
@click.option('--filter', 'filters', multiple=True, help='Filter in format field=value')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
@click.option('--all', 'fetch_all', is_flag=True, help='Fetch all results')
def get(id: Optional[str] = None, search: Optional[str] = None, 
        ordering: Optional[str] = None, filters: tuple = (), 
        json_output: bool = False, fetch_all: bool = False):
    """Get simulation details or list all simulations"""
    client = get_client()
    
    if id:
        simulation = client.get(f'/api/fhiaims-simulation/{id}/')
        if json_output:
            click.echo(json.dumps(simulation, indent=2))
            return
        # Format single simulation output
        console = Console()
        console.print(f"ID: {simulation['id']}")
        console.print(f"Name: {simulation.get('name', 'N/A')}")
        console.print(f"Created: {format_datetime(simulation['created_at'])}")
        if simulation.get('task'):
            status = get_status_string(simulation['task'].get('status'))
            console.print(f"Status: {status}")
    else:
        params = {}
        if search:
            params['search'] = search
        if ordering:
            params['ordering'] = ordering
        
        # Add filter parameters
        for f in filters:
            try:
                field, value = f.split('=', 1)
                params[field] = value
            except ValueError:
                click.echo(f"Invalid filter format: {f}. Use field=value", err=True)
                return
            
        # Updated this section
        if fetch_all:
            # list
            results = client.get_all('/api/fhiaims-simulation/', params=params)
        else:
            # dict
            results = client.get('/api/fhiaims-simulation/', params=params)
        
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        # Get the data from the results section and a string with page details
        items, footer_string = get_pagination_info(results)
        
        # Display table
        console = Console()

        if not items:
            console.print(f"[white]No simulations found with the given criteria:[/white]\n[green]{params}")
            return

        table = table_0
        table.title = "FHI-aims Simulations"
        table.caption = footer_string

        # TODO: discuss presentation. what to show?
        table.add_column("ID")
        table.add_column("Name")
        table.add_column("Created")
        table.add_column("Status")

        for item in items:
            status_code = item.get('task', {}).get('status')  # nested task
            status = get_status_string(status_code)
            table.add_row(
                str(item['id']),
                item.get('name', 'N/A'),
                format_datetime(item['created_at']),
                status
            )
        console.print(table)


@fhiaims_group.command()
@click.option('--name', help='Simulation name')
@click.option('--description', help='Simulation description')
@click.option('--control-file', required=True, help='Control file content')
@click.option('--geometry-file', required=True, help='Geometry file content')
@click.option('--generate-finite-diff', is_flag=True, help='Generate finite difference displacements')
@click.option('--finite-diff-displacement', type=float, help='Finite difference displacement value')
def create(name: Optional[str], description: Optional[str],
          control_file: str, geometry_file: str,
          generate_finite_diff: bool = False,
          finite_diff_displacement: Optional[float] = None):
    """Create a new FHI-aims simulation"""
    client = get_client()
    data = {
        'control_file': control_file,
        'geometry_file': geometry_file,
        'generate_finite_diff_displacements': generate_finite_diff,
    }
    if name:
        data['name'] = name
    if description:
        data['description'] = description
    if finite_diff_displacement is not None:
        data['finite_diff_displacement'] = finite_diff_displacement
        
    simulation = client.post('/api/fhiaims-simulation/', data=data)
    click.echo(f"Created simulation with ID: {simulation['id']}")


@fhiaims_group.command()
@click.argument('id')
def delete(id: str):
    """Delete a FHI-aims simulation"""
    client = get_client()
    client.delete(f'/api/fhiaims-simulation/{id}/')
    click.echo(f"Deleted simulation {id}")
