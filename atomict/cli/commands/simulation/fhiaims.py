import click
from rich.console import Console
from rich.table import Table
from rich.box import ASCII_DOUBLE_HEAD
from typing import Optional
import json

from atomict.cli.commands.helpers import format_datetime
from atomict.cli.core.client import get_client


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
@click.option('--json-output', is_flag=True, help='Output in JSON format')
@click.option('--all', 'fetch_all', is_flag=True, help='Fetch all results')
def get(id: Optional[str] = None, search: Optional[str] = None, 
        ordering: Optional[str] = None, json_output: bool = False,
        fetch_all: bool = False):
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
        console.print(f"Description: {simulation.get('description', 'N/A')}")
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

        # Handle pagination
        if isinstance(results, dict):
            items = results.get('results', [])
            count = results.get('count')
            page_size = len(items)
            if count and page_size:
                console = Console()
                # console.print(f"\nShowing page 1 of {(count - 1) // page_size + 1}")
                # console.print(f"Total items: {count}")
                # if results.get('next'):
                #     console.print("Use --all to fetch all results")
        else:
            items = results

        if isinstance(results, dict):
            # paginated
            footer_string = f"Showing page 1 of {(count - 1) // page_size + 1}"
            footer_string += ", " + f"Total items: {count}"
            footer_string += ". " + "Use --all to fetch all results"
        else:
            footer_string = ""

        # Display table
        console = Console()
        table = Table(
            title="FHI-aims Simulations",
            title_style="bold",
            title_justify="center",
            box=ASCII_DOUBLE_HEAD,
            caption=footer_string,
            caption_justify="center",
            show_header=True,
            header_style="bold cyan",
            # show_footer=True,
            # footer_style="bold cyan",
            highlight=True,
        )
        table.add_column("ID")
        table.add_column("Name")
        table.add_column("Created")
        table.add_column("Status")
        
        for item in items:
            status_code = item.get('task', {}).get('status')
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
