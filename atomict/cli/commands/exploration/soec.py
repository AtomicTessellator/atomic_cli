import click
from rich.console import Console
from rich.table import Table
from typing import Optional, List

from atomict.cli.core.client import get_client
from atomict.cli.commands.helpers import format_datetime, get_status_string


@click.group(name='soecexploration')
def soecexploration_group():
    """Manage SOEC explorations"""
    pass


@soecexploration_group.command()
@click.argument('id', required=False)
@click.option('--search', help='Search term')
@click.option('--ordering', help='Field to order results by')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
@click.option('--limit', type=int, help='Number of results per page')
@click.option('--offset', type=int, help='Results offset')
def get(id: Optional[str] = None, search: Optional[str] = None,
        ordering: Optional[str] = None, json_output: bool = False,
        limit: Optional[int] = None, offset: Optional[int] = None):
    """Get SOEC exploration details or list all explorations"""
    client = get_client()
    
    if id:
        exploration = client.get(f'/api/ea-exploration/{id}/')
        if json_output:
            click.echo(exploration)
        else:
            console = Console()
            console.print(f"ID: {exploration['id']}")
            console.print(f"Name: {exploration.get('name', 'N/A')}")
            console.print(f"Description: {exploration.get('description', 'N/A')}")
            console.print(f"Created: {format_datetime(exploration['created_at'])}")
            console.print(f"Edited: {format_datetime(exploration.get('edited_at', 'N/A'))}")
            console.print(f"Stress Algorithm: {exploration.get('stress_algorithm', 'N/A')}")
            console.print(f"Stress Method: {exploration.get('stress_method', 'N/A')}")
            console.print(f"Number of Last Samples: {exploration.get('num_last_samples', 'N/A')}")
            console.print(f"Strains List: {exploration.get('strains_list', [])}")
            if exploration.get('task'):
                status = get_status_string(exploration['task'].get('status'))
                console.print(f"Status: {status}")
    else:
        params = {}
        if search:
            params['search'] = search
        if ordering:
            params['ordering'] = ordering
        if limit:
            params['limit'] = limit
        if offset:
            params['offset'] = offset
            
        response = client.get('/api/ea-exploration/', params=params)
        
        if json_output:
            click.echo(response)
        else:
            console = Console()
            table = Table(show_header=True)
            table.add_column("ID")
            table.add_column("Name")
            table.add_column("Created")
            table.add_column("Stress Method")
            table.add_column("Status")
            
            for exp in response['results']:
                status_code = exp.get('task', {}).get('status')
                status = get_status_string(status_code)
                table.add_row(
                    str(exp['id']),
                    exp.get('name', 'N/A'),
                    format_datetime(exp['created_at']),
                    exp.get('stress_method', 'N/A'),
                    status
                )
            console.print(table)


@soecexploration_group.command()
@click.option('--name', help='Exploration name')
@click.option('--description', help='Exploration description')
@click.option('--project', required=True, help='Project ID')
@click.option('--stress-algorithm', help='Stress algorithm to use')
@click.option('--stress-method', help='Stress method to use')
@click.option('--num-last-samples', type=int, help='Number of last samples')
@click.option('--strains', multiple=True, help='List of strains')
@click.option('--relaxed-structure', help='Relaxed structure ID')
def create(name: Optional[str], description: Optional[str], project: str,
           stress_algorithm: Optional[str] = None, stress_method: Optional[str] = None,
           num_last_samples: Optional[int] = None, strains: Optional[List[str]] = None,
           relaxed_structure: Optional[str] = None):
    """Create a new SOEC exploration"""
    client = get_client()
    data = {
        'project': project
    }
    if name:
        data['name'] = name
    if description:
        data['description'] = description
    if stress_algorithm:
        data['stress_algorithm'] = stress_algorithm
    if stress_method:
        data['stress_method'] = stress_method
    if num_last_samples is not None:
        data['num_last_samples'] = num_last_samples
    if strains:
        data['strains_list'] = list(strains)
    if relaxed_structure:
        data['relaxed_structure'] = relaxed_structure

    # the API docs say this takes a FHIAims simulation object
    # the required data will have to be ironed out on the server side
    # then this can be revisited
    exploration = client.post('/api/ea-exploration/', json=data)
    click.echo(f"Created exploration with ID: {exploration['id']}")


@soecexploration_group.command()
@click.argument('id')
def delete(id: str):
    """Delete a SOEC exploration"""
    client = get_client()
    client.delete(f'/api/ea-exploration/{id}/')
    click.echo(f"Deleted exploration {id}")
