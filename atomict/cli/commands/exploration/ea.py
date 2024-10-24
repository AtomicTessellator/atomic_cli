import click
from rich.console import Console
from rich.table import Table
from typing import Optional

from atomict.cli.core.client import get_client
from atomict.cli.commands.helpers import format_datetime


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

@click.group(name='ea')
def ea_group():
    """Manage Electrolyzer Analysis (EA) explorations"""
    pass

@ea_group.command()
@click.argument('id', required=False)
@click.option('--search', help='Search term')
@click.option('--ordering', help='Field to order results by')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
@click.option('--limit', type=int, help='Number of results per page')
@click.option('--offset', type=int, help='Results offset')
def get(id: Optional[str] = None, search: Optional[str] = None,
        ordering: Optional[str] = None, json_output: bool = False,
        limit: Optional[int] = None, offset: Optional[int] = None):
    """Get EA exploration details or list all explorations"""
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
            console.print(f"Updated: {format_datetime(exploration.get('updated_at', 'N/A'))}")
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
            table.add_column("Updated")
            table.add_column("Status")
            
            for exp in response['results']:
                status_code = exp.get('task', {}).get('status')
                status = get_status_string(status_code)
                table.add_row(
                    str(exp['id']),
                    exp.get('name', 'N/A'),
                    format_datetime(exp['created_at']),
                    format_datetime(exp.get('updated_at', 'N/A')),
                    status
                )
            console.print(table)


@ea_group.command()
@click.option('--name', help='Exploration name')
@click.option('--description', help='Exploration description')
@click.option('--project', required=False, help='Project ID')
def create(name: Optional[str], description: Optional[str], project: str):
    """Create a new EA exploration"""
    client = get_client()
    data = {}
    if project:
        # TODO: backend work on this. Project ID validation is done but the error doesn't
        # get back to the client
        # also IDs should be created serverside and returned on creation
        # ideally, the required fields are marked in the DRF serializers
        data['project'] = project
    if name:
        data['name'] = name
    if description:
        data['description'] = description
        
    exploration = client.post('/api/ea-exploration/', data=data)
    click.echo(f"Created exploration with ID: {exploration['id']}")


@ea_group.command()
@click.argument('id')
def delete(id: str):
    """Delete an EA exploration"""
    client = get_client()
    client.delete(f'/api/ea-exploration/{id}/')
    click.echo(f"Deleted exploration {id}")
