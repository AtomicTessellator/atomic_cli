import json

import click
from rich.console import Console
from rich.table import Table
from typing import Optional, List

from atomict.cli.core.client import get_client
from atomict.cli.commands.helpers import format_datetime, get_status_string


@click.group(name='sqs')
def sqs_group():
    """Manage Special Quasirandom Structure (SQS) explorations"""
    pass


@sqs_group.command()
@click.argument('id', required=False)
@click.option('--search', help='Search term')
@click.option('--ordering', help='Field to order results by')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get(id: Optional[str] = None, search: Optional[str] = None,
        ordering: Optional[str] = None, json_output: bool = False):
    """Get SQS exploration details or list all explorations"""
    client = get_client()
    
    if id:
        exploration = client.get(f'/api/sqs-exploration/{id}/')
        if json_output:
            click.echo(exploration)
        else:
            console = Console()
            console.print(f"ID: {exploration['id']}")
            console.print(f"Name: {exploration.get('name', 'N/A')}")
            console.print(f"Description: {exploration.get('description', 'N/A')}")
            console.print(f"Created: {format_datetime(exploration['created_at'])}")
            console.print(f"Updated: {format_datetime(exploration.get('updated_at', 'N/A'))}")
            console.print(f"Max Size: {'Auto' if exploration.get('auto_max_size') else exploration.get('max_size', 'N/A')}")
            console.print(f"Atom Count Limit: {exploration.get('atom_count_upper_limit', 'N/A')}")
            console.print(f"Cluster Cutoffs: {exploration.get('cluster_cutoffs', [])}")
            if exploration.get('task'):
                status = get_status_string(exploration['task'].get('status'))
                console.print(f"Status: {status}")
                console.print(f"Running on: {exploration['task'].get('running_on', 'N/A')}")
    else:
        params = {}
        if search:
            params['search'] = search
        if ordering:
            params['ordering'] = ordering

        response = client.get('/api/sqs-exploration/', params=params)
        if 'results' in response:
            # TBD: standardize paginated responses
            response = response['results']

        if json_output:
            click.echo(json.dumps(response, indent=2))
        else:
            console = Console()
            table = Table(show_header=True)
            table.add_column("ID")
            table.add_column("Name")
            table.add_column("Max Size")
            table.add_column("Created")
            table.add_column("Status")
            print(response)
            for exp in response:
                status_code = exp.get('task', {}).get('status')
                status = get_status_string(status_code)
                table.add_row(
                    exp['id'],
                    exp.get('name', 'N/A'),
                    'Auto' if exp.get('auto_max_size') else str(exp.get('max_size', 'N/A')),
                    format_datetime(exp['created_at']),
                    status
                )
            console.print(table)


@sqs_group.command()
@click.option('--name', help='Exploration name')
@click.option('--description', help='Exploration description')
@click.option('--project', required=True, help='Project ID')
@click.option('--auto-max-size', is_flag=True, help='Automatically determine max size')
@click.option('--max-size', type=int, help='Maximum structure size')
@click.option('--atom-count-limit', type=int, help='Upper limit for atom count')
@click.option('--cutoffs', type=float, multiple=True, help='Cluster cutoff values')
@click.option('--starting-structure', help='Starting structure ID')
def create(name: Optional[str], description: Optional[str], project: str,
           auto_max_size: bool = False, max_size: Optional[int] = None,
           atom_count_limit: Optional[int] = None,
           cutoffs: Optional[List[float]] = None,
           starting_structure: Optional[str] = None):
    """Create a new SQS exploration"""
    client = get_client()
    data = {
        'project': project,
        'auto_max_size': auto_max_size
    }
    if name:
        data['name'] = name
    if description:
        data['description'] = description
    if max_size is not None:
        data['max_size'] = max_size
    if atom_count_limit is not None:
        data['atom_count_upper_limit'] = atom_count_limit
    if cutoffs:
        data['cluster_cutoffs'] = list(cutoffs)
    if starting_structure:
        data['starting_structure'] = starting_structure
        
    exploration = client.post('/api/sqs-exploration/', json=data)
    click.echo(f"Created exploration with ID: {exploration['id']}")


@sqs_group.command()
@click.argument('id')
def delete(id: str):
    """Delete an SQS exploration"""
    client = get_client()
    client.delete(f'/api/sqs-exploration/{id}/')
    click.echo(f"Deleted exploration {id}")
