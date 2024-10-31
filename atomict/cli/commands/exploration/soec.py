import json
import click
from rich.console import Console
from rich.panel import Panel
from typing import Optional

from atomict.cli.core.client import get_client
from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.utils import get_pagination_info


@click.group(name='ea')
def soecexploration_group():
    """Manage EA explorations and related resources"""
    pass


@soecexploration_group.command()
@click.argument('id', required=False)
@click.option('--search', help='Search term')
@click.option('--ordering', help='Field to order results by')
@click.option('--filter', 'filters', multiple=True, help='Filter in format field=value')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
@click.option('--all', 'fetch_all', is_flag=True, help='Fetch all results')
def get(id: Optional[str] = None, search: Optional[str] = None,
        ordering: Optional[str] = None, filters: tuple = (), 
        json_output: bool = False, fetch_all: bool = False):
    """Get EA exploration details or list all explorations"""
    client = get_client()
    console = Console()
    # WIP
    if id:
        result = client.get(f'/api/ea-exploration/{id}/')
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel(f"[bold]EA Exploration Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")
        console.print(f"Project: {result['project'].get('name', 'N/A')} ({result['project']['id']})")
        console.print(f"Status: {get_status_string(result.get('status'))}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at'))}")
        
        if result.get('parameters'):
            console.print("\n[bold]Parameters[/bold]")
            console.print_json(data=result['parameters'])
    else:
        params = {}
        if search is not None:
            params['search'] = search
        if ordering:
            params['ordering'] = ordering
        
        for f in filters:
            try:
                field, value = f.split('=', 1)
                params[field] = value
            except ValueError:
                console.print(f"[red]Invalid filter format: {f}. Use field=value[/red]")
                return

        if fetch_all:
            results = client.get_all('/api/ea-exploration/', params=params)
        else:
            results = client.get('/api/ea-exploration/', params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Project", "project", lambda x: x.get('name') if isinstance(x, dict) else None),
            ("Status", "status", get_status_string),
            ("Created", "created_at", format_datetime),
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(f"[white]No EA explorations found with the given criteria:[/white]\n[green]{params}")
            return

        table = create_table(
            columns=columns,
            items=items,
            title="EA Explorations",
            caption=footer_string
        )
        
        console.print(table)


@soecexploration_group.command()
@click.option('--name', required=True, help='Exploration name')
@click.option('--project', required=True, help='Project ID')
@click.option('--parameters', required=True, help='Exploration parameters (JSON string)')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def create(name: str, project: str, parameters: str, json_output: bool = False):
    """Create a new EA exploration"""
    client = get_client()
    console = Console()
    #WIP
    try:
        params_data = json.loads(parameters)
    except json.JSONDecodeError:
        console.print("[red]Invalid JSON format for parameters[/red]")
        return
    
    data = {
        'name': name,
        'project': project,
        'parameters': params_data
    }
        
    result = client.post('/api/ea-exploration/', data=data)
    
    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"[green]Created EA exploration '{name}' with ID: {result['id']}[/green]")


@soecexploration_group.command()
@click.argument('id')
def delete(id: str):
    """Delete an EA exploration"""
    #WIP
    client = get_client()
    client.delete(f'/api/ea-exploration/{id}/')
    console = Console()
    console.print(f"[green]Deleted EA exploration {id}[/green]") 