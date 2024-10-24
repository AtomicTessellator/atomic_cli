# cli/commands/project.py
import click
from rich.table import Table
from rich.console import Console
from rich.panel import Panel
from rich.tree import Tree
import json
from typing import Optional

from atomict.cli.commands.helpers import format_datetime, get_status_string


console = Console()

@click.group(name='project')
def project():
    """Manage projects and their related resources"""
    pass

@project.command()
@click.argument('id', required=False)
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get(id: Optional[str] = None, json_output: bool = False):
    """
    Get project details. If no ID is provided, lists all projects.
    """
    from atomict.cli.core.client import get_client
    
    client = get_client()

    if id:
        # Get single project
        result = client.get(f'/api/project/{id}/')
        if json_output:
            click.echo(json.dumps(result, indent=2))
            return

        # Single project detailed view
        console.print(Panel(f"[bold]Project: {result.get('name', 'N/A')}[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Created: {format_datetime(result.get('created_at', 'N/A'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at', 'N/A'))}")
        console.print()

        # Related simulations table
        if 'simulations' in result:
            console.print("[bold]Related Simulations[/bold]")
            sim_table = Table(show_header=True)
            sim_table.add_column("ID")
            sim_table.add_column("Name")
            sim_table.add_column("Status")
            
            for sim in result['simulations']:
                sim_table.add_row(
                    str(sim['id']),
                    sim.get('name', 'N/A'),
                    get_status_string(sim.get('status', 'N/A'))
                )
            console.print(sim_table)
            console.print()

        # Project notes table
        if 'notes' in result:
            console.print("[bold]Project Notes[/bold]")
            notes_table = Table(show_header=True)
            notes_table.add_column("ID")
            notes_table.add_column("Title")
            notes_table.add_column("Created")
            
            for note in result['notes']:
                notes_table.add_row(
                    str(note['id']),
                    note.get('title', 'N/A'),
                    format_datetime(note.get('created_at', 'N/A'))
                )
            console.print(notes_table)
            console.print()

        # Show tags as a tree
        if 'tags' in result:
            tag_tree = Tree("[bold]Project Tags[/bold]")
            for tag in result['tags']:
                tag_tree.add(f"[{tag.get('color', 'white')}]{tag.get('tag', 'N/A')}[/{tag.get('color', 'white')}]")
            console.print(tag_tree)

    else:
        # List all projects
        results = client.get_all('/api/project/')
        if json_output:
            click.echo(json.dumps(results, indent=2))
            return

        # Projects table view
        table = Table(show_header=True)
        table.add_column("ID")
        table.add_column("Name")
        table.add_column("Created")
        table.add_column("Updated")

        for project in results:
            table.add_row(
                str(project['id']),
                project.get('name', 'N/A'),
                format_datetime(project.get('created_at', 'N/A')),
                format_datetime(project.get('updated_at', 'N/A'))
            )

        console.print(table)
