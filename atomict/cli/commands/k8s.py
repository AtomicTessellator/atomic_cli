# cli/commands/k8s.py
import json

import click
from rich.table import Table
from rich.console import Console
from rich.panel import Panel


console = Console()


@click.group(name='k8s')
def k8s():
    """Manage Kubernetes jobs and clusters"""
    pass

@k8s.command()
@click.argument('id')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def get_job(id: str, json_output: bool):
    """Get K8s job details including logs"""
    from atomict.cli.core.client import get_client
    
    client = get_client()
    result = client.get(f'/api/k8s-job/{id}/')

    if json_output:
        click.echo(json.dumps(result, indent=2))
        return

    # Main job info
    console.print(Panel(f"[bold]K8s Job Details[/bold]"))
    console.print(f"ID: {result['id']}")
    console.print(f"Name: {result.get('name', 'N/A')}")
    console.print(f"Status: {result.get('status', 'N/A')}")
    console.print()

    # Resource requests/limits
    resource_table = Table(show_header=True, title="Resource Allocation")
    resource_table.add_column("Type")
    resource_table.add_column("Request")
    resource_table.add_column("Limit")
    
    resource_table.add_row(
        "CPU",
        result.get('cpu_request', 'N/A'),
        result.get('cpu_limit', 'N/A')
    )
    resource_table.add_row(
        "Memory",
        result.get('memory_request', 'N/A'),
        result.get('memory_limit', 'N/A')
    )
    console.print(resource_table)
    console.print()

    # Pod events if present
    if 'events' in result:
        console.print("[bold]Pod Events[/bold]")
        events_table = Table(show_header=True)
        events_table.add_column("Time")
        events_table.add_column("Type")
        events_table.add_column("Message")
        
        for event in result['events']:
            events_table.add_row(
                event.get('time', 'N/A'),
                event.get('type', 'N/A'),
                event.get('message', 'N/A')
            )
        console.print(events_table)

@k8s.command()
@click.argument('id')
def logs(id: str):
    """Stream logs from a K8s job"""
    raise NotImplementedError("Not implemented")
    from atomict.cli.core.client import get_client
    
    client = get_client()
    
    try:
        # TODO:
        for log_line in client.stream(f'/api/k8s-job/{id}/logs/'):
            console.print(log_line.get('message', ''))
    except KeyboardInterrupt:
        console.print("\n[yellow]Stopped log streaming[/yellow]")
