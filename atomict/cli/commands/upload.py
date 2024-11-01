# cli/commands/upload.py
import json
from typing import Optional

import click
from rich.progress import Progress
from rich.console import Console
from rich.table import Table
from pathlib import Path

console = Console()


@click.group(name='upload')
def upload():
    """Manage file uploads"""
    pass


@upload.command()
@click.argument('file_path', type=click.Path(exists=True))
@click.option('--type', 'file_type', help='Type of file being uploaded')
@click.option('--description', help='Description of the file')
def create(file_path: str, file_type: Optional[str], description: Optional[str]):
    """Upload a new file"""
    from atomict.cli.core.client import get_client
    
    client = get_client()
    path = Path(file_path)
    
    with Progress() as progress:
        task = progress.add_task(f"Uploading {path.name}...", total=path.stat().st_size)
        
        # Implementation would depend on your API's upload endpoint
        # This is a placeholder for the upload logic
        with open(path, 'rb') as f:
            data = {
                'file': f,
                'type': file_type,
                'description': description
            }
            result = client.post('/api/user-upload/', data)
            
        progress.update(task, completed=path.stat().st_size)
    
    console.print(f"[green]Successfully uploaded {path.name}[/green]")
    console.print(f"Upload ID: {result['id']}")

@upload.command()
@click.option('--limit', type=int, help='Number of results to return')
@click.option('--type', 'file_type', help='Filter by file type')
@click.option('--json-output', is_flag=True, help='Output in JSON format')
def list(limit: Optional[int], file_type: Optional[str], json_output: bool):
    """List uploaded files"""
    from atomict.cli.core.client import get_client
    
    params = {}
    if limit:
        params['limit'] = limit
    if file_type:
        params['type'] = file_type

    client = get_client()
    results = client.get_all('/api/user-upload/', params)

    if json_output:
        click.echo(json.dumps(results, indent=2))
        return

    table = Table(show_header=True)
    table.add_column("ID")
    table.add_column("Filename")
    table.add_column("Type")
    table.add_column("Size")
    table.add_column("Uploaded")

    for upload in results:
        table.add_row(
            str(upload['id']),
            upload.get('filename', 'N/A'),
            upload.get('type', 'N/A'),
            str(upload.get('size', 'N/A')),
            upload.get('created_at', 'N/A')
        )

    console.print(table)
