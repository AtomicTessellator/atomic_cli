# cli/commands/upload.py
from typing import Optional

import click
from rich.console import Console

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="upload")
def upload_group():
    """Manage file uploads"""
    pass


# TBD: not currently supported via API
# @upload.command()
# @click.argument('file_path', type=click.Path(exists=True))
# @click.option('--type', 'file_type', help='Type of file being uploaded')
# @click.option('--description', help='Description of the file')
# def create(file_path: str, file_type: Optional[str], description: Optional[str]):
#     """Upload a new file"""
#     from atomict.cli.core.client import get_client

#     client = get_client()
#     path = Path(file_path)

#     with Progress() as progress:
#         task = progress.add_task(f"Uploading {path.name}...", total=path.stat().st_size)

#         # Implementation would depend on your API's upload endpoint
#         # This is a placeholder for the upload logic
#         with open(path, 'rb') as f:
#             data = {
#                 'file': f,
#                 'type': file_type,
#                 'description': description
#             }
#             result = client.post('/api/user-upload/', data)

#         progress.update(task, completed=path.stat().st_size)

#     console.print(f"[green]Successfully uploaded {path.name}[/green]")
#     console.print(f"Upload ID: {result['id']}")


@upload_group.command()
@click.option("--limit", type=int, help="Number of results to return")
@click.option("--type", "file_type", help="Filter by file type")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get(limit: Optional[int], file_type: Optional[str], json_output: bool = False, fetch_all: bool = False):
    """Get uploaded files"""
    params = {}
    if limit:
        params["limit"] = limit
    if file_type:
        params["type"] = file_type

    client = get_client()
    
    if fetch_all:
        results = client.get_all("/api/user-upload/", params)
    else:
        results = client.get("/api/user-upload/", params)

    if json_output:
        console.print_json(data=results)
        return

    columns = [
        ("ID", "id", None),
        ("Filename", "orig_name", None),
        ("Type", "type", None),
        ("Size", "size", None),
        ("Uploaded", "uploaded", format_datetime),
    ]
    items, footer_string = get_pagination_info(results)

    table = create_table(
        columns=columns, items=items, title="Uploads", caption=footer_string
    )

    console.print(table)
