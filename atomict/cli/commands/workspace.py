# cli/commands/workspace.py
import base64
import os
from pathlib import Path
from typing import Optional

import click
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime
from atomict.cli.core.client import get_client

console = Console()


@click.group(name="workspace")
def workspace_group():
    """Manage workspace simulations and files"""
    pass


@workspace_group.command()
@click.argument('uuid', type=str)
@click.argument('path', type=str, default=".", required=False)
def cp(uuid: str, path: str = "."):
    """Download simulation files by UUID to a directory
    
    Example:
        tess workspace cp 8f328887-62ae-4ace-b37b-2c31640baa47 .
    """
    client = get_client()
    
    # Fetch simulation metadata
    try:
        result = client.get(f"/simulation/lookup/{uuid}/")
    except Exception as e:
        console.print(f"[red]Error fetching simulation: {e}[/red]")
        return
    
    sim_type = result.get("type", "Unknown")
    simulation = result.get("simulation", {})
    files = result.get("files", [])
    
    # Display simulation info
    console.print(f"\n[bold cyan]Simulation Type:[/bold cyan] {sim_type}")
    console.print(f"[bold cyan]ID:[/bold cyan] {simulation.get('id', 'N/A')}")
    console.print(f"[bold cyan]Name:[/bold cyan] {simulation.get('name') or '(unnamed)'}")
    
    if not files:
        console.print("\n[yellow]No files found for this simulation[/yellow]")
        return
    
    # Create output directory
    output_path = Path(path)
    output_path.mkdir(parents=True, exist_ok=True)
    
    console.print(f"\n[bold green]Downloading {len(files)} file(s) to: {output_path.absolute()}[/bold green]\n")
    
    # Download files with progress bar
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:
        
        download_task = progress.add_task("[cyan]Downloading files...", total=len(files))
        download_ct = 0
        for file_info in files:
            # Extract user_upload info (handle nested structure)
            user_upload = file_info.get('user_upload', {})
            if isinstance(user_upload, dict):
                file_id = user_upload.get('id')
                filename = user_upload.get('orig_name', f"file-{file_id}")
            else:
                # If user_upload is an ID string
                file_id = user_upload
                filename = f"file-{file_id}"
            
            if not file_id:
                console.print(f"[yellow]Skipping file with no ID[/yellow]")
                progress.advance(download_task)
                continue
            
            progress.update(download_task, description=f"[cyan]Downloading: {filename}")
            
            try:
                # Fetch file with content
                file_data = client.get(
                    f"/api/user-upload/{file_id}/",
                    params={"include_content": "true"}
                )
                
                file_content_b64 = file_data.get("file_content_base64")
                if not file_content_b64:
                    console.print(f"[red]✗ No content received for: {filename}[/red]")
                    progress.advance(download_task)
                    continue
                
                # Decode and save file
                file_path = output_path / filename
                with open(file_path, "wb") as f:
                    f.write(base64.b64decode(file_content_b64))
                
                progress.advance(download_task)
                
            except Exception as e:
                console.print(f"[red]✗ Error downloading {filename}: {e}[/red]")
                progress.advance(download_task)

            download_ct += 1
    
    console.print(f"\n[bold {'green' if download_ct else 'red'}]Downloaded {download_ct} file(s)[/bold {'green' if download_ct else 'red'}]")

