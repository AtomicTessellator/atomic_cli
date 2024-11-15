from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="user")
def user_group():
    """Manage users and user uploads"""
    pass


@user_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get user details or list all users"""
    client = get_client()

    if id:
        result = client.get(f"/api/user/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]User Details[/bold]"))
        console.print(f"Username: {result.get('username', 'N/A')}")
        console.print(f"Email: {result.get('email', 'N/A')}")
        console.print(f"First Name: {result.get('first_name', 'N/A')}")
        console.print(f"Last Name: {result.get('last_name', 'N/A')}")
        console.print(f"Date Joined: {format_datetime(result.get('date_joined'))}")
        console.print(f"Last Login: {format_datetime(result.get('last_login'))}")
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        if fetch_all:
            results = client.get_all("/api/user/", params=params)
        else:
            results = client.get("/api/user/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("Username", "username", None),
            ("Email", "email", None),
            ("First Name", "first_name", None),
            ("Last Name", "last_name", None),
            ("OpenAPI Key", "openapi_key", None),
            ("Date Joined", "date_joined", format_datetime),
            ("Last Login", "last_login", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No users found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="Users", caption=footer_string
        )

        console.print(table)


@user_group.command()
@click.argument("id", required=False)
@click.option("--user", help="User ID to filter by")
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_upload(
    id: Optional[str] = None,
    user: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get user upload details or list all uploads

    When listing uploads, --user can be provided to filter by user.
    """
    client = get_client()

    if id:
        result = client.get(f"/api/user-upload/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel.fit(f"[bold]User Upload Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"File Name: {result.get('orig_name', 'N/A')}")
        console.print(f"Size: {result.get('size', 'N/A')} bytes")
        console.print(f"Type: {result.get('type', 'N/A')}")
        console.print(f"Description: {result.get('users_description', 'N/A')}")
        console.print(f"Uploaded: {format_datetime(result.get('uploaded'))}")
    else:
        params = {}
        if search:
            params["search"] = search
        if user:
            params["user"] = user
        if ordering:
            params["ordering"] = ordering

        if fetch_all:
            results = client.get_all("/api/user-upload/", params=params)
        else:
            results = client.get("/api/user-upload/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("File Name", "orig_name", None),
            ("Size", "size", lambda x: f"{x} bytes"),
            ("Type", "type", None),
            ("Description", "users_description", None),
            ("Uploaded", "uploaded", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No uploads found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="User Uploads", caption=footer_string
        )

        console.print(table)


@user_group.command()
@click.argument("id")
@click.option(
    "--output",
    "-o",
    help="Output file path (default: current directory with original filename)",
)
def download(id: str, output: Optional[str] = None):
    """Download a user upload file"""
    client = get_client()

    try:
        with console.status(f"[bold blue]Downloading file..."):
            result = client.get(
                f"/api/user-upload/{id}/", params={"include_content": "true"}
            )

        filename = result.get("orig_name", f"file-{id}")
        if not output:
            output = filename

        # The file content is base64 encoded in the response
        file_content = result.get("file_content_base64")
        if not file_content:
            console.print("[red]No file content received[/red]")
            return

        import base64

        with open(output, "wb") as f:
            f.write(base64.b64decode(file_content))

        console.print(f"[green]Successfully downloaded to: {output}[/green]")
    except Exception as e:
        click.echo(f"[red]Error downloading file: {str(e)}[/red]", err=True)
