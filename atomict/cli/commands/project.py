# cli/commands/project.py
from typing import Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.tree import Tree

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info

console = Console()


@click.group(name="project")
def project_group():
    """Manage projects and their related resources"""
    pass


@project_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get project details or list all projects"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/project/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        # Single project detailed view
        console.print(Panel(f"[bold]Project Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")
        console.print(f"Smiles: {result.get('thumbnail_smiles', 'N/A')}")
        console.print(f"Created: {format_datetime(result.get('created_at', 'N/A'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at', 'N/A'))}")
        if result.get("description_html"):
            console.print("\n[bold]Description[/bold]")
            console.print(result["description_html"])
        console.print()

        # Related simulations table
        if result.get("simulations"):
            console.print("[bold]Related Simulations[/bold]")
            sim_columns = [
                ("ID", "id", None),
                ("Name", "name", None),
                ("Status", "status", get_status_string),
            ]
            sim_table = create_table(
                columns=sim_columns, items=result["simulations"], title="Simulations"
            )
            console.print(sim_table)
            console.print()

        # Project notes table
        if result.get("notes"):
            console.print("[bold]Project Notes[/bold]")
            notes_columns = [
                ("ID", "id", None),
                ("Title", "title", None),
                ("Created", "created_at", format_datetime),
            ]
            notes_table = create_table(
                columns=notes_columns, items=result["notes"], title="Notes"
            )
            console.print(notes_table)

        # Show tags as a tree
        if result.get("tags"):
            tag_tree = Tree("[bold]Project Tags[/bold]")
            for tag in result["tags"]:
                tag_tree.add(
                    f"[{tag.get('color', 'white')}]{tag.get('tag', 'N/A')}[/{tag.get('color', 'white')}]"
                )
            console.print(tag_tree)

    else:
        params = {}
        if search:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        # Add filter parameters
        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/project/", params=params)
        else:
            results = client.get("/api/project/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Smiles", "thumbnail_smiles", None),
            ("Created", "created_at", format_datetime),
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No projects found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="Projects", caption=footer_string
        )

        console.print(table)


@project_group.command()
@click.option("--name", required=True, help="Project name")
@click.option("--description", help="Project description")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create(name: str, description: Optional[str] = None, json_output: bool = False):
    """Create a new project"""
    client = get_client()
    console = Console()

    data = {
        "name": name,
    }
    if description:
        data["description"] = description

    result = client.post("/api/project/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"Created project with ID: {result['id']}")


@project_group.command()
@click.argument("id")
def delete(id: str):
    """Delete a project"""
    client = get_client()
    client.delete(f"/api/project/{id}/")
    console = Console()
    console.print(f"[green]Deleted project {id}[/green]")


@project_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_note(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get project note details or list all notes"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/project-note/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel(f"[bold]Project Note Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Title: {result.get('title', 'N/A')}")
        console.print(f"Project: {result.get('project', 'N/A')}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
        console.print(f"Updated: {format_datetime(result.get('updated_at'))}")
        if result.get("content"):
            console.print("\n[bold]Content[/bold]")
            console.print(result["content"])
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/project-note/", params=params)
        else:
            results = client.get("/api/project-note/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Title", "title", None),
            ("Project", "project", None),
            ("Created", "created_at", format_datetime),
            ("Updated", "updated_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No notes found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="Project Notes", caption=footer_string
        )

        console.print(table)


@project_group.command()
@click.option("--project", required=True, help="Project ID")
@click.option("--title", required=True, help="Note title")
@click.option("--content", help="Note content")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_note(
    project: str, title: str, content: Optional[str] = None, json_output: bool = False
):
    """Create a new project note"""
    client = get_client()
    console = Console()

    data = {
        "project": project,
        "title": title,
    }
    if content:
        data["content"] = content

    result = client.post("/api/project-note/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"[green]Created note with ID: {result['id']}[/green]")


@project_group.command()
@click.argument("id")
def delete_note(id: str):
    """Delete a project note"""
    client = get_client()
    client.delete(f"/api/project-note/{id}/")
    console = Console()
    console.print(f"[green]Deleted note {id}[/green]")


@project_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_star(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get project star details or list all starred projects"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/project-star/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel(f"[bold]Project Star Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(
            f"Project: {result['project'].get('name', 'N/A')} ({result['project']['id']})"
        )
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/project-star/", params=params)
        else:
            results = client.get("/api/project-star/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            # TODO: when serializers have a depth configured, use project.name
            ("Project", "project", None),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No starred projects found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="Starred Projects",
            caption=footer_string,
        )

        console.print(table)


@project_group.command()
@click.option("--project", required=True, help="Project ID to star")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_star(project: str, json_output: bool = False):
    """Star a project"""
    client = get_client()
    console = Console()
    data = {
        "project": project,
    }

    result = client.post("/api/project-star/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"[green]Starred project {project}[/green]")


@project_group.command()
@click.argument("id")
def delete_star(id: str):
    """Unstar a project"""
    client = get_client()
    client.delete(f"/api/project-star/{id}/")
    console = Console()
    console.print(f"[green]Unstarred project {id}[/green]")


@project_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_tag(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get project tag details or list all tags"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/project-tag/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel(f"[bold]Project Tag Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(
            f"Tag: [{result.get('color', 'white')}]{result.get('tag', 'N/A')}[/{result.get('color', 'white')}]"
        )
        console.print(f"Color: {result.get('color', 'N/A')}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/project-tag/", params=params)
        else:
            results = client.get("/api/project-tag/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Tag", "tag", None),
            ("Color", "color", None),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No tags found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="Project Tags", caption=footer_string
        )

        console.print(table)


@project_group.command()
@click.option("--tag", required=True, help="Tag name")
@click.option("--color", default="white", help="Tag color")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_tag(tag: str, color: str = "white", json_output: bool = False):
    """Create a new project tag"""
    client = get_client()
    console = Console()

    data = {
        "tag": tag,
        "color": color,
    }

    result = client.post("/api/project-tag/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created tag [{color}]{tag}[/{color}] with ID: {result['id']}[/green]"
        )


@project_group.command()
@click.argument("id")
def delete_tag(id: str):
    """Delete a project tag"""
    client = get_client()
    client.delete(f"/api/project-tag/{id}/")
    console = Console()
    console.print(f"[green]Deleted tag {id}[/green]")


@project_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_tag_project(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get project tag assignment details or list all assignments"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/project-tag-project/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel(f"[bold]Project Tag Assignment Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(
            f"Project: {result['project'].get('name', 'N/A')} ({result['project']['id']})"
        )
        console.print(
            f"Tag: [{result['tag'].get('color', 'white')}]{result['tag'].get('tag', 'N/A')}[/{result['tag'].get('color', 'white')}]"
        )
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/project-tag-project/", params=params)
        else:
            results = client.get("/api/project-tag-project/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            (
                "Project",
                "project",
                lambda x: x.get("name") if isinstance(x, dict) else None,
            ),
            ("Tag", "tag", lambda x: x.get("tag") if isinstance(x, dict) else None),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No tag assignments found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="Project Tag Assignments",
            caption=footer_string,
        )

        console.print(table)


@project_group.command()
@click.option("--project", required=True, help="Project ID")
@click.option("--tag", required=True, help="Tag ID")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_tag_project(project: str, tag: str, json_output: bool = False):
    """Assign a tag to a project"""
    client = get_client()
    console = Console()

    data = {
        "project": project,
        "project_tag": tag,
    }

    result = client.post("/api/project-tag-project/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"[green]Assigned tag {tag} to project {project}[/green]")


@project_group.command()
@click.argument("id")
def delete_tag_project(id: str):
    """Remove a tag assignment from a project"""
    client = get_client()
    client.delete(f"/api/project-tag-project/{id}/")
    console = Console()
    console.print(f"[green]Removed tag assignment {id}[/green]")


# These methods seem to be web frontend specifc

# @project.command()
# @click.argument('id', required=False)
# @click.option('--search', help='Search term')
# @click.option('--ordering', help='Field to order results by')
# @click.option('--filter', 'filters', multiple=True, help='Filter in format field=value')
# @click.option('--json-output', is_flag=True, help='Output in JSON format')
# @click.option('--all', 'fetch_all', is_flag=True, help='Fetch all results')
# def get_workbench_layout(id: Optional[str] = None, search: Optional[str] = None,
#                         ordering: Optional[str] = None, filters: tuple = (),
#                         json_output: bool = False, fetch_all: bool = False):
#     """Get project workbench layout details or list all layouts"""
#     client = get_client()
#     console = Console()

#     if id:
#         result = client.get(f'/api/project-workbench-layout/{id}/')
#         if json_output:
#             console.print_json(data=result)
#             return

#         console.print(Panel(f"[bold]Project Workbench Layout Details[/bold]"))
#         console.print(f"ID: {result['id']}")
#         console.print(f"Project: {result['project'].get('name', 'N/A')} ({result['project']['id']})")
#         console.print(f"Name: {result.get('name', 'N/A')}")
#         console.print(f"Created: {format_datetime(result.get('created_at'))}")
#         console.print(f"Updated: {format_datetime(result.get('updated_at'))}")
#         if result.get('layout'):
#             console.print("\n[bold]Layout Data[/bold]")
#             console.print_json(data=result['layout'])
#     else:
#         params = {}
#         if search is not None:
#             params['search'] = search
#         if ordering:
#             params['ordering'] = ordering

#         for f in filters:
#             try:
#                 field, value = f.split('=', 1)
#                 params[field] = value
#             except ValueError:
#                 console.print(f"[red]Invalid filter format: {f}. Use field=value[/red]")
#                 return

#         if fetch_all:
#             results = client.get_all('/api/project-workbench-layout/', params=params)
#         else:
#             results = client.get('/api/project-workbench-layout/', params=params)

#         if json_output:
#             console.print_json(data=results)
#             return

#         columns = [
#             ("ID", "id", None),
#             ("Project", "project", lambda x: x.get('name') if isinstance(x, dict) else None),
#             ("Name", "name", None),
#             ("Created", "created_at", format_datetime),
#             ("Updated", "updated_at", format_datetime),
#         ]

#         items, footer_string = get_pagination_info(results)

#         if not items:
#             console.print(f"[white]No layouts found with the given criteria:[/white]\n[green]{params}")
#             return

#         table = create_table(
#             columns=columns,
#             items=items,
#             title="Project Workbench Layouts",
#             caption=footer_string
#         )

#         console.print(table)


# @project.command()
# @click.option('--project', required=True, help='Project ID')
# @click.option('--name', required=True, help='Layout name')
# @click.option('--layout', required=True, help='Layout data (JSON string)')
# @click.option('--json-output', is_flag=True, help='Output in JSON format')
# def create_workbench_layout(project: str, name: str, layout: str, json_output: bool = False):
#     """Create a new project workbench layout"""
#     client = get_client()
#     console = Console()

#     try:
#         layout_data = json.loads(layout)
#     except json.JSONDecodeError:
#         console.print("[red]Invalid JSON format for layout data[/red]")
#         return

#     data = {
#         'project': project,
#         'name': name,
#         'layout': layout_data
#     }

#     result = client.post('/api/project-workbench-layout/', data=data)

#     if json_output:
#         console.print_json(data=result)
#     else:
#         console.print(f"[green]Created workbench layout '{name}' for project {project}[/green]")


# @project.command()
# @click.argument('id')
# def delete_workbench_layout(id: str):
#     """Delete a project workbench layout"""
#     client = get_client()
#     client.delete(f'/api/project-workbench-layout/{id}/')
#     console = Console()
#     console.print(f"[green]Deleted workbench layout {id}[/green]")


@project_group.command()
@click.argument("id", required=False)
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def get_molecule(
    id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """Get project molecule details or list all molecules"""
    client = get_client()
    console = Console()

    if id:
        result = client.get(f"/api/project-molecule/{id}/")
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel(f"[bold]Project Molecule Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Heading: {result.get('heading', 'N/A')}")
        console.print(f"Project: {result['project']}")
        console.print(f"Created: {format_datetime(result.get('created_at'))}")
    else:
        params = {}
        if search is not None:
            params["search"] = search
        if ordering:
            params["ordering"] = ordering

        for f in filters:
            try:
                field, value = f.split("=", 1)
                params[field] = value
            except ValueError:
                click.echo(
                    f"[red]Invalid filter format: {f}. Use field=value[/red]", err=True
                )
                return

        if fetch_all:
            results = client.get_all("/api/project-molecule/", params=params)
        else:
            results = client.get("/api/project-molecule/", params=params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Heading", "heading", None),
            (
                "Project",
                "project",
                lambda x: x.get("name") if isinstance(x, dict) else None,
            ),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No molecules found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns,
            items=items,
            title="Project Molecules",
            caption=footer_string,
        )

        console.print(table)


@project_group.command()
@click.option("--project", required=True, help="Project ID")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create_molecule(project: str, json_output: bool = False):
    """Add a molecule to a project"""
    client = get_client()
    console = Console()

    data = {
        "project": project,
    }

    result = client.post("/api/project-molecule/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"[green]Added project-molecule {project}[/green]")


@project_group.command()
@click.argument("id")
def delete_molecule(id: str):
    """Remove a molecule from a project"""
    client = get_client()
    client.delete(f"/api/project-molecule/{id}/")
    console = Console()
    console.print(f"[green]Removed molecule {id} from project[/green]")
