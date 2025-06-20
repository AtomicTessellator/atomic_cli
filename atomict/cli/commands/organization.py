# cli/commands/organization.py
import json
import logging
from typing import Any, Dict, Optional

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from atomict.cli.commands.common import create_table
from atomict.cli.commands.helpers import format_datetime, get_status_string
from atomict.cli.core.client import get_client
from atomict.cli.core.utils import get_pagination_info
from atomict.organization.client import (
    list_organizations, 
    get_organization, 
    create_organization, 
    update_organization, 
    delete_organization,
    list_organization_users,
    add_user_to_organization,
    remove_user_from_organization,
    list_organization_invites,
    send_organization_invite,
    delete_organization_invite,
    get_active_organization,
    set_active_organization,
)

console = Console()
logger = logging.getLogger(__name__)


@click.group(name="org")
def org_group():
    """Manage organizations"""
    pass


@org_group.command()
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
    """Get organization details or list all organizations"""
    client = get_client()

    if id:
        result = get_organization(id)
        if json_output:
            console.print_json(data=result)
            return

        console.print(Panel(f"[bold]Organization Details[/bold]"))
        console.print(f"ID: {result['id']}")
        console.print(f"Name: {result.get('name', 'N/A')}")
        console.print(f"Description: {result.get('description', 'N/A')}")
        # Organization model doesn't have created_at/updated_at timestamps
        
        # Show billing emails - handle both old and new formats
        if result.get("billing_primary_email"):
            console.print(f"Billing Primary Email: {result['billing_primary_email']}")
        if result.get("chem_primary_email"):
            console.print(f"Chemistry Primary Email: {result['chem_primary_email']}")
        
        # Backward compatibility for billing_emails
        if result.get("billing_emails"):
            console.print(f"Billing Emails: {', '.join(result['billing_emails'])}")

        # Show users if present
        if "users" in result:
            console.print("\n[bold]Users[/bold]")
            user_table = Table(show_header=True)
            user_table.add_column("ID")
            user_table.add_column("Username")
            user_table.add_column("Email")
            user_table.add_column("Role")

            for user in result["users"]:
                user_table.add_row(
                    str(user["id"]),
                    user.get("username", "N/A"),
                    user.get("email", "N/A"),
                    user.get("role", "N/A"),
                )
            console.print(user_table)
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
                raise click.ClickException(f"Invalid filter format: {f}. Use field=value")

        results = list_organizations(**params)

        if json_output:
            console.print_json(data=results)
            return

        columns = [
            ("ID", "id", None),
            ("Name", "name", None),
            ("Users", "users", lambda x: str(len(x)) if x else "0"),
            ("Created", "created_at", format_datetime),
        ]

        items, footer_string = get_pagination_info(results)

        if not items:
            console.print(
                f"[white]No organizations found with the given criteria:[/white]\n[green]{params}"
            )
            return

        table = create_table(
            columns=columns, items=items, title="Organizations", caption=footer_string
        )

        console.print(table)


@org_group.command()
@click.option("--name", required=True, help="Organization name")
@click.option("--description", help="Organization description")
@click.option("--billing-emails", help="Comma-separated billing email addresses")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def create(
    name: str,
    description: Optional[str] = None,
    billing_emails: Optional[str] = None,
    json_output: bool = False,
):
    """Create a new organization"""
    
    # Prepare parameters for the organization client function
    kwargs: Dict[str, Any] = {"name": name}
    if description:
        kwargs["description"] = description
    if billing_emails:
        # Use the first email as the primary billing email
        email_list = [email.strip() for email in billing_emails.split(",")]
        if email_list:
            kwargs["billing_emails"] = email_list  # This will be converted to billing_primary_email in the function

    result = create_organization(**kwargs)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Created organization '{name}' with ID: {result['id']}[/green]"
        )


@org_group.command()
@click.argument("id")
@click.option("--name", help="Organization name")
@click.option("--description", help="Organization description")
@click.option("--billing-emails", help="Comma-separated billing email addresses")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def update(
    id: str,
    name: Optional[str] = None,
    description: Optional[str] = None,
    billing_emails: Optional[str] = None,
    json_output: bool = False,
):
    """Update an organization"""
    
    kwargs: Dict[str, Any] = {}
    if name:
        kwargs["name"] = name
    if description:
        kwargs["description"] = description
    if billing_emails:
        # Use the first email as the primary billing email
        email_list = [email.strip() for email in billing_emails.split(",")]
        if email_list:
            kwargs["billing_primary_email"] = email_list[0]

    if not kwargs:
        console.print(
            "[yellow]No fields to update. Use --name, --description, or --billing-emails[/yellow]"
        )
        return

    result = update_organization(org_id=id, **kwargs)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(f"[green]Updated organization {id}[/green]")


@org_group.command()
@click.argument("id")
@click.option("--yes", is_flag=True, help="Skip confirmation prompt")
def delete(id: str, yes: bool = False):
    """Delete an organization"""

    if not yes:
        if not click.confirm(f"Are you sure you want to delete organization {id}?"):
            console.print("[yellow]Operation cancelled[/yellow]")
            return

    delete_organization(id)
    console.print(f"[green]Deleted organization {id}[/green]")


@org_group.command()
@click.argument("id")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def switch(id: str, json_output: bool = False):
    """Switch to a different organization"""
    client = get_client()

    # First verify the organization exists
    try:
        org = get_organization(id)
    except Exception as e:
        console.print(
            f"[red]Error: Organization {id} not found or not accessible[/red]"
        )
        return

    # Set as active organization
    result = set_active_organization(id)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Switched to organization '{org.get('name', id)}'[/green]"
        )


@click.group(name="users")
def users_group():
    """Manage organization users"""
    pass


@users_group.command("list")
@click.argument("org_id")
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def list_users(
    org_id: str,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """List users in an organization"""
    client = get_client()

    params = {"organisation": org_id}
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
        results = client.get_all("/api/organisation-user/", params=params)
    else:
        results = client.get("/api/organisation-user/", params=params)

    if json_output:
        console.print_json(data=results)
        return

    columns = [
        ("User ID", "user", lambda x: x.get("id") if isinstance(x, dict) else x),
        (
            "Username",
            "user",
            lambda x: x.get("username") if isinstance(x, dict) else "N/A",
        ),
        ("Email", "user", lambda x: x.get("email") if isinstance(x, dict) else "N/A"),
        ("Role", "role", None),
        ("Joined", "created_at", format_datetime),
    ]

    items, footer_string = get_pagination_info(results)

    if not items:
        console.print(f"[white]No users found in organization {org_id}[/white]")
        return

    table = create_table(
        columns=columns,
        items=items,
        title=f"Organization {org_id} Users",
        caption=footer_string,
    )

    console.print(table)


@users_group.command()
@click.argument("org_id")
@click.argument("email")
@click.option("--admin", is_flag=True, help="Add user as admin")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def add(org_id: str, email: str, admin: bool = False, json_output: bool = False):
    """Add a user to an organization"""
    client = get_client()

    data = {
        "organisation": org_id,
        "user_email": email,
        "role": "admin" if admin else "member",
    }

    result = client.post("/api/organisation-user/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        role = "admin" if admin else "member"
        console.print(
            f"[green]Added user {email} to organization {org_id} as {role}[/green]"
        )


@users_group.command()
@click.argument("org_id")
@click.argument("user_id")
@click.option("--yes", is_flag=True, help="Skip confirmation prompt")
def remove(org_id: str, user_id: str, yes: bool = False):
    """Remove a user from an organization"""
    client = get_client()

    if not yes:
        if not click.confirm(
            f"Are you sure you want to remove user {user_id} from organization {org_id}?"
        ):
            console.print("[yellow]Operation cancelled[/yellow]")
            return

    # Find the organization-user relationship
    params = {"organisation": org_id, "user": user_id}
    relationships = client.get("/api/organisation-user/", params=params)

    items, _ = get_pagination_info(relationships)
    if not items:
        console.print(f"[red]User {user_id} not found in organization {org_id}[/red]")
        return

    # Delete the relationship
    rel_id = items[0]["id"]
    client.delete(f"/api/organisation-user/{rel_id}/")
    console.print(f"[green]Removed user {user_id} from organization {org_id}[/green]")


@users_group.command()
@click.argument("org_id")
@click.argument("user_id")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def promote(org_id: str, user_id: str, json_output: bool = False):
    """Promote a user to admin in an organization"""
    client = get_client()

    # Find the organization-user relationship
    params = {"organisation": org_id, "user": user_id}
    relationships = client.get("/api/organisation-user/", params=params)

    items, _ = get_pagination_info(relationships)
    if not items:
        console.print(f"[red]User {user_id} not found in organization {org_id}[/red]")
        return

    # Update the relationship
    rel_id = items[0]["id"]
    data = {"role": "admin"}
    result = client.patch(f"/api/organisation-user/{rel_id}/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Promoted user {user_id} to admin in organization {org_id}[/green]"
        )


@users_group.command()
@click.argument("org_id")
@click.argument("user_id")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def demote(org_id: str, user_id: str, json_output: bool = False):
    """Demote a user from admin to member in an organization"""
    client = get_client()

    # Find the organization-user relationship
    params = {"organisation": org_id, "user": user_id}
    relationships = client.get("/api/organisation-user/", params=params)

    items, _ = get_pagination_info(relationships)
    if not items:
        console.print(f"[red]User {user_id} not found in organization {org_id}[/red]")
        return

    # Update the relationship
    rel_id = items[0]["id"]
    data = {"role": "member"}
    result = client.patch(f"/api/organisation-user/{rel_id}/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        console.print(
            f"[green]Demoted user {user_id} to member in organization {org_id}[/green]"
        )


# Add users subcommands to org group
org_group.add_command(users_group)


@click.group(name="invites")
def invites_group():
    """Manage organization invitations"""
    pass


@invites_group.command("list")
@click.argument("org_id")
@click.option("--search", help="Search term")
@click.option("--ordering", help="Field to order results by")
@click.option("--filter", "filters", multiple=True, help="Filter in format field=value")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
@click.option("--all", "fetch_all", is_flag=True, help="Fetch all results")
def list_invites(
    org_id: str,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    filters: tuple = (),
    json_output: bool = False,
    fetch_all: bool = False,
):
    """List pending invitations for an organization"""
    client = get_client()

    params = {"organisation": org_id}
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
        results = client.get_all("/api/organisation-invite/", params=params)
    else:
        results = client.get("/api/organisation-invite/", params=params)

    if json_output:
        console.print_json(data=results)
        return

    columns = [
        ("ID", "id", None),
        ("Email", "email", None),
        ("Role", "role", None),
        ("Status", "status", None),
        ("Invited", "created_at", format_datetime),
        ("Expires", "expires_at", format_datetime),
    ]

    items, footer_string = get_pagination_info(results)

    if not items:
        console.print(
            f"[white]No pending invitations found for organization {org_id}[/white]"
        )
        return

    table = create_table(
        columns=columns,
        items=items,
        title=f"Organization {org_id} Invitations",
        caption=footer_string,
    )

    console.print(table)


@invites_group.command()
@click.argument("org_id")
@click.argument("email")
@click.option("--admin", is_flag=True, help="Invite user as admin")
@click.option("--json-output", is_flag=True, help="Output in JSON format")
def send(org_id: str, email: str, admin: bool = False, json_output: bool = False):
    """Send an invitation to join an organization"""
    client = get_client()

    data = {
        "organisation": org_id,
        "email": email,
        "role": "admin" if admin else "member",
    }

    result = client.post("/api/organisation-invite/", data=data)

    if json_output:
        console.print_json(data=result)
    else:
        role = "admin" if admin else "member"
        console.print(
            f"[green]Sent invitation to {email} to join organization {org_id} as {role}[/green]"
        )


@invites_group.command()
@click.argument("invite_id")
@click.option("--yes", is_flag=True, help="Skip confirmation prompt")
def cancel(invite_id: str, yes: bool = False):
    """Cancel a pending invitation"""
    client = get_client()

    if not yes:
        if not click.confirm(
            f"Are you sure you want to cancel invitation {invite_id}?"
        ):
            console.print("[yellow]Operation cancelled[/yellow]")
            return

    client.delete(f"/api/organisation-invite/{invite_id}/")
    console.print(f"[green]Cancelled invitation {invite_id}[/green]")


# Add invites subcommands to org group
org_group.add_command(invites_group)

org = org_group
