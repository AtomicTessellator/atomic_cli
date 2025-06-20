from typing import Any, Dict, List, Optional

from atomict.api import delete, get, patch, post
from atomict.exceptions import APIValidationError, PermissionDenied

from .exceptions import (
    InvalidOrganizationDataError,
    OrganizationInviteNotFoundError,
    OrganizationNotFoundError,
    OrganizationPermissionError,
    OrganizationUserNotFoundError,
)


def list_organizations(**params: Any) -> Any:
    """
    List all organizations accessible to the user

    Args:
        **params: Additional query parameters for filtering

    Returns:
        dict: List of organizations with pagination info
    """
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    url = (
        "api/organisation/" if not query_string else f"api/organisation/?{query_string}"
    )

    try:
        return get(url)
    except PermissionDenied as e:
        raise OrganizationPermissionError("Access denied to organizations") from e


def get_organization(org_id: str) -> Any:
    """
    Get organization details by ID

    Args:
        org_id: str - UUID of the organization

    Returns:
        dict: Organization data

    Raises:
        OrganizationNotFoundError: If organization doesn't exist
        OrganizationPermissionError: If user lacks access
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")

    try:
        return get(f"api/organisation/{org_id}/")
    except APIValidationError as e:
        if "not found" in str(e).lower():
            raise OrganizationNotFoundError(f"Organization {org_id} not found") from e
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to organization {org_id}"
        ) from e


def create_organization(
    name: str,
    description: Optional[str] = None,
    billing_primary_email: Optional[str] = None,
    chem_primary_email: Optional[str] = None,
    # Keep billing_emails for backward compatibility but map to billing_primary_email
    billing_emails: Optional[List[str]] = None,
) -> Any:
    """
    Create a new organization

    Args:
        name: str - Name of the organization
        description: str - Optional description
        billing_primary_email: str - Optional primary billing email
        chem_primary_email: str - Optional primary chemistry email
        billing_emails: List[str] - Deprecated, use billing_primary_email instead

    Returns:
        dict: Created organization data

    Raises:
        InvalidOrganizationDataError: If required data is missing or invalid
        OrganizationPermissionError: If user lacks permission to create organizations
    """
    if not name:
        raise InvalidOrganizationDataError("Organization name is required")

    payload: Dict[str, Any] = {"name": name}

    if description:
        payload["description"] = description
    if billing_primary_email:
        payload["billing_primary_email"] = billing_primary_email
    if chem_primary_email:
        payload["chem_primary_email"] = chem_primary_email
    # Handle backward compatibility - use first email as primary billing email
    if billing_emails and billing_emails[0] and not billing_primary_email:
        payload["billing_primary_email"] = billing_emails[0]

    try:
        return post("api/organisation/", payload)
    except APIValidationError as e:
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError("Access denied to create organization") from e


def update_organization(org_id: str, **kwargs: Any) -> Any:
    """
    Update organization details

    Args:
        org_id: str - UUID of the organization
        **kwargs: Fields to update (name, description, billing_primary_email, chem_primary_email)

    Returns:
        dict: Updated organization data

    Raises:
        OrganizationNotFoundError: If organization doesn't exist
        InvalidOrganizationDataError: If update data is invalid
        OrganizationPermissionError: If user lacks permission
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")

    if not kwargs:
        raise InvalidOrganizationDataError(
            "At least one field must be provided for update"
        )

    # Filter valid fields - support new field names and backward compatibility
    valid_fields = {"name", "description", "billing_primary_email", "chem_primary_email", "billing_emails"}
    payload = {}
    
    for k, v in kwargs.items():
        if k in valid_fields:
            if k == "billing_emails" and v:
                # Handle backward compatibility - convert to billing_primary_email
                if isinstance(v, list) and v[0]:
                    payload["billing_primary_email"] = v[0]
            else:
                payload[k] = v

    if not payload:
        raise InvalidOrganizationDataError(
            f"No valid fields provided. Valid fields: {valid_fields}"
        )

    try:
        return patch(f"api/organisation/{org_id}/", payload)
    except APIValidationError as e:
        if "not found" in str(e).lower():
            raise OrganizationNotFoundError(f"Organization {org_id} not found") from e
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to update organization {org_id}"
        ) from e


def delete_organization(org_id: str) -> Any:
    """
    Delete an organization

    Args:
        org_id: str - UUID of the organization

    Returns:
        dict: Deletion response

    Raises:
        OrganizationNotFoundError: If organization doesn't exist
        OrganizationPermissionError: If user lacks permission
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")

    try:
        return delete(f"api/organisation/{org_id}/")
    except APIValidationError as e:
        if "not found" in str(e).lower():
            raise OrganizationNotFoundError(f"Organization {org_id} not found") from e
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to delete organization {org_id}"
        ) from e


def list_organization_users(org_id: str, **params: Any) -> Any:
    """
    List users in an organization

    Args:
        org_id: str - UUID of the organization
        **params: Additional query parameters

    Returns:
        dict: List of organization users

    Raises:
        OrganizationNotFoundError: If organization doesn't exist
        OrganizationPermissionError: If user lacks access
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")

    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    url = f"api/organisation/{org_id}/users"
    if query_string:
        url = f"{url}?{query_string}"

    try:
        return get(url)
    except APIValidationError as e:
        if "not found" in str(e).lower():
            raise OrganizationNotFoundError(f"Organization {org_id} not found") from e
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to organization {org_id} users"
        ) from e


def add_user_to_organization(
    org_id: str, user_email: str, is_admin: bool = False
) -> Any:
    """
    Add a user to an organization

    Args:
        org_id: str - UUID of the organization
        user_email: str - Email address of the user to add
        is_admin: bool - Whether the user should have admin privileges

    Returns:
        dict: Organization user data

    Raises:
        InvalidOrganizationDataError: If required data is missing
        OrganizationPermissionError: If user lacks permission
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")
    if not user_email:
        raise InvalidOrganizationDataError("User email is required")

    payload = {
        "organisation": org_id,
        "user_email": user_email,
        "is_admin": is_admin,
    }

    try:
        return post("api/organisation-user/", payload)
    except APIValidationError as e:
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to add user to organization {org_id}"
        ) from e


def remove_user_from_organization(org_user_id: str) -> Any:
    """
    Remove a user from an organization

    Args:
        org_user_id: str - ID of the organization user record

    Returns:
        dict: Deletion response

    Raises:
        OrganizationUserNotFoundError: If organization user doesn't exist
        OrganizationPermissionError: If user lacks permission
    """
    if not org_user_id:
        raise InvalidOrganizationDataError("Organization user ID is required")

    try:
        return delete(f"api/organisation-user/{org_user_id}/")
    except APIValidationError as e:
        if "not found" in str(e).lower():
            raise OrganizationUserNotFoundError(
                f"Organization user {org_user_id} not found"
            ) from e
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to remove organization user {org_user_id}"
        ) from e


def list_organization_invites(org_id: str, **params: Any) -> Any:
    """
    List pending invites for an organization

    Args:
        org_id: str - UUID of the organization
        **params: Additional query parameters

    Returns:
        dict: List of organization invites

    Raises:
        OrganizationNotFoundError: If organization doesn't exist
        OrganizationPermissionError: If user lacks access
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")

    # Use the organisation query parameter to filter invites by organization
    all_params = {"organisation": org_id, **params}
    query_string = "&".join(f"{k}={v}" for k, v in all_params.items())
    url = f"api/organisation-invite/?{query_string}"

    try:
        return get(url)
    except APIValidationError as e:
        if "not found" in str(e).lower():
            raise OrganizationNotFoundError(f"Organization {org_id} not found") from e
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to organization {org_id} invites"
        ) from e


def send_organization_invite(org_id: str, email: str, is_admin: bool = False) -> Any:
    """
    Send an invite to join an organization

    Args:
        org_id: str - UUID of the organization
        email: str - Email address to send invite to
        is_admin: bool - Currently ignored, admin status set when invite is accepted

    Returns:
        dict: Organization invite data

    Raises:
        InvalidOrganizationDataError: If required data is missing
        OrganizationPermissionError: If user lacks permission
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")
    if not email:
        raise InvalidOrganizationDataError("Email is required")

    payload = {
        "organisation": org_id,
        "email": email,
        # Note: is_admin not included - backend doesn't support it for invites yet
    }

    try:
        return post("api/organisation-invite/", payload)
    except APIValidationError as e:
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to invite user to organization {org_id}"
        ) from e


def delete_organization_invite(invite_id: str) -> Any:
    """
    Delete/cancel an organization invite

    Args:
        invite_id: str - ID of the organization invite

    Returns:
        dict: Deletion response

    Raises:
        OrganizationInviteNotFoundError: If invite doesn't exist
        OrganizationPermissionError: If user lacks permission
    """
    if not invite_id:
        raise InvalidOrganizationDataError("Invite ID is required")

    try:
        return delete(f"api/organisation-invite/{invite_id}/")
    except APIValidationError as e:
        if "not found" in str(e).lower():
            raise OrganizationInviteNotFoundError(
                f"Organization invite {invite_id} not found"
            ) from e
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to delete organization invite {invite_id}"
        ) from e


def get_active_organization() -> Any:
    """
    Get the user's currently active organization

    Returns:
        dict: Active organization data

    Raises:
        OrganizationPermissionError: If user lacks access
    """
    try:
        return get("api/preferences/get-active-org/")
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            "Access denied to get active organization"
        ) from e


def set_active_organization(org_id: str) -> Any:
    """
    Set the user's active organization

    Args:
        org_id: str - UUID of the organization to set as active

    Returns:
        dict: Response confirming the active organization change

    Raises:
        InvalidOrganizationDataError: If organization ID is missing
        OrganizationPermissionError: If user lacks access
    """
    if not org_id:
        raise InvalidOrganizationDataError("Organization ID is required")

    payload = {"organisation": org_id}

    try:
        return post("api/preferences/set-active-org/", payload)
    except APIValidationError as e:
        raise InvalidOrganizationDataError(str(e)) from e
    except PermissionDenied as e:
        raise OrganizationPermissionError(
            f"Access denied to set active organization {org_id}"
        ) from e
