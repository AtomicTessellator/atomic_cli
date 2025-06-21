from unittest.mock import Mock, patch

import pytest

from atomict.exceptions import APIValidationError, PermissionDenied
from atomict.organization.client import (
    add_user_to_organization,
    create_organization,
    delete_organization,
    delete_organization_invite,
    get_active_organization,
    get_organization,
    list_organization_invites,
    list_organization_users,
    list_organizations,
    remove_user_from_organization,
    send_organization_invite,
    set_active_organization,
    update_organization,
)
from atomict.organization.exceptions import (
    InvalidOrganizationDataError,
    OrganizationInviteNotFoundError,
    OrganizationNotFoundError,
    OrganizationPermissionError,
    OrganizationUserNotFoundError,
)


@pytest.fixture
def mock_get():
    with patch("atomict.organization.client.get") as mock:
        yield mock


@pytest.fixture
def mock_post():
    with patch("atomict.organization.client.post") as mock:
        yield mock


@pytest.fixture
def mock_patch():
    with patch("atomict.organization.client.patch") as mock:
        yield mock


@pytest.fixture
def mock_delete():
    with patch("atomict.organization.client.delete") as mock:
        yield mock


class TestListOrganizations:
    def test_list_organizations_without_params(self, mock_get):
        """Test listing organizations without additional parameters"""
        mock_get.return_value = {"results": [{"id": "123", "name": "Test Org"}]}

        result = list_organizations()

        assert result == {"results": [{"id": "123", "name": "Test Org"}]}
        mock_get.assert_called_once_with("api/organisation/")

    def test_list_organizations_with_params(self, mock_get):
        """Test listing organizations with query parameters"""
        mock_get.return_value = {"results": []}

        result = list_organizations(search="test", ordering="name")

        assert result == {"results": []}
        mock_get.assert_called_once_with("api/organisation/?search=test&ordering=name")

    def test_list_organizations_permission_denied(self, mock_get):
        """Test permission denied error handling"""
        mock_get.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError, match="Access denied to organizations"
        ):
            list_organizations()


class TestGetOrganization:
    def test_get_organization_success(self, mock_get):
        """Test successful organization retrieval"""
        expected_org = {
            "id": "123",
            "name": "Test Org",
            "description": "Test description",
        }
        mock_get.return_value = expected_org

        result = get_organization("123")

        assert result == expected_org
        mock_get.assert_called_once_with("api/organisation/123/")

    def test_get_organization_empty_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            get_organization("")

    def test_get_organization_not_found(self, mock_get):
        """Test organization not found error"""
        mock_get.side_effect = APIValidationError("Organization not found")

        with pytest.raises(
            OrganizationNotFoundError, match="Organization 123 not found"
        ):
            get_organization("123")

    def test_get_organization_permission_denied(self, mock_get):
        """Test permission denied error"""
        mock_get.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError, match="Access denied to organization 123"
        ):
            get_organization("123")

    def test_get_organization_api_validation_error(self, mock_get):
        """Test general API validation error"""
        mock_get.side_effect = APIValidationError("Invalid request")

        with pytest.raises(InvalidOrganizationDataError, match="Invalid request"):
            get_organization("123")


class TestCreateOrganization:
    def test_create_organization_minimal(self, mock_post):
        """Test creating organization with minimal required data"""
        expected_org = {"id": "123", "name": "Test Org"}
        mock_post.return_value = expected_org

        result = create_organization("Test Org")

        assert result == expected_org
        mock_post.assert_called_once_with("api/organisation/", {"name": "Test Org"})

    def test_create_organization_full_data(self, mock_post):
        """Test creating organization with all optional fields"""
        expected_org = {"id": "123", "name": "Test Org"}
        mock_post.return_value = expected_org

        result = create_organization(
            "Test Org",
            description="Test description",
            billing_primary_email="billing@test.com",
            chem_primary_email="chem@test.com",
        )

        assert result == expected_org
        mock_post.assert_called_once_with(
            "api/organisation/",
            {
                "name": "Test Org",
                "description": "Test description",
                "billing_primary_email": "billing@test.com",
                "chem_primary_email": "chem@test.com",
            },
        )

    def test_create_organization_billing_emails_backward_compatibility(self, mock_post):
        """Test backward compatibility with billing_emails parameter"""
        expected_org = {"id": "123", "name": "Test Org"}
        mock_post.return_value = expected_org
        billing_emails = ["billing@test.com", "admin@test.com"]

        result = create_organization(
            "Test Org", description="Test description", billing_emails=billing_emails
        )

        assert result == expected_org
        mock_post.assert_called_once_with(
            "api/organisation/",
            {
                "name": "Test Org",
                "description": "Test description",
                "billing_primary_email": "billing@test.com",  # First email used as primary
            },
        )

    def test_create_organization_empty_name(self):
        """Test error handling for empty organization name"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization name is required"
        ):
            create_organization("")

    def test_create_organization_api_validation_error(self, mock_post):
        """Test API validation error handling"""
        mock_post.side_effect = APIValidationError("Name already exists")

        with pytest.raises(InvalidOrganizationDataError, match="Name already exists"):
            create_organization("Test Org")

    def test_create_organization_permission_denied(self, mock_post):
        """Test permission denied error"""
        mock_post.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError, match="Access denied to create organization"
        ):
            create_organization("Test Org")


class TestUpdateOrganization:
    def test_update_organization_success(self, mock_patch):
        """Test successful organization update"""
        expected_org = {"id": "123", "name": "Updated Org"}
        mock_patch.return_value = expected_org

        result = update_organization("123", name="Updated Org")

        assert result == expected_org
        mock_patch.assert_called_once_with(
            "api/organisation/123/", {"name": "Updated Org"}
        )

    def test_update_organization_multiple_fields(self, mock_patch):
        """Test updating multiple fields"""
        expected_org = {"id": "123", "name": "Updated Org"}
        mock_patch.return_value = expected_org

        result = update_organization(
            "123",
            name="Updated Org",
            description="New description",
            billing_primary_email="new@test.com",
        )

        assert result == expected_org
        mock_patch.assert_called_once_with(
            "api/organisation/123/",
            {
                "name": "Updated Org",
                "description": "New description",
                "billing_primary_email": "new@test.com",
            },
        )

    def test_update_organization_billing_emails_backward_compatibility(
        self, mock_patch
    ):
        """Test backward compatibility with billing_emails parameter"""
        expected_org = {"id": "123", "name": "Updated Org"}
        mock_patch.return_value = expected_org

        result = update_organization(
            "123", name="Updated Org", billing_emails=["new@test.com", "admin@test.com"]
        )

        assert result == expected_org
        mock_patch.assert_called_once_with(
            "api/organisation/123/",
            {
                "name": "Updated Org",
                "billing_primary_email": "new@test.com",  # First email used as primary
            },
        )

    def test_update_organization_empty_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            update_organization("", name="Test")

    def test_update_organization_no_fields(self):
        """Test error handling when no fields provided"""
        with pytest.raises(
            InvalidOrganizationDataError, match="At least one field must be provided"
        ):
            update_organization("123")

    def test_update_organization_invalid_fields(self):
        """Test error handling for invalid fields"""
        with pytest.raises(
            InvalidOrganizationDataError, match="No valid fields provided"
        ):
            update_organization("123", invalid_field="value")

    def test_update_organization_not_found(self, mock_patch):
        """Test organization not found error"""
        mock_patch.side_effect = APIValidationError("Organization not found")

        with pytest.raises(
            OrganizationNotFoundError, match="Organization 123 not found"
        ):
            update_organization("123", name="Updated")

    def test_update_organization_permission_denied(self, mock_patch):
        """Test permission denied error"""
        mock_patch.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to update organization 123",
        ):
            update_organization("123", name="Updated")


class TestDeleteOrganization:
    def test_delete_organization_success(self, mock_delete):
        """Test successful organization deletion"""
        mock_delete.return_value = {"message": "Deleted successfully"}

        result = delete_organization("123")

        assert result == {"message": "Deleted successfully"}
        mock_delete.assert_called_once_with("api/organisation/123/")

    def test_delete_organization_empty_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            delete_organization("")

    def test_delete_organization_not_found(self, mock_delete):
        """Test organization not found error"""
        mock_delete.side_effect = APIValidationError("Organization not found")

        with pytest.raises(
            OrganizationNotFoundError, match="Organization 123 not found"
        ):
            delete_organization("123")

    def test_delete_organization_permission_denied(self, mock_delete):
        """Test permission denied error"""
        mock_delete.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to delete organization 123",
        ):
            delete_organization("123")


class TestListOrganizationUsers:
    def test_list_organization_users_success(self, mock_get):
        """Test successful listing of organization users"""
        expected_users = {"results": [{"id": "1", "email": "user@test.com"}]}
        mock_get.return_value = expected_users

        result = list_organization_users("123")

        assert result == expected_users
        mock_get.assert_called_once_with("api/organisation/123/users")

    def test_list_organization_users_with_params(self, mock_get):
        """Test listing users with query parameters"""
        expected_users = {"results": []}
        mock_get.return_value = expected_users

        result = list_organization_users("123", search="test", limit=10)

        assert result == expected_users
        mock_get.assert_called_once_with(
            "api/organisation/123/users?search=test&limit=10"
        )

    def test_list_organization_users_empty_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            list_organization_users("")

    def test_list_organization_users_not_found(self, mock_get):
        """Test organization not found error"""
        mock_get.side_effect = APIValidationError("Organization not found")

        with pytest.raises(
            OrganizationNotFoundError, match="Organization 123 not found"
        ):
            list_organization_users("123")

    def test_list_organization_users_permission_denied(self, mock_get):
        """Test permission denied error"""
        mock_get.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError, match="Access denied to organization 123 users"
        ):
            list_organization_users("123")


class TestAddUserToOrganization:
    def test_add_user_to_organization_success(self, mock_post):
        """Test successful user addition to organization"""
        expected_user = {"id": "1", "user_email": "user@test.com", "is_admin": False}
        mock_post.return_value = expected_user

        result = add_user_to_organization("123", "user@test.com")

        assert result == expected_user
        mock_post.assert_called_once_with(
            "api/organisation-user/",
            {"organisation": "123", "user_email": "user@test.com", "is_admin": False},
        )

    def test_add_user_to_organization_as_admin(self, mock_post):
        """Test adding user as admin"""
        expected_user = {"id": "1", "user_email": "admin@test.com", "is_admin": True}
        mock_post.return_value = expected_user

        result = add_user_to_organization("123", "admin@test.com", is_admin=True)

        assert result == expected_user
        mock_post.assert_called_once_with(
            "api/organisation-user/",
            {"organisation": "123", "user_email": "admin@test.com", "is_admin": True},
        )

    def test_add_user_to_organization_empty_org_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            add_user_to_organization("", "user@test.com")

    def test_add_user_to_organization_empty_email(self):
        """Test error handling for empty user email"""
        with pytest.raises(
            InvalidOrganizationDataError, match="User email is required"
        ):
            add_user_to_organization("123", "")

    def test_add_user_to_organization_api_validation_error(self, mock_post):
        """Test API validation error handling"""
        mock_post.side_effect = APIValidationError("User already exists")

        with pytest.raises(InvalidOrganizationDataError, match="User already exists"):
            add_user_to_organization("123", "user@test.com")

    def test_add_user_to_organization_permission_denied(self, mock_post):
        """Test permission denied error"""
        mock_post.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to add user to organization 123",
        ):
            add_user_to_organization("123", "user@test.com")


class TestRemoveUserFromOrganization:
    def test_remove_user_from_organization_success(self, mock_delete):
        """Test successful user removal from organization"""
        mock_delete.return_value = {"message": "User removed successfully"}

        result = remove_user_from_organization("user123")

        assert result == {"message": "User removed successfully"}
        mock_delete.assert_called_once_with("api/organisation-user/user123/")

    def test_remove_user_from_organization_empty_id(self):
        """Test error handling for empty organization user ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization user ID is required"
        ):
            remove_user_from_organization("")

    def test_remove_user_from_organization_not_found(self, mock_delete):
        """Test organization user not found error"""
        mock_delete.side_effect = APIValidationError("Organization user not found")

        with pytest.raises(
            OrganizationUserNotFoundError, match="Organization user user123 not found"
        ):
            remove_user_from_organization("user123")

    def test_remove_user_from_organization_permission_denied(self, mock_delete):
        """Test permission denied error"""
        mock_delete.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to remove organization user user123",
        ):
            remove_user_from_organization("user123")


class TestListOrganizationInvites:
    def test_list_organization_invites_success(self, mock_get):
        """Test successful listing of organization invites"""
        expected_invites = {"results": [{"id": "1", "email": "invite@test.com"}]}
        mock_get.return_value = expected_invites

        result = list_organization_invites("123")

        assert result == expected_invites
        mock_get.assert_called_once_with("api/organisation-invite/?organisation=123")

    def test_list_organization_invites_with_params(self, mock_get):
        """Test listing invites with query parameters"""
        expected_invites = {"results": []}
        mock_get.return_value = expected_invites

        result = list_organization_invites("123", search="test", limit=5)

        assert result == expected_invites
        mock_get.assert_called_once_with(
            "api/organisation-invite/?organisation=123&search=test&limit=5"
        )

    def test_list_organization_invites_empty_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            list_organization_invites("")

    def test_list_organization_invites_not_found(self, mock_get):
        """Test organization not found error"""
        mock_get.side_effect = APIValidationError("Organization not found")

        with pytest.raises(
            OrganizationNotFoundError, match="Organization 123 not found"
        ):
            list_organization_invites("123")

    def test_list_organization_invites_permission_denied(self, mock_get):
        """Test permission denied error"""
        mock_get.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to organization 123 invites",
        ):
            list_organization_invites("123")


class TestSendOrganizationInvite:
    def test_send_organization_invite_success(self, mock_post):
        """Test successful organization invite sending"""
        expected_invite = {"id": "1", "email": "invite@test.com"}
        mock_post.return_value = expected_invite

        result = send_organization_invite("123", "invite@test.com")

        assert result == expected_invite
        mock_post.assert_called_once_with(
            "api/organisation-invite/",
            {"organisation": "123", "email": "invite@test.com"},
        )

    def test_send_organization_invite_as_admin(self, mock_post):
        """Test sending admin invite (is_admin parameter is ignored)"""
        expected_invite = {"id": "1", "email": "admin@test.com"}
        mock_post.return_value = expected_invite

        result = send_organization_invite("123", "admin@test.com", is_admin=True)

        assert result == expected_invite
        mock_post.assert_called_once_with(
            "api/organisation-invite/",
            {"organisation": "123", "email": "admin@test.com"},
        )

    def test_send_organization_invite_empty_org_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            send_organization_invite("", "invite@test.com")

    def test_send_organization_invite_empty_email(self):
        """Test error handling for empty email"""
        with pytest.raises(InvalidOrganizationDataError, match="Email is required"):
            send_organization_invite("123", "")

    def test_send_organization_invite_api_validation_error(self, mock_post):
        """Test API validation error handling"""
        mock_post.side_effect = APIValidationError("Invalid email address")

        with pytest.raises(InvalidOrganizationDataError, match="Invalid email address"):
            send_organization_invite("123", "invalid-email")

    def test_send_organization_invite_permission_denied(self, mock_post):
        """Test permission denied error"""
        mock_post.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to invite user to organization 123",
        ):
            send_organization_invite("123", "invite@test.com")


class TestDeleteOrganizationInvite:
    def test_delete_organization_invite_success(self, mock_delete):
        """Test successful organization invite deletion"""
        mock_delete.return_value = {"message": "Invite deleted successfully"}

        result = delete_organization_invite("invite123")

        assert result == {"message": "Invite deleted successfully"}
        mock_delete.assert_called_once_with("api/organisation-invite/invite123/")

    def test_delete_organization_invite_empty_id(self):
        """Test error handling for empty invite ID"""
        with pytest.raises(InvalidOrganizationDataError, match="Invite ID is required"):
            delete_organization_invite("")

    def test_delete_organization_invite_not_found(self, mock_delete):
        """Test organization invite not found error"""
        mock_delete.side_effect = APIValidationError("Organization invite not found")

        with pytest.raises(
            OrganizationInviteNotFoundError,
            match="Organization invite invite123 not found",
        ):
            delete_organization_invite("invite123")

    def test_delete_organization_invite_permission_denied(self, mock_delete):
        """Test permission denied error"""
        mock_delete.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to delete organization invite invite123",
        ):
            delete_organization_invite("invite123")


class TestGetActiveOrganization:
    def test_get_active_organization_success(self, mock_get):
        """Test successful retrieval of active organization"""
        expected_org = {"id": "123", "name": "Active Org"}
        mock_get.return_value = expected_org

        result = get_active_organization()

        assert result == expected_org
        mock_get.assert_called_once_with("api/preferences/get-active-org/")

    def test_get_active_organization_permission_denied(self, mock_get):
        """Test permission denied error"""
        mock_get.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to get active organization",
        ):
            get_active_organization()


class TestSetActiveOrganization:
    def test_set_active_organization_success(self, mock_post):
        """Test successful setting of active organization"""
        expected_response = {"message": "Active organization set successfully"}
        mock_post.return_value = expected_response

        result = set_active_organization("123")

        assert result == expected_response
        mock_post.assert_called_once_with(
            "api/preferences/set-active-org/", {"organisation": "123"}
        )

    def test_set_active_organization_empty_id(self):
        """Test error handling for empty organization ID"""
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            set_active_organization("")

    def test_set_active_organization_api_validation_error(self, mock_post):
        """Test API validation error handling"""
        mock_post.side_effect = APIValidationError("Invalid organization ID")

        with pytest.raises(
            InvalidOrganizationDataError, match="Invalid organization ID"
        ):
            set_active_organization("invalid")

    def test_set_active_organization_permission_denied(self, mock_post):
        """Test permission denied error"""
        mock_post.side_effect = PermissionDenied("Access denied")

        with pytest.raises(
            OrganizationPermissionError,
            match="Access denied to set active organization 123",
        ):
            set_active_organization("123")


class TestOrganizationClientIntegration:
    """Integration-style tests for the organization client module"""

    @patch("atomict.organization.client.get")
    @patch("atomict.organization.client.post")
    def test_organization_workflow_integration(self, mock_post, mock_get):
        """Test a typical workflow using multiple organization functions"""
        # Mock responses for workflow
        mock_post.side_effect = [
            {"id": "org123", "name": "Test Organization"},  # create_organization
            {
                "id": "user456",
                "user_email": "user@test.com",
            },  # add_user_to_organization
            {"id": "invite789", "email": "invite@test.com"},  # send_organization_invite
        ]
        mock_get.side_effect = [
            {
                "id": "org123",
                "name": "Test Organization",
                "users": [],
            },  # get_organization
            {
                "results": [{"id": "user456", "email": "user@test.com"}]
            },  # list_organization_users
            {
                "results": [{"id": "invite789", "email": "invite@test.com"}]
            },  # list_organization_invites
        ]

        # Create organization
        org_result = create_organization("Test Organization", description="Test org")
        assert org_result["id"] == "org123"

        # Get the created organization
        org_data = get_organization("org123")
        assert org_data["name"] == "Test Organization"

        # Add user to organization
        user_result = add_user_to_organization("org123", "user@test.com")
        assert user_result["user_email"] == "user@test.com"

        # List organization users
        users = list_organization_users("org123")
        assert len(users["results"]) == 1

        # Send organization invite
        invite_result = send_organization_invite("org123", "invite@test.com")
        assert invite_result["email"] == "invite@test.com"

        # List organization invites
        invites = list_organization_invites("org123")
        assert len(invites["results"]) == 1
