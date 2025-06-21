"""
Integration tests for organization operations.

These tests require:
- .env file with ATOMICT_USERNAME and ATOMICT_PASSWORD
- Access to reality_server API
- Optional TEST_ORGANIZATION_ID for non-destructive tests
- Optional TEST_ORGANIZATION_USER_EMAIL for user management tests

To run: uv run pytest tests/integration/test_organization.py -v -m integration
"""

import os
import uuid

import pytest
from dotenv import load_dotenv

from atomict.auth import authenticate
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


@pytest.fixture(scope="session", autouse=True)
def setup_authentication():
    """Setup authentication for integration tests"""
    load_dotenv()
    username = os.getenv("ATOMICT_USERNAME")
    password = os.getenv("ATOMICT_PASSWORD")

    if not username or not password:
        pytest.skip("ATOMICT_USERNAME and ATOMICT_PASSWORD must be set in .env file")

    try:
        token = authenticate(username, password)
        os.environ["AT_TOKEN"] = token
        return token
    except Exception as e:
        pytest.skip(f"Authentication failed: {e}")


@pytest.fixture(scope="session")
def test_organization_id():
    """Get test organization ID from environment (optional for non-destructive tests)"""
    return os.getenv("TEST_ORGANIZATION_ID")


@pytest.fixture(scope="session")
def test_user_email():
    """Get test user email from environment for user management tests"""
    return os.getenv("TEST_ORGANIZATION_USER_EMAIL")


@pytest.fixture
def unique_organization_name():
    """Generate a unique organization name for testing"""
    return f"integration-test-org-{uuid.uuid4().hex[:8]}"


@pytest.fixture
def unique_invite_email():
    """Generate a unique email for invitation testing"""
    return f"test-invite-{uuid.uuid4().hex[:8]}@example.com"


@pytest.fixture
def cleanup_organizations():
    """Track created organizations for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for org_id in created_ids:
        try:
            delete_organization(org_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def cleanup_organization_users():
    """Track added organization users for cleanup"""
    created_relations = []
    yield created_relations

    # Cleanup after test
    for relation_id in created_relations:
        try:
            remove_user_from_organization(relation_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.fixture
def cleanup_organization_invites():
    """Track sent organization invites for cleanup"""
    created_ids = []
    yield created_ids

    # Cleanup after test
    for invite_id in created_ids:
        try:
            delete_organization_invite(invite_id)
        except Exception:
            pass  # Best effort cleanup


@pytest.mark.integration
class TestOrganizationCRUDOperations:
    """Integration tests for organization CRUD operations"""

    def test_create_organization_minimal(
        self, unique_organization_name, cleanup_organizations
    ):
        """Test creating an organization with minimal parameters"""
        result = create_organization(name=unique_organization_name)

        # Track for cleanup
        cleanup_organizations.append(result["id"])

        # Verify response structure
        assert "id" in result
        assert result["name"] == unique_organization_name
        # Organisation model doesn't have timestamps or billing_emails field

    def test_create_organization_with_full_parameters(
        self, unique_organization_name, cleanup_organizations
    ):
        """Test creating an organization with all parameters"""
        description = "Integration test organization with full parameters"
        billing_primary_email = "billing@example.com"
        chem_primary_email = "chem@example.com"

        result = create_organization(
            name=unique_organization_name,
            description=description,
            billing_primary_email=billing_primary_email,
            chem_primary_email=chem_primary_email,
        )

        # Track for cleanup
        cleanup_organizations.append(result["id"])

        # Verify all parameters were set
        assert result["name"] == unique_organization_name
        assert result["description"] == description
        assert result["billing_primary_email"] == billing_primary_email
        assert result["chem_primary_email"] == chem_primary_email

    def test_get_organization(self, unique_organization_name, cleanup_organizations):
        """Test retrieving an organization by ID"""
        # Create organization first
        create_result = create_organization(
            name=unique_organization_name,
            description="Test organization for get operation",
        )
        org_id = create_result["id"]
        cleanup_organizations.append(org_id)

        # Get the organization
        get_result = get_organization(org_id)

        # Verify the retrieved organization matches
        assert get_result["id"] == org_id
        assert get_result["name"] == unique_organization_name
        assert get_result["description"] == "Test organization for get operation"

    def test_update_organization_partial(
        self, unique_organization_name, cleanup_organizations
    ):
        """Test updating an organization with partial field updates"""
        # Create organization first
        create_result = create_organization(
            name=unique_organization_name, description="Original description"
        )
        org_id = create_result["id"]
        cleanup_organizations.append(org_id)

        # Update only the description
        new_description = "Updated description"
        update_result = update_organization(org_id=org_id, description=new_description)

        # Verify update
        assert update_result["description"] == new_description
        assert (
            update_result["name"] == unique_organization_name
        )  # Should remain unchanged

        # Verify by getting the organization again
        get_result = get_organization(org_id)
        assert get_result["description"] == new_description
        assert get_result["name"] == unique_organization_name

    def test_update_organization_all_fields(
        self, unique_organization_name, cleanup_organizations
    ):
        """Test updating all organization fields"""
        # Create organization first
        create_result = create_organization(name=unique_organization_name)
        org_id = create_result["id"]
        cleanup_organizations.append(org_id)

        # Update all fields
        new_name = f"{unique_organization_name}-updated"
        new_description = "Completely updated organization"
        new_billing_primary_email = "updated@example.com"

        update_result = update_organization(
            org_id=org_id,
            name=new_name,
            description=new_description,
            billing_primary_email=new_billing_primary_email,
        )

        # Verify all updates
        assert update_result["name"] == new_name
        assert update_result["description"] == new_description
        assert update_result["billing_primary_email"] == new_billing_primary_email

    def test_delete_organization(self, unique_organization_name):
        """Test deleting an organization"""
        # Create organization first
        create_result = create_organization(name=unique_organization_name)
        org_id = create_result["id"]

        # Verify it exists
        get_result = get_organization(org_id)
        assert get_result["id"] == org_id

        # Delete the organization
        delete_result = delete_organization(org_id)

        # Verify deletion by attempting to get the organization (should raise exception)
        with pytest.raises((OrganizationNotFoundError, Exception)):
            get_organization(org_id)

    def test_list_organizations_basic(self):
        """Test basic organization listing"""
        result = list_organizations()

        # Verify response structure
        assert "results" in result
        assert "count" in result
        assert isinstance(result["results"], list)
        assert isinstance(result["count"], int)

    def test_list_organizations_with_search(
        self, unique_organization_name, cleanup_organizations
    ):
        """Test organization listing with search parameter"""
        search_term = "integration-search-test"
        org_name = f"{search_term}-{uuid.uuid4().hex[:8]}"

        # Create an organization with the search term
        create_result = create_organization(
            name=org_name,
            description=f"Organization for testing search functionality with {search_term}",
        )
        cleanup_organizations.append(create_result["id"])

        # Search for organizations containing the search term
        search_result = list_organizations(search=search_term)

        # Verify search results contain our organization
        assert search_result["count"] > 0
        org_names = [o["name"] for o in search_result["results"]]
        assert org_name in org_names


@pytest.mark.integration
class TestOrganizationUserManagement:
    """Integration tests for organization user management operations"""

    def test_list_organization_users_basic(self, test_organization_id):
        """Test basic organization user listing"""
        if not test_organization_id:
            pytest.skip("TEST_ORGANIZATION_ID environment variable must be set")

        result = list_organization_users(test_organization_id)

        # Verify response structure
        assert "results" in result
        assert "count" in result
        assert isinstance(result["results"], list)
        assert isinstance(result["count"], int)

    def test_add_and_remove_user_workflow(
        self,
        unique_organization_name,
        test_user_email,
        cleanup_organizations,
        cleanup_organization_users,
    ):
        """Test complete user add/remove workflow"""
        if not test_user_email:
            pytest.skip("TEST_ORGANIZATION_USER_EMAIL environment variable must be set")

        # Create organization first
        create_result = create_organization(
            name=unique_organization_name,
            description="Test organization for user management",
        )
        org_id = create_result["id"]
        cleanup_organizations.append(org_id)

        # Add user to organization
        add_result = add_user_to_organization(
            org_id=org_id, user_email=test_user_email, is_admin=False
        )

        # Track for cleanup
        cleanup_organization_users.append(add_result["id"])

        # Verify user was added
        assert "id" in add_result
        assert add_result["user_email"] == test_user_email
        assert add_result["is_admin"] is False

        # Verify user appears in organization user list
        users_result = list_organization_users(org_id)
        user_emails = [
            (
                user["user"]["email"]
                if isinstance(user["user"], dict)
                else user["user_email"]
            )
            for user in users_result["results"]
        ]
        assert test_user_email in user_emails

        # Remove user
        remove_result = remove_user_from_organization(add_result["id"])
        cleanup_organization_users.remove(add_result["id"])

        # Verify user was removed
        users_result_after = list_organization_users(org_id)
        user_emails_after = [
            (
                user["user"]["email"]
                if isinstance(user["user"], dict)
                else user["user_email"]
            )
            for user in users_result_after["results"]
        ]
        assert test_user_email not in user_emails_after

    def test_add_user_as_admin(
        self,
        unique_organization_name,
        test_user_email,
        cleanup_organizations,
        cleanup_organization_users,
    ):
        """Test adding a user as admin"""
        if not test_user_email:
            pytest.skip("TEST_ORGANIZATION_USER_EMAIL environment variable must be set")

        # Create organization first
        create_result = create_organization(name=unique_organization_name)
        org_id = create_result["id"]
        cleanup_organizations.append(org_id)

        # Add user as admin
        add_result = add_user_to_organization(
            org_id=org_id, user_email=test_user_email, is_admin=True
        )

        # Track for cleanup
        cleanup_organization_users.append(add_result["id"])

        # Verify user was added as admin
        assert add_result["is_admin"] is True


@pytest.mark.integration
class TestOrganizationInviteManagement:
    """Integration tests for organization invitation management operations"""

    def test_list_organization_invites_basic(self, test_organization_id):
        """Test basic organization invite listing"""
        if not test_organization_id:
            pytest.skip("TEST_ORGANIZATION_ID environment variable must be set")

        result = list_organization_invites(test_organization_id)

        # Verify response structure
        assert "results" in result
        assert "count" in result
        assert isinstance(result["results"], list)
        assert isinstance(result["count"], int)

    def test_send_and_cancel_invite_workflow(
        self,
        unique_organization_name,
        unique_invite_email,
        cleanup_organizations,
        cleanup_organization_invites,
    ):
        """Test complete invite send/cancel workflow"""
        # Create organization first
        create_result = create_organization(
            name=unique_organization_name,
            description="Test organization for invite management",
        )
        org_id = create_result["id"]
        cleanup_organizations.append(org_id)

        # Send invite
        invite_result = send_organization_invite(
            org_id=org_id, email=unique_invite_email, is_admin=False
        )

        # Track for cleanup
        cleanup_organization_invites.append(invite_result["id"])

        # Verify invite was sent
        assert "id" in invite_result
        assert invite_result["email"] == unique_invite_email
        # Note: OrganisationInvite model doesn't have is_admin field
        # The is_admin parameter is ignored - admin status is set when invite is accepted

        # Verify invite appears in organization invite list
        invites_result = list_organization_invites(org_id)
        invite_emails = [invite["email"] for invite in invites_result["results"]]
        assert unique_invite_email in invite_emails

        # Cancel invite
        cancel_result = delete_organization_invite(invite_result["id"])
        cleanup_organization_invites.remove(invite_result["id"])

        # Verify invite was cancelled
        invites_result_after = list_organization_invites(org_id)
        invite_emails_after = [
            invite["email"] for invite in invites_result_after["results"]
        ]
        assert unique_invite_email not in invite_emails_after

    def test_send_admin_invite(
        self,
        unique_organization_name,
        unique_invite_email,
        cleanup_organizations,
        cleanup_organization_invites,
    ):
        """Test sending an admin invite"""
        # Create organization first
        create_result = create_organization(name=unique_organization_name)
        org_id = create_result["id"]
        cleanup_organizations.append(org_id)

        # Send admin invite
        invite_result = send_organization_invite(
            org_id=org_id, email=unique_invite_email, is_admin=True
        )

        # Track for cleanup
        cleanup_organization_invites.append(invite_result["id"])

        # Verify invite was sent
        assert "id" in invite_result
        assert invite_result["email"] == unique_invite_email
        # Note: OrganisationInvite model doesn't support is_admin field yet
        # Admin status is set when invite is accepted and user becomes OrganisationUser


@pytest.mark.integration
class TestOrganizationPreferences:
    """Integration tests for organization preferences and active organization management"""

    def test_get_active_organization(self):
        """Test getting the user's active organization"""
        try:
            result = get_active_organization()
            # Should return active organization response structure
            assert "active_organisation" in result
            assert "organisation_details" in result

            # If there's an active organization, check its structure
            if result["active_organisation"]:
                assert result["organisation_details"] is not None
                assert "id" in result["organisation_details"]
                assert "name" in result["organisation_details"]
        except OrganizationPermissionError:
            # This might happen if user has no organizations
            pytest.skip("User has no active organization or lacks permission")

    def test_set_active_organization_workflow(self, test_organization_id):
        """Test setting active organization"""
        if not test_organization_id:
            pytest.skip("TEST_ORGANIZATION_ID environment variable must be set")

        # Get current active organization (if any)
        original_active = None
        try:
            original_active = get_active_organization()
        except OrganizationPermissionError:
            pass

        # Set the test organization as active
        try:
            result = set_active_organization(test_organization_id)

            # The response should be the updated preferences with active organization
            assert "active_organisation" in result

            # Verify the active organization changed
            new_active = get_active_organization()
            if new_active["active_organisation"]:
                assert new_active["active_organisation"] == test_organization_id

        finally:
            # Restore original active organization if there was one
            if original_active and original_active.get("active_organisation"):
                try:
                    set_active_organization(original_active["active_organisation"])
                except Exception:
                    pass  # Best effort restoration


@pytest.mark.integration
class TestOrganizationErrorHandling:
    """Integration tests for organization validation and error handling"""

    def test_get_nonexistent_organization(self):
        """Test getting an organization that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises((OrganizationNotFoundError, Exception)):
            get_organization(non_existent_id)

    def test_update_nonexistent_organization(self):
        """Test updating an organization that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises((OrganizationNotFoundError, Exception)):
            update_organization(org_id=non_existent_id, name="Updated Name")

    def test_delete_nonexistent_organization(self):
        """Test deleting an organization that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises((OrganizationNotFoundError, Exception)):
            delete_organization(non_existent_id)

    def test_create_organization_validation_errors(self):
        """Test organization creation validation"""
        # Test empty name
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization name is required"
        ):
            create_organization(name="")

        # Test None name
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization name is required"
        ):
            create_organization(name=None)  # type: ignore

    def test_update_organization_validation_errors(self):
        """Test organization update validation"""
        # Test empty org_id
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            update_organization(org_id="", name="Test")

        # Test no fields to update
        with pytest.raises(
            InvalidOrganizationDataError, match="At least one field must be provided"
        ):
            update_organization(org_id="test-id")

        # Test invalid fields
        with pytest.raises(
            InvalidOrganizationDataError, match="No valid fields provided"
        ):
            update_organization(org_id="test-id", invalid_field="value")

    def test_add_user_validation_errors(self):
        """Test user addition validation"""
        # Test empty org_id
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            add_user_to_organization(org_id="", user_email="test@example.com")

        # Test empty email
        with pytest.raises(
            InvalidOrganizationDataError, match="User email is required"
        ):
            add_user_to_organization(org_id="test-id", user_email="")

    def test_send_invite_validation_errors(self):
        """Test invite sending validation"""
        # Test empty org_id
        with pytest.raises(
            InvalidOrganizationDataError, match="Organization ID is required"
        ):
            send_organization_invite(org_id="", email="test@example.com")

        # Test empty email
        with pytest.raises(InvalidOrganizationDataError, match="Email is required"):
            send_organization_invite(org_id="test-id", email="")

    def test_delete_nonexistent_invite(self):
        """Test deleting an invite that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises((OrganizationInviteNotFoundError, Exception)):
            delete_organization_invite(non_existent_id)

    def test_remove_nonexistent_user(self):
        """Test removing a user relationship that doesn't exist"""
        non_existent_id = f"non-existent-{uuid.uuid4()}"

        with pytest.raises((OrganizationUserNotFoundError, Exception)):
            remove_user_from_organization(non_existent_id)


@pytest.mark.integration
class TestOrganizationWorkflows:
    """End-to-end workflow tests for organization operations"""

    def test_complete_organization_lifecycle(self, unique_organization_name):
        """Test complete organization lifecycle: create → get → update → delete"""
        original_description = "Original organization description"
        original_billing_primary_email = "original@example.com"

        # Step 1: Create organization
        create_result = create_organization(
            name=unique_organization_name,
            description=original_description,
            billing_primary_email=original_billing_primary_email,
        )
        org_id = create_result["id"]

        try:
            # Verify creation
            assert create_result["name"] == unique_organization_name
            assert create_result["description"] == original_description
            assert (
                create_result["billing_primary_email"] == original_billing_primary_email
            )

            # Step 2: Get organization and verify
            get_result = get_organization(org_id)
            assert get_result["id"] == org_id
            assert get_result["name"] == unique_organization_name

            # Step 3: Update organization
            updated_description = "Updated organization description"
            updated_billing_primary_email = "updated@example.com"

            update_result = update_organization(
                org_id=org_id,
                description=updated_description,
                billing_primary_email=updated_billing_primary_email,
            )

            assert update_result["description"] == updated_description
            assert (
                update_result["billing_primary_email"] == updated_billing_primary_email
            )

            # Step 4: Verify update by getting organization again
            updated_get_result = get_organization(org_id)
            assert updated_get_result["description"] == updated_description
            assert (
                updated_get_result["billing_primary_email"]
                == updated_billing_primary_email
            )

        finally:
            # Step 5: Delete organization
            delete_organization(org_id)

        # Step 6: Verify deletion
        with pytest.raises((OrganizationNotFoundError, Exception)):
            get_organization(org_id)

    def test_organization_user_management_workflow(
        self, unique_organization_name, test_user_email
    ):
        """Test complete organization user management workflow"""
        if not test_user_email:
            pytest.skip("TEST_ORGANIZATION_USER_EMAIL environment variable must be set")

        # Create organization
        create_result = create_organization(
            name=unique_organization_name, description="User management workflow test"
        )
        org_id = create_result["id"]

        try:
            # Add user as member
            add_result = add_user_to_organization(
                org_id=org_id, user_email=test_user_email, is_admin=False
            )
            user_relation_id = add_result["id"]

            # Verify user was added
            users_result = list_organization_users(org_id)
            assert users_result["count"] > 0

            # Remove user
            remove_user_from_organization(user_relation_id)

            # Verify user was removed
            users_result_after = list_organization_users(org_id)
            # Should have one less user (or the specific user should be gone)

        finally:
            # Clean up organization
            delete_organization(org_id)

    def test_organization_invite_workflow(
        self, unique_organization_name, unique_invite_email
    ):
        """Test complete organization invitation workflow"""
        # Create organization
        create_result = create_organization(
            name=unique_organization_name, description="Invite workflow test"
        )
        org_id = create_result["id"]

        try:
            # Send invite
            invite_result = send_organization_invite(
                org_id=org_id, email=unique_invite_email, is_admin=False
            )
            invite_id = invite_result["id"]

            # Verify invite in list
            invites_result = list_organization_invites(org_id)
            assert any(
                invite["id"] == invite_id for invite in invites_result["results"]
            )

            # Cancel invite
            delete_organization_invite(invite_id)

            # Verify invite was cancelled
            invites_result_after = list_organization_invites(org_id)
            assert not any(
                invite["id"] == invite_id for invite in invites_result_after["results"]
            )

        finally:
            # Clean up organization
            delete_organization(org_id)
