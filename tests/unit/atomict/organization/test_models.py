from datetime import datetime

import pytest

from atomict.organization.models import (
    Organization,
    OrganizationInvite,
    OrganizationUser,
)


class TestOrganization:
    """Test Organization data model"""

    def test_organization_creation_minimal(self):
        """Test creating organization with minimal required fields"""
        org = Organization(id="123", name="Test Org")

        assert org.id == "123"
        assert org.name == "Test Org"
        assert org.description is None
        assert org.billing_emails is None
        assert org.created_at is None
        assert org.updated_at is None

    def test_organization_creation_full(self):
        """Test creating organization with all fields"""
        created_at = datetime.now()
        updated_at = datetime.now()
        billing_emails = ["billing@test.com", "admin@test.com"]

        org = Organization(
            id="123",
            name="Test Org",
            description="Test description",
            billing_emails=billing_emails,
            created_at=created_at,
            updated_at=updated_at,
        )

        assert org.id == "123"
        assert org.name == "Test Org"
        assert org.description == "Test description"
        assert org.billing_emails == billing_emails
        assert org.created_at == created_at
        assert org.updated_at == updated_at

    def test_organization_from_dict_minimal(self):
        """Test creating organization from dict with minimal data"""
        data = {"id": "123", "name": "Test Org"}

        org = Organization.from_dict(data)

        assert org.id == "123"
        assert org.name == "Test Org"
        assert org.description is None
        assert org.billing_emails is None
        assert org.created_at is None
        assert org.updated_at is None

    def test_organization_from_dict_full(self):
        """Test creating organization from dict with all fields"""
        created_at = "2023-01-01T00:00:00Z"
        updated_at = "2023-01-02T00:00:00Z"

        data = {
            "id": "123",
            "name": "Test Org",
            "description": "Test description",
            "billing_emails": ["billing@test.com", "admin@test.com"],
            "created_at": created_at,
            "updated_at": updated_at,
        }

        org = Organization.from_dict(data)

        assert org.id == "123"
        assert org.name == "Test Org"
        assert org.description == "Test description"
        assert org.billing_emails == ["billing@test.com", "admin@test.com"]
        assert org.created_at == created_at
        assert org.updated_at == updated_at

    def test_organization_from_dict_empty(self):
        """Test creating organization from empty dict"""
        data = {}

        org = Organization.from_dict(data)

        assert org.id == ""
        assert org.name == ""
        assert org.description is None
        assert org.billing_emails is None
        assert org.created_at is None
        assert org.updated_at is None

    def test_organization_from_dict_missing_fields(self):
        """Test creating organization from dict with missing fields"""
        data = {"id": "123", "description": "Test description"}

        org = Organization.from_dict(data)

        assert org.id == "123"
        assert org.name == ""  # Default for missing field
        assert org.description == "Test description"
        assert org.billing_emails is None
        assert org.created_at is None
        assert org.updated_at is None


class TestOrganizationUser:
    """Test OrganizationUser data model"""

    def test_organization_user_creation_minimal(self):
        """Test creating organization user with minimal required fields"""
        org_user = OrganizationUser(
            id="user123", user_id="456", organization_id="org789", is_admin=False
        )

        assert org_user.id == "user123"
        assert org_user.user_id == "456"
        assert org_user.organization_id == "org789"
        assert org_user.is_admin is False
        assert org_user.user_email is None
        assert org_user.created_at is None

    def test_organization_user_creation_full(self):
        """Test creating organization user with all fields"""
        created_at = datetime.now()

        org_user = OrganizationUser(
            id="user123",
            user_id="456",
            organization_id="org789",
            is_admin=True,
            user_email="user@test.com",
            created_at=created_at,
        )

        assert org_user.id == "user123"
        assert org_user.user_id == "456"
        assert org_user.organization_id == "org789"
        assert org_user.is_admin is True
        assert org_user.user_email == "user@test.com"
        assert org_user.created_at == created_at

    def test_organization_user_from_dict_minimal(self):
        """Test creating organization user from dict with minimal data"""
        data = {"id": "user123", "user_id": "456", "organization_id": "org789"}

        org_user = OrganizationUser.from_dict(data)

        assert org_user.id == "user123"
        assert org_user.user_id == "456"
        assert org_user.organization_id == "org789"
        assert org_user.is_admin is False  # Default value
        assert org_user.user_email is None
        assert org_user.created_at is None

    def test_organization_user_from_dict_full(self):
        """Test creating organization user from dict with all fields"""
        created_at = "2023-01-01T00:00:00Z"

        data = {
            "id": "user123",
            "user_id": "456",
            "organization_id": "org789",
            "is_admin": True,
            "user_email": "admin@test.com",
            "created_at": created_at,
        }

        org_user = OrganizationUser.from_dict(data)

        assert org_user.id == "user123"
        assert org_user.user_id == "456"
        assert org_user.organization_id == "org789"
        assert org_user.is_admin is True
        assert org_user.user_email == "admin@test.com"
        assert org_user.created_at == created_at

    def test_organization_user_from_dict_empty(self):
        """Test creating organization user from empty dict"""
        data = {}

        org_user = OrganizationUser.from_dict(data)

        assert org_user.id == ""
        assert org_user.user_id == ""
        assert org_user.organization_id == ""
        assert org_user.is_admin is False
        assert org_user.user_email is None
        assert org_user.created_at is None

    def test_organization_user_boolean_handling(self):
        """Test proper boolean handling for is_admin field"""
        # Test various truthy/falsy values
        test_cases = [
            (True, True),
            (False, False),
            (1, 1),  # Should preserve original value
            (0, 0),  # Should preserve original value
            ("true", "true"),  # Should preserve original value
            ("false", "false"),  # Should preserve original value
        ]

        for input_value, expected in test_cases:
            data = {
                "id": "user123",
                "user_id": "456",
                "organization_id": "org789",
                "is_admin": input_value,
            }

            org_user = OrganizationUser.from_dict(data)
            assert org_user.is_admin == expected


class TestOrganizationInvite:
    """Test OrganizationInvite data model"""

    def test_organization_invite_creation_minimal(self):
        """Test creating organization invite with minimal required fields"""
        invite = OrganizationInvite(
            id="invite123",
            organization_id="org456",
            email="invite@test.com",
            is_admin=False,
        )

        assert invite.id == "invite123"
        assert invite.organization_id == "org456"
        assert invite.email == "invite@test.com"
        assert invite.is_admin is False
        assert invite.created_at is None
        assert invite.expires_at is None

    def test_organization_invite_creation_full(self):
        """Test creating organization invite with all fields"""
        created_at = datetime.now()
        expires_at = datetime.now()

        invite = OrganizationInvite(
            id="invite123",
            organization_id="org456",
            email="admin@test.com",
            is_admin=True,
            created_at=created_at,
            expires_at=expires_at,
        )

        assert invite.id == "invite123"
        assert invite.organization_id == "org456"
        assert invite.email == "admin@test.com"
        assert invite.is_admin is True
        assert invite.created_at == created_at
        assert invite.expires_at == expires_at

    def test_organization_invite_from_dict_minimal(self):
        """Test creating organization invite from dict with minimal data"""
        data = {
            "id": "invite123",
            "organization_id": "org456",
            "email": "invite@test.com",
        }

        invite = OrganizationInvite.from_dict(data)

        assert invite.id == "invite123"
        assert invite.organization_id == "org456"
        assert invite.email == "invite@test.com"
        assert invite.is_admin is False  # Default value
        assert invite.created_at is None
        assert invite.expires_at is None

    def test_organization_invite_from_dict_full(self):
        """Test creating organization invite from dict with all fields"""
        created_at = "2023-01-01T00:00:00Z"
        expires_at = "2023-01-08T00:00:00Z"

        data = {
            "id": "invite123",
            "organization_id": "org456",
            "email": "admin@test.com",
            "is_admin": True,
            "created_at": created_at,
            "expires_at": expires_at,
        }

        invite = OrganizationInvite.from_dict(data)

        assert invite.id == "invite123"
        assert invite.organization_id == "org456"
        assert invite.email == "admin@test.com"
        assert invite.is_admin is True
        assert invite.created_at == created_at
        assert invite.expires_at == expires_at

    def test_organization_invite_from_dict_empty(self):
        """Test creating organization invite from empty dict"""
        data = {}

        invite = OrganizationInvite.from_dict(data)

        assert invite.id == ""
        assert invite.organization_id == ""
        assert invite.email == ""
        assert invite.is_admin is False
        assert invite.created_at is None
        assert invite.expires_at is None

    def test_organization_invite_email_validation(self):
        """Test that email field is preserved as-is (no validation in model)"""
        # Test various email formats - model should preserve them as-is
        test_emails = [
            "valid@test.com",
            "invalid-email",
            "",
            "multiple@emails@test.com",
            "user+tag@domain.co.uk",
        ]

        for email in test_emails:
            data = {"id": "invite123", "organization_id": "org456", "email": email}

            invite = OrganizationInvite.from_dict(data)
            assert invite.email == email


class TestModelsComparison:
    """Test model comparison and equality"""

    def test_organization_equality(self):
        """Test organization equality comparison"""
        org1 = Organization(id="123", name="Test Org")
        org2 = Organization(id="123", name="Test Org")
        org3 = Organization(id="456", name="Test Org")

        assert org1 == org2
        assert org1 != org3

    def test_organization_user_equality(self):
        """Test organization user equality comparison"""
        user1 = OrganizationUser(
            id="1", user_id="2", organization_id="3", is_admin=False
        )
        user2 = OrganizationUser(
            id="1", user_id="2", organization_id="3", is_admin=False
        )
        user3 = OrganizationUser(
            id="1", user_id="2", organization_id="3", is_admin=True
        )

        assert user1 == user2
        assert user1 != user3

    def test_organization_invite_equality(self):
        """Test organization invite equality comparison"""
        invite1 = OrganizationInvite(
            id="1", organization_id="2", email="test@example.com", is_admin=False
        )
        invite2 = OrganizationInvite(
            id="1", organization_id="2", email="test@example.com", is_admin=False
        )
        invite3 = OrganizationInvite(
            id="1", organization_id="2", email="other@example.com", is_admin=False
        )

        assert invite1 == invite2
        assert invite1 != invite3


class TestModelsStringRepresentation:
    """Test string representation of models"""

    def test_organization_repr(self):
        """Test organization string representation"""
        org = Organization(id="123", name="Test Org")
        repr_str = repr(org)

        assert "Organization" in repr_str
        assert "123" in repr_str
        assert "Test Org" in repr_str

    def test_organization_user_repr(self):
        """Test organization user string representation"""
        user = OrganizationUser(
            id="1", user_id="2", organization_id="3", is_admin=False
        )
        repr_str = repr(user)

        assert "OrganizationUser" in repr_str
        assert "1" in repr_str
        assert "2" in repr_str
        assert "3" in repr_str

    def test_organization_invite_repr(self):
        """Test organization invite string representation"""
        invite = OrganizationInvite(
            id="1", organization_id="2", email="test@example.com", is_admin=False
        )
        repr_str = repr(invite)

        assert "OrganizationInvite" in repr_str
        assert "1" in repr_str
        assert "2" in repr_str
        assert "test@example.com" in repr_str
