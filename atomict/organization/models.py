from dataclasses import dataclass
from datetime import datetime
from typing import Any, Dict, List, Optional


@dataclass
class Organization:
    """Organization data model"""

    id: str
    name: str
    description: Optional[str] = None
    billing_emails: Optional[List[str]] = None
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Organization":
        """Create Organization from API response"""
        return cls(
            id=data.get("id", ""),
            name=data.get("name", ""),
            description=data.get("description"),
            billing_emails=data.get("billing_emails"),
            created_at=data.get("created_at"),
            updated_at=data.get("updated_at"),
        )


@dataclass
class OrganizationUser:
    """Organization user data model"""

    id: str
    user_id: str
    organization_id: str
    is_admin: bool
    user_email: Optional[str] = None
    created_at: Optional[datetime] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "OrganizationUser":
        """Create OrganizationUser from API response"""
        return cls(
            id=data.get("id", ""),
            user_id=data.get("user_id", ""),
            organization_id=data.get("organization_id", ""),
            is_admin=data.get("is_admin", False),
            user_email=data.get("user_email"),
            created_at=data.get("created_at"),
        )


@dataclass
class OrganizationInvite:
    """Organization invite data model"""

    id: str
    organization_id: str
    email: str
    is_admin: bool
    created_at: Optional[datetime] = None
    expires_at: Optional[datetime] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "OrganizationInvite":
        """Create OrganizationInvite from API response"""
        return cls(
            id=data.get("id", ""),
            organization_id=data.get("organization_id", ""),
            email=data.get("email", ""),
            is_admin=data.get("is_admin", False),
            created_at=data.get("created_at"),
            expires_at=data.get("expires_at"),
        )
