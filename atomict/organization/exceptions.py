from atomict.exceptions import APIValidationError, PermissionDenied


class OrganizationNotFoundError(APIValidationError):
    """Raised when an organization is not found"""

    pass


class OrganizationUserNotFoundError(APIValidationError):
    """Raised when an organization user is not found"""

    pass


class OrganizationInviteNotFoundError(APIValidationError):
    """Raised when an organization invite is not found"""

    pass


class OrganizationPermissionError(PermissionDenied):
    """Raised when user lacks permission for organization operation"""

    pass


class InvalidOrganizationDataError(APIValidationError):
    """Raised when organization data is invalid"""

    pass
