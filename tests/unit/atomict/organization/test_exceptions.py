import pytest
from atomict.organization.exceptions import (
    OrganizationNotFoundError,
    OrganizationUserNotFoundError,
    OrganizationInviteNotFoundError,
    OrganizationPermissionError,
    InvalidOrganizationDataError,
)
from atomict.exceptions import APIValidationError, PermissionDenied


class TestOrganizationNotFoundError:
    """Test OrganizationNotFoundError exception"""

    def test_exception_inheritance(self):
        """Test that OrganizationNotFoundError inherits from APIValidationError"""
        assert issubclass(OrganizationNotFoundError, APIValidationError)

    def test_exception_creation_with_message(self):
        """Test creating exception with custom message"""
        message = "Organization 123 not found"
        error = OrganizationNotFoundError(message)
        
        assert str(error) == message
        assert isinstance(error, APIValidationError)

    def test_exception_creation_without_message(self):
        """Test creating exception without message"""
        error = OrganizationNotFoundError()
        
        assert isinstance(error, APIValidationError)

    def test_exception_raising(self):
        """Test raising OrganizationNotFoundError"""
        with pytest.raises(OrganizationNotFoundError) as exc_info:
            raise OrganizationNotFoundError("Test organization not found")
        
        assert str(exc_info.value) == "Test organization not found"

    def test_exception_catching_as_parent(self):
        """Test that exception can be caught as parent class"""
        with pytest.raises(APIValidationError):
            raise OrganizationNotFoundError("Test message")


class TestOrganizationUserNotFoundError:
    """Test OrganizationUserNotFoundError exception"""

    def test_exception_inheritance(self):
        """Test that OrganizationUserNotFoundError inherits from APIValidationError"""
        assert issubclass(OrganizationUserNotFoundError, APIValidationError)

    def test_exception_creation_with_message(self):
        """Test creating exception with custom message"""
        message = "Organization user 456 not found"
        error = OrganizationUserNotFoundError(message)
        
        assert str(error) == message
        assert isinstance(error, APIValidationError)

    def test_exception_creation_without_message(self):
        """Test creating exception without message"""
        error = OrganizationUserNotFoundError()
        
        assert isinstance(error, APIValidationError)

    def test_exception_raising(self):
        """Test raising OrganizationUserNotFoundError"""
        with pytest.raises(OrganizationUserNotFoundError) as exc_info:
            raise OrganizationUserNotFoundError("Test user not found")
        
        assert str(exc_info.value) == "Test user not found"

    def test_exception_catching_as_parent(self):
        """Test that exception can be caught as parent class"""
        with pytest.raises(APIValidationError):
            raise OrganizationUserNotFoundError("Test message")


class TestOrganizationInviteNotFoundError:
    """Test OrganizationInviteNotFoundError exception"""

    def test_exception_inheritance(self):
        """Test that OrganizationInviteNotFoundError inherits from APIValidationError"""
        assert issubclass(OrganizationInviteNotFoundError, APIValidationError)

    def test_exception_creation_with_message(self):
        """Test creating exception with custom message"""
        message = "Organization invite 789 not found"
        error = OrganizationInviteNotFoundError(message)
        
        assert str(error) == message
        assert isinstance(error, APIValidationError)

    def test_exception_creation_without_message(self):
        """Test creating exception without message"""
        error = OrganizationInviteNotFoundError()
        
        assert isinstance(error, APIValidationError)

    def test_exception_raising(self):
        """Test raising OrganizationInviteNotFoundError"""
        with pytest.raises(OrganizationInviteNotFoundError) as exc_info:
            raise OrganizationInviteNotFoundError("Test invite not found")
        
        assert str(exc_info.value) == "Test invite not found"

    def test_exception_catching_as_parent(self):
        """Test that exception can be caught as parent class"""
        with pytest.raises(APIValidationError):
            raise OrganizationInviteNotFoundError("Test message")


class TestOrganizationPermissionError:
    """Test OrganizationPermissionError exception"""

    def test_exception_inheritance(self):
        """Test that OrganizationPermissionError inherits from PermissionDenied"""
        assert issubclass(OrganizationPermissionError, PermissionDenied)

    def test_exception_creation_with_message(self):
        """Test creating exception with custom message"""
        message = "Access denied to organization 123"
        error = OrganizationPermissionError(message)
        
        assert str(error) == message
        assert isinstance(error, PermissionDenied)

    def test_exception_creation_without_message(self):
        """Test creating exception without message"""
        error = OrganizationPermissionError()
        
        assert isinstance(error, PermissionDenied)

    def test_exception_raising(self):
        """Test raising OrganizationPermissionError"""
        with pytest.raises(OrganizationPermissionError) as exc_info:
            raise OrganizationPermissionError("Test permission denied")
        
        assert str(exc_info.value) == "Test permission denied"

    def test_exception_catching_as_parent(self):
        """Test that exception can be caught as parent class"""
        with pytest.raises(PermissionDenied):
            raise OrganizationPermissionError("Test message")


class TestInvalidOrganizationDataError:
    """Test InvalidOrganizationDataError exception"""

    def test_exception_inheritance(self):
        """Test that InvalidOrganizationDataError inherits from APIValidationError"""
        assert issubclass(InvalidOrganizationDataError, APIValidationError)

    def test_exception_creation_with_message(self):
        """Test creating exception with custom message"""
        message = "Organization name is required"
        error = InvalidOrganizationDataError(message)
        
        assert str(error) == message
        assert isinstance(error, APIValidationError)

    def test_exception_creation_without_message(self):
        """Test creating exception without message"""
        error = InvalidOrganizationDataError()
        
        assert isinstance(error, APIValidationError)

    def test_exception_raising(self):
        """Test raising InvalidOrganizationDataError"""
        with pytest.raises(InvalidOrganizationDataError) as exc_info:
            raise InvalidOrganizationDataError("Test data validation error")
        
        assert str(exc_info.value) == "Test data validation error"

    def test_exception_catching_as_parent(self):
        """Test that exception can be caught as parent class"""
        with pytest.raises(APIValidationError):
            raise InvalidOrganizationDataError("Test message")


class TestExceptionHierarchy:
    """Test the overall exception hierarchy and relationships"""

    def test_all_exceptions_inherit_from_base_classes(self):
        """Test that all organization exceptions inherit from appropriate base classes"""
        # APIValidationError descendants
        api_validation_exceptions = [
            OrganizationNotFoundError,
            OrganizationUserNotFoundError,
            OrganizationInviteNotFoundError,
            InvalidOrganizationDataError,
        ]
        
        for exc_class in api_validation_exceptions:
            assert issubclass(exc_class, APIValidationError)
        
        # PermissionDenied descendants
        permission_exceptions = [
            OrganizationPermissionError,
        ]
        
        for exc_class in permission_exceptions:
            assert issubclass(exc_class, PermissionDenied)

    def test_exception_instances_are_unique(self):
        """Test that different exception instances are distinct"""
        exc1 = OrganizationNotFoundError("Test 1")
        exc2 = OrganizationNotFoundError("Test 2")
        exc3 = OrganizationUserNotFoundError("Test 1")
        
        assert exc1 != exc2
        assert exc1 != exc3
        assert exc2 != exc3

    def test_exception_chaining(self):
        """Test exception chaining with 'from' clause"""
        original_error = APIValidationError("Original error")
        
        with pytest.raises(OrganizationNotFoundError) as exc_info:
            try:
                raise original_error
            except APIValidationError as e:
                raise OrganizationNotFoundError("Organization not found") from e
        
        assert exc_info.value.__cause__ is original_error

    def test_multiple_exception_catching(self):
        """Test catching multiple exception types"""
        exceptions_to_test = [
            OrganizationNotFoundError("Not found"),
            OrganizationUserNotFoundError("User not found"),
            OrganizationInviteNotFoundError("Invite not found"),
            InvalidOrganizationDataError("Invalid data"),
        ]
        
        for exception in exceptions_to_test:
            with pytest.raises((
                OrganizationNotFoundError,
                OrganizationUserNotFoundError,
                OrganizationInviteNotFoundError,
                InvalidOrganizationDataError,
            )):
                raise exception

    def test_permission_error_distinction(self):
        """Test that permission errors are distinct from validation errors"""
        validation_error = InvalidOrganizationDataError("Invalid data")
        permission_error = OrganizationPermissionError("Access denied")
        
        assert isinstance(validation_error, APIValidationError)
        assert not isinstance(validation_error, PermissionDenied)
        
        assert isinstance(permission_error, PermissionDenied)
        assert not isinstance(permission_error, APIValidationError)


class TestExceptionMessages:
    """Test exception message handling and formatting"""

    def test_exception_with_formatted_messages(self):
        """Test exceptions with formatted string messages"""
        org_id = "123"
        user_id = "456"
        invite_id = "789"
        
        exceptions = [
            (OrganizationNotFoundError, f"Organization {org_id} not found"),
            (OrganizationUserNotFoundError, f"Organization user {user_id} not found"),
            (OrganizationInviteNotFoundError, f"Organization invite {invite_id} not found"),
            (OrganizationPermissionError, f"Access denied to organization {org_id}"),
            (InvalidOrganizationDataError, f"Invalid data for organization {org_id}"),
        ]
        
        for exc_class, message in exceptions:
            error = exc_class(message)
            assert str(error) == message

    def test_exception_with_empty_messages(self):
        """Test exceptions with empty messages"""
        exceptions = [
            OrganizationNotFoundError,
            OrganizationUserNotFoundError,
            OrganizationInviteNotFoundError,
            OrganizationPermissionError,
            InvalidOrganizationDataError,
        ]
        
        for exc_class in exceptions:
            error = exc_class("")
            assert str(error) == ""

    def test_exception_with_none_messages(self):
        """Test exceptions with None messages"""
        exceptions = [
            OrganizationNotFoundError,
            OrganizationUserNotFoundError,
            OrganizationInviteNotFoundError,
            OrganizationPermissionError,
            InvalidOrganizationDataError,
        ]
        
        for exc_class in exceptions:
            # Most exception classes handle None by converting to empty string or default
            error = exc_class(None)
            # The exact behavior depends on the base exception class implementation
            assert isinstance(error, exc_class)
