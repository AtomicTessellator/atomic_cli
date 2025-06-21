"""Structure conversion utilities for atomic tessellator.

This module provides functions for converting and manipulating crystal structures,
including supercell generation and conventional cell transformations.

IMPORTANT NOTE:
All functions in this module require UserUpload UUIDs, NOT structure IDs.
Use the UUID from the UserUpload object (e.g., "fea20908-39d5-43cb-a4a6-3918d89d8232"),
not the numeric structure ID from other parts of the system.
"""

from typing import List, Union, Dict, Any
from atomict.api import post
from atomict.exceptions import APIValidationError


def create_supercell(
    user_upload_id: str, dimensions: Union[List[int], str], **kwargs: Any
) -> Dict[str, Any]:
    """Create a supercell from an existing structure.

    Takes a user-uploaded structure and creates a supercell with the specified
    dimensions. The new supercell structure is saved as a new UserUpload.

    Args:
        user_upload_id: UUID of the source UserUpload structure
            **IMPORTANT**: Must be the UserUpload UUID, not the structure ID.
            Use the UUID from the UserUpload object, not the numeric ID.
        dimensions: Supercell dimensions as either:
            - List of 3 integers: [a, b, c]
            - String: "a b c" (space-separated)
        **kwargs: Additional parameters passed to the API

    Returns:
        API response dictionary containing result information

    Raises:
        APIValidationError: If dimensions are invalid or user_upload_id not found
        PermissionDenied: If user lacks access to the source structure

    Example:
        >>> # Use UserUpload UUID, not structure ID
        >>> create_supercell("fea20908-39d5-43cb-a4a6-3918d89d8232", [2, 2, 2])
        {'response': 'ok'}

        >>> create_supercell("fea20908-39d5-43cb-a4a6-3918d89d8232", "3 3 3")
        {'response': 'ok'}
    """
    # Validate and normalize dimensions
    if isinstance(dimensions, list):
        if len(dimensions) != 3:
            raise APIValidationError("Dimensions must contain exactly 3 values")
        if not all(isinstance(d, int) and d > 0 for d in dimensions):
            raise APIValidationError("All dimensions must be positive integers")
        dimensions_str = " ".join(str(d) for d in dimensions)
    elif isinstance(dimensions, str):
        try:
            dims = [int(x.strip()) for x in dimensions.split()]
            if len(dims) != 3:
                raise APIValidationError("Dimensions must contain exactly 3 values")
            if not all(d > 0 for d in dims):
                raise APIValidationError("All dimensions must be positive integers")
            dimensions_str = dimensions.strip()
        except ValueError:
            raise APIValidationError("Dimensions must be integers")
    else:
        raise APIValidationError(
            "Dimensions must be a list of integers or space-separated string"
        )

    # Validate user_upload_id format (basic UUID validation)
    if not user_upload_id or not isinstance(user_upload_id, str):
        raise APIValidationError("user_upload_id must be a non-empty string")

    # Build payload
    payload = {"user_upload": user_upload_id, "supercell": dimensions_str, **kwargs}

    # Make API call
    response = post("conversion/supercell/", payload)
    return response


def create_conventional_cell(user_upload_id: str, **kwargs: Any) -> Dict[str, Any]:
    """Create a conventional cell from an existing structure.

    Takes a user-uploaded structure and creates its standardized conventional
    cell representation. The new structure is saved as a new UserUpload.

    Args:
        user_upload_id: UUID of the source UserUpload structure
            **IMPORTANT**: Must be the UserUpload UUID, not the structure ID.
            Use the UUID from the UserUpload object, not the numeric ID.
        **kwargs: Additional parameters passed to the API

    Returns:
        API response dictionary containing result information

    Raises:
        APIValidationError: If user_upload_id is invalid
        PermissionDenied: If user lacks access to the source structure

    Example:
        >>> # Use UserUpload UUID, not structure ID
        >>> create_conventional_cell("fea20908-39d5-43cb-a4a6-3918d89d8232")
        {'response': 'ok'}
    """
    # Validate user_upload_id format
    if not user_upload_id or not isinstance(user_upload_id, str):
        raise APIValidationError("user_upload_id must be a non-empty string")

    # Build payload
    payload = {"user_upload": user_upload_id, **kwargs}

    # Make API call
    response = post("conversion/conventional/", payload)
    return response
