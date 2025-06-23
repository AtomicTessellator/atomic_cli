from atomict.api import get, post
from atomict.exceptions import APIValidationError
from typing import Optional, Dict, Any

# Constants
DEFAULT_BATCH_SIZE = 16
DEFAULT_NUM_BATCHES = 1


def create_mattergen(
    project_id: str,
    name: str,
    description: str = "",
    batch_size: int = DEFAULT_BATCH_SIZE,
    num_batches: int = DEFAULT_NUM_BATCHES,
    diffusion_guidance_factor: Optional[int] = None,
    action: str = "DRAFT",
    extra_kwargs: Optional[Dict[str, Any]] = None,
):
    """
    Create a MatterGen exploration.

    Args:
        project_id: Project ID to create the exploration in
        name: Name of the exploration
        description: Description of the exploration
        batch_size: Number of samples per batch (default: 16)
        num_batches: Number of batches to generate (default: 1)
        diffusion_guidance_factor: Guidance factor for diffusion process (optional)
        action: Action to perform - "DRAFT" or "LAUNCH" (default: "DRAFT")
        extra_kwargs: Additional parameters (e.g., cluster configuration)

    Returns:
        API response containing the created exploration

    Raises:
        APIValidationError: If parameters are invalid
    """
    # Validate action
    if action not in ["DRAFT", "LAUNCH"]:
        raise APIValidationError(f"Action must be 'DRAFT' or 'LAUNCH', got '{action}'")

    # Validate batch parameters
    if batch_size <= 0:
        raise APIValidationError(f"batch_size must be positive, got {batch_size}")

    if num_batches <= 0:
        raise APIValidationError(f"num_batches must be positive, got {num_batches}")

    if diffusion_guidance_factor is not None and diffusion_guidance_factor <= 0:
        raise APIValidationError(
            f"diffusion_guidance_factor must be positive, got {diffusion_guidance_factor}"
        )

    # Build payload with field mapping
    payload = {
        "project": project_id,  # LaunchableViewSet expects "project" not "project_id"
        "name": name,
        "description": description,
        "batch_size": batch_size,
        "num_batches": num_batches,
        "action": action,
    }

    # Add optional parameters
    if diffusion_guidance_factor is not None:
        payload["diffusion_guidance_factor"] = diffusion_guidance_factor

    # Add extra parameters (e.g., cluster configuration)
    if extra_kwargs:
        payload.update(extra_kwargs)

    return post("api/mattergen-exploration/", payload=payload)


def get_mattergen(id: str):
    """
    Get MatterGen
    """
    result = get(f"api/mattergen-exploration/{id}/")
    return result


def associate_user_upload_with_mattergen(user_upload_id: str, exploration_id: str):
    """
    Associate a user upload with a MatterGen
    """
    result = post(
        "api/mattergen-exploration-file/",
        payload={"user_upload_id": user_upload_id, "exploration_id": exploration_id},
    )
    return result
