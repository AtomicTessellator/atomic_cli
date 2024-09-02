import os
from typing import List


def environment_check(required_env_keys: List[str] = []):
    """
    Check if the environment is set up correctly.

    Args:
        required_env_keys (List[str]): List of environment keys that are required.

    Raises:
        ValueError: If any of the required environment keys are not set.
        ValueError: If the AT_USER and AT_PASS or AT_TOKEN are not set.
        ValueError: If the AT_SERVER is not set.
    """
    expected_env_keys = ["AT_SERVER", *required_env_keys]

    for k in expected_env_keys:
        if os.environ.get(k) in ["", None]:
            raise ValueError(f"{k} is not set")

    auth_creds = [os.environ.get("AT_USER"), os.environ.get("AT_PASS")]
    if None in auth_creds:
        if os.environ.get("AT_TOKEN") in ["", None]:
            raise ValueError("No authentication credentials found, please set AT_USER and AT_PASS or AT_TOKEN")
