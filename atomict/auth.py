import os

from atomict.api import post
from atomict.env import get, store


def authenticate(username: str, password: str) -> str:
    """Authenticate a user."""

    response = post("api-auth/", {"username": username, "password": password})

    store("token", response["token"])
    return response["token"]


def resolve_token() -> str:
    """
        Resolve the token from the environment or authenticate.

        Order of presedence:
        - Environment variable
        - Store
        - Authenticate


    """

    token = os.environ.get("AT_TOKEN")
    if token:
        return token

    try:
        token = get("token")
        if token:
            return token
    except FileNotFoundError:
        pass

    user = os.environ.get("AT_USER")
    passw = os.environ.get("AT_PASS")

    if None in [user, passw]:
        raise ValueError(
            "AT_TOKEN not set and user and password also not set AT_USER, AT_PASS"
        )

    return authenticate(user, passw)
