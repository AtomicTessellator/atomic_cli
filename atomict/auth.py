from atomict.api import post
from atomict.env import store


def authenticate(username: str, password: str) -> str:
    """Authenticate a user."""

    response = post("api-auth/", {"username": username, "password": password})

    store("token", response["token"])
    return response["token"]
