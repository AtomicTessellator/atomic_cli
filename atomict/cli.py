import os

import dotenv
import fire
import requests

from atomict.auth import authenticate
from atomict.env import check_environment, clear
from atomict.exceptions import APIValidationError, PermissionDenied


class AtomicT:
    """Atomic Tessellator CLI."""

    def auth(self):
        authenticate("alain@atomictessellator.com", "1181020")

    def clear(self):
        clear()


def main():
    dotenv.load_dotenv()
    check_environment()

    try:
        fire.Fire(AtomicT, name="atomict")
    except requests.exceptions.ConnectionError as e:
        print(
            f"{e}\n\nCould not connect to the Atomic Tessellator API. atomict is configured to use {os.environ.get('ATOMICT_API_ROOT')}"
        )
    except PermissionDenied as e:
        print(f"{e}\n\nPermission denied.")
    except APIValidationError as e:
        print(f"{e}")


if __name__ == "__main__":
    main()
