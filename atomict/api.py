import json
import os

import requests

from atomict.exceptions import APIValidationError, PermissionDenied


def get(path: str):
    api_root = os.environ.get("ATOMICT_API_ROOT")

    response = requests.get(f"{api_root}/{path}")

    if response.status_code == requests.codes.ok:
        resp = response.json()

        if "error" in resp:
            raise PermissionDenied(resp["error"])
        else:
            return resp

    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()


def post(path: str, payload: dict):
    payload_enc = json.dumps(payload)
    headers = {"Content-Type": "application/json"}
    api_root = os.environ.get("ATOMICT_API_ROOT")

    response = requests.post(f"{api_root}/{path}", data=payload_enc, headers=headers)

    if response.status_code == requests.codes.ok:
        resp = response.json()

        if "error" in resp:
            raise PermissionDenied(resp["error"])
        else:
            return resp

    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()
