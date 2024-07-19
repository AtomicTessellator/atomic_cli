import json
import os

import requests

from atomict.exceptions import APIValidationError, PermissionDenied


def get(path: str):
    api_root = os.environ.get("AT_SERVER", "https://api.atomictessellator.com")
    headers = {"Accept": "application/json", "Content-Type": "application/json"}

    if os.environ.get("AT_TOKEN"):
        headers["Authorization"] = f"Token {os.environ.get('AT_TOKEN')}"

    response = requests.get(f"{api_root}/{path}", headers=headers)

    content_type = response.headers.get("Content-Type")

    if response.status_code == requests.codes.ok and content_type == "application/json":
        resp = response.json()

        if "error" in resp and resp["error"] is not None:
            raise PermissionDenied(resp["error"])
        else:
            return resp
    elif (
        response.status_code == requests.codes.ok and content_type != "application/json"
    ):
        return response.content
    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()


def post(path: str, payload: dict, files=None, extra_headers={}):
    # Jesus christ this logic needs cleaning up
    if not files and "Content-Type" not in extra_headers:
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
    else:
        headers = {}

    if extra_headers:
        if "Content-Type" in extra_headers:
            headers["Content-Type"] = extra_headers["Content-Type"]
            payload = json.dumps(payload)
        else:
            headers.update(extra_headers)

    if os.environ.get("AT_TOKEN"):
        headers["Authorization"] = f"Token {os.environ.get('AT_TOKEN')}"

    api_root = os.environ.get("AT_SERVER")

    if files is not None:
        response = requests.post(
            f"{api_root}/{path}", data=payload, headers=headers, files=files
        )
    else:
        response = requests.post(f"{api_root}/{path}", data=payload, headers=headers)

    if response.status_code in [requests.codes.ok, requests.codes.created]:
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


def patch(path: str, payload: dict):
    payload_enc = json.dumps(payload)
    headers = {"Content-Type": "application/json"}

    if os.environ.get("AT_TOKEN"):
        headers["Authorization"] = f"Token {os.environ.get('AT_TOKEN')}"

    api_root = os.environ.get("AT_SERVER")
    response = requests.patch(f"{api_root}/{path}", data=payload_enc, headers=headers)

    if response.status_code == requests.codes.ok:
        resp = response.json()

        if resp.get("error") is not None:
            raise Exception(resp["error"])
        else:
            return resp

    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()


def delete(path: str):
    headers = {}

    if os.environ.get("AT_TOKEN"):
        headers["Authorization"] = f"Token {os.environ.get('AT_TOKEN')}"

    api_root = os.environ.get("AT_SERVER")
    response = requests.delete(f"{api_root}/{path}", headers=headers)

    if response.status_code == requests.codes.ok:
        resp = response.json()

        if resp.get("error") is not None:
            raise Exception(resp["error"])
        else:
            return resp

    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()
