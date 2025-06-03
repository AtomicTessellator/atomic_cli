import json
import os
import logging

import requests
from requests.exceptions import (
    HTTPError,
    ConnectionError,
    Timeout,
)
from tenacity import retry, stop_after_attempt, wait_exponential, after_log, retry_if_exception_type, retry_if_exception, before_sleep_log

from atomict.exceptions import APIValidationError, PermissionDenied
from atomict.utils.tenacity import before_log


logger = logging.getLogger(__name__)


def is_http_5xx_error(exception):
    """These are the HTTP errors that we should attempt to retry on."""
    return isinstance(exception, HTTPError) and 500 <= exception.response.status_code < 600


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=1, max=30),
    retry=(retry_if_exception_type((ConnectionError, Timeout)) | retry_if_exception(is_http_5xx_error)),
    before=before_log(logger, logging.INFO),
    after=after_log(logger, logging.INFO),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def get(path: str):
    api_root = os.environ.get("AT_SERVER", "https://api.atomictessellator.com")
    headers = {"Accept": "application/json", "Content-Type": "application/json"}

    if os.environ.get("AT_TOKEN"):
        headers["Authorization"] = f"Token {os.environ.get('AT_TOKEN')}"

    response = requests.get(f"{api_root}/{path}", headers=headers)

    content_type = response.headers.get("Content-Type")

    if response.status_code == requests.codes.ok and content_type == "application/json":
        resp = response.json()

        # This is truely hideous, an error can take many forms, not just permission denied
        # This really needs to be cleaned up
        # This hack about the paths is because the task object is special case where it has a legitimate
        # field called "error" - 
        # We need to standardize a response type from the server, likely the django rest framework response
        if "error" in resp and resp["error"] is not None and "api/tasks/" not in path:
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


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=1, max=30),
    retry=(retry_if_exception_type((ConnectionError, Timeout)) | retry_if_exception(is_http_5xx_error)),
    before=before_log(logger, logging.INFO),
    after=after_log(logger, logging.INFO),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
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

    api_root = os.environ.get("AT_SERVER", "https://api.atomictessellator.com")

    if files is not None:
        response = requests.post(
            f"{api_root}/{path}", data=payload, headers=headers, files=files
        )
    else:
        response = requests.post(f"{api_root}/{path}", data=payload, headers=headers)

    if response.status_code in [requests.codes.ok, requests.codes.created]:
        resp = response.json()

        if "error" in resp and "api/tasks/" not in path:
            raise PermissionDenied(resp["error"])
        else:
            return resp

    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=1, max=30),
    retry=(retry_if_exception_type((ConnectionError, Timeout)) | retry_if_exception(is_http_5xx_error)),
    before=before_log(logger, logging.INFO),
    after=after_log(logger, logging.INFO),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def patch(path: str, payload: dict):
    payload_enc = json.dumps(payload)
    headers = {"Content-Type": "application/json"}

    if os.environ.get("AT_TOKEN"):
        headers["Authorization"] = f"Token {os.environ.get('AT_TOKEN')}"

    api_root = os.environ.get("AT_SERVER", "https://api.atomictessellator.com")
    response = requests.patch(f"{api_root}/{path}", data=payload_enc, headers=headers)

    if response.status_code == requests.codes.ok:
        resp = response.json()

        if resp.get("error") is not None and "api/tasks/" not in path:
            raise Exception(resp["error"])
        else:
            return resp

    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=1, max=30),
    retry=(retry_if_exception_type((ConnectionError, Timeout)) | retry_if_exception(is_http_5xx_error)),
    before=before_log(logger, logging.INFO),
    after=after_log(logger, logging.INFO),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def delete(path: str):
    headers = {}

    if os.environ.get("AT_TOKEN"):
        headers["Authorization"] = f"Token {os.environ.get('AT_TOKEN')}"

    api_root = os.environ.get("AT_SERVER", "https://api.atomictessellator.com")
    response = requests.delete(f"{api_root}/{path}", headers=headers)

    if response.status_code == requests.codes.ok:
        resp = response.json()

        if resp.get("error") is not None and "api/tasks/" not in path:
            raise Exception(resp["error"])
        else:
            return resp

    elif response.status_code == requests.codes.bad_request:
        raise APIValidationError(response.json())
    elif response.status_code == requests.codes.forbidden:
        raise PermissionDenied(response.json())
    else:
        response.raise_for_status()
