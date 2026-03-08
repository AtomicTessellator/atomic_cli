from __future__ import annotations

from typing import Any, Mapping
from urllib.parse import urlencode

from atomict.api import delete, get, patch, post, put


def _clean_params(params: Mapping[str, Any] | None = None) -> dict[str, Any]:
    if not params:
        return {}
    return {key: value for key, value in params.items() if value is not None}


def build_path(
    resource_path: str,
    resource_id: str | None = None,
    params: Mapping[str, Any] | None = None,
) -> str:
    base_path = resource_path.strip("/")
    if resource_id is not None:
        base_path = f"{base_path}/{resource_id}"
    base_path = f"{base_path}/"

    query_params = _clean_params(params)
    if not query_params:
        return base_path

    return f"{base_path}?{urlencode(query_params, doseq=True)}"


def retrieve_resource(
    resource_path: str, resource_id: str, **params: Any
) -> dict[str, Any]:
    return get(build_path(resource_path, resource_id, params=params))


def list_resources(resource_path: str, **params: Any) -> dict[str, Any]:
    return get(build_path(resource_path, params=params))


def create_resource(
    resource_path: str,
    payload: Mapping[str, Any],
    *,
    json_payload: bool = True,
) -> dict[str, Any]:
    extra_headers = {"Content-Type": "application/json"} if json_payload else {}
    return post(build_path(resource_path), dict(payload), extra_headers=extra_headers)


def update_resource(
    resource_path: str,
    resource_id: str,
    fields: Mapping[str, Any],
    *,
    partial: bool = True,
) -> dict[str, Any]:
    path = build_path(resource_path, resource_id)
    payload = dict(fields)
    if partial:
        return patch(path, payload=payload)
    return put(path, payload=payload)


def delete_resource(resource_path: str, resource_id: str) -> dict[str, Any]:
    return delete(build_path(resource_path, resource_id))
