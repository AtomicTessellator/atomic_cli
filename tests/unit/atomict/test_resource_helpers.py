from unittest.mock import patch

from atomict.resource_helpers import (
    build_path,
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def test_build_path_without_params():
    assert build_path("api/project", "proj-1") == "api/project/proj-1/"


def test_build_path_with_params_and_lists():
    path = build_path(
        "api/project",
        params={"search": "alpha", "tag": ["one", "two"], "skip": None},
    )
    assert path == "api/project/?search=alpha&tag=one&tag=two"


@patch("atomict.resource_helpers.get")
def test_retrieve_resource(mock_get):
    mock_get.return_value = {"id": "proj-1"}

    result = retrieve_resource("api/project", "proj-1", include_ht=True)

    assert result == {"id": "proj-1"}
    mock_get.assert_called_once_with("api/project/proj-1/?include_ht=True")


@patch("atomict.resource_helpers.get")
def test_list_resources(mock_get):
    mock_get.return_value = {"results": []}

    result = list_resources("api/tasks", depth=2, status=3)

    assert result == {"results": []}
    mock_get.assert_called_once_with("api/tasks/?depth=2&status=3")


@patch("atomict.resource_helpers.post")
def test_create_resource_uses_json_headers(mock_post):
    mock_post.return_value = {"id": "note-1"}

    result = create_resource("api/project-note", {"title": "Example"})

    assert result == {"id": "note-1"}
    mock_post.assert_called_once_with(
        "api/project-note/",
        {"title": "Example"},
        extra_headers={"Content-Type": "application/json"},
    )


@patch("atomict.resource_helpers.patch")
def test_update_resource_defaults_to_patch(mock_patch):
    mock_patch.return_value = {"id": "task-1", "status": 6}

    result = update_resource("api/tasks", "task-1", {"status": 6})

    assert result == {"id": "task-1", "status": 6}
    mock_patch.assert_called_once_with("api/tasks/task-1/", payload={"status": 6})


@patch("atomict.resource_helpers.put")
def test_update_resource_can_use_put(mock_put):
    mock_put.return_value = {"id": "cluster-1"}

    result = update_resource(
        "api/k8s-cluster", "cluster-1", {"name": "new"}, partial=False
    )

    assert result == {"id": "cluster-1"}
    mock_put.assert_called_once_with(
        "api/k8s-cluster/cluster-1/", payload={"name": "new"}
    )


@patch("atomict.resource_helpers.delete")
def test_delete_resource(mock_delete):
    mock_delete.return_value = {"status": "ok"}

    result = delete_resource("api/object-link", "link-1")

    assert result == {"status": "ok"}
    mock_delete.assert_called_once_with("api/object-link/link-1/")
