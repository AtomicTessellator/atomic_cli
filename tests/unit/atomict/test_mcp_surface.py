from atomict.mcp_surface import (
    MCP_RESOURCE_REGISTRY,
    MCP_SOURCE_PARITY_SURFACE,
    resolve_operation,
)


def test_registry_contains_expected_resources():
    for resource in [
        "project",
        "project-note",
        "user-upload",
        "task",
        "fhiaims-simulation",
        "kpoint-exploration",
        "sqs-exploration",
        "ea-exploration",
        "mlrelax",
    ]:
        assert resource in MCP_RESOURCE_REGISTRY


def test_expected_core_operations_are_registered():
    assert "get" in MCP_RESOURCE_REGISTRY["project"]
    assert "list" in MCP_RESOURCE_REGISTRY["project"]
    assert "create" in MCP_RESOURCE_REGISTRY["project-note"]
    assert "download_content" in MCP_RESOURCE_REGISTRY["user-upload"]
    assert "cancel" in MCP_RESOURCE_REGISTRY["task"]
    assert "create" in MCP_RESOURCE_REGISTRY["sqs-target-concentration"]


def test_all_surface_operations_have_notes():
    for operation in MCP_SOURCE_PARITY_SURFACE:
        assert operation.notes, (
            f"Expected notes for {operation.resource}.{operation.operation}"
        )


def test_selected_operations_have_expected_notes():
    project_notes = MCP_RESOURCE_REGISTRY["project"]["get"].notes
    assert project_notes.startswith("Get project details.")
    assert "project_id (str): The project identifier." in project_notes
    assert "Returns:\n    dict: The API response payload." in project_notes

    user_upload_notes = MCP_RESOURCE_REGISTRY["user-upload"]["download_content"].notes
    assert user_upload_notes.startswith("Download the content of a user upload.")
    assert "upload_id (str): The user upload identifier." in user_upload_notes

    workspace_notes = MCP_RESOURCE_REGISTRY["workspace"]["lookup_simulation"].notes
    assert workspace_notes.startswith("Look up the workspace for a simulation.")
    assert "simulation_uuid (str): The simulation UUID." in workspace_notes

    task_notes = MCP_RESOURCE_REGISTRY["task"]["cancel"].notes
    assert task_notes.startswith("Cancel a running task.")
    assert "task_id (str): The task identifier." in task_notes

    sqs_notes = MCP_RESOURCE_REGISTRY["sqs-target-concentration"]["list"].notes
    assert sqs_notes.startswith("List target concentrations for an SQS exploration.")
    assert "**params (Any): Additional query parameters to include in the request." in sqs_notes


def test_selected_registered_operations_resolve_to_callables():
    selected = [
        MCP_RESOURCE_REGISTRY["project"]["get"],
        MCP_RESOURCE_REGISTRY["project-note"]["create"],
        MCP_RESOURCE_REGISTRY["user-upload"]["download_content"],
        MCP_RESOURCE_REGISTRY["task"]["cancel"],
        MCP_RESOURCE_REGISTRY["fhiaims-simulation"]["list"],
        MCP_RESOURCE_REGISTRY["sqs-exploration"]["create"],
        MCP_RESOURCE_REGISTRY["ea-exploration"]["update"],
        MCP_RESOURCE_REGISTRY["mlrelax"]["get"],
    ]

    for operation in selected:
        resolved = resolve_operation(operation)
        assert callable(resolved), (
            f"Failed to resolve {operation.module}.{operation.function}"
        )
