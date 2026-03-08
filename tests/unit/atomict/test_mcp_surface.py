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
        "qe-simulation",
    ]:
        assert resource in MCP_RESOURCE_REGISTRY


def test_expected_core_operations_are_registered():
    assert "get" in MCP_RESOURCE_REGISTRY["project"]
    assert "list" in MCP_RESOURCE_REGISTRY["project"]
    assert "create" in MCP_RESOURCE_REGISTRY["project-note"]
    assert "download_content" in MCP_RESOURCE_REGISTRY["user-upload"]
    assert "cancel" in MCP_RESOURCE_REGISTRY["task"]
    assert "create" in MCP_RESOURCE_REGISTRY["sqs-target-concentration"]


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
