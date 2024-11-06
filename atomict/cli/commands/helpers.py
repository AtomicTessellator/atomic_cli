# cli/commands/simulation/helpers.py
from datetime import datetime
from typing import Any, Dict, Optional

from rich.table import Table


# TODO: move to common?
def get_status_string(status_code: Optional[int]) -> str:
    """Convert status code to human readable string"""
    status_map = {
        0: "Draft",
        1: "Ready",
        2: "Running",
        3: "Completed",
        4: "Error",
        5: "Paused",
        6: "User Aborted",
    }
    if status_code is None:
        return "N/A"
    return status_map.get(status_code, "N/A")


def format_datetime(dt_str: str) -> str:
    """Format datetime string for display"""
    try:
        dt = datetime.fromisoformat(dt_str.replace("Z", "+00:00"))
        return dt.strftime("%Y-%m-%d %H:%M:%S")
    except (ValueError, AttributeError):
        return "N/A"


def create_detail_table(
    title: str, data: Dict[str, Any], fields: Dict[str, str]
) -> Table:
    """Create a detail table with specified fields"""
    table = Table(show_header=True, title=title)
    table.add_column("Property")
    table.add_column("Value")

    for label, field in fields.items():
        value = data.get(field, "N/A")
        table.add_row(label, str(value))

    return table


def format_task_status(task: Dict[str, Any]) -> Table:
    """Format task information into a table"""
    if not task:
        return None

    task_info = Table(show_header=True, title="Task Status")
    task_info.add_column("Property")
    task_info.add_column("Value")

    fields = {
        "Status": "status",
        "Progress": lambda t: f"{t.get('progress', 0)}%",
        "Running On": "running_on",
        "Created": lambda t: format_datetime(t.get("created_at", "")),
        "Updated": lambda t: format_datetime(t.get("updated_at", "")),
        "Error": "error",
    }

    for label, field in fields.items():
        if callable(field):
            value = field(task)
        else:
            value = task.get(field, "N/A")
        task_info.add_row(label, str(value))

    return task_info
