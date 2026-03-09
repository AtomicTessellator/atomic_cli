import logging
import os
import shutil
from collections.abc import Callable
from typing import Any

from atomict.infra.distwork.task import TaskStatus, update_task_status
from atomict.io.utils import human_filesize
from atomict.user.files import download_file, upload_single_file


def display_name(user_upload):
    if "users_name" in user_upload and user_upload["users_name"] not in ["", None]:
        return user_upload["users_name"]
    return user_upload["orig_name"]


def clear_workspace(sim: dict[str, Any], base_path: str = "./workspace") -> None:

    """Clear a local workspace directory.

    Args:
        sim (dict[str, Any]): The simulation payload that owns the workspace.
        base_path (str): The base directory containing workspace folders.

    Returns:
        None: This helper removes local files in place and does not return a value.
    """
    target_dir = os.path.join(base_path, sim["id"])

    if os.path.exists(target_dir):
        logging.warning(f"Removing existing workspace folder {target_dir}")
        shutil.rmtree(target_dir)


def download_workspace(
    workspace_files: list[dict[str, Any]], target_directory: str
) -> None:
    """Download workspace files to a local directory.

    Args:
        workspace_files (list[dict[str, Any]]): The workspace file payloads to
            download.
        target_directory (str): The local directory where the files should be
            written.

    Returns:
        None: This helper downloads files to disk and does not return a value.
    """

    os.makedirs(target_directory, exist_ok=True)

    total_bytes = sum([f["user_upload"]["size"] for f in workspace_files])
    finished_bytes = 0
    for sim_file in workspace_files:
        logging.info(
            f"Downloading file {display_name(sim_file['user_upload'])} ({human_filesize(sim_file['user_upload']['size'])})"
        )

        download_file(
            sim_file["user_upload"]["id"],
            f"{target_directory}/{sim_file['user_upload']['users_name']}",
        )

        finished_bytes += sim_file["user_upload"]["size"]
        logging.info(
            f"Downloaded {human_filesize(finished_bytes)} of {human_filesize(total_bytes)} ({finished_bytes/total_bytes*100:.1f}%)"
        )


def upload_workspace(
    sim: dict[str, Any],
    associate_function: Callable[[str, str], object],
    workspace_folder: str,
    starting_percent: int = 80,
) -> None:
    """Upload a local workspace directory.

    Args:
        sim (dict[str, Any]): The simulation payload that owns the workspace.
        associate_function (Callable[[str, str], object]): The callback used to
            associate an uploaded file with the simulation.
        workspace_folder (str): The local workspace directory to upload.
        starting_percent (int): The initial task progress percentage to report.

    Returns:
        None: This helper uploads files and updates task progress in place.
    """

    simulation_id = sim["id"]

    total_size = sum(
        os.path.getsize(os.path.join(root, file))
        for root, _, files in os.walk(workspace_folder)
        for file in files
    )
    uploaded_size = 0
    last_update_percent = starting_percent

    update_task_status(sim["task"]["id"],
        percent=last_update_percent,
        progress_indeterminate=False,
    )

    for root, _, files in os.walk(workspace_folder):
        for file in files:
            file_path = os.path.join(root, file)
            inner_workspace = file_path.replace(workspace_folder, "")
            file_size = os.path.getsize(file_path)

            try:
                result = upload_single_file(file_path, inner_workspace)
                if result["status"] != "OK":
                    logging.error(f"Failed to upload {inner_workspace}")
                    logging.error(result)
                    raise Exception(f"Failed to upload {inner_workspace} {result}")
                else:
                    logging.info(f"Uploaded {inner_workspace} OK")
            except Exception as e:
                logging.error(f"Failed to upload {inner_workspace}")
                logging.error(e)
                raise

            associate_function(result["UserUpload"]["id"], simulation_id)

            uploaded_size += file_size
            current_percent = 80 + int((uploaded_size / total_size) * 20)

            # Update status if at least 2% has changed
            if current_percent - last_update_percent >= 2:
                update_task_status(sim["task"]["id"], percent=current_percent)
                last_update_percent = current_percent
