import logging
import os
import shutil
from typing import Any, Dict, List, Optional

from atomict.infra.distwork.task import TaskStatus, update_task_status
from atomict.io.utils import human_filesize
from atomict.user.files import (
    delete_user_upload,
    download_file,
    get_user_uploads,
    upload_single_file,
)


def display_name(user_upload):
    if "users_name" in user_upload and user_upload["users_name"] not in ["", None]:
        return user_upload["users_name"]
    return user_upload["orig_name"]


def clear_workspace(sim, base_path: str = "./workspace"):

    target_dir = os.path.join(base_path, sim["id"])

    if os.path.exists(target_dir):
        logging.warning(f"Removing existing workspace folder {target_dir}")
        shutil.rmtree(target_dir)


def download_workspace(workspace_files, target_directory: str):
    """Download a workspace to a target directory

    Args:
        workspace_files list of UserUpload objects: List of files to download
        target_directory (str): The directory to download the files to
    """

    os.makedirs(target_directory, exist_ok=True)

    total_bytes = sum([f["user_upload"]["size"] for f in workspace_files])
    finished_bytes = 0
    for sim_file in workspace_files:
        logging.info(
            f"Downloading file {display_name(sim_file['user_upload'])} ({human_filesize(sim_file['user_upload']['size'])})"
        )

        download_file(
            sim_file["user_upload"]["uuid"],
            f"{target_directory}/{sim_file['user_upload']['users_name']}",
        )

        finished_bytes += sim_file["user_upload"]["size"]
        logging.info(
            f"Downloaded {human_filesize(finished_bytes)} of {human_filesize(total_bytes)} ({finished_bytes/total_bytes*100:.1f}%)"
        )


def upload_workspace(
    sim, associate_function, workspace_folder: str, starting_percent: int = 80
):
    """
    Uploads a workspace folder to the Atomic platform, and associates the uploaded files with the given simulation.

    Args:
        sim: The simulation object to associate the uploaded files with
        associate_function: The function to associate the uploaded files with the simulation (e.g. associate_user_upload_with_qe_simulation)
        workspace_folder: The folder to upload
        starting_percent: The starting percentage to use when updating the task status
    """

    simulation_id = sim["id"]

    total_size = sum(
        os.path.getsize(os.path.join(root, file))
        for root, _, files in os.walk(workspace_folder)
        for file in files
    )
    uploaded_size = 0
    last_update_percent = starting_percent

    update_task_status(
        sim["task"]["id"],
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
                if result is None:
                    raise Exception(
                        f"Failed to upload {inner_workspace} - no result returned"
                    )
                elif result.get("status") != "OK":
                    logging.error(f"Failed to upload {inner_workspace}")
                    logging.error(result)
                    raise Exception(f"Failed to upload {inner_workspace} {result}")
                else:
                    logging.info(f"Uploaded {inner_workspace} OK")
            except Exception as e:
                logging.error(f"Failed to upload {inner_workspace}")
                logging.error(e)
                raise

            if result and "UserUpload" in result:
                associate_function(result["UserUpload"]["id"], simulation_id)

            uploaded_size += file_size
            current_percent = 80 + int((uploaded_size / total_size) * 20)

            # Update status if at least 2% has changed
            if current_percent - last_update_percent >= 2:
                update_task_status(sim["task"]["id"], percent=current_percent)
                last_update_percent = current_percent


def delete_user_workspace(workspace_type: str = "all") -> Dict[str, Any]:
    """Delete user workspace files.

    Args:
        workspace_type: Type of workspace to delete ("all", "uploads", "temp")

    Returns:
        Dictionary with deletion results including counts and any errors
    """
    results = {"deleted_count": 0, "failed_count": 0, "errors": []}

    try:
        # Get all user uploads
        uploads_response = get_user_uploads()
        if uploads_response is None:
            uploads = []
        elif isinstance(uploads_response, dict):
            uploads = uploads_response.get("results", [])
        else:
            uploads = uploads_response

        for upload in uploads:
            try:
                upload_uuid = upload.get("uuid") or upload.get("id")
                if upload_uuid:
                    delete_user_upload(upload_uuid)
                    results["deleted_count"] += 1
                    logging.info(
                        f"Deleted upload: {upload.get('orig_name', upload_uuid)}"
                    )
            except Exception as e:
                results["failed_count"] += 1
                upload_name = upload.get(
                    "orig_name", upload.get("uuid", upload.get("id", "unknown"))
                )
                error_msg = f"Failed to delete {upload_name}: {str(e)}"
                results["errors"].append(error_msg)
                logging.error(error_msg)

    except Exception as e:
        error_msg = f"Failed to retrieve workspace files: {str(e)}"
        results["errors"].append(error_msg)
        logging.error(error_msg)

    return results


def get_workspace_summary() -> Dict[str, Any]:
    """Get workspace usage statistics and file counts.

    Returns:
        Dictionary containing workspace summary statistics
    """
    summary = {"total_files": 0, "total_size": 0, "file_types": {}, "recent_files": []}

    try:
        uploads_response = get_user_uploads(limit=1000)  # Get up to 1000 files
        if uploads_response is None:
            uploads = []
        elif isinstance(uploads_response, dict):
            uploads = uploads_response.get("results", [])
        else:
            uploads = uploads_response

        summary["total_files"] = len(uploads)

        for upload in uploads:
            file_size = upload.get("size", 0)
            summary["total_size"] += file_size

            # Track file types
            file_type = upload.get("type", "unknown")
            summary["file_types"][file_type] = (
                summary["file_types"].get(file_type, 0) + 1
            )

        # Get recent files (first 10)
        summary["recent_files"] = uploads[:10]

    except Exception as e:
        logging.error(f"Failed to get workspace summary: {str(e)}")
        summary["error"] = str(e)

    return summary


def clean_workspace(
    older_than_days: int = 30,
    file_types: Optional[List[str]] = None,
    max_files: int = 100,
) -> Dict[str, Any]:
    """Clean workspace of old files by criteria.

    Args:
        older_than_days: Delete files older than this many days
        file_types: List of file types to target (None = all types)
        max_files: Maximum number of files to delete in one operation

    Returns:
        Dictionary with cleaning results
    """
    from datetime import datetime, timedelta

    results = {"deleted_count": 0, "failed_count": 0, "skipped_count": 0, "errors": []}

    cutoff_date = datetime.now() - timedelta(days=older_than_days)

    try:
        uploads_response = get_user_uploads(ordering="-uploaded")
        if uploads_response is None:
            uploads = []
        elif isinstance(uploads_response, dict):
            uploads = uploads_response.get("results", [])
        else:
            uploads = uploads_response

        deleted_count = 0
        for upload in uploads:
            if deleted_count >= max_files:
                break

            # Check file age
            uploaded_str = upload.get("uploaded")
            if uploaded_str:
                try:
                    # Parse the uploaded date (assumes ISO format)
                    uploaded_date = datetime.fromisoformat(
                        uploaded_str.replace("Z", "+00:00")
                    )
                    if uploaded_date.replace(tzinfo=None) > cutoff_date:
                        results["skipped_count"] += 1
                        continue
                except Exception:
                    # Skip if we can't parse the date
                    results["skipped_count"] += 1
                    continue

            # Check file type filter
            if file_types and upload.get("type") not in file_types:
                results["skipped_count"] += 1
                continue

            # Delete the file
            try:
                upload_uuid = upload.get("uuid") or upload.get("id")
                if upload_uuid:
                    delete_user_upload(upload_uuid)
                    results["deleted_count"] += 1
                    deleted_count += 1
                    logging.info(
                        f"Cleaned old file: {upload.get('orig_name', upload_uuid)}"
                    )
            except Exception as e:
                results["failed_count"] += 1
                upload_name = upload.get(
                    "orig_name", upload.get("uuid", upload.get("id", "unknown"))
                )
                error_msg = f"Failed to clean {upload_name}: {str(e)}"
                results["errors"].append(error_msg)
                logging.error(error_msg)

    except Exception as e:
        error_msg = f"Failed to clean workspace: {str(e)}"
        results["errors"].append(error_msg)
        logging.error(error_msg)

    return results
