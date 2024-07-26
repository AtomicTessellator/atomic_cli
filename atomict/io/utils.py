

def human_filesize(size):
    """
    Convert a file size in bytes to a human readable format
    """
    for unit in ["bytes", "KB", "MB", "GB", "TB"]:
        if size < 1024.0:
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} PB"
