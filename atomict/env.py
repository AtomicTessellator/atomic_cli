import json
import os
import platform

DEFAULT_ATOMICT_API_ROOT = "https://api.atomictessellator.com"


def store(key, value):
    cfg_path = get_config_path("config.json")

    if os.path.exists(cfg_path):
        with open(get_config_path("config.json"), 'r') as f:
            config = json.load(f)
            config[key] = value
    else:
        config = {key: value}

    with open(get_config_path("config.json"), 'w') as f:
        json.dump(config, f)


def get(key):
    with open(get_config_path("config.json")) as f:
        config = json.load(f)
        return config.get(key)


def delete(key):
    with open(get_config_path("config.json")) as f:
        config = json.load(f)
        del config[key]

    with open(get_config_path("config.json"), "w") as f:
        json.dump(config, f)


def clear():
    path = get_config_path("config.json")

    if os.path.exists(path):
        os.remove(path)


def get_config_path(filename):
    system = platform.system()
    if system == "Darwin":  # macOS
        base_path = os.path.join(
            os.path.expanduser("~/Library/Application Support"), "atomict"
        )
    elif system == "Windows":
        base_path = os.path.join(os.environ["APPDATA"], "atomict")
    else:  # Linux and other Unix-like systems
        config_path = os.path.join(os.path.expanduser("~"), ".config")
        base_path = os.path.join(config_path, "atomict")

    # Create directory if it doesn't exist
    if not os.path.exists(base_path):
        os.makedirs(base_path)

    return os.path.join(base_path, filename)


def check_api_root():
    api_root = os.environ.get("ATOMICT_API_ROOT")

    if api_root is not None:
        # Make sure it's a valid URL without a trailing slash
        if not api_root.startswith("http"):
            raise ValueError(
                "ATOMICT_API_ROOT must be a valid URL, currently: {api_root}"
            )
        if api_root.endswith("/"):
            raise ValueError(
                "ATOMICT_API_ROOT must not end with a trailing slash, currently: {api_root}"
            )


def check_environment():
    check_api_root()
