# Atomic Tessellator - CLI package

## Installation
```
pip install atomict

# install utils package with scientific extensions (ase, spglib) 
pip install atomict[utils]

# ", with a git url
pip install git+https://github.com/AtomicTessellator/atomic_cli#egg=atomict[utils]

# Troubleshooting
If pip install is not working, try:
pip install --upgrade setuptools pip wheel
```

## Installation for devs
```
pip install -e ".[dev]"
```

Enable verbose logging for debugging:

```
export AT_DEBUG=enabled
```

# CLI Usage

## Get a list of available commands

```$ tess```

![Alt text](img/at.png?raw=true "at")

You can get help for each command with `tess <command> --help`. This will print the command-specific help and options.
## Log in and store authentication token

```$ tess login```

This will prompt you for your username and password.


## Get a list of available projects

```$ tess project get```

![Alt text](img/at_project_get.png?raw=true "tess project get")


## Get a list of adsorbates

```$ tess adsorbate get```

![Alt text](img/at_adsorbate_get.png?raw=true "tess adsorbate get")

## Search objects for substrings

```$ tess adsorbate get --search NH3```

![Alt text](img/at_adsorbate_get_search.png?raw=true "tess adsorbate get --search")

## Get a list of tasks

```$ tess task get```

![Alt text](img/at_task_get.png?raw=true "tess task get")

## Get tasks with a specific status

```$ tess task get --status completed```

![Alt text](img/at_task_get_completed.png?raw=true "tess task get --status completed")

## Configuration

Tab completion is available. Run the hidden command:

```tess completion```

This will print out the instructions for enabling tab completion for your shell.
