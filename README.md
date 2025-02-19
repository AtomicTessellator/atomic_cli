# Atomic Tessellator - CLI package

## Installation
```
pip install atomict

# install utils package with scientific extensions (ase, spglib) 
pip install atomict[utils]
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

```$ at```

![Alt text](img/at.png?raw=true "at")

You can get help for each command with `at <command> --help`. This will print the command-specific help and options.
## Log in and store authentication token

```$ at login```

This will prompt you for your username and password.


## Get a list of available projects

```$ at project get```

![Alt text](img/at_project_get.png?raw=true "at project get")


## Get a list of adsorbates

```$ at adsorbate get```

![Alt text](img/at_adsorbate_get.png?raw=true "at adsorbate get")

## Search objects for substrings

```$ at adsorbate get --search NH3```

![Alt text](img/at_adsorbate_get_search.png?raw=true "at adsorbate get --search")

## Get a list of tasks

```$ at task get```

![Alt text](img/at_task_get.png?raw=true "at task get")

## Get tasks with a specific status

```$ at task get --status completed```

![Alt text](img/at_task_get_completed.png?raw=true "at task get --status completed")

## Configuration

Tab completion is available. Run the hidden command:

```at completion```

This will print out the instructions for enabling tab completion for your shell.
