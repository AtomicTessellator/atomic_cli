# Atomic Tessellator - CLI package

## Installation
```
pip install atomict
```

## Installation for devs
```
pip install -e ".[dev]"
```

Enable verbose logging for debugging:

```
export AT_DEBUG=enabled
```

## CLI Usage
```
Usage: at [OPTIONS] COMMAND [ARGS]...

  Atomic Tessellator CLI - Manage simulations and computational resources

Options:
  -v, --verbose  Enable verbose output
  --version      Show the version and exit.
  --help         Show this message and exit.

Commands:
  adsorbate  Manage adsorbates
  catalysis  Manage catalysis explorations
  ea         Manage EA / SOEC explorations and related resources
  fhiaims    Manage FHI-aims simulations
  k8s        Manage Kubernetes jobs and clusters
  kpoint     Manage K-point simulations
  project    Manage projects and their related resources
  sqs        Manage Special Quasirandom Structure (SQS) explorations
  task       Manage tasks and their status
  upload     Manage file uploads
  user       Manage users and user uploads

```

# Configuration

Tab completion is available. Run the hidden command:

```at completion```

This will print out the instructions for enabling tab completion for your shell.

