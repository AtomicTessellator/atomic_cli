# Atomic Tessellator - CLI package

## Installation
```
pip install atomict
```

## Development installation
```
pip install -e ".[dev]"
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
Tab completion is available for all commands. To generate the script to enable this, run the hidden command:
```at completion```
This will print out the script with instructions.

## Roadmap
- [x] CLI
- [x] Download Public Datasets
- [ ] Local dataset caching
- [x] Authentication
- [ ] Reality Server API endpoints - Simulations

