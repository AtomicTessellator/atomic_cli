# cli/main.py
import click
from rich.console import Console
import logging
import os

# Import command groups
from .commands import (
    task,
    upload,
    project,
    k8s,
    adsorbate,
)

# from .commands.simulation import simulation
from .commands.simulation import fhiaims, kpoint
# from .commands.exploration import exploration_group
from .commands.exploration import (
    sqs,
    soec,
)


console = Console()


def setup_logging(verbose: bool):
    """Configure logging based on verbose flag and AT_DEBUG env var"""
    if os.getenv('AT_DEBUG'):
        # Most verbose logging when AT_DEBUG is set
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler('atomict.log')
            ]
        )
        # Also enable HTTP library debugging
        logging.getLogger('httpx').setLevel(logging.DEBUG)
        logging.getLogger('httpcore').setLevel(logging.DEBUG)
  
        # Log some debug info
        logging.debug('Debug mode enabled via AT_DEBUG')
        logging.debug(f'Python path: {os.getenv("PYTHONPATH")}')
        logging.debug(f'Working directory: {os.getcwd()}')
    else:
        # Normal logging based on verbose flag
        level = logging.DEBUG if verbose else logging.ERROR
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )


@click.group()
@click.option('-v', '--verbose', is_flag=True, default=False, help='Enable verbose output')
@click.version_option(version='0.1.0')
def cli(verbose: bool):
    """Atomic Tessellator CLI - Manage simulations and computational resources"""
    setup_logging(verbose)


# Add a completion command
@cli.command(hidden=True)
@click.argument('shell', type=click.Choice(['bash', 'zsh', 'fish']), required=False)
def completion(shell):
    """Generate shell completion script"""
    if shell is None:
        shell = os.environ.get('SHELL', '')
        shell = shell.split('/')[-1]
        if shell not in ['bash', 'zsh', 'fish']:
            shell = 'bash'  # default to bash if shell not detected

    completion_script = None
    if shell == 'bash':
        completion_script = '''
            # Add to ~/.bashrc:
eval "$(_AT_COMPLETE=bash_source at)"
            '''
    elif shell == 'zsh':
        completion_script = '''
            # Add to ~/.zshrc:
eval "$(_AT_COMPLETE=zsh_source at)"
            '''
    elif shell == 'fish':
        completion_script = '''
            # Add to ~/.config/fish/config.fish:
eval "$(_AT_COMPLETE=fish_source at)"
'''
    click.echo(f"# Shell completion for {shell}")
    click.echo(completion_script.strip())
    click.echo("# Don't forget to source your rc file! `source ~/.bashrc` or `source ~/.zshrc` ...")


cli.add_command(completion)
# cli.add_command(simulation)
cli.add_command(task.task)
cli.add_command(upload.upload)
cli.add_command(project.project)
cli.add_command(k8s.k8s)
cli.add_command(adsorbate.adsorbate)

# raise commands to top-level
cli.add_command(fhiaims.fhiaims_group)
cli.add_command(kpoint.kpoint_group)
cli.add_command(sqs.sqs_group)
cli.add_command(soec.soecexploration_group)

# cli.add_command(exploration.exploration)  # move this
# standardize this later
# cli.add_command(exploration_group)


def main():
    cli()


if __name__ == '__main__':
    main()
