# cli/main.py
import logging
import os
import sys

import click
from rich.console import Console

from atomict.__version__ import __version__
from atomict.cli.commands import login, token, user, workspace
from .commands import k8s, project, task, traj, upload, convert
from .commands.exploration import soec, sqs
from .commands.simulation import fhiaims, kpoint


console = Console()


def setup_logging(verbose: bool):
    """Configure logging based on verbose flag and AT_DEBUG env var"""
    if os.getenv("AT_DEBUG") == "enabled":
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[logging.StreamHandler(), logging.FileHandler("atomict.log")],
        )
        logging.getLogger("httpx").setLevel(logging.DEBUG)
        logging.getLogger("httpcore").setLevel(logging.DEBUG)
        logging.debug("Debug mode enabled via AT_DEBUG")
        logging.debug(f'Python path: {os.getenv("PYTHONPATH")}')
        logging.debug(f"Working directory: {os.getcwd()}")
    else:
        level = logging.DEBUG if verbose else logging.ERROR
        logging.basicConfig(
            level=level, format="%(asctime)s - %(levelname)s - %(message)s"
        )


@click.group(invoke_without_command=True)
@click.option(
    "-v", "--verbose", is_flag=True, default=False, help="Enable verbose output"
)
@click.version_option(prog_name="tess", version=__version__)
@click.pass_context
def cli(ctx, verbose: bool):
    """Atomic Tessellator CLI - Manage simulations and computational resources
    """
    setup_logging(verbose)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command(hidden=True)
@click.argument("shell", type=click.Choice(["bash", "zsh", "fish"]), required=False)
def completion(shell):
    """Generate shell completion script"""
    if shell is None:
        shell = os.environ.get("SHELL", "")
        shell = shell.split("/")[-1]
        if shell not in ["bash", "zsh", "fish"]:
            shell = "bash"

    completion_script = None
    if shell == "bash":
        completion_script = """
            # Add to ~/.bashrc:
if tess --version > /dev/null 2>&1; then
    eval "$(_TESS_COMPLETE=bash_source tess)"
fi
            """
    elif shell == "zsh":
        completion_script = """
            # Add to ~/.zshrc:
if tess --version > /dev/null 2>&1; then
    eval "$(_TESS_COMPLETE=zsh_source tess)"
fi
            """
    elif shell == "fish":
        completion_script = """
            # Add to ~/.config/fish/config.fish:
if type -q tess
    eval (env _TESS_COMPLETE=fish_source tess)
end
"""
    click.echo(f"# Shell completion for {shell}")
    click.echo((completion_script or "").strip())
    click.echo(
        "# Don't forget to source your rc file! `source ~/.bashrc` or `source ~/.zshrc` ..."
    )


cli.add_command(completion)
cli.add_command(convert.convert)
cli.add_command(task.task_group)
cli.add_command(upload.upload_group)
cli.add_command(project.project_group)
cli.add_command(k8s.k8s_group)
cli.add_command(fhiaims.fhiaims_group)
cli.add_command(kpoint.kpoint_group)
cli.add_command(sqs.sqs_group)
cli.add_command(soec.soecexploration_group)
cli.add_command(traj.traj)
cli.add_command(user.user_group)
cli.add_command(workspace.workspace_group)
cli.add_command(login._login)
cli.add_command(token._token)


def main():
    try:
        cli()
    except Exception as exc:
        Console().print(f"[red]Error: {str(exc)}. Exiting...[/red]")
        sys.exit(1)


if __name__ == "__main__":
    main()
