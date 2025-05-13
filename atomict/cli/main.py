# cli/main.py
import logging
import os
import sys

import click
from rich.console import Console

from atomict.__version__ import __version__
from atomict.cli.commands import login, user
from atomict.io.msgpack import load_msgpack

# Import command groups
from .commands import adsorbate, catalysis, k8s, project, task, traj, upload
from .commands.exploration import soec, sqs
from .commands.simulation import fhiaims, kpoint, vibes

console = Console()


def setup_logging(verbose: bool):
    """Configure logging based on verbose flag and AT_DEBUG env var"""
    if os.getenv("AT_DEBUG") == "enabled":
        # Most verbose logging when AT_DEBUG is set
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            handlers=[logging.StreamHandler(), logging.FileHandler("atomict.log")],
        )
        # Also enable HTTP library debugging
        logging.getLogger("httpx").setLevel(logging.DEBUG)
        logging.getLogger("httpcore").setLevel(logging.DEBUG)

        # Log some debug info
        logging.debug("Debug mode enabled via AT_DEBUG")
        logging.debug(f'Python path: {os.getenv("PYTHONPATH")}')
        logging.debug(f"Working directory: {os.getcwd()}")
    else:
        # Normal logging based on verbose flag
        level = logging.DEBUG if verbose else logging.ERROR
        logging.basicConfig(
            level=level, format="%(asctime)s - %(levelname)s - %(message)s"
        )

@click.command(name="convert")
@click.argument("input_file")
@click.argument("output_file")
@click.option("--allow-int-keys", is_flag=True, default=False, help="Allow integer keys in msgpack files (sets strict_map_key=False)")
def convert_command(input_file, output_file, allow_int_keys):
    """Convert between atomic structure file formats using ASE

    Supports all formats that ASE can read/write, with special handling for .atm and .atraj files.
    Usage examples:
      tess input.cif output.xyz
      tess input.xyz output.atm
      tess input.traj output.atraj
      
    Options:
      --allow-int-keys  Allow integer keys in msgpack files (default: False)
    """
    try:
        import os.path
        from ase.io import read, write
        from ase.io.formats import UnknownFileTypeError
        from atomict.io.msgpack import save_msgpack, save_msgpack_trajectory, load_msgpack, load_msgpack_trajectory
    except ImportError:
        console.print("[red]Error: ASE (Atomic Simulation Environment) is required for file conversion.[/red]")
        console.print("[yellow]Install it with: pip install ase[/yellow]")
        return

    RW_FORMATS = [
        'abinit-in', 'aims', 'bundletrajectory', 'castep-cell', 'castep-geom', 
        'castep-md', 'cfg', 'cif', 'crystal', 'cube', 'db', 'dftb', 'dlp4', 
        'dmol-arc', 'dmol-car', 'dmol-incoor', 'eon', 'espresso-in', 'extxyz', 
        'gaussian-in', 'gen', 'gpumd', 'gromacs', 'gromos', 'json', 'jsv', 
        'lammps-data', 'magres', 'mustem', 'mysql', 'netcdftrajectory', 'nwchem-in', 
        'onetep-in', 'postgresql', 'prismatic', 'proteindatabank', 'res', 
        'rmc6f', 'struct', 'sys', 'traj', 'turbomole', 'v-sim', 'vasp', 
        'vasp-xdatcar', 'xsd', 'xsf', 'xtd', 'xyz'
    ]

    try:
        input_ext = os.path.splitext(input_file)[1].lower()[1:]
        output_ext = os.path.splitext(output_file)[1].lower()[1:]

        if not os.path.exists(input_file):
            console.print(f"[red]Error: Input file '{input_file}' not found.[/red]")
            return

        # Consider .atm and .atraj as special msgpack formats
        msgpack_formats = ["atm"]
        traj_msgpack_formats = ["atraj"]  # Separate trajectory msgpack format
        
        if input_ext not in RW_FORMATS and input_ext not in msgpack_formats and input_ext not in traj_msgpack_formats:
            console.print(f"[red]Error: Format '{input_ext}' is not supported for reading.[/red]")
            console.print("[yellow]Supported read/write formats include:[/yellow]")
            # Display formats in multiple columns for better readability
            for i in range(0, len(RW_FORMATS), 5):
                console.print("[yellow]  " + ", ".join(RW_FORMATS[i:i+5]) + "[/yellow]")
            console.print("[yellow]Special formats: atm (msgpack), atraj (msgpack trajectory)[/yellow]")
            return
            
        if output_ext not in RW_FORMATS and output_ext not in msgpack_formats and output_ext not in traj_msgpack_formats:
            console.print(f"[red]Error: Format '{output_ext}' is not supported for writing.[/red]")
            console.print("[yellow]Supported read/write formats include:[/yellow]")
            # Display formats in multiple columns for better readability
            for i in range(0, len(RW_FORMATS), 5):
                console.print("[yellow]  " + ", ".join(RW_FORMATS[i:i+5]) + "[/yellow]")
            console.print("[yellow]Special formats: atm (msgpack), atraj (msgpack trajectory)[/yellow]")
            return

        # The strict_map_key parameter is the opposite of allow_int_keys
        strict_map_key = not allow_int_keys
        
        try:
            if input_ext in msgpack_formats:
                atoms = load_msgpack(input_file, strict_map_key=strict_map_key)
            elif input_ext in traj_msgpack_formats:
                atoms, _ = load_msgpack_trajectory(input_file, strict_map_key=strict_map_key)
            else:
                atoms = read(input_file)
        except UnknownFileTypeError:
            console.print(f"[red]Error: Unknown file type for input file '{input_file}'[/red]")
            console.print(f"[yellow]The file extension '{input_ext}' is not recognized.[/yellow]")
            console.print("[yellow]Make sure the file has the correct extension for its format.[/yellow]")
            return
        except Exception as e:
            console.print(f"[red]Error reading input file '{input_file}': {str(e)}. If loading .atraj, try adding the --allow-int-keys flag.[/red]")
            console.print(f"[yellow]Make sure '{input_ext}' is a valid format and the file is not corrupted.[/yellow]")
            if "using a non-string key" in str(e) and strict_map_key:
                console.print("[yellow]Try using --allow-int-keys to allow integer keys in the msgpack file.[/yellow]")
            return
        
        try:
            if output_ext in msgpack_formats:
                save_msgpack(atoms, output_file)
                console.print(f"[green]Successfully converted {input_file} to {output_file} (MSGPACK format)[/green]")
            elif output_ext in traj_msgpack_formats:
                save_msgpack_trajectory(atoms, output_file)
                console.print(f"[green]Successfully converted {input_file} to {output_file} (MSGPACK trajectory format)[/green]")
            else:
                write(output_file, atoms)
                console.print(f"[green]Successfully converted {input_file} to {output_file}[/green]")

        except UnknownFileTypeError:
            console.print(f"[red]Error: Unknown file type for output file '{output_file}'[/red]")
            console.print(f"[yellow]The file extension '{output_ext}' is not recognized.[/yellow]")
            console.print("[yellow]Make sure the file has the correct extension for its format.[/yellow]")
            return
        except Exception as e:
            console.print(f"[red]Error writing output file '{output_file}': {str(e)}[/red]")
            console.print(f"[yellow]Make sure '{output_ext}' is a valid format and you have write permissions.[/yellow]")
            return
            
    except Exception as e:
        logging.debug(f"Conversion failed with error: {str(e)}", exc_info=True)
        console.print(f"[red]Error during conversion: {str(e)}[/red]")
        console.print("[yellow]Try running with --verbose for more detailed error information.[/yellow]")


# Create main CLI group
@click.group(invoke_without_command=True)
@click.option(
    "-v", "--verbose", is_flag=True, default=False, help="Enable verbose output"
)
@click.version_option(prog_name="tess", version=__version__)
@click.pass_context
def cli(ctx, verbose: bool):
    """Atomic Tessellator CLI - Manage simulations and computational resources
    
    Dynamic commands (takes positional arguments only):\n
      tess [file_a] [file_b] Optional[--allow-int-keys] Filetype conversion from a -> b
    """
    setup_logging(verbose)
    
    # Show help if no subcommand was provided
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


# Add a completion command
@cli.command(hidden=True)
@click.argument("shell", type=click.Choice(["bash", "zsh", "fish"]), required=False)
def completion(shell):
    """Generate shell completion script"""
    if shell is None:
        shell = os.environ.get("SHELL", "")
        shell = shell.split("/")[-1]
        if shell not in ["bash", "zsh", "fish"]:
            shell = "bash"  # default to bash if shell not detected

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
    click.echo(completion_script.strip())
    click.echo(
        "# Don't forget to source your rc file! `source ~/.bashrc` or `source ~/.zshrc` ..."
    )


cli.add_command(completion)

# TODO: rename these to group
cli.add_command(task.task)
cli.add_command(upload.upload)
cli.add_command(project.project)
cli.add_command(k8s.k8s)
cli.add_command(adsorbate.adsorbate)

# raise commands to top-level
cli.add_command(fhiaims.fhiaims_group)
cli.add_command(kpoint.kpoint_group)
cli.add_command(catalysis.catalysis_group)  # WIP
cli.add_command(sqs.sqs_group)
cli.add_command(soec.soecexploration_group)
cli.add_command(traj.traj)
cli.add_command(user.user_group)
cli.add_command(login._login)
cli.add_command(vibes.vibes_group)
# from .commands.exploration import exploration_group
# cli.add_command(exploration.exploration)  # move this
# TBD: decide on how to group or put all commands at top level
# standardize this later
# cli.add_command(exploration_group)

# we could do `at [exploration/simulation/project/user/etc] [get/create/delete] [id]`
# OR
# `at [get/create/delete] [exploration/simulation/project/user/etc] [id]`
# OR
# some other grouping that reflects a user's workflow


def main():
    try:
        # Check if we have a simple file conversion case
        if len(sys.argv) >= 3:
            arg1 = sys.argv[1]
            arg2 = sys.argv[2]
            if '.' in arg1 and '.' in arg2:
                # Default is False (don't allow int keys)
                allow_int_keys = False
                
                # Check if --allow-int-keys flag is present in the arguments
                if "--allow-int-keys" in sys.argv:
                    allow_int_keys = True
                    
                convert_command.callback(arg1, arg2, allow_int_keys)
                return

        cli()
    except Exception as exc:
        Console().print(f"[red]Error: {str(exc)}. Exiting...[/red]")
        sys.exit(1)


if __name__ == "__main__":
    main()
