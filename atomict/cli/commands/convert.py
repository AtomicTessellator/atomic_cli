import os
import logging
import click
from rich.console import Console
from ase.io import read, write
from ase.io.formats import UnknownFileTypeError
from atomict.io.formats.atraj import read_atraj, write_atraj
from atomict.io.formats.tess import read_tess, write_tess

console = Console()

@click.command()
@click.argument("input_file", required=True)
@click.argument("output_file", required=True)
def convert(input_file: str, output_file: str):
    """Convert atomic structure files between supported formats (e.g. .xyz, .cif, .traj).""" 
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

        traj_msgpack_formats = ["atraj", "tess"]
        
        def _is_supported(ext: str) -> bool:
            return (
                ext in RW_FORMATS
                or ext in traj_msgpack_formats
            )
        
        def _print_supported_error(ext: str, action: str) -> None:
            console.print(f"[red]Error: Format '{ext}' is not supported for {action}.[/red]")
            console.print("[yellow]Supported read/write formats include:[/yellow]")
            for i in range(0, len(RW_FORMATS), 5):
                console.print("[yellow]  " + ", ".join(RW_FORMATS[i:i+5]) + "[/yellow]")
            console.print("[yellow]Special formats: atm (msgpack), atraj (msgpack trajectory), tess (msgpack trajectory)[/yellow]")
        
        if not _is_supported(input_ext):
            _print_supported_error(input_ext, "reading")
            return
        
        if not _is_supported(output_ext):
            _print_supported_error(output_ext, "writing")
            return

        try:
            if input_ext == "atraj":
                atoms, _ = read_atraj(input_file)
            elif input_ext == "tess":
                atoms, _ = read_tess(input_file)
            else:
                atoms = read(input_file, index=":")
        except UnknownFileTypeError:
            console.print(f"[red]Error: Unknown file type for input file '{input_file}'[/red]")
            console.print(f"[yellow]The file extension '{input_ext}' is not recognized.[/yellow]")
            console.print("[yellow]Make sure the file has the correct extension for its format.[/yellow]")
            return
        except Exception as e:
            console.print(f"[red]Error reading input file '{input_file}': {str(e)}.[/red]")
            console.print(f"[yellow]Make sure '{input_ext}' is a valid format and the file is not corrupted.[/yellow]")
            return
        
        try:
            if output_ext == "atraj":
                write_atraj(atoms, output_file)
            elif output_ext == "tess":
                write_tess(atoms, output_file)
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

