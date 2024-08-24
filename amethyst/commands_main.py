from typing import Optional
from loguru import logger

import click

from amethyst.amethyst import enumerate



@click.group()
def main():
    pass


@main.command()
@click.argument("core_smi", nargs=1, required=True, type=str)
@click.argument(
    "r_file", nargs=1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option("-d", "--delimiter", default=",", show_default=True, type=str, help="Symbol separating R-groups.") # set to , so CSV by default
@click.option("-e", "--enantiomers", is_flag=True, default=False, help="Option to generate possible enantiomers.")
@click.option("-o", "--output-file", type=click.File(mode="x"), help="File output path. Provided filepath must be empty.")
@click.option("--debug", is_flag=True, default=False, help="Enable debug logs. Funnily enough currently broken so it always debug logs.")
def generate(
    core_smi: str,
    r_file: str,
    delimiter: Optional[str] = ",",
    enantiomers: bool = False,
    output_file: Optional[str] = None,
    debug: Optional[bool] = False
) -> None:
    """Generate combinations of provided CORE_SMI and R_FILE."""
    
    if debug:
        logger.add("amethyst.log", level="DEBUG")
    else:
        logger.add("amethyst.log")
    if output_file is None:
        output_file = "output.txt"

    with open(r_file, "r") as f:
        r_num: int = len(f.readlines())
        if r_num > 1:
            multiple_rs = True
        elif r_num == 1:
            multiple_rs = False
        else:
            raise EOFError("File is empty")
    
    logger.debug(f"Core: {core_smi}")
    logger.debug(f"R file: {r_file}")
    logger.debug(f"Delimiter: {delimiter}") # FIXME - it works only with -d";" for ; delim????
    logger.debug(f"Enantiomers: {enantiomers}")
    logger.debug(f"Output file: {output_file}")

    output: str = ",".join(
        enumerate(
            core_smi, delimiter=delimiter, subs_path=r_file, multiple_rs=multiple_rs, enantiomers=enantiomers, output_smi=True
        )
    )

    with open(output_file, "w") as f:
        f.write(output)
        f.close()
