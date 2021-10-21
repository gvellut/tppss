import logging
import sys

import click
from colorama import Fore, Style

logger = logging.getLogger(__package__)

# specify colors for different logging levels
LOG_COLORS = {logging.ERROR: Fore.RED, logging.WARNING: Fore.YELLOW}


class ColorFormatter(logging.Formatter):
    def format(self, record, *args, **kwargs):
        if record.levelno in LOG_COLORS:
            record.msg = "{color_begin}{message}{color_end}".format(
                message=record.msg,
                color_begin=LOG_COLORS[record.levelno],
                color_end=Style.RESET_ALL,
            )
        return super().format(record, *args, **kwargs)


def setup_logging(is_debug):
    global logger
    if is_debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = ColorFormatter("%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


@click.group()
@click.option(
    "--debug",
    "is_debug",
    is_flag=True,
    help=("Flag to activate debug mode"),
    required=False,
)
@click.pass_context
def main(ctx, is_debug):
    """Computes sunset / sunrise time taking into account local topography"""
    setup_logging(is_debug)
    # special attribute of context
    ctx.obj = {"DEBUG": is_debug}


@main.command("day")
@click.pass_context
def tpss_day(ctx):
    """Computes sunset / sunrise time for a single day"""
    # TODO multiple days at once ?
    pass


@main.command("year")
@click.pass_context
def tpss_year(ctx):
    """Computes sunset / sunrise time for a whole year"""
    pass


if __name__ == "__main__":
    main()
