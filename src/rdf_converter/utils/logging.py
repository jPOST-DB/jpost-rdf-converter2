from rich.console import Console
from rich.logging import RichHandler
import logging
import sys

_console = Console(file=sys.stdout)

def get_logger(name: str = 'jpconverter') -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        datefmt='[%X]',
        handlers=[RichHandler(console=_console, markup=True)],
    )
    logger = logging.getLogger(name)
    return logger

def console() -> Console:
    return _console
