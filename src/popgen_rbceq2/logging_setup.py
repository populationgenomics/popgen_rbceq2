"""Logging configuration shared by the workflow entry point and the job scripts."""

import logging


def setup_logging(level: int = logging.INFO, force: bool = False) -> None:
    """Configure standard logging for the pipeline.

    Args:
        level: Minimum level to emit.
        force: Replace handlers already installed on the root logger. Needed in a Hail job,
            where an imported library may have configured logging first.
    """
    logging.basicConfig(
        format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=level,
        force=force,
    )
