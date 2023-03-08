"""
Run this code by using the following wrapper script:
tools/mksurdat/mksurdat

The wrapper script and its README include a full description and instructions.
"""

import os
import logging
import argparse
from configparser import ConfigParser

from slim.utils import abort, write_output
from slim.config_utils import get_config_value
from slim.slim_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)
from slim.mksurdat.history_average import HistoryAverage

logger = logging.getLogger(__name__)


def main():
    """
    Description
    -----------
    Calls function that generates a SLIM surface dataset (surdat file)
    """

    # set up logging allowing user control
    setup_logging_pre_config()

    # read the command line argument to obtain the path to the .cfg file
    parser = argparse.ArgumentParser()
    parser.add_argument("cfg_path", help="/path/name.cfg of input file, eg ./mksurdat.cfg")
    add_logging_args(parser)
    args = parser.parse_args()
    process_logging_args(args)
    mksurdat(args.cfg_path)


def mksurdat(cfg_path):
    """Implementation of mksurdat command"""
    # read the .cfg (config) file
    config = ConfigParser()
    config.read(cfg_path)
    section = config.sections()[0]  # name of the first section

    # required: user must set these in the .cfg file
    history_in = get_config_value(
        config=config, section=section, item="history_in", file_path=cfg_path
    )
    surdat_out = get_config_value(
        config=config, section=section, item="surdat_out", file_path=cfg_path
    )

    # required but fallback values available for variables omitted
    # entirely from the .cfg file
    num_of_hist_yrs = get_config_value(
        config=config,
        section=section,
        item="num_of_hist_yrs",
        file_path=cfg_path,
        allowed_values=range(1, 100),  # 1 to 99
        convert_to_type=int,
    )

    # Create HistoryAverage object
    # TODO Currently this reads the first history file only
    history_average = HistoryAverage.init_from_file(
        history_in,
    )

    # If output file exists, abort before starting work
    if os.path.exists(surdat_out):
        errmsg = "Output file already exists: " + surdat_out
        abort(errmsg)

    # --------------------------------
    # generate surface data properties
    # --------------------------------

    history_average.set_surdat_vars()
    logger.info("set_surdat_vars complete")

    # -------------------------------------
    # Output the new SLIM surface data file
    # -------------------------------------
    write_output(history_average.file, history_in, surdat_out, "surdat")
