"""
Run this code by using the following wrapper script:
tools/modify_input_files/surdat_modifier

The wrapper script includes a full description and instructions.
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
from slim.modify_input_files.modify_surdat import ModifySurdat

logger = logging.getLogger(__name__)


def main():
    """
    Description
    -----------
    Calls function that modifies a surdat file (surface dataset)
    """

    # set up logging allowing user control
    setup_logging_pre_config()

    # read the command line argument to obtain the path to the .cfg file
    parser = argparse.ArgumentParser()
    parser.add_argument("cfg_path", help="/path/name.cfg of input file, eg ./modify.cfg")
    add_logging_args(parser)
    args = parser.parse_args()
    process_logging_args(args)
    surdat_modifier(args.cfg_path)


def surdat_modifier(cfg_path):
    """Implementation of surdat_modifier command"""
    # read the .cfg (config) file
    config = ConfigParser()
    config.read(cfg_path)
    section = config.sections()[0]  # name of the first section

    # required: user must set these in the .cfg file
    surdat_in = get_config_value(
        config=config, section=section, item="surdat_in", file_path=cfg_path
    )
    surdat_out = get_config_value(
        config=config, section=section, item="surdat_out", file_path=cfg_path
    )

    # required but fallback values available for variables omitted
    # entirely from the .cfg file
    idealized = get_config_value(
        config=config,
        section=section,
        item="idealized",
        file_path=cfg_path,
        convert_to_type=bool,
    )
    lnd_lat_1 = get_config_value(
        config=config,
        section=section,
        item="lnd_lat_1",
        file_path=cfg_path,
        convert_to_type=float,
    )
    lnd_lat_2 = get_config_value(
        config=config,
        section=section,
        item="lnd_lat_2",
        file_path=cfg_path,
        convert_to_type=float,
    )
    lnd_lon_1 = get_config_value(
        config=config,
        section=section,
        item="lnd_lon_1",
        file_path=cfg_path,
        convert_to_type=float,
    )
    lnd_lon_2 = get_config_value(
        config=config,
        section=section,
        item="lnd_lon_2",
        file_path=cfg_path,
        convert_to_type=float,
    )

    landmask_file = get_config_value(
        config=config,
        section=section,
        item="landmask_file",
        file_path=cfg_path,
        can_be_unset=True,
    )

    lat_dimname = get_config_value(
        config=config, section=section, item="lat_dimname", file_path=cfg_path, can_be_unset=True
    )
    lon_dimname = get_config_value(
        config=config, section=section, item="lon_dimname", file_path=cfg_path, can_be_unset=True
    )

    # Create ModifySurdat object
    modify_surdat = ModifySurdat.init_from_file(
        surdat_in,
        lnd_lon_1,
        lnd_lon_2,
        lnd_lat_1,
        lnd_lat_2,
        landmask_file,
        lat_dimname,
        lon_dimname,
    )

    # If output file exists, abort before starting work
    if os.path.exists(surdat_out):
        errmsg = "Output file already exists: " + surdat_out
        abort(errmsg)

    # not required: user may set these in the .cfg file
#   max_soil_type = int(modify_surdat.file.mxsoil_type)
    soil_type = get_config_value(
        config=config,
        section=section,
        item="soil_type",
        file_path=cfg_path,
#       allowed_values=range(1, max_soil_type + 1),  # 1 to max_soil_type
        convert_to_type=int,
        can_be_unset=True,
    )

    soil_tk_1d = get_config_value(
        config=config,
        section=section,
        item="soil_tk_1d",
        file_path=cfg_path,
        convert_to_type=float,
        can_be_unset=True,
    )
    soil_cv_1d = get_config_value(
        config=config,
        section=section,
        item="soil_cv_1d",
        file_path=cfg_path,
        convert_to_type=float,
        can_be_unset=True,
    )

    # ------------------------------
    # modify surface data properties
    # ------------------------------

    # Set surdat variables in a rectangle that could be global (default).
    # Note that the land/ocean mask gets specified in the domain file for
    # MCT or the ocean mesh files for NUOPC. Here the user may specify
    # surdat variables inside a box but cannot change which points will
    # run as land and which as ocean.
    if idealized:
        modify_surdat.set_idealized()  # set 3D variables
        logger.info("idealized complete")

    # TODO slevis
    # Move vars_3d before if idealized and pass as arg into set_idealized
    # Repeat var, val loop from modify_surdat.py here

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    write_output(modify_surdat.file, surdat_in, surdat_out, "surdat")
