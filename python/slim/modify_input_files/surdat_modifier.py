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
    defaults = get_config_value(
        config=config,
        section=section,
        item="defaults",
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

    # dictionary of entries to loop over
    # "variable name": [type, allowed_values, index]
    vars_3d = {
        "glc_mask": [int, "from_file", 0],
        "alb_gvd": [float, None, 1],
        "alb_svd": [float, None, 2],
        "alb_gnd": [float, None, 3],
        "alb_snd": [float, None, 4],
        "alb_gvf": [float, None, 5],
        "alb_svf": [float, None, 6],
        "alb_gnf": [float, None, 7],
        "alb_snf": [float, None, 8],
        "bucketdepth": [float, None, 9],
        "emissivity": [float, None, 10],
        "snowmask": [float, None, 11],
        "roughness": [float, None, 12],
        "evap_res": [float, None, 13],
        "soil_type": [int, "from_file", 14],
        "soil_tk_1d": [float, None, 15],
        "soil_cv_1d": [float, None, 16],
        "glc_tk_1d": [float, None, 17],
        "glc_cv_1d": [float, None, 18],
    }
    # initialize entry
    entry = [None, None, None, None, None, None, None, None, None, None, None, None] * len(vars_3d)
    # not required: user may set these in the .cfg file
    for var, val in vars_3d.items():
        # obtain allowed values from surdat_in variable directly
        # TODO prefer to obtain from surdat_in variable's metadata which will
        #      contain more accurate information
        if val[1] is not None:
            allowed = modify_surdat.file[var]
        else:
            allowed = None
        # obtain user-defined values from the configure file
        entry[val[2]] = get_config_value(
            config=config,
            section=section,
            item=var,
            file_path=cfg_path,
            allowed_values=allowed,
            is_list=True,
            convert_to_type=val[0],
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
    if defaults:
        modify_surdat.set_defaults()  # set 3D variables
        logger.info("defaults complete")

    # User-selected values will overwrite either
    # - set_default's values if defaults = True or
    # - the input surdat's values if defaults = False

    for var, val in vars_3d.items():
        if entry[val[2]] is not None:
            modify_surdat.set_monthly_values(var=var, val=entry[val[2]])

    # ----------------------------------------------
    # Output the now modified SLIM surface data file
    # ----------------------------------------------
    write_output(modify_surdat.file, surdat_in, surdat_out, "surdat")
