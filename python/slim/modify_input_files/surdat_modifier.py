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
    glc_mask = get_config_value(
        config=config,
        section=section,
        item="glc_mask",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_gvd = get_config_value(
        config=config,
        section=section,
        item="alb_gvd",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_svd = get_config_value(
        config=config,
        section=section,
        item="alb_svd",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_gnd = get_config_value(
        config=config,
        section=section,
        item="alb_gnd",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_snd = get_config_value(
        config=config,
        section=section,
        item="alb_snd",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_gvf = get_config_value(
        config=config,
        section=section,
        item="alb_gvf",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_svf = get_config_value(
        config=config,
        section=section,
        item="alb_svf",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_gnf = get_config_value(
        config=config,
        section=section,
        item="alb_gnf",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    alb_snf = get_config_value(
        config=config,
        section=section,
        item="alb_snf",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    bucketdepth = get_config_value(
        config=config,
        section=section,
        item="bucketdepth",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    emissivity = get_config_value(
        config=config,
        section=section,
        item="emissivity",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    snowmask = get_config_value(
        config=config,
        section=section,
        item="snowmask",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    roughness = get_config_value(
        config=config,
        section=section,
        item="roughness",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    evap_res = get_config_value(
        config=config,
        section=section,
        item="evap_res",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
#   max_soil_type = int(modify_surdat.file.mxsoil_type)
    soil_type = get_config_value(
        config=config,
        section=section,
        item="soil_type",
        file_path=cfg_path,
#       allowed_values=range(1, max_soil_type + 1),  # 1 to max_soil_type
        is_list=True,
        convert_to_type=int,
        can_be_unset=True,
    )
    soil_tk_1d = get_config_value(
        config=config,
        section=section,
        item="soil_tk_1d",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    soil_cv_1d = get_config_value(
        config=config,
        section=section,
        item="soil_cv_1d",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    glc_tk_1d = get_config_value(
        config=config,
        section=section,
        item="glc_tk_1d",
        file_path=cfg_path,
        is_list=True,
        convert_to_type=float,
        can_be_unset=True,
    )
    glc_cv_1d = get_config_value(
        config=config,
        section=section,
        item="glc_cv_1d",
        file_path=cfg_path,
        is_list=True,
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

    # User-selected values will overwrite either
    # - set_idealized's default values if idealized = True or
    # - the input surdat's values if idealized = False
    # Dictionary of 3d variables to loop over
    vars_3d = {
        "glc_mask": glc_mask,
        "alb_gvd": alb_gvd,
        "alb_svd": alb_svd,
        "alb_gnd": alb_gnd,
        "alb_snd": alb_snd,
        "alb_gvf": alb_gvf,
        "alb_svf": alb_svf,
        "alb_gnf": alb_gnf,
        "alb_snf": alb_snf,
        "bucketdepth": bucketdepth,
        "emissivity": emissivity,
        "snowmask": snowmask,
        "roughness": roughness,
        "evap_res": evap_res,
        "soil_type": soil_type,
        "soil_tk_1d": soil_tk_1d,
        "soil_cv_1d": soil_cv_1d,
        "glc_tk_1d": glc_tk_1d,
        "glc_cv_1d": glc_cv_1d,
    }

    for var, val in vars_3d.items():
        if val is not None:
            modify_surdat.set_monthly_values(var=var, val=val)

    # ----------------------------------------------
    # Output the now modified CTSM surface data file
    # ----------------------------------------------
    write_output(modify_surdat.file, surdat_in, surdat_out, "surdat")
