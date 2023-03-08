"""
Run this code by using the following wrapper script:
/tools/mksurdat/mksurdat

The wrapper script and its README include a full description and instructions.
"""

import os
import logging

import xarray as xr

logger = logging.getLogger(__name__)


class HistoryAverage:
    """
    Description
    -----------
    """

    def __init__(self, my_data):

        self.file = my_data

        self.months = 12  # number of months (always 12)

    @classmethod
    def init_from_file(cls, history_in):
        """Initialize a HistoryAverage object from file history_in"""
        logger.info("Opening history_in file to be averaged: %s", history_in)
        my_file = xr.open_dataset(history_in)
        return cls(my_file)

    @staticmethod
    def setvar(self, var, val):
        """
        Sets variable var to variable val
        """
        self.file[var] = self.file[val]

    def set_surdat_vars(self):
        """
        Description
        -----------
        Set surdat variables
        """

        # Default values of 3d variables. For guidance in selecting values, see
        # /glade/p/cesmdata/cseg/inputdata/lnd/slim/surdat/
        # globalconst_alpha0.2_soilcv2e6_hc0.1_rs100.0_glc_hc0.01_f19_cdf5_20211105.nc
        # Dimensions are time,lsmlat,lsmlon
        # Dictionary of the variables that we will loop over
        # Same as the list found in modify_input_files/surdat_modifier.py but
        # with four l2xavg_Fall_flxdst variables added
        # TODO Currently assigning history vars to surdat vars as a test
        # TODO Later will need a calculation for each surdat var
        vars_out = {
            "glc_mask": "TSA",
#           "alb_gvd": 
#           "alb_svd": 
#           "alb_gnd": 
#           "alb_snd": 
#           "alb_gvf": 
#           "alb_svf": 
#           "alb_gnf": 
#           "alb_snf": 
#           "bucketdepth": 
#           "emissivity": 
#           "snowmask": 
#           "roughness": 
#           "evap_res": 
#           "l2xavg_Fall_flxdst1": 
#           "l2xavg_Fall_flxdst2": 
#           "l2xavg_Fall_flxdst3": 
#           "l2xavg_Fall_flxdst4": 
#           "soil_type": 
#           "soil_tk_1d": 
#           "soil_cv_1d": 
#           "glc_tk_1d": 
#           "glc_cv_1d": 
            "lsmlon": "lon",
            "lsmlat": "lat",
#           "time": ,
        }

        for var, val in vars_out.items():
            logger.info('processing %s, %s', var, val)
            # TODO Why self needed twice for this to work?
            self.setvar(self, var, val)

        # TODO Clear unnecessary variables from file
