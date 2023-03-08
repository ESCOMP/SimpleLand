"""
Run this code by using the following wrapper script:
/tools/mksurdat/mksurdat

The wrapper script and its README include a full description and instructions.
"""

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
        # TODO Later will need a calculation for each surdat var and will want
        #      to assign the 3d vars that include the time dimension
        vars_out = {
            "glc_mask": "TSA",
            "alb_gvd": "TSA",
            "alb_svd": "TSA",
            "alb_gnd": "TSA",
            "alb_snd": "TSA",
            "alb_gvf": "TSA",
            "alb_svf": "TSA",
            "alb_gnf": "TSA",
            "alb_snf": "TSA",
            "bucketdepth": "TSA",
            "emissivity": "TSA",
            "snowmask": "TSA",
            "roughness": "TSA",
            "evap_res": "TSA",
            "l2xavg_Fall_flxdst1": "TSA",
            "l2xavg_Fall_flxdst2": "TSA",
            "l2xavg_Fall_flxdst3": "TSA",
            "l2xavg_Fall_flxdst4": "TSA",
            "soil_type": "TSA",
            "soil_tk_1d": "TSA",
            "soil_cv_1d": "TSA",
            "glc_tk_1d": "TSA",
            "glc_cv_1d": "TSA",
            "lsmlon": "lon",
            "lsmlat": "lat",
            "time": "time",  # TODO Replace with months 1 to 12
        }

        for var, val in vars_out.items():
            logger.info('processing %s, %s', var, val)
            # TODO Why self needed twice for this to work?
            self.setvar(self, var, val)

        # TODO Clear unnecessary variables from file
