"""
Run this code by using the following wrapper script:
/tools/modify_input_files/surdat_modifier

The wrapper script includes a full description and instructions.
"""

import logging

from math import isclose
import numpy as np
import xarray as xr

from slim.utils import abort
from slim.config_utils import lon_range_0_to_360

logger = logging.getLogger(__name__)


class ModifySurdat:
    """
    Description
    -----------
    """

    def __init__(
        self, my_data, lon_1, lon_2, lat_1, lat_2, landmask_file, lat_dimname, lon_dimname
    ):

        self.file = my_data

        self.rectangle = self._get_rectangle(
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            longxy=self.file.lsmlon,
            latixy=self.file.lsmlat,
        )

        if landmask_file is not None:
            # overwrite self.not_rectangle with data from
            # user-specified .nc file in the .cfg file
            landmask_ds = xr.open_dataset(landmask_file)
            self.rectangle = landmask_ds.mod_lnd_props.data
            # CF convention has dimension and coordinate variable names the same
            if lat_dimname is None:  # set to default
                lat_dimname = "lsmlat"
            if lon_dimname is None:  # set to default
                lon_dimname = "lsmlon"
            lsmlat = landmask_ds.dims[lat_dimname]
            lsmlon = landmask_ds.dims[lon_dimname]

            for row in range(lsmlat):  # rows from landmask file
                for col in range(lsmlon):  # cols from landmask file
                    errmsg = (
                        "landmask_ds.mod_lnd_props not 0 or 1 at "
                        + f"row, col, value = {row} {col} {self.rectangle[row, col]}"
                    )
                    assert isclose(self.rectangle[row, col], 0, abs_tol=1e-9) or isclose(
                        self.rectangle[row, col], 1, abs_tol=1e-9
                    ), errmsg

        self.not_rectangle = np.logical_not(self.rectangle)
        self.months = int(max(self.file.time))  # number of months (typically 12)

    @classmethod
    def init_from_file(
        cls, surdat_in, lon_1, lon_2, lat_1, lat_2, landmask_file, lat_dimname, lon_dimname
    ):
        """Initialize a ModifySurdat object from file surdat_in"""
        logger.info("Opening surdat_in file to be modified: %s", surdat_in)
        my_file = xr.open_dataset(surdat_in)
        return cls(my_file, lon_1, lon_2, lat_1, lat_2, landmask_file, lat_dimname, lon_dimname)

    @staticmethod
    def _get_rectangle(lon_1, lon_2, lat_1, lat_2, longxy, latixy):
        """
        Description
        -----------
        """

        # ensure that lon ranges 0-360 in case user entered -180 to 180
        lon_1 = lon_range_0_to_360(lon_1)
        lon_2 = lon_range_0_to_360(lon_2)

        # determine the rectangle(s)
        # TODO This is not really "nearest" for the edges but isel didn't work
        rectangle_1 = longxy >= lon_1
        rectangle_2 = longxy <= lon_2
        eps = np.finfo(np.float32).eps  # to avoid roundoff issue
        rectangle_3 = latixy >= (lat_1 - eps)
        rectangle_4 = latixy <= (lat_2 + eps)

        if lon_1 <= lon_2:
            # rectangles overlap
            union_1 = np.logical_and(rectangle_1, rectangle_2)
        else:
            # rectangles don't overlap: stradling the 0-degree meridian
            union_1 = np.logical_or(rectangle_1, rectangle_2)

        if lat_1 < -90 or lat_1 > 90 or lat_2 < -90 or lat_2 > 90:
            errmsg = "lat_1 and lat_2 need to be in the range -90 to 90"
            abort(errmsg)
        elif lat_1 <= lat_2:
            # rectangles overlap
            union_2 = np.logical_and(rectangle_3, rectangle_4)
        else:
            # rectangles don't overlap: one in the north, one in the south
            union_2 = np.logical_or(rectangle_3, rectangle_4)

        # union rectangles overlap
        rectangle = np.logical_and(union_1, union_2)

        return rectangle

    def set_monthly_values(self, var, val):
        """
        Description
        -----------
        If user has specified monthly values, use them. Else do nothing.
        """
        if len(val) != self.months:
            errmsg = (
                "Error: Variable should have exactly "
                + str(self.months)
                + " entries in the configure file: "
                + var
            )
            abort(errmsg)
        for mon in self.file.time - 1:  # loop over the months
            # set 3D variable
            self.setvar_lev1(var, val[int(mon)], lev1_dim=mon)

    def setvar_lev1(self, var, val, lev1_dim):
        """
        Sets 3d variable var to value val in user-defined rectangle,
        defined as "other" in the function

        HINT for working with 2d or 4d variables instead:
        See ctsm subdirectory /python/ctsm/modify_input_files,
        file modify_fsurdat.py for templates of the corresponding functions
        """
        self.file[var][lev1_dim, ...] = self.file[var][lev1_dim, ...].where(
            self.not_rectangle, other=val
        )

    def set_defaults(self):
        """
        Description
        -----------
        Set default surdat values in a rectangle defined by lon/lat limits
        """

        # Overwrite in rectangle(s)
        # ------------------------
        # If defaults, then user makes changes to variables as follows.
        # Values in the user-defined rectangle are replaced.
        # Values outside the rectangle are preserved.
        # ------------------------

        # Default values of 3d variables. For guidance in selecting values, see
        # /glade/p/cesmdata/cseg/inputdata/lnd/slim/surdat/
        # globalconst_alpha0.2_soilcv2e6_hc0.1_rs100.0_glc_hc0.01_f19_cdf5_20211105.nc
        # Dimensions are time,lsmlat,lsmlon
        # Dictionary of the variables that we will loop over
        vars_3d = {
            "glc_mask": [0] * self.months,
            "alb_gvd": [0.2] * self.months,
            "alb_svd": [0.8] * self.months,
            "alb_gnd": [0.3] * self.months,
            "alb_snd": [0.6] * self.months,
            "alb_gvf": [0.2] * self.months,
            "alb_svf": [0.8] * self.months,
            "alb_gnf": [0.3] * self.months,
            "alb_snf": [0.6] * self.months,
            "bucketdepth": [200] * self.months,
            "emissivity": [1] * self.months,
            "snowmask": [50] * self.months,
            "roughness": [0.1] * self.months,
            "evap_res": [100] * self.months,
            "soil_type": [0] * self.months,
            "soil_tk_1d": [1.5] * self.months,
            "soil_cv_1d": [2e6] * self.months,
            "glc_tk_1d": [2.4] * self.months,
            "glc_cv_1d": [1.9e6] * self.months,
        }

        for var, val in vars_3d.items():
            self.set_monthly_values(var=var, val=val)
