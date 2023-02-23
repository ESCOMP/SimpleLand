"""
Run this code by using the following wrapper script:
/tools/modify_input_files/surdat_modifier

The wrapper script includes a full description and instructions.
"""

import os
import logging

from math import isclose
import numpy as np
import xarray as xr

from slim.utils import abort, update_metadata
from slim.git_utils import get_slim_git_short_hash
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

    def write_output(self, surdat_in, surdat_out):
        """
        Description
        -----------
        Write output file

        Arguments
        ---------
        surdat_in:
            (str) Command line entry of input surface dataset
        surdat_out:
            (str) Command line entry of output surface dataset
        """

        # update attributes
        title = "Modified SLIM surdat file"
        summary = "Modified SLIM surdat file"
        contact = "N/A"
        data_script = os.path.abspath(__file__) + " -- " + get_slim_git_short_hash()
        description = "Modified this file: " + surdat_in
        update_metadata(
            self.file,
            title=title,
            summary=summary,
            contact=contact,
            data_script=data_script,
            description=description,
        )

        # abort if output file already exists
        file_exists = os.path.exists(surdat_out)
        if file_exists:
            errmsg = "Output file already exists: " + surdat_out
            abort(errmsg)

        # mode 'w' overwrites file if it exists
        self.file.to_netcdf(path=surdat_out, mode="w", format="NETCDF3_64BIT")
        logger.info("Successfully created surdat_out: %s", surdat_out)
        self.file.close()

    def set_monthly_values(self, var, val):
        """
        Description
        -----------
        If user has specified monthly values, use them. Else do nothing.
        """
        if len(val) != self.months:
            errmsg = (
                "Error: Variable should have exactly "
                + self.months
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

        See ctsm subdirectory /python/ctsm/modify_input_files,
        file modify_fsurdat.py for templates of the corresponding
        functions that set 2d and 4d variables in the user-defined rectangle
        """
        self.file[var][lev1_dim, ...] = self.file[var][lev1_dim, ...].where(
            self.not_rectangle, other=val
        )

    def set_idealized(self):
        """
        Description
        -----------
        Set surdat variables in a rectangle defined by lon/lat limits
        """

        # Overwrite in rectangle(s)
        # ------------------------
        # If idealized, the user makes changes to variables as follows.
        # Values in the user-defined rectangle are replaced.
        # Values outside the rectangle are preserved.
        # ------------------------

        # Default values
        glc_mask = [0] * self.months
        alb_gvd = [0] * self.months
        alb_svd = [0] * self.months
        alb_gnd = [0] * self.months
        alb_snd = [0] * self.months
        alb_gvf = [0] * self.months
        alb_svf = [0] * self.months
        alb_gnf = [0] * self.months
        alb_snf = [0] * self.months
        bucketdepth = [0] * self.months
        emissivity = [0] * self.months
        snowmask = [0] * self.months
        roughness = [0] * self.months
        evap_res = [0] * self.months
        l2xavg_Fall_flxdst1 = [0] * self.months
        l2xavg_Fall_flxdst2 = [0] * self.months
        l2xavg_Fall_flxdst3 = [0] * self.months
        l2xavg_Fall_flxdst4 = [0] * self.months
        soil_type = [0] * self.months
        soil_tk_1d = [0] * self.months
        soil_cv_1d = [0] * self.months
        glc_tk_1d = [0] * self.months
        glc_cv_1d = [0] * self.months

        # dictionary of 3d variables to loop over
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
            "l2xavg_Fall_flxdst1": l2xavg_Fall_flxdst1,
            "l2xavg_Fall_flxdst2": l2xavg_Fall_flxdst2,
            "l2xavg_Fall_flxdst3": l2xavg_Fall_flxdst3,
            "l2xavg_Fall_flxdst4": l2xavg_Fall_flxdst4,
            "soil_type": soil_type,
            "soil_tk_1d": soil_tk_1d,
            "soil_cv_1d": soil_cv_1d,
            "glc_tk_1d": glc_tk_1d,
            "glc_cv_1d": glc_cv_1d,
        }
        for var, val in vars_3d.items():
            if val is not None:
                self.set_monthly_values(var=var, val=val)
