#!/usr/bin/env python3

"""System tests for surdat_modifier

"""

import os
import re

import unittest
import tempfile
import shutil

import numpy as np
import xarray as xr

from slim.path_utils import path_to_slim_root
from slim.config_utils import lon_range_0_to_360
from slim.utils import write_output
from slim import unit_testing
from slim.modify_input_files.surdat_modifier import surdat_modifier

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysSurdatModifier(unittest.TestCase):
    """System tests for surdat_modifier"""

    def setUp(self):
        """
        Obtain path to the existing modify_surdat_template.cfg file
        Make /_tempdir for use by these tests
        Obtain path and names for the files being created in /_tempdir:
        - modify_surdat.cfg
        - surdat_out.nc
        - surdat_in.nc
        Generate dummy surdat_in file and save
        Come up with modifications to be introduced to surdat_in
        """
        self._cfg_template_path = os.path.join(
            path_to_slim_root(), "tools/modify_input_files/modify_surdat_template.cfg"
        )
        self._tempdir = tempfile.mkdtemp()
        self._cfg_file_path = os.path.join(self._tempdir, "modify_surdat.cfg")
        self._surdat_out = os.path.join(self._tempdir, "surdat_out.nc")
        self._surdat_in = os.path.join(self._tempdir, "surdat_in.nc")
        months = 12

        # -----------------------------------------------------------
        # create dummy SLIM surdat file
        # -----------------------------------------------------------
        # get lon/lat that would normally come from a surdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        self._lon_range = [2, 10]  # expected in ascending order: [min, max]
        self._lat_range = [3, 12]  # expected in ascending order: [min, max]
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min(self._lon_range),
            _max_lon=max(self._lon_range),
            _min_lat=min(self._lat_range),
            _max_lat=max(self._lat_range),
        )
        lon_1d = longxy[0, :]
        lat_1d = latixy[:, 0]
        # create xarray dataset containing lev1 variables;
        # the surdat_modify tool reads variables like this from a surdat file
        var_1d = np.arange(cols)
        ones_3d = np.ones((months, rows, cols))
        var_lev1 = var_1d * ones_3d
        self._surdat_in_data = xr.Dataset(
            data_vars=dict(
                time=(["time"], np.arange(months) + 1),
                lsmlon=(["lsmlon"], lon_1d),
                lsmlat=(["lsmlat"], lat_1d),
                glc_mask=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_gvd=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_svd=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_gnd=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_snd=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_gvf=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_svf=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_gnf=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_snf=(["time", "lsmlat", "lsmlon"], var_lev1),
                bucketdepth=(["time", "lsmlat", "lsmlon"], var_lev1),
                emissivity=(["time", "lsmlat", "lsmlon"], var_lev1),
                snowmask=(["time", "lsmlat", "lsmlon"], var_lev1),
                roughness=(["time", "lsmlat", "lsmlon"], var_lev1),
                evap_res=(["time", "lsmlat", "lsmlon"], var_lev1),
                l2xavg_Fall_flxdst1=(["time", "lsmlat", "lsmlon"], var_lev1),
                l2xavg_Fall_flxdst2=(["time", "lsmlat", "lsmlon"], var_lev1),
                l2xavg_Fall_flxdst3=(["time", "lsmlat", "lsmlon"], var_lev1),
                l2xavg_Fall_flxdst4=(["time", "lsmlat", "lsmlon"], var_lev1),
                soil_type=(["time", "lsmlat", "lsmlon"], var_lev1),
                soil_tk_1d=(["time", "lsmlat", "lsmlon"], var_lev1),
                soil_cv_1d=(["time", "lsmlat", "lsmlon"], var_lev1),
                glc_tk_1d=(["time", "lsmlat", "lsmlon"], var_lev1),
                glc_cv_1d=(["time", "lsmlat", "lsmlon"], var_lev1),
            )
        )
        # Add attributes to all the variables
        attr_map = {
            "glc_mask": ["unitless", 1e36, "Glacier/ice sheet mask", [0, 1]],
            "alb_gvd": ["unitless", 1e36, "Visible direct albedo for bare ground", []],
            "alb_svd": ["unitless", 1e36, "Visible direct albedo for deep snow", []],
            "alb_gnd": ["unitless", 1e36, "NIR direct albedo for bare ground", []],
            "alb_snd": ["unitless", 1e36, "NIR direct albedo for deep snow", []],
            "alb_gvf": ["unitless", 1e36, "Visible diffuse albedo for bare ground", []],
            "alb_svf": ["unitless", 1e36, "Visible diffuse albedo for deep snow", []],
            "alb_gnf": ["unitless", 1e36, "NIR diffuse albedo for bare ground", []],
            "alb_snf": ["unitless", 1e36, "NIR diffuse albedo for deep snow", []],
            "bucketdepth": ["kg/m2", 1e36, "Bucket capacity", []],
            "emissivity": ["unitless", 1e36, "Surface emissivity for longwave radiation", []],
            "snowmask": ["kg/m2", 1e36, "Snow-masking depth", []],
            "roughness": ["m", 1e36, "Vegetation height", []],
            "evap_res": ["s/m", 1e36, "Evaporative resistance", []],
            "l2xavg_Fall_flxdst1": ["unknown", 1e36, "Dust flux", []],
            "l2xavg_Fall_flxdst2": ["unknown", 1e36, "Dust flux", []],
            "l2xavg_Fall_flxdst3": ["unknown", 1e36, "Dust flux", []],
            "l2xavg_Fall_flxdst4": ["unknown", 1e36, "Dust flux", []],
            "soil_type": ["unitless", 1e36, "Soil type (unused)", [0]],
            "soil_tk_1d": ["W/m/K", 1e36, "Soil thermal conductivity", []],
            "soil_cv_1d": ["J/m3/K", 1e36, "Soil heat capacity", []],
            "glc_tk_1d": ["W/m/K", 1e36, "Ice thermal conductivity", []],
            "glc_cv_1d": ["J/m3/K", 1e36, "Ice heat capacity", []],
            "lsmlat": ["degrees north", False, "Coordinate latitude", []],
            "lsmlon": ["degrees east", False, "Coordinate longitude", []],
            "time": ["month", False, "", []],
        }
        for var, val in attr_map.items():
            self._surdat_in_data[var].attrs["Units"] = val[0]
            self._surdat_in_data[var].attrs["_FillValue"] = val[1]
            self._surdat_in_data[var].attrs["long_name"] = val[2]
            self._surdat_in_data[var].attrs["valid_range"] = val[3]

        # save in tempdir; _in and _out files are the same file in this case
        write_output(self._surdat_in_data, self._surdat_in, self._surdat_in, "surdat")
        # come up with modifications to be introduced to surdat_in
        self._modified_1 = ones_3d.astype(int)
        self._modified_2 = 0 * self._modified_1
        self._modified_3 = 0.5 * self._modified_1
        self._modified_4 = 195 * self._modified_1

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_minimalInfo(self):
        """
        This test specifies a minimal amount of information
        Create .cfg file, run the tool, compare surdat_in to surdat_out
        """

        self._create_config_file_minimal()

        # run the surdat_modifier tool
        surdat_modifier(self._cfg_file_path)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        surdat_out_data = xr.open_dataset(self._surdat_out)
        # assert that surdat_out equals surdat_in
        self.assertTrue(surdat_out_data.equals(self._surdat_in_data))

    def test_allInfo(self):
        """
        This version specifies all possible information
        Create .cfg file, run the tool, compare surdat_in to surdat_out
        Here also compare surdat_out to surdat_out_baseline
        """

        self._create_config_file_complete()

        # run the surdat_modifier tool
        surdat_modifier(self._cfg_file_path)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # compare surdat_out to surdat_in
        surdat_out_data = xr.open_dataset(self._surdat_out)
        # assert that surdat_out does not equal surdat_in
        self.assertFalse(surdat_out_data.equals(self._surdat_in_data))

        # -----------------------------------------------------------
        # compare surdat_out to surdat_out_baseline
        # -----------------------------------------------------------
        # generate surdat_out_baseline by merging surdat_in into the
        # modified dataset and compare to surdat_out
        modified_1_through_4 = xr.Dataset(
            data_vars=dict(
                glc_mask=(["time", "lsmlat", "lsmlon"], self._modified_1),
                alb_gvd=(["time", "lsmlat", "lsmlon"], self._modified_2),
                alb_svd=(["time", "lsmlat", "lsmlon"], self._modified_3),
                bucketdepth=(["time", "lsmlat", "lsmlon"], self._modified_4),
            )
        )
        surdat_out_base_data = modified_1_through_4.merge(self._surdat_in_data, compat="override")

        # assert that surdat_out equals surdat_out_baseline
        self.assertTrue(surdat_out_data.equals(surdat_out_base_data))

    def _create_config_file_minimal(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open(self._cfg_file_path, "w", encoding="utf-8") as cfg_out:
            with open(self._cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *surdat_in *=", line):
                        line = f"surdat_in = {self._surdat_in}"
                    elif re.match(r" *surdat_out *=", line):
                        line = f"surdat_out = {self._surdat_out}"
                    cfg_out.write(line)

    def _create_config_file_complete(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open(self._cfg_file_path, "w", encoding="utf-8") as cfg_out:
            with open(self._cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *surdat_in *=", line):
                        line = f"surdat_in = {self._surdat_in}"
                    elif re.match(r" *surdat_out *=", line):
                        line = f"surdat_out = {self._surdat_out}"
                    elif re.match(r" *defaults *=", line):
                        line = "defaults = False"
                    elif re.match(r" *lnd_lat_1 *=", line):
                        line = "lnd_lat_1 = " + str(min(self._lat_range)) + "\n"
                    elif re.match(r" *lnd_lat_2 *=", line):
                        line = "lnd_lat_2 = " + str(max(self._lat_range)) + "\n"
                    elif re.match(r" *lnd_lon_1 *=", line):
                        line = "lnd_lon_1 = " + str(min(self._lon_range)) + "\n"
                    elif re.match(r" *lnd_lon_2 *=", line):
                        line = "lnd_lon_2 = " + str(max(self._lon_range)) + "\n"
                    elif re.match(r" *glc_mask *=", line):
                        # in .cfg file user enters list of monthly (i.e. 12)
                        # values without punctuation (e.g. brackets or commas)
                        line = "glc_mask = " + str(self._modified_1[:, 0, 0])[1:-1] + "\n"
                    elif re.match(r" *alb_gvd *=", line):
                        # in .cfg file user enters list of monthly (i.e. 12)
                        # values without punctuation (e.g. brackets or commas)
                        line = "alb_gvd = " + str(self._modified_2[:, 0, 0])[1:-1] + "\n"
                    elif re.match(r" *alb_svd *=", line):
                        # in .cfg file user enters list of monthly (i.e. 12)
                        # values without punctuation (e.g. brackets or commas)
                        line = "alb_svd = " + str(self._modified_3[:, 0, 0])[1:-1] + "\n"
                    elif re.match(r" *bucketdepth *=", line):
                        # in .cfg file user enters list of monthly (i.e. 12)
                        # values without punctuation (e.g. brackets or commas)
                        line = "bucketdepth = " + str(self._modified_4[:, 0, 0])[1:-1] + "\n"
                    cfg_out.write(line)

    def _get_longxy_latixy(self, _min_lon, _max_lon, _min_lat, _max_lat):
        """
        Return longxy, latixy, cols, rows
        Function copied from test_unit_modify_surdat.py
        TODO Move to a separate file of test utilities?
        """
        cols = _max_lon - _min_lon + 1
        rows = _max_lat - _min_lat + 1

        long = np.arange(_min_lon, _max_lon + 1)
        long = [lon_range_0_to_360(longitude) for longitude in long]
        longxy = long * np.ones((rows, cols))
        compare = np.repeat([long], rows, axis=0)  # alternative way to form
        # assert this to confirm intuitive understanding of these matrices
        np.testing.assert_array_equal(longxy, compare)

        lati = np.arange(_min_lat, _max_lat + 1)
        self.assertEqual(min(lati), _min_lat)
        self.assertEqual(max(lati), _max_lat)
        latixy_transp = lati * np.ones((cols, rows))
        compare = np.repeat([lati], cols, axis=0)  # alternative way to form
        # assert this to confirm intuitive understanding of these matrices
        np.testing.assert_array_equal(latixy_transp, compare)
        latixy = np.transpose(latixy_transp)

        return longxy, latixy, cols, rows


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
