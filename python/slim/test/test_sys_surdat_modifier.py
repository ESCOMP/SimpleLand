#!/usr/bin/env python3

"""System tests for surdat_modifier

"""

import os
import re

import unittest
import tempfile

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
        Obtain path to the existing:
        - modify_surdat_template.cfg file
#       - /testinputs directory and surdat_in, located in /testinputs
        Make /_tempdir for use by these tests.
        Obtain path and names for the files being created in /_tempdir:
        - modify_surdat.cfg
        - surdat_out.nc
        """
        self._cfg_template_path = os.path.join(
            path_to_slim_root(), "tools/modify_input_files/modify_surdat_template.cfg"
        )
        self._tempdir = tempfile.mkdtemp()
        print(self._tempdir)
        self._cfg_file_path = os.path.join(self._tempdir, "modify_surdat.cfg")
        self._surdat_out = os.path.join(self._tempdir, "surdat_out.nc")
        self._surdat_in = os.path.join(self._tempdir, "surdat_in.nc")
        self._months = 12

        # -----------------------------------------------------------
        # create dummy SLIM surdat file
        # TODO Functionalize?

        # get lon/lat that would normally come from a surdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = 2  # expects min_lon < max_lon
        min_lat = 3  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=10, _min_lat=min_lat, _max_lat=12
        )
        # create xarray dataset containing lev1 variables;
        # the surdat_modify tool reads variables like this from a surdat file
        var_1d = np.arange(cols)
        var_lev1 = var_1d * np.ones((self._months, rows, cols))
        self._surdat_in_data = xr.Dataset(
            data_vars=dict(
                time=(["time"], np.arange(self._months) + 1),
                lsmlon=(["lsmlon"], longxy[0,:]),
                lsmlat=(["lsmlat"], latixy[:,0]),
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
        # save in tempdir; _in and _out files are the same file in this case
        write_output(self._surdat_in_data, self._surdat_in, self._surdat_in, 'surdat')
        # -----------------------------------------------------------

    def tearDown(self):
        """
        Remove temporary directory
        """
#       shutil.rmtree(self._tempdir, ignore_errors=True)

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

        # compare surdat_out to surdat_out_baseline
        # -----------------------------------------------------------
        # generate surdat_out_baseline as was done above for surdat_in
        # no need to output surdat_out_baseline; just compare
        # TODO Functionalize?

        # get lon/lat that would normally come from a surdat file
        # self._get_longxy_latixy will convert -180 to 180 to 0-360 longitudes
        # get cols, rows also
        min_lon = 2  # expects min_lon < max_lon
        min_lat = 3  # expects min_lat < max_lat
        longxy, latixy, cols, rows = self._get_longxy_latixy(
            _min_lon=min_lon, _max_lon=10, _min_lat=min_lat, _max_lat=12
        )
        # create xarray dataset containing lev1 variables;
        # the surdat_modify tool reads variables like this from a surdat file
        var_1d = np.arange(cols)
        var_lev1 = var_1d * np.ones((self._months, rows, cols))
        # the next 4 lines must agree with the modifications in
        # def _create_config_file_complete
        modified_1 = np.ones((self._months, rows, cols))
        modified_2 = 0 * np.ones((self._months, rows, cols))
        modified_3 = 0.5 * np.ones((self._months, rows, cols))
        modified_4 = 195 * np.ones((self._months, rows, cols))
        surdat_out_base_data = xr.Dataset(
            data_vars=dict(
                time=(["time"], np.arange(self._months) + 1),
                lsmlon=(["lsmlon"], longxy[0,:]),
                lsmlat=(["lsmlat"], latixy[:,0]),
                glc_mask=(["time", "lsmlat", "lsmlon"], modified_1),
                alb_gvd=(["time", "lsmlat", "lsmlon"], modified_2),
                alb_svd=(["time", "lsmlat", "lsmlon"], modified_3),
                alb_gnd=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_snd=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_gvf=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_svf=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_gnf=(["time", "lsmlat", "lsmlon"], var_lev1),
                alb_snf=(["time", "lsmlat", "lsmlon"], var_lev1),
                bucketdepth=(["time", "lsmlat", "lsmlon"], modified_4),
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
        # -----------------------------------------------------------

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
                    elif re.match(r" *idealized *=", line):
                        line = "idealized = False"
                    elif re.match(r" *lnd_lat_1 *=", line):
                        line = "lnd_lat_1 = 3\n"
                    elif re.match(r" *lnd_lat_2 *=", line):
                        line = "lnd_lat_2 = 12\n"
                    elif re.match(r" *lnd_lon_1 *=", line):
                        line = "lnd_lon_1 = 2\n"
                    elif re.match(r" *lnd_lon_2 *=", line):
                        line = "lnd_lon_2 = 10\n"
                    elif re.match(r" *glc_mask *=", line):
                        line = "glc_mask = 1 1 1 1 1 1 1 1 1 1 1 1\n"
                    elif re.match(r" *alb_gvd *=", line):
                        line = "alb_gvd = 0 0 0 0 0 0 0 0 0 0 0 0\n"
                    elif re.match(r" *alb_svd *=", line):
                        line = "alb_svd = 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n"
                    elif re.match(r" *bucketdepth *=", line):
                        line = "bucketdepth = 195 195 195 195 195 195 195 195 195 195 195 195\n"
                    cfg_out.write(line)

    def _get_longxy_latixy(self, _min_lon, _max_lon, _min_lat, _max_lat):
        """
        Return longxy, latixy, cols, rows
        TODO: This function is copied from test_unit_modify_surdat.py
              May wish to separate into a separate file
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
