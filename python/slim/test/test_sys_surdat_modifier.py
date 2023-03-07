#!/usr/bin/env python3

"""System tests for surdat_modifier

"""

import os
import re

import unittest
import tempfile
import shutil

import xarray as xr

from slim.path_utils import path_to_slim_root
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
        - /testinputs directory and surdat_in, located in /testinputs
        Make /_tempdir for use by these tests.
        Obtain path and names for the files being created in /_tempdir:
        - modify_surdat.cfg
        - surdat_out.nc
        """
        self._cfg_template_path = os.path.join(
            path_to_slim_root(), "tools/modify_input_files/modify_surdat_template.cfg"
        )
        # TODO Generate simple surdat_in (see test_unit_modify_surdat.py)
        # rather than pointing to a /testinputs directory?
        # If so, then need to output this to _tempdir and read back in
        # as _surdat_in. For call to write_output see surdat_modifier.py.
        # Now see next TODO.
        testinputs_path = os.path.join(path_to_slim_root(), "python/slim/test/testinputs")
        self._surdat_in = os.path.join(
            testinputs_path,
            "slim_realistic_fromCLM5_alb2000_hc2000_rs2000annual_f19_20190621.nc",
        )
        self._tempdir = tempfile.mkdtemp()
        self._cfg_file_path = os.path.join(self._tempdir, "modify_surdat.cfg")
        self._surdat_out = os.path.join(self._tempdir, "surdat_out.nc")

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

        surdat_in_data = xr.open_dataset(self._surdat_in)
        surdat_out_data = xr.open_dataset(self._surdat_out)
        # assert that surdat_out equals surdat_in
        self.assertTrue(surdat_out_data.equals(surdat_in_data))

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
        surdat_in_data = xr.open_dataset(self._surdat_in)
        surdat_out_data = xr.open_dataset(self._surdat_out)
        # assert that surdat_out does not equal surdat_in
        self.assertFalse(surdat_out_data.equals(surdat_in_data))

        # compare surdat_out to surdat_out_baseline located in /testinputs
        # TODO Generate surdat_out_baseline as was done above for _surdat_in.
        # No need to output surdat_out_baseline; just compare.
        surdat_out_baseline = self._surdat_in[:-3] + "_modified" + self._surdat_in[-3:]
        surdat_out_base_data = xr.open_dataset(surdat_out_baseline)
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
                        line = "idealized = True"
                    elif re.match(r" *lnd_lat_1 *=", line):
                        line = "lnd_lat_1 = -10\n"
                    elif re.match(r" *lnd_lat_2 *=", line):
                        line = "lnd_lat_2 = 45\n"
                    elif re.match(r" *lnd_lon_1 *=", line):
                        line = "lnd_lon_1 = 10\n"
                    elif re.match(r" *lnd_lon_2 *=", line):
                        line = "lnd_lon_2 = 341\n"
                    elif re.match(r" *glc_mask *=", line):
                        line = "glc_mask = 1 1 1 1 1 1 1 0 1 1 1 1\n"
                    elif re.match(r" *alb_gvd *=", line):
                        line = "alb_gvd = 0 0.1 0.2 0.3 0.4 0.5 0.5 0.4 0.3 0.2 0.1 0\n"
                    elif re.match(r" *alb_svd *=", line):
                        line = "alb_svd = 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n"
                    elif re.match(r" *bucketdepth *=", line):
                        line = "bucketdepth = 195 205 195 195 195 195 195 195 195 195 195 195\n"
                    cfg_out.write(line)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
