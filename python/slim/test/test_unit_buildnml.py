#!/usr/bin/env python3

"""Unit tests for buildnml
"""

import unittest
import tempfile
import shutil
import os
import logging

from CIME.nmlgen import NamelistGenerator
from CIME.BuildTools.configure import FakeCase
from CIME.utils import expect

# pylint: disable=wrong-import-order,unused-import
from slim import add_slim_cime_py_to_path
from slim import unit_testing

from slim_cime_py.buildnml import (
    check_nml_dtime,
    check_nml_performance,
    check_nml_general,
    check_nml_history,
    check_nml_data,
    check_nml_initial_conditions,
)

logger = logging.getLogger(__name__)

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestPathUtils(unittest.TestCase):
    """Tests of buildnml"""

    def setUp(self):
        """Initialize"""
        self._testdir = tempfile.mkdtemp()
        # namelist definition file
        lnd_root = os.path.normpath(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir
            )
        )
        namelist_xml_dir = os.path.join(lnd_root, "cime_config")
        definition_file = [os.path.join(namelist_xml_dir, "namelist_definition_slim.xml")]
        for file_ in definition_file:
            expect(os.path.isfile(file_), "Namelist XML file %s not found!" % file_)

        self.case = FakeCase(compiler=None, mpilib=None, debug=None)
        self.case.set_value("DIN_LOC_ROOT", ".")
        self.case.set_value("LND_DOMAIN_PATH", ".")
        self.case.set_value("LND_DOMAIN_FILE", "domain.nc")
        self.case.set_value("SLIM_START_TYPE", "cold")
        self.case.set_value("LND_GRID", "0.9x1.25")
        # Create the namelist generator object - independent of instance
        self.nmlgen = NamelistGenerator(self.case, definition_file)
        # ------------------------------------------------------
        # Create config dictionary
        # ------------------------------------------------------
        self.config = {}
        self.config["lnd_grid"] = self.case.get_value("LND_GRID")
        self.config["slim_scenario"] = "1850_control"

        # ------------------------------------------------------
        # Initialize namelist defaults
        # ------------------------------------------------------
        self.nmlgen.init_defaults(infiles=[], config=self.config)

    def tearDown(self):
        """Finalize"""
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_check_nml_performance(self):
        """Test the check nml performance subroutine"""
        self.nmlgen.set_value("nsegspc", 20)
        check_nml_performance(self.nmlgen)
        self.nmlgen.set_value("nsegspc", -1)
        with self.assertRaisesRegex(SystemExit, "nsegspc must be positive"):
            check_nml_performance(self.nmlgen)

    def test_check_nml_general(self):
        """Test the check nml general subroutine"""
        self.nmlgen.set_value("res", self.case.get_value("LND_GRID"))
        self.nmlgen.set_value("slim_start_type", "cold")
        check_nml_general(self.nmlgen)
        for var in ("slim_start_type", "res"):
            val = self.nmlgen.get_value(var)
            self.nmlgen.set_value(var, None)
            with self.assertRaisesRegex(SystemExit, var + " must be set"):
                check_nml_general(self.nmlgen)

            self.nmlgen.set_value(var, val)

    def test_check_nml_data(self):
        """Test the check nml data subroutine"""
        self.nmlgen.set_value("res", self.case.get_value("LND_GRID"))
        # Loop through the scenarios that we have datasets for
        self.config["lnd_grid"] = "1.9x2.5"
        for scen in ("global_uniform", "realistic_from_1850", "realistic_from_2000"):
            self.config["slim_scenario"] = scen
            self.nmlgen.init_defaults(infiles=[], config=self.config)
            check_nml_data(self.nmlgen)

        # make the dataset unset and make sure it fails
        var = "mml_surdat"
        self.nmlgen.set_value(var, "UNSET")
        with self.assertRaisesRegex(SystemExit, var + " file is NOT set and is required"):
            check_nml_data(self.nmlgen)

    def test_check_init_data(self):
        """Test the check nml initial data subroutine"""
        check_nml_initial_conditions(self.nmlgen, self.case)
        # A cold start should have a blank finidat file
        self.case.set_value("SLIM_START_TYPE", "cold")
        self.nmlgen.init_defaults(infiles=[], config=self.config)
        check_nml_initial_conditions(self.nmlgen, self.case)
        val = self.nmlgen.get_value("finidat")
        expect("' '", val)
        # If you have a cold-start and explicitly set the finidat file that should error out
        self.nmlgen.set_value("finidat", "file_is_set.nc")
        with self.assertRaisesRegex(
            SystemExit, "finidat is set but SLIM_START_TYPE is cold which is a contradiction"
        ):
            check_nml_initial_conditions(self.nmlgen, self.case)
        # Set to startup so that finidat is required
        self.case.set_value("SLIM_START_TYPE", "startup")
        self.nmlgen.init_defaults(infiles=[], config=self.config)
        self.nmlgen.set_value("finidat", "file_is_set.nc")
        check_nml_initial_conditions(self.nmlgen, self.case)
        # Don't set the IC file, so make sure it aborts
        self.nmlgen.set_value("finidat", "UNSET")
        check_nml_initial_conditions(self.nmlgen, self.case)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
