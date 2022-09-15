#!/usr/bin/env python3

"""Unit tests for buildnml
"""

import unittest
import tempfile
import shutil
import os
import logging

from pathlib import Path

from CIME.nmlgen import NamelistGenerator
from CIME.BuildTools.configure import FakeCase
from CIME.utils import expect

# pylint: disable=wrong-import-order,unused-import
from slim import add_slim_cime_py_to_path
from slim import unit_testing
from slim.slim_logging import setup_logging

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
        self._pwd = os.getcwd()
        # namelist definition file
        lnd_root = os.path.normpath(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir
            )
        )
        namelist_xml_dir = os.path.join(lnd_root, "cime_config")
        self._definition_file = [os.path.join(namelist_xml_dir, "namelist_definition_slim.xml")]
        for file_ in self._definition_file:
            expect(os.path.isfile(file_), "Namelist XML file %s not found!" % file_)

        setup_logging(logging.DEBUG)
        os.chdir(self._testdir)
        self.case = FakeCase(compiler=None, mpilib=None, debug=None)
        self.case.set_value("RUNDIR", self._testdir)
        self.case.set_value("RUN_TYPE", "startup")
        self.case.set_value("RUN_STARTDATE", "2000-01-01")
        self.case.set_value("RUN_REFCASE", "case.std")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        self.case.set_value("RUN_REFTOD", "00000")
        self.case.set_value("RUN_REFDIR", "cesm2_init")
        self.case.set_value("DIN_LOC_ROOT", ".")
        self.case.set_value("LND_DOMAIN_PATH", ".")
        self.case.set_value("LND_DOMAIN_FILE", "domain.nc")
        self.case.set_value("SLIM_SCENARIO", "global_uniform")
        self.case.set_value("SLIM_START_TYPE", "any")
        self.case.set_value("LND_GRID", "1.9x2.5")

        self.InitNML()

    def InitNML(self):
        """Re initialize the Namelist"""

        # Recreate the object so that it will start empty
        self.nmlgen = NamelistGenerator(self.case, self._definition_file)
        # ------------------------------------------------------
        # Create config dictionary
        # ------------------------------------------------------
        self.config = {}
        self.config["lnd_grid"] = self.case.get_value("LND_GRID")
        self.config["slim_scenario"] = self.case.get_value("SLIM_SCENARIO")
        self.config["slim_start_type"] = self.case.get_value("SLIM_START_TYPE")

        # ------------------------------------------------------
        # Initialize namelist defaults
        # ------------------------------------------------------
        self.nmlgen.init_defaults(infiles=[], config=self.config)

    def tearDown(self):
        """Finalize"""
        shutil.rmtree(self._testdir, ignore_errors=True)
        os.chdir(self._pwd)

    def test_check_nml_performance(self):
        """Test the check nml performance subroutine"""
        self.nmlgen.set_value("nsegspc", 20)
        check_nml_performance(self.nmlgen)
        self.nmlgen.set_value("nsegspc", -1)
        with self.assertRaisesRegex(SystemExit, "nsegspc must be positive"):
            check_nml_performance(self.nmlgen)

    def test_check_nml_general(self):
        """Test the check nml general subroutine"""
        self.case.set_value("SLIM_START_TYPE", "cold")
        self.InitNML()
        self.nmlgen.set_value("res", self.case.get_value("LND_GRID"))
        check_nml_general(self.nmlgen)
        for var in ("slim_start_type", "res"):
            val = self.nmlgen.get_value(var)
            self.nmlgen.set_value(var, None)
            with self.assertRaisesRegex(SystemExit, var + " must be set"):
                check_nml_general(self.nmlgen)

            self.nmlgen.set_value(var, val)

    def test_check_nml_data(self):
        """Test the check nml data subroutine"""
        self.case.set_value("LND_GRID", "1.9x2.5")
        self.nmlgen.set_value("res", self.case.get_value("LND_GRID"))
        # Loop through the scenarios that we have datasets for
        for scen in ("global_uniform", "realistic_from_1850", "realistic_from_2000"):
            self.case.set_value("SLIM_SCENARIO", scen)
            print("SLIM_SCENARIO = " + scen)
            self.InitNML()
            check_nml_data(self.nmlgen)
        # 1-degree has one dataset
        self.case.set_value("LND_GRID", "0.9x1.25")
        self.case.set_value("SLIM_SCENARIO", "realistic_from_2000")
        self.InitNML()
        check_nml_data(self.nmlgen)
        # Check that 1-degree for another scenario fails
        self.case.set_value("SLIM_SCENARIO", "global_uniform")
        self.InitNML()
        with self.assertRaisesRegex(SystemExit, " file is NOT set and is required"):
            check_nml_data(self.nmlgen)

        # make the dataset unset and make sure it fails
        var = "mml_surdat"
        self.nmlgen.set_value(var, "UNSET")
        with self.assertRaisesRegex(SystemExit, var + " file is NOT set and is required"):
            check_nml_data(self.nmlgen)

    def test_check_init_data(self):
        """Test the check nml initial data subroutine"""
        check_nml_initial_conditions(self.nmlgen, self.case)

        #
        # Check that normal startup acts as expected
        #
        # A cold start should have a blank finidat file
        self.case.set_value("SLIM_START_TYPE", "cold")
        self.InitNML()
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
        self.InitNML()
        self.nmlgen.set_value("finidat", "file_is_set.nc")
        Path("file_is_set.nc").touch()
        check_nml_initial_conditions(self.nmlgen, self.case)
        # Don't set the IC file, so make sure it aborts
        self.nmlgen.set_value("finidat", "UNSET")
        check_nml_initial_conditions(self.nmlgen, self.case)
        # Make sure nrevsn can't be set
        self.nmlgen.set_value("nrevsn", "file_is_set.nc")
        with self.assertRaisesRegex(
            SystemExit, "nrevsn can NOT be set except when RUN_TYPE is a branch"
        ):
            check_nml_initial_conditions(self.nmlgen, self.case)

    def test_check_init_data_branch(self):
        """Test the check nml initial data subroutine for branch"""
        #
        # Check that a branch works correctly
        #
        self.case.set_value("RUN_TYPE", "branch")
        self.case.set_value("SLIM_START_TYPE", "required")
        self.case.set_value("RUN_REFCASE", "TESTCASE")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        nrevsn = "TESTCASE.slim.r.0001-01-01-00000.nc"
        Path(nrevsn).touch()
        self.InitNML()
        check_nml_initial_conditions(self.nmlgen, self.case)
        expect(
            self.nmlgen.get_value("nrevsn"),
            nrevsn,
            "finidat should be set correctly for hybrid case",
        )
        # Make sure finidat can't be set
        self.nmlgen.set_value("finidat", "file_is_set.nc")
        with self.assertRaisesRegex(SystemExit, "finidat can NOT be set when RUN_TYPE is a branch"):
            check_nml_initial_conditions(self.nmlgen, self.case)

    def test_check_init_data_hybrid(self):
        """Test the check nml initial data subroutine for hybrid"""
        #
        # Check that a hybrid works correctly
        #
        self.case.set_value("RUN_TYPE", "hybrid")
        self.case.set_value("SLIM_START_TYPE", "required")
        self.case.set_value("RUN_REFCASE", "TESTCASE")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        finidat = "TESTCASE.slim.r.0001-01-01-00000.nc"
        Path(finidat).touch()
        self.InitNML()
        check_nml_initial_conditions(self.nmlgen, self.case)
        expect(
            self.nmlgen.get_value("finidat"),
            finidat,
            "finidat should be set correctly for hybrid case",
        )


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
