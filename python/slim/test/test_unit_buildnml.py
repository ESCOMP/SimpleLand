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
    # pylint: disable=too-many-public-methods
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
        self.case = FakeCase(compiler=None, comp_interface="nuopc", mpilib=None, debug=None)
        self.case.set_value("RUNDIR", self._testdir)
        self.case.set_value("RUN_TYPE", "startup")
        self.case.set_value("RUN_STARTDATE", "2000-01-01")
        self.case.set_value("RUN_REFCASE", "case.std")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        self.case.set_value("RUN_REFTOD", "00000")
        self.case.set_value("RUN_REFDIR", "cesm2_init")
        self.case.set_value("DIN_LOC_ROOT", ".")
        self.case.set_value("CALENDAR", "NO_LEAP")
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

    def test_check_nml_history(self):
        """Test the check nml history subroutine"""
        self.nmlgen.set_value("hist_mfilt", [20])
        check_nml_history(self.nmlgen)
        self.nmlgen.set_value("hist_mfilt", [0])
        with self.assertRaisesRegex(SystemExit, "hist_mfilt must be 1 or larger"):
            check_nml_history(self.nmlgen)

        self.InitNML()
        # Make a list of settings to a list of history tape streams
        self.nmlgen.set_value("use_noio", ".false.")
        self.nmlgen.set_value("hist_empty_htapes", ".true.")
        self.nmlgen.set_value("hist_mfilt", [1, 1, 2, 3, 4, 5])
        self.nmlgen.set_value("hist_ndens", [1, 1, 2, 1, 1, 1])
        self.nmlgen.set_value("hist_nhtfrq", [1, 1, -24, 2, 2, 2])
        self.nmlgen.set_value("hist_avgflag_pertape", ["A", "I", "X", "M", "A", "A"])
        self.nmlgen.set_value("hist_fincl1", ["A:I", "B:A", "C:X", "D:M", "E:A", "F:A"])
        self.nmlgen.set_value("hist_fincl2", ["A", "B", "C", "D", "E", "F"])
        self.nmlgen.set_value("hist_fincl3", ["A", "B", "C", "D", "E", "F"])
        self.nmlgen.set_value("hist_fincl4", ["A", "B", "C", "D", "E", "F"])
        self.nmlgen.set_value("hist_fincl5", ["A", "B", "C", "D", "E", "F"])
        self.nmlgen.set_value("hist_fincl6", ["A", "B", "C", "D", "E", "F"])
        check_nml_history(self.nmlgen)
        # and another complex case
        self.InitNML()
        self.nmlgen.set_value("hist_empty_htapes", ".false.")
        self.nmlgen.set_value("hist_fexcl1", ["A", "B", "C", "D", "E", "F"])
        check_nml_history(self.nmlgen)
        # Check that use_noio works if you don't set any hist_* options
        self.InitNML()
        self.nmlgen.set_value("use_noio", ".true.")
        check_nml_history(self.nmlgen)

    def test_check_nml_history_simple_fails_bad_timeavg(self):
        """Test the check nml history subroutine for simple fails bad time avg"""
        self.nmlgen.set_value("hist_fincl1", ["A:Z"])
        with self.assertRaisesRegex(
            SystemExit, "History averaging option Z is not valid in hist_fincl1"
        ):
            check_nml_history(self.nmlgen)

    def test_check_nml_history_simple_fails_bad_characters(self):
        """Test the check nml history subroutine for simple fails bad characters in field"""
        self.nmlgen.set_value("hist_fincl1", ["A%$#@!~"])
        with self.assertRaisesRegex(
            SystemExit,
            r"History field name hist_fincl1 has invalid characters or whitespace in it=",
        ):
            check_nml_history(self.nmlgen)

    def test_check_nml_history_simple_fails_white_space_in_field(self):
        """Test the check nml history subroutine for simple fails bad characters in field"""
        self.nmlgen.set_value("hist_fincl1", [" "])
        with self.assertRaisesRegex(
            SystemExit,
            r"History field name hist_fincl1 has invalid characters or whitespace in it=",
        ):
            check_nml_history(self.nmlgen)

    def test_check_nml_history_complex_fails_array_size_not_consistent(self):
        """Test the check nml history subroutine for complex fails array size not consistent"""
        self.nmlgen.set_value("hist_mfilt", [1])
        self.nmlgen.set_value("hist_ndens", [1, 1])
        self.nmlgen.set_value("hist_fincl2", ["A"])
        with self.assertRaisesRegex(
            SystemExit, r"hist_mfilt array size does not agree with the expected size of 2"
        ):
            check_nml_history(self.nmlgen)

    def test_check_nml_history_complex_fails_array_size(self):
        """Test the check nml history subroutine for complex fails array size"""
        self.nmlgen.set_value("hist_mfilt", [1, 1, 1, 1, 1, 1, 1])
        self.nmlgen.set_value("hist_fincl2", ["A"])
        self.nmlgen.set_value("hist_fincl3", ["A"])
        self.nmlgen.set_value("hist_fincl4", ["A"])
        self.nmlgen.set_value("hist_fincl5", ["A"])
        self.nmlgen.set_value("hist_fincl6", ["A"])
        with self.assertRaisesRegex(
            SystemExit, r"hist_mfilt array size does not agree with the expected size of 6"
        ):
            check_nml_history(self.nmlgen)

    def test_check_nml_history_complex_fails_excl_with_empty(self):
        """Test the check nml history subroutine for complex fails exclude used with
        empty history
        """
        self.nmlgen.set_value("hist_empty_htapes", ".true.")
        self.nmlgen.set_value("hist_fexcl1", ["A"])
        with self.assertRaisesRegex(
            SystemExit, "hist_fexcl1 can not be set if hist_empty_htapes is set to true"
        ):
            check_nml_history(self.nmlgen)

    def test_check_nml_history_fails_excl_use_noio(self):
        """Test the check nml history subroutine for fails use_noio used with other
        history settings
        """
        self.nmlgen.set_value("use_noio", ".true.")
        self.nmlgen.set_value("hist_fincl1", ["A"])
        with self.assertRaisesRegex(
            SystemExit,
            r"use_noio turns off all history output"
            + ", so no hist_ namelist option should also be set.*",
        ):
            check_nml_history(self.nmlgen)

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
            check_nml_data(self.nmlgen, self.case)
        # 1-degree has one dataset
        self.case.set_value("LND_GRID", "0.9x1.25")
        self.case.set_value("SLIM_SCENARIO", "realistic_from_2000")
        self.InitNML()
        check_nml_data(self.nmlgen, self.case)
        # Check that 1-degree for another scenario fails
        self.case.set_value("SLIM_SCENARIO", "global_uniform")
        self.InitNML()
        with self.assertRaisesRegex(SystemExit, " file is NOT set and is required"):
            check_nml_data(self.nmlgen, self.case)

        # make the dataset unset and make sure it fails
        var = "mml_surdat"
        self.nmlgen.set_value(var, "UNSET")
        with self.assertRaisesRegex(SystemExit, var + " file is NOT set and is required"):
            check_nml_data(self.nmlgen, self.case)

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

    def test_check_set_user_defined(self):
        """Test the check nml initial data subroutine for user_defined """
        self.case.set_value("SLIM_SCENARIO", "user_defined")
        self.InitNML()
        with self.assertRaisesRegex(
            SystemExit, "When SLIM_SCENARIO is set to user_defined, you must provide the mml_surdat"
        ):
            check_nml_data(self.nmlgen, self.case)


    def test_check_use_init_interp(self):
        """Test the check nml initial data subroutine for use_init_interp options"""
        self.case.set_value("SLIM_START_TYPE", "startup")
        finidat = "finidat_file_to_interpolate_from.nc"
        self.nmlgen.set_value("finidat", finidat)
        Path(finidat).touch()
        self.nmlgen.set_value("use_init_interp", ".true.")
        finidat_dest = "finidat_file_to_create.nc"
        self.nmlgen.set_value("finidat_interp_dest", finidat_dest)
        check_nml_initial_conditions(self.nmlgen, self.case)
        check_nml_data(self.nmlgen, self.case)

    def test_check_use_init_interp_fails_cold(self):
        """Test the check nml initial data subroutine for use_init_interp options that fail 1"""
        # Cold start with use_init_interp should fail
        self.case.set_value("SLIM_START_TYPE", "cold")
        self.nmlgen.set_value("use_init_interp", ".true.")
        check_nml_initial_conditions(self.nmlgen, self.case)
        with self.assertRaisesRegex(
            SystemExit, "use_init_interp can not be set to TRUE for a cold start"
        ):
            check_nml_data(self.nmlgen, self.case)

    def test_check_use_init_interp_fails_setdest(self):
        """Test the check nml initial data subroutine for use_init_interp options that fail 2"""
        # Setting finidat_interp_dest without use_init_interp should fail
        self.nmlgen.set_value("use_init_interp", ".false.")
        finidat_dest = "finidat_file_to_create.nc"
        self.nmlgen.set_value("finidat_interp_dest", finidat_dest)
        check_nml_initial_conditions(self.nmlgen, self.case)
        with self.assertRaisesRegex(
            SystemExit, "finidat_interp_dest can NOT be set if use_init_interp is not on"
        ):
            check_nml_data(self.nmlgen, self.case)

    def test_check_use_init_interp_fails_branch(self):
        """Test the check nml initial data subroutine for use_init_interp options that fail 3"""
        # Branch type should fail for use_init_interp set
        self.case.set_value("RUN_TYPE", "branch")
        self.case.set_value("SLIM_START_TYPE", "required")
        self.case.set_value("RUN_REFCASE", "TESTCASE")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        nrevsn = "TESTCASE.slim.r.0001-01-01-00000.nc"
        self.nmlgen.set_value("use_init_interp", ".true.")
        Path(nrevsn).touch()
        check_nml_initial_conditions(self.nmlgen, self.case)
        with self.assertRaisesRegex(
            SystemExit, "use_init_interp can NOT be set to TRUE for a branch run type"
        ):
            check_nml_data(self.nmlgen, self.case)

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

    def test_check_dtime(self):
        """Test the check nml dtime"""
        self.case.set_value("NCPL_BASE_PERIOD", "hour")
        self.case.set_value("CALENDAR", "GREGORIAN")
        self.case.set_value("LND_NCPL", 1)
        self.InitNML()
        check_nml_dtime(self.nmlgen, self.case)
        expect(
            self.nmlgen.get_value("dtime"),
            3600,
            "dtime should be 3600 seconds for an hour",
        )
        self.case.set_value("NCPL_BASE_PERIOD", "day")
        self.case.set_value("CALENDAR", "NO_LEAP")
        self.case.set_value("LND_NCPL", 48)
        self.InitNML()
        check_nml_dtime(self.nmlgen, self.case)
        expect(
            self.nmlgen.get_value("dtime"),
            1800,
            "dtime should be 1800 seconds for 48 cycles per day",
        )
        self.case.set_value("NCPL_BASE_PERIOD", "year")
        self.case.set_value("CALENDAR", "NO_LEAP")
        self.case.set_value("LND_NCPL", 365)
        self.InitNML()
        check_nml_dtime(self.nmlgen, self.case)
        expect(
            self.nmlgen.get_value("dtime"),
            86400,
            "dtime should be 86400 seconds for 365 cycles per year",
        )
        self.case.set_value("NCPL_BASE_PERIOD", "year")
        self.case.set_value("CALENDAR", "NO_LEAP")
        self.case.set_value("LND_NCPL", 36500)
        self.InitNML()
        check_nml_dtime(self.nmlgen, self.case)
        expect(
            self.nmlgen.get_value("dtime"),
            8640,
            "dtime should be 8640 seconds for 36500 cycles per decade",
        )

    def test_check_dtime_fail_invalid_cal_year(self):
        """Test the check nml dtime fail test for invalid calendar year"""
        self.case.set_value("NCPL_BASE_PERIOD", "year")
        self.case.set_value("CALENDAR", "GREGORIAN")
        self.case.set_value("LND_NCPL", 1)
        with self.assertRaisesRegex(
            SystemExit, "ERROR: Invalid CALENDAR for NCPL_BASE_PERIOD year"
        ):
            check_nml_dtime(self.nmlgen, self.case)

    def test_check_dtime_fail_invalid_cal_decade(self):
        """Test the check nml dtime fail test for invalid calendar decade"""
        self.case.set_value("NCPL_BASE_PERIOD", "decade")
        self.case.set_value("CALENDAR", "GREGORIAN")
        self.case.set_value("LND_NCPL", 1)
        with self.assertRaisesRegex(
            SystemExit, "ERROR: Invalid CALENDAR for NCPL_BASE_PERIOD decade"
        ):
            check_nml_dtime(self.nmlgen, self.case)

    def test_check_dtime_fail_invalid_base_period(self):
        """Test the check nml dtime fail test for invalid base period"""
        self.case.set_value("NCPL_BASE_PERIOD", "minute")
        self.case.set_value("CALENDAR", "GREGORIAN")
        self.case.set_value("LND_NCPL", 1)
        with self.assertRaisesRegex(SystemExit, "ERROR: Invalid NCPL_BASE_PERIOD "):
            check_nml_dtime(self.nmlgen, self.case)

    def test_check_dtime_fail_invalid_division(self):
        """Test the check nml dtime fail test for invalid coupling division"""
        self.case.set_value("NCPL_BASE_PERIOD", "day")
        self.case.set_value("CALENDAR", "GREGORIAN")
        self.case.set_value("LND_NCPL", 47)
        with self.assertRaisesRegex(
            SystemExit, "ERROR: LND_NCPL=47 doesn't divide evenly into NCPL_BASE_PERIOD day"
        ):
            check_nml_dtime(self.nmlgen, self.case)

    def test_check_dtime_fail_too_short(self):
        """Test the check nml dtime fail test for too short"""
        self.case.set_value("NCPL_BASE_PERIOD", "hour")
        self.case.set_value("LND_NCPL", 3600)
        with self.assertRaisesRegex(
            SystemExit,
            "ERROR: LND_NCPL=3600 is too frequent which gives a time step that is too short",
        ):
            check_nml_dtime(self.nmlgen, self.case)

    def test_check_dtime_fail_too_long(self):
        """Test the check nml dtime fail test for too long"""
        self.case.set_value("NCPL_BASE_PERIOD", "year")
        self.case.set_value("CALENDAR", "NO_LEAP")
        self.case.set_value("LND_NCPL", 5)
        with self.assertRaisesRegex(
            SystemExit,
            "ERROR: LND_NCPL=5 is too infrequent which gives a time step that is too long",
        ):
            check_nml_dtime(self.nmlgen, self.case)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
