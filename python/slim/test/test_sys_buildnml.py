#!/usr/bin/env python3

"""System tests for buildnml
"""

import unittest
import tempfile
import shutil
import os
import re
import logging

from pathlib import Path

from CIME.BuildTools.configure import FakeCase
from CIME.utils import expect

# pylint: disable=wrong-import-order,unused-import
from slim import add_slim_cime_py_to_path
from slim import unit_testing

from slim_cime_py.buildnml import buildnml

logger = logging.getLogger(__name__)

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


def getVariableFromNML(nmlfile, variable):
    """Get a variable from the namelist file"""
    with open(nmlfile, "r") as nfile:
        for line in nfile:
            if variable in line:
                print("lnd_in:" + line)
                match = re.search('= ["]*([ a-zA-Z0-9._//-]+)["]*', line)
                if match is not None:
                    return match.group(1)
                match = re.search("= [']*([ a-zA-Z0-9._//-]+)[']*", line)
                if match is not None:
                    return match.group(1)
    return None


def addLinesToUserNL(user_nl_file, lines):
    """Add some lines to the user_nl_file for SLIM"""
    if os.path.exists(user_nl_file):
        os.remove(user_nl_file)
    print("Add lines to " + user_nl_file)
    with open(user_nl_file, "x") as userfile:
        for line in lines:
            userfile.write(line)
            print(line)
    userfile.close()
    print("Close file")


class TestBuildNML(unittest.TestCase):
    """System Tests of buildnml"""

    def setUp(self):
        """Initialize"""
        self._testdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self._testdir)
        # namelist definition file
        lnd_root = os.path.normpath(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir
            )
        )
        self.case = FakeCase(compiler=None, mpilib=None, debug=None)
        self.case.set_value("CASEROOT", self._testdir)
        self.case.set_value("RUN_TYPE", "any")
        self.case.set_value("RUN_STARTDATE", "2000-01-01")
        self.case.set_value("RUN_REFCASE", "case.std")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        self.case.set_value("RUN_REFTOD", "00000")
        self.case.set_value("RUN_REFDIR", "cesm2_init")
        self.case.set_value("RUNDIR", ".")
        self.case.set_value("NINST_LND", 1)
        self.case.set_value("NCPL_BASE_PERIOD", "day")
        self.case.set_value("LND_NCPL", 48)
        self.case.set_value("SLIM_FORCE_COLDSTART", "off")
        self.case.set_value("SLIM_SCENARIO", "global_uniform")
        self.case.set_value("COMP_ROOT_DIR_LND", lnd_root)
        self.case.set_value("DIN_LOC_ROOT", ".")
        self.case.set_value("LND_DOMAIN_PATH", ".")
        self.case.set_value("LND_DOMAIN_FILE", "domain.nc")
        self.case.set_value("SLIM_START_TYPE", "cold")
        self.case.set_value("LND_GRID", "1.9x2.5")

    def tearDown(self):
        """Finalize"""
        os.chdir(self.curdir)
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_simple(self):
        """Test a simple call of buildnml"""
        for scenario in ("global_uniform", "realistic_from_1850", "realistic_from_2000"):
            self.case.set_value("SLIM_SCENARIO", scenario)
            buildnml(self.case, self._testdir, "slim")
            expect(
                os.path.isfile("Buildconf/slimconf/lnd_in"),
                "Namelist file lnd_in should exist in Buildconf after running buildnml",
            )
            expect(
                os.path.isfile("lnd_in"), "Namelist file lnd_in should exist after running buildnml"
            )
            expect(
                os.path.isfile("Buildconf/slim.input_data_list"),
                "Input data list file should exist after running buildnml",
            )

    def test_hybrid_start(self):
        """Test a hybrid startup call of buildnml"""
        self.case.set_value("SLIM_START_TYPE", "required")
        self.case.set_value("RUN_TYPE", "hybrid")
        self.case.set_value("RUN_REFCASE", "TESTCASE")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        self.case.set_value("RUN_REFTOD", "00000")
        Path("TESTCASE.slim.r.0001-01-01-00000.nc").touch()
        buildnml(self.case, self._testdir, "slim")
        expect(
            os.path.isfile("Buildconf/slimconf/lnd_in"),
            "Namelist file lnd_in should exist in Buildconf after running buildnml",
        )
        expect(os.path.isfile("lnd_in"), "Namelist file lnd_in should exist after running buildnml")
        expect(
            os.path.isfile("Buildconf/slim.input_data_list"),
            "Input data list file should exist after running buildnml",
        )
        value = getVariableFromNML("lnd_in", "finidat")
        self.assertEqual(
            value, "TESTCASE.slim.r.0001-01-01-00000.nc", msg="finidat not set as expected"
        )

    def test_branch_start(self):
        """Test a branch startup call of buildnml"""
        self.case.set_value("SLIM_START_TYPE", "required")
        self.case.set_value("RUN_TYPE", "branch")
        self.case.set_value("RUN_REFCASE", "TESTCASE")
        self.case.set_value("RUN_REFDATE", "0001-01-01")
        self.case.set_value("RUN_REFTOD", "00000")
        Path("TESTCASE.slim.r.0001-01-01-00000.nc").touch()
        buildnml(self.case, self._testdir, "slim")
        expect(
            os.path.isfile("Buildconf/slimconf/lnd_in"),
            "Namelist file lnd_in should exist in Buildconf after running buildnml",
        )
        expect(os.path.isfile("lnd_in"), "Namelist file lnd_in should exist after running buildnml")
        expect(
            os.path.isfile("Buildconf/slim.input_data_list"),
            "Input data list file should exist after running buildnml",
        )
        value = getVariableFromNML("lnd_in", "nrevsn")
        self.assertEqual(
            value, "TESTCASE.slim.r.0001-01-01-00000.nc", msg="nrevsn not set as expected"
        )

    def test_start_types(self):
        """Test start types of buildnml"""
        # Cold start types
        finidat = " "
        lines = []
        for stype in ("cold", "any", "required"):
            print("Type: " + stype)
            addLinesToUserNL("user_nl_slim", lines)

            self.case.set_value("SLIM_START_TYPE", stype)
            buildnml(self.case, self._testdir, "slim")
            expect(
                os.path.isfile("Buildconf/slimconf/lnd_in"),
                "Namelist file lnd_in should exist in Buildconf after running buildnml",
            )
            expect(
                os.path.isfile("lnd_in"), "Namelist file lnd_in should exist after running buildnml"
            )
            expect(
                os.path.isfile("Buildconf/slim.input_data_list"),
                "Input data list file should exist after running buildnml",
            )
            value = getVariableFromNML("lnd_in", "finidat")
            self.assertEqual(value, finidat, msg="finidat not set as expected")
        stype = "required"
        finidat = "TESTFINIDATFILENAME.nc"
        Path(finidat).touch()
        lines = ["finidat = '" + finidat + "'\n"]
        print("Type: " + stype)
        addLinesToUserNL("user_nl_slim", lines)

        self.case.set_value("SLIM_START_TYPE", stype)
        buildnml(self.case, self._testdir, "slim")
        expect(
            os.path.isfile("Buildconf/slimconf/lnd_in"),
            "Namelist file lnd_in should exist in Buildconf after running buildnml",
        )
        expect(os.path.isfile("lnd_in"), "Namelist file lnd_in should exist after running buildnml")
        expect(
            os.path.isfile("Buildconf/slim.input_data_list"),
            "Input data list file should exist after running buildnml",
        )
        value = getVariableFromNML("lnd_in", "finidat")
        self.assertEqual(value, finidat, msg="finidat not set as expected")
        os.remove("user_nl_slim")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
