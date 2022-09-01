#!/usr/bin/env python3

"""System tests for buildnml
"""

import unittest
import tempfile
import shutil
import os
import logging

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


class TestBuildNML(unittest.TestCase):
    """System Tests of buildnml"""

    def setUp(self):
        """Initialize"""
        self._testdir = tempfile.mkdtemp()
        # namelist definition file
        lnd_root = os.path.normpath(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir
            )
        )
        self.case = FakeCase(compiler=None, mpilib=None, debug=None)
        self.case.set_value("CASEROOT", self._testdir)
        self.case.set_value("RUN_TYPE", "startup")
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
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_simple(self):
        """Test a simple call of buildnml"""
        buildnml(self.case, self._testdir, "slim")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
