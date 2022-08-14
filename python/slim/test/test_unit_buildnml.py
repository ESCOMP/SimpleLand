#!/usr/bin/env python3

"""Unit tests for buildnml
"""

import unittest
import tempfile
import shutil
import os

from unittest import mock
from slim import unit_testing
from slim import add_slim_cime_py_to_path
from slim_cime_py import buildnml

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestPathUtils(unittest.TestCase):
    """Tests of buildnml"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
