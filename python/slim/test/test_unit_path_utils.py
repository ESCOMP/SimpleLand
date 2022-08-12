#!/usr/bin/env python3

"""Unit tests for path_utils
"""

import unittest
import tempfile
import shutil
import os

from unittest import mock
from slim import unit_testing
from slim import path_utils

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestPathUtils(unittest.TestCase):
    """Tests of path_utils"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def _slim_path_in_cesm(self):
        """Returns the path to a slim directory nested inside a typical cesm
        directory structure, where self._testdir is the root of the cesm
        checkout
        """
        return os.path.join(self._testdir, "components", "slim")

    def _cime_path_in_cesm(self):
        """Returns the path to a cime directory nested inside a typical
        cesm directory structure, where self._testdir is the root of the
        cesm checkout
        """
        return os.path.join(self._testdir, "cime")

    def _make_cesm_dirs(self):
        """Makes a directory structure for a typical CESM layout, where
        self._testdir is the root of the CESM checkout.

        This makes the slim directory and the cime directory.

        Returns a tuple, (slim_path, cime_path)
        """
        slim_path = self._slim_path_in_cesm()
        cime_path = self._cime_path_in_cesm()
        os.makedirs(slim_path)
        os.makedirs(cime_path)
        return (slim_path, cime_path)

    def test_pathToCime_standaloneOnlyWithCime(self):
        """Test path_to_cime with standalone_only, where cime is in the location
        it should be with a standalone checkout
        """
        slim_path = os.path.join(self._testdir, "slim")
        actual_path_to_cime = os.path.join(slim_path, "cime")
        os.makedirs(actual_path_to_cime)

        with mock.patch("slim.path_utils.path_to_slim_root", return_value=slim_path):
            path_to_cime = path_utils.path_to_cime(standalone_only=True)

        self.assertEqual(path_to_cime, actual_path_to_cime)

    def test_pathToCime_standaloneOnlyWithoutCime(self):
        """Test path_to_cime with standalone_only, where cime is missing"""
        slim_path = os.path.join(self._testdir, "slim")
        os.makedirs(slim_path)

        with mock.patch("slim.path_utils.path_to_slim_root", return_value=slim_path):
            with self.assertRaisesRegex(RuntimeError, "Cannot find cime"):
                _ = path_utils.path_to_cime(standalone_only=True)

    def test_pathToCime_standaloneOnlyWithCimeInCesm(self):
        """Test path_to_cime with standalone_only, where cime is missing from
        the standalone structure, but cime is present in the CESM
        directory structure: should raise an exception rather than
        finding that cime
        """
        slim_path, _ = self._make_cesm_dirs()

        with mock.patch("slim.path_utils.path_to_slim_root", return_value=slim_path):
            with self.assertRaisesRegex(RuntimeError, "Cannot find cime"):
                _ = path_utils.path_to_cime(standalone_only=True)

    def test_pathToCime_cimeInCesm(self):
        """Test path_to_cime, where cime is not in the standalone directory but
        is present in the CESM directory structure
        """
        slim_path, actual_path_to_cime = self._make_cesm_dirs()

        with mock.patch("slim.path_utils.path_to_slim_root", return_value=slim_path):
            path_to_cime = path_utils.path_to_cime()

        self.assertEqual(path_to_cime, actual_path_to_cime)

    def test_pathToCime_notInCesmCheckout(self):
        """Test path_to_cime, where cime is not in the standalone directory, and
        we don't appear to be in a CESM checkout
        """
        slim_path = os.path.join(self._testdir, "components", "foo")
        os.makedirs(slim_path)
        os.makedirs(self._cime_path_in_cesm())

        with mock.patch("slim.path_utils.path_to_slim_root", return_value=slim_path):
            with self.assertRaisesRegex(
                RuntimeError,
                "Cannot find cime.*don't seem to be within a CESM checkout",
            ):
                _ = path_utils.path_to_cime()

    def test_pathToCime_noCimeInCesm(self):
        """Test path_to_cime, where we appear to be within a CESM checkout, but
        there is no cime directory"""
        slim_path = self._slim_path_in_cesm()
        os.makedirs(slim_path)

        with mock.patch("slim.path_utils.path_to_slim_root", return_value=slim_path):
            with self.assertRaisesRegex(RuntimeError, "Cannot find cime.*or within CESM checkout"):
                _ = path_utils.path_to_cime()

    def test_pathToCime_cimeInStandaloneAndCesm(self):
        """Test path_to_cime, where there is a cime directory both in the
        standalone checkout and in the enclosing CESM checkout. Should
        give us the cime in the standalone checkout.
        """
        slim_path, _ = self._make_cesm_dirs()
        actual_path_to_cime = os.path.join(slim_path, "cime")
        os.makedirs(actual_path_to_cime)

        with mock.patch("slim.path_utils.path_to_slim_root", return_value=slim_path):
            path_to_cime = path_utils.path_to_cime()

        self.assertEqual(path_to_cime, actual_path_to_cime)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
