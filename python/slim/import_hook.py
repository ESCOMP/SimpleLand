# pylint: disable=missing-module-docstring
import sys
import os

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "cime")
sys.path.append(_CIMEROOT)
