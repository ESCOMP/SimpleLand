#pylint: disable=missing-module-docstring
import sys
import os
_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..", "cime")
_LIB_DIR = os.path.join(_CIMEROOT, "scripts", "lib")
sys.path.append(_LIB_DIR)
