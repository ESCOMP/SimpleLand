"""Adds slim_cime_py lib to path

Any file that can potentially be run as a top-level script (with an if __name__ ==
'__main__' block) that needs slim_cime_py should import this. This includes unit test
modules.

This should be the very first slim module imported. That way, slim_cime_py will be added to your
path before other imports. That is, your script should have:

# Standard library imports go here
# Then:
from slim import add_slim_cime_py_to_path
"""

from slim.path_utils import add_slim_cime_pylib_to_path

_ = add_slim_cime_pylib_to_path()
