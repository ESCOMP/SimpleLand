"""General-purpose utility functions"""

import logging
import os
import sys
import string
import pdb

from datetime import date
from getpass import getuser

logger = logging.getLogger(__name__)

def abort(errmsg):
    """Abort the program with the given error message

    No traceback is given, but if the logging level is DEBUG, then we'll enter pdb
    """
    if logger.isEnabledFor(logging.DEBUG):
        pdb.set_trace()

    sys.exit("ERROR: {}".format(errmsg))
