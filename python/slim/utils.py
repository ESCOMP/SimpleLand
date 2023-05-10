"""General-purpose utility functions"""

import logging
import os
import sys
import pdb

from datetime import date
from getpass import getuser

from slim.git_utils import get_slim_git_short_hash

logger = logging.getLogger(__name__)


def abort(errmsg):
    """Abort the program with the given error message

    No traceback is given, but if the logging level is DEBUG, then we'll enter pdb
    """
    if logger.isEnabledFor(logging.DEBUG):
        pdb.set_trace()

    sys.exit("ERROR: {}".format(errmsg))


def update_metadata(file, title, summary, contact, data_script, description):
    """
    Description
    -----------
    Update netcdf file's metadata
    Arguments
    ---------
    title: No more than short one-sentence explanation.
    summary: No more than two-sentence explanation.
    contact: E.g. CAM bulletin board at https://bb.cgd.ucar.edu
    data_script: Script or instructions used to generate the dataset.
    description: Anything else that's relevant. Capturing the command-line
                 would be good (sys.argv) here or in data_script.
    """

    # update attributes
    today = date.today()
    today_string = today.strftime("%Y-%m-%d")

    # This is the required metadata for inputdata files
    file.attrs["title"] = title
    file.attrs["summary"] = summary
    file.attrs["creator"] = getuser()
    file.attrs["contact"] = contact
    file.attrs["creation_date"] = today_string
    file.attrs["data_script"] = data_script
    file.attrs["description"] = description

    # delete unrelated attributes if they exist
    del_attrs = [
        "source_code",
        "SVN_url",
        "hostname",
        "history",
        "History_Log",
        "Logname",
        "Host",
        "Version",
        "Compiler_Optimized",
    ]
    attr_list = file.attrs

    for attr in del_attrs:
        if attr in attr_list:
            del file.attrs[attr]


def write_output(file, file_in, file_out, file_type):
    """
    Description
    -----------
    Write output file
    Arguments
    ---------
    file_in:
        (str) User-defined entry of input file
    file_out:
        (str) User-defined entry of output file
    file_type:
        (str) examples: mesh, surdat
    """

    # update attributes
    title = "Modified " + file_type + " file"
    summary = "Modified " + file_type + " file"
    contact = "N/A"
    data_script = os.path.abspath(__file__) + " -- " + get_slim_git_short_hash()
    description = "Modified this file: " + file_in
    update_metadata(
        file,
        title=title,
        summary=summary,
        contact=contact,
        data_script=data_script,
        description=description,
    )

    # mode 'w' overwrites file if it exists
    file.to_netcdf(path=file_out, mode="w", format="NETCDF3_64BIT")
    logger.info("Successfully created: %s", file_out)
    file.close()
