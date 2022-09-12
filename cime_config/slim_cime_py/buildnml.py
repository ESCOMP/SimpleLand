"""
SLIM namelist creator
"""
import os
import shutil
import logging

from CIME.buildnml import create_namelist_infile
from CIME.nmlgen import NamelistGenerator
from CIME.utils import expect

logger = logging.getLogger(__name__)

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_dtime(nmlgen, case):
    ####################################################################################
    """Set the namelist settings for time-step"""
    # pylint: disable=global-statement
    global logger
    # ------------------------------------------------------
    logger.info(" check_nml_dtime")
    ncpl_base_period = case.get_value("NCPL_BASE_PERIOD")
    if ncpl_base_period == "hour":
        basedt = 3600
    elif ncpl_base_period == "day":
        basedt = 3600 * 24
    elif ncpl_base_period == "year":
        if case.get_value("CALENDAR") == "NO_LEAP":
            basedt = 3600 * 24 * 365
        else:
            expect(False, "Invalid CALENDAR for NCPL_BASE_PERIOD %s " % ncpl_base_period)
    elif ncpl_base_period == "decade":
        if case.get_value("CALENDAR") == "NO_LEAP":
            basedt = 3600 * 24 * 365 * 10
        else:
            expect(
                False,
                "invalid NCPL_BASE_PERIOD NCPL_BASE_PERIOD %s " % ncpl_base_period,
            )
    else:
        expect(False, "invalid NCPL_BASE_PERIOD NCPL_BASE_PERIOD %s " % ncpl_base_period)

    logger.info(" basedt = %s", str(basedt))

    if basedt < 0:
        expect(False, "basedt invalid overflow for NCPL_BASE_PERIOD %s " % ncpl_base_period)

    lnd_ncpl = int(case.get_value("LND_NCPL"))
    if basedt % lnd_ncpl != 0:
        expect(
            False,
            "lnd_ncpl %s doesn't divide evenly into basedt %s\n" % (lnd_ncpl, basedt),
        )
    else:
        dtime = basedt // lnd_ncpl
    nmlgen.set_value("dtime", value=dtime)


# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_general(nmlgen):
    ####################################################################################
    """Set the namelist settings for general settings"""
    # pylint: disable=global-statement
    global logger
    # ------------------------------------------------------
    logger.info(" check_nml_general")
    for var in ("slim_start_type", "res"):
        expect(nmlgen.get_value(var) is not None, var + " must be set")


####################################################################################
def check_file(filename, case):
    ####################################################################################
    """Check that the file exists"""
    if os.path.isabs(filename):
        expect(os.path.isfile(filename), "filename must exist:" + filename)
    else:
        rundir = case.get_value("RUNDIR")
        fname = os.path.normpath(os.path.join(rundir, filename))
        expect(os.path.isfile(fname), "filename must exist:" + fname)


# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_performance(nmlgen):
    ####################################################################################
    """Set the namelist settings for performance"""
    # pylint: disable=global-statement
    global logger
    # ------------------------------------------------------
    logger.info(" check_nml_performance")
    expect(int(nmlgen.get_value("nsegspc")) > 0, "nsegspc must be positive")


# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_history(nmlgen):
    ####################################################################################
    """Set the namelist settings for history"""
    # pylint: disable=global-statement
    global logger
    # ------------------------------------------------------
    logger.info(" check_nml_history")

    hist_mfilt = nmlgen.get_value("hist_mfilt")
    for mfilt in hist_mfilt:
        logger.info(" hist_mfilt = %d", int(mfilt))
        if int(mfilt) <= 0:
            raise SystemExit("hist_mfilt must be 1 or larger")


# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_initial_conditions(nmlgen, case, inst_string=""):
    ####################################################################################
    """Set the namelist settings for initial conditions"""
    # pylint: disable=global-statement
    global logger
    # ------------------------------------------------------
    logger.info(" check_nml_initial_conditions")
    start_type = case.get_value("SLIM_START_TYPE")
    run_type = case.get_value("RUN_TYPE")
    run_refcase = case.get_value("RUN_REFCASE")
    run_refdate = case.get_value("RUN_REFDATE")
    run_reftod = case.get_value("RUN_REFTOD")
    rundir = case.get_value("RUNDIR")
    if run_type in ("hybrid", "branch"):
        slim_startfile = "%s.slim%s.r.%s-%s.nc" % (
            run_refcase,
            inst_string,
            run_refdate,
            run_reftod,
        )
        if not os.path.exists(os.path.join(rundir, slim_startfile)):
            slim_startfile = "%s.slim.r.%s-%s.nc" % (
                run_refcase,
                run_refdate,
                run_reftod,
            )

    nrevsn = nmlgen.get_value("nrevsn")
    finidat = nmlgen.get_value("finidat")
    #
    # Non branch types
    #
    if run_type != "branch":
        # Handle a cold start
        if start_type == "cold":
            logger.info(" finidat = %s", finidat)
            if finidat != " " and finidat != "UNSET" and finidat is not None:
                print(" finidat = '" + finidat + "'")
                raise SystemExit(
                    "finidat is set but SLIM_START_TYPE is cold which is a contradiction"
                )
            nmlgen.set_value("finidat", value=" ")

        # Set to blank meaning a cold start if still UNSET
        if finidat == " " or finidat == "UNSET" or finidat is None:
            if run_type == "hybrid" and start_type != "cold":
                finidat = slim_startfile
                check_file(finidat, case)
                nmlgen.set_value("finidat", value=finidat)
            else:
                nmlgen.set_value("finidat", value=" ")
                logger.warning("WARNING: SLIM is starting up from a cold state")

        else:
            check_file(finidat, case)

        if nrevsn is not None:
            raise SystemExit("nrevsn can NOT be set except when RUN_TYPE is a branch")
    #
    # branch types
    #
    else:
        if nrevsn is None:
            nrevsn = slim_startfile
            nmlgen.set_value("nrevsn", value=nrevsn)

        check_file(nrevsn, case)
        if finidat is not None:
            print(finidat)
            raise SystemExit("finidat can NOT be set when RUN_TYPE is a branch")


# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_data(nmlgen):
    ####################################################################################
    """Set the namelist settings for data"""
    # pylint: disable=global-statement
    global logger
    # ------------------------------------------------------
    logger.info(" check_nml_data")

    mml_surdat = nmlgen.get_value("mml_surdat")
    if mml_surdat == "UNSET":
        raise SystemExit("mml_surdat file is NOT set and is required")


# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
# Turn off unused-argument for inst_string, since isn't in place right now
# pylint: disable=unused-argument
####################################################################################
def _create_namelists(case, confdir, inst_string, infile, nmlgen, data_list_path):
    ####################################################################################
    """Write out the namelist for this component.

    Most arguments are the same as those for `NamelistGenerator`. The
    `inst_string` argument is used as a suffix to distinguish files for
    different instances. The `confdir` argument is used to specify the directory
    in which output files will be placed.
    """
    # pylint: disable=global-statement
    global logger
    # ------------------------------------------------------
    # Create config dictionary
    # ------------------------------------------------------
    config = {}
    config["lnd_grid"] = case.get_value("LND_GRID")
    config["slim_scenario"] = case.get_value("SLIM_SCENARIO")
    config["slim_start_type"] = case.get_value("SLIM_START_TYPE")

    logger.info(" SLIM lnd grid is %s", config["lnd_grid"])

    # ------------------------------------------------------
    # Initialize namelist defaults
    # ------------------------------------------------------
    nmlgen.init_defaults(infile, config)

    # ------------------------------------------------------
    #  Process different namelists and parts of the namelist
    # ------------------------------------------------------
    check_nml_dtime(nmlgen, case)
    check_nml_general(nmlgen)
    check_nml_performance(nmlgen)
    check_nml_history(nmlgen)
    check_nml_initial_conditions(nmlgen, case, inst_string)
    check_nml_data(nmlgen)

    # ----------------------------------------------------
    # Write output namelist
    # ----------------------------------------------------
    logger.info("Write namelists")
    namelist_file = os.path.join(confdir, "lnd_in")
    nmlgen.write_output_file(
        namelist_file,
        data_list_path,
        groups=["slim_inparm", "slim_data_and_initial", "slim_history", "slim_perf"],
    )


###############################################################################
def buildnml(case, caseroot, compname):
    ###############################################################################
    """Build the slim namelist"""
    # pylint: disable=global-statement
    global logger

    # Build the component namelist
    if compname != "slim":
        raise AttributeError

    lnd_root = case.get_value("COMP_ROOT_DIR_LND")

    # -----------------------------------------------------
    # Clear out old data
    # -----------------------------------------------------

    input_data_list = os.path.join(caseroot, "Buildconf", "slim.input_data_list")
    if os.path.exists(input_data_list):
        os.remove(input_data_list)

    # -----------------------------------------------------
    # Set confdir
    # -----------------------------------------------------

    confdir = os.path.join(caseroot, "Buildconf", "slimconf")
    if not os.path.isdir(confdir):
        os.makedirs(confdir)

    # namelist definition file
    namelist_xml_dir = os.path.join(lnd_root, "cime_config")
    definition_file = [os.path.join(namelist_xml_dir, "namelist_definition_slim.xml")]
    for file_ in definition_file:
        expect(os.path.isfile(file_), "Namelist XML file %s not found!" % file_)

    # Create the namelist generator object - independent of instance
    nmlgen = NamelistGenerator(case, definition_file)

    # ----------------------------------------------------
    # Clear out old data list
    # ----------------------------------------------------
    data_list_path = os.path.join(case.get_case_root(), "Buildconf", "slim.input_data_list")
    if os.path.exists(data_list_path):
        os.remove(data_list_path)

    ### Independent of instance...
    rundir = case.get_value("RUNDIR")

    # -----------------------------------------------------
    # loop over instances
    # -----------------------------------------------------

    ninst_lnd = case.get_value("NINST_LND")
    ninst = int(ninst_lnd)
    for inst_counter in range(1, ninst + 1):

        # determine instance string
        inst_string = ""
        if ninst > 1:
            inst_string = "_" + "%04d" % inst_counter

        # If multi-instance case does not have restart file, use
        # single-case restart for each instance
        rpointer = "rpointer.lnd"
        if os.path.isfile(os.path.join(rundir, rpointer)) and (
            not os.path.isfile(os.path.join(rundir, rpointer + inst_string))
        ):
            shutil.copy(
                os.path.join(rundir, rpointer),
                os.path.join(rundir, rpointer + inst_string),
            )
        ###
        ### Namelist infile
        ###
        infile_lines = []

        user_nl_file = os.path.join(caseroot, "user_nl_slim" + inst_string)
        infile = os.path.join(confdir, "namelist_infile")

        create_namelist_infile(case, user_nl_file, infile, "\n".join(infile_lines))
        namelist_infile = [infile]

        # create namelist
        _create_namelists(case, confdir, inst_string, namelist_infile, nmlgen, data_list_path)
        # -----------------------------------------------------
        # copy resolved namelist to rundir
        # -----------------------------------------------------
        if os.path.isdir(rundir):
            file1 = os.path.join(confdir, "lnd_in")
            file2 = os.path.join(rundir, "lnd_in")
            if ninst > 1:
                file2 += inst_string
            logger.info("SLIM namelist copy: file1 %s file2 %s ", file1, file2)
            shutil.copy(file1, file2)
