"""
SLIM namelist creator
"""
import sys, os, shutil

_CIMEROOT = os.environ.get("CIMEROOT")
_LIBDIR = os.path.join(_CIMEROOT, "scripts", "Tools")
sys.path.append(_LIBDIR)

from standard_script_setup          import *
from CIME.buildnml                  import create_namelist_infile, parse_input
from CIME.nmlgen                    import NamelistGenerator
from CIME.case                      import Case
from CIME.utils                     import expect, run_cmd

logger = logging.getLogger(__name__)

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_dtime( nmlgen ):
####################################################################################
    """ Set the namelist settings for time-step
    """
    global logger
    #------------------------------------------------------
    logger.info( " check_nml_dtime" )

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_general( nmlgen ):
####################################################################################
    """ Set the namelist settings for general settings
    """
    global logger
    #------------------------------------------------------
    logger.info( " check_nml_general" )

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_performance( nmlgen ):
####################################################################################
    """ Set the namelist settings for performance
    """
    global logger
    #------------------------------------------------------
    logger.info( " check_nml_performance" )

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_history( nmlgen ):
####################################################################################
    """ Set the namelist settings for history
    """
    global logger
    #------------------------------------------------------
    logger.info( " check_nml_history" )

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_initial_conditions( nmlgen ):
####################################################################################
    """ Set the namelist settings for initial conditions
    """
    global logger
    #------------------------------------------------------
    logger.info( " check_nml_initial_conditions" )

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def check_nml_data( nmlgen ):
####################################################################################
    """ Set the namelist settings for data
    """
    global logger
    #------------------------------------------------------
    logger.info( " check_nml_data" )

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches,too-many-statements
####################################################################################
def _create_namelists(case, confdir, inst_string, infile, nmlgen, data_list_path):
####################################################################################
    """Write out the namelist for this component.

    Most arguments are the same as those for `NamelistGenerator`. The
    `inst_string` argument is used as a suffix to distinguish files for
    different instances. The `confdir` argument is used to specify the directory
    in which output files will be placed.
    """
    global logger
    #------------------------------------------------------
    # Create config dictionary
    #------------------------------------------------------
    config = {}
    config['lnd_grid']      = case.get_value("LND_GRID")
    config['slim_scenario'] = case.get_value("SLIM_SCENARIO")

    logger.info( " SLIM lnd grid is %s", config['lnd_grid'] )

    #------------------------------------------------------
    # Initialize namelist defaults
    #------------------------------------------------------
    nmlgen.init_defaults(infile, config)

    #------------------------------------------------------
    #  Process different namelists and parts of the namelist
    #------------------------------------------------------
    check_nml_dtime( nmlgen ) 
    check_nml_general( nmlgen ) 
    check_nml_performance( nmlgen )
    check_nml_history( nmlgen )
    check_nml_initial_conditions( nmlgen )
    check_nml_data( nmlgen )

    #----------------------------------------------------
    # Write output namelist
    #----------------------------------------------------
    namelist_file = os.path.join(confdir, "lnd_in")
    nmlgen.write_output_file(namelist_file, data_list_path, \
                             groups=['slim_inparm', 'slim_data_and_initial', 'slim_history', 'slim_perf'])


###############################################################################
def buildnml(case, caseroot, compname):
###############################################################################
    """Build the slim namelist """
    global logger

    # Build the component namelist 
    if compname != "slim":
        raise AttributeError

    lnd_root = case.get_value("COMP_ROOT_DIR_LND")

    # -----------------------------------------------------
    # Clear out old data
    # -----------------------------------------------------

    input_data_list = os.path.join(caseroot,"Buildconf","slim.input_data_list")
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
        expect(os.path.isfile(file_), "Namelist XML file %s not found!" % file_ )

    # Create the namelist generator object - independent of instance
    nmlgen = NamelistGenerator(case, definition_file)

    #----------------------------------------------------
    # Clear out old data list
    #----------------------------------------------------
    data_list_path = os.path.join(case.get_case_root(), "Buildconf", "slim.input_data_list")
    if os.path.exists(data_list_path):
        os.remove(data_list_path)

    ### Independent of instance...
    startfile_type = "finidat"
    start_type = "default"
    run_type = case.get_value("RUN_TYPE")
    if run_type == "hybrid":
        start_type = "startup"
    elif run_type != "startup":
        start_type = run_type

    slim_force_coldstart = case.get_value("SLIM_FORCE_COLDSTART")
    if run_type == "branch":
        startfile_type = "nrevsn"
        if slim_force_coldstart == "on":
           slim_force_coldstart = "off"
           logger.warning( "WARNING: You've turned on SLIM_FORCE_COLDSTART for a branch run_type, which is a contradiction, the coldstart will be ignored\n" +
                           "  turn off SLIM_FORCE_COLDSTART, or set RUN_TYPE=hybrid to get rid of this warning"
                         )

    if (slim_force_coldstart == "on"):
        logger.warning( "WARNING: SLIM is starting up from a cold state" )
        start_type = "cold"


    run_startdate = case.get_value("RUN_STARTDATE")
    start_ymd = run_startdate.replace('-','')

    inputdata_file = os.path.join(caseroot,"Buildconf","slim.input_data_list")
    
    rundir = case.get_value("RUNDIR")
    
    # -----------------------------------------------------
    # loop over instances
    # -----------------------------------------------------

    ninst_lnd = case.get_value("NINST_LND")
    ninst = int(ninst_lnd)
    for inst_counter in range(1, ninst+1):

        # determine instance string
        inst_string = ""
        if ninst > 1:
            inst_string = '_' + '%04d' % inst_counter

        # If multi-instance case does not have restart file, use
        # single-case restart for each instance
        rpointer = "rpointer.lnd" 
        if (os.path.isfile(os.path.join(rundir,rpointer)) and
            (not os.path.isfile(os.path.join(rundir,rpointer + inst_string)))):
            shutil.copy(os.path.join(rundir, rpointer),
                        os.path.join(rundir, rpointer + inst_string))

        ###
        ### instance dependent information...
        ###
        run_refcase = case.get_value("RUN_REFCASE")
        run_refdate = case.get_value("RUN_REFDATE")
        run_reftod = case.get_value("RUN_REFTOD")
        rundir     = case.get_value("RUNDIR")
        if run_type == "hybrid" or run_type == "branch":
            slim_startfile = "%s.slim%s.r.%s-%s.nc"%(run_refcase,inst_string,run_refdate,run_reftod)
            if not os.path.exists(os.path.join(rundir, slim_startfile)):
                slim_startfile = "%s.slim.r.%s-%s.nc"%(run_refcase,run_refdate,run_reftod)
            slim_icfile = "%s = \'%s\'"%(startfile_type, slim_startfile)
        else:
            slim_icfile = ""
        ###

        infile_lines = []
        infile_lines.append(slim_icfile)

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
            logger.info("SLIM namelist copy: file1 %s file2 %s " %(file1, file2))
            shutil.copy(file1,file2)
