#!/usr/bin/env python2
"""
jobsub: a tool for EUTelescope job submission to Marlin 

Steering files are generated based on steering templates by substituting
job and run specific variables. The variable information can be
provided to jobsub by command line argument, by a config file or by a
text file with run parameters (in comma-separated value/csv format).

Run 
python jobsub.py --help
to see the list of command line options.

"""
import sys
import logging

def parseIntegerString(nputstr=""):
    """
    return a list of selected values when a string in the form:
    1-4,6
    would return:
    1,2,3,4,6
    as expected...
    (from http://thoughtsbyclayg.blogspot.de/2008/10/parsing-list-of-numbers-in-python.html)

    """
    selection = list()
    # tokens are comma seperated values
    tokens = [substring.strip() for substring in nputstr.split(',')]
    for i in tokens:
        try:
            # typically tokens are plain old integers
            selection.append(int(i))
        except ValueError:
            # if not, then it might be a range
            token = [int(k.strip()) for k in i.split('-')]
            if len(token) > 1:
                token.sort()
                # we have items seperated by a dash
                # try to build a valid range
                first = token[0]
                last = token[len(token)-1]
                for value in range(first, last+1):
                    selection.append(value)
    return selection # end parseIntegerString

def ireplace(old, new, text):
    """ 
    case insensitive search and replace function searching through string and returning the filtered string
    (based on http://stackoverflow.com/a/4773614)

    """
    idx = 0
    occur = 0
    while idx < len(text):
        index_l = text.lower().find(old.lower(), idx)
        if index_l == -1:
            if occur == 0:
                raise EOFError("Could not find string "+old)
            return text
        text = text[:index_l] + new + text[index_l + len(old):]
        idx = index_l + len(new)
        occur = occur+1
    if occur == 0:
        raise EOFError("Could not find string "+old)
    return text


def loadparamsfromcsv(csvfilename, runs):
    """ Load and parse the csv file for the given set of runs and
    return nested dictionary: a collection of dictionaries, one for
    each csv row matching a run number.

    """
    import csv
    import os.path
    from sys import exit # use sys.exit instead of built-in exit (latter raises exception)

    class CommentedFile:
        """ Decorator for text files: filters out comments (i.e. first char of line #)
        Based on http://www.mfasold.net/blog/2010/02/python-recipe-read-csvtsv-textfiles-and-ignore-comment-lines/
        
        """
        def __init__(self, f, commentstring="#"):
            self.f = f
            self.commentstring = commentstring
            self.linecount = 0
        def rewind(self):
            self.f.seek(0)
            self.linecount = 0
        def next(self):
            line = self.f.next()
            self.linecount += 1
            while line.startswith(self.commentstring) or not line.strip(): # test if line commented or empty
                line = self.f.next()
                self.linecount += 1
            return line
        def __iter__(self):
            return self

    log = logging.getLogger('jobsub')
    parameters_csv = {} # store all information needed from the csv file
    if csvfilename is None: 
        return parameters_csv # if no file name given, return empty collection here
    if not os.path.isfile(csvfilename): # check if file exists
        log.error("Could not find the specified csv file '"+csvfilename+"'!")
        exit(1)
    try:
        log.debug("Opening csv file '"+csvfilename+"'.")
        csvfile = open(csvfilename, 'rb')
        filteredfile = CommentedFile(csvfile)
        try:
            # contruct a sample for the csv format sniffer:
            sample = ""
            try:
                while (len(sample)<1024):
                    sample += filteredfile.next()
            except StopIteration:
                log.debug("End of csv file reached, sample limited to " + str(len(sample))+ " bytes")
            dialect = csv.Sniffer().sniff(sample) # test csv file format details
            log.debug("Determined the CSV dialect as follows: delimiter=%s, doublequote=%s, escapechar=%s, lineterminator=%s, quotechar=%s , quoting=%s, skipinitialspace=%s", dialect.delimiter, dialect.doublequote, dialect.escapechar, list(ord(c) for c in dialect.lineterminator), dialect.quotechar, dialect.quoting, dialect.skipinitialspace)
            filteredfile.rewind() # back to beginning of file
            reader = csv.DictReader(filteredfile, dialect=dialect) # now process CSV file contents here and load them into memory
            reader.next() # python < 2.6 requires an actual read access before filling 'DictReader.fieldnames'
            log.debug("CSV file contains the header info: %s", reader.fieldnames)
            try:
                reader.fieldnames = [field.lower() for field in reader.fieldnames] # convert to lower case keys to avoid confusion
                reader.fieldnames = [field.strip() for field in reader.fieldnames] # remove leading and trailing white space
            except TypeError:
                log.error("Could not process the CSV file header information. csv.DictReader returned fieldnames: %s", reader.fieldnames)
                exit(1)
            if not "runnumber" in reader.fieldnames: # verify that we have a column "runnumber"
                log.error("Could not find a column with header label 'RunNumber' in file '"+csvfilename+"'!")
                exit(1)
            if "" in reader.fieldnames:
                log.warning("Column without header label encountered in csv file '"+csvfilename+"'!")
            log.info("Successfully loaded csv file'"+csvfilename+"'.")
            # first: search through csv file to find corresponding runnumber entry line for every run
            filteredfile.rewind() # back to beginning of file
            reader.next()   # .. and skip the header line
            missingRuns = list(runs) # list of runs to look for in csv file
            for row in reader: # loop over all rows once
                try:
                    for run in missingRuns: # check all runs if runnumber matches
                        if int(row["runnumber"]) == run:
                            log.debug("Found entry in csv file for run "+str(run)+" on line "+ str(filteredfile.linecount))
                            parameters_csv[run] = {}
                            parameters_csv[run].update(row)
                            missingRuns.remove(run)
                            break
                except ValueError: # int conversion error
                    log.warn("Could not interpret run number on line "+str(filteredfile.linecount)+" in file '"+csvfilename+"'.")
                    continue
                if len(missingRuns)==0:
                    log.debug("Finished search for runs in csv file before reaching end of file")
                    break
            log.debug("Searched over "+str(filteredfile.linecount)+" lines in file '"+csvfilename+"'.")
            if not len(missingRuns)==0:
                log.error("Could not find an entry for the following run numbers in '"+csvfilename+"': "+', '.join(map(str, missingRuns)))
        finally:
            csvfile.close()
    except csv.Error, e:
        log.error("Problem loading the csv file '"+csvfilename+"'(%s): %s"%(e.errno, e.strerror))
        exit(1)
    return parameters_csv

def checkSteer(sstring):
    """ Check string for any occurance of @.*@ and return boolean. """
    log = logging.getLogger('jobsub')
    import re
    hits = re.findall("@.*@", sstring)
    if hits:
        log.error ("Missing configuration parameters: "+', '.join(map(str, hits)))
        return False
    else:
        return True

def check_program(name):
    """ Searches PATH environment variable for executable given by parameter """
    import os
    for dir in os.environ['PATH'].split(os.pathsep):
        prog = os.path.join(dir, name)
        if os.path.exists(prog): return prog

def runMarlin(filenamebase, jobtask, silent):
    """ Runs Marlin and stores log of output """
    from sys import exit # use sys.exit instead of built-in exit (latter raises exception)
    log = logging.getLogger('jobsub.' + jobtask)

    # check for Marlin executable
    cmd = check_program("Marlin")
    if cmd:
        log.debug("Found Marlin executable: " + cmd)
    else:
        log.error("Marlin executable not found in PATH!")
        exit(1)

    # search for stdbuf command: adjust stdout buffering
    stdbuf = check_program("stdbuf")
    if stdbuf:
        log.debug("Found stdbuf, will use line buffered output.")
        # -oL: adjust standard output stream buffering to line buffered
        cmd = stdbuf + " -oL " + cmd

    # need some addtional libraries for process interaction
    from subprocess import Popen, PIPE
    from threading  import Thread # threading used for non-blocking process output parsing
    from time import sleep
    try:
        from Queue import Queue, Empty # python 2.x
    except ImportError:
        from queue import Queue, Empty  # python 3.x

    import datetime
    import shlex        

    # parsing process output using threads
    # (approach from http://stackoverflow.com/a/4896288)
    def enqueue_output(out, queue):
        """ feed queue with readline output """
        for line in iter(out.readline, ''):
            queue.put(line)
        out.close()
    ON_POSIX = 'posix' in sys.builtin_module_names
    cmd = cmd+" "+filenamebase+".xml"
    rcode = None # the return code that will be set by a later subprocess method
    try:
        # run process
        log.info ("Now running Marlin on "+filenamebase+".xml")
        log.debug ("Executing: "+cmd)
        p = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE, bufsize=1, close_fds=ON_POSIX)
        # setup output queues and threads
        qout = Queue()
        tout = Thread(target=enqueue_output, args=(p.stdout, qout))
        qerr = Queue()
        terr = Thread(target=enqueue_output, args=(p.stderr, qerr))
        # threads die with the program
        tout.daemon = True
        terr.daemon = True 
        tout.start()
        terr.start()
        # open log file
        log_file = open(filenamebase+".log", "w")
        # print timestamp to log file
        log_file.write("---=== Analysis started on " + datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p") + " ===---\n\n")
        try:
            while p.poll() is None:
                # read line without blocking
                try:  
                    line = qout.get_nowait() # or q.get(timeout=.1)
                    if not silent:
                        if "WARNING" in line.strip():
                            log.warning(line.strip())
                        elif "ERROR" in line.strip():
                            log.error(line.strip())
                        else:
                            log.info(line.strip())
                    log_file.write(line)
                except Empty:
                    pass
                
                try:  
                    line = qerr.get_nowait() # or q.get(timeout=.1)
                    log.error(line.strip())
                    log_file.write(line)                     
                except Empty:
                    sleep(0.005) # sleep for 5 ms to avoid excessive CPU load

            # process done
            tout.join() # finish stdout thread; wait for remaining buffer to be read
            terr.join() # finish stderr thread
            # process the remainder of the buffers now stored in our queues
            while not qout.empty() or not qerr.empty():
                # read line without blocking
                try:  
                    line = qout.get_nowait() # or q.get(timeout=.1)
                    if not silent:
                        if "WARNING" in line.strip():
                            log.warning(line.strip())
                        elif "ERROR" in line.strip():
                            log.error(line.strip())
                        else:
                            log.info(line.strip())
                    log_file.write(line)
                except Empty:
                    pass
                
                try:  
                    line = qerr.get_nowait() # or q.get(timeout=.1)
                    log.error(line.strip())
                    log_file.write(line)                     
                except Empty:
                    pass
        finally:
            log_file.close()
        rcode = p.returncode # get the return code
    except OSError, e:
        log.critical("Problem with Marlin execution: Command '%s' resulted in error #%s, %s", cmd, e.errno, e.strerror)
        exit(1)
    return rcode

def submitNAF(filenamebase, jobtask, qsubfile, runnr):
    """ Submits the Marlin job to NAF """
    import os
    from sys import exit # use sys.exit instead of built-in exit (latter raises exception)
    log = logging.getLogger('jobsub.' + jobtask)
    # We are running on NAF.

    # check for qsub executable
    cmd = check_program("qsub")
    if cmd:
        log.debug("Found qsub executable: " + cmd)
    else:
        log.error("qsub executable not found in PATH!")
        exit(1)

    # Add qsub parameters:
    #qsub -@ qsubParams.txt BIN
    cmd = cmd+" -@ "+qsubfile+" -N \"Run"+runnr+"\" "
    
    # check for Marlin executable
    marlin = check_program("Marlin")
    if marlin:
        log.debug("Found Marlin executable: " + marlin)
        cmd = cmd+" "+marlin
    else:
        log.error("Marlin executable not found in PATH!")
        exit(1)

    cmd = cmd+" "+filenamebase+".xml"
    rcode = None # the return code that will be set by a later subprocess method
    try:
        # run process
        log.info ("Now submitting Marlin job: "+filenamebase+".xml to NAF")
        log.debug ("Executing: "+cmd)
        os.popen(cmd)
    except OSError, e:
        log.critical("Problem with NAF submission: Command '%s' resulted in error #%s, %s", cmd, e.errno, e.strerror)
        exit(1)
    return 0

def submitLXPLUS(filenamebase, jobtask, bsubfile, runnr):
    """ Submits the Marlin job to LXPLUS """
    import os
    from sys import exit # use sys.exit instead of built-in exit (latter raises exception)
    log = logging.getLogger('jobsub.' + jobtask)
    # We are running on LXPLUS.

    # check for bsub executable
    cmd = check_program("bsub")
    if cmd:
        log.debug("Found bsub executable: " + cmd)
    else:
        log.error("bsub executable not found in PATH!")
        exit(1)

    # Add bsub parameters:
    #bsub < bsubparams.txt BIN
    cmd = cmd+" < "+bsubfile+" -J \"Run"+runnr+"\" "
    
    # check for Marlin executable
    marlin = check_program("Marlin")
    if marlin:
        log.debug("Found Marlin executable: " + marlin)
        cmd = cmd+" "+marlin
    else:
        log.error("Marlin executable not found in PATH!")
        exit(1)

    filename = os.path.abspath(filenamebase+".xml")
    cmd = cmd+" "+filename
    rcode = None # the return code that will be set by a later subprocess method
    try:
        # run process
        log.info ("Now submitting Marlin job: "+filenamebase+".xml to LXPLUS")
        log.debug ("Executing: "+cmd)
        os.popen(cmd)
    except OSError, e:
        log.critical("Problem with LXPLUS submission: Command '%s' resulted in error #%s, %s", cmd, e.errno, e.strerror)
        exit(1)
    return 0

def zipLogs(path, filename):
    """  stores output from Marlin in zip file; enables compression if necessary module is available """
    import zipfile
    import os.path
    log = logging.getLogger('jobsub')
    try:     # compression module might not be available, therefore try import here
        import zlib
        compression = zipfile.ZIP_DEFLATED
        log.debug("Creating *compressed* log archive")
    except ImportError: # no compression module available, use flat files
        compression = zipfile.ZIP_STORED
        log.debug("Creating flat log archive")
    try:
        zf = zipfile.ZipFile(os.path.join(path, filename)+".zip", mode='w') # create new zip file
        try:
            zf.write(os.path.join("./", filename)+".xml", compress_type=compression) # store in zip file
            zf.write(os.path.join("./", filename)+".log", compress_type=compression) # store in zip file
            os.remove(os.path.join("./", filename)+".xml") # delete file
            os.remove(os.path.join("./", filename)+".log") # delete file
            log.info("Logs written to "+os.path.join(path, filename)+".zip")
        finally:
            log.debug("Closing log archive file")
            zf.close()
    except IOError: # could not create zip file - path non-existant?!
        log.error("Input/Output error: Could not create log and steering file archive ("+os.path.join(path, filename)+".zip"+")!")


def main(argv=None):
    """  main routine of jobsub: a tool for EUTelescope job submission to Marlin """
    log = logging.getLogger('jobsub') # set up logging
    formatter = logging.Formatter('%(asctime)s %(name)s(%(levelname)s): %(message)s',"%H:%M:%S")
    handler_stream = logging.StreamHandler()
    handler_stream.setFormatter(formatter)
    log.addHandler(handler_stream)
    # using this decorator, we can count the number of error messages
    class callcounted(object):
        """Decorator to determine number of calls for a method"""
        def __init__(self,method):
            self.method=method
            self.counter=0
        def __call__(self,*args,**kwargs):
            self.counter+=1
            return self.method(*args,**kwargs)
    log.error=callcounted(log.error)

    import os.path
    import ConfigParser
    try:
        import argparse
    except ImportError:
        log.debug("No locally installed argparse module found; trying the package provided with jobsub.")
        # argparse is not installed; use (old) version provided with jobsub
        # determine path to subdirectory
        libdir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(__file__))),"pymodules","argparse")
        if libdir not in sys.path:
            sys.path.append(libdir)
        # try again loading the module
        try:
            import argparse
        except ImportError:
            # nothing we can do now
            log.critical("Could not load argparse module. For python versions prior to 2.7, please install it from http://code.google.com/p/argparse")
            return 1

    if argv is None:
        argv = sys.argv
        progName = os.path.basename(argv.pop(0))

    # command line argument parsing
    parser = argparse.ArgumentParser(prog=progName, description="A tool for the convenient run-specific modification of Marlin steering files and their execution through the Marlin processor")
    parser.add_argument('--version', action='version', version='Revision: $Revision$, $LastChangedDate$')
    parser.add_argument('--option', '-o', action='append', metavar="NAME=VALUE", help="Specify further options such as 'beamenergy=5.3'. This switch be specified several times for multiple options or can parse a comma-separated list of options. This switch overrides any config file options.")
    parser.add_argument("-c", "--conf-file", "--config", help="Load specified config file with global and task specific variables", metavar="FILE")
    parser.add_argument("-n", "--naf-file", "--naf", help="Specify qsub parameter file for NAF submission. Run NAF submission via qsub instead of calling Marlin directly", metavar="FILE")
    parser.add_argument("-lx", "--lxplus-file", "--lxplus", help="Specify bsub parameter file for LXPLUS submission. Run LXPLUS submission via bsub instead of calling Marlin directly", metavar="FILE")
    parser.add_argument("--concatenate", action="store_true", default=False, help="Modifies run range treatment: concatenate all runs into first run (e.g. to combine runs for alignment) by combining every options that includes the string '@RunRange@' multiple times, once for each run of the range specified.")
    parser.add_argument("-csv", "--csv-file", help="Load additional run-specific variables from table (text file in csv format)", metavar="FILE")
    parser.add_argument("--log-file", help="Save submission log to specified file", metavar="FILE")
    parser.add_argument("-l", "--log", default="info", help="Sets the verbosity of log messages during job submission where LEVEL is either debug, info, warning or error", metavar="LEVEL")
    parser.add_argument("-s", "--silent", action="store_true", default=False, help="Suppress non-error (stdout) Marlin output to console")
    parser.add_argument("--dry-run", action="store_true", default=False, help="Write steering files but skip actual Marlin execution")
    parser.add_argument("--subdir", action="store_true", default=False, help="Execute every job in its own subdirectory instead of all in the base path")
    parser.add_argument("--plain", action="store_true", default=False, help="Output written to stdout/stderr and log file in prefix-less format i.e. without time stamping")
    parser.add_argument("jobtask", help="Which task to submit (e.g. convert, hitmaker, align); task names are arbitrary and can be set up by the user; they determine e.g. the config section and default steering file names.")
    parser.add_argument("runs", help="The runs to be analyzed; can be a list of single runs and/or a range, e.g. 1056-1060.", nargs='*')
    parser.add_argument("-g", "--graphic", action="store_true", default=False)
    args = parser.parse_args(argv)

    #if desired, import the colorer module
    if args.graphic:
        import Colorer

    # set the logging level
    numeric_level = getattr(logging, "INFO", None) # default: INFO messages and above
    if args.log:
        # Convert log level to upper case to allow the user to specify --log=DEBUG or --log=debug
        numeric_level = getattr(logging, args.log.upper(), None)
        if not isinstance(numeric_level, int):
            log.error('Invalid log level: %s' % args.log)
            return 2
    handler_stream.setLevel(numeric_level)
    log.setLevel(numeric_level)

    if args.plain:
        formatter = logging.Formatter('%(message)s')
        handler_stream.setFormatter(formatter)

    # set up submission log file if requested on command line
    if args.log_file:
        handler_file = logging.FileHandler([args.log_file])
        handler_file.setFormatter(formatter)
        handler_file.setLevel(numeric_level)
        log.addHandler(handler_file) 

    log.debug( "Command line arguments used: %s ", args )

    runs = list()
    for runnum in args.runs:
        try:
            log.debug("Parsing run-range argument: '%s'", runnum)
            runs = runs + parseIntegerString(runnum)
        except ValueError:
            log.error("The list of runs contains non-integer and non-range values: '%s'", runnum)
            return 2

    if not runs:
        log.error("No run numbers were specified. Please see '"+progName+" --help' for details.")
        return 2

    if len(runs) > len(set(runs)): # sets items are unique
        log.error("At least one run is specified multiple times!")
        return 2

    # dictionary keeping our parameters
    # here you can set some minimal default config values that will (possibly) be overwritten by the config file
    parameters = {"templatepath":".", "templatefile":args.jobtask+"-tmp.xml", "logpath":"."}

    # read in config file if specified on command line
    if args.conf_file:
        config = ConfigParser.SafeConfigParser()
        # local variables useful in the context of the config; access using %(NAME)s in config
        config.set("DEFAULT", "HOME",str(os.environ.get('HOME')))
        if not os.environ.get('EUTELESCOPE') is None:
            config.set("DEFAULT", "EUTelescopePath", str(os.environ.get('EUTELESCOPE')))
        else:
            log.debug("Environment variable EUTELESCOPE not found; will not be set for steering/config file parsing")
        try:
            if not config.read([args.conf_file]): # loads global defaults
                log.error("Could not read config file '%s'!", args.conf_file)
                return 1
            # merge with defaults and create final set of configuration parameters
            if config.has_section(args.jobtask):
                parameters.update(dict(config.items(args.jobtask)))
            else:
                log.warning("Config file '%s' is missing a section [%s]!", args.conf_file, args.jobtask)
            log.info("Loaded config file %s", args.conf_file)
        except ConfigParser.InterpolationMissingOptionError, err: # if interpolation during config parsing fails
            log.error('Bad value substitution in config file '+str(args.conf_file)+ ": missing '%s' key in section [%s] for option '%s'."%(err.reference, err.section, err.option))
            if err.reference == "eutelescopepath":
                log.error('EUTELESCOPE environment variable not set but required in config through "EUTelescopePath" key - please source build_env.sh in EUTelescope top directory or set variable manually.')
            return 1
    else:
        log.warn("No config file specified")
    
    # Parse option part of the  argument here -> overwriting config options
    if args.option is None:
        log.debug("Nothing to parse: No additional config options specified through command line arguments. ")
    else:
        try:
            # now parse any options given through the -o cmd line switch
            cmdoptions = dict(opt.strip().split('=', 1) for optlist in args.option for opt in optlist.split(',')) # args.option is a list of lists of strings we need to split at every '='
        except ValueError:
            log.error( "Command line error: cannot parse --option argument(s). Please use a '--option name=value' format. ")
            return 2
        for key in cmdoptions: # and overwrite our current config settings
            log.debug( "Parsing cmd line: Setting "+key+" to value '"+cmdoptions[key]+"', possibly overwriting corresponding config file option")
            parameters[key.lower()] = cmdoptions[key]

    log.debug( "Our final config:")
    for key, value in parameters.items():
        log.debug ( "     "+key+" = "+value)

    steeringTmpFileName = os.path.join(parameters["templatepath"], parameters["templatefile"])
    if not os.path.isfile(steeringTmpFileName):
        log.critical("Steering file template '"+steeringTmpFileName+"' not found!")
        return 1

    log.debug( "Opening steering file template "+steeringTmpFileName)
    steeringStringBase = open(steeringTmpFileName, "r").read()

    #Query replace steering template with our parameter set
    log.debug ("Generating base steering file")
    for key in parameters.keys():
        # check if we actually find all parameters from the config in the steering file
        try:
            # need not to search for config variables only concerning submission control
            if (not key == "templatefile" and not key == "templatepath"):
                # if using concatenation, we have a modified behavior in case the key contains "@RunRange@": then the key is replaced for every run
                if args.concatenate and parameters[key].lower().find("@runrange@")>-1:
                    log.info("Concatenation: Option '" + key + "' contains string '@RunRange@', will fill for all runs of specified range")
                    runiter = iter(runs)
                    firstrun = runiter.next() # skip first run to leave one instance of @RunNumber@ in the file for main loop
                    log.debug("Concatenation: doing substitution for first run: "+str(firstrun))
                    runkey = ireplace("@runrange@","@RunNumber@",parameters[key]) # insert run placeholder into first key, will be later replaced in main loop
                    # replace key and add the same key again for next run-through
                    steeringStringBase = ireplace("@" + key + "@", runkey+" "+"@"+key+"@", steeringStringBase) 
                    thisrun=runiter.next() 
                    for nextrun in runiter:
                        runkey = ireplace("@runrange@",str(thisrun).zfill(6),parameters[key]) # insert run number into key
                        log.debug("Concatenation: doing substitution for run "+str(thisrun)+" using " + runkey)
                        steeringStringBase = ireplace("@" + key + "@", runkey+" "+"@"+key+"@", steeringStringBase)
                        thisrun=nextrun # effectively skipping last, will need special treatment again
                    # last run: do not add key again
                    runkey = ireplace("@runrange@",str(thisrun).zfill(6),parameters[key]) # insert run number into key
                    log.debug("Concatenation: doing substitution for last run "+str(thisrun)+" using " + runkey)
                    steeringStringBase = ireplace("@" + key + "@", runkey, steeringStringBase)
                else: # the common case when not concatenating: just replace keyword
                    steeringStringBase = ireplace("@" + key + "@", parameters[key], steeringStringBase)
        except EOFError:
            if (not key == "eutelescopepath" and not key == "home" and not key == "logpath"): # do not warn about default content of config
                log.warn("Parameter '" + key + "' was not found in template file "+parameters["templatefile"])

    if args.concatenate:
        # replace list of runs with first run only
        log.info("Concatenating runs into first run")
        runs = runs[0:1] # slice run list down to first item

    # CSV table
    log.debug ("Loading csv file (if requested)")
    parameters_csv = loadparamsfromcsv(args.csv_file, runs) # store all information needed from the csv file

    # setup mechanism to deal with user pressing ctrl-c in a safe way while we execute marlin later
    import signal
    keepRunning = {'Sigint':'no'}
    def signal_handler(signal, frame):
        """ log if SIGINT detected, set variable to indicate status """
        log.critical ('You pressed Ctrl+C!')
        keepRunning['Sigint'] = 'seen'
    prevINTHandler = signal.signal(signal.SIGINT, signal_handler)

    log.info("Will now start processing the following runs: "+', '.join(map(str, runs)))
    # now loop over all runs
    for run in runs:
        if keepRunning['Sigint'] == 'seen':
            log.critical("Stopping to process remaining runs now")
            break  # if we received ctrl-c (SIGINT) we stop processing here

        runnr = str(run).zfill(6)
        log.info ("Now generating steering file for run number "+runnr+"..")

        # make a copy of the preprocessed steering file content
        steeringString = steeringStringBase

        # if we have a csv file we can parse, we will lookup the runnumber and replace any
        # variables identified by the csv header by the run specific value
        if parameters_csv:
            try:
                for field in parameters_csv[run].keys():
                    # check if we actually find all parameters from the csv file in the steering file - warn if not
                    log.debug("Parsing steering file for csv field name '%s'", field)
                    try:
                        # check that the field name is not empty and do not yet replace the runnumber
                        if not field == "" and not field == "runnumber":                    
                            steeringString = ireplace("@" + field + "@", parameters_csv[run][field], steeringString)
                    except EOFError:
                        log.warn("Parameter '" + field + "' from the csv file was not found in the template file (already overwritten by config file parameters?)")
            except KeyError:
                log.warning("Run #" + runnr + " was not found in the specified CSV file - will skip this run! ")
                continue

        try:
            steeringString = ireplace("@RunNumber@", runnr, steeringString)
        except EOFError:
            log.error("No reference to run number ('@RunNumber@') found in template file "+steeringTmpFileName)
            return 1
                
        if not checkSteer(steeringString):
            return 1

        if args.naf_file and args.lxplus_file:
            log.critical("Not possible to submit to both NAF and LXPLUS at the same time!")
            return 1

        if args.naf_file:
            args.naf_file = os.path.abspath(args.naf_file)
            if not os.path.isfile(args.naf_file):
                log.critical("NAF submission parameters file '"+args.naf_file+"' not found!")
                return 1
        elif args.lxplus_file:
            args.lxplus_file = os.path.abspath(args.lxplus_file)
            if not os.path.isfile(args.lxplus_file):
                log.critical("LXPLUS submission parameters file '"+args.lxplus_file+"' not found!")
                return 1

        log.debug ("Writing steering file for run number "+runnr)
        # When  running in subdirectories for every job, create it:
        if args.subdir:
            basedirectory = "run"+runnr
            if not os.path.exists(basedirectory):
                os.makedirs(basedirectory)

            # Decend into subdirectory:
            savedPath = os.getcwd()
            os.chdir(basedirectory)
        
        # Write the steering file:
        basefilename = args.jobtask+"-"+runnr
        steeringFile = open(basefilename+".xml", "w")

        try:
            steeringFile.write(steeringString)
        finally:
            steeringFile.close()

        # bail out if running a dry run
        if args.dry_run:
            log.info("Dry run: skipping Marlin execution. Steering file written to "+basefilename+'.xml')
        elif args.naf_file:
            rcode = submitNAF(basefilename, args.jobtask, args.naf_file, runnr) # start NAF submission
            if rcode == 0:
                log.info("NAF job submitted")
            else:
                log.error("NAF submission returned with error code "+str(rcode))
        elif args.lxplus_file:
            rcode = submitLXPLUS(basefilename, args.jobtask, args.lxplus_file, runnr) # start LXPLUS submission
            if rcode == 0:
                log.info("LXPLUS job submitted")
            else:
                log.error("LXPLUS submission returned with error code "+str(rcode))
        else:
            rcode = runMarlin(basefilename, args.jobtask, args.silent) # start Marlin execution
            if rcode == 0:
                log.info("Marlin execution done")
            else:
                log.error("Marlin returned with error code "+str(rcode))
            zipLogs(parameters["logpath"], basefilename)

        # Return to old directory:
        if args.subdir:
            os.chdir(savedPath)
        
    # return to the previous signal handler
    signal.signal(signal.SIGINT, prevINTHandler)
    if log.error.counter>0:
        log.warning("There were "+str(log.error.counter)+" error messages reported")

    return 0

if __name__ == "__main__":
    sys.exit(main())
