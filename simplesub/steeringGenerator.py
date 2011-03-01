import subprocess

def parseArgs(args):
    #remove program name
    import sys
    progName = args.pop(0)
    dryRun = False
    runs = []
    try:
        argindex = args.index("-a") + 1
    except:
        print("You must supply a mode with -a. Try \'python", progName, "-a help\' for a list of modes.")
        sys.exit(1)
    mode = args[ argindex ]
    args.pop( argindex )
    args.pop( argindex - 1 )
    if( args.count("-d") > 0 ):
        args.pop( args.index("-d") )
        dryRun = True
    modes = []
    if( mode == "all"): modes = ["rawtohit", "aligntel", "aligndut", "aligndut2", "fitter"]
    else: modes = [mode]
    try:
        for runnum in args:
            runs.append(int(runnum))
    except:
        print("All these values are not numbers:", args, "\nThey should be.")
        sys.exit(1)
    return modes, args, dryRun

def checkSteer(sstring):
    import re, sys
    m = re.search("@.*@", sstring)
    if( m != None):
        print ("Missing param with key: " , m.group(0))
        sys.exit(1)

def generateSteeringFile(run, alignruns, fun, siteConfig, runConfig, jobs):
    params = {}
    params["RunNumber"] = run
    params["AlignInputFiles"] = alignruns
    siteConfig(params)
    runConfig(params)
    fun(params)
    #Slurp template file 
    steeringString = open(params["TemplatePath"] + "/" + params["TemplateFile"], "r").read()
    #Query replace
    for key in params.keys():
        steeringString = steeringString.replace("@" + key + "@", params[key])
    checkSteer(steeringString)
    jobs.append( [run, steeringString] )

def startJob(mode, job, procs):
    fname = mode + "-steer-" + job[0] + ".xml"
    steeringOut = open( fname , "w")
    steeringOut.write(job[1])
    steeringOut.close()
    fname2 = mode + "-log-" + job[0] + ".txt"
    logStream = open(fname2, "w")
    print("starting job: Marlin " + fname)
    job = subprocess.Popen(["Marlin", fname], stdout=logStream, stderr=logStream )
    procs[job.pid] = [job, fname, fname2]

def cleanUpJob(procl):
    import os
    log = open( procl[2], "r").read()
    print(log)
    os.remove(procl[1])
    os.remove(procl[2])

def jobMaker( functions, runConfig):
    import sys, time, os
    from siteConfig import siteConfig 
    modes, runs, dryRun = parseArgs(sys.argv)
    globParams = {}
    siteConfig(globParams)
    runConfig(globParams)
    if( len(runs) == 0):
        runs = eval(globParams["DataSet"])
    for mode in modes:
        jobs = []
        func = None
        try:
            func = functions[mode]
        except:
            print( mode , " does not appear to be a valid mode. Should be one of ", functions.keys())
            sys.exit()
        if (mode == "aligndut") or (mode == "aligntel") or (mode == "daffittersingle"):
            runstring = ""
            for run in runs:
                runstring = runstring + globParams["ResultPath"] + "/run" + str(run).zfill(6) + "-hit-p.slcio\n"
            generateSteeringFile(str(runs[0]).zfill(6), runstring, func, siteConfig, runConfig, jobs)
        else:
            for run in runs:
                generateSteeringFile(str(run).zfill(6), "", functions[mode], siteConfig, runConfig, jobs )
        procs = {}
        if dryRun:
            for steer in jobs:
                print( steer[1])
        else:
            while len(procs) or len(jobs):
                if (len(procs) >= int(globParams["NParallel"])) or len(jobs) == 0:
                    print( len(procs), "job(s) running, waiting for one to finish...")
                    (pid,status) = os.wait()
                    cleanUpJob(procs[pid])
                    del procs[pid]
                if len(jobs):
                    startJob(mode, jobs[0], procs)
                    del jobs[0]
        print( "Done")
