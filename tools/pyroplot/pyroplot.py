#!/usr/bin/env python2
"""
Python ROOT plotter - A tool for selecting and assembling histogram
plots and comparision plots from multiple ROOT files at once

Run 
python pyroplot.py --help
to see the list of command line options.

"""

# TODO:
# optimize ratio plots for 2D type histograms
# check: if compare: need at least two files!

import sys
from sys import exit # use sys.exit instead of built-in exit (latter raises exception)

import logging
import os


def timeStampCanvas ( canvas ):
    import rootpy
    from rootpy.plotting import Canvas
    from ROOT import TText
    import datetime
    # apply time stamp:
    text=TText()
    canvas.cd(0) #select canvas
    text.SetTextSize(0.015)
    text.DrawTextNDC(.7,0.008,datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))

def markCanvas ( canvas, text, x = 0.14, y = 0.007, size = 0.013, color = 1 ):
    markPad(text, x, y, size, color, canvas, 0)

def markPad ( text, x = 0.14, y = 0.007, size = 0.013, color = 1, canvas = None, npad = -1 ):
    """
    Puts a line of text on a specific canvas and pad. If none are specified, the current pad is used
    """
    import rootpy
    from rootpy.plotting import Canvas
    from ROOT import TText
    t=TText()
    if npad >= 0 and canvas:
        canvas.cd(npad) #select pad (0 = whole canvas)
    t.SetTextSize(size)
    t.SetTextColor(color)
    t.DrawTextNDC(x,y,text)

def printPlots(cans,outPath,separateFiles):
    """
    Prints a list of canvases to a single pdf file or separate files.
    """
    import rootpy
    from rootpy.plotting import Canvas
    log = logging.getLogger('pyroplot')
    if separateFiles:
        log.info( "Saving %d canvases to separate files in folder '%s'"%(len(cans),outPath))
    else:
        log.info( "Saving %d canvases to single file '%s'"%(len(cans),outPath))
    pdfOpen = False
    # loop over list of canvases and save them to file(s)
    for idx, cname in enumerate(sorted(cans)):
        can = cans[cname]
        timeStampCanvas(can)
        if not separateFiles:
            if not pdfOpen:
                log.debug( "Saving canvas #%i to '%s'"%(idx, outPath))
                # if we have multiple canvases we need a special operation 
                # to open the pdf for appending
                if len(cans) > 1:
                    can.Print(outPath+"(") # open the pdf file for appending
                else:
                    can.Print(outPath) # just write the file if handling single canvas
                pdfOpen = True
            else:
                log.debug( "Storing canvas #"+str(idx))
                # check if we are on the last canvas: then close the pdf
                if idx == (len(cans)-1):
                    log.debug( ".. and closing pdf file")
                    can.Print(outPath+")") # close the pdf file
                else:
                    can.Print(outPath) # print in the same file
                log.debug( "Done with canvas #"+str(idx))
        else:
            # set up file name to print to
            filename = os.path.join(outPath,
                                    # in case of subdirs in TFile:
                                        os.path.dirname(cname), 
                                    os.path.basename(cname)+".pdf")
            log.info( "Saving canvas #"+str(idx)+" to '" + filename +"'")
            # might have to create directory path first (if root file had subfolders)
            if not os.path.exists(os.path.dirname(filename)):
                        os.makedirs(os.path.dirname(filename))
            can.Print(filename)


def makePlotCollection( histodicts, files, filesdescr, outPath, docompare, doseparateFiles, logscale):
    log = logging.getLogger('pyroplot')
    canvas = {}
    if docompare:
        canvas.update( makeComparisionPage(histodicts, files, filesdescr, doseparateFiles))
    else:
        for idx, histos in enumerate(histodicts):
            canvas.update( makePage(histos, files[idx], filesdescr[files[idx]], doseparateFiles, logscale))
    # now print the canvases
    printPlots(canvas,outPath,doseparateFiles)

def plotHistos ( histos, text = "", option = "", statbox = True):
    """ 
    Plots a list of histograms
    """
    log = logging.getLogger('pyroplot')
    import rootpy
    from rootpy.plotting import Hist, HistStack
    from ROOT import kRed,gPad,TPaveStats
    # only for 1D and 2D type histograms:
    # determine and set maximum for all histograms
    stack = HistStack()
    for hist in histos:
        if not hist.__class__.__name__ == "Hist3D":
            stack.Add(hist)
    # determine maximum value and set it for all histograms
    if stack.GetHists():
        maxplot = stack.GetMaximum()
        minplot = stack.GetMinimum()
        for hist in stack.GetHists():
            hist.SetMaximum(maxplot)
            # special treatment for log scale Y
            if gPad.GetLogy():
                hist.SetMinimum(1.)
            else:
                # if histogram minimum is positive, set to 0.
                if minplot > 0:
                    hist.SetMinimum(0.)
    for idx, hist in enumerate(histos):
        try:
            thisopt = option
            # if not first histo, add "same" to options so previous ones are not overwritten
            if idx:
                if statbox:
                    thisopt += "sames"
                else:
                    thisopt += "same"
            histcopy = hist.DrawCopy(thisopt)
            # on first go: print identifying text on pad
            if not idx and text:
                markPad(text=text, x=.14, y=.85, size=0.041)
                if statbox:
                    thisopt += "sames"
                else:
                    thisopt += "same"
                histcopy.Draw(thisopt)
            gPad.Update()
            if statbox:
                try:
                    statBox = histcopy.GetListOfFunctions().FindObject("stats")
                    # offset from last statbox
                    offset = .18
                    # needs to be larger for Profile & 2D histos
                    if (hist.__class__.__name__=="Profile"
                        or hist.__class__.__name__=="Hist2D"
                        or hist.__class__.__name__=="Profile2D"
                        or hist.__class__.__name__=="Hist3D"):
                        offset = .26
                    statBox.SetY1NDC(statBox.GetY1NDC()-offset*(idx))
                    statBox.SetY2NDC(statBox.GetY2NDC()-offset*(idx))
                    statBox.SetTextColor(hist.GetLineColor())
                    statBox.SetBorderSize(2)
                except AttributeError:
                    log.debug("Could not get statbox for histogram " + hist.GetName())
            
        except rootpy.ROOTError, e:
            log.error("Drawing histogram %s caused ROOTError exception ('%s')"%(hist.GetName(),e.msg))
            gPad.Clear() # otherwise this happens again when drawing..
            markPad ( text="Could not draw %s ('%s')"%(hist.GetName(),e.msg), x = 0.14, y = 0.4, size = 0.023, color = kRed)
            return # give up!


def makePage( histos, fileName, fileDescr, separateFiles, logscale):
    """
    Prepares a canvas with one histogram per pad
    """
    import rootpy
    from rootpy.plotting import Hist, Canvas
    from ROOT import kBlue,gPad
    log = logging.getLogger('pyroplot')
    cans = {}
    log.info( "Drawing histograms .." )
    for idx, name in enumerate(sorted(histos.keys())):
        if separateFiles:
            log.debug( "Creating new canvas with index %d."%(idx))
            c=Canvas( 600, 800)
            cans[name]=c
            markCanvas(c, fileName, 0.05, y = 0.009, size = 0.025, color = 2 )
        if not separateFiles and (idx)%6 == 0:
            log.debug( "Creating new canvas with index %d."%(idx/6))
            # start a new canvas
            c=Canvas( 600, 800)
            cans[fileName+'/'+name]=c
            c.Divide(2,3)
            markCanvas(c, fileName, 0.05, y = 0.009, size = 0.025, color = 2 )
        # draw the histogram
        hist = histos[name]
        log.debug( "Drawing histogram #" + str(idx%6+1) +": " + hist.GetName() + " (" + hist.__class__.__name__ + ") in canvas #" + str(int(idx/6) ))
        hist.color = kBlue
        if not separateFiles:
            c.cd(idx%6+1)
        if logscale:
            gPad.SetLogy()
        plotHistos(histos=[hist],text=name)
        if fileDescr:
            markPad(text=fileDescr, x=.14, y=.8, size=0.041, color = 2)
    return cans


def makeComparisionPage( histodicts , fileNames, fileDescr, separateFiles):
    """
    Prepares a canvas comparing multiple histograms: plotting all in one pad and their ratios in the second
    """
    import rootpy
    from rootpy.plotting import Hist, Canvas, Legend
    import ROOT
    from ROOT import gPad
    log = logging.getLogger('pyroplot')
    cans = {}
    colors = [ROOT.kBlue, ROOT.kRed+1,ROOT.kViolet-1, ROOT.kOrange+7,ROOT.kGreen-7,ROOT.kOrange-6,
              ROOT.kPink-9,ROOT.kTeal-6,ROOT.kBlue+4,ROOT.kAzure+2]
    log.info( "Drawing histograms .." )
    # prepare set of histograms to compare to the reference on (the first)
    # loop over the reference set of histos (sorted by key):
    for hidx, refname in enumerate(sorted(histodicts[0].keys())):
        # prepare canvas
        if separateFiles:
            log.debug( "Creating new canvas with index %d."%(hidx))
            c=Canvas( 600, 270)
            cans[refname] = c
            c.Divide(3,1)
            c.cd(1)
        if not separateFiles and (hidx)%4 == 0:
            log.debug( "Creating new canvas with index %d."%(hidx/3))
            # start a new canvas
            c=Canvas( 600, 800)
            cans[refname] = c
            c.Divide(3,4)
        # prepare histograms for drawing
        log.debug( "Drawing histogram #" + str(hidx+1) +" (" + refname + ") on canvas #" + str(len(cans)) )
        hists = []
        ratiohists = []
        hiter = iter (histodicts)
        # special treatment for tprofile: prepare the reference projection for the ratio
        if histodicts[0][refname].__class__.__name__=="Profile":
            refProj = histodicts[0][refname].ProjectionX()
            refProj.SetName("reference_proj")
        for idx, h in enumerate(hiter):
            # make sure we have this histogram loaded:
            if not refname in h:
                continue
            # access the corresponding histogram of the other files at the same hidx as used for ref
            h[refname].color = colors[idx]
            h[refname].linestyle = idx
            hists.append (h[refname])
            # prepare ratio is this is not the first (i.e. reference) histogram
            if idx:
                # special treatment for TProfile:
                if h[refname].__class__.__name__=="Profile":
                    myratio = Hist(h[refname].nbins(), h[refname].lowerbound(), h[refname].upperbound()) #dummy hist
                    myratioproj = h[refname].ProjectionX()
                    myratioproj.SetName("cmp_hist_proj"+str(idx))
                    try:
                        myratio.divide(myratioproj,refProj)
                    except rootpy.ROOTError, e:
                        log.error("Calculation of ratio for histogram %s caused ROOTError exception ('%s')"%(h[refname].GetName(),e.msg))
                        break
                    myratio.color = colors[idx]
                else:
                    myratio = h[refname].clone() # make sure that the ratio has the right type
                    try:
                        myratio.Divide(h[refname], histodicts[0][refname]) # divide by reference hist
                    except rootpy.ROOTError, e:
                        log.error("Calculation of ratio for histogram %s caused ROOTError exception ('%s')"%(h[refname].GetName(),e.msg))
                        break
                myratio.yaxis.SetTitle("(h_{cmp} - h_{ref})/h_{ref}")
                myratio.SetTitle("ratio to reference")
                myratio.SetMaximum(2)
                myratio.SetMinimum(0)
                myratio.SetStats(0)
                ratiohists.append(myratio)
        # go to right pad to draw main histos
        if not separateFiles:
            c.cd(3*(hidx%4)+1) # 12 pads (4x3) but we always skip one for the ratio and one for the log plot
        plotHistos(histos=hists,text=refname)
        # go to right pad to draw log plots
        if not separateFiles:
            c.cd(3*(hidx%4)+2)
        else:
            c.cd(2)
        gPad.SetLogy()
        plotHistos(histos=hists, text="log scale")
        # go to last pad to draw ratios
        if not separateFiles:
            c.cd(3*(hidx%4)+3)
        else:
            c.cd(3)
        log.debug( "Plotting ratios for #" + str(hidx%3+1) +": " + refname + " in canvas #" + str(int(hidx/3) ))
        plotHistos(histos=ratiohists, statbox=False, text="ratio")
        #create legend
        legend = Legend(nentries=len(histodicts)+len(fileDescr)+1, leftmargin=0.25, 
                        topmargin=0.05, rightmargin=0.25, entryheight=0.05)
        legend.AddEntry(histodicts[0][refname],label=os.path.basename(fileNames[0]), style="l")
        if fileNames[0] in fileDescr:
            if fileDescr[fileNames[0]]:
                legend.AddEntry(None,label=fileDescr[fileNames[0]],style="")
        legend.AddEntry(None,label="(reference)",style="")
        for idx,ratiohist in enumerate(ratiohists):
            legend.AddEntry(ratiohist,
                            label=os.path.basename(fileNames[idx+1]), 
                            style="l")
            # add additional info if specified by the user
            if fileNames[idx+1] in fileDescr:
                if fileDescr[fileNames[idx+1]]:
                    legend.AddEntry(None,label=fileDescr[fileNames[idx+1]],style="")
        legend.SetBorderSize(0)
        legend.SetMargin(0.3)
        legend.Draw()
        # we need to keep the legend around in memory or it will be deleted before the plot is drawn
        ROOT.SetOwnership( legend, False ) # True to own, False to release
    return cans

def findHistogramsInFile(filename, regexlist, strict, verbose = False):
    """ searches the ROOT file given by filename and returns a list of files matching 
    one of the regular expressions in regexlist
    """
    log = logging.getLogger('pyroplot')
    import re
    import rootpy
    from rootpy.io import root_open
    # speed up search routine by explicitly caching the regular expressions
    cRegEx = [] # stores the compiled reg ex
    for thisSelection in regexlist:
        if strict:
            # compile with implied beginning and ending of line markers (^ and $)
            cRegEx.append( re.compile("^%s$"%thisSelection) )
        else:
            # if only searching for partial matches: do not need ^ or $
            cRegEx.append( re.compile(thisSelection) )
    selectedHistos = []
    nobj = 0
    if verbose:
        # some header output only needed when listing the contents
        print "=================================="
        print " File: %s"%(filename)
        print " (matched by regex) path/and/objectname"
        print "=================================="
    # open file and loop over it
    f = root_open(filename)
    try:
        # recursively walk through the file
        for path, dirs, objects in f.walk():
            nobj+=len(objects) # sum up every object we encounter
            for obj in objects:
                fullpath = path + '/' + obj
                # try to find a match in our list of compiled regular expressions
                foundMatch = False
                for idx, thisRegEx in enumerate(cRegEx):
                    log.debug("Matching regex %s against %s"%(regexlist[idx],fullpath))
                    matches = False
                    if strict: 
                        matches = thisRegEx.match(fullpath)
                    else:
                        matches = thisRegEx.search(fullpath)
                    if matches:
                        foundMatch = True
                if foundMatch:
                    selectedHistos.append( fullpath )
                    log.debug('Found match: %s '%(fullpath))
                if verbose:
                    if foundMatch:
                        sys.stdout.write("(*) ")
                    else:
                        sys.stdout.write("( ) ")
                    print fullpath
    finally:
        f.close()
    log.info("%s: %d/%d objects selected."%(filename, len(selectedHistos),nobj))
    return selectedHistos

def loadHistogramsFromFile(filename, histonames, with2d, with3d):
    """ 
    loads specified histograms from the ROOT file given by filename and returns them.
    The histograms will no longer be associated with the file.
    """
    log = logging.getLogger('pyroplot')
    import rootpy
    from rootpy.io import root_open
    from rootpy.plotting import Hist
    histos = {}
    f =  root_open(filename);
    nignored = 0
    for h in histonames:
        try:
            histo = f.Get(h)
        except rootpy.io.DoesNotExist:
            # this can happen if the reference file contains more histos than the others
            log.warn("%s not found in file %s"%(h,filename))
            continue
        # might want to ignore Hist2D etc: VERY SLOW and plot processing not optimally suited yet
        if ((histo.__class__.__name__=="Hist" 
             or histo.__class__.__name__=="Profile")
            or ((with3d or with2d) and histo.__class__.__name__=="Hist2D")
            or (with3d and (histo.__class__.__name__=="Profile2D" 
                or histo.__class__.__name__=="Hist3D"))):
            histo.SetDirectory(0) # remove association with file
            histos[h] = histo
        else:
            log.debug("IGNORING %s as it is of class '%s'"%(h,histo.__class__.__name__))
            nignored += 1
    f.close()
    log.info("Loaded %d histograms from file %s"%(len(histos),filename))
    if nignored:
        log.info("IGNORED %d matching 2D/3D histograms: to see these use the --with-2D or --with-3D switches."%(nignored))
    return histos

def run( argv = sys.argv ):
    # Must be done before :py:mod:`rootpy` logs any messages.
    import logging;
    log = logging.getLogger('pyroplot') # set up logging

    try:
        import ROOT
    except ImportError:
        # module failed to load - maybe PYTHONPATH is not set correctly?
        # guess the right path, but that is only possible if ROOTSYS is set:
        if os.environ.get('ROOTSYS') is None:
            print "ERROR: Could not load the Python ROOT module. Please make sure that your ROOT installation is compiled with Python support and that your PYTHONPATH is set correctly and includes libPyROOT.so"
            exit(1)
        sys.path.append(os.path.join(os.environ.get('ROOTSYS'),"lib"))
        sys.path.append(os.path.join(os.environ.get('ROOTSYS'),"lib","root"))
        # try again:
        try:
            import ROOT
        except ImportError:
            print "ERROR: Could not load the Python ROOT module. Please make sure that your ROOT installation is compiled with Python support and that your PYTHONPATH is set correctly and includes libPyROOT.so"
            exit(1)

    try:
        import rootpy
    except ImportError:
        # rootpy is not installed; use (old) version provided with EUTelescope
        # determine (real) path to subdirectory pymodules (relative to current path)
        libdir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(__file__))),"pymodules","rootpy")
        # search for any rootpy folders
        import glob
        rootpydirs = glob.glob(libdir+"*")
        if not rootpydirs:
            print "Error: Could not find the rootpy module provided with EUTelescope in %s!"%(libdir)
        else:
            # add last entry to python search path (subfolder rootpy where the modules are located)
            sys.path.append(rootpydirs[-1])
        # try again loading the module
        try:
            import rootpy
        except ImportError:
            print "Error: Could not load the rootpy modules. Please install them from http://www.rootpy.org/install.html"
            exit(1)
        except SyntaxError:
            req_version = (2,5)
            cur_version = sys.version_info
            if cur_version < req_version:
                print "Error: Python version too old: due to its dependency on rootpy, this script requires a Python interpreter version 2.6 or later (installed: %s.%s.%s)!"%(cur_version[:3])
                exit(1)
            print "Error: Failed to load rootpy module! Possibly incompatible with installed Python version (%s.%s.%s)?"%(cur_version[:3])
            exit(1)

    from rootpy import log; log = log["/pyroplot"]
    rootpy.log.basic_config_colorized()
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = 1001

    import argparse
    # command line argument parsing
    parser = argparse.ArgumentParser(description="Python ROOT plotter - A tool for selecting and assembling histogram plots and comparision plots from multiple ROOT files at once")
    parser.add_argument('--version', action='version', version='Revision: $Revision$, $LastChangedDate$')
    parser.add_argument("-l", "--log-level", default="info", help="Sets the verbosity of log messages where LEVEL is either debug, info, warning or error", metavar="LEVEL")
    parser.add_argument("--compare", action="store_true", default=False, help="Compare the selected histograms between files (ratio plots, chi2) where the first file provides the reference.")
    parser.add_argument("-log", "--log-scale", action="store_true", default=False, help="Uses a logarithmic scale for the y axis; only relevant when not using '--compare'.")
    parser.add_argument('--select', '-s', action='append', help="Specify regular expression(s) for histogram selection.")
    parser.add_argument("--selection-from-file", help="Load list of regular expressions for histogram selection from file (plain text file, one reg ex per line).", metavar="FILE")
    parser.add_argument("--one-file-per-histogram", action="store_true", default=False, help="Writes one file per histogram instead of storing all plots in one single file.")
    parser.add_argument("-o","--output", default="./overview.pdf", help="Output path and file name. If the file does not end in '.pdf' it will be assumed to be a path and created if needed. If --one-file-per-histogram is set, this will be the output directory for the plots.", metavar="FILE/PATH")
    parser.add_argument("--with-2D","-2D", action="store_true", default=False, help="Also loads TH2-type histograms.")    
    parser.add_argument("--with-3D","-3D", action="store_true", default=False, help="Also loads TH3-type and Profile2D-type histograms, implies --with-2D.")    
    parser.add_argument("--list-only", "--list", action="store_true", default=False, help="Do not generate plots but only list objects in ROOT file(s) and indicate which ones would be selected.")
    parser.add_argument("--strict", action="store_true", default=False, help="Require the selection to match the full histogram path and name (with implied '^' and '$') instead of only a partial match.")
    parser.add_argument("files", help="The files to be processed; additional info STRING to be included in the plot legend can be added by specifiying FILE:STRING", nargs='+')
    # parse the arguments
    args = parser.parse_args(argv)
    # set the logging level
    numeric_level = getattr(logging, "INFO", None) # default: INFO messages and above
    if args.log_level:
        # Convert log level to upper case to allow the user to specify --log-level=DEBUG or --log-level=debug
        numeric_level = getattr(logging, args.log_level.upper(), None)
        if not isinstance(numeric_level, int):
            log.error('Invalid log level: %s' % args.log_level)
            exit(2)
    log.setLevel(numeric_level)
    log.debug( "Command line arguments used: %s ", args )

    log.debug("Using rootpy %s from %s"%(rootpy.__version__,rootpy.__file__))

    # laod and combine all specified reg ex
    regexs = []
    # first from file
    if args.selection_from_file:
        f = open(args.selection_from_file, 'r')
        try:
            lines = f.read().splitlines()
            for line in lines:
                if line: # test if line is not empty (would match anything)
                    log.debug("Loading reg ex from file " + args.selection_from_file 
                              + ": '" + line +"'")
                    regexs.append(line)
        finally:
            f.close()
    if args.select:
        for arg in args.select:
            log.debug("Using reg ex from command line: " + arg)
            regexs.append(arg)
    # still nothing to select? use default
    if not regexs:
        import inspect
        filepath = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(__file__))),"default.sel")
        try:
            f = open(filepath, 'r')
            try:
                lines = f.read().splitlines()
                for line in lines:
                    if line: # test if line is not empty (would match anything)
                        log.debug("Loading reg ex from file " + filepath + ": '" + line +"'")
                        regexs.append(line)
            finally:
                f.close()
        except IOError:
            log.warn("Could not find the file with the default selection ('"+filepath+"'), will use default of '.*' (select all)")
            regexs.append('.*')

    # parse output file name and verify that it ends in '.pdf'
    outputFilePath = ""
    fileName, fileExtension = os.path.splitext(args.output)
    if not fileExtension == '.pdf':
        log.debug("Output argument does not end in '.pdf': '%s'. Assuming it's meant to be a path"%args.output)
        if not args.one_file_per_histogram:
            # append default name for single histogram file
            outputFilePath = os.path.join(args.output,"overview.pdf")
        else:
            outputFilePath = args.output
    else:
        if args.one_file_per_histogram:
            # all we need is the path, strip the file and append a slash
            outputFilePath = os.path.dirname(args.output)
        else:
            outputFilePath = args.output

    # parse file names and extract additionally provided info
    fileNames = []
    fileDescr = {}
    for thisFile in args.files:
        s = thisFile.strip().split(':', 1) # try to split the string
        if (len(s)==1):
            # didn't work, only have one entry
            fileNames.append(s[0])
            fileDescr[s[0]] = ""
        else:
            fileNames.append(s[0])
            fileDescr[s[0]] = s[1]
                         
    histoDicts = [] # our histograms: each element will store a dict()
                    # of histogram objects with its full path in the
                    # root file as key
    selectedHistos = []
    # loop over all files
    for idx, thisFile in enumerate( fileNames ):
        # only search for matching histo names on first iteration if doing a comparison between files
        if not idx or not args.compare:
            selectedHistos = findHistogramsInFile(thisFile, regexs, args.strict, args.list_only)
        if args.list_only:
            continue
        h = loadHistogramsFromFile( thisFile , selectedHistos, args.with_2D, args.with_3D)
        histoDicts.append(h) # append to main histo list
    if histoDicts:
        log.info("Input file(s) read. %d histograms matched selection criteria and were loaded"%(sum(len(histos) for histos in histoDicts)))
        makePlotCollection(histoDicts, fileNames, fileDescr, outputFilePath, args.compare, args.one_file_per_histogram,args.log_scale)

    log.info("done")    
    
if __name__ == "__main__":
    run( sys.argv[1:])
