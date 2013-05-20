#!/usr/bin/env python
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

def markCanvas ( canvas, textstring, x = 0.14, y = 0.007, size = 0.013, color = 1 ):
    markPad(textstring, x, y, size, color, canvas, 0)

def markPad ( textstring, x = 0.14, y = 0.007, size = 0.013, color = 1, canvas = None, npad = -1 ):
    import rootpy
    from rootpy.plotting import Canvas
    from ROOT import TText
    text=TText()
    if npad >= 0 and canvas:
        canvas.cd(npad) #select pad (0 = whole canvas)
    text.SetTextSize(size)
    text.SetTextColor(color)
    text.DrawTextNDC(x,y,textstring)

def printPlots(cans,outPath,separateFiles):
    import rootpy
    from rootpy.plotting import Canvas
    log = logging.getLogger('pyroplot')
    if separateFiles:
        log.info( "Saving %d canvases to separate files in folder '%s'"%(len(cans),outPath))
    else:
        log.info( "Saving %d canvases to single file '%s'"%(len(cans),outPath))
    pdfOpen = False
    # loop over list of canvases and save them to file(s)
    for idx, cname in enumerate(cans):
        can = cans[cname]
        timeStampCanvas(can)
        if not separateFiles:
            if not pdfOpen:
                log.debug( "Saving canvas to '" + outPath +"'")
                can.Print(outPath+"(") # open the pdf file for appending
                pdfOpen = True
            else:
                log.debug( "Storing canvas #"+str(idx))
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
    if pdfOpen:
        can.Print(outPath+")") # close the pdf file


def makePlotCollection( histodicts, files, filesdescr, outPath, docompare, doseparateFiles):
    log = logging.getLogger('pyroplot')
    canvas = {}
    if docompare:
        canvas.update( makeComparisionPage(histodicts, files, filesdescr, doseparateFiles))
    else:
        for idx, histos in enumerate(histodicts):
            canvas.update( makePage(histos, files[idx], doseparateFiles))
    # now print the canvases
    printPlots(canvas,outPath,doseparateFiles)


def makePage( histos, fileName, fileDescr, separateFiles):
    import rootpy
    from rootpy.plotting import Hist, Canvas
    from ROOT import kBlue,kRed,kGreen,kYellow,gPad,TPaveStats
    log = logging.getLogger('pyroplot')
    cans = {}
    for idx, name in enumerate(histos):
        if separateFiles:
            log.debug( "Creating new canvas with index %d."%(idx))
            c=Canvas( 600, 800)
            cans[name]=c
            markCanvas(c, fileName+" "+fileDescr, 0.05, y = 0.009, size = 0.025, color = 2 )
        if not separateFiles and (idx)%6 == 0:
            log.debug( "Creating new canvas with index %d."%(idx/6))
            # start a new canvas
            c=Canvas( 600, 800)
            cans[name]=c
            c.Divide(2,3)
            markCanvas(c, fileName+" "+fileDescr, 0.05, y = 0.009, size = 0.025, color = 2 )
        # draw the histogram
        hist = histos[name]
        log.debug( "Drawing histogram #" + str(idx%6+1) +": " + hist.GetName() + " in canvas #" + str(int(idx/6) ))
        hist.color = kBlue
        if not separateFiles:
            c.cd(idx%6+1)
        h = hist.DrawCopy("")
        gPad.Update()
        statBox = h.GetListOfFunctions().FindObject("stats")
        statBox.SetName('new_stat')
        statBox.SetTextColor(kBlue)
        statBox.SetBorderSize(2)
    return cans


def makeComparisionPage( histodicts , fileNames, fileDescr, separateFiles):
    import rootpy
    from rootpy.plotting import Hist, Canvas, Legend
    import ROOT
    from ROOT import gPad,TPaveStats
    log = logging.getLogger('pyroplot')
    cans = {}
    colors = [ROOT.kBlue, ROOT.kRed+1,ROOT.kViolet-1, ROOT.kOrange+7,ROOT.kGreen-7,ROOT.kOrange-6,
              ROOT.kPink-9,ROOT.kTeal-6,ROOT.kBlue+4,ROOT.kAzure+2]
    # prepare set of histograms to compare to the reference on (the first)
    # loop over the reference set of histos:
    for hidx, refname in enumerate(histodicts[0]):
        # prepare canvas
        if separateFiles:
            log.debug( "Creating new canvas with index %d."%(hidx))
            c=Canvas( 600, 270)
            cans[refname] = c
            c.Divide(2,1)
            c.cd(1)
        if not separateFiles and (hidx)%3 == 0:
            log.debug( "Creating new canvas with index %d."%(hidx/3))
            # start a new canvas
            c=Canvas( 600, 800)
            cans[refname] = c
            c.Divide(2,3)
        # draw the reference histogram
        refhist = histodicts[0][refname]
        log.info( "Drawing reference histogram #" + str(hidx%3+1) +": " + refhist.GetName() + " in canvas #" + str(int(hidx/3) ))
        refhist.color = colors[0]
        # go to right pad
        if not separateFiles:
            c.cd(2*(hidx%3)+1) # six pads (2x3) but we always skip one for the ratio
        refcopy = refhist.DrawCopy("")
        gPad.Update()
        statBox = refcopy.GetListOfFunctions().FindObject("stats")
        statBox.SetTextColor(colors[0])
        statBox.SetBorderSize(2)
        # draw remaining histograms and prepare ratio plots
        # prepare histos for residual/ratio plots
        ratiohists = []
        riter = iter (histodicts)
        next (riter) # skip first run for ratio since it's the reference histo
        for idx, r in enumerate(riter):
            # make sure we have this histogram loaded:
            if not r[refname]:
                continue
            # access the corresponding histogram of the other files at the same hidx as used for ref
            r[refname].color = colors[idx+1]
            r[refname].linestyle = idx
            h = r[refname].DrawCopy("sames")
            gPad.Update()
            statBox = h.GetListOfFunctions().FindObject("stats")
            statBox.SetY1NDC(statBox.GetY1NDC()-.18*(idx+1))
            statBox.SetY2NDC(statBox.GetY2NDC()-.18*(idx+1))
            statBox.SetTextColor(colors[idx+1])
            statBox.SetBorderSize(2)
            myratio = r[refname].Clone()
            myratio.Add(refhist,-1.) # subtract reference hist
            myratio.Divide(refhist)
            myratio.yaxis.SetTitle("(h_{cmp} - h_{ref})/h_{ref}")
            myratio.SetTitle("relative difference to reference")
            myratio.SetMaximum(1)
            myratio.SetMinimum(-1)
            myratio.SetStats(0)
            ratiohists.append(myratio)
        if not separateFiles:
            c.cd(2*(hidx%3)+2)
        else:
            c.cd(2)
        for idx,ratiohist in enumerate(ratiohists):
            arg = "same"
            if not idx:
                arg = ""
            ratiohist.DrawCopy(arg)
        #create legend
        legend = Legend(nentries=len(histodicts)+len(fileDescr)+1, leftmargin=0.25, 
                        topmargin=0.05, rightmargin=0.25, entryheight=0.05)
        legend.AddEntry(refhist,label=os.path.basename(fileNames[0]), legendstyle="l")
        if fileNames[0] in fileDescr:
            legend.AddEntry(None,label=fileDescr[fileNames[0]],legendstyle="")
        legend.AddEntry(None,label="(reference)",legendstyle="")
        for idx,ratiohist in enumerate(ratiohists):
            legend.AddEntry(ratiohist,
                            label=os.path.basename(fileNames[idx+1]), 
                            legendstyle="l")
            # add additional info if specified by the user
            if fileNames[idx+1] in fileDescr:
                legend.AddEntry(None,label=fileDescr[fileNames[idx+1]],legendstyle="")
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
        print "=================================="
    # open file and loop over it
    with root_open(filename) as f:
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
    for h in histonames:
        try:
            histo = f.Get(h)
        except rootpy.io.DoesNotExist:
            # this can happen if the reference file contains more histos than the others
            log.warn("%s not found in file %s"%(h,filename))
            continue
        if issubclass(type(histo), Hist):
            print "subclass checks out for %s: "%(h,type(histo))
        # might want to ignore Hist2D etc: VERY SLOW and plot processing not optimally suited yet
        if ((histo.__class__.__name__=="Hist" 
             or histo.__class__.__name__=="Profile")
            or ((with3d or with2d) and histo.__class__.__name__=="Hist2D")
            or (with3d and histo.__class__.__name__=="Profile2D" 
                or histo.__class__.__name__=="Hist3D")):
            log.debug("%s inherits from Hist/TH1"%h)
            histo.SetDirectory(0) # remove association with file
            histos[h] = histo
        else:
            log.warn("IGNORING %s as it is of class '%s'. Enable with --with-2D or --with-3D"%(h,histo.__class__.__name__))
    f.close()
    log.info("Loaded %d histograms from file %s"%(len(histos),filename))
    return histos

def run( argv = sys.argv ):
    # Must be done before :py:mod:`rootpy` logs any messages.
    import logging;
    log = logging.getLogger('pyroplot') # set up logging
    try:
        import rootpy
        from rootpy import log; log = log["/pyroplot"]
        rootpy.log.basic_config_colorized()
    except ImportError:
        log.error("Could not load the rootpy modules. Please install from http://www.rootpy.org/install.html")
        exit(1)
    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        ROOT.gErrorIgnoreLevel = 1001
    except ImportError:
        log.error("Could not load the Python ROOT module. Please make sure that your ROOT installation is correctly set up and includes libPyROOT.so")
        exit(1)
    import argparse
    # command line argument parsing
    parser = argparse.ArgumentParser(description="Python ROOT plotter - A tool for selecting and assembling histogram plots and comparision plots from multiple ROOT files at once")
    parser.add_argument('--version', action='version', version='Revision: $Revision$, $LastChangedDate$')
    parser.add_argument("-l", "--log", default="info", help="Sets the verbosity of log messages where LEVEL is either debug, info, warning or error", metavar="LEVEL")
    parser.add_argument("--compare", action="store_true", default=False, help="Compare the selected histograms between files (ratio plots, chi2) where the first file provides the reference.")
    parser.add_argument('--select', '-s', action='append', help="Specify regular expression(s) for histogram selection.")
    parser.add_argument("--selection-from-file", help="Load list of regular expressions for histogram selection from file (plain text file, one reg ex per line).", metavar="FILE")
    parser.add_argument("--one-file-per-histogram", action="store_true", default=False, help="Writes one file per histogram instead of storing all plots in one single file.")
    parser.add_argument("-o","--output", default="./overview.pdf", help="Output path and file name. If the file does not end in '.pdf' it will be assumed to be a path and created if needed. If --one-file-per-histogram is set, this will be the output directory for the plots.", metavar="FILE/PATH")
    parser.add_argument("--with-2D","-2D", action="store_true", default=False, help="Also loads TH2-type histograms.")    
    parser.add_argument("--with-3D","-3D", action="store_true", default=False, help="Also loads TH3-type and Profile2D-type histograms, implies --with-2D.")    
    parser.add_argument("--list-only", action="store_true", default=False, help="Do not generate plots but only list objects in ROOT file(s) and indicate which ones would be selected.")
    parser.add_argument("--strict", action="store_true", default=False, help="Require the selection to match the full histogram path and name (with implied '^' and '$') instead of only a partial match.")
    parser.add_argument("files", help="The files to be processed; additional info STRING to be included in the plot legend can be added by specifiying FILE:STRING", nargs='+')
    # parse the arguments
    args = parser.parse_args(argv)
    # set the logging level
    numeric_level = getattr(logging, "INFO", None) # default: INFO messages and above
    if args.log:
        # Convert log level to upper case to allow the user to specify --log=DEBUG or --log=debug
        numeric_level = getattr(logging, args.log.upper(), None)
        if not isinstance(numeric_level, int):
            log.error('Invalid log level: %s' % args.log)
            exit(2)
    log.setLevel(numeric_level)
    log.debug( "Command line arguments used: %s ", args )

    # laod and combine all specified reg ex
    regexs = []
    # first from file
    if args.selection_from_file:
        with open(args.selection_from_file, 'r') as f:
            lines = f.read().splitlines()
            for line in lines:
                if line: # test if line is not empty (would match anything)
                    log.debug("Loading reg ex from file " + args.selection_from_file 
                              + ": '" + line +"'")
                    regexs.append(line)
    if args.select:
        for arg in args.select:
            log.debug("Using reg ex from command line: " + arg)
            regexs.append(arg)
    # still nothing to select? use default
    if not regexs:
        import inspect
        filepath = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"default.sel")))
        try:
            with open(filepath, 'r') as f:
                lines = f.read().splitlines()
                for line in lines:
                    if line: # test if line is not empty (would match anything)
                        log.debug("Loading reg ex from file " + filepath + ": '" + line +"'")
                        regexs.append(line)
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
        makePlotCollection(histoDicts, fileNames, fileDescr, outputFilePath, args.compare, args.one_file_per_histogram)

    log.info("done")    
    
if __name__ == "__main__":
    run( sys.argv[1:])
