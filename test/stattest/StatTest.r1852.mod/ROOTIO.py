"""
ROOT interfaces and utilities
"""
from Utils import WrongDataType

def stripPathName( aPath ):
    """
    Remove file name identifier and trailing / from ROOT
    TDirectory path
    """
    return aPath.GetPath()[aPath.GetPath().find(":")+1:].rstrip('/')

def readDirectory( tdir ):
    """
    Recursively read content of ROOT's TDirectory tdir
    and add full pathnames to outputlist
    """
    outputlist=[]
    from ROOT import TIter
    nextkey = TIter( tdir.GetListOfKeys() )
    key = nextkey()
    while key:
        obj = key.ReadObj()
        if obj.IsA().InheritsFrom("TH1") or obj.IsA().InheritsFrom("TTree"):
            curdirname = stripPathName(tdir)
            if len(curdirname)>0:
                objpath = "%s/%s"%(stripPathName(tdir),obj.GetName())
            else:
                objpath = obj.GetName()
                #print obj.GetName()," :-> ",objpath
            outputlist.append( objpath )
        elif obj.IsA().InheritsFrom("TDirectory"):
            from ROOT import gDirectory
            #print gDirectory.GetPath(),obj.GetName(),obj.GetPath()
            outputlist += readDirectory( obj )
        key=nextkey()
    return outputlist

def checkTreeHasBranch( fileobj , treename , branchnane ):
    """
    Check if Tree with name treename contained in TFile fileobj,
    has a TBranch with name matching the regular expression branchname.
    If found returns the name of the branch, otherwise None.
    Note: Returns the first occurance of a branch matching the re
    """
    thetree=fileobj.Get(treename)
    if not thetree.IsA().InheritsFrom("TTree"):
        return None
    import re
    mm=re.compile("^%s$"%branchnane)
    for bridx in range(thetree.GetNbranches()):
        candidate = thetree.GetListOfBranches()[bridx].GetName()
        if mm.match( candidate ):
            return candidate
    return None

def buildUnbinnedInputs(tfile,paths):
    """
    Creates and returns a list of StatTest.IO.Branch objects
    corresponding to the specified list of tuples (treename, branchname,conf)
    """
    inputs=[]
    from IO import Branch
    from Interface import Configuration
    for treename,branchname,conf in paths:
        #print treename,branchname,conf #########################
        #print 'Creo Branch in file: ',tfile,' Da Tree.Branch: ',treename,branchname,
        #print 'Typo container: ',conf[Configuration.TYPEKEY],' Di dimensione: ',conf[Configuration.SIZEKEY]
        if conf.has_key(Configuration.ELEMENTKEY):          
            un = Branch( tfile , treename , branchname , 
                        Branch.BranchType.stringToValue( conf[Configuration.TYPEKEY] ),
                        conf[Configuration.SIZEKEY],
                        conf[Configuration.ELEMENTKEY])
        else:
            un = Branch( tfile , treename , branchname , 
                        Branch.BranchType.stringToValue( conf[Configuration.TYPEKEY] ),
                        conf[Configuration.SIZEKEY])
        un.name = "%s:%s"%(treename,branchname)
        inputs.append( un )
    return inputs

def buildHistogramInputs( tfile, paths ):
    """
    Creates and returns a list of StatTest.IO.Histogram objects
    correponding to the specified list of paths
    """
    inputs = []
    from IO import Histogram
    for objectname in paths:
        try:
            hh = Histogram(tfile,objectname)
            hh.name = objectname
            inputs.append( hh )
        except WrongDataType:
            #Not an histogram, skip
            pass
    return inputs

def makePage( algorithm , pagename , cdash = False, prefix=""):
    from ROOT import TCanvas,kBlue,kRed,gROOT,kGreen,kYellow
    gROOT.SetBatch(True)
    c=TCanvas( algorithm.output.name , algorithm.output.name )
    c.Divide(1,2)
    from Interface import Result
    aColor = None
    if algorithm.output.result == Result.FAILED:
        aColor = kRed
    if algorithm.output.result == Result.NOTPASSED:
        aColor = kYellow
    if algorithm.output.result == Result.SUCCESS:
        aColor = kGreen
    if aColor:
        c.SetFillColor( aColor )
    aPad = c.cd(1)
    from Utils import draw
    lims = ()
    if "TH1" not in algorithm.test.dataset1.__class__.__name__:
        lims = ( 100, 
                min( algorithm.test.dataset1.tolist() + algorithm.test.dataset2.tolist() ),
                max( algorithm.test.dataset1.tolist() + algorithm.test.dataset2.tolist() )
               )
    h1=draw( algorithm.test.dataset1 , kBlue , ""    , lims , algorithm.output.name )
    h1.SetName(h1.GetName()+"_new")
    aPad.Update()
    from ROOT import TPaveStats
    statBox = h1.GetListOfFunctions().FindObject("stats")
    statBox.SetName('new_stat')
    statBox.SetY1NDC(statBox.GetY1NDC()-.18)
    statBox.SetY2NDC(statBox.GetY2NDC()-.18)
    statBox.SetTextColor(kBlue)
    statBox.SetBorderSize(2)
    h2=draw( algorithm.test.dataset2 , kRed  , "sames", lims , algorithm.output.name+"ref")
    h2.SetName(h2.GetName()+"_ref")
    aPad.Update()
    statBox2 = h2.GetListOfFunctions().FindObject("stats")
    statBox2.SetName('ref_stat')
    statBox2.SetTextColor(kRed)
    from ROOT import TPaveText
    pave=TPaveText(0.02,0.85,0.35,0.99,"NDC")
    pave.SetTextColor(aColor)
    pave.SetFillColor(1)
    pave.AddText(" %s "%algorithm.output.result)
    pave.AddText("(p-val: %s Test: %s)"%(algorithm.output.value,
                                         algorithm.test.__class__.__name__))
    pave.Draw()
    c.cd(2)
    if 'residuals' in algorithm.test.__dict__:
        algorithm.test.residuals.Draw()
    else:
        from Utils import makeResiduals
        algorithm.test.residuals = makeResiduals( h1 , h2 )
        algorithm.test.residuals.Draw()
    c.Print(pagename+prefix)
    # only print CDASH info if the test has failed (saves space on webserver)
    if cdash == True and not (algorithm.output.result == Result.SUCCESS):
        import os
        c.Print(os.path.dirname(pagename)+"/"+algorithm.output.name.replace('/','_')+".png")
        print " <DartMeasurementFile name=\"" + algorithm.output.name +  "\" type=\"image/png\"> " + os.path.dirname(pagename) + "/" + algorithm.output.name.replace('/','_') + ".png" + " </DartMeasurementFile>"

def testme( filename="AtlasECAL_pi-_100_QGSP_BERT_95ref02.root" ):
    """
    Test function
    """
    from ROOT import TFile
    f=TFile.Open(filename)
    output=readDirectory( f )
    for name in output:
        print "Full path name: %s for object of name: %s and type: %s"%(name,f.Get(name).GetName(),f.Get(name).IsA().GetName())
        #print output
    return buildHistogramInputs(f, output )


