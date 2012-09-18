"""
Created on March 2012

@author: Andrea Dotti (andrea.dotti@cern.ch)

This module contains definition of tests to perform
regression testing. There are tests  to
perform a statistical comparison between binned
and unbinned distributions.
Tests input are ROOT's histograms or TTrees.
When possible the ROOT functionalities are used.

See testme function for exanmple
"""
 
from Utils import DataType,BaseClass,Error
from Utils import getFrom1DHistogram,getFromTree,makeROOTgof

class SkipTest(Error):
    pass
class ZeroEntries(SkipTest):
    def __init__(self,hname):
        SkipTest.__init__(self,"No data: zero statistics, cannot perform test (%s)"%hname)

class StatTest:
    """
    Base Class definint a test object
    """
    def __init__(self,data1,data2,dataType):
        """
        Create a new test for the specified iunputs of
        type data type.
        Values are cached.
        """
        self.dataType = dataType
        self.dataset1 = data1
        self.dataset2 = data2
        self.cached = False
    def calc( self ):
        """
        If needed perform the
        calculations for tests.
        """
        if not self.cached:
            self.statvalue = self._stat()
            self.pvalue = self._pval()
            self.ndfvalue = self._ndf()
            self.cached = True
    def _stat(self ):
        """
        To be implemented by derived class:
        perform calculation of statistics value
        """
        raise BaseClass()
    def _pval(self):
        """
        To be implemented by derived class:
        perform calculation of p-value
        """
        raise BaseClass()
    def _ndf(self):
        """
        To be implemented by derived class:
        perform calculation of degrees of freedom
        """
        raise BaseClass()

    def value(self):
        """
        Return value of statistics
        """
        raise BaseClass()

    def pval(self):
        """
        Return p-value
        """
        if not self.cached:
            self.calc()
        return self.pvalue
    def stat(self):
        """
        Return value of statistics
        """
        if not self.cached:
            self.calc()
        return self.statvalue
    def ndf(self):
        """
        Return degrees-of-freedom
        """
        if not self.cached:
            self.calc()
        return self.ndfvalue

class Binned1DChi2Test(StatTest):
    """
    Chi2 Test for Binned distributions.
    Inputs: Two 1D histograms
    """
    def __init__(self,h1,h2,option="UU"):
        StatTest.__init__(self,h1,h2,DataType.BINNED1D)
        self.option = option
        self.cached = False
        from ROOT import TH1F
        theaxis = h1.GetXaxis()
        self.residuals = TH1F("%s-Residuals"%h1.GetName(),"%s Residuals"%h1.GetName(),
                              theaxis.GetNbins(),theaxis.GetXmin(),theaxis.GetXmax())
        self.localpval = None
        self.localchi2 = None
        self.localndf = None
    def value(self):
        return self.pval()
    def _pval(self):
        return self.localpval
    def _stat(self):
        #This is the first method to be called...
        #from ctypes import c_double
        mysize = self.residuals.GetXaxis().GetNbins()
        #myarray =  mysize * c_double
        #myres = myarray()
        #pval = self.dataset1.Chi2Test( self.dataset2 ,"UU",myres)
        import ROOT
        from array import array
        myres = array('d',mysize*[0])
        chi2 = ROOT.Double(0)
        ndf = ROOT.Long(0)
        igood=ROOT.Long(0)
        #Chi2 test does not work if there is no statistics        
        if self.dataset1.GetEntries() == 0:
            raise ZeroEntries(self.dataset1.GetName())
        if self.dataset2.GetEntries() == 0:
            raise ZeroEntries(self.dataset1.GetName())
        self.localpval = self.dataset1.Chi2TestX( self.dataset2,chi2,ndf,igood,self.option,myres)
        self.localchi2 = chi2
        self.localndf  = ndf
        for binidx in range(0,mysize):
            self.residuals.SetBinContent( binidx+1, myres[binidx] )
        return self.localpval
    def _ndf(self):
        return self.localndf

class BinnedWeighted1DChi2Test(Binned1DChi2Test):
    """
    Chi2 Test for Binned distributions, distributions are weighted.
    """
    def __init__(self,h1,h2):
        Binned1DChi2Test.__init__(self,h1,h2,"WW")

class AndersonDarlingTest(StatTest):
    """
    Anderson-Darling Test.
    Inputs: Two samples. Both binned or unbinned are accepted
    """
    def __init__(self,h1,h2,dataType=DataType.UNBINNED):
        if dataType == DataType.BINNED1D:
            first = getFrom1DHistogram( h1 )
            second= getFrom1DHistogram( h2 )
        elif dataType == DataType.UNBINNED:
            first = getFromTree( h1 )
            second = getFromTree( h2 )
            #print first[0:10]
            #print second[0:10]
        else:
            raise WrongDataType()
        StatTest.__init__(self,first,second,dataType)
        from Utils import makeROOTgof
        self.gof=makeROOTgof(self.dataset1,self.dataset2)
    def value(self):
        return self.pval()
    def _pval(self):
        return self.gof.AndersonDarling2SamplesTest("p")
    def _stat(self):
        return self.gof.AndersonDarling2SamplesTest("t")
    def _ndf(self):
        return 0

def BinnedAndersonDarlingTest(AndersonDarlingTest):
    def __init__(self,d1,d2):
        AndersonDarlingTest.__init__(self,d1,d2,DataType.BINNED1D)

class KolmogorovSmirnovTest(StatTest):
    """
    Kolmogorov-Smirnov Test.
    Inputs: Two samples. Both binned or unbinned are accepted (see ROOT's documentation
    for TH1::KolmogorovSmirnov test for binned distributions)
    """
    def __init__(self,d1,d2,dataType=DataType.UNBINNED):
        if dataType == DataType.BINNED1D:
            print 'WARNING: Using test on binned distributions'
            first = d1 #getFrom1DHistogram( d1 )
            second= d2 #getFrom1DHistogram( d2 )
        elif dataType == DataType.UNBINNED:
            first = getFromTree( d1 )
            second = getFromTree( d2 )
        else:
            raise WringDataType()
        StatTest.__init__(self,first,second,dataType)
        if dataType == DataType.UNBINNED:
            self.gof=makeROOTgof(self.dataset1,self.dataset2)
    def value(self):
        return self.pval()
    def _ndf(self):
        return 0
    def _pval(self):
        if self.dataType == DataType.UNBINNED:
            return self.gof.KolmogorovSmirnov2SamplesTest("p")
        else:
            return self.dataset1.KolmogorovTest( self.dataset2 ,"D")
    def _stat(self):
        if self.dataType == DataType.UNBINNED:
            return self.gof.KolmogorovSmirnov2SamplesTest("t")
        else:
            return self.dataset1.KolmogorovTest( self.dataset2 ,"M")
        
class BinnedKolmogorovSmirnovTest(KolmogorovSmirnovTest):
    def __init__(self,h1,h2):
        KolmogorovSmirnovTest.__init__(self,h1,h2,DataType.BINNED1D)
        
def testme():
    """
    Example function
    """
    from ROOT import TH1F
    h1 = TH1F("h1","h",100,-10,10)
    h2 = TH1F("h2","h",100,-10,10)
    h1.FillRandom("gaus")
    h2.FillRandom("gaus")
    h1.SetLineColor(2)
    h1.DrawCopy()
    h2.DrawCopy("same")
    binned = Binned1DChi2Test( h1 , h2 )
    print 'For Chi2Test:',binned.pval()
    binned = AndersonDarlingTest( h1 , h2 , DataType.BINNED1D)
    print 'For AD: ',binned.pval()
    binned = KolmogorovSmirnovTest( h1 , h2, DataType.BINNED1D)
    print 'For KS: ',binned.pval()


_testsMap = {
        "BinnedAndersonDarlingTest" : BinnedAndersonDarlingTest,
        "AndersonDarlingTest" : AndersonDarlingTest ,
        "Binned1DChi2Test" : Binned1DChi2Test,
        "BinnedWeighted1DChi2Test" : BinnedWeighted1DChi2Test,
        "BinnedKolmogorovSmirnovTest" : BinnedKolmogorovSmirnovTest,
        "KolmogorovSmirnovTest" : KolmogorovSmirnovTest
        }

def getTestByName( name ):
    """
    Returns the test class by name
    """
    if _testsMap.has_key(name):
        return _testsMap[name] 
    return StatTest
