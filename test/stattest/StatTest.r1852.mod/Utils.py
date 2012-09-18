import array
def getFrom1DHistogram( histo ):
    """
    Get histogram content as list.
    Excluding under-flow and over-flow bins
    """
    data=[]
    for bin in range(histo.GetXaxis().GetFirst(),histo.GetXaxis().GetLast()+1):
        binctr = histo.GetXaxis().GetBinCenter(bin)
        bincon = histo.GetBinContent(bin)
        if ( bincon > 0 ):
            data+=( [ binctr for e in xrange(0,int(bincon)) ] )
    data.sort()
    return array.array('d',data)

def getFromTree( branch ):
    """
    Get data from Branch input as array of doubles
    """
    return array.array('d',branch())

def draw( obj , linecol , opt="" , limits=None , title=None):
    """
    Draw an histogram of object obj, with specified line color.
    obj can be a histogram or a dataset.
    """
    if "TH1" in obj.__class__.__name__:
        obj.SetLineColor( linecol )
        obj.Draw(opt)
        hh = obj
    else:
        from ROOT import TH1F
        hh=TH1F( title , title , limits[0],limits[1],limits[2])
        [ hh.Fill( e ) for e in obj ]
        hh.SetLineColor( linecol )
        hh.DrawCopy(opt)
    return hh
               
def makeResiduals( h1 , h2 ):
    """
    Create residuals histogram
    Copyied from TH1::Chi2TestX for UU 1D distributions
    """
    hres = h1.Clone()
    hres.SetName( h1.GetName() + "-Residuals" )
    hres.SetTitle( "Residuals" )
    hres.Reset()
    sum1 = 0
    sum2 = 0
    for bin in range( h1.GetXaxis().GetFirst() , h1.GetXaxis().GetLast() ):
        sum1 += h1.GetBinContent( bin )
        sum2 += h2.GetBinContent( bin )
    sum = sum1+sum2
    from math import sqrt
    for bin in range( h1.GetXaxis().GetFirst() , h1.GetXaxis().GetLast() ):
        bin1 = h1.GetBinContent( bin )
        bin2 = h2.GetBinContent( bin )
        binsum= bin1+bin2
        if binsum > 0:
            nexp1 = binsum*sum1/sum
            res = (bin1-nexp1)/sqrt(nexp1)
            #Habermann correction for residuals
            correc = (1.-sum1/sum)*(1.-binsum/sum)
            res /= sqrt(correc)
            hres.SetBinContent( bin, res )
    return hres
    
    
                
def makeROOTgof( data1, data2):
    """
    Utility function, creates a ROOT's GoFTest class instance
    """
    from ROOT import Math
    return Math.GoFTest( len(data1) , data1, len(data2), data2)


import logging
class logger:
    """ logging class
    this is a simple interface to log messages
    """
    def __init__(self):
        self.FORMAT = '$(asctime)s %(name)s %(levelname)s: %(message)s'
        self.LEVEL = logging.WARNING
        try:
            logging.basicConfig(level=self.LEVEL, format = self.FORMAT)
        except TypeError: #< 2.5
            logging.basicConfig()
        l=logging.getLogger()
        l.setLevel(self.LEVEL)
    def getLogger(self,name):
        return logging.getLogger()

class DataType:
    BINNED1D = 1
    UNBINNED = 2

class Error(Exception):
    def __init__(self,msg):
        Exception.__init__(self,"StatTest",msg)
        self.message=msg
        self.logger = logger().getLogger('StatTest')
        self.logger.critical(str(self.message))

class WrongDataType(Error):
    def __init__(self):
        Error.__init__(self,"DataType not recognized")
class BaseClass(Error):
    def __init__(self):
        Error.__init__(self,"Use derived class")
        
class NotYet(Error):
    def __init__(self):
        Error.__init__(self,"Not yet implemented")

        #import Interface
        #import Tests
