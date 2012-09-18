"""
IO Module
This module defines the inout/output interfaces
"""

from Utils import Error,logger,NotYet,WrongDataType
_logger=logger().getLogger('IO')
        
class Input:
    """
    Base Class for input
    """
    def __init__(self):
        self.input = None
        self.name = None
    def getObject(self):
        """
        Returns the input object
        """
        return self.input
    def __str__(self):
        try:
            return 'Input: %s'%self.name
        except:
            return 'Input: %s'%self.name
    def __call__(self):
        raise BaseClass()
    
class FileInput(Input):
    """
    Read input from File.
    Currently ROOT's TFile inputs are supported
    """
    class FileType:
        ROOT = 1
    def __init__(self,filename,afiletype=FileType.ROOT):
        self.filetype = afiletype
        if filename != None:
            if self.filetype == FileInput.FileType.ROOT:
                from ROOT import TFile
                self.tfile = TFile.Open( filename )
            else:
                _logger.error("Functionality not yet implemented")
                raise NotYet()
        else:
            self.tfile = None
        Input.__init__(self)
    def _readHisto(self,name):
        self.input = self.tfile.Get( name )
        if not self.input.IsA().InheritsFrom("TH1"):
            raise WrongDataType()
    def _readBranch(self,treename,branchname):
        tree=self.tfile.Get(treename)
        if not tree.IsA().InheritsFrom("TTree"):
            raise WrongDataType()
        self.input = tree.GetBranch(branchname)
    def getFileObject(self):
        return self.tfile
    def __str__(self):
        fn = str(self.tfile)
        if self.tfile.__class__.__name__ == 'TFile':
            fn = self.tfile.GetName()
        return "%s (from file:%s)"%(Input.__str__(self),fn)

class Histogram(FileInput):
    """
    Get a histogram from a file
    """
    def __init__(self,afile,hname):
        """
        Create an Input histogram contained in file 'file'
        """
        #Check if file is a name (string) or object
        if type(afile)==type(str()):
            if afile.find('.root')>0:
                FileInput.__init__(self,afile,FileInput.FileType.ROOT)
            else:
                _logger.error("Cannot recognize filen: %s"%file)
                raise NotYet()
        else:
            FileInput.__init__(self,None)
            self.tfile=afile
        self._readHisto(hname)
    def __call__(self):
        return self.getObject()
    
class Branch(FileInput):
    """
    Get data from a Branch
    """
    class BranchType:
        """
        Supported data types: float, double, int, vector<double>
        """
        from ctypes import c_float,c_double,c_int
        FLOAT = 1 
        DOUBLE = 2
        INT = 3
        VECTORDOUBLE = 4 #Special case for std::vector<double>
        _typesMap={ FLOAT : c_float ,
                   DOUBLE : c_double ,
                   INT : c_int ,
                   VECTORDOUBLE : None , #Need special care
                   }
        @staticmethod
        def stringToValue( str ):
            if str == "FLOAT":
                return Branch.BranchType.FLOAT
            elif str == "DOUBLE":
                return Branch.BranchType.DOUBLE
            elif str == "INT":
                return Branch.BranchType.INT
            elif str == "VECTORDOUBLE":
                return Branch.BranchType.VECTORDOUBLE
            print "ERROR: Cannot recognize %s"%str
            raise WrongDataType()
    def __init__(self,afile,treename,brname,brtype=BranchType.DOUBLE,nelements=1,element=-1):
        """
        Create a Branch input
        @input afile: the file from where to read
        @input treename: the tree name
        @input brname: the name of the branch to read from
        @input brtype: the data type contained in the branch
        @input nelements: number of elements in the branch (i.e. array size of the branch)
        @input element: which element of the array to use (-1: all)
        Note: data are cached, use resetCache function to reset cache.
        """
        self._cache=None
        self.branchType = brtype
        self.nelements = nelements
        self.element = element
        #Check if file is a name (string) or object
        if type(afile)==type(str()):
            if afile.find('.root')>0:
                FileInput.__init__(self,afile,FileInput.FileType.ROOT)
            else:
                _logger.error("Cannot recognize file name: %s"%afile)
                raise NotYet()
        else:
            FileInput.__init__(self,None)
            self.tfile=afile
        self._readBranch(treename,brname)
    def resetCache(self):
        self._cache=None
    def __call__(self):
        if self.filetype == FileInput.FileType.ROOT:
            return self.readFromROOT()
        else:
            _logger.error("Input File Type not supported")
            raise NotYet()
    def readFromROOT(self):
        #Check if data have already been read,
        #if not read in branch data
        if self._cache == None:
            branch = self.getObject()
            if self.branchType == Branch.BranchType.VECTORDOUBLE:
                #Use some ROOT's stuff to read vector<double>
                #Note we need to use ROOT.Double to have correct data type in vectors
                from ROOT import std,Double,AddressOf
                #Create a vector to contain doubles
                vect = std.vector(Double)(self.nelements,0)
                data = AddressOf( vect )
            else:
                #Create an array to store data for each entry
                datatype = self.nelements * Branch.BranchType._typesMap[ self.branchType ] 
                data = datatype()
            #Set address of branch to the array
            branch.SetAddress( data )
            #Loop on events
            self._cache = []
            for entry in xrange( branch.GetEntries() ):
                branch.GetEntry( entry )
                #First simple case of simple variable
                if self.nelements == 1 and self.branchType != Branch.BranchType.VECTORDOUBLE:
                    self._cache.append( data[0] )
                else:
                    #Check if specific element is needed or all elements
                    if self.element == -1: 
                        #Loop on elements of arrays
                        d = []
                        for elem in xrange( self.nelements ):
                            if self.branchType == Branch.BranchType.VECTORDOUBLE:
                                d.append( vect.at(elem) )
                            else:
                                d.append( data[elem] )
                    else:
                        #Get array/vector element self.element
                        if self.branchType == Branch.BranchType.VECTORDOUBLE:
                            d = vect.at( self.element )
                        else:
                            d = data[self.element]
                    self._cache.append( d )
        return self._cache

    
def testme(fn="AtlasECAL_pi-_100_QGSP_BERT_95ref02.root"):
    input1 = Histogram( fn , "Spectra/hProtonSpectrum" )
    histogram1= input1()
    print input1,histogram1
    #Re-use open file
    input2 = Histogram( input1.getFileObject() , "Spectra/hNeutronSpectrum" )
    print input2,input2()

    input3 = Branch( "AtlasECAL_pi-_100_QGSP_BERT_95p01.root", "SimplifiedCalorimeter","EDEP_ACT" , Branch.BranchType.FLOAT )
    branch=input3()
    print input3,len(branch)

    input4 = Branch( input3.getFileObject() ,"SimplifiedCalorimeter","EDEP_CAL" , Branch.BranchType.FLOAT )
    br2=input4()
    print input4,len(input4()),br2[0:3] #Pring using cached values

    #Reading a std::vector<double>
    input5 = Branch( input3.getFileObject() , "SimplifiedCalorimeter","R",Branch.BranchType.VECTORDOUBLE,10)
    br3 = input5()
    print input4,len(br3),br3[0]
