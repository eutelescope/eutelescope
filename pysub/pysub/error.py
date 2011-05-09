## @package pysub.error
# This package contains all the errors defined and used
# within the pysub @package.

## Base class for errors in Pysub
#
class PysubError( Exception ):
    pass

## Error with Marlin execution
#
class MarlinError( PysubError ):

    def __init__ (self, message, errno ):
        self._message = message
        self._errno   = errno

## Problem with the GRID submission
class GRIDSubmissionError( PysubError ):

    def __init__ (self, message):
        self._message = message

## Problem with file verification
class VerificationError( PysubError ):
    pass


## Missing generic file
#
class MissingFileError ( PysubError ):
    def __init__ (self, filename ):
        self._filename = filename

## Missing configuration file
#
class MissingConfigurationFileError( MissingFileError ):
    pass

## Missing input file
#
class MissingInputFileError( MissingFileError ):
    pass

## Missing grid library file
class MissingLibraryFileError( MissingFileError ):
    pass

## Missing GEAR file
#
class MissingGEARFileError( MissingFileError ):
    pass

## Missing Steering template
#
class MissingSteeringTemplateError( MissingFileError ):
    pass

## Missing hotpixel file
class MissingHotPixelFileError( MissingFileError ):
    pass

## Missing offset file
class MissingOffsetFileError( MissingFileError ):
    pass

## Missing alignment file
class MissingAlignmentFileError( MissingFileError ):
    pass

## Missing pedestal file
class MissingPedestalFileError( MissingFileError ):
    pass

## Missing output file
#
class MissingOutputFileError( MissingFileError ):
    pass

## Missing histogram file
#
class MissingHistogramFileError( MissingFileError ):
    pass

## Missing joboutput file
#
class MissingJoboutputFileError( MissingFileError ):
    pass

## Missing GRID folder
#
class MissingGRIDFolderError( MissingFileError ) :
    pass

## Missing file on the GRID
#
class MissingFileOnGRIDError( MissingFileError ):
    pass

## Missing input file on the GRID
#
class MissingInputFileOnGRIDError ( MissingFileOnGRIDError ):
    pass

## Missing hotpixel file on the GRID
#
class MissingHotPixelFileOnGRIDError ( MissingFileOnGRIDError ) :
    pass

## Missing pedestal file on the GRID
#
class MissingPedestalFileOnGRIDError ( MissingFileOnGRIDError ) :
    pass

## Missing offset file on the GRID
#
class MissingOffsetFileOnGRIDError ( MissingFileOnGRIDError ) :
    pass

## Missing alignment file on the GRID
#
class MissingAlignmentFileOnGRIDError ( MissingFileOnGRIDError ) :
    pass

## File already on the GRID
#
class FileAlreadyOnGRIDError( MissingFileError ):
    pass


## Output file already on the GRID
#
class OutputFileAlreadyOnGRIDError( FileAlreadyOnGRIDError ):
    pass

## Histogram file already on the GRID
#
class HistogramFileAlreadyOnGRIDError( FileAlreadyOnGRIDError ):
    pass

## Joboutput file already on the GRID
#
class JoboutputFileAlreadyOnGRIDError( FileAlreadyOnGRIDError ) :
    pass

## Not enough file for continue 
#
class NotEnoughFilesError( PysubError ):
    pass

## Problem copying file from GRID
#
class GRID_LCG_CPError( MissingFileError ):
    pass

## Problem copying / registering file to GRID
#
class GRID_LCG_CRError( MissingFileError ):
    pass

## Stop execution error
#
# This error is so critical that the execution has to
# be terminated.
#
class StopExecutionError( PysubError ):
    def __init__( self, message ):
        self._message = message
