// Author: Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#include "EUTELESCOPE.h"

// system includes
#include <algorithm>
#include <string>

using namespace eutelescope;

const std::string EUTELESCOPE::GEOFILENAME		= "telescope_geometry.root";
const bool EUTELESCOPE::DUMPGEOROOT				= true;

const char *   EUTELESCOPE::HEADERVERSION       = "HeaderVersion";
const char *   EUTELESCOPE::NOOFEVENT           = "NoOfEvent";
const char *   EUTELESCOPE::DATATYPE            = "DataType";
const char *   EUTELESCOPE::DATETIME            = "DateTime";
const char *   EUTELESCOPE::DAQHWNAME           = "DAQHWName";
const char *   EUTELESCOPE::DAQHWVERSION        = "DAQHWVersion";
const char *   EUTELESCOPE::DAQSWNAME           = "DAQSWName";
const char *   EUTELESCOPE::DAQSWVERSION        = "DAQSWVersion";
const char *   EUTELESCOPE::SIMULSWNAME         = "SimulSWName";
const char *   EUTELESCOPE::SIMULSWVERSION      = "SimulSWVersion";
const char *   EUTELESCOPE::GEOID               = "GeoID";
const char *   EUTELESCOPE::EUDRBMODE           = "EUDRBMode";
const char *   EUTELESCOPE::EUDRBDET            = "EUDRBDet";
const char *   EUTELESCOPE::BEAMLOCATION        = "BeamLocation";
const char *   EUTELESCOPE::BEAMTYPE            = "BeamType";
const char *   EUTELESCOPE::BEAMENERGY          = "BeamEnergy";
const char *   EUTELESCOPE::NOOFDETECTOR        = "NoOfDetector";
const char *   EUTELESCOPE::MINX                = "MinX";
const char *   EUTELESCOPE::MAXX                = "MaxX";
const char *   EUTELESCOPE::MINY                = "MinY";
const char *   EUTELESCOPE::MAXY                = "MaxY";
const char *   EUTELESCOPE::INTERMEDIATEFILE    = "IntermediateFile";
const char *   EUTELESCOPE::USERCOMMENT         = "UserComment";
const char *   EUTELESCOPE::APPLIEDPROCESSOR    = "AppliedProcessor";
const char *   EUTELESCOPE::DAQDATA             = "DAQData";
const char *   EUTELESCOPE::SIMULDATA           = "SimulData";
const char *   EUTELESCOPE::CONVDATA            = "ConvData";
const char *   EUTELESCOPE::EUDRB               = "EUDRB";
const char *   EUTELESCOPE::IPHCIMAGER          = "IPHCIMAGER";
const char *   EUTELESCOPE::SUCIMAIMAGER        = "SUCIMAIMAGER";
const char *   EUTELESCOPE::EUDAQ               = "EUDAQ";
const int      EUTELESCOPE::GOODPIXEL           =  0;
const char *   EUTELESCOPE::EVENTTYPE           = "EventType";
const int      EUTELESCOPE::BADPIXEL            =  1;
const int      EUTELESCOPE::HITPIXEL            = -1;
const int      EUTELESCOPE::MISSINGPIXEL        =  2;
const int      EUTELESCOPE::FIRINGPIXEL         =  3;
const char *   EUTELESCOPE::ABSOLUTENOISEVALUE  = "AbsoluteNoiseValue";
const char *   EUTELESCOPE::NOISEDISTRIBUTION   = "NoiseDistribution";
const char *   EUTELESCOPE::DEADPIXEL           = "DeadPixel";
const char *   EUTELESCOPE::ABSOLUTEPEDEVALUE   = "AbsolutePedeValue";
const char *   EUTELESCOPE::FULLFRAME           = "FullFrame";
const char *   EUTELESCOPE::ROWWISE             = "RowWise";
const char *   EUTELESCOPE::MEANRMS             = "MeanRMS";
const char *   EUTELESCOPE::AIDAPROFILE         = "AIDAProfile";
const char *   EUTELESCOPE::FIXEDFRAME          = "FixedFrame";
const char *   EUTELESCOPE::DFIXEDFRAME         = "DFixedFrame";
const char *   EUTELESCOPE::SPARSECLUSTER       = "SparseCluster";
const char *   EUTELESCOPE::SPARSECLUSTER2      = "SparseCluster2";
const char *   EUTELESCOPE::SPARSECLUSTER3      = "SparseCluster3";
const char *   EUTELESCOPE::BRICKEDCLUSTER      = "BrickedCluster";
const char *   EUTELESCOPE::ANALOG              = "Analog";
const char *   EUTELESCOPE::DIGITAL             = "Digital";
const char *   EUTELESCOPE::BINARY              = "Binary";
const char *   EUTELESCOPE::FLAGONLY            = "FlagOnly";
const char *   EUTELESCOPE::MATRIXDEFAULTENCODING    = "sensorID:7,xMin:12,xMax:12,yMin:12,yMax:12";
const char *   EUTELESCOPE::CLUSTERDEFAULTENCODING   = "sensorID:7,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,quality:7";

const char *   EUTELESCOPE::PULSEDEFAULTENCODING     = "sensorID:7,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5,type:5,quality:5";
const char *   EUTELESCOPE::ZSDATADEFAULTENCODING    = "sensorID:7,sparsePixelType:5";
const char *   EUTELESCOPE::ZSCLUSTERDEFAULTENCODING = "sensorID:7,sparsePixelType:5,quality:5";
const char *   EUTELESCOPE::HITENCODING              = "sensorID:7,properties:7";
const char *   EUTELESCOPE::FIXEDWEIGHT              = "FixedWeight";


namespace eutelescope {

  ClusterQuality operator&(ClusterQuality a, ClusterQuality b) {
    return ClusterQuality( static_cast<int>(a) & static_cast<int>(b) );
  }

  ClusterQuality operator|(ClusterQuality a, ClusterQuality b) {
    return ClusterQuality( static_cast<int>(a) | static_cast<int>(b) );
  }

  ClusterQuality operator|=(ClusterQuality a, ClusterQuality b) {
    return ClusterQuality( static_cast<int>(a) | static_cast<int>(b) );
  }

  std::ostream& operator<<(std::ostream& os, const ClusterQuality& quality) {
    if ( quality == kGoodCluster ) {
      os << "kGoodCluster (" << static_cast<int > (quality) << ")";
      return os;
    }
    bool moreThanOne = false;
    if ( quality & kIncompleteCluster ) {
      if ( moreThanOne ) os << ", ";
      os << "kIncompleteCluster";
      moreThanOne = true;
    }
    if ( quality & kBorderCluster ) {
      if ( moreThanOne ) os << ", ";
      os << "kBorderCluster";
      moreThanOne = true;
    }
    if ( quality & kMergedCluster ) {
      if ( moreThanOne ) os << ", ";
      os << "kMergedCluster";
      moreThanOne = true;
    }
    os << " (" << static_cast<int> (quality) << ")";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const SparsePixelType& type) {
    if ( type == kEUTelBaseSparsePixel ) os << "kEUTelBaseSparsePixel" ;
    else if ( type == kEUTelSimpleSparsePixel ) os << "kEUTelSimpleSparsePixel";
    else if ( type == kEUTelGenericSparsePixel ) os << "kEUTelGenericSparsePixel";
    else if ( type == kEUTelGeometricPixel ) os << "kEUTelGeometricPixel";
    // add here your type
    else if ( type == kUnknownPixelType ) os << "kUnknownPixelType";
    os << " (" << static_cast<int> (type ) << ")";
    return os;
  }

  std::string ucase( const std::string & inputString ) {

    std::string result( inputString );
    transform( result.begin(), result.end(), result.begin(), ::toupper);
    return result;

  }

  std::string lcase( const std::string & inputString ) {

    std::string result( inputString );
    transform( result.begin(), result.end(), result.begin(), ::tolower);
    return result;

  }


  std::string trim( const std::string & inputString ) {

    // define all possible space types you want to remove
    static const std::string spaces = "\t\n\r\v";
    // find the first character that is not a space
    size_t b = inputString.find_first_not_of( spaces );
    // find the last character that is not a space
    size_t e = inputString.find_first_not_of( spaces );
    if ( b == std::string::npos || e == std::string::npos ) {
      // the string is empty!
      return "";
    }
    return std::string(inputString, b, e - b + 1 );
  }

  std::string escape( const std::string & inputString ) {

    std::ostringstream ret;
    ret << std::setfill('0') << std::hex;
    for ( size_t iPos = 0; iPos < inputString.length(); ++iPos ) {
      if ( inputString[iPos] == '\\' ) ret << "\\\\";
      else if ( inputString[iPos] < 32 ) ret << "\\x" << std::setw(2) << int( inputString[ iPos ] );
      else ret << inputString[ iPos ] ;
    }
    return ret.str();

  }

  std::string firstline( const std::string & inputString ) {

    size_t iPos = inputString.find( '\n' );
    return inputString.substr( 0, iPos );

  }

  std::vector<std::string > split( const std::string & inputString, const std::string & delim ) {
    std::string copy( inputString );
    std::vector<std::string > result;
    size_t iPos;
    while ( ( iPos = copy.find_first_of( delim ) ) != std::string::npos ) {
      result.push_back( copy.substr( 0, iPos ) );
      copy = copy.substr( iPos + 1 );
    }
    result.push_back( copy );
    return result;
  }
}
