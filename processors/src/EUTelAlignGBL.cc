// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Authors:
// Philipp Roloff, DESY <mailto:philipp.roloff@desy.de>
// Joerg Behr, Hamburg Uni/DESY <joerg.behr@desy.de>
// Slava Libov, DESY <mailto:vladyslav.libov@desy.de>
// Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
// Daniel Pitzl, DESY <mailto:daniel.pitzl@desy.de>
//
// Version: $Id: EUTelAlignGBL.cc,v 1.48 2009-08-01 10:49:46 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if GEAR and MARLINUTIL are used
#if defined(USE_GEAR) && defined(USE_MARLINUTIL)

// eutelescope includes ".h"
#include "EUTelAlignGBL.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelExceptions.h"
#include "EUTelPStream.h" // process streams redi::ipstream
#include "EUTelAlignmentConstant.h"

// GBL:
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// marlin util includes
#include "Mille.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

// ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TRandom.h>
//#include <TMinuit.h>
#include <TSystem.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

AIDA::IHistogram1D * nAllHitHistoGBLAlign;

AIDA::IHistogram1D * dx01HistGBLAlign;
AIDA::IHistogram1D * dy01HistGBLAlign;
AIDA::IHistogram1D * dx02HistGBLAlign;
AIDA::IHistogram1D * dy02HistGBLAlign;
AIDA::IHistogram1D * dx03HistGBLAlign;
AIDA::IHistogram1D * dy03HistGBLAlign;
AIDA::IHistogram1D * dx04HistGBLAlign;
AIDA::IHistogram1D * dy04HistGBLAlign;
AIDA::IHistogram1D * dx05HistGBLAlign;
AIDA::IHistogram1D * dy05HistGBLAlign;

AIDA::IHistogram1D * tridxHistGBLAlign;
AIDA::IHistogram1D * tridyHistGBLAlign;
AIDA::IHistogram1D * tridxcHistGBLAlign;
AIDA::IHistogram1D * tridycHistGBLAlign;
AIDA::IHistogram1D * ntriHistGBLAlign;

AIDA::IHistogram1D * dridxHistGBLAlign;
AIDA::IHistogram1D * dridyHistGBLAlign;
AIDA::IHistogram1D * dridxcHistGBLAlign;
AIDA::IHistogram1D * dridycHistGBLAlign;
AIDA::IHistogram1D * ndriHistGBLAlign;

AIDA::IHistogram1D * sixdxHistGBLAlign;
AIDA::IHistogram1D * sixdyHistGBLAlign;
AIDA::IHistogram1D * sixkxHistGBLAlign;
AIDA::IHistogram1D * sixkyHistGBLAlign;

AIDA::IHistogram1D * selxHistGBLAlign;
AIDA::IHistogram1D * selyHistGBLAlign;
AIDA::IHistogram1D * selaxHistGBLAlign;
AIDA::IHistogram1D * selayHistGBLAlign;
AIDA::IHistogram1D * seldxHistGBLAlign;
AIDA::IHistogram1D * seldyHistGBLAlign;
AIDA::IHistogram1D * selkxHistGBLAlign;
AIDA::IHistogram1D * selkyHistGBLAlign;

AIDA::IHistogram1D * seldx1HistGBLAlign;
AIDA::IHistogram1D * seldy1HistGBLAlign;
AIDA::IHistogram1D * seldx3HistGBLAlign;
AIDA::IHistogram1D * seldy3HistGBLAlign;
AIDA::IHistogram1D * seldx4HistGBLAlign;
AIDA::IHistogram1D * seldy4HistGBLAlign;
AIDA::IHistogram1D * seldx5HistGBLAlign;
AIDA::IHistogram1D * seldy5HistGBLAlign;
AIDA::IHistogram1D * seldx6HistGBLAlign;
AIDA::IHistogram1D * seldy6HistGBLAlign;

AIDA::IHistogram1D * gblndfHistGBLAlign;
AIDA::IHistogram1D * gblchi2HistGBLAlign;
AIDA::IHistogram1D * gblprbHistGBLAlign;

AIDA::IHistogram1D * badxHistGBLAlign;
AIDA::IHistogram1D * badyHistGBLAlign;
AIDA::IHistogram1D * badaxHistGBLAlign;
AIDA::IHistogram1D * badayHistGBLAlign;
AIDA::IHistogram1D * baddxHistGBLAlign;
AIDA::IHistogram1D * baddyHistGBLAlign;
AIDA::IHistogram1D * badkxHistGBLAlign;
AIDA::IHistogram1D * badkyHistGBLAlign;

AIDA::IHistogram1D * baddx1HistGBLAlign;
AIDA::IHistogram1D * baddy1HistGBLAlign;
AIDA::IHistogram1D * baddx3HistGBLAlign;
AIDA::IHistogram1D * baddy3HistGBLAlign;
AIDA::IHistogram1D * baddx4HistGBLAlign;
AIDA::IHistogram1D * baddy4HistGBLAlign;
AIDA::IHistogram1D * baddx5HistGBLAlign;
AIDA::IHistogram1D * baddy5HistGBLAlign;
AIDA::IHistogram1D * baddx6HistGBLAlign;
AIDA::IHistogram1D * baddy6HistGBLAlign;

AIDA::IHistogram1D * goodx1HistGBLAlign;
AIDA::IHistogram1D * goody1HistGBLAlign;
AIDA::IHistogram1D * goodx6HistGBLAlign;
AIDA::IHistogram1D * goody6HistGBLAlign;

AIDA::IHistogram1D * gblax0HistGBLAlign;
AIDA::IHistogram1D * gbldx0HistGBLAlign;
AIDA::IHistogram1D * gblrx0HistGBLAlign;
AIDA::IHistogram1D * gblry0HistGBLAlign;
AIDA::IHistogram1D * gblpx0HistGBLAlign;
AIDA::IHistogram1D * gblpy0HistGBLAlign;

AIDA::IHistogram1D * gblax1HistGBLAlign;
AIDA::IHistogram1D * gbldx1HistGBLAlign;
AIDA::IHistogram1D * gblrx1HistGBLAlign;
AIDA::IHistogram1D * gblry1HistGBLAlign;

AIDA::IHistogram1D * gblax2HistGBLAlign;
AIDA::IHistogram1D * gbldx2HistGBLAlign;
AIDA::IHistogram1D * gblrx2HistGBLAlign;
AIDA::IHistogram1D * gblry2HistGBLAlign;

AIDA::IHistogram1D * gblax3HistGBLAlign;
AIDA::IHistogram1D * gbldx3HistGBLAlign;
AIDA::IHistogram1D * gblrx3HistGBLAlign;
AIDA::IHistogram1D * gblry3HistGBLAlign;

AIDA::IHistogram1D * gblax4HistGBLAlign;
AIDA::IHistogram1D * gbldx4HistGBLAlign;
AIDA::IHistogram1D * gblrx4HistGBLAlign;
AIDA::IHistogram1D * gblry4HistGBLAlign;

AIDA::IHistogram1D * gblax5HistGBLAlign;
AIDA::IHistogram1D * gbldx5HistGBLAlign;
AIDA::IHistogram1D * gblrx5HistGBLAlign;
AIDA::IHistogram1D * gblry5HistGBLAlign;

AIDA::IHistogram1D * gblax6HistGBLAlign;
AIDA::IHistogram1D * gbldx6HistGBLAlign;
AIDA::IHistogram1D * gbldy6HistGBLAlign;
AIDA::IHistogram1D * gblrx6HistGBLAlign;
AIDA::IHistogram1D * gblry6HistGBLAlign;

AIDA::IHistogram1D * gblkx1HistGBLAlign;
AIDA::IHistogram1D * gblkx2HistGBLAlign;
AIDA::IHistogram1D * gblkx3HistGBLAlign;
AIDA::IHistogram1D * gblkx4HistGBLAlign;
AIDA::IHistogram1D * gblkx5HistGBLAlign;
AIDA::IHistogram1D * gblkx6HistGBLAlign;

AIDA::IHistogram1D * nmHistGBLAlign;

#endif

gbl::MilleBinary * milleAlignGBL; // for producing MillePede-II binary file

//------------------------------------------------------------------------------
EUTelAlignGBL::EUTelAlignGBL(): Processor("EUTelAlignGBL") {

  //some default values
  FloatVec MinimalResidualsX;
  FloatVec MinimalResidualsY;
  FloatVec MaximalResidualsX;
  FloatVec MaximalResidualsY;

  FloatVec PedeUserStartValuesX;
  FloatVec PedeUserStartValuesY;

  FloatVec PedeUserStartValuesGamma;

  FloatVec SensorZPositions;

  FloatVec SensorXShifts;
  FloatVec SensorYShifts;

  FloatVec SensorGamma;

  FloatVec SensorAlpha;

  FloatVec SensorBeta;

  //maybe one has to chose a larger value than 6?

  for( int i = 0; i < 6; i++ ) {

    MinimalResidualsX.push_back(0.0);
    MinimalResidualsY.push_back(0.0);
    MaximalResidualsX.push_back(0.0);
    MaximalResidualsY.push_back(0.0);

    PedeUserStartValuesX.push_back(0.0);
    PedeUserStartValuesY.push_back(0.0);

    PedeUserStartValuesGamma.push_back(0.0);

    float zpos = 20000.0 +  20000.0 * (float)i; // [um]
    SensorZPositions.push_back(zpos);

    SensorXShifts.push_back(0.0);
    SensorYShifts.push_back(0.0);

    SensorGamma.push_back(0.0);
    SensorAlpha.push_back(0.0);
    SensorBeta.push_back(0.0);
  }

  // modify processor description
  _description =
    "EUTelAlignGBL uses the MILLE program to write data files for MILLEPEDE II.";

  // input collections
  std::vector<std::string > HitCollectionNameVecExample;
  HitCollectionNameVecExample.push_back("corrhits");

  registerInputCollections(LCIO::TRACKERHIT,"HitCollectionName",
      "Hit collections name",
      _hitCollectionName,HitCollectionNameVecExample);

  registerProcessorParameter( "Ebeam",
      "Beam energy [GeV]",
      _eBeam, static_cast < double >( 4.0));

  registerProcessorParameter( "IsFirstAlignStep",
      "Bool: yes/no",
      _IsFirstAlignStep, static_cast < int >( 0));

  registerOptionalParameter("ExcludePlanes","Exclude planes from fit according to their sensor ids.",_excludePlanes_sensorIDs ,std::vector<int>());

  registerOptionalParameter("FixedPlanes","Fix sensor planes in the fit according to their sensor ids.",_FixedPlanes_sensorIDs ,std::vector<int>());


  registerOptionalParameter("MaxTrackCandidatesTotal","Maximal number of track candidates (Total).",_maxTrackCandidatesTotal, static_cast <int> (10000000));
  registerOptionalParameter("MaxTrackCandidates","Maximal number of track candidates.",_maxTrackCandidates, static_cast <int> (2000));

  registerOptionalParameter("BinaryFilename","Name of the Millepede binary file.",_binaryFilename, string ("mille.bin"));

  registerOptionalParameter("TelescopeResolution","Resolution of the telescope for Millepede (sigma_x=sigma_y.",_telescopeResolution, static_cast <float> (0.010));

  registerOptionalParameter("AlignMode","Number of alignment constants used. Available mode are: "
      "\nXYZShifts - shifts in X and Y"
      "\nXYShiftsRotZ - shifts in X and Y and rotation around the Z axis,"
      "\nXYZShiftsRotZ - shifts in X,Y and Z and rotation around the Z axis",
      _alignModeString, std::string("XYShiftsRotZ"));

  registerOptionalParameter("AlignmentConstantLCIOFile","This is the name of the LCIO file name with the output alignment"
      "constants (add .slcio)",_alignmentConstantLCIOFile, static_cast< string > ( "alignment.slcio" ) );

  registerOptionalParameter("AlignmentConstantCollectionName", "This is the name of the alignment collection to be saved into the slcio file",
      _alignmentConstantCollectionName, static_cast< string > ( "alignment" ));

  registerOptionalParameter( "triCut", "Upstream triplet residual cut [um]", _triCut, 0.30 );
  registerOptionalParameter( "driCut", "Downstream triplet residual cut [um]", _driCut, 0.40 );
  registerOptionalParameter( "sixCut", "Upstream-Downstream Track matching cut [um]", _sixCut, 0.60 );
  registerOptionalParameter( "slopeCut", "t(d)riplet slope cut [radian]", _slopeCut, 0.01 );

  registerOptionalParameter("GeneratePedeSteerfile","Generate a steering file for the pede program.",_generatePedeSteerfile, static_cast <int> (0));

  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt"));

  registerProcessorParameter( "kappa",
      "global factor to Highland formula",
      _kappa, static_cast <double>(1.0)); // 1.0 means HL as is, 1.2 means 20% additional scattering

  registerProcessorParameter( "targetthick",
      "target thickness in um",
      _targetthick, static_cast <double>(0.0)); 
}


//------------------------------------------------------------------------------
void EUTelAlignGBL::init() {
  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  streamlog_out( ERROR2 ) << "Marlin was not built with GEAR support." << endl;
  streamlog_out( ERROR2 ) << "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if( Global::GEAR == 0x0 ) {
    streamlog_out( ERROR2) << "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _histogramSwitch = true;

  // Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

  _planePosition = new double[_nTelPlanes]; // z pos

  for(int i = 0; i < _nTelPlanes; i++) {

    //_planeID[i]=_siPlanesLayerLayout->getID(i);
    _planePosition[i]=_siPlanesLayerLayout->getLayerPositionZ(i);
    //_planeThickness[i]=_siPlanesLayerLayout->getLayerThickness(i);
    //_planeX0[i]=_siPlanesLayerLayout->getLayerRadLength(i);
    //_planeResolution[i] = _siPlanesLayerLayout->getSensitiveResolution(i);
  }

  streamlog_out(MESSAGE4) << "assumed beam energy " << _eBeam << " GeV" <<  endl;
  //lets guess the number of planes

    // the number of planes is got from the GEAR description and is
    // the sum of the telescope reference planes and the DUT (if
    // any)
    _nPlanes = _siPlanesParameters->getSiPlanesNumber();
    if( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) {
      ++_nPlanes;
    }

  // an associative map for getting also the sensorID ordered
  map< double, int > sensorIDMap;

  //lets create an array with the z positions of each layer
  for(  int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    _siPlaneZPosition.push_back(_siPlanesLayerLayout->getLayerPositionZ(iPlane));
    sensorIDMap.insert( make_pair( _siPlanesLayerLayout->getLayerPositionZ(iPlane), _siPlanesLayerLayout->getID(iPlane) ) );
  }

  if  ( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) {
    _siPlaneZPosition.push_back(_siPlanesLayerLayout->getDUTPositionZ());
    sensorIDMap.insert( make_pair( _siPlanesLayerLayout->getDUTPositionZ(),  _siPlanesLayerLayout->getDUTID() ) ) ;
  }

  //lets sort the array with increasing z
  sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());

  cout << "planes sorted along z:" << endl;
  for( size_t i = 0; i < _siPlaneZPosition.size(); i++ ) {
    cout << i << "  z " << _siPlaneZPosition[i] << endl;
  }
  // the user is giving sensor ids for the planes to be excluded.
  // These sensor ids have to be converted to a local index
  // according to the planes positions along the z axis.

  cout << "FixedPlanes " << _FixedPlanes_sensorIDs.size() << endl;

  for( size_t i = 0; i < _FixedPlanes_sensorIDs.size(); i++ ) {
    cout << "plane " << _FixedPlanes_sensorIDs[i] << " fixed\n";
    map< double, int >::iterator iter = sensorIDMap.begin();
    int counter = 0;
    while ( iter != sensorIDMap.end() ) {
      if( iter->second == _FixedPlanes_sensorIDs[i] ) {
	_FixedPlanes.push_back(counter);
	break;
      }
      ++iter;
      ++counter;
    }
  }

  cout << "ExcludedPlanes " << _excludePlanes_sensorIDs.size() << endl;
  for( size_t i = 0; i < _excludePlanes_sensorIDs.size(); i++ ) {
    cout << "plane " << _excludePlanes_sensorIDs[i] << " excluded\n";
    map< double, int >::iterator iter = sensorIDMap.begin();
    int counter = 0;
    while ( iter != sensorIDMap.end() ) {
      if( iter->second == _excludePlanes_sensorIDs[i]) {
	_excludePlanes.push_back(counter);
	break;
      }
      ++iter;
      ++counter;
    }
  }

  // strip from the map the sensor id already sorted.
  map< double, int >::iterator iter = sensorIDMap.begin();
  int counter = 0;
  while ( iter != sensorIDMap.end() ) {
    bool excluded = false;
    for( size_t i = 0; i < _excludePlanes.size(); i++) {
      if(_excludePlanes[i] == counter) {
	// printf("excludePlanes %2d of %2d (%2d) \n", i, _excludePlanes_sensorIDs.size(), counter );
	excluded = true;
	break;
      }
    }
    if(!excluded)
      _orderedSensorID_wo_excluded.push_back( iter->second );

    _orderedSensorID.push_back( iter->second );

    ++iter;
    ++counter;
  }

  //consistency
  if( (int)_siPlaneZPosition.size() != _nPlanes ) {
    streamlog_out( ERROR2 ) << "the number of detected planes is " << _nPlanes << " but only " << _siPlaneZPosition.size() << " layer z positions were found!"  << endl;
    exit(-1);
  }

#endif

  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;
  _printEventCounter = 0;

  // Initialize number of excluded planes
  _nExcludePlanes = _excludePlanes.size();

  streamlog_out( MESSAGE2 ) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << endl;

  // Initialise Mille statistics:
  _nMilleTracks = 0;

  _telescopeResolX = new double[_nPlanes];
  _telescopeResolY = new double[_nPlanes];
  _telescopeResolZ = new double[_nPlanes];

  //  printf("print resolution  X: %2d, Y: %2d, Z: %2d \n", _resolutionX.size(), _resolutionY.size(), _resolutionZ.size() );

  // booking histograms
  bookHistos();

  streamlog_out( MESSAGE2 ) << "Initialising Mille..." << endl;

  unsigned int reserveSize = 8000;
  milleAlignGBL = new gbl::MilleBinary( _binaryFilename, reserveSize );

  streamlog_out( MESSAGE2 ) << "The filename for the binary file is: " << _binaryFilename.c_str() << endl;


  /*if( _distanceMaxVec.size() > 0 ) {
    if(_distanceMaxVec.size() !=  static_cast<unsigned int>(_nPlanes )-1 ) {
    streamlog_out( WARNING2 ) << "Consistency check of the DistanceMaxVec array failed. Its size is different compared to the number of planes! Will now use _distanceMax for each pair of planes." << endl;
    _distanceMaxVec.clear();
    for(int i = 0; i < _nPlanes-1; i++) {
    _distanceMaxVec.push_back(_distanceMax);
    }
    }
    }
    else {
    _distanceMaxVec.clear();
    for(int i = 0; i < _nPlanes-1; i++) {
    _distanceMaxVec.push_back(_distanceMax);
    }
    }
    */

  // apply correction to cut if is not first alignment step
  if(!_IsFirstAlignStep){
    _triCut = _triCut*6./_eBeam;
    _driCut = _driCut*6./_eBeam;
    _sixCut = _sixCut*6./_eBeam;
  }

  if(_alignModeString.compare("XYShiftsRotZ") == 0 ) {
	_alignMode = Utility::alignMode::XYShiftsRotZ;
  } else if( _alignModeString.compare("XYShifts") == 0 ) {
	_alignMode = Utility::alignMode::XYShifts;
  } else if( _alignModeString.compare("XYZShiftsRotZ") == 0 ) {
	_alignMode = Utility::alignMode::XYZShiftsRotZ;
  } else {
	streamlog_out(ERROR) << "The chosen AlignMode: '" << _alignModeString << "' is invalid. Please correct your steering template and retry!" << std::endl;
	throw InvalidParameterException("AlignMode");
  }


  streamlog_out( MESSAGE2 ) << "end of init" << endl;
}


//------------------------------------------------------------------------------
void EUTelAlignGBL::processRunHeader( LCRunHeader * rdr ) {

  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() ) ;

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.

  if( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) {
    streamlog_out( ERROR2 ) << "Error during the geometry consistency check: " << endl;
    streamlog_out( ERROR2 ) << "The run header says the GeoID is " << header->getGeoID() << endl;
    streamlog_out( ERROR2 ) << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() << endl;

#ifdef EUTEL_INTERACTIVE
    string answer;
    while (true) {
      streamlog_out( ERROR2 ) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << endl;
      cin >> answer;
      // put the answer in lower case before making the comparison.
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if( answer == "q" ) {
	exit(-1);
      } else if( answer == "c" ) {
	break;
      }
    }
#endif

  }

  // increment the run counter
  ++_iRun;
}

//------------------------------------------------------------------------------
Eigen::Matrix<double,5,5> Jac55new( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  Eigen::Matrix<double,5,5> jac = Eigen::Matrix<double,5,5>::Identity();
  //jac.UnitMatrix();
  jac(3,1) = ds; // x = x0 + xp * ds
  jac(4,2) = ds; // y = y0 + yp * ds
  return jac;
}


//------------------------------------------------------------------------------
void EUTelAlignGBL::processEvent( LCEvent * event ) {

  if( _iEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE2 ) << "Processing event "
      << setw(6) << setiosflags(ios::right)
      << event->getEventNumber() << " in run "
      << setw(6) << setiosflags(ios::right)
      << event->getRunNumber()
      << ", currently having "
      << _nMilleTracks << " tracks "
      << endl;
  }

  if( _nMilleTracks > _maxTrackCandidatesTotal ) {
    throw StopProcessingException(this);
  }

  // fill resolution arrays
  for( int help = 0; help < _nPlanes; help++ ) {
    _telescopeResolX[help] = _telescopeResolution;
    _telescopeResolY[help] = _telescopeResolution;
  }

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

  if( evt->getEventType() == kEORE ) {
    streamlog_out( DEBUG2 ) << "EORE found: nothing else to do." << endl;
    return;
  }

  std::vector<std::vector<EUTelAlignGBL::HitsInPlane> > _hitsArray(_nPlanes - _nExcludePlanes, std::vector<EUTelAlignGBL::HitsInPlane>() );
  std::vector<int> indexconverter (_nPlanes,-1);

  //  if( _nExcludePlanes > 0 )
{
	int icounter = 0;
	for( int i = 0; i < _nPlanes; i++ ) {
		int excluded = 0; //0 - not excluded, 1 - excluded
		if( _nExcludePlanes > 0 ) {
			for( int helphelp = 0; helphelp < _nExcludePlanes; helphelp++ ) {
				if( i == _excludePlanes[helphelp] ) {
					excluded = 1;
					break;//leave the for loop
				}
			}
		}
		if( excluded == 1 )
			indexconverter[i] = -1;
		else {
			indexconverter[i] = icounter;
			icounter++;
		}
	}
}


    for( size_t i = 0; i < _hitCollectionName.size(); i++ ) {

      LCCollection* collection;
      try {
	collection = event->getCollection(_hitCollectionName[i]);
      } catch (DataNotAvailableException& e) {
	//	streamlog_out( WARNING2 ) << "No input collection "
	//			  << _hitCollectionName[i] << " found for event "
	//			  << event->getEventNumber()
	//			  << " in run " << event->getRunNumber() << endl;
	throw SkipEventException(this);
      }
      int layerIndex = -1;
      HitsInPlane hitsInPlane;

	// loop over all hits in collection:

        if(_printEventCounter < 100) streamlog_out(DEBUG2) << "Event " << event->getEventNumber() << " contains " 
	  << collection->getNumberOfElements() << " hits" << endl;
        nAllHitHistoGBLAlign->fill(collection->getNumberOfElements());

	for( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ ) {


	  TrackerHitImpl * hit = static_cast<TrackerHitImpl*> ( collection->getElementAt(iHit) );

	  LCObjectVec clusterVector = hit->getRawHits();

	  double minDistance =  numeric_limits< double >::max();
	  double * hitPosition = const_cast<double * > (hit->getPosition());

	  for(  int i = 0; i < (int)_siPlaneZPosition.size(); i++ ) {
	    double distance = std::abs( hitPosition[2] - _siPlaneZPosition[i] );
	    if( distance < minDistance ) {
	      minDistance = distance;
	      layerIndex = i;
	    }
	  }
	  if( minDistance > 30 /* mm */ ) {
	    // advise the user that the guessing wasn't successful
	    streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
	      "Please check the consistency of the data with the GEAR file" << endl;
	  }

	  // Getting positions of the hits.
	  // ------------------------------
	  hitsInPlane.measuredX = hit->getPosition()[0]; // mm !! (ARGH, PLEASE DO NOT USE [um])
	  hitsInPlane.measuredY = hit->getPosition()[1]; // 
	  hitsInPlane.measuredZ = hit->getPosition()[2]; // 

	  if(_iEvt < 5) streamlog_out( DEBUG3 ) << "hit x = " << hit->getPosition()[0] << endl;
	  if(_iEvt < 5) streamlog_out( DEBUG3 ) << "hit y = " << hit->getPosition()[1] << endl;
	  if(_iEvt < 5) streamlog_out( DEBUG3 ) << "hit z = " << hit->getPosition()[2] << endl;


	  //printf("hit %5d of %5d , at %-8.3f %-8.3f %-8.3f, %5d %5d \n", iHit , collection->getNumberOfElements(), hitsInPlane.measuredX*1E-3, hitsInPlane.measuredY*1E-3, hitsInPlane.measuredZ*1E-3, indexconverter[layerIndex], layerIndex );

	  if( indexconverter[layerIndex] != -1 )
	    _hitsArray[indexconverter[layerIndex]].push_back( hitsInPlane );

	} // end loop over all hits in collection

    }//loop over hits

  /*
     streamlog_out( MESSAGE2 ) << "hits per plane: ";
     for( size_t i = 0; i < _hitsArray.size(); i++ )
     streamlog_out( MESSAGE2 ) << _hitsArray[i].size() << "  ";
     streamlog_out( MESSAGE2 ) << endl;
     */
  // mode 1 or 3: tracks are read in

  // DP: correlate telescope hits from different planes

  int iA = indexconverter[0]; // plane 0

  if( iA >= 0 ) { // not excluded
    for( size_t jA = 0; jA < _hitsArray[iA].size(); jA++ ) { // hits in plane

      int iB = indexconverter[1]; // plane 1
      if( iB >= 0 ) { // not excluded
	for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
	  double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
	  double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
	  dx01HistGBLAlign->fill( dx*1e3 );
	  dy01HistGBLAlign->fill( dy*1e3 );
	}//loop hits jB
      }//iB valid

      iB = indexconverter[2]; // plane 2
      if( iB >= 0 ) { // not excluded
	for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
	  double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
	  double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
	  dx02HistGBLAlign->fill( dx*1e3 );
	  dy02HistGBLAlign->fill( dy*1e3 );
	}//loop hits jB
      }//iB valid

      iB = indexconverter[3]; // plane 3
      if( iB >= 0 ) { // not excluded
	for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
	  double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
	  double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
	  dx03HistGBLAlign->fill( dx*1e3 );
	  dy03HistGBLAlign->fill( dy*1e3 );
	}//loop hits jB
      }//iB valid

      iB = indexconverter[4]; // plane 4
      if( iB >= 0 ) { // not excluded
	for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
	  double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
	  double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
	  dx04HistGBLAlign->fill( dx*1e3 );
	  dy04HistGBLAlign->fill( dy*1e3 );
	}//loop hits jB
      }//iB valid

      iB = indexconverter[5]; // plane 5
      if( iB >= 0 ) { // not excluded
	for( size_t jB = 0; jB < _hitsArray[iB].size(); jB++ ) {
	  double dx = _hitsArray[iB][jB].measuredX - _hitsArray[iA][jA].measuredX;
	  double dy = _hitsArray[iB][jB].measuredY - _hitsArray[iA][jA].measuredY;
	  dx05HistGBLAlign->fill( dx*1e3 );
	  dy05HistGBLAlign->fill( dy*1e3 );
	}//loop hits jB
      }//iB valid

    }//loop hits jA
  }// iA valid

  // triplets 0-1-2:
  // driplets 3-4-5:

  // i is plane index
  // j is hit index
  // k is triplet index

  int i0 = indexconverter[0]; // plane 0
  int i1 = indexconverter[1]; // plane 1
  int i2 = indexconverter[2]; // plane 2

  int i3 = indexconverter[3]; // plane 3
  int i4 = indexconverter[4]; // plane 4
  int i5 = indexconverter[5]; // plane 5


  if( i0*i1*i2*i3*i4*i5 >= 0 ) { // not excluded


    int ntri = 0;
    double xmA[99];
    double ymA[99];
    double zmA[99];
    double sxA[99];
    double syA[99];
    int hts[6][99]; // 6 planes

    double p = _eBeam; // beam momentum

    for( size_t j0 = 0; j0 < _hitsArray[i0].size(); j0++ ) { // hits in plane

      for( size_t j2 = 0; j2 < _hitsArray[i2].size(); j2++ ) {

	double dx02 = _hitsArray[i2][j2].measuredX - _hitsArray[i0][j0].measuredX;
	double dy02 = _hitsArray[i2][j2].measuredY - _hitsArray[i0][j0].measuredY;
	double dz02 = _hitsArray[i2][j2].measuredZ - _hitsArray[i0][j0].measuredZ;

	double avx = 0.5 * ( _hitsArray[i0][j0].measuredX + _hitsArray[i2][j2].measuredX ); // average
	double avy = 0.5 * ( _hitsArray[i0][j0].measuredY + _hitsArray[i2][j2].measuredY ); // average
	double avz = 0.5 * ( _hitsArray[i0][j0].measuredZ + _hitsArray[i2][j2].measuredZ ); // average

	double tx = dx02 / dz02; // angle theta x
	double ty = dy02 / dz02; // angle theta x

	// middle plane 1 for triplets:

	for( size_t j1 = 0; j1 < _hitsArray[i1].size(); j1++ ) {

	  // triplet residual:

	  double zs = _hitsArray[i1][j1].measuredZ - avz;
	  double xs = avx + tx * zs; // doublet extrapolation at 1
	  double ys = avy + ty * zs;

	  double zDUT = _hitsArray[i0][j0].measuredZ + 2.0 * avz + 0.5*(_planePosition[3] - _planePosition[2]); // halfway b/wp2 and p3
	  double xA = _hitsArray[i0][j0].measuredX + tx * zDUT; // doublet extrapolation at zDUT
	  double yA = _hitsArray[i0][j0].measuredY + ty * zDUT; // doublet extrapolation at zDUT

	  double dx = _hitsArray[i1][j1].measuredX - xs; // residual 
	  double dy = _hitsArray[i1][j1].measuredY - ys;

          //if(_iEvt < 5) streamlog_out( WARNING2 ) << "dx triplet = " << dx << endl;

	  if( abs(dy) < _triCut ) tridxHistGBLAlign->fill( dx*1e3 );
	  if( abs(dx) < _triCut ) tridyHistGBLAlign->fill( dy*1e3 );

	  if( abs(dx) < _triCut  && abs(dy) < _triCut && abs(dx02) < _slopeCut*dz02 && abs(dy02) < _slopeCut*dz02 ) {

	    tridxcHistGBLAlign->fill( dx*1e3 );
	    tridycHistGBLAlign->fill( dy*1e3 );

            // poor man's isolation
            bool cont = false;
	    for( size_t j1a = j1+1; j1a < _hitsArray[i1].size(); j1a++ ) {
	      if( fabs(_hitsArray[i1][j1].measuredX - _hitsArray[i1][j1a].measuredX) < 0.1) cont = true;
	      if( fabs(_hitsArray[i1][j1].measuredY - _hitsArray[i1][j1a].measuredY) < 0.1) cont = true;
	    }
	    if(cont) continue;

	    // apply fiducial cut
            if ( fabs(xA) >  9.) continue;
            if (     -yA  < -4.) continue;

	    if( ntri < 99 ) {
	      xmA[ntri] = avx;
	      ymA[ntri] = avy;
	      zmA[ntri] = avz;
	      sxA[ntri] = tx;
	      syA[ntri] = ty;
	      hts[0][ntri] = j0;
	      hts[1][ntri] = j1;
	      hts[2][ntri] = j2;
	    }
	    ntri++;
	  }//valid triplet

	}//loop hits j1
      }//loop hits j2
    }//loop hits j0

    ntriHistGBLAlign->fill( ntri );

    if( ntri >= 99 ) {
      streamlog_out( WARNING2 ) << "Maximum number of track candidates reached. Maybe further tracks were skipped" << endl;
    }

    // driplets 3-4-5:

    int ndri = 0;
    double xmB[99];
    double ymB[99];
    double zmB[99];
    double sxB[99];
    double syB[99];


    for( size_t j3 = 0; j3 < _hitsArray[i3].size(); j3++ ) { // hits in plane

      for( size_t j5 = 0; j5 < _hitsArray[i5].size(); j5++ ) {

	double dx35 = _hitsArray[i5][j5].measuredX - _hitsArray[i3][j3].measuredX;
	double dy35 = _hitsArray[i5][j5].measuredY - _hitsArray[i3][j3].measuredY;
	double dz35 = _hitsArray[i5][j5].measuredZ - _hitsArray[i3][j3].measuredZ;

	double avx = 0.5 * ( _hitsArray[i3][j3].measuredX + _hitsArray[i5][j5].measuredX ); // average
	double avy = 0.5 * ( _hitsArray[i3][j3].measuredY + _hitsArray[i5][j5].measuredY ); // average
	double avz = 0.5 * ( _hitsArray[i3][j3].measuredZ + _hitsArray[i5][j5].measuredZ ); // average

	double tx = dx35 / dz35; // angle theta x
	double ty = dy35 / dz35; // angle theta x

	// middle plane 4 for triplets:

	for( size_t j4 = 0; j4 < _hitsArray[i4].size(); j4++ ) {

	  // triplet residual:

	  double zs = _hitsArray[i4][j4].measuredZ - avz;
	  double xs = avx + tx * zs; // track at 4
	  double ys = avy + ty * zs;

	  double zDUT = _hitsArray[i5][j5].measuredZ - 2.0 * avz - 0.5*(_planePosition[3] - _planePosition[2]);
	  double xA = _hitsArray[i5][j5].measuredX - tx * zDUT; // doublet extrapolation at zDUT
	  double yA = _hitsArray[i5][j5].measuredY - ty * zDUT; // doublet extrapolation at zDUT

	  double dx = _hitsArray[i4][j4].measuredX - xs;
	  double dy = _hitsArray[i4][j4].measuredY - ys;

          //if(_iEvt < 5) streamlog_out( WARNING2 ) << "dx driplet = " << dx << endl;

	  if( abs(dy) < _driCut ) dridxHistGBLAlign->fill( dx*1e3 );
	  if( abs(dx) < _driCut ) dridyHistGBLAlign->fill( dy*1e3 );
	  
	  // for larger targets, slopeCut on driplet is not reasonable, so relax cut. here by 5 times angular width of 10m Alu at 1 GeV
	  if( abs(dx) < _driCut  && abs(dy) < _driCut && abs(dx35) < (_slopeCut+0.006)*dz35 && abs(dy35) < (_slopeCut+0.006)*dz35 ) {

	    dridxcHistGBLAlign->fill( dx*1e3 );
	    dridycHistGBLAlign->fill( dy*1e3 );

	    // I guess that this is poor man's isolation
            bool cont = false;
	    for( size_t j4a = j4+1; j4a < _hitsArray[i4].size(); j4a++ ) {
	      if( fabs(_hitsArray[i4][j4].measuredX - _hitsArray[i4][j4a].measuredX) < 0.1) cont = true;
	      if( fabs(_hitsArray[i4][j4].measuredY - _hitsArray[i4][j4a].measuredY) < 0.1) cont = true;
	    }
	    if(cont) continue;

	    // apply fiducial cut
            if ( fabs(xA) >  9.0) continue;
            if (     -yA  < -4.0) continue;

	    if( ndri < 99 ) {
	      xmB[ndri] = avx;
	      ymB[ndri] = avy;
	      zmB[ndri] = avz;
	      sxB[ndri] = tx;
	      syB[ndri] = ty;
	      hts[3][ndri] = j3;
	      hts[4][ndri] = j4;
	      hts[5][ndri] = j5;
	    }
	    ndri++;

	  }//valid driplet

	}//loop hits j4
      }//loop hits j5
    }//loop hits j3

    ndriHistGBLAlign->fill( ndri );

    if( ndri >= 99 ) {
      streamlog_out( WARNING2 ) << "Maximum number of track candidates reached. Maybe further tracks were skipped" << endl;
    }

    // match triplets A to driplets B:

    int nm = 0;


    for( int kA = 0; kA < ntri && kA < 99; ++kA ) { // i = A

      int j2 = hts[2][kA];

      for( int kB = 0; kB < ndri && kB < 99; ++kB ) { // j = B = REF

	int j3 = hts[3][kB];

	double zmid = 0.5 * ( _hitsArray[i2][j2].measuredZ + _hitsArray[i3][j3].measuredZ );

	double xA = xmA[kA] + sxA[kA] * ( zmid - zmA[kA] ); // B at zmid
	double yA = ymA[kA] + syA[kA] * ( zmid - zmA[kA] );

	double xB = xmB[kB] + sxB[kB] * ( zmid - zmB[kB] ); // B at zmid
	double yB = ymB[kB] + syB[kB] * ( zmid - zmB[kB] );

	double dx = xB - xA; // matching residual
	double dy = yB - yA;

	if( abs(dy) < _sixCut ) sixdxHistGBLAlign->fill( dx*1e3 );
	if( abs(dx) < _sixCut ) sixdyHistGBLAlign->fill( dy*1e3 );

	if( abs(dx) < _sixCut  && abs(dy) < _sixCut ) { // triplet-driplet match


	  double kx = sxB[kB] - sxA[kA]; //kink
	  double ky = syB[kB] - syA[kA];

	  sixkxHistGBLAlign->fill( kx*1E3 );
	  sixkyHistGBLAlign->fill( ky*1E3 );

	  // GBL point vector for the trajectory, all in [mm] !!
	  // GBL with triplet A as seed
	  std::vector<gbl::GblPoint> traj_points;

	  // build up trajectory:
	  std::vector<unsigned int> ilab; // 0-5 = telescope, 6 = DUT, 7 = REF
	  vector<double> sPoint;

	  // the arc length at the first measurement plane is 0.
	  double s = 0;

	  Eigen::Matrix2d proL2m = Eigen::Matrix2d::Identity();

	  // this depends on threshold, energy, dz, and number of iterations done ...
	  // We make an accurate, heuristic  guess
	  double resx = -1.; // [mm] telescope initial resolution for ONLY PREALIGNED t'scope! this is not 3.24 !
	  double resy = -1.; // [mm] telescope initial resolution

	  double distplane = _planePosition[1] - _planePosition[0];
          if(_iEvt < 5) streamlog_out( MESSAGE2 ) << "distplane = " << distplane << endl;

	  if( distplane > 100. ) {
	    resx = 100. - p*8;
	    resy = resx; 
	    if (p > 11) {
	      resx = 8;
	      resy = resx;
	    }
	  }

	  if( distplane < 30. ) {
	    resx = 20. - p*1.5; // if only tracks with prob(chi2,ndf) > 0.001 are passed to Mille
	    resy = resx; 
	  }

	  if( distplane > 30. && distplane < 100. ) {
	    resx = 50 - p*3.8; 
	    resy = resx; 
	    if (p > 10) {
	      resx = 6.;
	      resy = resx;
	    }
	  }

	  if(!_IsFirstAlignStep){
	    // 2nd iteration 20
	    if( distplane < 30. ) {
	      resx = 4.4 - p*0.1;
	      resy = resx;
	    }
	    // 2nd iteration 150
	    if( distplane > 100. ) {
	      resx = 18 - p*1.5;
	      resy = resx;
	    }
	  }

          resx = resx/1000.; // finally convert to [mm]
          resy = resy/1000.;
	  if(_printEventCounter < 10) streamlog_out( MESSAGE2 ) << "res x = " << resx << endl;

	  Eigen::Vector2d measPrec; // precision = 1/resolution^2
	  measPrec[0] = 1.0 / resx / resx;
	  measPrec[1] = 1.0 / resy / resy;

	  // scatter:
	  Eigen::Vector2d scat = Eigen::Vector2d::Zero(); //mean is zero

	  double epsSi = 55e-3 / 93.66 + 0.050 / 286.6; // Si + Kapton
	  double epsAir = -1.; // define later when dz is known
	  double epsAlu = _targetthick/88.97; // Alu target
	  double sumeps = 0.0;


	  // loop over all scatterers first to calculate sumeps, needed for log correctionin Highland
	  for( int ipl = 0; ipl < 6; ++ipl ){
	    sumeps += epsSi;
	    if( ipl < 5) {
	      distplane = _planePosition[ipl+1] - _planePosition[ipl];
	      epsAir =   distplane  / 304200.; 
	      sumeps += epsAir;
	    }

	  }
	  sumeps += epsAlu;
	  // done with calculating sum eps

	  // Paper showed HL predicts too high angle, at least for biased measurement. 
	  double tetSi = _kappa*0.0136 * sqrt(epsSi) / p * ( 1 + 0.038*std::log(sumeps) );
	  double tetAir = -1.;


	  Eigen::Vector2d wscatSi;
	  wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
	  wscatSi[1] = 1.0 / ( tetSi * tetSi );

	  Eigen::Vector2d wscatAir;


	  Eigen::Matrix2d alDer2; // alignment derivatives
	  alDer2(0,0) = 1.0; // dx/dx GBL sign convetion
	  alDer2(1,0) = 0.0; // dy/dx
	  alDer2(0,1) = 0.0; // dx/dy
	  alDer2(1,1) = 1.0; // dy/dy

	  Eigen::Matrix<double,2,3> alDer3; // alignment derivatives
	  alDer3(0,0) = 1.0; // dx/dx
	  alDer3(1,0) = 0.0; // dy/dx
	  alDer3(0,1) = 0.0; // dx/dy
	  alDer3(1,1) = 1.0; // dy/dy

	  Eigen::Matrix<double, 2,4> alDer4; // alignment derivatives
	  alDer4(0,0) = 1.0; // dx/dx
	  alDer4(1,0) = 0.0; // dy/dx
	  alDer4(0,1) = 0.0; // dx/dy
	  alDer4(1,1) = 1.0; // dy/dy
	  alDer4(0,3) = sxA[kA]; // dx/dz
	  alDer4(1,3) = syA[kA]; // dx/dz

	  // telescope planes 0-5:

	  double rx[6];
	  double ry[6];

	  int jhit = hts[0][kA];
	  //double zprev = _hitsArray[0][jhit].measuredZ; // first plane, including any pre-alignment
	  double step = .0;
	  unsigned int iLabel;

	  for( int ipl = 0; ipl < 6; ++ipl ) {


	    if( ipl < 3 )
	      jhit = hts[ipl][kA];
	    else
	      jhit = hts[ipl][kB];

	    double zz = _hitsArray[ipl][jhit].measuredZ; // [mm]

	    //double step = zz - zprev;
	    //
	  // transport matrix in (q/p, x', y', x, y) space
	    auto jacPointToPoint = Jac55new( step );
	    gbl::GblPoint *point = new gbl::GblPoint( jacPointToPoint );
	    s += step;
	    //zprev = zz;

	    // extrapolate triplet vector A to each plane:
	    double dz = zz - zmA[kA];
	    double xs = (xmA[kA] + sxA[kA] * dz); // Ax at plane
	    double ys = (ymA[kA] + syA[kA] * dz); // Ay at plane
	    if(_printEventCounter < 10) std::cout << "dz = " << dz << "   xs = " << xs << "   ys = " << ys << std::endl;

	    rx[ipl] = (_hitsArray[ipl][jhit].measuredX - xs); // resid hit-triplet, in micrometer ...
	    ry[ipl] = (_hitsArray[ipl][jhit].measuredY - ys); // resid
	  
		Eigen::Vector2d meas;
	    meas[0] = rx[ipl]; // fill meas vector for GBL
	    meas[1] = ry[ipl];

	    point->addMeasurement( proL2m, meas, measPrec );

	    point->addScatterer( scat, wscatSi );


	    if( _alignMode == Utility::alignMode::XYShifts ) { // only x and y shifts
	      // global labels for MP:
	      std::vector<int> globalLabels(2);
	      globalLabels[0] = 1 + 2*ipl;
	      globalLabels[1] = 2 + 2*ipl;
	      point->addGlobals( globalLabels, alDer2 ); // for MillePede alignment
	    }
	    else if( _alignMode == Utility::alignMode::XYShiftsRotZ ) { // with rot
	      std::vector<int> globalLabels(3);
	      globalLabels[0] = 1 + 3*ipl; // x
	      globalLabels[1] = 2 + 3*ipl; // y
	      globalLabels[2] = 3 + 3*ipl; // rot
	      //alDer3[0][2] = -_hitsArray[ipl][jhit].measuredY; // dx/dphi
	      //alDer3[1][2] =  _hitsArray[ipl][jhit].measuredX; // dy/dphi
	      alDer3(0,2) = -ys; // dx/dphi
	      alDer3(1,2) =  xs; // dy/dphi
	      point->addGlobals( globalLabels, alDer3 ); // for MillePede alignment
	    }
	    else if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) { // with rot and z shift
	      std::vector<int> globalLabels(4);
	      globalLabels[0] = 1 + 4*ipl;
	      globalLabels[1] = 2 + 4*ipl;
	      globalLabels[2] = 3 + 4*ipl;
	      globalLabels[3] = 4 + 4*ipl; // z
	      //alDer4[0][2] = -_hitsArray[ipl][jhit].measuredY; // dx/dphi
	      //alDer4[1][2] =  _hitsArray[ipl][jhit].measuredX; // dy/dphi
	      alDer4(0,2) = -ys; // dx/dphi
	      alDer4(1,2) =  xs; // dy/dphi
	      point->addGlobals( globalLabels, alDer4 ); // for MillePede alignment
	    }

	    sPoint.push_back( s );
	    iLabel = sPoint.size();
	    ilab.push_back(iLabel);
	    traj_points.push_back(*point);

	    delete point;


	    if( ipl < 5) {
	      distplane = _planePosition[ipl+1] - _planePosition[ipl];

	      //std::cout << " --- Add air --- " << std::endl;
	      step = 0.21*distplane; // in [mm]
	      epsAir =   0.5*distplane  / 304200.;  // in [mm]
	      tetAir = _kappa*0.0136 * sqrt(epsAir) / p * ( 1 + 0.038*std::log(sumeps) );

	      wscatAir[0] = 1.0 / ( tetAir * tetAir ); // weight
	      wscatAir[1] = 1.0 / ( tetAir * tetAir ); 

	      gbl::GblPoint * point = new gbl::GblPoint( Jac55new( step ) );
	      point->addScatterer( scat, wscatAir );

	      s += step;
	      traj_points.push_back(*point);
	      sPoint.push_back( s );
	      delete point;

	      step = 0.58*distplane; // in [mm]
	      if(ipl ==2){ // insert point at centre
		step = step / 2.;
		gbl::GblPoint * pointcentre = new gbl::GblPoint( Jac55new( step ) );

	        if(_targetthick > 0.001) { // at least 1 um target
		  double tetAlu = _kappa*0.0136 * sqrt(epsAlu) / p * ( 1 + 0.038*std::log(sumeps) );

	          Eigen::Vector2d wscatAlu;
	          wscatAlu[0] = 1.0 / ( tetAlu * tetAlu ); // weight
	          wscatAlu[1] = 1.0 / ( tetAlu * tetAlu ); 

		  pointcentre->addScatterer( scat, wscatAlu );
		}
		s += step;
		traj_points.push_back(*pointcentre);
		sPoint.push_back( s );
		//centre_label = sPoint.size();
		delete pointcentre;
		// now again step/2
	      }

	      gbl::GblPoint * point1 = new gbl::GblPoint( Jac55new( step ) );
	      point1->addScatterer( scat, wscatAir );

	      s += step;
	      traj_points.push_back(*point1);
	      sPoint.push_back( s );
	      delete point1;

	      step = 0.21*distplane; // remaing distance to next plane, in [mm]
	    }



	  } // loop over planes

	  // REF:

	  // monitor what we put into GBL:

	  selxHistGBLAlign->fill( xA ); // triplet at mid
	  selyHistGBLAlign->fill( yA );
	  selaxHistGBLAlign->fill( sxA[iA]*1E3 );//track slope
	  selayHistGBLAlign->fill( syA[iA]*1E3 );
	  seldxHistGBLAlign->fill( dx*1E3 ); // triplet-driplet match
	  seldyHistGBLAlign->fill( dy*1E3 );
	  selkxHistGBLAlign->fill( kx*1E3 ); // triplet-driplet kink
	  selkyHistGBLAlign->fill( ky*1E3 );

	  seldx1HistGBLAlign->fill( rx[1]*1E3 ); // triplet interpol
	  seldy1HistGBLAlign->fill( ry[1]*1E3 );
	  seldx3HistGBLAlign->fill( rx[3]*1E3 ); // triplet extrapol
	  seldy3HistGBLAlign->fill( ry[3]*1E3 );
	  seldx4HistGBLAlign->fill( rx[4]*1E3 );
	  seldy4HistGBLAlign->fill( ry[4]*1E3 );
	  seldx5HistGBLAlign->fill( rx[5]*1E3 );
	  seldy5HistGBLAlign->fill( ry[5]*1E3 );


	  double Chi2;
	  int Ndf;
	  double lostWeight;

	  gbl::GblTrajectory traj(traj_points, false); // curvature = false
	  traj.fit( Chi2, Ndf, lostWeight );
	  //traj.getLabels(ilab); // instead pushback sPoint.size() when adding plane

	  // debug:

	  if(_printEventCounter < 10){
	    streamlog_out(DEBUG4) << "traj with " << traj.getNumPoints() << " points:" << endl;
	    for( int ipl = 0; ipl < sPoint.size(); ++ipl ){
	      streamlog_out(DEBUG4) << "  GBL point " << ipl;
	      streamlog_out(DEBUG4) << "  z " << sPoint[ipl]; 
	      streamlog_out(DEBUG4) << endl;
	    }
	    for( int ipl = 0; ipl < 6; ++ipl ){
	      streamlog_out(DEBUG4) << " plane " << ipl << ", lab " << ilab[ipl];
	      streamlog_out(DEBUG4) << " z " << sPoint[ilab[ipl]-1];
	      streamlog_out(DEBUG4) << "  dx " << rx[ipl];
	      streamlog_out(DEBUG4) << "  dy " << ry[ipl];
	      streamlog_out(DEBUG4) << "  chi2 " << Chi2;
	      streamlog_out(DEBUG4) << "  ndf " << Ndf;
	      streamlog_out(DEBUG4) << endl;
	    }

	     streamlog_out(DEBUG2)  << " Is traj valid? " << traj.isValid() << std::endl;
	    traj.printPoints();
	    //traj.printTrajectory();
	    //traj.printData();
	    _printEventCounter++;
	  }

	  //cout << " chi2 " << Chi2 << ", ndf " << Ndf << endl;

	  gblndfHistGBLAlign->fill( Ndf );
	  gblchi2HistGBLAlign->fill( Chi2 );
	  double probchi = TMath::Prob( Chi2, Ndf );
	  gblprbHistGBLAlign->fill( probchi );

	  // bad fits:

	  if( probchi < 0.01 ) {

	    badxHistGBLAlign->fill( xA ); // triplet at DUT
	    badyHistGBLAlign->fill( yA );
	    badaxHistGBLAlign->fill( sxA[iA]*1E3 );
	    badayHistGBLAlign->fill( syA[iA]*1E3 );
	    baddxHistGBLAlign->fill( dx*1E3 ); // triplet-driplet match
	    baddyHistGBLAlign->fill( dy*1E3 );
	    badkxHistGBLAlign->fill( kx*1E3 ); // triplet-driplet kink
	    badkyHistGBLAlign->fill( ky*1E3 );

	    baddx1HistGBLAlign->fill( rx[1]*1E3 ); // triplet interpol
	    baddy1HistGBLAlign->fill( ry[1]*1E3 );
	    baddx3HistGBLAlign->fill( rx[3]*1E3 ); // triplet extrapol
	    baddy3HistGBLAlign->fill( ry[3]*1E3 );
	    baddx4HistGBLAlign->fill( rx[4]*1E3 );
	    baddy4HistGBLAlign->fill( ry[4]*1E3 );
	    baddx5HistGBLAlign->fill( rx[5]*1E3 );
	    baddy5HistGBLAlign->fill( ry[5]*1E3 );

	  }// bad fit

	  else {

	    goodx1HistGBLAlign->fill( rx[1]*1E3 ); // triplet interpol
	    goody1HistGBLAlign->fill( ry[1]*1E3 );

	  } // OK fit

	  // look at fit:

	  Eigen::VectorXd aCorrection;
	  Eigen::MatrixXd aCovariance;

	  double ax[8];
	  double ay[8];
	  unsigned int k = 0;

	  // at plane 0:
          unsigned int ndata = 2;
	  unsigned int ndim = 2;
	  Eigen::VectorXd aResiduals(ndim);
	  Eigen::VectorXd aMeasErrors(ndim);
	  Eigen::VectorXd aResErrors(ndim);
	  Eigen::VectorXd aDownWeights(ndim);

	  int ipos = ilab[0];
	  traj.getResults( ipos, aCorrection, aCovariance );
          traj.getMeasResults( ipos, ndata, aResiduals, aMeasErrors, aResErrors, aDownWeights );

	  //track = q/p, x', y', x, y
	  //        0,   1,  2,  3, 4

	  gblax0HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	  gbldx0HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
	  gblrx0HistGBLAlign->fill( (aResiduals[0])*1e3 ); // residual x [um]
	  gblry0HistGBLAlign->fill( (aResiduals[1])*1e3 ); // residual y [um]
	  gblpx0HistGBLAlign->fill( (aResiduals[0] / aResErrors[0]) );
	  gblpy0HistGBLAlign->fill( (aResiduals[1] / aResErrors[1]) ); 
	  ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	  ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	  k++;

	  ipos = ilab[1];
	  traj.getResults( ipos, aCorrection, aCovariance );
	  gblax1HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	  gbldx1HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
	  gblrx1HistGBLAlign->fill( (rx[1] - aCorrection[3])*1e3 ); // residual x [um]
	  gblry1HistGBLAlign->fill( (ry[1] - aCorrection[4])*1e3 ); // residual y [um]
	  ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	  ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	  k++;

	  ipos = ilab[2];
	  traj.getResults( ipos, aCorrection, aCovariance );
	  gblax2HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	  gbldx2HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
	  gblrx2HistGBLAlign->fill( (rx[2] - aCorrection[3])*1e3 ); // residual x [um]
	  gblry2HistGBLAlign->fill( (ry[2] - aCorrection[4])*1e3 ); // residual y [um]
	  ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	  ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	  k++;

	  ipos = ilab[3];
	  traj.getResults( ipos, aCorrection, aCovariance );
	  gblax3HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	  gbldx3HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
	  gblrx3HistGBLAlign->fill( (rx[3] - aCorrection[3]*1e3) ); // residual x [um]
	  gblry3HistGBLAlign->fill( (ry[3] - aCorrection[4])*1e3 ); // residual y [um]
	  ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	  ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	  k++;

	  ipos = ilab[4];
	  traj.getResults( ipos, aCorrection, aCovariance );
	  gblax4HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	  gbldx4HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
	  gblrx4HistGBLAlign->fill( (rx[4] - aCorrection[3])*1e3 ); // residual x [um]
	  gblry4HistGBLAlign->fill( (ry[4] - aCorrection[4])*1e3 ); // residual y [um]
	  ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	  ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	  k++;

	  ipos = ilab[5];
	  traj.getResults( ipos, aCorrection, aCovariance );
	  gblax5HistGBLAlign->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	  gbldx5HistGBLAlign->fill( aCorrection[3]*1e3 ); // shift x [um]
	  gblrx5HistGBLAlign->fill( (rx[5] - aCorrection[3])*1e3 ); // residual x [um]
	  gblry5HistGBLAlign->fill( (ry[5] - aCorrection[4])*1e3 ); // residual y [um]
	  ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	  ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	  k++;

	  // kinks: 1,2 = tele, 3 = DUT, 4,5 = tele

	  gblkx1HistGBLAlign->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
	  gblkx2HistGBLAlign->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
	  gblkx3HistGBLAlign->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
	  gblkx4HistGBLAlign->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
	  gblkx5HistGBLAlign->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad]

	  // do not pass very bad tracks to mille
	  if(probchi > 0.001) {
	    traj.milleOut( *milleAlignGBL );
	    nm++;
	  }

	} // match

      }//loop kB

    }//loop kA

    //streamlog_out( MESSAGE2 ) << "tracks found: " << _iEvt << ": " << nm << endl;

    nmHistGBLAlign->fill( nm );

    _nMilleTracks += nm;

  }//i3*i4*i5 valid

  // count events:

  ++_iEvt;
  if( isFirstEvent() ) _isFirstEvent = false;

}

//------------------------------------------------------------------------------
void EUTelAlignGBL::end() {

  delete [] _telescopeResolY;
  delete [] _telescopeResolX;
  delete [] _telescopeResolZ;

  // close the output file:
  delete milleAlignGBL;

  // if write the pede steering file
  if( _generatePedeSteerfile ) {

    streamlog_out( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

    ofstream steerFile;
    steerFile.open(_pedeSteerfileName.c_str());

    if( steerFile.is_open()) {

      // find first and last excluded plane
      int firstnotexcl = _nPlanes;
      int lastnotexcl = 0;

      // loop over all planes:

      for( int ipl = 0; ipl < _nPlanes; ipl++) {

	int excluded = 0;

	// loop over excluded planes:

	for( int jpl = 0; jpl < _nExcludePlanes; jpl++ ) {
	  if( ipl == _excludePlanes[jpl] ) excluded = 1;
	}

	if( excluded == 0 && firstnotexcl > ipl ) firstnotexcl = ipl;

	if( excluded == 0 && lastnotexcl < ipl ) lastnotexcl = ipl;

      } // end loop over all planes

      steerFile << "Cfiles" << endl;
      steerFile << _binaryFilename << endl;
      steerFile << endl;

      steerFile << "Parameter" << endl;

      int counter = 0;
      int nfix = 0;

      // loop over all planes:

      for( int ipl = 0; ipl < _nPlanes; ipl++) {

	int excluded = 0; // flag for excluded planes

	// check in list of excluded planes:

	for( int iex = 0; iex < _nExcludePlanes; iex++) {
	  if( ipl == _excludePlanes[iex] )
	    excluded = 1;
	}

	cout << "Plane " << ipl << " exclude = " << excluded << endl;

	if( excluded == 0 ) {

	  bool fixed = false;
	  for( size_t i = 0;i< _FixedPlanes.size(); i++ ) {
	    if( _FixedPlanes[i] == ipl )
	      fixed = true;
	  }

	  // if fixed planes

	  if( fixed || (_FixedPlanes.size() == 0 && (ipl == firstnotexcl || ipl == lastnotexcl) ) ) {
	    nfix++;
	    if( _alignMode == Utility::alignMode::XYShifts ) {
	      steerFile << (counter * 2 + 1) << "  0.0 -1.0" << endl;
	      steerFile << (counter * 2 + 2) << "  0.0 -1.0" << endl;
	    }
	    if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {
	      steerFile << (counter * 3 + 1) << "  0.0 -1.0" << endl; // fix x
	      steerFile << (counter * 3 + 2) << "  0.0 -1.0" << endl; // fix y
	      steerFile << (counter * 3 + 3) << "  0.0 -1.0" << endl; // fix rot
	    }
	    if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
	      steerFile << (counter * 4 + 1) << "  0.0 -1.0" << endl;
	      steerFile << (counter * 4 + 2) << "  0.0 -1.0" << endl;
	      steerFile << (counter * 4 + 3) << "  0.0 -1.0" << endl;
	    }
	  }

	  else {

	    if( _alignMode == Utility::alignMode::XYShifts ) {
	      steerFile << (counter * 2 + 1) << "  0.0  0.0" << endl;
	      steerFile << (counter * 2 + 2) << "  0.0  0.0" << endl;
	    }

	    if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {
	      steerFile << (counter * 3 + 1) << "  0.0  0.0" << endl;
	      steerFile << (counter * 3 + 2) << "  0.0  0.0" << endl;
	      steerFile << (counter * 3 + 3) << "  0.0  0.0" << endl;
	    }

	    if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
	      steerFile << (counter * 4 + 1) << "  0.0  0.0" << endl;
	      steerFile << (counter * 4 + 2) << "  0.0  0.0" << endl;
	      steerFile << (counter * 4 + 3) << "  0.0  0.0" << endl;
	    }

	  }// not fixed

	  // special for z shift:

	  if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {
	    if( ipl == 1 )
	      steerFile << (counter * 4 + 4) << "  0.0 -1.0" << endl;
	    else if( ipl == 4 )
	      steerFile << (counter * 4 + 4) << "  0.0 -1.0" << endl;
	    else
	      steerFile << (counter * 4 + 4) << "  0.0  0.0" << endl;
	  }

	  counter++;

	} // end if plane not excluded

      } // end loop over all planes

      if( nfix < 2 ) {

	if( _alignMode == Utility::alignMode::XYShifts ) {

	  steerFile << "Constraint 0 ! sum dx = 0" << endl;
	  steerFile << " 1  1.0" << endl;
	  steerFile << " 3  1.0" << endl;
	  steerFile << " 5  1.0" << endl;
	  steerFile << " 7  1.0" << endl;
	  steerFile << " 9  1.0" << endl;
	  steerFile << "11  1.0" << endl;

	  steerFile << "Constraint 0 ! sum dy = 0" << endl;
	  steerFile << " 2  1.0" << endl;
	  steerFile << " 4  1.0" << endl;
	  steerFile << " 6  1.0" << endl;
	  steerFile << " 8  1.0" << endl;
	  steerFile << "10  1.0" << endl;
	  steerFile << "12  1.0" << endl;
	}

	if( _alignMode == Utility::alignMode::XYShiftsRotZ ) {

	  steerFile << "Constraint 0 ! sum dx = 0" << endl;
	  steerFile << " 1  1.0" << endl;
	  steerFile << " 4  1.0" << endl;
	  steerFile << " 7  1.0" << endl;
	  steerFile << "10  1.0" << endl;
	  steerFile << "13  1.0" << endl;
	  steerFile << "16  1.0" << endl;

	  steerFile << "Constraint 0 ! sum dy = 0" << endl;
	  steerFile << " 2  1.0" << endl;
	  steerFile << " 5  1.0" << endl;
	  steerFile << " 8  1.0" << endl;
	  steerFile << "11  1.0" << endl;
	  steerFile << "14  1.0" << endl;
	  steerFile << "17  1.0" << endl;

	  steerFile << "Constraint 0 ! sum dphi = 0" << endl;
	  steerFile << " 3  1.0" << endl;
	  steerFile << " 6  1.0" << endl;
	  steerFile << " 9  1.0" << endl;
	  steerFile << "12  1.0" << endl;
	  steerFile << "15  1.0" << endl;
	  steerFile << "18  1.0" << endl;
	}

	if( _alignMode == Utility::alignMode::XYZShiftsRotZ ) {

	  steerFile << "Constraint 0 ! sum dx = 0" << endl;
	  steerFile << " 1  1.0" << endl;
	  steerFile << " 5  1.0" << endl;
	  steerFile << " 9  1.0" << endl;
	  steerFile << "13  1.0" << endl;
	  steerFile << "17  1.0" << endl;
	  steerFile << "21  1.0" << endl;

	  steerFile << "Constraint 0 ! sum dy = 0" << endl;
	  steerFile << " 2  1.0" << endl;
	  steerFile << " 6  1.0" << endl;
	  steerFile << "10  1.0" << endl;
	  steerFile << "14  1.0" << endl;
	  steerFile << "18  1.0" << endl;
	  steerFile << "22  1.0" << endl;

	  steerFile << "Constraint 0 ! sum dphi = 0" << endl;
	  steerFile << " 3  1.0" << endl;
	  steerFile << " 7  1.0" << endl;
	  steerFile << "11  1.0" << endl;
	  steerFile << "15  1.0" << endl;
	  steerFile << "19  1.0" << endl;
	  steerFile << "23  1.0" << endl;

	  steerFile << "Constraint 0 ! sum dz = 0" << endl;
	  steerFile << " 4  1.0" << endl;
	  steerFile << " 8  1.0" << endl;
	  steerFile << "12  1.0" << endl;
	  steerFile << "16  1.0" << endl;
	  steerFile << "20  1.0" << endl;
	  steerFile << "24  1.0" << endl;
	}

      }//nfix < 2

      steerFile << endl;
      steerFile << "! chiscut 5.0 2.5" << endl;
      steerFile << "outlierdownweighting 4" << endl;
      steerFile << "dwfractioncut 0.1" << endl;
      steerFile << endl;
      steerFile << "method inversion 10  0.1" << endl;
      // Use 10 OpenMP threads to process the data:
      steerFile << "threads 10 1" << endl;
      steerFile << endl;
      steerFile << "histprint" << endl;
      steerFile << endl;
      steerFile << "end" << endl;

      steerFile.close();

      streamlog_out( MESSAGE2 ) << "File " << _pedeSteerfileName << " written." << endl;

    }
    else {
      streamlog_out( ERROR2 ) << "Could not open steering file." << endl;
    }

  } // end if write the pede steering file

  streamlog_out( MESSAGE2 ) << endl;
  streamlog_out( MESSAGE2 ) << "Number of tracks used: " << _nMilleTracks << endl;

  streamlog_out( MESSAGE2 ) << endl;
  streamlog_out( MESSAGE2 ) << "Successfully finished" << endl;

}//end

//------------------------------------------------------------------------------
void EUTelAlignGBL::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    streamlog_out( MESSAGE2 ) << "Booking histograms..." << endl;

    //DP:
    
    nAllHitHistoGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "nallhit", 201, -0.5, 200.5 );
    nAllHitHistoGBLAlign->setTitle( "Telescope hits/event;telescope hits;events" );

    dx01HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dx01", 100, -5000, 5000 );
    dx01HistGBLAlign->setTitle( "dx01;x_{1}-x_{0} [um];hit pairs" );

    dy01HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dy01", 100, -5000, 5000 );
    dy01HistGBLAlign->setTitle( "dy01;y_{1}-y_{0} [um];hit pairs" );

    dx02HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dx02", 100, -5000, 5000 );
    dx02HistGBLAlign->setTitle( "dx02;x_{2}-x_{0} [um];hit pairs" );

    dy02HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dy02", 100, -5000, 5000 );
    dy02HistGBLAlign->setTitle( "dy02;y_{2}-y_{0} [um];hit pairs" );

    dx03HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dx03", 100, -5000, 5000 );
    dx03HistGBLAlign->setTitle( "dx03;x_{3}-x_{0} [um];hit pairs" );

    dy03HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dy03", 100, -5000, 5000 );
    dy03HistGBLAlign->setTitle( "dy03;y_{3}-y_{0} [um];hit pairs" );

    dx04HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dx04", 100, -5000, 5000 );
    dx04HistGBLAlign->setTitle( "dx04;x_{4}-x_{0} [um];hit pairs" );

    dy04HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dy04", 100, -5000, 5000 );
    dy04HistGBLAlign->setTitle( "dy04;y_{4}-y_{0} [um];hit pairs" );

    dx05HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dx05", 100, -5000, 5000 );
    dx05HistGBLAlign->setTitle( "dx05;x_{5}-x_{0} [um];hit pairs" );

    dy05HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dy05", 100, -5000, 5000 );
    dy05HistGBLAlign->setTitle( "dy05;y_{5}-y_{0} [um];hit pairs" );

    tridxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridx", 100, -1000, 1000 );
    tridxHistGBLAlign->setTitle( "tridx;x_{1}-x_{02} [um];triplet" );

    tridyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridy", 100, -1000, 1000 );
    tridyHistGBLAlign->setTitle( "tridy;y_{1}-y_{02} [um];triplet" );

    tridxcHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridxc", 100, -100, 100 );
    tridxcHistGBLAlign->setTitle( "tridxc;x_{1}-x_{02} [um];triplet" );

    tridycHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "tridyc", 100, -100, 100 );
    tridycHistGBLAlign->setTitle( "tridyc;y_{1}-y_{02} [um];triplet" );

    ntriHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "ntri", 21, -0.5, 20.5 );
    ntriHistGBLAlign->setTitle( "ntri;triplets;events" );

    dridxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridx", 100, -1000, 1000 );
    dridxHistGBLAlign->setTitle( "dridx;x_{4}-x_{35} [um];triplet" );

    dridyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridy", 100, -1000, 1000 );
    dridyHistGBLAlign->setTitle( "dridy;y_{4}-y_{35} [um];triplet" );

    dridxcHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridxc", 100, -100, 100 );
    dridxcHistGBLAlign->setTitle( "dridxc;x_{4}-x_{35} [um];triplet" );

    dridycHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "dridyc", 100, -100, 100 );
    dridycHistGBLAlign->setTitle( "dridyc;y_{4}-y_{35} [um];triplet" );

    ndriHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "ndri", 21, -0.5, 20.5 );
    ndriHistGBLAlign->setTitle( "ndri;driplets;events" );

    sixdxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixdx", 250, -1000, 1000 );
    sixdxHistGBLAlign->setTitle( "sixdx;x_{A}-x_{B} [um];triplet pairs" );

    sixdyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixdy", 250, -1000, 1000 );
    sixdyHistGBLAlign->setTitle( "sixdy;y_{A}-y_{B} [um];triplet pairs" );

    sixkxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixkx", 100, -25, 25 );
    sixkxHistGBLAlign->setTitle( "sixkx;x kink angle [mrad];triplet pairs" );

    sixkyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sixky", 100, -25, 25 );
    sixkyHistGBLAlign->setTitle( "sixky;y kink angle [mrad];triplet pairs" );

    // GBL:

    selxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selx", 150, -15, 15 );
    selxHistGBLAlign->setTitle( "x at DUT, sel GBL;x [mm];tracks" );

    selyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "sely", 100, -10, 10 );
    selyHistGBLAlign->setTitle( "y at DUT, sel GBL;y [mm];tracks" );

    selaxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selax", 100, -25, 25 );
    selaxHistGBLAlign->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

    selayHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selay", 100, -25, 25 );
    selayHistGBLAlign->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

    seldxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx", 100, -5000, 5000 );
    seldxHistGBLAlign->setTitle( "track match x, sel GBL;#Deltax [#mum];tracks" );

    seldyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy", 100, -5000, 5000 );
    seldyHistGBLAlign->setTitle( "track match y, sel GBL;#Deltay [#mum];tracks" );

    selkxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selkx", 100, -25, 25 );
    selkxHistGBLAlign->setTitle( "kink x, sel GBL;kink x [mrad];tracks" );

    selkyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "selky", 100, -25, 25 );
    selkyHistGBLAlign->setTitle( "kink y, sel GBL;kink y [mrad];tracks" );

    seldx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx1", 100, -1000, 1000 );
    seldx1HistGBLAlign->setTitle( "triplet resid x at 1, sel GBL;#Deltax [#mum];tracks" );

    seldy1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy1", 100, -1000, 1000 );
    seldy1HistGBLAlign->setTitle( "triplet resid y at 1, sel GBL;#Deltay [#mum];tracks" );

    seldx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx3", 100, -1000, 1000 );
    seldx3HistGBLAlign->setTitle( "triplet resid x at 3, sel GBL;#Deltax [#mum];tracks" );

    seldy3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy3", 100, -1000, 1000 );
    seldy3HistGBLAlign->setTitle( "triplet resid y at 3, sel GBL;#Deltay [#mum];tracks" );

    seldx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx4", 100, -1000, 1000 );
    seldx4HistGBLAlign->setTitle( "triplet resid x at 4, sel GBL;#Deltax [#mum];tracks" );

    seldy4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy4", 100, -1000, 1000 );
    seldy4HistGBLAlign->setTitle( "triplet resid y at 4, sel GBL;#Deltay [#mum];tracks" );

    seldx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx5", 100, -5000, 5000 );
    seldx5HistGBLAlign->setTitle( "triplet resid x at 5, sel GBL;#Deltax [#mum];tracks" );

    seldy5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy5", 100, -5000, 5000 );
    seldy5HistGBLAlign->setTitle( "triplet resid y at 5, sel GBL;#Deltay [#mum];tracks" );

    seldx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldx6", 100, -5000, 5000 );
    seldx6HistGBLAlign->setTitle( "triplet resid x at DUT, sel GBL;#Deltax [#mum];tracks" );

    seldy6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "seldy6", 100, -5000, 5000 );
    seldy6HistGBLAlign->setTitle( "triplet resid y at DUT, sel GBL;#Deltay [#mum];tracks" );

    gblndfHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblndf", 16, -0.5, 15.5 );
    gblndfHistGBLAlign->setTitle( "GBL fit NDF;GBL NDF;tracks" );

    gblchi2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblchi2", 100, 0, 100 );
    gblchi2HistGBLAlign->setTitle( "GBL fit chi2;GBL chi2;tracks" );

    gblprbHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblprb", 100, 0, 1 );
    gblprbHistGBLAlign->setTitle( "GBL fit probability;GBL fit probability;tracks" );

    // bad fits:

    badxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badx", 150, -15, 15 );
    badxHistGBLAlign->setTitle( "x at DUT, bad GBL;x [mm];tracks" );

    badyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "bady", 100, -10, 10 );
    badyHistGBLAlign->setTitle( "y at DUT, bad GBL;y [mm];tracks" );

    badaxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badax", 100, -25, 25 );
    badaxHistGBLAlign->setTitle( "track angle x, bad GBL;x angle [mrad];tracks" );

    badayHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baday", 100, -25, 25 );
    badayHistGBLAlign->setTitle( "track angle y, bad GBL;y angle [mrad];tracks" );

    baddxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx", 100, -5000, 5000 );
    baddxHistGBLAlign->setTitle( "track match x, bad GBL;#Deltax [#mum];tracks" );

    baddyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy", 100, -5000, 5000 );
    baddyHistGBLAlign->setTitle( "track match y, bad GBL;#Deltay [#mum];tracks" );

    badkxHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badkx", 100, -25, 25 );
    badkxHistGBLAlign->setTitle( "kink x, bad GBL;kink x [mrad];tracks" );

    badkyHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "badky", 100, -25, 25 );
    badkyHistGBLAlign->setTitle( "kink y, bad GBL;kink y [mrad];tracks" );

    baddx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx1", 100, -1000, 1000 );
    baddx1HistGBLAlign->setTitle( "triplet resid x at 1, bad GBL;#Deltax [#mum];tracks" );

    baddy1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy1", 100, -1000, 1000 );
    baddy1HistGBLAlign->setTitle( "triplet resid y at 1, bad GBL;#Deltay [#mum];tracks" );

    baddx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx3", 100, -1000, 1000 );
    baddx3HistGBLAlign->setTitle( "triplet resid x at 3, bad GBL;#Deltax [#mum];tracks" );

    baddy3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy3", 100, -1000, 1000 );
    baddy3HistGBLAlign->setTitle( "triplet resid y at 3, bad GBL;#Deltay [#mum];tracks" );

    baddx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx4", 100, -1500, 1500 );
    baddx4HistGBLAlign->setTitle( "triplet resid x at 4, bad GBL;#Deltax [#mum];tracks" );

    baddy4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy4", 100, -1500, 1500 );
    baddy4HistGBLAlign->setTitle( "triplet resid y at 4, bad GBL;#Deltay [#mum];tracks" );

    baddx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx5", 100, -3000, 3000 );
    baddx5HistGBLAlign->setTitle( "triplet resid x at 5, bad GBL;#Deltax [#mum];tracks" );

    baddy5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy5", 100, -3000, 3000 );
    baddy5HistGBLAlign->setTitle( "triplet resid y at 5, bad GBL;#Deltay [#mum];tracks" );

    baddx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddx6", 100, -250, 250 );
    baddx6HistGBLAlign->setTitle( "triplet resid x at DUT, bad GBL;#Deltax [#mum];tracks" );

    baddy6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "baddy6", 100, -250, 250 );
    baddy6HistGBLAlign->setTitle( "triplet resid y at DUT, bad GBL;#Deltay [#mum];tracks" );

    // good fits:

    goodx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goodx1", 100, -1000, 1000 );
    goodx1HistGBLAlign->setTitle( "triplet resid x at 1, good GBL;#Deltax [#mum];tracks" );

    goody1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goody1", 100, -1000, 1000 );
    goody1HistGBLAlign->setTitle( "triplet resid y at 1, good GBL;#Deltay [#mum];tracks" );

    goodx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goodx6", 100, -250, 250 );
    goodx6HistGBLAlign->setTitle( "triplet resid x at 6, good GBL;#Deltax [#mum];tracks" );

    goody6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "goody6", 100, -250, 250 );
    goody6HistGBLAlign->setTitle( "triplet resid y at 6, good GBL;#Deltay [#mum];tracks" );

    // look at fit:

    gblax0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax0", 100, -5, 5 );
    gblax0HistGBLAlign->setTitle( "GBL angle at plane 0;x angle at plane 0 [mrad];tracks" );

    gbldx0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx0", 500, -250, 250 );
    gbldx0HistGBLAlign->setTitle( "GBL shift at plane 0;x shift at plane 0 [#mum];tracks" );

    gblrx0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx0", 500, -250, 250 );
    gblrx0HistGBLAlign->setTitle( "GBL resid at plane 0;x resid at plane 0 [#mum];tracks" );

    gblry0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry0", 500, -250, 250 );
    gblry0HistGBLAlign->setTitle( "GBL resid at plane 0;y resid at plane 0 [#mum];tracks" );

    gblpx0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblpx0", 100, -5, 5 );
    gblpx0HistGBLAlign->setTitle( "GBL pull at plane 0;x pull at plane 0;tracks" );

    gblpy0HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblpy0", 100, -5, 5 );
    gblpy0HistGBLAlign->setTitle( "GBL pull at plane 0;y pull at plane 0;tracks" );


    gblax1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax1", 100, -5, 5 );
    gblax1HistGBLAlign->setTitle( "GBL angle at plane 1;x angle at plane 1 [mrad];tracks" );

    gbldx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx1", 500, -250, 250 );
    gbldx1HistGBLAlign->setTitle( "GBL shift at plane 1;x shift at plane 1 [#mum];tracks" );

    gblrx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx1", 500, -250, 250 );
    gblrx1HistGBLAlign->setTitle( "GBL resid at plane 1;x resid at plane 1 [#mum];tracks" );

    gblry1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry1", 500, -250, 250 );
    gblry1HistGBLAlign->setTitle( "GBL resid at plane 1;y resid at plane 1 [#mum];tracks" );


    gblax2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax2", 100, -5, 5 );
    gblax2HistGBLAlign->setTitle( "GBL angle at plane 2;x angle at plane 2 [mrad];tracks" );

    gbldx2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx2", 500, -250, 250 );
    gbldx2HistGBLAlign->setTitle( "GBL shift at plane 2;x shift at plane 2 [#mum];tracks" );

    gblrx2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx2", 500, -250, 250 );
    gblrx2HistGBLAlign->setTitle( "GBL resid at plane 2;x resid at plane 2 [#mum];tracks" );

    gblry2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry2", 500, -250, 250 );
    gblry2HistGBLAlign->setTitle( "GBL resid at plane 2;y resid at plane 2 [#mum];tracks" );


    gblax3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax3", 100, -5, 5 );
    gblax3HistGBLAlign->setTitle( "GBL angle at plane 3;x angle at plane 3 [mrad];tracks" );

    gbldx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx3", 500, -250, 250 );
    gbldx3HistGBLAlign->setTitle( "GBL shift at plane 3;x shift at plane 3 [#mum];tracks" );

    gblrx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx3", 500, -250, 250 );
    gblrx3HistGBLAlign->setTitle( "GBL resid at plane 3;x resid at plane 3 [#mum];tracks" );

    gblry3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry3", 500, -250, 250 );
    gblry3HistGBLAlign->setTitle( "GBL resid at plane 3;y resid at plane 3 [#mum];tracks" );


    gblax4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax4", 100, -5, 5 );
    gblax4HistGBLAlign->setTitle( "GBL angle at plane 4;x angle at plane 4 [mrad];tracks" );

    gbldx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx4", 500, -250, 250 );
    gbldx4HistGBLAlign->setTitle( "GBL shift at plane 4;x shift at plane 4 [#mum];tracks" );

    gblrx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx4", 500, -250, 250 );
    gblrx4HistGBLAlign->setTitle( "GBL resid at plane 4;x resid at plane 4 [#mum];tracks" );

    gblry4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry4", 500, -250, 250 );
    gblry4HistGBLAlign->setTitle( "GBL resid at plane 4;y resid at plane 4 [#mum];tracks" );


    gblax5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax5", 100, -5, 5 );
    gblax5HistGBLAlign->setTitle( "GBL angle at plane 5;x angle at plane 5 [mrad];tracks" );

    gbldx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx5", 500, -250, 250 );
    gbldx5HistGBLAlign->setTitle( "GBL shift at plane 5;x shift at plane 5 [#mum];tracks" );

    gblrx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx5", 500, -250, 250 );
    gblrx5HistGBLAlign->setTitle( "GBL resid at plane 5;x resid at plane 5 [#mum];tracks" );

    gblry5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry5", 500, -250, 250 );
    gblry5HistGBLAlign->setTitle( "GBL resid at plane 5;y resid at plane 5 [#mum];tracks" );


    gblax6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblax6", 100, -5, 5 );
    gblax6HistGBLAlign->setTitle( "GBL angle at DUT;x angle at DUT [mrad];tracks" );

    gbldx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldx6", 100, -250, 250 );
    gbldx6HistGBLAlign->setTitle( "GBL shift at DUT;x shift at DUT [#mum];tracks" );

    gbldy6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gbldy6", 100, -1000, 1000 );
    gbldy6HistGBLAlign->setTitle( "GBL shift at DUT;y shift at DUT [#mum];tracks" );

    gblrx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblrx6", 100, -250, 250 );
    gblrx6HistGBLAlign->setTitle( "GBL resid at DUT;x resid at DUT [#mum];tracks" );

    gblry6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblry6", 100, -250, 250 );
    gblry6HistGBLAlign->setTitle( "GBL resid at DUT;y resid at DUT [#mum];tracks" );


    gblkx1HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx1", 100, -5, 5 );
    gblkx1HistGBLAlign->setTitle( "GBL kink angle at plane 1;plane 1 kink [mrad];tracks" );

    gblkx2HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx2", 100, -5, 5 );
    gblkx2HistGBLAlign->setTitle( "GBL kink angle at plane 2;plane 2 kink [mrad];tracks" );

    gblkx3HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx3", 100, -5, 5 );
    gblkx3HistGBLAlign->setTitle( "GBL kink angle at plane 3;plane 3 kink [mrad];tracks" );

    gblkx4HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx4", 100, -5, 5 );
    gblkx4HistGBLAlign->setTitle( "GBL kink angle at plane 4;plane 4 kink [mrad];tracks" );

    gblkx5HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx5", 100, -5, 5 );
    gblkx5HistGBLAlign->setTitle( "GBL kink angle at plane 5;plane 5 kink [mrad];tracks" );

    gblkx6HistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "gblkx6", 100, -5, 5 );
    gblkx6HistGBLAlign->setTitle( "GBL kink angle at plane 6;plane 6 kink [mrad];tracks" );

    nmHistGBLAlign = AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "nm", 21, -0.5, 20.5 );
    nmHistGBLAlign->setTitle( "track matches;track matches;events" );

  }//try
  catch( lcio::Exception& e ) {

#ifdef EUTEL_INTERACTIVE
    streamlog_out( ERROR2 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
    string answer;
    while ( true ) {
      streamlog_out( ERROR2 ) << "[q]/[c]" << endl;
      cin >> answer;
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if( answer == "q" ) {
	exit(-1);
      }
      else if( answer == "c" )
	_histogramSwitch = false;
      break;
    }
#else
    streamlog_out( WARNING2 ) << "No AIDAProcessor initialized. Continue without histogramming" << endl;

#endif

  }
#endif

}

#endif // USE_GEAR
