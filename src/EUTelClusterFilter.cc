// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
 
// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelClusterFilter.h"
#include "EUTelExceptions.h"
#include "EUTelROI.h"
#include "EUTelMatrixDecoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelClusterFilter::EUTelClusterFilter () :Processor("EUTelClusterFilter") {

  // modify processor description
  _description = "EUTelClusterFilter is a very powerful tool. It allows to select among an input collection of TrackerPulse\n"
    "only the clusters fulfilling a certain set of selection criteria.\n"
    "The user can modify the switch on and off each selection cut and set the proper value for that via the processor parameter.";


  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERPULSE, "InputPulseCollectionName",
                           "This is the input Tracker Pulse collection that should be filtered",
                           _inputPulseCollectionName, string ("cluster"));


  // since v00-00-09 both the status and noise collections become compulsory
  registerInputCollection(LCIO::TRACKERDATA, "NoiseCollectionName",
                          "This is the name of the noise collection.\n"
                          "The presence of this collection in the event is allowing all the noise based selection cuts",
                          _noiseCollectionName, string( "noiseDB" ) );

  registerInputCollection(LCIO::TRACKERRAWDATA,
                          "StatusCollectionName","This is the name of the status collection.\n"
                          "The presence of this collection in the event is allowing all the noise based selection cuts",
                          _statusCollectionName, string( "statusDB" ) );

  // this is the output collection
  registerOutputCollection (LCIO::TRACKERPULSE,"OutputPulseCollectionName",
                            "This is the output Tracker Pulse collection containing the filtered clusters",
                            _outputPulseCollectionName, string ("filteredcluster"));

  // assuming 3 as the number of detectors
  FloatVec minTotalChargeVecExample;
  minTotalChargeVecExample.push_back(90);
  minTotalChargeVecExample.push_back(75);
  minTotalChargeVecExample.push_back(80);

  registerProcessorParameter("ClusterMinTotalCharge", "This is the minimum allowed total charge in to a cluster.\n"
                             "One floating point number for each sensor in the telescope",
                             _minTotalChargeVec, minTotalChargeVecExample);

  // again assuming 3 detectors
  FloatVec minNChargeVecExample;
  minNChargeVecExample.push_back(9); // this is the number of pixels in the cluster
  minNChargeVecExample.push_back(90);
  minNChargeVecExample.push_back(75);
  minNChargeVecExample.push_back(80);

  registerProcessorParameter("ClusterNMinCharge", "This is the minimum charge that a cluster of N pixels has to have.\n"
                             "The first figure has to be the number of pixels to consider in the cluster, \n"
                             "then one float number for each sensor.",
                             _minNChargeVec, minNChargeVecExample);

  FloatVec minNSNRVecExample;
  minNSNRVecExample.push_back(9); // this is the number of pixels in
                                  // the cluster
  minNSNRVecExample.push_back(25.0);
  minNSNRVecExample.push_back(21.0);
  minNSNRVecExample.push_back(20.0);

  registerProcessorParameter("ClusterNMinSNR", "This is the minimum SNR that a cluster of N pixels has to have.\n"
                             "The first figure has to be the number of pixels to consider in the cluster, \n"
                             "then one float number for each sensor. Setting N = 0 is enough to disable the cut.",
                             _minNSNRVec, minNSNRVecExample);

  // again for 3 planes
  FloatVec minNxNChargeVecExample;
  minNxNChargeVecExample.push_back(3);
  minNxNChargeVecExample.push_back(0);
  minNxNChargeVecExample.push_back(0);
  minNxNChargeVecExample.push_back(0);

  registerProcessorParameter("ClusterNxNMinCharge","This is the minimum charge that a cluster of N times N pixels has to have.\n"
                             "The first figure is the subcluster size in pixels (odd number), then one floating number for each \n"
                             "planes. To switch this selection off, set all numbers to zero.",
                             _minNxNChargeVec, minNxNChargeVecExample);

  // again for 3 planes
  FloatVec minNxNSNRVecExample;
  minNxNSNRVecExample.push_back(3);
  minNxNSNRVecExample.push_back(0);
  minNxNSNRVecExample.push_back(0);
  minNxNSNRVecExample.push_back(0);

  registerProcessorParameter("ClusterNxNMinSNR","This is the minimum SNR that a cluster of N times N pixels has to have.\n"
                             "The first figure is the subcluster size in pixels (odd number), then one floating number for each \n"
                             "planes. To switch this selection off, set at least the first number to zero.",
                             _minNxNSNRVec, minNxNSNRVecExample);



  FloatVec minSeedChargeVecExample;
  minSeedChargeVecExample.push_back(20);
  minSeedChargeVecExample.push_back(25);
  minSeedChargeVecExample.push_back(21);

  registerProcessorParameter("SeedMinCharge", "This is the minimum allowed charge that the seed pixel of a cluster has to have.\n"
                             "One floating number for each detector",
                             _minSeedChargeVec, minSeedChargeVecExample);

  FloatVec minSeedSNRVecExample;
  minSeedSNRVecExample.push_back(10);
  minSeedSNRVecExample.push_back(12);
  minSeedSNRVecExample.push_back(14);

  registerProcessorParameter("SeedMinSNR", "This is the minimum allowed SNR that the seed pixel of a cluster has to have.\n"
                             "One floating number for each detector. Set to 0 to disable",
                             _minSeedSNRVec, minSeedSNRVecExample);


  IntVec clusterQualityVecExample;
  clusterQualityVecExample.push_back(0);
  clusterQualityVecExample.push_back(0);
  clusterQualityVecExample.push_back(0);

  registerProcessorParameter("ClusterQuality", "This is the required quality for the cluster.\n"
                             "One integer number for each detector according to ClusterQuality.\n"
                             "Put a negative number to disable the cut",
                             _clusterQualityVec, clusterQualityVecExample);

  FloatVec roiExample;
  roiExample.push_back( -1 );
  roiExample.push_back( 10. );
  roiExample.push_back( 10. );
  roiExample.push_back( 40. );
  roiExample.push_back( 40. );

  registerProcessorParameter("InsideRegion", "Define here ROI's. The first number (integer) is the detector ID.\n"
                             "The other four float are xBotLeft  yBotLeft xTopRight yTopRight.\n"
                             "To disable it, just put a negative number as detector ID.",
                             _tempInsideROI, roiExample, roiExample.size() );

  registerProcessorParameter("OutsideRegion","Define here ROI's. The first number (integer) is the detector ID.\n"
                             "The other four float are xBotLeft  yBotLeft xTopRight yTopRight.\n"
                             "To disable it, just put a negative number as detector ID.",
                             _tempOutsideROI, roiExample, roiExample.size() );


  IntVec minClusterNoVecExample;
  minClusterNoVecExample.push_back(0);
  minClusterNoVecExample.push_back(0);
  minClusterNoVecExample.push_back(0);

  registerProcessorParameter("MinClusterPerPlane", "This is the minimum required number of cluster per plane.\n"
                             "One integer number for each detector. Write 0 to disable the cut",
                             _minClusterNoVec, minClusterNoVecExample);

  IntVec maxClusterNoVecExample;
  maxClusterNoVecExample.push_back(-1);
  maxClusterNoVecExample.push_back(-1);
  maxClusterNoVecExample.push_back(-1);

  registerProcessorParameter("MaxClusterPerPlane", "This is the maximum allowed number of cluster per plane.\n "
                             "One integer number for each detector. Write a negative number to disable the cut",
                             _maxClusterNoVec, maxClusterNoVecExample);

  FloatVec maxClusterNoiseVecExample;
  maxClusterNoiseVecExample.push_back(50.);
  maxClusterNoiseVecExample.push_back(45.);
  maxClusterNoiseVecExample.push_back(48.);

  registerProcessorParameter("MaxClusterNoise", "This is maximum allowed cluster noise.\n"
                             "One floating number for each detector. Write a negative number to disable the cut",
                             _maxClusterNoiseVec, maxClusterNoiseVecExample );

  FloatVec minTotalSNRVecExample;
  minTotalSNRVecExample.push_back(10.);
  minTotalSNRVecExample.push_back(12.);
  minTotalSNRVecExample.push_back(11.);

  registerProcessorParameter("MinTotalClusterSNR", "This is the minimum allow total cluster SNR\n"
                             "One floating number for each detector. Write 0 to disable the cut",
                             _minTotalSNRVec, minTotalSNRVecExample );


  registerProcessorParameter("SameNumberOfHits", "Setting this to true will select only events having the same number \n"
                             "of hits for each plane.",
                             _sameNumberOfHitSwitch, static_cast<bool > (false) );


  registerOptionalParameter("SkipEmptyEvent","If true, a SkipEventException is thrown if after selection\n"
                            "there are no cluster left.",
                            _skipEmptyEvent, static_cast<bool > ( false ) );

  // set the global noise switch to on
  _noiseRelatedCuts = true;




  registerProcessorParameter("DFFNumberOfHits","This is a cut on the number of hit pixels inside the digital fixed frame\n"
                             "cluster algorithm. One cut for each sensor plane.\n",

                             _DFFNHitsCuts, std::vector<int>(9, 0) );




}

void EUTelClusterFilter::init () {

  printParameters ();

  // check and set properly the switches
  // total cluster charge
  if (count_if( _minTotalChargeVec.begin(),  _minTotalChargeVec.end(),  bind2nd(greater<float>(), 0) ) != 0 ) {
    // ok there is at least one sensor for which this cut is on
    _minTotalChargeSwitch = true;
  } else {
    _minTotalChargeSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinTotalChargeSwitch " << _minTotalChargeSwitch << endl;

  // N cluster charge
  if ( !_minNChargeVec.empty() ) {
    if ( _minNChargeVec[0] != 0 ) {
      _minNChargeSwitch = true;
    } else {
      _minNChargeSwitch = false;
    }
  } else {
    _minNChargeSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinNChargeSwitch " << _minNChargeSwitch << endl;

  // N cluster SNR
  if ( (!_minNSNRVec.empty()) && (_minNSNRVec[0] != 0) ) {
    _minNSNRSwitch = true;
  } else {
    _minNSNRSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinNSNRSwitch " << _minNSNRSwitch << endl;

  // NxN cluster charge
  if (  ( count_if( _minNxNChargeVec.begin(), _minNxNChargeVec.end(), bind2nd(greater<float>(),0)) != 0  ) &&
        ( _minNxNChargeVec[0] != 0 ) ) {
    _minNxNChargeSwitch = true;
  } else {
    _minNxNChargeSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinNxNChargeSwitch " << _minNxNChargeSwitch << endl;

  // NxN cluster SNR
  if (  ( count_if( _minNxNSNRVec.begin(), _minNxNSNRVec.end(), bind2nd(greater<float>(),0)) != 0  ) &&
        ( _minNxNSNRVec[0] != 0 ) ) {
    _minNxNSNRSwitch = true;
  } else {
    _minNxNSNRSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinNxNSNRSwitch " << _minNxNSNRSwitch << endl;


  // Seed charge
  if ( count_if( _minSeedChargeVec.begin(), _minSeedChargeVec.end(), bind2nd(greater<float>(), 0) ) != 0 ) {
    _minSeedChargeSwitch = true;
  } else {
    _minSeedChargeSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinSeedChargeSwitch " << _minSeedChargeSwitch << endl;

  // Seed SNR
  if ( count_if( _minSeedSNRVec.begin(), _minSeedSNRVec.end(), bind2nd(greater<float >(), 0) ) != 0 ) {
    _minSeedSNRSwitch = true;
  } else {
    _minSeedSNRSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinSeedSNRSwitch " << _minSeedChargeSwitch << endl;

  // Cluster quality
  if ( count_if(_clusterQualityVec.begin(), _clusterQualityVec.end(), bind2nd(greater_equal<int>(), 0) ) != 0 ) {
    _clusterQualitySwitch = true;
  } else {
    _clusterQualitySwitch = false;
  }

  streamlog_out( DEBUG2 ) << "ClusterQualitySwitch " << _clusterQualitySwitch << endl;

  // Minimum cluster SNR
  if ( count_if(_minTotalSNRVec.begin(), _minTotalSNRVec.end(), bind2nd(greater<float >(), 0) ) != 0 ) {
    _minTotalSNRSwitch = true;
  } else {
    _minTotalSNRSwitch = false;
  }

  streamlog_out( DEBUG2 ) << "MinTotalSNRSwitch " << _minTotalSNRSwitch << endl;

  // Minimum cluster number
  if ( count_if(_minClusterNoVec.begin(), _minClusterNoVec.end(), bind2nd(greater<int>(), 0) ) != 0 ) {
    _minClusterNoSwitch = true;
  } else {
    _minClusterNoSwitch = false;
  }

  streamlog_out ( DEBUG2 ) << "MinClusterNoSwitch " << _minClusterNoSwitch << endl;

  // Maximum cluster number
  if ( count_if(_maxClusterNoVec.begin(), _maxClusterNoVec.end(), bind2nd(greater<int>(), 0) ) != 0 ) {
    _maxClusterNoSwitch = true;
    vector<int >::iterator start = _maxClusterNoVec.begin();
    while ( true ) {
      start = find_if ( start, _maxClusterNoVec.end(), bind2nd(less<int>(), 0));
      if ( start == _maxClusterNoVec.end() ) {
        break;
      } else {
        (*start) = numeric_limits<int>::max();
      }
    }
  } else {
    _maxClusterNoSwitch = false;
  }

  // number of hit pixel inside a cluster for DFF cluser
  if ( count_if(_DFFNHitsCuts.begin(), _DFFNHitsCuts.end(), bind2nd(greater<int>(), 0) ) != 0 ) {
    _dffnhitsswitch = true;
  } else {
    _dffnhitsswitch = false;
  }


  streamlog_out( DEBUG2 )  << "MaxClusterNoSwitch " << _maxClusterNoSwitch << endl;


  // inside ROI
  vector<float>::iterator roiIter = _tempInsideROI.begin();
  while ( roiIter != _tempInsideROI.end() ) {
    int detectorID = static_cast<int> ( *roiIter++ );
    float xBL = ( *roiIter++ );
    float yBL = ( *roiIter++ );
    float xTR = ( *roiIter++ );
    float yTR = ( *roiIter++ );
    if ( detectorID >= 0 ) {
      EUTelROI roi( detectorID, xBL, yBL, xTR, yTR);
      _insideROIVec.push_back( roi );
      streamlog_out( DEBUG1 )  << roi << endl;
    }
  }

  _insideROISwitch = (  !_insideROIVec.empty() );
  streamlog_out( DEBUG2 )  << "InsideROISwitch " << _insideROISwitch << endl;

  // outside ROI
  roiIter = _tempOutsideROI.begin();
  while ( roiIter != _tempOutsideROI.end() ) {
    int detectorID = static_cast<int> ( *roiIter++ );
    float xBL = ( *roiIter++ );
    float yBL = ( *roiIter++ );
    float xTR = ( *roiIter++ );
    float yTR = ( *roiIter++ );
    if ( detectorID >= 0 ) {
      EUTelROI roi( detectorID, xBL, yBL, xTR, yTR);
      _outsideROIVec.push_back( roi );
      streamlog_out( DEBUG1 )  << roi << endl;
    }
  }

  _outsideROISwitch = (  !_outsideROIVec.empty() );
  streamlog_out( DEBUG2 ) << "OutsideROISwitch " << _outsideROISwitch << endl;

  // max cluster noise
  if ( count_if(_maxClusterNoiseVec.begin(), _maxClusterNoiseVec.end(), bind2nd(greater_equal<float >(), 0)) != 0 ) {
    _maxClusterNoiseSwitch = true;
  } else {
    _maxClusterNoiseSwitch = false;
  }
  streamlog_out( DEBUG2 )  << "MaxClusterNoiseSwitch " << _maxClusterNoiseSwitch << endl;


  streamlog_out( DEBUG2 ) << "SameNumberofHitSwitch " << _sameNumberOfHitSwitch  << endl;

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  _isFirstEvent = true;


}


void EUTelClusterFilter::initializeGeometry(LCEvent * event) {

  // try the noise collection
  try {

    // get the noise collection
    LCCollectionVec * noiseCollection = dynamic_cast< LCCollectionVec * > ( event->getCollection( _noiseCollectionName ) ) ;

    // prepare a CellIDDecoder
    CellIDDecoder< TrackerDataImpl > noiseDecoder( noiseCollection ) ;

    // this is the size
    _noOfDetectors  = noiseCollection->size();

    // clear the ancillaryIndexMap
    _ancillaryIndexMap.clear();

    for ( size_t iDetector = 0; iDetector < noiseCollection->size(); ++iDetector ) {
      TrackerDataImpl * noise = dynamic_cast< TrackerDataImpl * > ( noiseCollection->getElementAt( iDetector ) );
      int sensorID =  noiseDecoder( noise )["sensorID"];
      _ancillaryIndexMap.insert( make_pair( sensorID, iDetector ) );
    }
  } catch ( lcio::DataNotAvailableException ) {
    // just catch it!
    _noOfDetectors  = 0;
  }
}

void EUTelClusterFilter::checkCriteria() {

  if ( _noOfDetectors == 0 ) {
    // it means that the information is not available and the criteria
    // can't be verified. Returning immediately
    return ;
  }

  // reset the cluster counter
  _totalClusterCounter.clear();
  _acceptedClusterCounter.clear();
  for ( size_t iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) {
    _totalClusterCounter.push_back(0);
    _acceptedClusterCounter.push_back(0);
  }

  _rejectionMap.clear();


  // check the consistency of selection thresholds
  if ( _minTotalChargeSwitch ) {
    if (  _minTotalChargeVec.size() != static_cast< unsigned >( _noOfDetectors ) ) {
      streamlog_out ( ERROR1 ) << "The threshold vector on the total cluster charge did not match the right size \n"
                               <<  "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minTotalChargeVec.size() <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _minTotalChargeSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Total cluster charge criterion verified and switched on" << endl;
      vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinTotalChargeCut", rejectedCounter));
    }
  }

  if ( _minTotalSNRSwitch ) {
    if ( _minTotalSNRVec.size() != static_cast< unsigned >( _noOfDetectors ) ) {
      streamlog_out ( ERROR1 ) << "The threshold vector on the total cluster SNR did not match the right size \n"
                               <<  "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minTotalSNRVec.size() <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _minTotalSNRSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Total cluster SNR criterion verified and switched on" << endl;
      vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinTotalSNRCut", rejectedCounter));
    }
  }

  if ( _minNChargeSwitch ) {
    unsigned int module = _noOfDetectors + 1;
    if ( _minNChargeVec.size() % module != 0 ) {
      streamlog_out ( ERROR1 ) << "The threshold vector for the N pixels charge did not match the right size \n "
                               <<  "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minNChargeVec.size() <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _minNChargeSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) <<"N pixel charge criterion verified and switched on" << endl;
      vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinNChargeCut", rejectedCounter));
    }
  }

  if ( _minNSNRSwitch ) {
    unsigned int module = _noOfDetectors + 1;
    if ( _minNSNRVec.size() % module != 0 ) {
      streamlog_out ( ERROR1 ) << "The threshold vector for the N pixels SNR did not match the right size \n "
                               <<  "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minNSNRVec.size()  <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _minNSNRSwitch =  false;
    } else {
      streamlog_out ( DEBUG1 ) <<"N pixel SNR criterion verified and switched on" << endl;
      vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinNSNRCut", rejectedCounter));
    }
  }

  if ( _minNxNChargeSwitch ) {
    unsigned int module = _noOfDetectors + 1;
    if ( _minNxNChargeVec.size() % module != 0 ) {
      streamlog_out ( ERROR1 ) << "The threshold vector for the N x N pixels charge did not match the right size \n "
                               <<  "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minNxNChargeVec.size()  << "\n"
                               << "Disabling the selection criterion and continue without" << endl;
      _minNxNChargeSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) <<"NxN pixel charge criterion verified and switched on" << endl;
      vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinNxNChargeCut", rejectedCounter));
    }
  }

  if ( _minNxNSNRSwitch ) {
    unsigned int module = _noOfDetectors + 1;
    if ( _minNxNSNRVec.size() % module != 0 ) {
      streamlog_out ( ERROR1 ) << "The threshold vector for the N x N pixels SNR did not match the right size \n "
                               <<  "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minNxNSNRVec.size() << "\n"
                               << "Disabling the selection criterion and continue without" << endl;
      _minNxNSNRSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) <<"NxN pixel SNR criterion verified and switched on" << endl;
      vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinNxNSNRCut", rejectedCounter));
    }
  }

  if ( _minSeedChargeSwitch ) {
    if (  _minSeedChargeVec.size() != static_cast< unsigned >( _noOfDetectors ) ) {
      streamlog_out ( ERROR1 ) << "The threshold vector on the seed charge did not match the right size \n"
                               <<  "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minSeedChargeVec.size()
                               <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _minSeedChargeSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Seed charge criterion verified and switched on" << endl;
      vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinSeedChargeCut", rejectedCounter));
    }
  }

  if ( _minSeedSNRSwitch ) {
    if ( _minSeedSNRVec.size() != static_cast< unsigned >( _noOfDetectors )) {
      streamlog_out ( ERROR1 ) << "The threshold vector on the seed SNR did not match the right size \n"
                               <<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minSeedSNRVec.size()
                               <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _minSeedSNRSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Seed SNR criterion verified and switched on" << endl;
      vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinSeedSNRCut", rejectedCounter));
    }
  }

  if ( _clusterQualitySwitch ) {
    if ( _clusterQualityVec.size() != static_cast< unsigned >( _noOfDetectors )) {
      streamlog_out ( ERROR1 ) << "The cluster quality vector did not match the right size \n"
                               << "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _minSeedChargeVec.size()
                               <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _clusterQualitySwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Cluster quality criterion verified and switched on" << endl;
      vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("ClusterQualityCut", rejectedCounter ));
    }
  }

  if ( _minClusterNoSwitch ) {
    if ( _minClusterNoVec.size() != static_cast< unsigned >( _noOfDetectors )) {
      streamlog_out ( ERROR1 ) << "The minimum cluster number vector did not match the right size \n"
                               << "The number of planes is " << _noOfDetectors << " while the thresholds are "
                               << _minClusterNoVec.size() <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _minClusterNoSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Minimum cluster number criterion verified and switched on" << endl;
      vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinClusterNoCut", rejectedCounter ));
    }
  }

  if ( _maxClusterNoSwitch ) {
    if ( _maxClusterNoVec.size() != static_cast< unsigned >( _noOfDetectors )) {
      streamlog_out ( ERROR1 ) << "The maximum cluster number vector did not match the right size \n"
                               << "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _maxClusterNoVec.size() <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _maxClusterNoSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Maximum cluster number criterion verified and switched on" << endl;
      vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MaxClusterNoCut", rejectedCounter ));
    }
  }

  if ( _maxClusterNoiseSwitch ) {
    if ( _maxClusterNoiseVec.size() != static_cast< unsigned >( _noOfDetectors )) {
      streamlog_out ( ERROR1 ) << "The maximum cluster noise vector did not match the right size \n"
                               << "The number of planes is " << _noOfDetectors
                               << " while the thresholds are " << _maxClusterNoiseVec.size()   << "\n"
                               << "Disabling the selection criterion and continue without" << endl;
      _maxClusterNoSwitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Maximum cluster noise criterion verified and switched on" << endl;
      vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MaxClusterNoiseCut", rejectedCounter ));
    }
  }

  if ( _dffnhitsswitch ) {
    if ( _DFFNHitsCuts.size() != static_cast< unsigned >( _noOfDetectors )) {
      streamlog_out ( ERROR1 ) << "The minimum hit pixel vector did not match the right size \n"
                               << "The number of planes is " << _noOfDetectors << " while the thresholds are "
                               << _DFFNHitsCuts.size() <<  "\n"
                               <<  "Disabling the selection criterion and continue without" << endl;
      _dffnhitsswitch = false;
    } else {
      streamlog_out ( DEBUG1 ) << "Minimum hit pixel number criterion verified and switched on" << endl;
      vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
      _rejectionMap.insert( make_pair("MinHitPixel", rejectedCounter ));
    }
  }



  if ( _insideROISwitch ) {
    streamlog_out ( DEBUG1 ) << "Inside ROI criterion verified and switched on" << endl;
    vector<unsigned int > rejectedCounter( _noOfDetectors, 0 );
    _rejectionMap.insert( make_pair("InsideROICut", rejectedCounter ));
  }

  if ( _outsideROISwitch ) {
    streamlog_out ( DEBUG1 ) << "Outside ROI criterion verified and switched on" << endl;
    vector<unsigned int > rejectedCounter( _noOfDetectors, 0 );
    _rejectionMap.insert( make_pair("OutsideROICut", rejectedCounter ));
  }

  if ( _sameNumberOfHitSwitch ) {
    streamlog_out ( DEBUG1 ) << "Same number of hits criterion verified and switched on" << endl;
    vector<unsigned int > rejectedCounter( _noOfDetectors, 0 );
    _rejectionMap.insert( make_pair("SameNumberOfHitCut", rejectedCounter ));
  }

}


void EUTelClusterFilter::processRunHeader (LCRunHeader * rdr) {

  // increment the run counter
  ++_iRun;

  if ( isFirstEvent() ) {
    unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr) ;
    runHeader->addProcessor(type());
  }

}


void EUTelClusterFilter::processEvent (LCEvent * event) {

    ++_iEvt;

    if ( isFirstEvent() )
    {
        // try to guess the total number of sensors
        initializeGeometry( event );
        checkCriteria();
        _isFirstEvent = false;
    }


    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> ( event );
    if ( evt->getEventType() == kEORE )
    {
        streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
        return ;
    }
    else if ( evt->getEventType() == kUNKNOWN )
    {
        streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }

    try
    {
        LCCollectionVec * pulseCollectionVec    =   dynamic_cast <LCCollectionVec *> (evt->getCollection(_inputPulseCollectionName));
        LCCollectionVec * filteredCollectionVec =   new LCCollectionVec(LCIO::TRACKERPULSE);
        CellIDEncoder<TrackerPulseImpl> outputEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, filteredCollectionVec);
        CellIDDecoder<TrackerPulseImpl> inputDecoder(pulseCollectionVec);

        vector<int > acceptedClusterVec;
        vector<int > clusterNoVec(_noOfDetectors, 0);

        // CLUSTER BASED CUTS
        for ( int iPulse = 0; iPulse < pulseCollectionVec->getNumberOfElements(); iPulse++ )
        {
            streamlog_out ( DEBUG1 ) << "Filtering cluster " << iPulse + 1  << " / " << pulseCollectionVec->getNumberOfElements() << endl;
            TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl* > (pulseCollectionVec->getElementAt(iPulse));
            ClusterType type         = static_cast<ClusterType> (static_cast<int> ( inputDecoder(pulse)["type"] ));
            EUTelVirtualCluster * cluster;
            SparsePixelType       pixelType;

            if ( type == kEUTelDFFClusterImpl )
            {
                cluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData() ) );
            }
            else if ( type == kEUTelBrickedClusterImpl )
            {
                cluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData() ) );

                if ( _noiseRelatedCuts )
                {
                    // the EUTelBrickedClusterImpl doesn't contain the noise and status
                    // information in the TrackerData object. So this is the right
                    // place to attach to the cluster the noise information.
                    //! ---
                    //! ((NOTE TAKI)): actually i think each cluster does contain its own noise values already!
                    //! each candidate, that was created in clusearch, already had its noise set properly!
                    //! is this information kept when storing the clusters inside the corresponding collection?
                    //! not sure what happened to the status, as well!
                    //! well... this routine here seems to set the noise again but will set noise = 0, if a pixel is a bad one.
                    //! this might cause the bricked cluster to see some pixels with noise = 0 again!
                    //! ---

                    try
                    {
                        LCCollectionVec * noiseCollectionVec  = dynamic_cast<LCCollectionVec * > ( evt->getCollection( _noiseCollectionName )) ;
                        LCCollectionVec * statusCollectionVec = dynamic_cast<LCCollectionVec * > ( evt->getCollection( _statusCollectionName )) ;
                        CellIDDecoder<TrackerDataImpl> noiseDecoder(noiseCollectionVec);

                        int detectorID  = cluster->getDetectorID();
                        int detectorPos = _ancillaryIndexMap[ detectorID ];
                        TrackerDataImpl    * noiseMatrix  = dynamic_cast<TrackerDataImpl    *> ( noiseCollectionVec->getElementAt(detectorPos) );
                        TrackerRawDataImpl * statusMatrix = dynamic_cast<TrackerRawDataImpl *> ( statusCollectionVec->getElementAt(detectorPos) );
                        EUTelMatrixDecoder   noiseMatrixDecoder(noiseDecoder, noiseMatrix);

                        int xSeed, ySeed, xClusterSize, yClusterSize;
                        cluster->getCenterCoord(xSeed, ySeed);
                        cluster->getClusterSize(xClusterSize, yClusterSize);
                        vector<float > noiseValues;
                        for ( int yPixel = ySeed - ( yClusterSize / 2 ); yPixel <= ySeed + ( yClusterSize / 2 ); yPixel++ )
                        {
                            for ( int xPixel = xSeed - ( xClusterSize / 2 ); xPixel <= xSeed + ( xClusterSize / 2 ); xPixel++ )
                            {

                                // always check we are still within the sensor!!!
                                if ( ( xPixel >= noiseMatrixDecoder.getMinX() )  &&  ( xPixel <= noiseMatrixDecoder.getMaxX() ) &&
                                     ( yPixel >= noiseMatrixDecoder.getMinY() )  &&  ( yPixel <= noiseMatrixDecoder.getMaxY() ) )
                                {
                                    int index = noiseMatrixDecoder.getIndexFromXY(xPixel, yPixel);

                                    // the corresponding position in the status matrix has to be HITPIXEL
                                    // in the EUTelClusteringProcessor, we verify also that
                                    // the pixel isHit, but this cannot be done in this
                                    // processor, since the status matrix could have been reset
                                    //
                                    // bool isHit  = ( statusMatrix->getADCValues()[index] ==
                                    // EUTELESCOPE::HITPIXEL );
                                    //
                                    if( static_cast< int >( statusMatrix->getADCValues().size() ) > index )
                                    {
                                    bool isBad  = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::BADPIXEL );
                                    if ( !isBad )
                                    {
                                        noiseValues.push_back( noiseMatrix->getChargeValues()[index] );
                                    }
                                    else
                                    {
                                        noiseValues.push_back( 0. );
                                    }
                                    }

                                }
                                else
                                {
                                    noiseValues.push_back( 0. );
                                }
                            }
                        }
                        cluster->setNoiseValues( noiseValues );
                    }
                    catch ( lcio::Exception& e )
                    {
                        streamlog_out ( ERROR1 ) << e.what() << endl << "Continuing w/o noise based cuts" << endl;
                        _noiseRelatedCuts = false;
                    }
                }
            }
            else if ( type == kEUTelFFClusterImpl )
            {
                cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData() ) );

                if ( _noiseRelatedCuts )
                {
                    // the EUTelFFClusterImpl doesn't contain the noise and status
                    // information in the TrackerData object. So this is the right
                    // place to attach to the cluster the noise information.
                    try
                    {
                        LCCollectionVec * noiseCollectionVec  = dynamic_cast<LCCollectionVec * > ( evt->getCollection( _noiseCollectionName )) ;
                        LCCollectionVec * statusCollectionVec = dynamic_cast<LCCollectionVec * > ( evt->getCollection( _statusCollectionName )) ;
                        CellIDDecoder<TrackerDataImpl> noiseDecoder(noiseCollectionVec);

                        int detectorID  = cluster->getDetectorID();
                        int detectorPos = _ancillaryIndexMap[ detectorID ];
                        TrackerDataImpl    * noiseMatrix  = dynamic_cast<TrackerDataImpl    *> ( noiseCollectionVec->getElementAt(detectorPos) );
                        TrackerRawDataImpl * statusMatrix = dynamic_cast<TrackerRawDataImpl *> ( statusCollectionVec->getElementAt(detectorPos) );
                        EUTelMatrixDecoder   noiseMatrixDecoder(noiseDecoder, noiseMatrix);

                        int xSeed, ySeed, xClusterSize, yClusterSize;
                        cluster->getCenterCoord(xSeed, ySeed);
                        cluster->getClusterSize(xClusterSize, yClusterSize);
                        vector<float > noiseValues;
                        for ( int yPixel = ySeed - ( yClusterSize / 2 ); yPixel <= ySeed + ( yClusterSize / 2 ); yPixel++ )
                        {
                            for ( int xPixel = xSeed - ( xClusterSize / 2 ); xPixel <= xSeed + ( xClusterSize / 2 ); xPixel++ )
                            {

                                // always check we are still within the sensor!!!
                                if ( ( xPixel >= noiseMatrixDecoder.getMinX() )  &&  ( xPixel <= noiseMatrixDecoder.getMaxX() ) &&
                                     ( yPixel >= noiseMatrixDecoder.getMinY() )  &&  ( yPixel <= noiseMatrixDecoder.getMaxY() ) )
                                {
                                    int index = noiseMatrixDecoder.getIndexFromXY(xPixel, yPixel);

                                    // the corresponding position in the status matrix has to be HITPIXEL
                                    // in the EUTelClusteringProcessor, we verify also that
                                    // the pixel isHit, but this cannot be done in this
                                    // processor, since the status matrix could have been reset
                                    //
                                    // bool isHit  = ( statusMatrix->getADCValues()[index] ==
                                    // EUTELESCOPE::HITPIXEL );
                                    // 
                                    ///
                                    if( static_cast< int >(statusMatrix->getADCValues().size()) > index )
                                    {

                                    bool isBad  = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::BADPIXEL );
                                    if ( !isBad )
                                    {
                                        noiseValues.push_back( noiseMatrix->getChargeValues()[index] );
                                    }
                                    else
                                    {
                                        noiseValues.push_back( 0. );
                                    }
                                    }

                                }
                                else
                                {
                                    noiseValues.push_back( 0. );
                                }
                            }
                        }
                        cluster->setNoiseValues( noiseValues );
                    }
                    catch ( lcio::Exception& e )
                    {
                        streamlog_out ( ERROR1 ) << e.what() << endl << "Continuing w/o noise based cuts" << endl;
                        _noiseRelatedCuts = false;
                    }
                }
            }
            else if ( type == kEUTelSparseClusterImpl )
            {
                // knowing that is a sparse cluster is not enough we need also
                // to know also the sparse pixel type. This information is
                // available in the "original_zsdata" collection. Let's get it!
                LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
                TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
                CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
                pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));

                if ( pixelType == kEUTelGenericSparsePixel )
                {
                    cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel > ( static_cast<TrackerDataImpl* > (pulse->getTrackerData() ));

                    auto recasted = dynamic_cast<EUTelSparseClusterImpl<EUTelGenericSparsePixel>*>( cluster );

                    if ( _noiseRelatedCuts )
                    {
                        // the EUTelSparseClusterImpl<EUTelGenericSparsePixel>
                        // doesn't contain any intrinsic noise information. So we
                        // need to get them from the input noise collection.
                        try
                        {
                            LCCollectionVec * noiseCollectionVec = dynamic_cast<LCCollectionVec *> ( evt->getCollection( _noiseCollectionName ));
                            CellIDDecoder<TrackerDataImpl > noiseDecoder( noiseCollectionVec ) ;

                            int detectorID = cluster->getDetectorID();
                            int detectorPos = _ancillaryIndexMap[ detectorID ];
                            TrackerDataImpl    * noiseMatrix = dynamic_cast<TrackerDataImpl *> ( noiseCollectionVec->getElementAt( detectorPos ));
                            EUTelMatrixDecoder   noiseMatrixDecoder( noiseDecoder, noiseMatrix ) ;

			    auto& pixelVec = recasted->getPixels();

                            vector<float > noiseValues;
                            for(auto& sparsePixel: pixelVec) {
                                int index = noiseMatrixDecoder.getIndexFromXY( sparsePixel.getXCoord(), sparsePixel.getYCoord() );
                                noiseValues.push_back( noiseMatrix->getChargeValues()[ index ] );
                            }
                            cluster->setNoiseValues( noiseValues ) ;
                        } catch ( lcio::Exception&  e ) {
                            streamlog_out ( ERROR1 )  << e.what() << "\n" << "Continuing without noise based cuts" << endl;
                            _noiseRelatedCuts = false;
                        }
                        streamlog_out ( DEBUG1 ) << "Noise related cuts may be used" << endl;
                    }
                }
                else
                {
                  streamlog_out ( ERROR4 ) << "Unknown pixel type. Sorry for quitting" << endl;
                  throw UnknownDataTypeException("Pixel type unknown");
                }
		}
            else
            {
                streamlog_out ( ERROR4 ) << "Unknown cluster type. Sorry for quitting" << endl;
                throw UnknownDataTypeException("Cluster type unknown");
            }

            // increment the event counter
            _totalClusterCounter[ _ancillaryIndexMap[ cluster->getDetectorID() ] ]++;

            bool isAccepted = true;

            if ( type == kEUTelDFFClusterImpl )
            {
                isAccepted &= isAboveNumberOfHitPixel(cluster);
            }
            else
            {
                isAccepted &= isAboveMinTotalCharge(cluster);
                isAccepted &= isAboveMinTotalSNR(cluster);
                isAccepted &= isAboveNMinCharge(cluster);
                isAccepted &= isAboveNMinSNR(cluster);
                isAccepted &= isAboveNxNMinCharge(cluster);
                isAccepted &= isAboveNxNMinSNR(cluster);
                isAccepted &= isAboveMinSeedCharge(cluster);
                isAccepted &= isAboveMinSeedSNR(cluster);
                isAccepted &= isBelowMaxClusterNoise(cluster);
            }
            isAccepted &= hasQuality(cluster);
            isAccepted &= isInsideROI(cluster);
            isAccepted &= isOutsideROI(cluster);

            if ( isAccepted )  acceptedClusterVec.push_back(iPulse);

            delete cluster;

        }

        vector<int >::iterator cluIter = acceptedClusterVec.begin();
        while ( cluIter != acceptedClusterVec.end() )
        {
            TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl* > (pulseCollectionVec->getElementAt(*cluIter));
            int detectorID  = inputDecoder(pulse)["sensorID"];
            clusterNoVec[ _ancillaryIndexMap[ detectorID ] ]++;
            ++cluIter;
        }

        bool areClusterEnoughTemp     = areClusterEnough(clusterNoVec);
        bool areClusterTooManyTemp    = areClusterTooMany(clusterNoVec);
        bool hasSameNumberOfHitTemp   = hasSameNumberOfHit(clusterNoVec);
        bool isEventAccepted = areClusterEnoughTemp && !areClusterTooManyTemp;
        isEventAccepted &= hasSameNumberOfHitTemp;

        if ( ! isEventAccepted ) acceptedClusterVec.clear();

        if ( acceptedClusterVec.empty() )
        {
            delete filteredCollectionVec;
            if ( ! areClusterEnoughTemp )
            {
                streamlog_out ( DEBUG1 ) << "Not enough clusters passed the selection " << endl;
            }
            else if ( areClusterTooManyTemp )
            {
               streamlog_out ( DEBUG1 ) << "Too many clusters passed the selection " << endl;
            }
            else
            {
                streamlog_out ( DEBUG1 ) << "No clusters passed the selection" << endl;
            }
            if ( _skipEmptyEvent )
            {
                if ( ! areClusterEnoughTemp )
                {
                    streamlog_out (WARNING0 ) << "Skipping event because too few clusters passed the selection" << endl;
                }
                else if ( areClusterTooManyTemp )
                {
                    streamlog_out (WARNING0 ) << "Skipping event because too many clusters passed the selection" << endl;
                }
                else
                {
                    streamlog_out ( WARNING0 ) << "Skipping event because no clusters passed the selection" << endl;
                }
                throw SkipEventException(this);
            }
            else
            {
                return;
            }
        }
        else
        {
            vector<int >::iterator iter = acceptedClusterVec.begin();
            while ( iter != acceptedClusterVec.end() )
            {
                TrackerPulseImpl * pulse     = dynamic_cast<TrackerPulseImpl *> ( pulseCollectionVec->getElementAt( *iter ) );
                TrackerPulseImpl * accepted  = new TrackerPulseImpl;
                accepted->setCellID0( pulse->getCellID0() );
                accepted->setCellID1( pulse->getCellID1() );
                accepted->setTime(    pulse->getTime()    );
                accepted->setCharge(  pulse->getCharge()  );
                accepted->setQuality( pulse->getQuality() );
                accepted->setTrackerData( pulse->getTrackerData() );
                filteredCollectionVec->push_back(accepted);
                _acceptedClusterCounter[ _ancillaryIndexMap[ inputDecoder(pulse)["sensorID"] ] ]++;
                ++iter;
            }
            evt->addCollection(filteredCollectionVec, _outputPulseCollectionName);
        }
    }
    catch (DataNotAvailableException& e )
    {
        streamlog_out ( MESSAGE2 )  << "Input collection not found in the current event." << endl;
        return;
    }
}


bool EUTelClusterFilter::areClusterEnough(std::vector<int > clusterNoVec) const {

  if ( ! _minClusterNoSwitch ) return true;

  bool areEnough = true;

  for ( size_t iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) {

    if ( clusterNoVec[iDetector] < _minClusterNoVec[iDetector] ) {
      streamlog_out ( DEBUG2 )  << "Rejected event because on detector " << iDetector
                                << " there are " << clusterNoVec[iDetector] << "\nwhile "
                                << "at least " << _minClusterNoVec[iDetector] << " were requested" << endl;
      _rejectionMap["MinClusterNoCut"][iDetector]++;
      areEnough = false;
    }
  }

  return areEnough;

}

bool EUTelClusterFilter::areClusterTooMany(std::vector<int > clusterNoVec) const {

  if ( ! _maxClusterNoSwitch ) return false;

  bool areTooMany = false;

  for ( size_t iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) {

    if ( clusterNoVec[iDetector] > _maxClusterNoVec[iDetector] ) {
      streamlog_out ( DEBUG2 )  << "Rejected event because on detector " << iDetector
                                << " there are " << clusterNoVec[iDetector] << "\nwhile "
                                << "only " << _maxClusterNoVec[iDetector] << " were allowed" << endl;
      _rejectionMap["MaxClusterNoCut"][iDetector]++;
      areTooMany = true;
    }
  }

  return areTooMany;
}

bool EUTelClusterFilter::hasSameNumberOfHit(std::vector<int > clusterNoVec) const {

  if ( ! _sameNumberOfHitSwitch ) return true;

  bool hasSameNumber = true;

  vector<int >::iterator iter = clusterNoVec.begin() + 1;
  while ( iter != clusterNoVec.end() ) {
    if ( *( iter - 1) !=  (*iter ) ) {
      streamlog_out ( DEBUG2 )  << "Rejected event because the number of hit is different from one sensor to the other"  << endl;
      hasSameNumber = false;
      break;
    }
    ++iter;
  }

  if ( ! hasSameNumber ) {
    vector<int >::iterator iter = clusterNoVec.begin();
    vector<unsigned int >::iterator iter2 = _rejectionMap["SameNumberOfHitCut"].begin();
    while ( ( iter  != clusterNoVec.end() ) ||
            ( iter2 != _rejectionMap["SameNumberOfHitCut"].end() ) ) {
      (*iter2) += (*iter);
      ++iter;
      ++iter2;
    }
  }

  return hasSameNumber;
}

bool EUTelClusterFilter::isAboveNumberOfHitPixel(EUTelVirtualCluster * cluster) const {
  if ( !_dffnhitsswitch ) {
    return true;
  }
  streamlog_out ( DEBUG1 ) << "Filtering against number of hit pixel inside a cluster " << endl;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap[ detectorID ];

  if ( static_cast< int >(cluster->getTotalCharge()) >= _DFFNHitsCuts[detectorPos] ) return true;
  else {
    streamlog_out ( DEBUG2 )  << "Rejected cluster because the number of hit pixel is " << static_cast< int >(cluster->getTotalCharge())
                              << " and the threshold is " << _DFFNHitsCuts[detectorPos] << endl;
    _rejectionMap["MinHitPixel"][detectorPos]++;
    return false;
  }
}



bool EUTelClusterFilter::isAboveMinTotalCharge(EUTelVirtualCluster * cluster) const {

  if ( !_minTotalChargeSwitch ) {
    return true;
  }
  streamlog_out ( DEBUG1 ) << "Filtering against the total charge " << endl;

  int detectorID  = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap[ detectorID ];

  if ( cluster->getTotalCharge() > _minTotalChargeVec[detectorPos] ) return true;
  else {
    streamlog_out ( DEBUG2 )  << "Rejected cluster because its charge is " << cluster->getTotalCharge()
                              << " and the threshold is " << _minTotalChargeVec[detectorPos] << endl;
    _rejectionMap["MinTotalChargeCut"][detectorPos]++;
    return false;
  }
}

bool EUTelClusterFilter::isAboveMinTotalSNR(EUTelVirtualCluster * cluster) const {

  if ( !_noiseRelatedCuts   ) return true;
  if ( !_minTotalSNRSwitch  ) return true;

  int detectorID  = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];

  streamlog_out ( DEBUG1 ) << "Filtering against the minimum total SNR " << endl;
  if  ( cluster->getClusterSNR() > _minTotalSNRVec[ detectorPos ] ) return true;
  else {
    streamlog_out ( DEBUG2 )  << "Rejected cluster because its SNR is " << cluster->getClusterSNR()
                              << " and the threshold is " << _minTotalSNRVec[ detectorPos ] << endl;
    _rejectionMap["MinTotalSNRCut"][detectorPos]++;
    return false;
  }
}

bool EUTelClusterFilter::isAboveNMinCharge(EUTelVirtualCluster * cluster) const {

  if ( !_minNChargeSwitch ) return true;

  streamlog_out ( DEBUG1 ) << "Filtering against the N Pixel charge " << endl;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  vector<float >::const_iterator iter = _minNChargeVec.begin();
  while ( iter != _minNChargeVec.end() ) {
    int nPixel      = static_cast<int > (*iter);
    float charge    = cluster->getClusterCharge(nPixel);
    float threshold = (* (iter + detectorPos + 1) );
    if ( charge > threshold ) {
      iter += _noOfDetectors + 1;
    } else {
      streamlog_out ( DEBUG2 ) << "Rejected cluster because its charge over " << (*iter) << " is " << charge
                               << " and the threshold is " << threshold << endl;
      _rejectionMap["MinNChargeCut"][detectorPos]++;
      return false;
    }
  }
  return true;
}


bool EUTelClusterFilter::isAboveNMinSNR(EUTelVirtualCluster * cluster) const {

  if ( !_noiseRelatedCuts ) return true;
  if ( !_minNSNRSwitch    ) return true;

  streamlog_out ( DEBUG1 ) << "Filtering against the N pixel SNR " << endl;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  vector<float >::const_iterator iter = _minNSNRVec.begin();
  while ( iter !=  _minNSNRVec.end() ) {
    int nPixel      = static_cast<int > (*iter);
    float SNR       = cluster->getClusterSNR(nPixel);
    float threshold = (* (iter + detectorPos + 1 ) );
    if ( SNR > threshold ) {
      iter += _noOfDetectors + 1;
    } else {
      streamlog_out ( DEBUG2 )  << "Rejected cluster because its SNR over " << (*iter) << " is " << SNR
                                << " and the threshold is " << threshold  << endl;
      _rejectionMap["MinNSNRCut"][detectorPos]++;
      return false;
    }
  }
  return true;
}




bool EUTelClusterFilter::isAboveNxNMinCharge(EUTelVirtualCluster * cluster) const {

  if ( !_minNxNChargeSwitch ) return true;

  streamlog_out ( DEBUG1 ) << "Filtering against the N x N pixel charge" << endl;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  vector<float >::const_iterator iter = _minNxNChargeVec.begin();
  while ( iter != _minNxNChargeVec.end() ) {
    int nxnPixel    = static_cast<int > ( *iter ) ;
    float charge    = cluster->getClusterCharge(nxnPixel, nxnPixel);
    float threshold = (* ( iter + detectorPos + 1 )) ;
    if ( ( threshold <= 0) || (charge > threshold) ) {
      iter += _noOfDetectors + 1;
    } else {
      streamlog_out ( DEBUG2 ) << "Rejected cluster because its charge within a " << (*iter) << " x " << (*iter)
                               << " subcluster is " << charge << " and the threshold is " << threshold << endl;
      _rejectionMap["MinNxNChargeCut"][detectorPos]++;
      return false;
    }
  }
  return true;
}


bool EUTelClusterFilter::isAboveNxNMinSNR(EUTelVirtualCluster * cluster) const {

  if ( !_noiseRelatedCuts  ) return true;
  if ( !_minNxNSNRSwitch   ) return true;

  streamlog_out ( DEBUG1 ) << "Filtering against the N x N pixel charge" << endl;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  vector<float >::const_iterator iter = _minNxNSNRVec.begin();
  while ( iter != _minNxNSNRVec.end() ) {
    int nxnPixel    = static_cast<int > ( *iter ) ;
    float snr       = cluster->getClusterSNR(nxnPixel, nxnPixel);
    float threshold = (* ( iter + detectorPos + 1 )) ;
    if ( ( threshold <= 0) || (snr > threshold) ) {
      iter += _noOfDetectors + 1;
    } else {
      streamlog_out ( DEBUG2 )  << "Rejected cluster because its SNR within a " << (*iter) << " x " << (*iter)
                                << " subcluster is " << snr << " and the threshold is " << threshold << endl;
      _rejectionMap["MinNxNSNRCut"][detectorPos]++;
      return false;
    }
  }
  return true;

}

bool EUTelClusterFilter::isAboveMinSeedCharge(EUTelVirtualCluster * cluster) const {

  if ( !_minSeedChargeSwitch ) return true;

  streamlog_out ( DEBUG1 ) << "Filtering against the seed charge " << endl;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  if ( cluster->getSeedCharge() > _minSeedChargeVec[detectorPos] ) return true;
  else {
    streamlog_out ( DEBUG2 )  << "Rejected cluster because its seed charge is " << cluster->getSeedCharge()
                              << " and the threshold is " <<  _minSeedChargeVec[detectorPos] << endl;
    _rejectionMap["MinSeedChargeCut"][detectorPos]++;
    return false;
  }
}

bool EUTelClusterFilter::isAboveMinSeedSNR(EUTelVirtualCluster * cluster) const {

  if ( !_noiseRelatedCuts  ) return true;
  if ( !_minSeedSNRSwitch  ) return true;

  streamlog_out ( DEBUG1 ) << "Filtering against the seed SNR " << endl;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  if ( cluster->getSeedSNR() > _minSeedSNRVec[detectorPos] ) return true;
  else {
    streamlog_out ( DEBUG2 ) << "Rejected cluster because its seed charge is " << cluster->getSeedSNR()
                             << " and the threshold is " <<  _minSeedSNRVec[detectorPos] << endl;
    _rejectionMap["MinSeedSNRCut"][detectorPos]++;
    return false;
  }
}



bool EUTelClusterFilter::hasQuality(EUTelVirtualCluster * cluster) const {

  if ( !_clusterQualitySwitch ) return true;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  if ( _clusterQualityVec[detectorID] < 0 ) return true;

  ClusterQuality actual = cluster->getClusterQuality();
  ClusterQuality needed = static_cast<ClusterQuality> ( _clusterQualityVec[detectorPos] );

  if ( actual == needed ) return true;
  else {
    streamlog_out ( DEBUG2 ) <<  "Rejected cluster because its quality " << static_cast<int> (actual)
                             << " is not " << _clusterQualityVec[detectorPos] << endl;
    _rejectionMap["ClusterQualityCut"][detectorPos]++;
    return false;
  }
}

bool EUTelClusterFilter::isBelowMaxClusterNoise(EUTelVirtualCluster * cluster) const {

  if ( !_noiseRelatedCuts       ) return true;
  if ( !_maxClusterNoiseSwitch  ) return true;

  streamlog_out ( DEBUG1 ) << "Filtering against the maximum cluster noise"  << endl;
  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  if (  ( cluster->getClusterNoise() < _maxClusterNoiseVec[detectorPos] ) ||
        ( _maxClusterNoiseVec[detectorID] < 0 ) ) return true;
  else {
    streamlog_out ( DEBUG2 )  << "Rejected cluster because its noise is " << cluster->getClusterNoise()
                              << " and the threshold is " <<  _maxClusterNoiseVec[detectorPos] << endl;
    _rejectionMap["MaxClusterNoiseCut"][detectorPos]++;
    return false;
  }
}


bool EUTelClusterFilter::isInsideROI(EUTelVirtualCluster * cluster) const {

  if ( !_insideROISwitch ) return true;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  float x, y;
  cluster->getCenterOfGravity(x, y);

  bool tempAccepted = true;
  vector<EUTelROI>::const_iterator iter = _insideROIVec.begin();
  vector<EUTelROI>::const_iterator end  = _insideROIVec.end();
  while ( true ) {

    // the HasSameID requires the sensorID, so don't replace it with
    // detectorPos
    iter = find_if( iter, end, HasSameID(detectorID));
    if ( iter != end ) {
      if ( (*iter).isInside(x,y) ) {
        tempAccepted &= true;
      }  else {
        _rejectionMap["InsideROICut"][detectorPos]++;
        tempAccepted &= false;
      }
      ++iter;
    } else break;
  }
  return tempAccepted;

}

bool EUTelClusterFilter::isOutsideROI(EUTelVirtualCluster * cluster) const {

  if ( !_outsideROISwitch ) return true;

  int detectorID = cluster->getDetectorID();
  int detectorPos = _ancillaryIndexMap [ detectorID ];
  float x, y;
  cluster->getCenterOfGravity(x, y);

  bool tempAccepted = true;
  vector<EUTelROI>::const_iterator iter = _outsideROIVec.begin();
  vector<EUTelROI>::const_iterator end  = _outsideROIVec.end();
  while ( true ) {
    iter = find_if( iter, end, HasSameID(detectorID));
    if ( iter != end ) {
      if ( !(*iter).isInside(x,y) ) {
        tempAccepted &= true;
      }  else {
        _rejectionMap["OutsideROICut"][detectorPos]++;
        tempAccepted &= false;
      }
      ++iter;
    } else break;
  }
  return tempAccepted;

}

void EUTelClusterFilter::check (LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelClusterFilter::end() {
  streamlog_out ( MESSAGE4 )  << printSummary() << endl;

}

string EUTelClusterFilter::printSummary() const {

  stringstream  ss;

  int bigSpacer    = 20;
  int smallSpacer  = 12;

  int fullLineWidth = bigSpacer + _noOfDetectors * smallSpacer + 1;
  stringstream doubleLine;
  stringstream singleLine;
  for ( int i = 0; i < fullLineWidth ; i++ ) {
    doubleLine << "=";
    singleLine << "-";
  }


  ss << doubleLine.str() << endl
     << " Rejection summary " << endl
     << doubleLine.str() << endl;

  map<string, vector<unsigned int> >::iterator iter = _rejectionMap.begin();
  while ( iter != _rejectionMap.end() ) {
    ss << " " << setiosflags(ios::left) << setw(bigSpacer) << (*iter).first << resetiosflags(ios::left) ;
    vector<unsigned int>::iterator iter2 = (*iter).second.begin();
    while ( iter2 != (*iter).second.end() ) {
      ss << setw(smallSpacer) << ( *iter2) ;
      ++iter2;
    }
    ss << "\n" << singleLine.str() <<  endl;
    ++iter;
  }
  ss << singleLine.str() << "\n" << " " << setiosflags(ios::left)
     << setw(bigSpacer)  << "Accepted clusters " << resetiosflags(ios::left);
  vector<unsigned int>::const_iterator iter3 = _acceptedClusterCounter.begin();
  while ( iter3 != _acceptedClusterCounter.end() ) {
    ss <<  setw(smallSpacer) << ( *iter3 );
    ++iter3;
  }

  ss << "\n" << singleLine.str() << "\n" << " "
     << setw(bigSpacer) << "Total no of clusters" << resetiosflags(ios::left);
  vector<unsigned int>::const_iterator iter4 = _totalClusterCounter.begin();
  while ( iter4 != _totalClusterCounter.end() ) {
    ss <<  setw(smallSpacer) << ( *iter4 );
    ++iter4;
  }
  ss << "\n" << doubleLine.str() << endl;

  return ss.str();

}
