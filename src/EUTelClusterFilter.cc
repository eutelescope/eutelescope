// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelClusterFilter.cc,v 1.9 2007-07-31 14:45:50 bulgheroni Exp $
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

  registerOptionalParameter("NoiseCollectionName","This is the name of the noise collection.\n"
			    "The presence of this collection in the event is allowing all the noise based selection cuts",
			    _noiseCollectionName, string( "noiseDB" ) );
  
  registerOptionalParameter("StatusCollectionName","This is the name of the status collection.\n"
			    "The presence of this collection in the event is allowing all the noise based selection cuts",
			    _statusCollectionName, string( "statusDB" ) );


  registerOptionalParameter("SkipEmptyEvent","If true, a SkipEventException is thrown if after selection\n"
			    "there are no cluster left.",
			    _skipEmptyEvent, static_cast<bool > ( false ) );

  // set the global noise switch to on
  _noiseRelatedCuts = true;

 
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
  
  message<DEBUG> ( log () << "MinTotalChargeSwitch " << _minTotalChargeSwitch );
  
  // N cluster charge
  if ( _minNChargeVec.size() != 0 ) {
    if ( _minNChargeVec[0] != 0 ) {
      _minNChargeSwitch = true;
    } else {
      _minNChargeSwitch = false;
    }
  } else {
    _minNChargeSwitch = false;
  }

  message<DEBUG> ( log () << "MinNChargeSwitch " << _minNChargeSwitch );

  // N cluster SNR
  if ( (_minNSNRVec.size() != 0) && (_minNSNRVec[0] != 0) ) {
    _minNSNRSwitch = true;
  } else {
    _minNSNRSwitch = false;
  }

  message<DEBUG> ( log () << "MinNSNRSwitch " << _minNSNRSwitch );

  // NxN cluster charge
  if (  ( count_if( _minNxNChargeVec.begin(), _minNxNChargeVec.end(), bind2nd(greater<float>(),0)) != 0  ) &&
	( _minNxNChargeVec[0] != 0 ) ) {
    _minNxNChargeSwitch = true;
  } else {
    _minNxNChargeSwitch = false;
  }

  message<DEBUG> ( log () << "MinNxNChargeSwitch " << _minNxNChargeSwitch );	 

  // NxN cluster SNR
  if (  ( count_if( _minNxNSNRVec.begin(), _minNxNSNRVec.end(), bind2nd(greater<float>(),0)) != 0  ) &&
	( _minNxNSNRVec[0] != 0 ) ) {
    _minNxNSNRSwitch = true;
  } else {
    _minNxNSNRSwitch = false;
  }

  message<DEBUG> ( log () << "MinNxNSNRSwitch " << _minNxNSNRSwitch );


  // Seed charge
  if ( count_if( _minSeedChargeVec.begin(), _minSeedChargeVec.end(), bind2nd(greater<float>(), 0) ) != 0 ) {
    _minSeedChargeSwitch = true;
  } else {
    _minSeedChargeSwitch = false;
  }

  message<DEBUG> ( log () << "MinSeedChargeSwitch " << _minSeedChargeSwitch );

  // Seed SNR
  if ( count_if( _minSeedSNRVec.begin(), _minSeedSNRVec.end(), bind2nd(greater<float >(), 0) ) != 0 ) {
    _minSeedSNRSwitch = true;
  } else {
    _minSeedSNRSwitch = false;
  }

  message<DEBUG> ( log () << "MinSeedSNRSwitch " << _minSeedChargeSwitch );

  // Cluster quality
  if ( count_if(_clusterQualityVec.begin(), _clusterQualityVec.end(), bind2nd(greater_equal<int>(), 0) ) != 0 ) {
    _clusterQualitySwitch = true;
  } else {
    _clusterQualitySwitch = false;
  }

  message<DEBUG> ( log () << "ClusterQualitySwitch " << _clusterQualitySwitch );

  // Minimum cluster SNR
  if ( count_if(_minTotalSNRVec.begin(), _minTotalSNRVec.end(), bind2nd(greater<float >(), 0) ) != 0 ) {
    _minTotalSNRSwitch = true;
  } else {
    _minTotalSNRSwitch = false;
  }

  message<DEBUG> ( log () << "MinTotalSNRSwitch " << _minTotalSNRSwitch );

  // Minimum cluster number
  if ( count_if(_minClusterNoVec.begin(), _minClusterNoVec.end(), bind2nd(greater<int>(), 0) ) != 0 ) {
    _minClusterNoSwitch = true;
  } else {
    _minClusterNoSwitch = false;
  }
  
  message<DEBUG> ( log() << "MinClusterNoSwitch " << _minClusterNoSwitch ) ;

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
  
  message<DEBUG> ( log() << "MaxClusterNoSwitch " << _maxClusterNoSwitch ) ;


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
      message<DEBUG> ( log() << roi );
    }
  }
  
  _insideROISwitch = (  _insideROIVec.size() > 0 );
  message<DEBUG> ( log() << "InsideROISwitch " << _insideROISwitch );

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
      message<DEBUG> ( log() << roi );
    }
  }
  
  _outsideROISwitch = (  _outsideROIVec.size() > 0 );
  message<DEBUG> ( log() << "OutsideROISwitch " << _outsideROISwitch );
 
  // max cluster noise
  if ( count_if(_maxClusterNoiseVec.begin(), _maxClusterNoiseVec.end(), bind2nd(greater_equal<float >(), 0)) != 0 ) {
    _maxClusterNoiseSwitch = true;
  } else {
    _maxClusterNoiseSwitch = false;
  }
  message<DEBUG> ( log() << "MaxClusterNoiseSwitch " << _maxClusterNoiseSwitch );


  message<DEBUG> ( log() << "SameNumberofHitSwitch " << _sameNumberOfHitSwitch );

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  _isFirstEvent = true;


}

void EUTelClusterFilter::processRunHeader (LCRunHeader * rdr) {

  // increment the run counter
  ++_iRun;

  if ( isFirstEvent() ) {
    EUTelRunHeaderImpl * runHeader = static_cast<EUTelRunHeaderImpl *> (rdr);
    
    _noOfDetectors = runHeader->getNoOfDetector();
    
    // reset the cluster counter
    _totalClusterCounter.clear();
    _acceptedClusterCounter.clear();
    for ( int iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) {
      _totalClusterCounter.push_back(0);
      _acceptedClusterCounter.push_back(0);
    }

    // check the consistency of selection thresholds
    _rejectionMap.clear();
    
    if ( _minTotalChargeSwitch ) {
      if (  _minTotalChargeVec.size() != (unsigned) _noOfDetectors ) {
	message<ERROR> (log() << "The threshold vector on the total cluster charge did not match the right size \n"
			<<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minTotalChargeVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minTotalChargeSwitch = false;
      } else {
	message<DEBUG> ( "Total cluster charge criterion verified and switched on" );
	vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinTotalChargeCut", rejectedCounter));
      }
    }

    if ( _minTotalSNRSwitch ) {
      if ( _minTotalSNRVec.size() != ( unsigned ) _noOfDetectors ) {
	message<ERROR> (log() << "The threshold vector on the total cluster SNR did not match the right size \n"
			<<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minTotalSNRVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minTotalSNRSwitch = false;
      } else {
	message<DEBUG> ( "Total cluster SNR criterion verified and switched on" );
	vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinTotalSNRCut", rejectedCounter));
      }
    }	
    
    if ( _minNChargeSwitch ) {
      unsigned int module = _noOfDetectors + 1;
      if ( _minNChargeVec.size() % module != 0 ) {
	message<ERROR> (log() << "The threshold vector for the N pixels charge did not match the right size \n "
			<<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minNChargeVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minNChargeSwitch = false;
      } else {
	message<DEBUG> ("N pixel charge criterion verified and switched on");
	vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinNChargeCut", rejectedCounter));
      }
    }
    
    if ( _minNSNRSwitch ) {
      unsigned int module = _noOfDetectors + 1;
      if ( _minNSNRVec.size() % module != 0 ) {
	message<ERROR> (log() << "The threshold vector for the N pixels SNR did not match the right size \n "
			<<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minNSNRVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minNSNRSwitch =  false;
      } else {
	message<DEBUG> ("N pixel SNR criterion verified and switched on");
	vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinNSNRCut", rejectedCounter));
      }
    }	

    if ( _minNxNChargeSwitch ) {
      unsigned int module = _noOfDetectors + 1;
      if ( _minNxNChargeVec.size() % module != 0 ) {
	message<ERROR> ( log() << "The threshold vector for the N x N pixels charge did not match the right size \n "
			 <<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minNxNChargeVec.size() 
			 << "\n"
			 << "Disabling the selection criterion and continue without");
	_minNxNChargeSwitch = false; 
      } else {
	message<DEBUG> ("NxN pixel charge criterion verified and switched on");
	vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinNxNChargeCut", rejectedCounter));
      }
    }

    if ( _minNxNSNRSwitch ) {
      unsigned int module = _noOfDetectors + 1;
      if ( _minNxNSNRVec.size() % module != 0 ) {
	message<ERROR> ( log() << "The threshold vector for the N x N pixels SNR did not match the right size \n "
			 <<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minNxNSNRVec.size() 
			 << "\n"
			 << "Disabling the selection criterion and continue without");
	_minNxNSNRSwitch = false; 
      } else {
	message<DEBUG> ("NxN pixel SNR criterion verified and switched on");
	vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinNxNSNRCut", rejectedCounter));
      }
    }
    
    if ( _minSeedChargeSwitch ) {
      if (  _minSeedChargeVec.size() != (unsigned) _noOfDetectors ) {
	message<ERROR> (log() << "The threshold vector on the seed charge did not match the right size \n"
			<<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minSeedChargeVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minSeedChargeSwitch = false;
      } else {
	message<DEBUG> ( "Seed charge criterion verified and switched on" );
	vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinSeedChargeCut", rejectedCounter));
      }
    }

    if ( _minSeedSNRSwitch ) {
      if ( _minSeedSNRVec.size() != (unsigned) _noOfDetectors ) {
	message<ERROR> (log() << "The threshold vector on the seed SNR did not match the right size \n"
			<<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minSeedSNRVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minSeedSNRSwitch = false;	
      } else {
	message<DEBUG> ( "Seed SNR criterion verified and switched on" );
	vector<unsigned int >  rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinSeedSNRCut", rejectedCounter));
      }
    }

    if ( _clusterQualitySwitch ) {
      if ( _clusterQualityVec.size() != (unsigned) _noOfDetectors ) {
	message<ERROR> (log() << "The cluster quality vector did not match the right size \n"
			<< "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minSeedChargeVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_clusterQualitySwitch = false;
      } else {
	message<DEBUG> ( "Cluster quality criterion verified and switched on");
	vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("ClusterQualityCut", rejectedCounter ));
      }
    }

    if ( _minClusterNoSwitch ) {
      if ( _minClusterNoVec.size() != (unsigned) _noOfDetectors ) {
	message<ERROR> (log() << "The minimum cluster number vector did not match the right size \n"
			<< "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minClusterNoVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minClusterNoSwitch = false;
      } else {
	message<DEBUG> ( "Minimum cluster number criterion verified and switched on");
	vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinClusterNoCut", rejectedCounter ));
      }
    }	

    if ( _maxClusterNoSwitch ) {
      if ( _maxClusterNoVec.size() != (unsigned) _noOfDetectors ) {
	message<ERROR> (log() << "The maximum cluster number vector did not match the right size \n"
			<< "The number of planes is " << _noOfDetectors << " while the thresholds are " << _maxClusterNoVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_maxClusterNoSwitch = false;
      } else {
	message<DEBUG> ( "Maximum cluster number criterion verified and switched on");
	vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MaxClusterNoCut", rejectedCounter ));
      }
    }

    if ( _maxClusterNoiseSwitch ) {
      if ( _maxClusterNoiseVec.size() != (unsigned) _noOfDetectors ) {
	message<ERROR> ( log() << "The maximum cluster noise vector did not match the right size \n"
			 << "The number of planes is " << _noOfDetectors << " while the thresholds are " << _maxClusterNoiseVec.size()
			 << "\n"
			 << "Disabling the selection criterion and continue without");
	_maxClusterNoSwitch = false;
      } else {
	message<DEBUG> ( "Maximum cluster noise criterion verified and switched on");
	vector<unsigned int > rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MaxClusterNoiseCut", rejectedCounter ));
      }
    }
    
    if ( _insideROISwitch ) {
      message<DEBUG> ( "Inside ROI criterion verified and switched on");
      vector<unsigned int > rejectedCounter( _noOfDetectors, 0 );
      _rejectionMap.insert( make_pair("InsideROICut", rejectedCounter ));
    }

    if ( _outsideROISwitch ) {
      message<DEBUG> ( "Outside ROI criterion verified and switched on");
      vector<unsigned int > rejectedCounter( _noOfDetectors, 0 );
      _rejectionMap.insert( make_pair("OutsideROICut", rejectedCounter ));
    }

    if ( _sameNumberOfHitSwitch ) {
      message<DEBUG> ( "Same number of hits criterion verified and switched on");
      vector<unsigned int > rejectedCounter( _noOfDetectors, 0 );
      _rejectionMap.insert( make_pair("SameNumberOfHitCut", rejectedCounter ));
    }

  }
}


void EUTelClusterFilter::processEvent (LCEvent * event) {

  ++_iEvt;

  if ( isFirstEvent() ) _isFirstEvent = false;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> ( event );
  if ( evt->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do.");
    return ;
  }
  
  if (_iEvt % 10 == 0) 
    message<MESSAGE> ( log() << "Filtering clusters on event " << _iEvt ) ;

  try {

    LCCollectionVec * pulseCollectionVec    =   dynamic_cast <LCCollectionVec *> (evt->getCollection(_inputPulseCollectionName));
    LCCollectionVec * filteredCollectionVec =   new LCCollectionVec(LCIO::TRACKERPULSE);
    CellIDEncoder<TrackerPulseImpl> outputEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, filteredCollectionVec);
    CellIDDecoder<TrackerPulseImpl> inputDecoder(pulseCollectionVec);
    
    vector<int > acceptedClusterVec;
    vector<int > clusterNoVec(_noOfDetectors, 0);
    
    // CLUSTER BASED CUTS
    for ( int iPulse = 0; iPulse < pulseCollectionVec->getNumberOfElements(); iPulse++ ) {
      message<DEBUG> ( log() << "Filtering cluster " << iPulse + 1  << " / " << pulseCollectionVec->getNumberOfElements() ) ;
      TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl* > (pulseCollectionVec->getElementAt(iPulse));
      ClusterType type = static_cast<ClusterType> (static_cast<int> ( inputDecoder(pulse)["type"] ));
      EUTelVirtualCluster * cluster;
      
      if ( type == kEUTelFFClusterImpl )  {
	cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData() ) );
	
	if ( _noiseRelatedCuts ) {
	  
	  // the EUTelFFClusterImpl doesn't contain the noise and status
	  // information in the TrackerData object. So this is the right
	  // place to attach to the cluster the noise information.
	  try {
	    LCCollectionVec * noiseCollectionVec  = dynamic_cast<LCCollectionVec * > ( evt->getCollection( _noiseCollectionName )) ;
	    LCCollectionVec * statusCollectionVec = dynamic_cast<LCCollectionVec * > ( evt->getCollection( _statusCollectionName )) ;
	    CellIDDecoder<TrackerDataImpl> noiseDecoder(noiseCollectionVec);
	    
	    int detectorID = cluster->getDetectorID();
	    TrackerDataImpl    * noiseMatrix  = dynamic_cast<TrackerDataImpl    *> ( noiseCollectionVec->getElementAt(detectorID) );
	    TrackerRawDataImpl * statusMatrix = dynamic_cast<TrackerRawDataImpl *> ( statusCollectionVec->getElementAt(detectorID) );
	    EUTelMatrixDecoder   noiseMatrixDecoder(noiseDecoder, noiseMatrix);
	    
	    int xSeed, ySeed, xClusterSize, yClusterSize;
	    cluster->getSeedCoord(xSeed, ySeed);
	    cluster->getClusterSize(xClusterSize, yClusterSize);
	    vector<float > noiseValues;
	    for ( int yPixel = ySeed - ( yClusterSize / 2 ); yPixel <= ySeed + ( yClusterSize / 2 ); yPixel++ ) {
	      for ( int xPixel = xSeed - ( xClusterSize / 2 ); xPixel <= xSeed + ( xClusterSize / 2 ); xPixel++ ) {
		
		// always check we are still within the sensor!!!
		if ( ( xPixel >= noiseMatrixDecoder.getMinX() )  &&  ( xPixel <= noiseMatrixDecoder.getMaxX() ) &&
		     ( yPixel >= noiseMatrixDecoder.getMinY() )  &&  ( yPixel <= noiseMatrixDecoder.getMaxY() ) ) {
		  int index = noiseMatrixDecoder.getIndexFromXY(xPixel, yPixel);
		  
		  // the corresponding position in the status matrix has to be HITPIXEL
		  // in the EUTelClusteringProcessor, we verify also that
		  // the pixel isHit, but this cannot be done in this
		  // processor, since the status matrix could have been reset
		  // 
		  // bool isHit  = ( statusMatrix->getADCValues()[index] ==
		  // EUTELESCOPE::HITPIXEL );
		  //
		  bool isBad  = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::BADPIXEL );
		  if ( !isBad ) {
		    noiseValues.push_back( noiseMatrix->getChargeValues()[index] );
		  } else {
		    noiseValues.push_back( 0. );
		  }
		} else {
		  noiseValues.push_back( 0. );
		}
	      }
	    }
	    try {
	      cluster->setNoiseValues( noiseValues );
	    } catch ( IncompatibleDataSetException& e ) {
	      message<ERROR> ( log() << e.what() << "\n" << "Continuing without noise based cuts" );
	      _noiseRelatedCuts = false;
	    }
	  } catch ( DataNotAvailableException& e ) {
	    message<ERROR> ( log() << e.what() << "\n" << "Continuing without noise based cuts" );
	    _noiseRelatedCuts = false;
	  }
	  message<DEBUG> ( "Noise related cuts may be used" );
	}
	      
      } else {
	message<ERROR> ("Unknown cluster type. Sorry for quitting");
	throw UnknownDataTypeException("Cluster type unknown");
      }
    
      // increment the event counter 
      _totalClusterCounter[cluster->getDetectorID()]++;
    
      bool isAccepted = true;
   
      isAccepted &= isAboveMinTotalCharge(cluster);
      isAccepted &= isAboveMinTotalSNR(cluster);
      isAccepted &= isAboveNMinCharge(cluster);
      isAccepted &= isAboveNMinSNR(cluster);
      isAccepted &= isAboveNxNMinCharge(cluster);
      isAccepted &= isAboveNxNMinSNR(cluster);
      isAccepted &= isAboveMinSeedCharge(cluster);
      isAccepted &= isAboveMinSeedSNR(cluster);
      isAccepted &= hasQuality(cluster);
      isAccepted &= isBelowMaxClusterNoise(cluster);
      isAccepted &= isInsideROI(cluster); 
      isAccepted &= isOutsideROI(cluster);

      if ( isAccepted )  acceptedClusterVec.push_back(iPulse); 

      delete cluster;

    }

    vector<int >::iterator cluIter = acceptedClusterVec.begin(); 
    while ( cluIter != acceptedClusterVec.end() ) {
      TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl* > (pulseCollectionVec->getElementAt(*cluIter));
      int detectorID = inputDecoder(pulse)["sensorID"];
      clusterNoVec[detectorID]++;
      ++cluIter;
    }

    bool isEventAccepted = areClusterEnough(clusterNoVec) && !areClusterTooMany(clusterNoVec);
    isEventAccepted &= hasSameNumberOfHit(clusterNoVec);

    if ( ! isEventAccepted ) acceptedClusterVec.clear();

    if ( acceptedClusterVec.empty() ) {
      delete filteredCollectionVec;
      message<DEBUG> ( "No cluster passed the selection" );
      if ( _skipEmptyEvent ) throw SkipEventException(this);
      else return;
    } else {
      vector<int >::iterator iter = acceptedClusterVec.begin();
      while ( iter != acceptedClusterVec.end() ) {
	TrackerPulseImpl * pulse     = dynamic_cast<TrackerPulseImpl *> ( pulseCollectionVec->getElementAt( *iter ) );
	TrackerPulseImpl * accepted  = new TrackerPulseImpl;
	accepted->setCellID0( pulse->getCellID0() );
	accepted->setCellID1( pulse->getCellID1() );
	accepted->setTime(    pulse->getTime()    );
	accepted->setCharge(  pulse->getCharge()  );
	accepted->setQuality( pulse->getQuality() );
	accepted->setTrackerData( pulse->getTrackerData() );
	filteredCollectionVec->push_back(accepted);
	_acceptedClusterCounter[ inputDecoder(pulse)["sensorID"] ]++;
	++iter;
      }
      evt->addCollection(filteredCollectionVec, _outputPulseCollectionName);
    }
  } catch (DataNotAvailableException& e ) {
    message<WARNING> ( log() << "Input collection not found in the current event. Skipping..." );
    return;
  }
}  


bool EUTelClusterFilter::areClusterEnough(std::vector<int > clusterNoVec) const {
  
  if ( ! _minClusterNoSwitch ) return true;

  bool areEnough = true;

  for ( int iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) {

    if ( clusterNoVec[iDetector] < _minClusterNoVec[iDetector] ) {
      message<DEBUG> ( log() << "Rejected event because on detector " << iDetector 
		       << " there are " << clusterNoVec[iDetector] << "\nwhile " 
		       << "at least " << _minClusterNoVec[iDetector] << " were requested" );
      _rejectionMap["MinClusterNoCut"][iDetector]++;
      areEnough = false;
    }
  }
  
  return areEnough;

}

bool EUTelClusterFilter::areClusterTooMany(std::vector<int > clusterNoVec) const {
  
  if ( ! _maxClusterNoSwitch ) return false;
  
  bool areTooMany = false;

  for ( int iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) {
    
    if ( clusterNoVec[iDetector] > _maxClusterNoVec[iDetector] ) {
      message<DEBUG> ( log() << "Rejected event because on detector " << iDetector 
		       << " there are " << clusterNoVec[iDetector] << "\nwhile " 
		       << "only " << _maxClusterNoVec[iDetector] << " were allowed" );
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
      message<DEBUG> ( log() << "Rejected event because the number of hit is different from one sensor to the other"  );
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

bool EUTelClusterFilter::isAboveMinTotalCharge(EUTelVirtualCluster * cluster) const {
  
  if ( !_minTotalChargeSwitch ) {
    return true;
  }
  message<DEBUG> ( "Filtering against the total charge " ) ;


  int detectorID = cluster->getDetectorID();
  
  if ( cluster->getTotalCharge() > _minTotalChargeVec[detectorID] ) return true;
  else {
    message<DEBUG> ( log() << "Rejected cluster because its charge is " << cluster->getTotalCharge() 
		     << " and the threshold is " << _minTotalChargeVec[detectorID] ) ;
    _rejectionMap["MinTotalChargeCut"][detectorID]++;
    return false;
  }
}

bool EUTelClusterFilter::isAboveMinTotalSNR(EUTelVirtualCluster * cluster) const {
  
  if ( !_noiseRelatedCuts   ) return true;
  if ( !_minTotalSNRSwitch  ) return true;
  
  int detectorID = cluster->getDetectorID();

  message<DEBUG> ( "Filtering against the minimum total SNR " ) ;
  if  ( cluster->getClusterSNR() > _minTotalSNRVec[detectorID] ) return true;
  else {
    message<DEBUG> ( log() << "Rejected cluster because its SNR is " << cluster->getClusterSNR() 
		     << " and the threshold is " << _minTotalSNRVec[detectorID] ) ;
    _rejectionMap["MinTotalSNRCut"][detectorID]++;
    return false;
  }
}
  
bool EUTelClusterFilter::isAboveNMinCharge(EUTelVirtualCluster * cluster) const {
  
  if ( !_minNChargeSwitch ) return true;
  
  message<DEBUG> ( "Filtering against the N Pixel charge " ) ;

  int detectorID = cluster->getDetectorID();

  vector<float >::const_iterator iter = _minNChargeVec.begin();
  while ( iter != _minNChargeVec.end() ) {
    int nPixel      = static_cast<int > (*iter);
    float charge    = cluster->getClusterCharge(nPixel);
    float threshold = (* (iter + detectorID + 1) );
    if ( charge > threshold ) {
      iter += _noOfDetectors + 1;
    } else {
      message<DEBUG> ( log() << "Rejected cluster because its charge over " << (*iter) << " is " << charge 
		       << " and the threshold is " << threshold );
      _rejectionMap["MinNChargeCut"][detectorID]++;
      return false;
    }
  }
  return true;
}
    

bool EUTelClusterFilter::isAboveNMinSNR(EUTelVirtualCluster * cluster) const {

  if ( !_noiseRelatedCuts ) return true;
  if ( !_minNSNRSwitch    ) return true;
  
  message<DEBUG> ( "Filtering against the N pixel SNR " );
  
  int detectorID = cluster->getDetectorID();
  
  vector<float >::const_iterator iter = _minNSNRVec.begin();
  while ( iter !=  _minNSNRVec.end() ) {
    int nPixel      = static_cast<int > (*iter);
    float SNR       = cluster->getClusterSNR(nPixel);
    float threshold = (* (iter + detectorID + 1 ) );
    if ( SNR > threshold ) {
      iter += _noOfDetectors + 1;
    } else {
      message<DEBUG> ( log() << "Rejected cluster because its SNR over " << (*iter) << " is " << SNR 
		       << " and the threshold is " << threshold );
      _rejectionMap["MinNSNRCut"][detectorID]++;
      return false;
    }
  }
  return true;
}




bool EUTelClusterFilter::isAboveNxNMinCharge(EUTelVirtualCluster * cluster) const {
  
  if ( !_minNxNChargeSwitch ) return true;
  
  message<DEBUG> ( "Filtering against the N x N pixel charge" );
  
  int detectorID = cluster->getDetectorID();
  vector<float >::const_iterator iter = _minNxNChargeVec.begin();
  while ( iter != _minNxNChargeVec.end() ) {
    int nxnPixel    = static_cast<int > ( *iter ) ;
    float charge    = cluster->getClusterCharge(nxnPixel, nxnPixel);
    float threshold = (* ( iter + detectorID + 1 )) ;
    if ( ( threshold <= 0) || (charge > threshold) ) {
      iter += _noOfDetectors + 1;
    } else {
      message<DEBUG> ( log() << "Rejected cluster because its charge within a " << (*iter) << " x " << (*iter) 
		       << " subcluster is " << charge << " and the threshold is " << threshold ) ;
      _rejectionMap["MinNxNChargeCut"][detectorID]++;
      return false;
    }
  }
  return true;
}


bool EUTelClusterFilter::isAboveNxNMinSNR(EUTelVirtualCluster * cluster) const {
  
  if ( !_noiseRelatedCuts  ) return true;
  if ( !_minNxNSNRSwitch   ) return true;

  message<DEBUG> ( "Filtering against the N x N pixel charge" );
  
  int detectorID = cluster->getDetectorID();
  vector<float >::const_iterator iter = _minNxNSNRVec.begin();
  while ( iter != _minNxNSNRVec.end() ) {
    int nxnPixel    = static_cast<int > ( *iter ) ;
    float snr       = cluster->getClusterSNR(nxnPixel, nxnPixel);
    float threshold = (* ( iter + detectorID + 1 )) ;
    if ( ( threshold <= 0) || (snr > threshold) ) {
      iter += _noOfDetectors + 1;
    } else {
      message<DEBUG> ( log() << "Rejected cluster because its SNR within a " << (*iter) << " x " << (*iter) 
		       << " subcluster is " << snr << " and the threshold is " << threshold ) ;
      _rejectionMap["MinNxNSNRCut"][detectorID]++;
      return false;
    }
  }
  return true;

}

bool EUTelClusterFilter::isAboveMinSeedCharge(EUTelVirtualCluster * cluster) const {
  
  if ( !_minSeedChargeSwitch ) return true;

  message<DEBUG> ( "Filtering against the seed charge " ) ;
  
  int detectorID = cluster->getDetectorID();
  
  if ( cluster->getSeedCharge() > _minSeedChargeVec[detectorID] ) return true;
  else {
    message<DEBUG> (  log() << "Rejected cluster because its seed charge is " << cluster->getSeedCharge()
		      << " and the threshold is " <<  _minSeedChargeVec[detectorID] );
    _rejectionMap["MinSeedChargeCut"][detectorID]++;
    return false;
  }
}

bool EUTelClusterFilter::isAboveMinSeedSNR(EUTelVirtualCluster * cluster) const {

  if ( !_noiseRelatedCuts  ) return true;
  if ( !_minSeedSNRSwitch  ) return true;

  message<DEBUG> ( "Filtering against the seed SNR " ); 

  int detectorID = cluster->getDetectorID();
  if ( cluster->getSeedSNR() > _minSeedSNRVec[detectorID] ) return true;
  else {
    message<DEBUG> ( log() << "Rejected cluster because its seed charge is " << cluster->getSeedSNR()
		     << " and the threshold is " <<  _minSeedSNRVec[detectorID] );
    _rejectionMap["MinSeedSNRCut"][detectorID]++;
    return false;
  }
}



bool EUTelClusterFilter::hasQuality(EUTelVirtualCluster * cluster) const {
  
  if ( !_clusterQualitySwitch ) return true;
  
  int detectorID = cluster->getDetectorID();
  
  if ( _clusterQualityVec[detectorID] < 0 ) return true;
  
  ClusterQuality actual = cluster->getClusterQuality();
  ClusterQuality needed = static_cast<ClusterQuality> ( _clusterQualityVec[detectorID] );
  
  if ( actual == needed ) return true;
  else {
    message<DEBUG> ( log() << "Rejected cluster because its quality " << static_cast<int> (actual)
		     << " is not " << _clusterQualityVec[detectorID] );
    _rejectionMap["ClusterQualityCut"][detectorID]++;
    return false;
  }
}    

bool EUTelClusterFilter::isBelowMaxClusterNoise(EUTelVirtualCluster * cluster) const {
  
  if ( !_noiseRelatedCuts       ) return true;
  if ( !_maxClusterNoiseSwitch  ) return true;

  message<DEBUG> ( "Filtering against the maximum cluster noise"  ) ;
  int detectorID = cluster->getDetectorID();
  if (  ( cluster->getClusterNoise() < _maxClusterNoiseVec[detectorID] ) ||
	( _maxClusterNoiseVec[detectorID] < 0 ) ) return true;
  else {
    message<DEBUG> ( log() << "Rejected cluster because its noise is " << cluster->getClusterNoise()
		     << " and the threshold is " <<  _maxClusterNoiseVec[detectorID] );
    _rejectionMap["MaxClusterNoiseCut"][detectorID]++;
    return false;
  }
}


bool EUTelClusterFilter::isInsideROI(EUTelVirtualCluster * cluster) const {
  
  if ( !_insideROISwitch ) return true;
  
  int detectorID = cluster->getDetectorID();
  float x, y;
  cluster->getCenterOfGravity(x, y);

  bool tempAccepted = true;
  vector<EUTelROI>::const_iterator iter = _insideROIVec.begin();
  vector<EUTelROI>::const_iterator end  = _insideROIVec.end();
  while ( true ) {
    iter = find_if( iter, end, HasSameID(detectorID));
    if ( iter != end ) {
      if ( (*iter).isInside(x,y) ) {
	tempAccepted &= true;
      }  else {
	_rejectionMap["InsideROICut"][detectorID]++;
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
	_rejectionMap["OutsideROICut"][detectorID]++;
	tempAccepted &= false;
      }
      ++iter;
    } else break;
  }
  return tempAccepted;
  
}

void EUTelClusterFilter::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelClusterFilter::end() {
  message<MESSAGE> ( log() << printSummary() );
  
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
