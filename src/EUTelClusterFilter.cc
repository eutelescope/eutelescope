// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelClusterFilter.cc,v 1.6 2007-07-10 14:31:15 bulgheroni Exp $
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

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h> 
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <limits>

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


  FloatVec minSeedChargeVecExample;
  minSeedChargeVecExample.push_back(20);
  minSeedChargeVecExample.push_back(25);
  minSeedChargeVecExample.push_back(21);
  
  registerProcessorParameter("SeedMinCharge", "This is the minimum allowed charge that the seed pixel of a cluster has to have.\n"
			     "One floating number for each detector",
			     _minSeedChargeVec, minSeedChargeVecExample);

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

  // NxN cluster charge
  if (  ( count_if( _minNxNChargeVec.begin(), _minNxNChargeVec.end(), bind2nd(greater<float>(),0)) != 0  ) &&
	( _minNxNChargeVec[0] != 0 ) ) {
    _minNxNChargeSwitch = true;
  } else {
    _minNxNChargeSwitch = false;
  }
	 

  // Seed charge
  if ( count_if( _minSeedChargeVec.begin(), _minSeedChargeVec.end(), bind2nd(greater<float>(), 0) ) != 0 ) {
    _minSeedChargeSwitch = true;
  } else {
    _minSeedChargeSwitch = false;
  }

  message<DEBUG> ( log () << "MinSeedChargeSwitch " << _minSeedChargeSwitch );

  // Cluster quality
  if ( count_if(_clusterQualityVec.begin(), _clusterQualityVec.end(), bind2nd(greater_equal<int>(), 0) ) != 0 ) {
    _clusterQualitySwitch = true;
  } else {
    _clusterQualitySwitch = false;
  }

  message<DEBUG> ( log () << "ClusterQualitySwitch " << _clusterQualitySwitch );

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
    
    if ( _minNChargeSwitch ) {
      unsigned int module = _noOfDetectors + 1;
      if ( _minNChargeVec.size() % module != 0 ) {
	message<ERROR> (log() << "The threshold vector for the N pixels charge did not match the right size \n "
			<<  "The number of planes is " << _noOfDetectors << " while the thresholds are " << _minSeedChargeVec.size() 
			<<  "\n"
			<<  "Disabling the selection criterion and continue without");
	_minNChargeSwitch = false;
      } else {
	message<DEBUG> ("N pixel charge criterion verified and switched on");
	vector<unsigned int> rejectedCounter(_noOfDetectors, 0);
	_rejectionMap.insert( make_pair("MinNChargeCut", rejectedCounter));
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
    
    if ( type == kEUTelFFClusterImpl ) 
      cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData() ) );
    else {
      message<ERROR> ("Unknown cluster type. Sorry for quitting");
      throw UnknownDataTypeException("Cluster type unknown");
    }
    
    bool isAccepted = true;
    
    isAccepted &= isAboveMinTotalCharge(cluster);
    isAccepted &= isAboveNMinCharge(cluster);
    isAccepted &= isAboveNxNMinCharge(cluster);
    isAccepted &= isAboveMinSeedCharge(cluster);
    isAccepted &= hasQuality(cluster);
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

  if ( ! isEventAccepted ) acceptedClusterVec.clear();

  if ( acceptedClusterVec.empty() ) {
    delete filteredCollectionVec;
    message<DEBUG> ( "No cluster passed the selection" );
    throw SkipEventException(this);
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
      ++iter;
    }
    evt->addCollection(filteredCollectionVec, _outputPulseCollectionName);
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
    unsigned int current = _rejectionMap["MinTotalChargeCut"][detectorID];
    _rejectionMap["MinTotalChargeCut"][detectorID] = current + 1;
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

string  EUTelClusterFilter::printSummary() const {
  
  stringstream  ss;

  ss << "========================================================\n" 
     << " Rejection summary \n"
     << "========================================================\n"  ;
  
  map<string, vector<unsigned int> >::iterator iter = _rejectionMap.begin();
  while ( iter != _rejectionMap.end() ) {
    ss << " " << (*iter).first << "\t";
    
    vector<unsigned int>::iterator iter2 = (*iter).second.begin();
    while ( iter2 != (*iter).second.end() ) {
      ss << ( *iter2)  << "  ";
      ++iter2;
    }
    ss << "\n"
       << "--------------------------------------------------------\n";
    
    ++iter;
  }
  ss << "\n"
     << "========================================================\n" ;

  return ss.str();
    
}
