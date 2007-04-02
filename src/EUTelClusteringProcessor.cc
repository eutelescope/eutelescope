// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelClusteringProcessor.cc,v 1.6 2007-04-02 14:21:10 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

//
//   Because of the complexity of this processor some debug features
//   have been implemented and left within the code, so that other
//   users can test the processor functionality. To active the debug
//   mode replace all "///" with "" and re-build.  To shut down the
//   debug option, just replace the string "/*DEBUG*/" with "///
//   /*DEBUG*/"
//   Emacs Alt+Shift+5 does it very easily!
//


// eutelescope includes ".h" 
#include "EUTELESCOPE.h"
#include "EUTelClusteringProcessor.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h> 
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
/// /*DEBUG*/ #include <fstream> 

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

  /// /*DEBUG*/ ofstream logfile;

EUTelClusteringProcessor::EUTelClusteringProcessor () :Processor("EUTelClusteringProcessor") {

  // modify processor description
  _description =
    "EUTelClusteringProcessor subtract the pedestal value from the input data";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "DataCollectionName",
			   "Input calibrated data collection name",
			   _dataCollectionName, string ("data"));

  registerInputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
			   "Noise (input) collection name",
			   _noiseCollectionName, string("noise"));

  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
			   "Pixel status (input) collection name",
			   _statusCollectionName, string("status"));

  registerOutputCollection(LCIO::TRACKERDATA, "ClusterCollectionName",
			   "Cluster (output) collection name",
			   _clusterCollectionName, string("cluster"));


  // now the optional parameters
  registerProcessorParameter ("ClusteringAlgo",
			      "Select here which algorithm should be used for clustering",
			      _clusteringAlgo, string(EUTELESCOPE::FIXEDFRAME));
  
  registerProcessorParameter ("ClusterSizeX",
			      "Maximum allowed cluster size along x (only odd numbers)",
			      _xClusterSize, static_cast<int> (5));

  registerProcessorParameter ("ClusterSizeY",
			      "Maximum allowed cluster size along y (only odd numbers)",
			      _yClusterSize, static_cast<int> (5));

  registerProcessorParameter ("SeedPixelCut",
			      "Threshold in SNR for seed pixel identification",
			      _seedPixelCut, static_cast<float> (4.5));

  registerProcessorParameter ("ClusterCut",
			      "Threshold in SNR for cluster identification",
			      _clusterCut, static_cast<float> (3.0));

}


void EUTelClusteringProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // in the case the FIXEDFRAME algorithm is selected, the check if
  // the _xClusterSize and the _yClusterSize are odd numbers
  if ( _clusteringAlgo == EUTELESCOPE::FIXEDFRAME ) {
    bool isZero = ( _xClusterSize <= 0 );
    bool isEven = ( _xClusterSize % 2 == 0 );
    if ( isZero || isEven ) {
      throw InvalidParameterException("_xClusterSize has to be positive and odd");
    }
    isZero = ( _yClusterSize <= 0 );
    isEven = ( _yClusterSize % 2 == 0 );
    if ( isZero || isEven ) {
      throw InvalidParameterException("_yClusterSize has to be positive and odd");
    }
  }

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // reset the content of the total cluster vector
  _totCluster.clear();

}

void EUTelClusteringProcessor::processRunHeader (LCRunHeader * rdr) {

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  EUTelRunHeaderImpl *  runHeader = static_cast<EUTelRunHeaderImpl*>(rdr);

  // the four vectors containing the first and the last pixel
  // along both the directions
  _minX = runHeader->getMinX();
  _maxX = runHeader->getMaxX();
  _minY = runHeader->getMinY();
  _maxY = runHeader->getMaxY();

  // increment the run counter
  ++_iRun;

}


void EUTelClusteringProcessor::processEvent (LCEvent * evt) {

  if (_iEvt % 10 == 0) 
    cout << "[" << name() << "] Clustering event " << _iEvt << endl;
  
  if ( _clusteringAlgo == EUTELESCOPE::FIXEDFRAME ) fixedFrameClustering(evt);
  
}

void EUTelClusteringProcessor::fixedFrameClustering(LCEvent * evt) {
  
  LCCollectionVec * inputCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_dataCollectionName));
  LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
  LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
  
  if (isFirstEvent()) {
    
      /// /*DEBUG*/ logfile.open("clustering.log");
    
    // this is the right place to cross check whether the pedestal and
    // the input data are at least compatible. I mean the same number
    // of detectors and the same number of pixels in each place.
    
    if  ( inputCollectionVec->getNumberOfElements() != noiseCollectionVec->getNumberOfElements()) {
      stringstream ss;
      ss << "Input data and pedestal are incompatible" << endl
	 << "Input collection has    " << inputCollectionVec->getNumberOfElements()    << " detectors," << endl
	 << "Noise collection has    " << noiseCollectionVec->getNumberOfElements() << " detectors," << endl;
      throw IncompatibleDataSetException(ss.str());
    }
    
    for (_iDetector = 0; _iDetector < inputCollectionVec->getNumberOfElements(); _iDetector++) 
      _totCluster.push_back(0);
    
    _isFirstEvent = false;
  }
  
    /// /*DEBUG*/ logfile << "Event " << _iEvt << endl;
  
  LCCollectionVec * clusterCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  
  for (_iDetector = 0; _iDetector < inputCollectionVec->getNumberOfElements(); _iDetector++) {
    
      /// /*DEBUG*/ logfile << "  Working on detector " << _iDetector << endl;
    
    // get the calibrated data 
    TrackerDataImpl    * data   = dynamic_cast<TrackerDataImpl*>   (inputCollectionVec->getElementAt(_iDetector));
    TrackerDataImpl    * noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt(_iDetector));
    TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt(_iDetector));

    // reset the status
    resetStatus(status);

    // initialize the cluster counter 
    short clusterCounter = 0;

    _seedCandidateMap.clear();
    
    for (unsigned int iPixel = 0; iPixel < data->getChargeValues().size(); iPixel++) {
      if (status->getADCValues()[iPixel] == EUTELESCOPE::GOODPIXEL) {
	if (data->getChargeValues()[iPixel] > _seedPixelCut * noise->getChargeValues()[iPixel]) {
	  _seedCandidateMap.insert(make_pair(data->getChargeValues()[iPixel], iPixel));
	}
      }
    }

    // continue only if seed candidate map is not empty!
    if ( _seedCandidateMap.size() != 0 ) {

        /// /*DEBUG*/ logfile << "  Seed candidates " << _seedCandidateMap.size() << endl;

      // now built up a cluster for each seed candidate 
      map<float, unsigned int>::iterator mapIter = _seedCandidateMap.end();     
      while ( mapIter != _seedCandidateMap.begin() ) {
	--mapIter;	
	// check if this seed candidate has not been already added to a
	// cluster
	if ( status->adcValues()[(*mapIter).second] == EUTELESCOPE::GOODPIXEL ) {
	  // if we enter here, this means that at least the seed pixel
	  // wasn't added yet to another cluster.  Note that now we need
	  // to build a candidate cluster that has to pass the
	  // clusterCut to be considered a good cluster
	  double clusterCandidateSignal    = 0.;
	  double clusterCandidateNoise2    = 0.;
	  FloatVec clusterCandidateCharges;
	  IntVec   clusterCandidateIndeces;
	  int seedX, seedY;
	  getXYFromIndex((*mapIter).second, seedX, seedY);

	  // start looping around the seed pixel. Remember that the seed
	  // pixel has to stay in the center of cluster
	  ClusterQuality cluQuality = kGoodCluster;
	  for (int yPixel = seedY - (_yClusterSize / 2); yPixel <= seedY + (_yClusterSize / 2); yPixel++) {
	    for (int xPixel =  seedX - (_xClusterSize / 2); xPixel <= seedX + (_xClusterSize / 2); xPixel++) {
	      // always check we are still within the sensor!!!
	      if ( ( xPixel >= _minX[_iDetector] )  &&  ( xPixel <= _maxX[_iDetector] ) &&
		   ( yPixel >= _minY[_iDetector] )  &&  ( yPixel <= _maxY[_iDetector] ) ) {
		int index = getIndexFromXY(xPixel, yPixel);
		bool isHit  = ( status->getADCValues()[index] == EUTELESCOPE::HITPIXEL  );
		bool isGood = ( status->getADCValues()[index] == EUTELESCOPE::GOODPIXEL );
		if ( isGood && !isHit ) {
		  clusterCandidateSignal += data->getChargeValues()[index];
		  clusterCandidateNoise2 += noise->getChargeValues()[index];
		  clusterCandidateCharges.push_back(data->getChargeValues()[index]);
		  clusterCandidateIndeces.push_back(index);
		} else if (isHit) {
		  // this can be a good place to flag the current
		  // cluster as kMergedCluster, but it would introduce
		  // a bias since the at least another cluster (the
		  // one which this pixel belong to) is not flagged.
		  //
		  // In order to flag all merged clusters and possibly
		  // try to separate the different contributions use
		  // the EUTelSeparateClusterProcessor
		  clusterCandidateCharges.push_back(0.);
		} else if (!isGood) {
		  cluQuality = cluQuality | kIncompleteCluster;
		  clusterCandidateCharges.push_back(0.);
		}
	      } else {
		cluQuality = cluQuality | kBorderCluster;
		clusterCandidateCharges.push_back(0.);
	      }
	    }
	  }
	  // at this point we have built the cluster candidate,
	  // we need to validate it
	  if ( clusterCandidateSignal > _clusterCut * sqrt(clusterCandidateNoise2) ) {
	    // the cluster candidate is a good cluster
	    // mark all pixels belonging to the cluster as hit
	    IntVec::iterator indexIter = clusterCandidateIndeces.begin();
	    TrackerDataImpl * cluster = new TrackerDataImpl;
	    CellIDEncoder<TrackerDataImpl> idClusterEncoder(EUTELESCOPE::CLUSTERDEFAULTENCODING, clusterCollection);
	    idClusterEncoder["sensorID"]      = _iDetector;
	    idClusterEncoder["clusterID"]     = clusterCounter;
	    idClusterEncoder["xSeed"]         = seedX;
	    idClusterEncoder["ySeed"]         = seedY;
	    idClusterEncoder["xCluSize"]      = _xClusterSize;
	    idClusterEncoder["yCluSize"]      = _yClusterSize;
	    idClusterEncoder["quality"]       = static_cast<int>(cluQuality);
	    idClusterEncoder.setCellID(cluster);
	    
	      /// /*DEBUG*/ logfile << "  Cluster no " <<  clusterCounter << " seedX " << seedX << " seedY " << seedY << endl;
	    
	    while ( indexIter != clusterCandidateIndeces.end() ) {
	      status->adcValues()[(*indexIter)] = EUTELESCOPE::HITPIXEL;
	      ++indexIter;
	    }

	      /// /*DEBUG*/ for (unsigned int iPixel = 0; iPixel < clusterCandidateIndeces.size(); iPixel++) {
	      /// /*DEBUG*/  logfile << "  x " <<  getXFromIndex(clusterCandidateIndeces[iPixel])
	      /// /*DEBUG*/	      << "  y " <<  getYFromIndex(clusterCandidateIndeces[iPixel])
	      /// /*DEBUG*/              << "  s " <<  clusterCandidateCharges[iPixel] << endl;
	      /// /*DEBUG*/ }

	    // copy the candidate charges inside the cluster
	    cluster->setChargeValues(clusterCandidateCharges);
	    clusterCollection->push_back(cluster);
	    _totCluster[_iDetector] += 1;
	    ++clusterCounter;

	  } else {
	    // the cluster has not passed the cut!
	  }
	}
      }
    }
  }
  evt->addCollection(clusterCollection,_clusterCollectionName);
  
  ++_iEvt;
  
}



void EUTelClusteringProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}
 

void EUTelClusteringProcessor::end() {
  cout << "[" << name() <<"] Successfully finished" << endl;

  for (_iDetector = 0; _iDetector < (signed) _totCluster.size() ; _iDetector++) {
    cout << "Found " << _totCluster[_iDetector] << " clusters on detector " << _iDetector << endl;
     /// /*DEBUG*/ logfile << "Found " << _totCluster[_iDetector] << " clusters on detector " << _iDetector << endl;
  }
    /// /*DEBUG*/  logfile.close();
}


void EUTelClusteringProcessor::resetStatus(IMPL::TrackerRawDataImpl * status) {
  
  ShortVec::iterator iter = status->adcValues().begin();
  while ( iter != status->adcValues().end() ) {
    if ( *iter == EUTELESCOPE::HITPIXEL ) {
      *iter = EUTELESCOPE::GOODPIXEL;
    }
    ++iter; 
  }

}
