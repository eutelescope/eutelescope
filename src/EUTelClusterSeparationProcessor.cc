// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelClusterSeparationProcessor.cc,v 1.4 2007-05-21 11:46:24 bulgheroni Exp $
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
#include "EUTelFFClusterImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelClusterSeparationProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <string>
#include <set>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelClusterSeparationProcessor::EUTelClusterSeparationProcessor () :Processor("EUTelClusterSeparationProcessor") {

  // modify processor description
  _description = "EUTelClusterSeparationProcessor separates merging clusters";


  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "ClusterCollectionName",
			   "Cluster collection name ",
			   _clusterCollectionName, string ("cluster"));

  // now the optional parameters
  registerProcessorParameter ("SeparationAlgorithm",
			      "Select which algorithm to use for cluster separation",
			      _separationAlgo, string(EUTELESCOPE::FLAGONLY));

  registerProcessorParameter ("MinimumDistance",
			      "Minimum distance allowed between separated clusters (0 == only touching clusters)",
			      _minimumDistance, static_cast<float> (0));
    
}


void EUTelClusterSeparationProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

}

void EUTelClusterSeparationProcessor::processRunHeader (LCRunHeader * rdr) {

  // increment the run counter
  ++_iRun;

}


void EUTelClusterSeparationProcessor::processEvent (LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> ( event );
  if ( evt->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do.");
    return ;
  }
  
  if (_iEvt % 10 == 0) 
    message<MESSAGE> ( log() << "Separating clusters on event " << _iEvt ) ;

  LCCollectionVec       *   clusterCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_clusterCollectionName));
  vector< pair<int, int > > mergingPairVector;

  for ( int iCluster = 0 ; iCluster < clusterCollectionVec->getNumberOfElements() ; iCluster++) {

    EUTelFFClusterImpl * cluster = static_cast< EUTelFFClusterImpl *> (clusterCollectionVec->getElementAt(iCluster));
    
    int  iOtherCluster     = iCluster + 1;
    bool isExisisting      = (iOtherCluster < clusterCollectionVec->getNumberOfElements() );
    bool isOnSameDetector  = true;

    
    while ( isOnSameDetector && isExisisting ) {

      // get the next cluster in the collection
      EUTelFFClusterImpl * otherCluster = static_cast< EUTelFFClusterImpl *> (clusterCollectionVec->getElementAt(iOtherCluster));
      
      // check if the two are on the same detector
      if ( cluster->getDetectorID() == otherCluster->getDetectorID() ) {
	
	if ( _minimumDistance == 0 ) {
	  // ok we need to calculate the touching distance
	  float radius      = cluster->getExternalRadius();
	  float otherRadius = otherCluster->getExternalRadius();
	  _minimumDistance  = radius + otherRadius;
	}

	// ok they are on the same plane, so it makes sense check it
	// they are merging
	float distance = cluster->getDistance(otherCluster);
	
	if ( distance < _minimumDistance ) {
	  // they are merging! we need to apply the separation
	  // algorithm
	  mergingPairVector.push_back( make_pair(iCluster, iOtherCluster) );
	}

      } else {
	isOnSameDetector = false;
      }
      isExisisting = (++iOtherCluster <  clusterCollectionVec->getNumberOfElements() );
    }  
  }
  
  // at this point we have inserted into the mergingPairVector all the
  // pairs of merging clusters. we can try to put together all groups
  // of clusters, but only in the case the mergingPairVector has a non
  // null size
  if ( mergingPairVector.empty() ) {
    ++_iEvt;
    return ;
  }
  
  // all merging clusters are collected into a vector of set. Each set
  // is a group of clusters all merging.
  vector< set<int > > mergingSetVector;
  groupingMergingPairs(mergingPairVector, &mergingSetVector) ;

  applySeparationAlgorithm(mergingSetVector, clusterCollectionVec);

  ++_iEvt;
  
}
  

bool EUTelClusterSeparationProcessor::applySeparationAlgorithm(std::vector<std::set <int > > setVector, 
							       LCCollectionVec * collectionVec) const {

  //  message<DEBUG> ( log() << "Applying cluster separation algorithm " << _separationAlgo );
  message<DEBUG> ( log () <<  "Found " ); // << setVector.size() << " group(s) of merging clusters on event " << _iEvt );
  if ( _separationAlgo == EUTELESCOPE::FLAGONLY ) {

#ifdef MARLINDEBUG
    int iCounter = 0;    
#endif

    vector<set <int > >::iterator vectorIterator = setVector.begin();    
    while ( vectorIterator != setVector.end() ) {

#ifdef MARLINDEBUG
      message<DEBUG> ( log() <<  "     Group " << (iCounter++) << " with the following clusters " << endl;
#endif

      set <int >::iterator setIterator = (*vectorIterator).begin();
      while ( setIterator != (*vectorIterator).end() ) {
	EUTelFFClusterImpl * cluster = static_cast<EUTelFFClusterImpl *> (collectionVec->getElementAt( *setIterator )) ;

#ifdef MARLINDEBUG
	int xSeed, ySeed;
	int detectorID = cluster->getDetectorID();
	int clusterID  = cluster->getClusterID();
	int xSize, ySize;
	ClusterQuality quality = cluster->getClusterQuality();
	cluster->getClusterSize(xSize, ySize);
	cluster->getSeedCoord(xSeed, ySeed);
	message<DEBUG> ( log()  << "         Cluster " << (*setIterator)  << " (" << detectorID << ":" << clusterID << ":" << xSeed 
			 << ", " << ySeed << ":" << xSize << "," << ySize << ":" << static_cast<int>(quality) << ") " );
#endif

	cluster->setClusterQuality ( cluster->getClusterQuality() | kMergedCluster );
	++setIterator;
      }
      ++vectorIterator;
    }
    return true;
  }

  return false;


}

void EUTelClusterSeparationProcessor::groupingMergingPairs(std::vector< std::pair<int , int> > pairVector, 
							   std::vector< std::set<int > > * setVector) const {

  vector< pair<int, int> >::iterator iter = pairVector.begin();
  while ( iter != pairVector.end() ) {

    set<int > tempSet;
    // add the pair to the tempSet
    tempSet.insert(pairVector.front().first);
    tempSet.insert(pairVector.front().second);

    vector<vector< pair<int, int> >::iterator > tempIterVec;
    tempIterVec.push_back(iter);
     
    if ( iter + 1 != pairVector.end() ) {
      vector< pair<int, int> >::iterator otherIter = iter + 1;
      while (otherIter != pairVector.end() ) {
 	if ( ( tempSet.find(otherIter->first)  != tempSet.end() ) ||
 	     ( tempSet.find(otherIter->second) != tempSet.end() ) ) {
 	  tempSet.insert(otherIter->first);
 	  tempSet.insert(otherIter->second);
 	  tempIterVec.push_back(otherIter);
	}
	++otherIter;
      }
    }
    setVector->push_back(tempSet);

    // remove the pairs already grouped starting from the last
    vector<vector< pair<int, int> >::iterator >::reverse_iterator iterIter = tempIterVec.rbegin();
    while ( iterIter != tempIterVec.rend() ) {
      pairVector.erase(*iterIter);
      ++iterIter;
    }
  }
}

void EUTelClusterSeparationProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelClusterSeparationProcessor::end() {
  message<MESSAGE> ( "Successfully finished" );
}

