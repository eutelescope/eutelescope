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
#include "EUTelFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelClusterSeparationProcessor.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <string>
#include <set>
#include <memory>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelClusterSeparationProcessor::EUTelClusterSeparationProcessor () :Processor("EUTelClusterSeparationProcessor") {

  // modify processor description
  _description = "EUTelClusterSeparationProcessor separates merging clusters";


  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERPULSE, "ClusterCollectionName",
                           "Cluster collection name ",
                           _clusterCollectionName, string ("cluster"));

  // and the output collection
  registerOutputCollection (LCIO::TRACKERPULSE, "ClusterOutputCollectionName",
                            "Cluster output collection name",
                            _clusterOutputCollectionName, string ("splitcluster" ));

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

  // set to zero the run counter
  _iRun = 0;

}

void EUTelClusterSeparationProcessor::processRunHeader (LCRunHeader * rdr) {

  auto_ptr<EUTelRunHeaderImpl > runHeader ( new EUTelRunHeaderImpl( rdr ) ) ;
  runHeader->addProcessor( type() );

  // increment the run counter
  ++_iRun;

}


void EUTelClusterSeparationProcessor::processEvent (LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> ( event );

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  LCCollectionVec * clusterCollectionVec;

  try {
    clusterCollectionVec  =  dynamic_cast <LCCollectionVec *> (evt->getCollection(_clusterCollectionName));
  } catch (lcio::DataNotAvailableException& e ) {
    streamlog_out ( DEBUG2 ) << "No cluster collection found on this event. Returnin now" << endl;
    return ;
  }


  LCCollectionVec * outputCollectionVec  =  new LCCollectionVec(LCIO::TRACKERPULSE);
  CellIDEncoder<TrackerPulseImpl> outputEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, outputCollectionVec);
  CellIDDecoder<TrackerPulseImpl> cellDecoder(clusterCollectionVec);

  vector< pair<int, int > >       mergingPairVector;

  for ( int iCluster = 0 ; iCluster < clusterCollectionVec->getNumberOfElements() ; iCluster++) {

    TrackerPulseImpl   * pulse   = dynamic_cast<TrackerPulseImpl *>   ( clusterCollectionVec->getElementAt(iCluster) );
    ClusterType          type    = static_cast<ClusterType> (static_cast<int>( cellDecoder(pulse)["type"] ) );

    // all clusters have to inherit from the virtual cluster (that is
    // a TrackerDataImpl with some utility methods).
    EUTelVirtualCluster    * cluster;

    if ( type == kEUTelFFClusterImpl )
    {
      cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) ) ;
    }
    else if ( type == kEUTelBrickedClusterImpl )
    {
      cluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) ) ;
    }
    else if ( type == kEUTelSparseClusterImpl ) {

      // ok the cluster is of sparse type, but we also need to know
      // the kind of pixel description used. This information is
      // stored in the corresponding original data collection.

      LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
      TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
      CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
      SparsePixelType pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));

      // now we know the pixel type. So we can properly create a new
      // instance of the sparse cluster
      if ( pixelType == kEUTelGenericSparsePixel ) {
        cluster = new EUTelSparseClusterImpl< EUTelGenericSparsePixel >
          ( static_cast<TrackerDataImpl *> ( pulse->getTrackerData()  ) );
      } else {
        streamlog_out ( ERROR4 ) << "Unknown pixel type. Sorry for quitting." << endl;
        throw UnknownDataTypeException("Pixel type unknown");
      }

    } else {
      streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
      throw UnknownDataTypeException("Cluster type unknown");
    }

    int  iOtherCluster     = iCluster + 1;
    bool isExisisting      = (iOtherCluster < clusterCollectionVec->getNumberOfElements() );
    bool isOnSameDetector  = true;


    while ( isOnSameDetector && isExisisting ) {

      // get the next cluster in the collection
      TrackerPulseImpl    * otherPulse   = dynamic_cast<TrackerPulseImpl *> (clusterCollectionVec->getElementAt(iOtherCluster)) ;
      EUTelVirtualCluster * otherCluster;

      if ( type == kEUTelFFClusterImpl )
        otherCluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (otherPulse->getTrackerData()) );
      else if ( type == kEUTelBrickedClusterImpl )
        otherCluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*> (otherPulse->getTrackerData()) );
      else
      {
        streamlog_out ( ERROR4 ) << "Unknown cluster type. Sorry for quitting" << endl;
        throw UnknownDataTypeException("Cluster type unknown");
      }

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
      delete otherCluster;
    }

    delete cluster;
  }

  // at this point we have inserted into the mergingPairVector all the
  // pairs of merging clusters. we can try to put together all groups
  // of clusters, but only in the case the mergingPairVector has a non
  // null size
  if ( mergingPairVector.empty() ) {
    // if the mergingPairVector is empty, then there are no merging
    // clusters, i.e. that the input and the output collections are
    // exactly the same
    for ( int iPulse = 0; iPulse < clusterCollectionVec->getNumberOfElements() ; iPulse++ ) {
      TrackerPulseImpl * pulse    = dynamic_cast<TrackerPulseImpl *> ( clusterCollectionVec->getElementAt( iPulse ) );
      TrackerPulseImpl * newPulse = new TrackerPulseImpl;
      newPulse->setCellID0( pulse->getCellID0() );
      newPulse->setCellID1( pulse->getCellID1() );
      newPulse->setTime   ( pulse->getTime() );
      newPulse->setCharge ( pulse->getCharge() );
      newPulse->setQuality( pulse->getQuality() );
      newPulse->setTrackerData( pulse->getTrackerData() );

      outputCollectionVec->push_back( newPulse );
    }
    evt->addCollection( outputCollectionVec, _clusterOutputCollectionName );
    return ;
  }

  // all merging clusters are collected into a vector of set. Each set
  // is a group of clusters all merging.
  vector< set<int > > mergingSetVector;
  groupingMergingPairs(mergingPairVector, &mergingSetVector) ;

  applySeparationAlgorithm(mergingSetVector, clusterCollectionVec, outputCollectionVec);
  evt->addCollection( outputCollectionVec, _clusterOutputCollectionName );


}


bool EUTelClusterSeparationProcessor::applySeparationAlgorithm(std::vector<std::set <int > > setVector,
                                                               LCCollectionVec * inputCollectionVec,
                                                               LCCollectionVec * outputCollectionVec) const {

  //  message<DEBUG5> ( log() << "Applying cluster separation algorithm
  //  " << _separationAlgo );
  streamlog_out ( DEBUG0 ) <<  "Found "  << setVector.size() << " group(s) of merging clusters " << endl;

  if ( _separationAlgo == EUTELESCOPE::FLAGONLY ) {

    // with this splitting algorithm, the input and output collections
    // have the same number of entries. The only difference is that
    // there maybe some cluster in the output collection being flagged
    // as merged. So the first thing is to copy the content of the
    // input collection to the output one. The second step is to
    // modify the content of the output collection according to what
    // follow.
    for ( int iPulse = 0; iPulse < inputCollectionVec->getNumberOfElements() ; iPulse++ ) {
      TrackerPulseImpl * pulse    = dynamic_cast<TrackerPulseImpl *> ( inputCollectionVec->getElementAt( iPulse ) );
      TrackerPulseImpl * newPulse = new TrackerPulseImpl;
      newPulse->setCellID0( pulse->getCellID0() );
      newPulse->setCellID1( pulse->getCellID1() );
      newPulse->setTime   ( pulse->getTime() );
      newPulse->setCharge ( pulse->getCharge() );
      newPulse->setQuality( pulse->getQuality() );
      newPulse->setTrackerData( pulse->getTrackerData() );

      outputCollectionVec->push_back( newPulse );
    }


    CellIDDecoder<TrackerPulseImpl> cellDecoder(outputCollectionVec);


    int iCounter = 0;

    vector<set <int > >::iterator vectorIterator = setVector.begin();
    while ( vectorIterator != setVector.end() ) {


      streamlog_out ( DEBUG4 )  <<  "     Group " << (iCounter++) << " with the following clusters " << endl;

      set <int >::iterator setIterator = (*vectorIterator).begin();
      while ( setIterator != (*vectorIterator).end() ) {
        TrackerPulseImpl    * pulse   = dynamic_cast<TrackerPulseImpl * > ( outputCollectionVec->getElementAt( *setIterator ) ) ;
        ClusterType           type    = static_cast<ClusterType> (static_cast<int>( cellDecoder(pulse)["type"] ) );

        EUTelVirtualCluster * cluster;

        if ( type == kEUTelFFClusterImpl )
          cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );
        else if ( type == kEUTelBrickedClusterImpl )
          cluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );
        else
        {
          message<ERROR5> ( "Unknown cluster type. Sorry for quitting" ) ;
          throw UnknownDataTypeException("Cluster type unknown");
        }

        streamlog_out ( DEBUG4 ) << ( * cluster ) << endl;

        try {
	//TODO: FIxMe//cluster->setClusterQuality ( cluster->getClusterQuality() | kMergedCluster );
        } catch (lcio::ReadOnlyException& e) {
          streamlog_out ( WARNING2 )  << "Attempt to change the cluster quality on the original data\n"
            "This is possible only when the " << name() << " is not applied to data already on tape\n"
            "In this case only the pulse containing this cluster will have the proper quality" << endl;
        }
        pulse->setQuality ( static_cast<int> (ClusterQuality( pulse->getQuality() | kMergedCluster )) );
        ++setIterator;
        delete cluster;
      }
      ++vectorIterator;
    }
    return true;
  }

  return false;


}

void EUTelClusterSeparationProcessor::groupingMergingPairs(std::vector< std::pair<int , int> > pairVector,
                                                           std::vector< std::set<int > > * setVector) const {

  streamlog_out ( DEBUG0 ) << "Grouping merging pairs of clusters " << endl;

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

void EUTelClusterSeparationProcessor::check (LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelClusterSeparationProcessor::end() {
  streamlog_out ( MESSAGE2 ) << "Successfully finished" << endl;
}

