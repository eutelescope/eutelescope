// Version: $Id$
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $ $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
 
// since v00-00-09 this processor is built only if Marlin has GEAR
// support. In theory, since that version EUTelescope require GEAR,
// but it is better being sure
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelClusteringProcessor.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelSparseClusterImpl.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h" // this is because we want to use GEAR.

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCCollectionVec.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// system includes <>
#include <algorithm>
#include <string>
#include <vector>
#include <memory>
#include <list>
#include <cstdio>
#include <stdio.h>
#include <iostream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

static const int  MAXCLUSTERSIZE = 4096;


EUTelClusteringProcessor::EUTelClusteringProcessor () 
: Processor("EUTelClusteringProcessor"), 
  _nzsDataCollectionName(""),
  _zsDataCollectionName(""),
  _noiseCollectionName(""),
  _statusCollectionName(""),
  _pulseCollectionName(""),
  _hotPixelCollectionName(""),
  _initialPulseCollectionSize(0),
  _dummyCollectionName(""),
  _iRun(0),
  _nzsClusteringAlgo(""),
  _zsClusteringAlgo(""),
  _dataFormatType(""),
  _ffXClusterSize(0),
  _ffYClusterSize(0),
  _ffSeedCut(0.0),
  _sparseSeedCut(0.0),
  _ffClusterCut(0.0),
  _sparseClusterCut(0.0),
  _sparseMinDistanceSquared(2),
  _sparseMinDistance(0.0),
  _iEvt(0),
  _fillHistos(false),
  _histoInfoFileName(""),
  _seedCandidateMap(),
  _totClusterMap(),
  _noOfDetector(0),
  _ExcludedPlanes(),
  _clusterSpectraNVector(),
  _clusterSpectraNxNVector(),
  _clusterSignalHistos(),
  _clusterSizeXHistos(),
  _clusterSizeYHistos(),
  _seedSignalHistos(),
  _hitMapHistos(),
  _seedSNRHistos(),
  _clusterNoiseHistos(),
  _clusterSNRHistos(),
  _cluster_vs_seedSNRHistos(),
  _eventMultiplicityHistos(),
  _clusterSignal_NHistos(),
  _clusterSNR_NHistos(),
  _clusterSignal_NxNHistos(),
  _clusterSNR_NxNHistos(),
  _siPlanesParameters(NULL),
  _siPlanesLayerLayout(NULL),
  _isGeometryReady(false),
  _layerIndexMap(),
  _dutLayerIndexMap(),
  _ancillaryIndexMap(),
  _orderedSensorIDVec(),
  _sensorIDVec(),
  zsInputDataCollectionVec(NULL),
  nzsInputDataCollectionVec(NULL),
  pulseCollectionVec(NULL),
  noiseCollectionVec(NULL),
  statusCollectionVec(NULL),
  hotPixelCollectionVec(NULL),
  hasNZSData(false),
  hasZSData(false),
  _hitIndexMapVec()
 {
  
  // modify processor description
  _description =
    "EUTelClusteringProcessor is looking for clusters into a calibrated pixel matrix.";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "NZSDataCollectionName",
                           "Input calibrated data not zero suppressed collection name",
                           _nzsDataCollectionName, string ("data"));

  registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                           "Input of Zero Suppressed data",
                           _zsDataCollectionName, string ("zsdata") );

  registerInputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
                           "Noise (input) collection name",
                           _noiseCollectionName, string("noise"));

  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                           "Pixel status (input) collection name",
                           _statusCollectionName, string("status"));

  registerOutputCollection(LCIO::TRACKERPULSE, "PulseCollectionName",
                           "Cluster (output) collection name",
                           _pulseCollectionName, string("cluster"));

  // I believe it is safer not allowing the dummyCollection to be
  // renamed by the user. I prefer to set it once for ever here and
  // eventually, only if really needed, in the future allow add
  // another registerOutputCollection.
  _dummyCollectionName = "original_data";


  // now the optional parameters
  registerProcessorParameter ("ClusteringAlgo",
                              "Select here which algorithm should be used for clustering.\n"
                              "Available algorithms are:\n"
                              "-> FixedFrame: for custer with a given size\n"
                              "-> BrickedCluster: for bricked clustering on raw data",
                              _nzsClusteringAlgo, string(EUTELESCOPE::FIXEDFRAME));

  registerProcessorParameter ("ZSClusteringAlgo",
                              "Select here which algorithm should be used for clustering.\n"
                              "Available algorithms are:\n"
                              "-> SparseCluster: for cluster in ZS frame\n"
                              "-> SparseCluster2: for cluster in ZS frame with better performance\n"
                              "-> FixedFrame: for cluster with a given size\n"
                              "-> DFixedFrame: for digital cluster with a given size\n"
                              "-> BrickedCluster: for bricked clustering on zs data\n",
                              _zsClusteringAlgo, string(EUTELESCOPE::SPARSECLUSTER));

  registerProcessorParameter ("DataFormatType",
                              "Select herewith the type of the data format you are expecting from the sensors.\n"
                              "Available types of the data format:\n"
                              "-> Analog: smooth distribution of pixel ADC values from Min to Max\n"
                              "-> Digital: descrete distribution of pixel ADC values from Min to Max\n"
                              "-> Binary: only two values of the signal - 0 and 1\n",
                              _dataFormatType, string(EUTELESCOPE::BINARY));


  registerProcessorParameter ("FFClusterSizeX",
                              "Maximum allowed cluster size along x (only odd numbers)",
                              _ffXClusterSize, static_cast<int> (5));

  registerProcessorParameter ("FFClusterSizeY",
                              "Maximum allowed cluster size along y (only odd numbers)",
                              _ffYClusterSize, static_cast<int> (5));

  registerProcessorParameter ("FFSeedCut",
                              "Threshold in SNR for seed pixel identification",
                              _ffSeedCut, static_cast<float> (4.5));

  registerProcessorParameter ("FFClusterCut",
                              "Threshold in SNR for cluster identification",
                              _ffClusterCut, static_cast<float> (3.0));

  registerProcessorParameter("HistoInfoFileName", "This is the name of the histogram information file",
                             _histoInfoFileName, string( "histoinfo.xml" ) );

  registerProcessorParameter("SparseSeedCut","Threshold in SNR for seed pixel contained in ZS data",
                             _sparseSeedCut, static_cast<float > (4.5));

  registerProcessorParameter("SparseClusterCut","Threshold in SNR for clusters contained in ZS data",
                             _sparseClusterCut, static_cast<float > (3.0) );

  registerProcessorParameter("SparseMinDistanceSquared","Minimum distance squared between sparsified pixel ( touching == 2) ",
                             _sparseMinDistanceSquared, static_cast<int>(2) );
  
  registerProcessorParameter("SparseMinDistance","Minimum distance between sparsified pixel ( touching == sqrt(2)) ",
                             _sparseMinDistance, static_cast<float > (0.0 ) );

//  registerOptionalParameter("HotPixelDBFile","This is the name of the LCIO file name with the output hotpixel db (add .slcio)",
//                             _hotPixelDBFile, static_cast< string > ( "hotpixel.slcio" ) );

  registerOptionalParameter("HotPixelCollectionName","This is the name of the hotpixel collection",
                             _hotPixelCollectionName, static_cast< string > ( "hotpixel_m26" ) );


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  IntVec clusterNxNExample;
  clusterNxNExample.push_back(3);
  clusterNxNExample.push_back(5);

  registerOptionalParameter("ClusterNxN", "The list of cluster NxN to be filled."
                            "For example 3 means filling the 3x3 histogram spectrum",
                            _clusterSpectraNxNVector, clusterNxNExample);

  IntVec clusterNExample;
  clusterNExample.push_back(4);
  clusterNExample.push_back(9);
  clusterNExample.push_back(14);
  clusterNExample.push_back(19);
  clusterNExample.push_back(25);

  registerOptionalParameter("ClusterN", "The list of cluster N to be filled."
                            "For example 7 means filling the cluster spectra with the 7 most significant pixels",
                            _clusterSpectraNVector, clusterNExample );
#endif

  registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling",
                             _fillHistos, static_cast< bool > ( true ) );

  registerOptionalParameter("ExcludedPlanes", "The list of sensor ids that have to be excluded from the clustering.",
                             _ExcludedPlanes, std::vector<int> () );
  _isFirstEvent = true;
}


void EUTelClusteringProcessor::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to

  printParameters ();

  // in the case the FIXEDFRAME algorithm is selected, the check if
  // the _ffXClusterSize and the _ffYClusterSize are odd numbers
  if ( 
          _nzsClusteringAlgo == EUTELESCOPE::FIXEDFRAME  || _zsClusteringAlgo == EUTELESCOPE::FIXEDFRAME
          ||
          _nzsClusteringAlgo == EUTELESCOPE::DFIXEDFRAME  || _nzsClusteringAlgo == EUTELESCOPE::DFIXEDFRAME
       ) {
    bool isZero = ( _ffXClusterSize <= 0 );
    bool isEven = ( _ffXClusterSize % 2 == 0 );
    if ( isZero || isEven ) {
      throw InvalidParameterException("_ffXClusterSize has to be positive and odd");
    }
    isZero = ( _ffYClusterSize <= 0 );
    isEven = ( _ffYClusterSize % 2 == 0 );
    if ( isZero || isEven ) {
      throw InvalidParameterException("_ffYClusterSize has to be positive and odd");
    }
  }


  if (  
          _nzsClusteringAlgo == EUTELESCOPE::BRICKEDCLUSTER  
          ||
          _zsClusteringAlgo == EUTELESCOPE::BRICKEDCLUSTER
          )
  {
      if ( ! (_ffXClusterSize == 3 ) && (_ffYClusterSize == 3 ) )
      {
          streamlog_out ( ERROR2 ) << "[init()] For bricked pixel clustering the cluster size has to be 3x3 at the moment(!). Sorry!";
          throw InvalidParameterException("Set cluster size to 3x3 for bricked clustering!");
      }
  }

  // reset hotpixel map vectors
  _hitIndexMapVec.clear();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // the geometry is not yet initialized, so set the corresponding
  // switch to false
  _isGeometryReady = false;
}



void EUTelClusteringProcessor::processRunHeader (LCRunHeader * rdr) {


  auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );

  runHeader->addProcessor( type() );

  // increment the run counter
  ++_iRun;

}

void EUTelClusteringProcessor::initializeGeometry( LCEvent * event ) throw ( marlin::SkipEventException ) {
  // set the total number of detector to zero. This number can be
  // different from the one written in the gear description because
  // the input collection can contain only a fraction of all the
  // sensors.
  //
  // we assume that the no of detectors is the sum of the elements in
  // the NZS input collection and the in ZS one.
  _noOfDetector = 0;
  _sensorIDVec.clear();


  streamlog_out( DEBUG5 ) << "Initializing geometry" << endl;

  try {
    nzsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > (event->getCollection( _nzsDataCollectionName ) );
    _noOfDetector += nzsInputDataCollectionVec->getNumberOfElements();

    CellIDDecoder<TrackerDataImpl > cellDecoder( nzsInputDataCollectionVec );
    for ( size_t i = 0 ; i < nzsInputDataCollectionVec->size(); ++i ) {
      TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( nzsInputDataCollectionVec->getElementAt( i ) );
      _sensorIDVec.push_back( cellDecoder( data ) ["sensorID"] );
    }

  } catch ( lcio::DataNotAvailableException ) {
    // do nothing
    streamlog_out( DEBUG5 ) << "_nzsDataCollectionName " << _nzsDataCollectionName.c_str() << " not found " << endl; 
  }

  try 
  {
    zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _zsDataCollectionName ) ) ;
    _noOfDetector += zsInputDataCollectionVec->getNumberOfElements();

    CellIDDecoder<TrackerDataImpl > cellDecoder( zsInputDataCollectionVec );
    for ( size_t i = 0; i < zsInputDataCollectionVec->size(); ++i ) 
    {
      TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( i ) ) ;
      _sensorIDVec.push_back( cellDecoder( data )[ "sensorID" ] );
      _totClusterMap.insert( make_pair( cellDecoder( data )[ "sensorID" ] , 0 ));
    }

  } catch ( lcio::DataNotAvailableException ) {
    // do nothing again
    streamlog_out( DEBUG5 ) << "_zsDataCollectionName " << _zsDataCollectionName.c_str() << " not found " << endl; 
  }


  _siPlanesParameters  = const_cast< gear::SiPlanesParameters*  > ( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast< gear::SiPlanesLayerLayout* > ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  // now let's build a map relating the position in the layerindex
  // with the sensorID.
  _layerIndexMap.clear();
  for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); ++iLayer ) 
  {
      streamlog_out( DEBUG1) << 
          "Telescope: iLayer " << iLayer <<
          "["<< _siPlanesLayerLayout->getNLayers() << 
          "], getID:" <<  _siPlanesLayerLayout->getID( iLayer ) << endl;
      _layerIndexMap.insert( make_pair( _siPlanesLayerLayout->getID( iLayer ), iLayer ) );
  }

  // check if there is a DUT section or not
  _dutLayerIndexMap.clear();
  if( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT )    {

    // for the time being this is quite useless since, if there is a
    // DUT this is just one, but anyway it will become useful in a
    // short time.
    streamlog_out( DEBUG1) << "DUT:  getID: " << _siPlanesLayerLayout->getDUTID() << endl;
    _dutLayerIndexMap.insert( make_pair( _siPlanesLayerLayout->getDUTID(), 0 ) );
  }


  // now another map relating the position in the ancillary
  // collections (noise, pedestal and status) with the sensorID
  _ancillaryIndexMap.clear();
  _orderedSensorIDVec.clear();

  try {
    // this is the exemplary ancillary collection
    LCCollectionVec * noiseCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _noiseCollectionName ) );

    // prepare also a cell decoder
    CellIDDecoder< TrackerDataImpl > noiseDecoder( noiseCollectionVec );

    for ( size_t iDetector = 0 ; iDetector < noiseCollectionVec->size(); ++iDetector ) {
      TrackerDataImpl * noise = dynamic_cast< TrackerDataImpl * > ( noiseCollectionVec->getElementAt ( iDetector ) );
      _ancillaryIndexMap.insert( make_pair( noiseDecoder( noise ) ["sensorID"], iDetector ) );
      _orderedSensorIDVec.push_back( noiseDecoder( noise ) ["sensorID"] );
    }
  } catch (  lcio::DataNotAvailableException ) {
    streamlog_out( WARNING2 ) << "Unable to initialize the geometry. Trying with the following event" << endl;
    _isGeometryReady = false;
    throw SkipEventException( this ) ;
  }

  if ( _noOfDetector == 0 ) {
    streamlog_out( WARNING2 ) << "Unable to initialize the geometry. Trying with the following event" << endl;
    _isGeometryReady = false;
    throw SkipEventException( this );
  } else {
    _isGeometryReady = true;
  }
  
}

void EUTelClusteringProcessor::modifyEvent( LCEvent * /* event */ )
{
    return;
}


void EUTelClusteringProcessor::initializeHotPixelMapVec(  )
{
    streamlog_out(DEBUG4) << "initializeHotPixelMapVec, hotPixelCollectionVec size = " << hotPixelCollectionVec->size() << endl;
    
    // prepare some decoders
    CellIDDecoder<TrackerDataImpl> cellDecoder( hotPixelCollectionVec );
    CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );

    for ( unsigned int iDetector = 0 ; iDetector < hotPixelCollectionVec->size(); iDetector++ )         
    {
        if( _hitIndexMapVec.size() < iDetector+1 )            _hitIndexMapVec.resize(iDetector+1);

        TrackerDataImpl * hotData = dynamic_cast< TrackerDataImpl * > ( hotPixelCollectionVec->getElementAt( iDetector ) );
        int sensorID            = static_cast<int > ( cellDecoder( hotData )["sensorID"] );
 
        //if this is an excluded sensor go to the next element

        bool foundexcludedsensor = false;
        for(size_t j = 0; j < _ExcludedPlanes.size(); ++j)
        {
            if(_ExcludedPlanes[j] == sensorID)
            {
                foundexcludedsensor = true;
            }
        }
        if(foundexcludedsensor)  continue;

 
        if ( _layerIndexMap.find( sensorID ) == _layerIndexMap.end()   )
        {
            streamlog_out( DEBUG5 ) << "Sensor " << sensorID << " not found in the present data, skipping its hotPixel information." << endl;
            continue;
        }
       
        // the noise map. we only need this map for decoding issues.
        TrackerDataImpl    *noise  = 0;
 
        // the noise map. we only need this map for decoding issues.
        noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));

        // prepare the matrix decoder
        EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );

        // now prepare the EUTelescope interface to sparsified data.  
        auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > >  sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( hotData ));

        streamlog_out ( DEBUG1 ) << "Processing sparse data on detector " << sensorID << " with "
                                 << sparseData->size() << " pixels " << endl;
        
        for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ ) 
        {
            // loop over all pixels in the sparseData object.      
            EUTelGenericSparsePixel *sparsePixel =  new EUTelGenericSparsePixel() ;

            sparseData->getSparsePixelAt( iPixel, sparsePixel );
            int decoded_XY_index = matrixDecoder.getIndexFromXY( sparsePixel->getXCoord(), sparsePixel->getYCoord() ); // unique pixel index !!

            streamlog_out ( DEBUG1 )   << 
                " iPixel " << iPixel <<
                " idet " << iDetector <<
                " decoded_XY_index " << decoded_XY_index << endl;
            
            if( _hitIndexMapVec[iDetector].find( decoded_XY_index ) == _hitIndexMapVec[iDetector].end() )
            {
                _hitIndexMapVec[iDetector].insert ( make_pair ( decoded_XY_index, EUTELESCOPE::FIRINGPIXEL ) );               
                streamlog_out ( DEBUG4 ) 
                    << "adding hot pixel ["<< iPixel <<"]"
                    << " idet " << iDetector 
                    << " decoded_XY_index " << decoded_XY_index  
                    << " [" << sparsePixel->getXCoord() 
                    << " "<< sparsePixel->getYCoord() << "]" 
                    << " status : " << EUTELESCOPE::FIRINGPIXEL << endl;
              }
              else
              {
                  streamlog_out ( ERROR5 ) << "hot pixel [index " << decoded_XY_index << "] reoccured ?!" << endl;
              }
        }
  }
    
}
 
void EUTelClusteringProcessor::initializeStatusCollection(  )
{
    // prepare some decoders
    CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
    CellIDDecoder<TrackerDataImpl> statusDecoder( statusCollectionVec );
    CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );


    for ( unsigned int iDetector = 0 ; iDetector < zsInputDataCollectionVec->size(); iDetector++ ) 
    {
        if( _hitIndexMapVec.size() < iDetector+1 )            _hitIndexMapVec.resize(iDetector+1);
    
        // get the TrackerData and guess which kind of sparsified data it
        // contains.

        TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( iDetector ) );
        int sensorID            = static_cast<int > ( cellDecoder( zsData )["sensorID"] );


        //if this is an excluded sensor go to the next element

        bool foundexcludedsensor = false;
        for(size_t j = 0; j < _ExcludedPlanes.size(); ++j)
        {
            if(_ExcludedPlanes[j] == sensorID)
            {
                foundexcludedsensor = true;
            }
        }
        if(foundexcludedsensor)  continue;

        // get the noise and the status matrix with the right detectorID
        TrackerRawDataImpl * status = 0;

        // the noise map. we only need this map for decoding issues.
        TrackerDataImpl    * noise  = 0;
   
        // get the noise and the status matrix with the right detectorID
        status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));        
        vector< short > statusVec = status->adcValues();

        //the noise map. we only need this map for decoding issues.
        noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));


        // prepare the matrix decoder
        EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );

        // now prepare the EUTelescope interface to sparsified data.  
        auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > >  sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ));

        streamlog_out ( DEBUG1 ) << "Processing sparse data on detector " << sensorID << " with "
                                 << sparseData->size() << " pixels " << endl;
        
        for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ ) 
        {
            // loop over all pixels in the sparseData object.      
            EUTelGenericSparsePixel *sparsePixel =  new EUTelGenericSparsePixel() ;

            sparseData->getSparsePixelAt( iPixel, sparsePixel );
            int decoded_XY_index = matrixDecoder.getIndexFromXY( sparsePixel->getXCoord(), sparsePixel->getYCoord() ); // unique pixel index !!

            if( _hitIndexMapVec[iDetector].find( decoded_XY_index ) == _hitIndexMapVec[iDetector].end() )
            {
                int old_size = status->adcValues().size();
                // increment
              
                int new_size     = old_size + 1 ;
                int last_element = old_size;
              
                status->adcValues().resize( new_size  );
                _hitIndexMapVec[iDetector].insert ( make_pair ( decoded_XY_index, last_element ) );

              
                status->adcValues()[ last_element ] = EUTELESCOPE::HITPIXEL ;  // adcValues is a vector, there fore must address the elements incrementally
            }
            else
            {
                status->adcValues()[ _hitIndexMapVec[iDetector][ decoded_XY_index]  ] = EUTELESCOPE::HITPIXEL ;
            }
        }
    }

    return;
}


void EUTelClusteringProcessor::readCollections (LCEvent * event)
{

    try {
        nzsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > (event->getCollection( _nzsDataCollectionName ) );
        streamlog_out ( DEBUG4 ) << "nzsInputDataCollectionVec: " << _nzsDataCollectionName.c_str() << " found " << endl;
    } catch ( lcio::DataNotAvailableException ) {
        // do nothing
        streamlog_out ( DEBUG4 ) << "nzsInputDataCollectionVec: " << _nzsDataCollectionName.c_str() << " not found " << endl;
    }

    
    try {
        zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _zsDataCollectionName ) ) ;
        streamlog_out ( DEBUG4 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << endl;
    } catch ( lcio::DataNotAvailableException ) {
        // do nothing again
        streamlog_out ( DEBUG4 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " not found " << endl;
    }

   
    //
    //    
    statusCollectionVec = 0;
    try 
    {
        statusCollectionVec  = dynamic_cast < LCCollectionVec * > (event->getCollection( _statusCollectionName ));
        streamlog_out ( DEBUG4 ) << "statusCollectionName: " << _statusCollectionName.c_str() << " found " << endl;
    }
    catch (lcio::DataNotAvailableException& e ) 
    {
        streamlog_out ( DEBUG4 ) << "No status collection found in the event" << endl;
    }
 
    noiseCollectionVec = 0;
    try 
    {
        noiseCollectionVec  = dynamic_cast < LCCollectionVec * > (event->getCollection( _noiseCollectionName ));
        streamlog_out ( DEBUG4 ) << "noiseCollectionName: " << _noiseCollectionName.c_str() << " found " << endl;
    }
    catch (lcio::DataNotAvailableException& e ) 
    {
        streamlog_out ( DEBUG4 ) << "No noise pixel DB collection found in the event" << endl;
    }
 
    // the hotpixel db file should be read only once, and 
    // only after the statusCollection has been created (opened) 
    // 
    if( isFirstEvent() && statusCollectionVec != 0 )
    {
        hotPixelCollectionVec = 0;
        try 
        {
            hotPixelCollectionVec = static_cast< LCCollectionVec* >  (event->getCollection( _hotPixelCollectionName )) ;
            initializeHotPixelMapVec();
            streamlog_out ( DEBUG5 ) << "hotPixelCollectionName: " << _hotPixelCollectionName.c_str() << " found " << endl;
        } 
        catch (lcio::DataNotAvailableException& e ) 
        {
            streamlog_out ( WARNING5 ) << "No hot pixel DB collection found in the event" << endl;
        }

    }
 
    // in the current event it is possible to have either full frame and
    // zs data. Here is the right place to guess what we have

    //hasNZSData = true;
    hasZSData  = true;

    try 
    {
        event->getCollection(_nzsDataCollectionName);
    } 
    catch (lcio::DataNotAvailableException& e) 
    {
        hasNZSData = false;
        streamlog_out ( DEBUG4 ) << "No NZS data found in the event" << endl;
    }

    try 
    {
        event->getCollection( _zsDataCollectionName ) ;
    } 
    catch (lcio::DataNotAvailableException& e ) 
    {
        hasZSData = false;
        streamlog_out ( DEBUG4 ) << "No ZS data found in the event" << endl;
    }

    if ( !hasNZSData && !hasZSData ) 
    {
        streamlog_out ( MESSAGE2 ) << "The current event doesn't contain neither ZS nor NZS data collections: skip # " << event->getEventNumber() << endl;
//            << "Leaving this event without any further processing" << endl;
        throw SkipEventException( this ) ;
    }

    return;
}

void EUTelClusteringProcessor::processEvent (LCEvent * event) 
{

  ++_iEvt;

  // first of all we need to be sure that the geometry is properly
  // initialized!
  // 
  if ( !_isGeometryReady ) 
  {
    initializeGeometry( event ) ;
  }

  // 
  // read noise, status, and hotpixel collections
  // 
  readCollections(event);




#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  // book the histograms now
  if ( _fillHistos && isFirstEvent() ) 
  {
    bookHistos();
  }
#endif


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  if ( evt->getEventType() == kEORE ) 
  {
    streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  endl;
    return;
  }
  else if ( evt->getEventType() == kUNKNOWN ) 
  {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  // prepare a pulse collection to add all clusters found
  // this can be either a new collection or already existing in the
  // event  
  LCCollectionVec * pulseCollection;
  bool pulseCollectionExists = false;
  _initialPulseCollectionSize = 0;
  try 
  {
    pulseCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _pulseCollectionName ) );
    pulseCollectionExists = true;
    _initialPulseCollectionSize = pulseCollection->size();
  } 
  catch ( lcio::DataNotAvailableException& e ) 
  {
    pulseCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
  }

  // 
  // non Zero Suppresed (RAW data)
  // 
  if ( hasNZSData ) 
  {
    // put here all the possible algorithm applicable to NZS data
    if ( _nzsClusteringAlgo == EUTELESCOPE::FIXEDFRAME )     fixedFrameClustering(evt, pulseCollection);
    if ( _nzsClusteringAlgo == EUTELESCOPE::BRICKEDCLUSTER ) nzsBrickedClustering(evt, pulseCollection);   // force 3x3 clusters         
  }

  // 
  // ZS data
  //
  if ( hasZSData ) 
  {
    // put here all the possible algorithm applicable to ZS data
    if      ( _zsClusteringAlgo == EUTELESCOPE::SPARSECLUSTER  ) sparseClustering(evt, pulseCollection);
    else if ( _zsClusteringAlgo == EUTELESCOPE::SPARSECLUSTER2 ) sparseClustering(evt, pulseCollection);
    else if ( _zsClusteringAlgo == EUTELESCOPE::SPARSECLUSTER3 ) sparseClustering(evt, pulseCollection);
    else if ( _zsClusteringAlgo == EUTELESCOPE::FIXEDFRAME     ) zsFixedFrameClustering(evt, pulseCollection);
    else if ( _zsClusteringAlgo == EUTELESCOPE::DFIXEDFRAME    ) digitalFixedFrameClustering(evt, pulseCollection);

    // the bricked clustering type is solely needed for TAKI sensors type 
    else if ( _zsClusteringAlgo == EUTELESCOPE::BRICKEDCLUSTER ) zsBrickedClustering(evt, pulseCollection); // force 3x3 clusters      
  }

  // if the pulseCollection is not empty add it to the event
  if ( ! pulseCollectionExists && ( pulseCollection->size() != _initialPulseCollectionSize )) 
  {
    evt->addCollection( pulseCollection, _pulseCollectionName );
  }
  
  if ( pulseCollection->size() != _initialPulseCollectionSize ) 
  {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    if ( _fillHistos ) fillHistos(event);
#endif
  }
  if ( ! pulseCollectionExists && ( pulseCollection->size() == _initialPulseCollectionSize ) ) 
  {
    delete pulseCollection;
  }

  _isFirstEvent = false;
}





void EUTelClusteringProcessor::digitalFixedFrameClustering(LCEvent * evt, LCCollectionVec * pulseCollection) 
{
  streamlog_out ( DEBUG4 ) << "Looking for clusters in the zs data with digital FixedFrame algorithm " << endl;

  // get the collections of interest from the event.
//  LCCollectionVec * zsInputDataCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
//  LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection( _statusCollectionName ));
//  LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection( _noiseCollectionName ));

  // prepare some decoders
  CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
  CellIDDecoder<TrackerDataImpl> statusDecoder( statusCollectionVec );
  CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );

  // this is the equivalent of the dummyCollection in the fixed frame
  // clustering. BTW we should consider changing that "meaningful"
  // name! This contains cluster and not yet pulses
  bool isDummyAlreadyExisting = false;
  LCCollectionVec * sparseClusterCollectionVec = NULL;
  try {
    sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata") );
    isDummyAlreadyExisting = true ;
  } catch (lcio::DataNotAvailableException& e) {
    sparseClusterCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
    isDummyAlreadyExisting = false;
  }
  size_t dummyCollectionInitialSize = sparseClusterCollectionVec->size();



  //auto_ptr<LCCollectionVec > sparseClusterCollectionVec ( new  LCCollectionVec(LCIO::TRACKERDATA) );
  CellIDEncoder<TrackerDataImpl> idZSClusterEncoder( EUTELESCOPE::CLUSTERDEFAULTENCODING, sparseClusterCollectionVec  );
  // prepare an encoder also for the pulse collection
  CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
  // utility
  short limitExceed    = 0;

  if ( isFirstEvent() ) 
  {
    // For the time being nothing to do specifically in the first
    // event.
  }
  
  
  dim2array<bool> pixelmatrix(_ffXClusterSize, _ffYClusterSize, false);
  dim2array<int> pixelchargematrix(_ffXClusterSize, _ffYClusterSize, _charge);

  for ( unsigned int i = 0 ; i < zsInputDataCollectionVec->size(); i++ ) 
  {
    // get the TrackerData and guess which kind of sparsified data it
    // contains.
    TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( i ) );
    SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );

    int _sensorID            = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
    int sensorID            = _sensorID;
        
    //if this is an excluded sensor go to the next element
    bool foundexcludedsensor = false;
    for(size_t j = 0; j < _ExcludedPlanes.size(); ++j)
      {
        if(_ExcludedPlanes[j] == _sensorID)
          {
            foundexcludedsensor = true;
          }
      }
    if(foundexcludedsensor)
      continue;

    // reset the cluster counter for the clusterID
    int clusterID = 0;

    // get the noise and the status matrix with the right detectorID
//    TrackerRawDataImpl * status = 0;
    // the noise map. we only need this map for decoding issues.
    TrackerDataImpl    * noise  = 0;
   
    // get the noise and the status matrix with the right detectorID
//    status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));        
    //the noise map. we only need this map for decoding issues.
    noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));


//    if( _dataFormatType == EUTELESCOPE::BINARY )
//    {
//        if ( isFirstEvent() ) 
//        {
//            status->adcValues().clear();
//        }
//    }

    // reset the status

    // now that we know which is the sensorID, we can ask to GEAR
    // which are the minX, minY, maxX and maxY.
    int _minX, _minY, _maxX, _maxY;
    _minX = 0;
    _minY = 0;

    // this sensorID can be either a reference plane or a DUT, do it
    // differently...
    if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() ){
      // this is a reference plane
      _maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
      _maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
      
    } else if ( _dutLayerIndexMap.find( sensorID ) != _dutLayerIndexMap.end() ) {
      // ok it is a DUT plane
      _maxX = _siPlanesLayerLayout->getDUTSensitiveNpixelX() - 1;
      _maxY = _siPlanesLayerLayout->getDUTSensitiveNpixelY() - 1;

    } else {
      // this is not a reference plane neither a DUT... what's that?
//      throw  InvalidGeometryException ("Unknown sensorID " + to_string( sensorID ));
//      exit(-1);
        streamlog_out( ERROR5 ) << "Unknown sensorID " << sensorID << ", perhaps your GEAR file is incomplete." << endl;
        continue;
    }

    //todo: declare a vector of sensormatrizes as a class member in
    //order to not allocate a new object for each loop iteration.
    //sensormatrix.push_back(dim2array<bool>((unsigned int)(_maxX+1 - _minX), (unsigned int)(_maxY+1 - _minY), false));
    //    dim2array<bool> sensormatrix((unsigned int)(_maxX+1 - _minX), (unsigned int)(_maxY+1 - _minY), false);

    std::map<unsigned int, std::map<unsigned int, bool> > sensormatrix;
    std::map<unsigned int, std::map<unsigned int, int> > chargematrix;


    // prepare the matrix decoder
    EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );


    //insert some noise pixel by hand. was only used for debugging.
    // for(int i = 0; i < 250; i++)
    //       {
    //         int   index  = matrixDecoder.getIndexFromXY( i,i );
    //         status->adcValues()[index] = EUTELESCOPE::HITPIXEL;
    //       }

    // prepare a data vector mimicking the TrackerData data of the
    // standard digitalFixedFrameClustering. Initialize all the entries to zero.
    //    vector<float > dataVec( status->getADCValues().size(), 0. );
    vector<float > dataVec( noise->getChargeValues().size(), 0. );


    //seed candidates
    list<seed> seedcandidates;

    const int xoffset = _minX;
    const int yoffset = _minY;

    //    bool firstfoundhitpixel = true;

    if ( type == kEUTelGenericSparsePixel ) 
    {
      // now prepare the EUTelescope interface to sparsified data.
      auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( zsData ));

      streamlog_out ( DEBUG1 ) << "Processing sparse data on detector " << _sensorID << " with "
                               << sparseData->size() << " pixels " << endl;

      // loop over all pixels in the sparseData object.
      auto_ptr<EUTelGenericSparsePixel > sparsePixel( new EUTelGenericSparsePixel );
      for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ ) 
      {

          sparseData->getSparsePixelAt( iPixel, sparsePixel.get() );
          int index = matrixDecoder.getIndexFromXY( sparsePixel->getXCoord(), sparsePixel->getYCoord() );
	  float signal = sparsePixel->getSignal();
	  dataVec[ index  ] = signal;

	  streamlog_out ( DEBUG1) << "pixel x " << sparsePixel->getXCoord() << " y " <<  sparsePixel->getYCoord()
				  << " charge " << signal << endl;

          if(static_cast<int>(_hitIndexMapVec.size()) > sensorID ){
              if( _hitIndexMapVec[sensorID].find( index ) != _hitIndexMapVec[sensorID].end() )
          {
              streamlog_out ( DEBUG1) <<
                  " iDetector " << sensorID << 
                  " iPixel " << iPixel << 
                  " unique index " << index <<
                  " at x = " << sparsePixel->getXCoord() <<
                  " y= " << sparsePixel->getYCoord() << 
		" charge = " << signal << endl;
              continue;
          }
          }

                  sensormatrix[sparsePixel->getXCoord()][sparsePixel->getYCoord()] = true;
                  chargematrix[sparsePixel->getXCoord()][sparsePixel->getYCoord()] = signal;
      }
    }
    else 
    {
      throw UnknownDataTypeException("Unknown sparsified pixel");
    }


    
    // ------------------------------------------------------------------------------------------------------------------
    // now the seed pixel finding !!
    // 

    // RIA: 
    // this is neibhours counting only!
    // 
    const int stepx = static_cast<int>(_ffXClusterSize / 2);
    const int stepy = static_cast<int>(_ffYClusterSize / 2);
  
///    const int stepx = 1; wrong or OK ?? could skipp this counting for reeeally sparse data 
///    i.e. when one knows that there are mostly 1-2 pixels per clusters expected.
///    const int stepy = 1;
    
    
    std::map<unsigned int, std::map<unsigned int, bool> >::iterator pos;
    for(pos = sensormatrix.begin(); pos != sensormatrix.end(); ++pos) 
    {
        std::map<unsigned int, bool>::iterator sec;
        for (sec = sensormatrix[(*pos).first].begin(); sec != sensormatrix[(*pos).first].end(); ++sec) 
        {
            //              if(sensormatrix.at(i,j))
            if(sensormatrix[pos->first][sec->first]) 
            {
                const unsigned int i = pos->first;
                const unsigned int j = sec->first;

                // number of neighbours
                int nb = 0; 
 
                // total number of pixels in a cluster around the seed candidate 
                // (also diagonal elements are counted)
                int npixel_cl = 0; 

 
                // first npixel_cl will be determined
                for(unsigned int index_x = i-stepx; index_x <= (i + stepx); index_x++)
                {
                    for(unsigned int index_y = j-stepy; index_y <= (j + stepy);index_y++)
                    {

                        if(index_x > 0 && index_y > 0)
                        {
                            std::map<unsigned int, std::map<unsigned int, bool> >::const_iterator z = sensormatrix.find(index_x);
                            if(z!=sensormatrix.end() && z->second.find(index_y)!=z->second.end())
                            {
                                if(sensormatrix[index_x][index_y])
                                {
                                    npixel_cl++;                                    
                                }
                            }
                        }
                    }
                }
 
                if(npixel_cl > 1)
                {
                    if(i>=1)
		      for(int index_x = static_cast<int>(i-1); index_x <= static_cast<int>(i + 1); index_x++)
                    {
                        if(index_x >= 0)
                        {
                            std::map<unsigned int, std::map<unsigned int, bool> >::const_iterator z = sensormatrix.find(index_x);
                            if(z!=sensormatrix.end() && z->second.find(j)!=z->second.end())
                            {
                                if(sensormatrix[index_x][j])
                                  nb++;
                            }
                        }
                    }
                    
                    if(j>=1)
		      for(int index_y = static_cast<int>(j-1); index_y <= static_cast<int>(j + 1); index_y++)
                    {
                        if(index_y >= 0)
                        {
                            std::map<unsigned int, std::map<unsigned int, bool> >::const_iterator z = sensormatrix.find(i);
                            if(z!=sensormatrix.end() && z->second.find(index_y)!=z->second.end())
                            {
                                if(sensormatrix[i][index_y])
                                  nb++;
                            }
                        }
                    }
                }         // could all this passage be skipped ?

                //fill this pixel into the list of found seed pixel candidates
                seedcandidates.push_back(seed(i,j,nb,npixel_cl));
              }
          }
      } 
    // sort the list of seed pixel candidates. the first criteria is
    // the number of neighbours without diagonal neighbours. then the
    // second criteria is the total number of neighbours
 
    // sorts the seed list according to the "operator<" definition in the seed class
    seedcandidates.sort();       

    //end of seed pixel finding!

    //if at least one seed pixel candidate was found, then ...
    if(!seedcandidates.empty())
    {
        list<seed>::iterator i;
        //loop over all found seed pixel candidates
        for( i = seedcandidates.begin(); i != seedcandidates.end(); ++i)
        {
            //check that this pixel was not used before.
            
            std::map<unsigned int, std::map<unsigned int, bool> >::const_iterator z = sensormatrix.find( i->x );
            if(z!=sensormatrix.end() && z->second.find( i->y )!=z->second.end())
            { 
                if(sensormatrix[i->x][i->y])
                {
                    std::vector<pixel> pix;
                    // select pixels around the seed pixel

                    if(i->x >= static_cast<unsigned int>(stepx) && i->y >= static_cast<unsigned int>(stepy))
                    {
		      for(int index_x = static_cast<int>(i->x - stepx); index_x <= static_cast<int>(i->x + stepx); index_x++)
                        {
			  for(int index_y = static_cast<int>(i->y - stepy); index_y <= static_cast<int>(i->y + stepy);index_y++)
                            {
                                if(index_x >= 0 && index_y >= 0)
                                {
                                    std::map<unsigned int, std::map<unsigned int, bool> >::const_iterator z = sensormatrix.find( index_x );
                                    if(z!=sensormatrix.end() && z->second.find( index_y )!=z->second.end())
                                    {
                                        if(sensormatrix[index_x][index_y])
                                        {
                                            pix.push_back(pixel(index_x, index_y));
                                        }
                                    }
                                }
                            }
                        }
                    }
                          
                    // pix is a vector with all found "good" pixel, that
                    // were not used before in a different cluster.

                    // cut on the number of pixel. dont
                    // apply this cut here, use it in the
                    // filtering processor?

                    if(!pix.empty()) 
                    {
                        // we found a cluster ...

                        IntVec   clusterCandidateIndeces;
                        FloatVec clusterCandidateCharges;
                        ClusterQuality cluQuality = kGoodCluster;

                        // the pixel coordinates of the seed pixels are
                        // needed later
                        int seedX = -1;
                        int seedY = -1;

                        // reset the pixel matrix
                        // a matrix of pixel for this cluster. it is needed
                        // for decoding issues.
                        pixelmatrix.pad(false);

                        // loop over all hit pixels inside this cluster
                        for(unsigned int j = 0; j < pix.size(); j++)
                        {
                            // remove pixels, that were assigned to this
                            // cluster from the dummy sensor map. this
                            // pixel will then not be used then in other clusters
                            // sensormatrix.set(pix[j].x,pix[j].y,false);

                            std::map<unsigned int, std::map<unsigned int, bool> >::const_iterator z = sensormatrix.find( pix[j].x );
                            if(z!=sensormatrix.end() && z->second.find( pix[j].y )!=z->second.end())
                            {
                                sensormatrix[pix[j].x][pix[j].y] = false;
                            }   

                            // dont forget to apply the offset correction!
//                            int index = matrixDecoder.getIndexFromXY(pix[j].x + xoffset, pix[j].y + yoffset);

                            if(pix[j].x == i->x  && pix[j].y == i->y)
                            {
                                // this is the seed pixel!
                                seedX = pix[j].x + xoffset;
                                seedY = pix[j].y + yoffset;
                            }
                            else
                            {
                                // this is a neighbour pixel!
                                // nothing to do?
                            }

                                clusterCandidateIndeces.push_back(-1);
                                cluQuality = cluQuality | kIncompleteCluster | kMergedCluster ;
                        }
                
                        // sanity check
                        if(seedX == -1 || seedY == -1)
                        {
			  streamlog_out(DEBUG5) << "a cluster was found but no seed pixel coordinates!" << endl;
			  streamlog_out(DEBUG5) << pix.size() << " " << i->x << " " << i->y << endl;
			  exit(-1);
                        }

                        // now lets fill the cluster pixel matrix, which is required
                        // by the decoding of the cluster into a 1d array (clusterCandidateCharges).
                        for(unsigned int j = 0; j < pix.size(); j++)
                        {
                            // set the hits. all other pixels are by
                            // default false. the seed pixel is in the
                            // center of this matrix.

                            pixelmatrix.set(
                                    pix[j].x + xoffset - seedX + static_cast<int>(_ffXClusterSize / 2),
                                    pix[j].y + yoffset - seedY + static_cast<int>(_ffYClusterSize / 2),
                                    true
                                    );
                            pixelchargematrix.set(
                                    pix[j].x + xoffset - seedX + static_cast<int>(_ffXClusterSize / 2),
                                    pix[j].y + yoffset - seedY + static_cast<int>(_ffYClusterSize / 2),				    
				    chargematrix[pix[j].x][pix[j].y]
				    );
                        }

                        // loop over the cluster pixels and fill them into
                        // the 1d array. The ordering of the two loops is
                        // copied from the CoG shift method of the class EUTelDFFClusterImpl
			int xsize=0;
			int ysize=0;
                           
                        for(int yPixel = 0; yPixel < _ffYClusterSize; yPixel++)
                        {
                            for(int xPixel = 0; xPixel < _ffXClusterSize; xPixel++)
                            {
                                if(pixelmatrix.at(xPixel,yPixel))
                                {
				  //clusterCandidateCharges.push_back(1.0);
				  clusterCandidateCharges.push_back(pixelchargematrix.at(xPixel,yPixel));				  
                                }
                                else
                                {
                                    clusterCandidateCharges.push_back(0.0);
                                }
                            }
                        }

                        for(int yPixel = 0; yPixel < _ffYClusterSize; yPixel++){
			  if( pixelmatrix.at(0,yPixel) || pixelmatrix.at(1,yPixel) || pixelmatrix.at(2,yPixel) ) ysize ++;
			}
			for(int xPixel = 0; xPixel < _ffXClusterSize; xPixel++){
			  if( pixelmatrix.at(xPixel,0) || pixelmatrix.at(xPixel,1) || pixelmatrix.at(xPixel,2) ) xsize ++;
			}


                        // check whether this cluster is partly outside
                        // the sensor matrix
                        if(
                                   (seedX - stepx ) < _minX
                                || (seedX + stepx ) > _maxX
                                || (seedY - stepy ) < _minY
                                || (seedY + stepy ) > _maxY
                                )
                        {
                            cluQuality = cluQuality | kBorderCluster;
                        }
     
                        // the final cluster creation

                        // the final result of the clustering will enter in a
                        // TrackerPulseImpl in order to be algorithm independent

                        TrackerPulseImpl * pulse = new TrackerPulseImpl;
                        CellIDEncoder<TrackerPulseImpl> idPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
                        idPulseEncoder["sensorID"]      = _sensorID;
                        idPulseEncoder["xSeed"]         = seedX;
                        idPulseEncoder["ySeed"]         = seedY;
                        idPulseEncoder["xCluSize"]      = _ffXClusterSize; //xsize; //_ffXClusterSize;
                        idPulseEncoder["yCluSize"]      = _ffYClusterSize; //ysize; //_ffYClusterSize;
                        idPulseEncoder["type"]          = static_cast<int>(kEUTelDFFClusterImpl);
                        idPulseEncoder.setCellID(pulse);

                       TrackerDataImpl * cluster = new TrackerDataImpl;
                        CellIDEncoder<TrackerDataImpl> idClusterEncoder(EUTELESCOPE::CLUSTERDEFAULTENCODING, sparseClusterCollectionVec);
                        idClusterEncoder["sensorID"]      = _sensorID;
                        idClusterEncoder["xSeed"]         = seedX;
                        idClusterEncoder["ySeed"]         = seedY;
                        idClusterEncoder["xCluSize"]      = _ffXClusterSize; //xsize; //_ffXClusterSize;
                        idClusterEncoder["yCluSize"]      = _ffYClusterSize; //ysize; //_ffYClusterSize;
                        idClusterEncoder["quality"]       = static_cast<int>(cluQuality);
                        idClusterEncoder.setCellID(cluster);

                        streamlog_out (DEBUG0) << "  Cluster no " <<  clusterID << " seedX " << seedX << " seedY " << seedY << endl;
/*
                        IntVec::iterator indexIter = clusterCandidateIndeces.begin();
                        while ( indexIter != clusterCandidateIndeces.end() ) 
                        {
                            if((*indexIter) != -1)
                            {
                                if( _dataFormatType == EUTELESCOPE::BINARY )
                                { 
                                    status->adcValues()[ _indexMap[(*indexIter)] ] = EUTELESCOPE::HITPIXEL;
                                }else{
                                    status->adcValues()[(*indexIter)] = EUTELESCOPE::HITPIXEL;
                                }
                            }
                            ++indexIter;
                        }
*/

                        // copy the candidate charges inside the cluster
                        cluster->setChargeValues(clusterCandidateCharges);
                        sparseClusterCollectionVec->push_back(cluster);


//continue;

                        EUTelDFFClusterImpl * eutelCluster = new EUTelDFFClusterImpl( cluster );
                        pulse->setCharge(eutelCluster->getTotalCharge());

                        delete eutelCluster;

                        pulse->setQuality(static_cast<int>(cluQuality));
                        pulse->setTrackerData(cluster);
                        pulseCollection->push_back(pulse);




                        // increment the cluster counters
                        _totClusterMap[ sensorID ] += 1;
                        ++clusterID;
                        if ( clusterID >= MAXCLUSTERSIZE ) {
                            ++limitExceed;
                            --clusterID;
                            streamlog_out ( WARNING2 ) << "Event " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                << " on detector " << _sensorID
                                << " contains more than " << MAXCLUSTERSIZE << " cluster (" << clusterID + limitExceed << ")" << endl;
                            
                        } 
                    }
                }
            }
        } //loop over all found seed pixel candidates :: END
    } // LOOP over all seedcandidates :: END
  } // for ( unsigned int i = 0 ; i < zsInputDataCollectionVec->size(); i++ ) :: END
  

  // if the sparseClusterCollectionVec isn't empty add it to the
  // current event. The pulse collection will be added afterwards

  if ( ! isDummyAlreadyExisting ) 
  {
      if ( sparseClusterCollectionVec->size() != dummyCollectionInitialSize ) 
      {
          evt->addCollection( sparseClusterCollectionVec, "original_zsdata" );
      } 
      else 
      {
          delete sparseClusterCollectionVec;
      }
  }

}



void EUTelClusteringProcessor::zsFixedFrameClustering(LCEvent * evt, LCCollectionVec * pulseCollection) {

  streamlog_out ( DEBUG4 ) << "Looking for clusters in the zs data with FixedFrame algorithm " << endl;

  // get the collections of interest from the event.
//  LCCollectionVec * zsInputDataCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
//  LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection( _noiseCollectionName ));
//  LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection( _statusCollectionName ));
  // prepare some decoders
  CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
  CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );

  // this is the equivalent of the dummyCollection in the fixed frame
  // clustering. BTW we should consider changing that "meaningful"
  // name! This contains cluster and not yet pulses
  bool isDummyAlreadyExisting = false;
  LCCollectionVec * sparseClusterCollectionVec = NULL;
  try {
    sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata") );
    isDummyAlreadyExisting = true ;
  } catch (lcio::DataNotAvailableException& e) {
    sparseClusterCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
    isDummyAlreadyExisting = false;
  }
  size_t dummyCollectionInitialSize = sparseClusterCollectionVec->size();

  // prepare an encoder also for the dummy collection.
  // even if it is a sparse data set, since the clustering is FF apply
  // the standard CLUSTERDEFAULTENCODING
  CellIDEncoder<TrackerDataImpl> idClusterEncoder( EUTELESCOPE::CLUSTERDEFAULTENCODING, sparseClusterCollectionVec  );

  // prepare an encoder also for the pulse collection
  CellIDEncoder<TrackerPulseImpl> idPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);

  // utility
  short limitExceed    = 0;

  if ( isFirstEvent() ) {

    // For the time being nothing to do specifically in the first
    // event.

  }

  for ( unsigned int i = 0 ; i < zsInputDataCollectionVec->size(); i++ ) {
    // get the TrackerData and guess which kind of sparsified data it
    // contains.
    TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( i ) );
    SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );

    int sensorID             = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
    //if this is an excluded sensor go to the next element
    bool foundexcludedsensor = false;
    for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
      {
        if(_ExcludedPlanes[i] == sensorID)
          {
            foundexcludedsensor = true;
          }
      }
    if(foundexcludedsensor)
      continue;
    // now that we know which is the sensorID, we can ask to GEAR
    // which are the minX, minY, maxX and maxY.
    int minX, minY, maxX, maxY;
    minX = 0;
    minY = 0;

    // this sensorID can be either a reference plane or a DUT, do it
    // differently...
    if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() ){
      // this is a reference plane
      maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
      maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
    } else if ( _dutLayerIndexMap.find( sensorID ) != _dutLayerIndexMap.end() ) {
      // ok it is a DUT plane
      maxX = _siPlanesLayerLayout->getDUTSensitiveNpixelX() - 1;
      maxY = _siPlanesLayerLayout->getDUTSensitiveNpixelY() - 1;
    } else {
      // this is not a reference plane neither a DUT... what's that?
//      throw  InvalidGeometryException ("Unknown sensorID " + to_string( sensorID ));
      streamlog_out( ERROR5 ) << "Unknown sensorID " << sensorID << ", perhaps your GEAR file is incomplete." << endl;
      continue;
    }

    // reset the cluster counter for the clusterID
    int clusterID = 0;

    // get the noise and the status matrix with the right detectorID
    TrackerDataImpl    * noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));
    TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));


    // reset the status
    resetStatus(status);

    // prepare the matrix decoder
    EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );

    // prepare a data vector mimicking the TrackerData data of the
    // standard FixedFrameClustering. Initialize all the entries to zero.
    vector<float > dataVec( noise->getChargeValues().size(), 0. );

    // prepare a multimap for the seed candidates
    multimap<float , int > seedCandidateMap;

    if ( type == kEUTelGenericSparsePixel ) {

      // now prepare the EUTelescope interface to sparsified data.
      auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( zsData ));

      streamlog_out ( DEBUG1 ) << "Processing sparse data on detector " << sensorID << " with "
                               << sparseData->size() << " pixels " << endl;

      // loop over all pixels in the sparseData object.
      auto_ptr<EUTelGenericSparsePixel > sparsePixel( new EUTelGenericSparsePixel );
      for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ ) {
        sparseData->getSparsePixelAt( iPixel, sparsePixel.get() );
        int   index  = matrixDecoder.getIndexFromXY( sparsePixel->getXCoord(), sparsePixel->getYCoord() );
        float signal = sparsePixel->getSignal();
        dataVec[ index  ] = signal;
        if( static_cast<int>(status->getADCValues().size()) < index )
        {
            status->adcValues().resize(index+1);
        }
        if (  ( signal  > _ffSeedCut * noise->getChargeValues()[ index ] ) &&
              ( status->getADCValues()[ index ] == EUTELESCOPE::GOODPIXEL ) ) {
          seedCandidateMap.insert ( make_pair ( signal, index ) );
          streamlog_out ( DEBUG1 ) << "Added pixel " << sparsePixel->getXCoord()
                                   << ", " << sparsePixel->getYCoord()
                                   << " with signal " << signal
                                   << " to the seedCandidateMap" << endl;
        }

      }
    } else {
      throw UnknownDataTypeException("Unknown sparsified pixel");
    }

    if ( !seedCandidateMap.empty() ) {

      streamlog_out ( DEBUG0 ) << "There are " << seedCandidateMap.size() << " seed candidates." << endl;

      // now build up a cluster for each seed candidate
      multimap<float, int >::reverse_iterator rMapIter = seedCandidateMap.rbegin();
      while ( rMapIter != seedCandidateMap.rend() ) {
      
        // Remove hot pixel:
        if( _hitIndexMapVec.size() > static_cast<unsigned int>(sensorID)) {
              if( _hitIndexMapVec[sensorID].find( (*rMapIter).second ) != _hitIndexMapVec[sensorID].end() )
              {
                int seedX, seedY;
                matrixDecoder.getXYFromIndex ( (*rMapIter).second, seedX, seedY );
                streamlog_out ( DEBUG5 ) << "Detector " << sensorID << " Pixel " << seedX << " " << seedY << " -- HOTPIXEL, skipping... " << endl;
                ++rMapIter;
                continue;
              }
        }
        if ( status->adcValues()[ (*rMapIter).second ] == EUTELESCOPE::GOODPIXEL ) {
          // if we enter here, this means that at least the seed pixel
          // wasn't added yet to another cluster.  Note that now we need
          // to build a candidate cluster that has to pass the
          // clusterCut to be considered a good cluster
          double clusterCandidateSignal    = 0.;
          double clusterCandidateNoise2    = 0.;
          FloatVec clusterCandidateCharges;
          IntVec   clusterCandidateIndeces;
          int seedX, seedY;
          matrixDecoder.getXYFromIndex ( (*rMapIter).second, seedX, seedY );

          // start looping around the seed pixel. Remember that the seed
          // pixel has to stay in the center of cluster
          ClusterQuality cluQuality = kGoodCluster;
          for (int yPixel = seedY - (_ffYClusterSize / 2); yPixel <= seedY + (_ffYClusterSize / 2); yPixel++) {
            for (int xPixel =  seedX - (_ffXClusterSize / 2); xPixel <= seedX + (_ffXClusterSize / 2); xPixel++) {
              // always check we are still within the sensor!!!
              if ( ( xPixel >= minX )  &&  ( xPixel <= maxX ) &&
                   ( yPixel >= minY )  &&  ( yPixel <= maxY ) ) {
                int index = matrixDecoder.getIndexFromXY(xPixel, yPixel);
                
                bool isHit  = ( status->getADCValues()[index] == EUTELESCOPE::HITPIXEL  );
                bool isGood = ( status->getADCValues()[index] == EUTELESCOPE::GOODPIXEL );

                if(isGood)
                  clusterCandidateIndeces.push_back(index);
                else
                  clusterCandidateIndeces.push_back(-1);

                if ( isGood && !isHit ) {
                  // if the pixel wasn't selected, then its signal
                  // will be 0.0. Mark it in the status
                  if ( dataVec[ index ] == 0.0 )
                    status->adcValues()[ index ] = EUTELESCOPE::MISSINGPIXEL ;
                  clusterCandidateSignal += dataVec[ index ] ;
                  clusterCandidateNoise2 += pow ( noise->getChargeValues() [ index ], 2 );
                  clusterCandidateCharges.push_back( dataVec[ index ] );
                } else if ( isHit ) {
                  // this can be a good place to flag the current
                  // cluster as kMergedCluster, but it would introduce
                  // a bias since the at least another cluster (the
                  // one which this pixel belong to) is not flagged.
                  //
                  // In order to flag all merged clusters and possibly
                  // try to separate the different contributions use
                  // the EUTelSeparateClusterProcessor. In this
                  // processor not all the merged clusters will be
                  // flagged as kMergedCluster | kIncompleteCluster
                  cluQuality = cluQuality | kIncompleteCluster | kMergedCluster ;
                  clusterCandidateCharges.push_back(0.);
                } else if ( !isGood ) {
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
          if ( clusterCandidateSignal > _ffClusterCut * sqrt( clusterCandidateNoise2 ) ) {
            // the cluster candidate is a good cluster
            // mark all pixels belonging to the cluster as hit
            IntVec::iterator indexIter = clusterCandidateIndeces.begin();

            // the final result of the clustering will enter in a
            // TrackerPulseImpl in order to be algorithm independent
            TrackerPulseImpl * pulse = new TrackerPulseImpl;
            idPulseEncoder["sensorID"]      = sensorID;
            idPulseEncoder["xSeed"]         = seedX;
            idPulseEncoder["ySeed"]         = seedY;
            idPulseEncoder["xCluSize"]      = _ffXClusterSize;
            idPulseEncoder["yCluSize"]      = _ffYClusterSize;
            idPulseEncoder["type"]          = static_cast<int>(kEUTelFFClusterImpl);
            idPulseEncoder.setCellID(pulse);

            TrackerDataImpl * cluster = new TrackerDataImpl;
            idClusterEncoder["sensorID"]      = sensorID;
            idClusterEncoder["xSeed"]         = seedX;
            idClusterEncoder["ySeed"]         = seedY;
            idClusterEncoder["xCluSize"]      = _ffXClusterSize;
            idClusterEncoder["yCluSize"]      = _ffYClusterSize;
            idClusterEncoder["quality"]       = static_cast<int>(cluQuality);
            idClusterEncoder.setCellID(cluster);


            streamlog_out (DEBUG0) << "  Cluster no " <<  clusterID << " seedX " << seedX << " seedY " << seedY << endl;


            while ( indexIter != clusterCandidateIndeces.end() ) {
              if((*indexIter) != -1)
                status->adcValues()[(*indexIter)] = EUTELESCOPE::HITPIXEL;
              ++indexIter;
            }


            // copy the candidate charges inside the cluster
            cluster->setChargeValues(clusterCandidateCharges);
            sparseClusterCollectionVec->push_back(cluster);

            EUTelFFClusterImpl * eutelCluster = new EUTelFFClusterImpl( cluster );
            pulse->setCharge(eutelCluster->getTotalCharge());
            delete eutelCluster;

            pulse->setQuality(static_cast<int>(cluQuality));
            pulse->setTrackerData(cluster);
            pulseCollection->push_back(pulse);

            // increment the cluster counters
            _totClusterMap[ sensorID ] += 1;
            ++clusterID;
            if ( clusterID >= MAXCLUSTERSIZE ) {
              ++limitExceed;
              --clusterID;
              streamlog_out ( WARNING2 ) << "Event " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                         << " on detector " << sensorID
                                         << " contains more than " << MAXCLUSTERSIZE << " cluster (" << clusterID + limitExceed << ")" << endl;
            }
          }
        }
        ++rMapIter;
      }
    }


  }

  // if the sparseClusterCollectionVec isn't empty add it to the
  // current event. The pulse collection will be added afterwards
  if ( ! isDummyAlreadyExisting ) {
    if ( sparseClusterCollectionVec->size() != dummyCollectionInitialSize ) {
      evt->addCollection( sparseClusterCollectionVec, "original_zsdata" );
    } else {
      delete sparseClusterCollectionVec;
    }
  }

}


void EUTelClusteringProcessor::zsBrickedClustering(LCEvent * evt, LCCollectionVec * pulseCollection) {

  streamlog_out ( DEBUG4 ) << "Looking for clusters in the zs data with zsBrickedClustering algorithm " << endl;

  // get the collections of interest from the event.
//  LCCollectionVec * zsInputDataCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
//  LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
//  LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
  
  // prepare some decoders
  CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
  CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );

  // this is the equivalent of the dummyCollection in the fixed frame
  // clustering. BTW we should consider changing that "meaningful"
  // name! This contains cluster and not yet pulses
  bool isDummyAlreadyExisting = false;
  LCCollectionVec * sparseClusterCollectionVec = NULL;
  try
    {
      sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata") );
      isDummyAlreadyExisting = true ;
    }
  catch (lcio::DataNotAvailableException& e)
    {
      sparseClusterCollectionVec =  new LCCollectionVec(LCIO::TRACKERDATA);
      isDummyAlreadyExisting = false;
    }
  CellIDEncoder<TrackerDataImpl> idZSClusterEncoder( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec  );

  // prepare an encoder also for the pulse collection
  CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);

  // utility
  short limitExceed    = 0;

  if ( isFirstEvent() )
    {

      // For the time being nothing to do specifically in the first
      // event.

    }

  for ( unsigned int i = 0 ; i < zsInputDataCollectionVec->size(); i++ )
    {
      // get the TrackerData and guess which kind of sparsified data it
      // contains.
      TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( i ) );
      SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
      int sensorID             = static_cast<int > ( cellDecoder( zsData )["sensorID"] );

      // now that we know which is the sensorID, we can ask to GEAR
      // which are the minX, minY, maxX and maxY.
      int minX, minY, maxX, maxY;
      minX = 0;
      minY = 0;

      // this sensorID can be either a reference plane or a DUT, do it
      // differently...
      if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() )
        {
          // this is a reference plane
          maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
          maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
        }
      else if ( _dutLayerIndexMap.find( sensorID ) != _dutLayerIndexMap.end() )
        {
          // ok it is a DUT plane
          maxX = _siPlanesLayerLayout->getDUTSensitiveNpixelX() - 1;
          maxY = _siPlanesLayerLayout->getDUTSensitiveNpixelY() - 1;
        }
      else
        {
          // this is not a reference plane neither a DUT... what's that?
//          throw  InvalidGeometryException ("Unknown sensorID " + to_string( sensorID ));
          streamlog_out( ERROR5 ) << "Unknown sensorID " << sensorID << ", perhaps your GEAR file is incomplete." << endl;
          continue;
        }

      // reset the cluster counter for the clusterID
      int clusterID = 0;

      // get the noise and the status matrix with the right detectorID
      TrackerDataImpl    * noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));
      TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));

      // reset the status
      resetStatus(status);

      // prepare the matrix decoder
      EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );

      // prepare a data vector mimicking the TrackerData data of the
      // standard FixedFrameClustering. Initialize all the entries to zero.
      // NOTE
      // TAKI 0.0001 instead of 0.0, because we have integers coming in from the DUT.
      // And these might very well be 0 -> 0.0 quite often! So 0.0001 is used here.
      // If the 0.0001 value is found here later on again, then we know that the corresponding pixel was not transmitted!
      vector<float > dataVec( status->getADCValues().size(), 0.0001 );

      // prepare a multimap for the seed candidates
      multimap<float , int > seedCandidateMap;

      if ( type == kEUTelGenericSparsePixel )
        {

          // now prepare the EUTelescope interface to sparsified data.
	  auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > sparseData(new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( zsData ));

          streamlog_out ( DEBUG1 ) << "Processing sparse data on detector " << sensorID << " with "
                                   << sparseData->size() << " pixels " << endl;

          // loop over all pixels in the sparseData object.
          auto_ptr<EUTelGenericSparsePixel > sparsePixel( new EUTelGenericSparsePixel );
          for ( unsigned int iPixel = 0; iPixel < sparseData->size(); iPixel++ )
            {
              sparseData->getSparsePixelAt( iPixel, sparsePixel.get() );
              int   index  = matrixDecoder.getIndexFromXY( sparsePixel->getXCoord(), sparsePixel->getYCoord() );
              float signal = sparsePixel->getSignal();
              dataVec[ index ] = signal;

              //! CUT 1
              if (  ( signal  > _ffSeedCut * noise->getChargeValues()[ index ] ) &&
                    ( status->getADCValues()[ index ] == EUTELESCOPE::GOODPIXEL ) )
                {
                  seedCandidateMap.insert ( make_pair ( signal, index ) );
                  streamlog_out ( DEBUG1 ) << "Added pixel " << sparsePixel->getXCoord()
                                           << ", " << sparsePixel->getYCoord()
                                           << " with signal " << signal
                                           << " to the seedCandidateMap" << endl;

                  if ( noise->getChargeValues()[ index ] < 0.01 )
                    {
                      streamlog_out ( ERROR2 ) << "ZERO NOISE SEED PIXEL ADDED!"
                                               << "\n x=" << sparsePixel->getXCoord()
                                               << "\n y=" << sparsePixel->getYCoord()
                                               << "\n amp=" << signal
                                               << "\n status=" << status->getADCValues()[ index ]
                                               <<    " GOODP   =  0,"
                                               <<    " BAD     =  1,"
                                               <<    " HIT     = -1,"
                                               <<    " MISSING =  2,"
                                               <<    " FIRING  =  3.";
                    }
                }

            }
        }
      else
        {
          throw UnknownDataTypeException("Unknown sparsified pixel");
        }

      if ( !seedCandidateMap.empty() )
        {

          streamlog_out ( DEBUG0 ) << "  Seed candidates " << seedCandidateMap.size() << endl;

          // now build up a cluster for each seed candidate
          multimap<float, int >::reverse_iterator rMapIter = seedCandidateMap.rbegin();
          while ( rMapIter != seedCandidateMap.rend() )
            {
              if ( status->adcValues()[ (*rMapIter).second ] == EUTELESCOPE::GOODPIXEL )
                {
                  // if we enter here, this means that at least the seed pixel
                  // wasn't added yet to another cluster.  Note that now we need
                  // to build a candidate cluster that has to pass the
                  // clusterCut to be considered a good cluster
                  ClusterQuality cluQuality = kGoodCluster; // (so far)

                  //gather signal and noise around seed pixel:
                  // prepare a vector to store the noise values
                  vector<float > noiseValueVec;
                  FloatVec clusterCandidateCharges;
                  IntVec   clusterCandidateIndeces;

                  // start looping around the seed pixel. Remember that the seed
                  // pixel has to stay in the center of cluster
                  int seedX, seedY;
                  matrixDecoder.getXYFromIndex ( (*rMapIter).second, seedX, seedY );

                  for (int yPixel = seedY - (_ffYClusterSize / 2); yPixel <= seedY + (_ffYClusterSize / 2); yPixel++)
                    {
                      //!HACK TAKI
                      //streamlog_out ( MESSAGE4 ) <<
                      //" y=" << setw(4) << yPixel << ", x= (" << seedX - (_ffXClusterSize / 2) << "-"<< seedX + (_ffXClusterSize / 2) << ") \t";
                      //!HACK TAKI

                      for (int xPixel =  seedX - (_ffXClusterSize / 2); xPixel <= seedX + (_ffXClusterSize / 2); xPixel++)
                        {
                          // always check we are still within the sensor!!!
                          if ( ( xPixel >= minX )  &&  ( xPixel <= maxX ) &&
                               ( yPixel >= minY )  &&  ( yPixel <= maxY ) )
                            {

                              //linear (1D) index of current (2D) pixel
                              int index = matrixDecoder.getIndexFromXY(xPixel, yPixel);

                              //get noise for each and every pixel!
                              noiseValueVec.push_back(noise->getChargeValues()[ index ]);

                              bool isHit  = ( status->getADCValues()[index] == EUTELESCOPE::HITPIXEL  ); //this is set for pixels already used for another cluster
                              bool isGood = ( status->getADCValues()[index] == EUTELESCOPE::GOODPIXEL );

                              if ( isGood ) //normal case, good means not marked as used by another cluster, yet
                                {
                                  clusterCandidateCharges.push_back( dataVec[ index ] );
                                  clusterCandidateIndeces.push_back( index ); //used to flag used pixels afterwards!

                                  // If the pixel wasn't selected (zs), then its signal
                                  // will still be 0.0001, because dataVec was initialized this way.
                                  // Mark this in the status!
                                  if ( dataVec[ index ] == 0.0001 )
                                    {
                                      status->adcValues()[ index ] = EUTELESCOPE::MISSINGPIXEL;
                                    }

                                  //!HACK TAKI
                                  //streamlog_out ( MESSAGE4 ) << setw(4) << dataVec[ index ];
                                  //!HACK TAKI

                                }
                              else if ( isHit ) //used by another cluster
                                {
                                  // this can be a good place to flag the current
                                  // cluster as kMergedCluster, but it would introduce
                                  // a bias since the at least another cluster (the
                                  // one which this pixel belong to) is not flagged.
                                  //
                                  // In order to flag all merged clusters and possibly
                                  // try to separate the different contributions use
                                  // the EUTelSeparateClusterProcessor. In this
                                  // processor not all the merged clusters will be
                                  // flagged as kMergedCluster | kIncompleteCluster
                                  cluQuality = cluQuality | kIncompleteCluster | kMergedCluster ;
                                  clusterCandidateCharges.push_back(0.0);
                                  clusterCandidateIndeces.push_back( -1 ); //used to flag used pixels afterwards!
                                }
                              else if ( !isGood ) //pixel is bad, firing or missing (when it was first checked) for some reason!
                                {
                                  cluQuality = cluQuality | kIncompleteCluster;
                                  clusterCandidateCharges.push_back(0.0);
                                  clusterCandidateIndeces.push_back( -1 ); //used to flag used pixels afterwards!
                                }
                            }
                          else //the pixel was outside the matrix!
                            {
                              //!no index to push back here, still have to push something in!
                              cluQuality = cluQuality | kBorderCluster;
                              clusterCandidateCharges.push_back(0.0);
                              clusterCandidateIndeces.push_back( -1 ); //used to flag used pixels afterwards!
                              noiseValueVec.push_back(0.0);
                            }
                        }
                      //!HACK TAKI
                      //streamlog_out ( MESSAGE4 ) << endl;
                      //!HACK TAKI
                    } //END  //gathering values  //for (int yPixel = seedY - (_ffYClusterSize / 2); yPixel <= seedY + (_ffYClusterSize / 2); yPixel++)

                  //! build a cluster candidate object from the values obtained
                  if ( ! ( (clusterCandidateCharges.size()==9)&&(noiseValueVec.size()==9) ) )
                    {
                      //the crucial one is clusterCandidateCharges!
                      streamlog_out ( ERROR2 ) << "In event " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                               << " on detector " << sensorID << ":" << endl
                                               << "NOT ENOUGH/TOO MUCH DATA GATHERED TO FORM A 3x3 CLUSTER!! SORRY" << endl
                                               << "There should be 9 noise values,            but there are" << noiseValueVec.size() << "."<< endl
                                               << "There should be 9 signal values,           but there are" << clusterCandidateCharges.size() << "."<< endl
                                               << "There are " << clusterCandidateIndeces.size() << " candidate pixel indeces." << endl;
                      throw IncompatibleDataSetException("NOT ENOUGH/TOO MUCH DATA GATHERED TO FORM A 3x3 CLUSTER");
                    }

                  // the final result of the clustering will enter in a
                  // TrackerPulseImpl in order to be algorithm independent
                  TrackerPulseImpl* pulse = new TrackerPulseImpl; //this will be deleted if the candidate does NOT make it through the cluster cut check, otherwise it will be added to a collection
                  CellIDEncoder<TrackerPulseImpl> idPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
                  idPulseEncoder["sensorID"]      = sensorID;
                  idPulseEncoder["xSeed"]         = seedX;
                  idPulseEncoder["ySeed"]         = seedY;
                  idPulseEncoder["xCluSize"]      = _ffXClusterSize;
                  idPulseEncoder["yCluSize"]      = _ffYClusterSize;
                  //streamlog_out (MESSAGE2) << "  idPulseEncoder:   setting: _ffXClusterSize=" <<  _ffXClusterSize << " _ffYClusterSize " << _ffYClusterSize << endl;
                  idPulseEncoder["type"]          = static_cast<int>(kEUTelBrickedClusterImpl);
                  idPulseEncoder.setCellID(pulse);

                  TrackerDataImpl* clusterData = new TrackerDataImpl; //this will be deleted if the candidate does NOT make it through the cluster cut check, otherwise it will be added to a collection
                  CellIDEncoder<TrackerDataImpl> idClusterEncoder(EUTELESCOPE::CLUSTERDEFAULTENCODING, sparseClusterCollectionVec );
                  idClusterEncoder["sensorID"]      = sensorID;
                  idClusterEncoder["xSeed"]         = seedX;
                  idClusterEncoder["ySeed"]         = seedY;
                  idClusterEncoder["xCluSize"]      = _ffXClusterSize;
                  idClusterEncoder["yCluSize"]      = _ffYClusterSize;
                  //streamlog_out (MESSAGE2) << "  idClusterEncoder:   setting: _ffXClusterSize=" <<  _ffXClusterSize << " _ffYClusterSize " << _ffYClusterSize << endl;
                  idClusterEncoder["quality"]       = static_cast<int>(cluQuality);
                  idClusterEncoder.setCellID(clusterData);


                  //!TAKI:
                  //!I THINK IT IS VERY IMPORTANT TO BE SURE,
                  //!THAT IN clusterCandidateCharges THE PIXELS ARE REALLY SORTED FROM TOP LEFT TO BOTTOM RIGHT (in a 1D way but still)!
                  //!IF THEY ARE NOT, THEN THE WHOLE CLUSTER WILL BE MESSED UP. The for-in-for loop above should ensure that tho.
                  clusterData->setChargeValues(clusterCandidateCharges); //copy data in
                  EUTelBrickedClusterImpl* brickedClusterCandidate = new EUTelBrickedClusterImpl(clusterData); //this will be deleted in any case
                  brickedClusterCandidate->setNoiseValues(noiseValueVec);
                  pulse->setCharge(brickedClusterCandidate->getTotalCharge());

                  //! CUT 2
                  // we need to validate the cluster candidate:
                  //if ( brickedClusterCandidate->getClusterSNR() > _ffClusterCut )
                  if ( brickedClusterCandidate->getClusterSNR(3) > _ffClusterCut ) //!HACK TAKI !! important
                    {
                      //! the cluster candidate is a good cluster
                      //! mark all pixels belonging to the cluster as hit

                      // need to leave a few pixels untouched!!
                      // (again this is only good for the 3x3 implementation with even rows skewed left!!)
                      int pixelNumberInsideClusterToSkip1, pixelNumberInsideClusterToSkip2;
                      if (seedY % 2 == 0) //!NOTE
                        //! IT IS ESSENTIAL HERE, THAT THE INDICES ARE PUSHED IN IN THE CORRECT ORDER!
                        //! (indices representing small y to big y, small x to big x - and y in the outer for-loop, x in the inner)
                        //! x=-1,y=-1 -> x=0,y=-1 -> x=+1,y=-1 -> x=-1,y=0 -> x=0,y=0 -> x=+1,y=0 -> x=-1,y=+1 -> x=0,y=+1 -> x=+1,y=+1
                        //! and there must not be any index missing!!
                        {
                          pixelNumberInsideClusterToSkip1 = 2; //3rd pixel = top right
                          pixelNumberInsideClusterToSkip2 = 8; //9th pixel = bottom right
                        }
                      else
                        {
                          pixelNumberInsideClusterToSkip1 = 0; //1st pixel = top left
                          pixelNumberInsideClusterToSkip2 = 6; //7th pixel = bottom left
                        }

                      IntVec::iterator indexIter = clusterCandidateIndeces.begin();
                      while ( indexIter != clusterCandidateIndeces.end() )
                        {
                          if (
                              (indexIter == clusterCandidateIndeces.begin()+pixelNumberInsideClusterToSkip1)
                              ||
                              (indexIter == clusterCandidateIndeces.begin()+pixelNumberInsideClusterToSkip2)
                              )
                            {
                              //nothing
                            }
                          else
                            {
                              if ( (*indexIter) != -1 )
                                {
                                  status->adcValues()[(*indexIter)] = EUTELESCOPE::HITPIXEL;
                                }
                            }
                          ++indexIter;
                        }

                      //!HACK TAKI DEBUG OUTPUT
                      /*
                        streamlog_out (MESSAGE4) << "=== CREATED A CLUSTER THAT MADE IT THROUGH THE CLUSTER_CUT_CHECK ===" << endl;
                        streamlog_out (MESSAGE4) << " == INFO BEGIN:" << endl;
                        streamlog_out (MESSAGE4) << "  = Run#: " << evt->getRunNumber() << " Evt#: " << evt->getEventNumber() << " Clu_ID: " << clusterID << endl;
                        brickedClusterCandidate->debugOutput();
                        streamlog_out (MESSAGE4) << " == INFO END." << endl;
                      */
                      //!HACK TAKI DEBUG OUTPUT

                      sparseClusterCollectionVec->push_back(clusterData); //taki: don't really understand, what this is good for
                      pulse->setQuality(static_cast<int>(cluQuality));
                      pulse->setTrackerData(clusterData);
                      pulseCollection->push_back(pulse);

                      // increment the cluster counters
                      _totClusterMap[ sensorID ] += 1;
                      ++clusterID;
                      if ( clusterID >= MAXCLUSTERSIZE )
                        {
                          ++limitExceed;
                          --clusterID;
                          streamlog_out ( WARNING2 ) << "Event " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                                     << " on detector " << sensorID
                                                     << " contains more than " << MAXCLUSTERSIZE << " cluster (" << clusterID + limitExceed << ")" << endl;
                        }
                    } //END: if ( brickedClusterCandidate->getClusterSNR() > _ffClusterCut )
                  else
                    {
                      delete clusterData;
                      delete pulse;
                    }

                  delete brickedClusterCandidate;

                } //END: if ( currentSeedpixelcandidate == EUTELESCOPE::GOODPIXEL )

              ++rMapIter;

            } //END: while (not all seed candidates in the map have been processed) ((while ( rMapIter != seedCandidateMap.rend() )))
        } //END: if ( seedCandidateMap.size() != 0 )
    } //for ( unsigned int i = 0 ; i < zsInputDataCollectionVec->size(); i++ )

  // if the sparseClusterCollectionVec isn't empty add it to the
  // current event. The pulse collection will be added afterwards
  if ( ! isDummyAlreadyExisting )
    {
      if ( sparseClusterCollectionVec->size() != 0 )
        {
          evt->addCollection( sparseClusterCollectionVec, "original_zsdata" );
        }
      else
        {
          delete sparseClusterCollectionVec;
        }
    }

}




void EUTelClusteringProcessor::sparseClustering(LCEvent* evt, LCCollectionVec* pulseCollection)
{
	//only for noise
	CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );



	// prepare some decoders
	CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );

	bool isDummyAlreadyExisting = false;
	LCCollectionVec* sparseClusterCollectionVec = NULL;

	try
	{
		sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "original_zsdata") );
		isDummyAlreadyExisting = true ;
	}
	catch (lcio::DataNotAvailableException& e)
	{
		sparseClusterCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
		isDummyAlreadyExisting = false;
	}

	CellIDEncoder<TrackerDataImpl> idZSClusterEncoder( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );

	// prepare an encoder also for the pulse collection
	CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);

	// in the zsInputDataCollectionVec we should have one TrackerData for each
	// detector working in ZS mode. We need to loop over all of them
	for ( unsigned int idetector = 0 ; idetector < zsInputDataCollectionVec->size(); idetector++ )
	{
		// get the TrackerData and guess which kind of sparsified data it contains.
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( idetector ) );
		SparsePixelType type = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
		int sensorID = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
	    

		//if this is an excluded sensor go to the next element
		bool foundexcludedsensor = false;
		for(size_t iexclude = 0; iexclude < _ExcludedPlanes.size(); ++iexclude)
		{
			if(_ExcludedPlanes[iexclude] == sensorID)
			{
				foundexcludedsensor = true;
			}
		}
		if(foundexcludedsensor)
		{	
			continue;
		}


		if ( type == kEUTelGenericSparsePixel )
		{

			// now prepare the EUTelescope interface to sparsified data.
			auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel > > sparseData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> ( zsData ) );

			streamlog_out ( DEBUG2 ) << "Processing sparse data on detector " << sensorID << " with " << sparseData->size() << " pixels " << endl;

			int hitPixelsInEvent = sparseData->size();
			std::vector<EUTelGenericSparsePixel> hitPixelVec;
			EUTelGenericSparsePixel* pixel = new EUTelGenericSparsePixel;

			//This for-loop loads all the hits of the given event and detector plane and stores them
			for(int i = 0; i < hitPixelsInEvent; ++i )
			{
				//Load the information of the hit pixel into genericPixel
				sparseData->getSparsePixelAt( i, pixel );
				EUTelGenericSparsePixel hitPixel( *pixel );

				//and push this pixel back
				hitPixelVec.push_back( hitPixel );
			}	

			std::vector<EUTelGenericSparsePixel> newlyAdded;
			//We now cluster those hits together
			while( !hitPixelVec.empty() )
			{
				// prepare a TrackerData to store the cluster candidate
				auto_ptr< TrackerDataImpl > zsCluster ( new TrackerDataImpl );
				// prepare a reimplementation of sparsified cluster
				auto_ptr<EUTelSparseClusterImpl<EUTelGenericSparsePixel > > sparseCluster ( new EUTelSparseClusterImpl<EUTelGenericSparsePixel>( zsCluster.get() ) );

				std::vector<EUTelGenericSparsePixel> cluCandidate;

				//First we need to take any pixel, so let's take the first one
				//Add it to the cluster as well as the newly added pixels
				newlyAdded.push_back( hitPixelVec.front() );
				//sparseCluster->addSparsePixel( &(hitPixelVec.front()) );
				cluCandidate.push_back( hitPixelVec.front() );
				//And remove it from the original collection
				hitPixelVec.erase( hitPixelVec.begin() );

				//Now process all newly added pixels, initially this is the just previously added one
				//but in the process of neighbour finding we continue to add new pixels
				while( !newlyAdded.empty() )
				{
					bool newlyDone = true;
					int  x1, x2, y1, y2, dX, dY;

					//check against all pixels in the hitPixelVec
					for( std::vector<EUTelGenericSparsePixel>::iterator hitVec = hitPixelVec.begin(); hitVec != hitPixelVec.end(); ++hitVec )
					{
						//get the relevant infos from the newly added pixel
						x1 = newlyAdded.front().getXCoord();
						y1 = newlyAdded.front().getYCoord();

						//and the pixel we test against
						x2 = hitVec->getXCoord();
						y2 = hitVec->getYCoord();

						dX = x1 - x2;
						dY = y1 - y2;
						int distance = dX*dX+dY*dY;
						//if they pass the spatial and temporal cuts, we add them
						if( distance <= _sparseMinDistanceSquared )
						{
							//add them to the cluster as well as to the newly added ones
							newlyAdded.push_back( *hitVec );
							cluCandidate.push_back( *hitVec );
						//	sparseCluster->addSparsePixel( &(*hitVec) );
							//and remove it from the original collection
							hitPixelVec.erase( hitVec );
							//for the pixel we test there might be other neighbours, we still have to check
							newlyDone = false;
							break;
						}
					}

					//if no neighbours are found, we can delete the pixel from the newly added
					//we tested against _ALL_ non cluster pixels, there are no other pixels
					//which could be neighbours
					if(newlyDone) newlyAdded.erase( newlyAdded.begin() );
				}
				
				// get the noise matrix with the right detectorID
				TrackerDataImpl* noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));
				// prepare the matrix decoder
				EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );
				// prepare a vector to store the noise values
				vector<float> noiseValueVec;

				//Hot pixel removement:
				while(!cluCandidate.empty())
				{
					EUTelGenericSparsePixel pixel = cluCandidate.front();
					cluCandidate.erase( cluCandidate.begin() );

          				int index = matrixDecoder.getIndexFromXY( pixel.getXCoord(), pixel.getYCoord() );
	  				if( _hitIndexMapVec[idetector].find( index ) != _hitIndexMapVec[idetector].end() )
					{
						// do nothing
					}
					else
					{
						sparseCluster->addSparsePixel( &pixel );
						noiseValueVec.push_back(noise->getChargeValues()[ index ]);
					}
				}

				sparseCluster->setNoiseValues( noiseValueVec );

				//Now we need to process the found cluster
				if ( (sparseCluster->size() > 0) && (sparseCluster->getSeedSNR() >= _sparseSeedCut) && (sparseCluster->getClusterSNR() >= _sparseClusterCut) )
				{
					// set the ID for this zsCluster
					idZSClusterEncoder["sensorID"] = sensorID;
					idZSClusterEncoder["sparsePixelType"] = static_cast<int>( type );
					idZSClusterEncoder["quality"] = 0;
					idZSClusterEncoder.setCellID( zsCluster.get() );

					// add it to the cluster collection
					sparseClusterCollectionVec->push_back( zsCluster.get() );

					// prepare a pulse for this cluster
					auto_ptr<TrackerPulseImpl> zsPulse ( new TrackerPulseImpl );
					idZSPulseEncoder["sensorID"] = sensorID;
					idZSPulseEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
					idZSPulseEncoder.setCellID( zsPulse.get() );

					//zsPulse->setCharge( sparseCluster->getTotalCharge() );
					zsPulse->setTrackerData( zsCluster.release() );
					pulseCollection->push_back( zsPulse.release() );

					// last but not least increment the totClusterMap
					_totClusterMap[ sensorID ] += 1;
				} //cluster processing if

				else
				{
					//in the case the cluster candidate is not passing the threshold ...
					//forget about them, the memory should be automatically cleaned by auto_ptr's
				}
			} //loop over all found clusters

			delete pixel;
		}
		else
		{
		     throw UnknownDataTypeException("Unknown sparsified pixel");
		}
	} // this is the end of the loop over all ZS detectors

	// if the sparseClusterCollectionVec isn't empty add it to the
	// current event. The pulse collection will be added afterwards
	if ( ! isDummyAlreadyExisting )
	{
		if ( sparseClusterCollectionVec->size() != 0 )
		{
			evt->addCollection( sparseClusterCollectionVec, "original_zsdata" );
		}
		else
		{
			delete sparseClusterCollectionVec;
		}
	}
}


void EUTelClusteringProcessor::fixedFrameClustering(LCEvent * evt, LCCollectionVec * pulseCollection) {

  streamlog_out ( DEBUG4 ) << "Looking for clusters in the RAW frame with FixedFrame algorithm " << endl;

//  LCCollectionVec * nzsInputDataCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_nzsDataCollectionName));
//  LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
//  LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));

  CellIDDecoder<TrackerDataImpl> cellDecoder( nzsInputDataCollectionVec );

  if (isFirstEvent()) {

    // check if for each TrackerData into the NZS corresponds one
    // TrackerData with noise information having the same number of
    // pixels.

    for ( unsigned int i = 0 ; i < nzsInputDataCollectionVec->size(); i++ ) {
      TrackerDataImpl    * nzsData = dynamic_cast<TrackerDataImpl* > ( nzsInputDataCollectionVec->getElementAt( i ) );
      int detectorID     = cellDecoder( nzsData ) ["sensorID"];
      //if this is an excluded sensor go to the next element
      bool foundexcludedsensor = false;
      for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
        {
          if(_ExcludedPlanes[i] == detectorID)
            {
              foundexcludedsensor = true;
            }
        }
      if(foundexcludedsensor)
        continue;
      TrackerDataImpl    * noise   = dynamic_cast<TrackerDataImpl* >    ( noiseCollectionVec->getElementAt( _ancillaryIndexMap[ detectorID ] ) );
      TrackerRawDataImpl * status  = dynamic_cast<TrackerRawDataImpl *> ( statusCollectionVec->getElementAt( _ancillaryIndexMap[ detectorID ] ) );


      if ( ( noise->chargeValues().size() != status->adcValues().size() ) ||
           ( noise->chargeValues().size() != nzsData->getChargeValues().size() ) ) {
        throw IncompatibleDataSetException("NZS data and noise/status size mismatch");
      }
    }


    /// /* DEBUG */    logfile.open("clustering.log");
  }

  streamlog_out ( DEBUG0 ) << "Event " << _iEvt << endl;

  bool isDummyAlreadyExisting = false;
  LCCollectionVec * dummyCollection = NULL;
  try {
    dummyCollection = dynamic_cast< LCCollectionVec* > ( evt->getCollection( _dummyCollectionName ) );
    isDummyAlreadyExisting = true;

  } catch (lcio::DataNotAvailableException& e) {
    dummyCollection =  new LCCollectionVec(LCIO::TRACKERDATA);
    isDummyAlreadyExisting = false;
  }

  for ( int i = 0; i < nzsInputDataCollectionVec->getNumberOfElements(); i++) {

    // get the calibrated data
    TrackerDataImpl    * nzsData = dynamic_cast<TrackerDataImpl*>  (nzsInputDataCollectionVec->getElementAt( i ) );
    int sensorID                 = cellDecoder( nzsData ) ["sensorID"];
    //if this is an excluded sensor go to the next element
    bool foundexcludedsensor = false;
    for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
      {
        if(_ExcludedPlanes[i] == sensorID)
          {
            foundexcludedsensor = true;
          }
      }
    if(foundexcludedsensor)
      continue;
    // now that we know which is the sensorID, we can ask to GEAR
    // which are the minX, minY, maxX and maxY.
    int minX, minY, maxX, maxY;
    minX = 0;
    minY = 0;

    // this sensorID can be either a reference plane or a DUT, do it
    // differently...
    if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() ){
      // this is a reference plane
      maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
      maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
    } else if ( _dutLayerIndexMap.find( sensorID ) != _dutLayerIndexMap.end() ) {
      // ok it is a DUT plane
      maxX = _siPlanesLayerLayout->getDUTSensitiveNpixelX() - 1;
      maxY = _siPlanesLayerLayout->getDUTSensitiveNpixelY() - 1;
    } else {
      // this is not a reference plane neither a DUT... what's that?
//      throw  InvalidGeometryException ("Unknown sensorID " + to_string( sensorID ));
      streamlog_out( ERROR5 ) << "Unknown sensorID " << sensorID << ", perhaps your GEAR file is incomplete." << endl;
      continue;
    }

    streamlog_out ( DEBUG0 ) << "  Working on detector " << sensorID << endl;

    TrackerDataImpl    * noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));
    TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));

    // prepare the matrix decoder
    EUTelMatrixDecoder matrixDecoder(cellDecoder, nzsData);

    // reset the status
    resetStatus(status);

    // initialize the cluster counter
    short clusterCounter = 0;
    short limitExceed    = 0;

    _seedCandidateMap.clear();

    for (unsigned int iPixel = 0; iPixel < nzsData->getChargeValues().size(); iPixel++) 
    {
        if (status->getADCValues()[iPixel] == EUTELESCOPE::GOODPIXEL) 
        {
            if ( nzsData->getChargeValues()[iPixel] > _ffSeedCut * noise->getChargeValues()[iPixel]) 
            {
                _seedCandidateMap.push_back(make_pair( nzsData->getChargeValues()[iPixel], iPixel));
            }
        }
    }

    // continue only if seed candidate map is not empty!
    if ( !_seedCandidateMap.empty() ) {

      streamlog_out ( DEBUG0 ) << "There are << " << _seedCandidateMap.size() << " seed candidates." << endl;

      // now built up a cluster for each seed candidate
      //Start by sorting the vector in order of smallest to largest seed signal
      std::sort(_seedCandidateMap.begin(),_seedCandidateMap.end());
      vector< pair< float, unsigned int > >::iterator mapIter = _seedCandidateMap.end();
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
          matrixDecoder.getXYFromIndex((*mapIter).second,seedX, seedY);

          // start looping around the seed pixel. Remember that the seed
          // pixel has to stay in the center of cluster
          ClusterQuality cluQuality = kGoodCluster;
          for (int yPixel = seedY - (_ffYClusterSize / 2); yPixel <= seedY + (_ffYClusterSize / 2); yPixel++) {
            for (int xPixel =  seedX - (_ffXClusterSize / 2); xPixel <= seedX + (_ffXClusterSize / 2); xPixel++) {
              // always check we are still within the sensor!!!
              if ( ( xPixel >= minX )  &&  ( xPixel <= maxX ) &&
                   ( yPixel >= minY )  &&  ( yPixel <= maxY ) ) {
                int index = matrixDecoder.getIndexFromXY(xPixel, yPixel);

                bool isHit  = ( status->getADCValues()[index] == EUTELESCOPE::HITPIXEL  );
                bool isGood = ( status->getADCValues()[index] == EUTELESCOPE::GOODPIXEL );
                
                if(isGood)
                  clusterCandidateIndeces.push_back(index);
                else
                  clusterCandidateIndeces.push_back(-1);
                
                if ( isGood && !isHit ) {
                  clusterCandidateSignal += nzsData->getChargeValues()[index];
                  clusterCandidateNoise2 += pow(noise->getChargeValues()[index] , 2);
                  clusterCandidateCharges.push_back(nzsData->getChargeValues()[index]);
                } else if (isHit) {
                  // this can be a good place to flag the current
                  // cluster as kMergedCluster, but it would introduce
                  // a bias since the at least another cluster (the
                  // one which this pixel belong to) is not flagged.
                  //
                  // In order to flag all merged clusters and possibly
                  // try to separate the different contributions use
                  // the EUTelSeparateClusterProcessor. In this
                  // processor not all the merged clusters will be
                  // flagged as kMergedCluster | kIncompleteCluster
                  cluQuality = cluQuality | kIncompleteCluster | kMergedCluster ;
                  clusterCandidateCharges.push_back(0.);
                } else if (!isGood) {
                  cluQuality = cluQuality | kIncompleteCluster;
                  clusterCandidateCharges.push_back(0.);
                }
              } else {
                cluQuality = cluQuality | kBorderCluster;
                clusterCandidateCharges.push_back(0.);
                clusterCandidateIndeces.push_back( -1 ) ;
              }
            }
          }

          // at this point we have built the cluster candidate,
          // we need to validate it
          if ( clusterCandidateSignal > _ffClusterCut * sqrt(clusterCandidateNoise2) ) {
            // the cluster candidate is a good cluster
            // mark all pixels belonging to the cluster as hit
            IntVec::iterator indexIter = clusterCandidateIndeces.begin();

            // the final result of the clustering will enter in a
            // TrackerPulseImpl in order to be algorithm independent
            TrackerPulseImpl * pulse = new TrackerPulseImpl;
            CellIDEncoder<TrackerPulseImpl> idPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
            idPulseEncoder["sensorID"]      = sensorID;
            idPulseEncoder["xSeed"]         = seedX;
            idPulseEncoder["ySeed"]         = seedY;
            idPulseEncoder["xCluSize"]      = _ffXClusterSize;
            idPulseEncoder["yCluSize"]      = _ffYClusterSize;
            idPulseEncoder["type"]          = static_cast<int>(kEUTelFFClusterImpl);
            idPulseEncoder.setCellID(pulse);


            TrackerDataImpl * cluster = new TrackerDataImpl;
            CellIDEncoder<TrackerDataImpl> idClusterEncoder(EUTELESCOPE::CLUSTERDEFAULTENCODING, dummyCollection);
            idClusterEncoder["sensorID"]      = sensorID;
            idClusterEncoder["xSeed"]         = seedX;
            idClusterEncoder["ySeed"]         = seedY;
            idClusterEncoder["xCluSize"]      = _ffXClusterSize;
            idClusterEncoder["yCluSize"]      = _ffYClusterSize;
            idClusterEncoder["quality"]       = static_cast<int>(cluQuality);
            idClusterEncoder.setCellID(cluster);

            streamlog_out (DEBUG0) << "  Cluster no " <<  clusterCounter << " seedX " << seedX << " seedY " << seedY << endl;

            while ( indexIter != clusterCandidateIndeces.end() ) {
              if (*indexIter != -1 ) {
                status->adcValues()[(*indexIter)] = EUTELESCOPE::HITPIXEL;
              }
              ++indexIter;
            }

            // copy the candidate charges inside the cluster
            cluster->setChargeValues(clusterCandidateCharges);
            dummyCollection->push_back(cluster);

            EUTelFFClusterImpl * eutelCluster = new EUTelFFClusterImpl( cluster );
            pulse->setCharge(eutelCluster->getTotalCharge());
            delete eutelCluster;

            pulse->setQuality(static_cast<int>(cluQuality));
            pulse->setTrackerData(cluster);
            pulseCollection->push_back(pulse);

            // increment the cluster counters
            _totClusterMap[ sensorID ] += 1;
            ++clusterCounter;
            if ( clusterCounter > MAXCLUSTERSIZE ) {
              ++limitExceed;
              --clusterCounter;
              streamlog_out ( WARNING2 ) << "Event " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                         << " on detector " << sensorID
                                         << " contains more than " << MAXCLUSTERSIZE << " cluster (" << clusterCounter + limitExceed << ")" << endl;
            }
          } else {
            // the cluster has not passed the cut!

          }
        }
      }
    }
  }

  if ( ! isDummyAlreadyExisting ) {
    if ( dummyCollection->size() != 0 ) {
      evt->addCollection(dummyCollection,_dummyCollectionName);
    } else {
      delete dummyCollection;
    }
  }

}

void EUTelClusteringProcessor::nzsBrickedClustering(LCEvent * evt, LCCollectionVec * pulseCollection)
{
  streamlog_out ( DEBUG4 ) << "Looking for clusters in the RAW data with nzsBrickedClustering algorithm " << endl;

  // get the collections of interest from the event.
//  LCCollectionVec * nzsInputDataCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _nzsDataCollectionName ));
//  LCCollectionVec * noiseCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
//  LCCollectionVec * statusCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));

  CellIDDecoder<TrackerDataImpl> cellDecoder( nzsInputDataCollectionVec );

  if (isFirstEvent())
    {
      // check if for each TrackerData into the NZS corresponds one
      // TrackerData with noise information having the same number of
      // pixels.

      for ( unsigned int i = 0 ; i < nzsInputDataCollectionVec->size(); i++ )
        {
          TrackerDataImpl    * nzsData = dynamic_cast<TrackerDataImpl* > ( nzsInputDataCollectionVec->getElementAt( i ) );
          int detectorID     = cellDecoder( nzsData ) ["sensorID"];

          TrackerDataImpl    * noise   = dynamic_cast<TrackerDataImpl* >    ( noiseCollectionVec->getElementAt( _ancillaryIndexMap[ detectorID ] ) );
          TrackerRawDataImpl * status  = dynamic_cast<TrackerRawDataImpl *> ( statusCollectionVec->getElementAt( _ancillaryIndexMap[ detectorID ] ) );

          if ( ( noise->chargeValues().size() != status->adcValues().size() ) ||
               ( noise->chargeValues().size() != nzsData->chargeValues().size() ) )
            {
              throw IncompatibleDataSetException("NZS data and noise/status size mismatch");
            }
        }
      //DEBUG:
      //logfile.open("clustering.log");
    }

  streamlog_out ( DEBUG0 ) << "nzsBrickedClustering: Event " << _iEvt << endl;


  bool isDummyAlreadyExisting = false;
  LCCollectionVec * dummyCollection = NULL;
  try
    {
      dummyCollection = dynamic_cast< LCCollectionVec* > ( evt->getCollection( _dummyCollectionName ) );
      isDummyAlreadyExisting = true;
    }
  catch (lcio::DataNotAvailableException& e)
    {
      dummyCollection =  new LCCollectionVec(LCIO::TRACKERDATA);
      isDummyAlreadyExisting = false;
    }

  // prepare an encoder also for the pulse collection
  CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);


  // utility
  short limitExceed    = 0;
  if ( isFirstEvent() )
    {
      // For the time being nothing to do specifically in the first
      // event.
    }


  for ( int i = 0; i < nzsInputDataCollectionVec->getNumberOfElements(); i++)
    {
      // get the calibrated data
      TrackerDataImpl    * nzsData = dynamic_cast<TrackerDataImpl*>  (nzsInputDataCollectionVec->getElementAt( i ) );
      EUTelMatrixDecoder matrixDecoder(cellDecoder, nzsData);
      int sensorID                 = static_cast<int > ( cellDecoder( nzsData )["sensorID"] );

      // now that we know which is the sensorID, we can ask to GEAR
      // which are the minX, minY, maxX and maxY.
      int minX, minY, maxX, maxY;
      minX = 0;
      minY = 0;

      // this sensorID can be either a reference plane or a DUT, do it
      // differently...
      if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() )
        {
          // this is a reference plane
          maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
          maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
        }
      else if ( _dutLayerIndexMap.find( sensorID ) != _dutLayerIndexMap.end() )
        {
          // ok it is a DUT plane
          maxX = _siPlanesLayerLayout->getDUTSensitiveNpixelX() - 1;
          maxY = _siPlanesLayerLayout->getDUTSensitiveNpixelY() - 1;
        }
      else
        {
          // this is not a reference plane neither a DUT... what's that?
//          throw  InvalidGeometryException ("Unknown sensorID " + to_string( sensorID ));
          streamlog_out( ERROR5 ) << "Unknown sensorID " << sensorID << ", perhaps your GEAR file is incomplete." << endl;
          continue;
        }

      // get the noise and the status matrix with the right detectorID
      TrackerDataImpl    * noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));
      TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));

      // reset the cluster counter for the clusterID
      int clusterID = 0;

      // reset the status
      resetStatus(status);

      // prepare a multimap for the seed candidates
      vector< pair <float, int> > seedCandidateMap;

      // fill the seed candidate map
      for (unsigned int iPixel = 0; iPixel < nzsData->getChargeValues().size(); iPixel++)
        {
          if (status->getADCValues()[iPixel] == EUTELESCOPE::GOODPIXEL)
            {
              //! CUT 1
              if ( nzsData->getChargeValues()[iPixel] > _ffSeedCut * noise->getChargeValues()[iPixel])
                {
                  seedCandidateMap.push_back(make_pair( nzsData->getChargeValues()[iPixel], iPixel));
                  streamlog_out ( MESSAGE2 )
                    << "Added pixel at (index=" << iPixel
                    << ") with signal " << nzsData->getChargeValues()[iPixel]
                    << " to the seedCandidateMap" << endl;

                  if ( noise->getChargeValues()[ iPixel ] < 0.01 )
                    {
                      streamlog_out ( ERROR2 )    << "ZERO NOISE SEED PIXEL ADDED (nszBrickedClustering)!"
                                                  << "\n index=" << iPixel
                                                  << "\n amp=" << nzsData->getChargeValues()[ iPixel ]
                                                  << "\n status=" << status->getADCValues()[ iPixel ]
                                                  <<    " GOODP   =  0,"
                                                  <<    " BAD     =  1,"
                                                  <<    " HIT     = -1,"
                                                  <<    " MISSING =  2,"
                                                  <<    " FIRING  =  3.";
                    }
                }
            }
        }

      streamlog_out ( DEBUG0 ) << "The number of seed candidates is: " << seedCandidateMap.size() << endl;
      std::sort(seedCandidateMap.begin(),seedCandidateMap.end());
      if ( !seedCandidateMap.empty() )
        {
          // now build up a cluster for each seed candidate
          vector< pair<float, int> >::iterator rMapIter = seedCandidateMap.end();
          while ( rMapIter != seedCandidateMap.begin() )
            {
              rMapIter--;
              if ( status->adcValues()[ (*rMapIter).second ] == EUTELESCOPE::GOODPIXEL )
                {
                  // if we enter here, this means that at least the seed pixel
                  // wasn't added yet to another cluster.  Note that now we need
                  // to build a candidate cluster that has to pass the
                  // clusterCut to be considered a good cluster
                  ClusterQuality cluQuality = kGoodCluster; // (so far)

                  //gather signal and noise around seed pixel:
                  // prepare a vector to store the noise values
                  vector<float > noiseValueVec;
                  FloatVec clusterCandidateCharges;
                  IntVec   clusterCandidateIndeces;

                  // start looping around the seed pixel. Remember that the seed
                  // pixel has to stay in the center of cluster
                  int seedX, seedY;
                  matrixDecoder.getXYFromIndex ( (*rMapIter).second, seedX, seedY );

                  for (int yPixel = seedY - (_ffYClusterSize / 2); yPixel <= seedY + (_ffYClusterSize / 2); yPixel++)
                    {
                      //!HACK TAKI
                      //streamlog_out ( MESSAGE4 ) <<
                      //" y=" << setw(4) << yPixel << ", x= (" << seedX - (_ffXClusterSize / 2) << "-"<< seedX + (_ffXClusterSize / 2) << ") \t";
                      //!HACK TAKI

                      for (int xPixel =  seedX - (_ffXClusterSize / 2); xPixel <= seedX + (_ffXClusterSize / 2); xPixel++)
                        {
                          // always check we are still within the sensor!!!
                          if ( ( xPixel >= minX )  &&  ( xPixel <= maxX ) &&
                               ( yPixel >= minY )  &&  ( yPixel <= maxY ) )
                            {

                              //linear (1D) index of current (2D) pixel
                              int index = matrixDecoder.getIndexFromXY(xPixel, yPixel);

                              //get noise for each and every pixel! (because it's available)
                              noiseValueVec.push_back(noise->getChargeValues()[ index ]);

                              bool isHit  = ( status->getADCValues()[index] == EUTELESCOPE::HITPIXEL  ); //this is set for pixels already used for another cluster
                              bool isGood = ( status->getADCValues()[index] == EUTELESCOPE::GOODPIXEL );

                              if ( isGood ) //normal case
                                {
                                  clusterCandidateCharges.push_back( nzsData->getChargeValues()[index] );
                                  clusterCandidateIndeces.push_back( index ); //used to flag used pixels afterwards!

                                  //no need to check the signal presence here
                                  //(compare to zs (bricked) clustering where we take a look, if dataVec[...] has changed here!)

                                  //!HACK TAKI
                                  //streamlog_out ( MESSAGE4 ) << setw(4) << dataVec[ index ];
                                  //!HACK TAKI
                                }
                              else if ( isHit )
                                {
                                  // this can be a good place to flag the current
                                  // cluster as kMergedCluster, but it would introduce
                                  // a bias since the at least another cluster (the
                                  // one which this pixel belong to) is not flagged.
                                  //
                                  // In order to flag all merged clusters and possibly
                                  // try to separate the different contributions use
                                  // the EUTelSeparateClusterProcessor. In this
                                  // processor not all the merged clusters will be
                                  // flagged as kMergedCluster | kIncompleteCluster
                                  clusterCandidateCharges.push_back(0.0);
                                  clusterCandidateIndeces.push_back( -1 ); //used to flag used pixels afterwards!
                                  cluQuality = cluQuality | kIncompleteCluster | kMergedCluster ;
                                }
                              else if ( !isGood ) //pixel is bad, firing or missing (when it was first checked) for some reason!
                                {
                                  clusterCandidateCharges.push_back(0.0);
                                  clusterCandidateIndeces.push_back( -1 ); //used to flag used pixels afterwards!
                                  cluQuality = cluQuality | kIncompleteCluster;
                                }
                            }
                          else //the pixel was outside the matrix!
                            {
                              //!no index to push back here, still have to push something in!
                              clusterCandidateCharges.push_back(0.0);
                              clusterCandidateIndeces.push_back( -1 ); //used to flag used pixels afterwards!
                              cluQuality = cluQuality | kBorderCluster;
                              noiseValueVec.push_back(0.0);
                            }
                        }
                      //!HACK TAKI
                      //streamlog_out ( MESSAGE4 ) << endl;
                      //!HACK TAKI
                    } //END  //gathering values  //for (int yPixel = seedY - (_ffYClusterSize / 2); yPixel <= seedY + (_ffYClusterSize / 2); yPixel++)

                  //! build a cluster candidate object from the values obtained
                  if ( ! ( (clusterCandidateCharges.size()==9)&&(noiseValueVec.size()==9) ) )
                    {
                      //the crucial one is clusterCandidateCharges!
                      streamlog_out ( ERROR2 ) << "In event " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                               << " on detector " << sensorID << ":" << endl
                                               << "NOT ENOUGH/TOO MUCH DATA GATHERED TO FORM A 3x3 CLUSTER!! SORRY" << endl
                                               << "There should be 9 noise values,            but there are" << noiseValueVec.size() << "."<< endl
                                               << "There should be 9 signal values,           but there are" << clusterCandidateCharges.size() << "."<< endl
                                               << "There are " << clusterCandidateIndeces.size() << " candidate pixel indeces." << endl;
                      throw IncompatibleDataSetException("NOT ENOUGH/TOO MUCH DATA GATHERED TO FORM A 3x3 CLUSTER");
                    }

                  // the final result of the clustering will enter in a
                  // TrackerPulseImpl in order to be algorithm independent
                  TrackerPulseImpl* pulse = new TrackerPulseImpl; //this will be deleted if the candidate does NOT make it through the cluster cut check, otherwise it will be added to a collection
                  CellIDEncoder<TrackerPulseImpl> idPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
                  idPulseEncoder["sensorID"]      = sensorID;
                  idPulseEncoder["xSeed"]         = seedX;
                  idPulseEncoder["ySeed"]         = seedY;
                  idPulseEncoder["xCluSize"]      = _ffXClusterSize;
                  idPulseEncoder["yCluSize"]      = _ffYClusterSize;
                  //streamlog_out (MESSAGE2) << "  idPulseEncoder:   setting: _ffXClusterSize=" <<  _ffXClusterSize << " _ffYClusterSize " << _ffYClusterSize << endl;
                  idPulseEncoder["type"]          = static_cast<int>(kEUTelBrickedClusterImpl);
                  idPulseEncoder.setCellID(pulse);

                  TrackerDataImpl* clusterData = new TrackerDataImpl; //this will be deleted if the candidate does NOT make it through the cluster cut check, otherwise it will be added to a collection
                  CellIDEncoder<TrackerDataImpl> idClusterEncoder(EUTELESCOPE::CLUSTERDEFAULTENCODING, dummyCollection ); //GOBACK
                  idClusterEncoder["sensorID"]      = sensorID;
                  idClusterEncoder["xSeed"]         = seedX;
                  idClusterEncoder["ySeed"]         = seedY;
                  idClusterEncoder["xCluSize"]      = _ffXClusterSize;
                  idClusterEncoder["yCluSize"]      = _ffYClusterSize;
                  //streamlog_out (MESSAGE2) << "  idClusterEncoder:   setting: _ffXClusterSize=" <<  _ffXClusterSize << " _ffYClusterSize " << _ffYClusterSize << endl;
                  idClusterEncoder["quality"]       = static_cast<int>(cluQuality);
                  idClusterEncoder.setCellID(clusterData);


                  //!TAKI:
                  //!I THINK IT IS VERY IMPORTANT TO BE SURE,
                  //!THAT IN clusterCandidateCharges THE PIXELS ARE REALLY SORTED FROM TOP LEFT TO BOTTOM RIGHT (in a 1D way but still)!
                  //!IF THEY ARE NOT, THEN THE WHOLE CLUSTER WILL BE MESSED UP. The for-in-for loop above should ensure that tho.
                  clusterData->setChargeValues(clusterCandidateCharges); //copy data in
                  EUTelBrickedClusterImpl* brickedClusterCandidate = new EUTelBrickedClusterImpl(clusterData); //this will be deleted in any case
                  brickedClusterCandidate->setNoiseValues(noiseValueVec);
                  pulse->setCharge(brickedClusterCandidate->getTotalCharge());

                  //! CUT 2
                  // we need to validate the cluster candidate:
                  //if ( brickedClusterCandidate->getClusterSNR() > _ffClusterCut )
                  if ( brickedClusterCandidate->getClusterSNR(3) > _ffClusterCut ) //!HACK TAKI !! important
                    {
                      //! the cluster candidate is a good cluster
                      //! mark all pixels belonging to the cluster as hit

                      // need to leave a few pixels untouched!!
                      // (again this is only good for the 3x3 implementation with even rows skewed left!!)
                      int pixelNumberInsideClusterToSkip1, pixelNumberInsideClusterToSkip2;
                      if (seedY % 2 == 0) //!NOTE
                        //! IT IS ESSENTIAL HERE, THAT THE INDICES ARE PUSHED IN IN THE CORRECT ORDER!
                        //! (indices representing small y to big y, small x to big x - and y in the outer for-loop, x in the inner)
                        //! x=-1,y=-1 -> x=0,y=-1 -> x=+1,y=-1 -> x=-1,y=0 -> x=0,y=0 -> x=+1,y=0 -> x=-1,y=+1 -> x=0,y=+1 -> x=+1,y=+1
                        //! and there must not be any index missing!!
                        {
                          pixelNumberInsideClusterToSkip1 = 2; //3rd pixel = top right
                          pixelNumberInsideClusterToSkip2 = 8; //9th pixel = bottom right
                        }
                      else
                        {
                          pixelNumberInsideClusterToSkip1 = 0; //1st pixel = top left
                          pixelNumberInsideClusterToSkip2 = 6; //7th pixel = bottom left
                        }

                      IntVec::iterator indexIter = clusterCandidateIndeces.begin();
                      while ( indexIter != clusterCandidateIndeces.end() )
                        {
                          if (
                              (indexIter == clusterCandidateIndeces.begin()+pixelNumberInsideClusterToSkip1)
                              ||
                              (indexIter == clusterCandidateIndeces.begin()+pixelNumberInsideClusterToSkip2)
                              )
                            {
                              //nothing
                            }
                          else
                            {
                              if ( (*indexIter) != -1 )
                                {
                                  status->adcValues()[(*indexIter)] = EUTELESCOPE::HITPIXEL;
                                }
                            }
                          ++indexIter;
                        }

                      //!HACK TAKI DEBUG OUTPUT
                      /*
                        streamlog_out (MESSAGE4) << "=== CREATED A CLUSTER THAT MADE IT THROUGH THE CLUSTER_CUT_CHECK ===" << endl;
                        streamlog_out (MESSAGE4) << " == INFO BEGIN:" << endl;
                        streamlog_out (MESSAGE4) << "  = Run#: " << evt->getRunNumber() << " Evt#: " << evt->getEventNumber() << " Clu_ID: " << clusterID << endl;
                        brickedClusterCandidate->debugOutput();
                        streamlog_out (MESSAGE4) << " == INFO END." << endl;
                      */
                      //!HACK TAKI DEBUG OUTPUT

                      dummyCollection->push_back(clusterData); //taki: don't really understand, what this is good for
                      pulse->setQuality(static_cast<int>(cluQuality));
                      pulse->setTrackerData(clusterData);
                      pulseCollection->push_back(pulse);

                      // increment the cluster counters
                      _totClusterMap[ sensorID ] += 1;
                      ++clusterID;
                      if ( clusterID >= MAXCLUSTERSIZE )
                        {
                          ++limitExceed;
                          --clusterID;
                          streamlog_out ( WARNING2 ) << "Event " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                                     << " on detector " << sensorID
                                                     << " contains more than " << MAXCLUSTERSIZE << " cluster (" << clusterID + limitExceed << ")" << endl;
                        }
                    } //END: if ( brickedClusterCandidate->getClusterSNR() > _ffClusterCut )
                  else
                    {
                      delete clusterData;
                      delete pulse;
                    }

                  delete brickedClusterCandidate;

                } //END: if ( currentSeedpixelcandidate == EUTELESCOPE::GOODPIXEL )
            } //END: while (not all seed candidates in the map have been processed) ((while ( rMapIter != seedCandidateMap.rend() )))
        } //END: if ( seedCandidateMap.size() != 0 )
    } //for ( unsigned int i = 0 ; i < zsInputDataCollectionVec->size(); i++ )

  // if the sparseClusterCollectionVec isn't empty add it to the
  // current event. The pulse collection will be added afterwards
  if ( ! isDummyAlreadyExisting )
    {
      if ( dummyCollection->size() != 0 )
        {
          evt->addCollection( dummyCollection, _dummyCollectionName );
        }
      else
        {
          delete dummyCollection;
        }
    }
}


void EUTelClusteringProcessor::check (LCEvent * /* evt */) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelClusteringProcessor::end() {

  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;

  map< int, int >::iterator iter = _totClusterMap.begin();
  while ( iter != _totClusterMap.end() ) {

    streamlog_out ( MESSAGE2 ) << "Found " << iter->second << " clusters on detector " << iter->first << endl;
    ++iter;

  }
}




void EUTelClusteringProcessor::resetStatus(IMPL::TrackerRawDataImpl * status) {

    int i = 0;
    ShortVec::iterator iter = status->adcValues().begin();
    while ( iter != status->adcValues().end() )
    {
        if ( *iter == EUTELESCOPE::HITPIXEL )
        {
            *iter = EUTELESCOPE::GOODPIXEL;
        }
        else if ( *iter == EUTELESCOPE::MISSINGPIXEL )
        {
            *iter = EUTELESCOPE::GOODPIXEL;
        }
        ++iter;
        i++;
    } 
}




#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void EUTelClusteringProcessor::fillHistos (LCEvent * evt) {


    
  EUTelEventImpl * eutelEvent = static_cast<EUTelEventImpl*> (evt);
  EventType type              = eutelEvent->getEventType();

  if ( type == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return ;
  } else if ( type == kUNKNOWN ) {
    // if it is unknown we had already issued a warning to the user at
    // the beginning of the processEvent. If we get to here, it means
    // that the assumption that the event was a data event was
    // correct, so no harm to continue...
  }


  try {

    LCCollectionVec * pulseCollectionVec = dynamic_cast<LCCollectionVec*>  (evt->getCollection(_pulseCollectionName));
    CellIDDecoder<TrackerPulseImpl > cellDecoder(pulseCollectionVec);

    // I also need the noise collection too fill in the SNR histograms
    LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
    LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
    CellIDDecoder<TrackerDataImpl > noiseDecoder(noiseCollectionVec);

    vector<unsigned short> eventCounterVec( _noOfDetector, 0 );

    for ( int iPulse = _initialPulseCollectionSize; iPulse < pulseCollectionVec->getNumberOfElements(); iPulse++ ) {
      TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl*> ( pulseCollectionVec->getElementAt(iPulse) );
      ClusterType        type  = static_cast<ClusterType> ( static_cast<int> ( cellDecoder(pulse)["type"] ));
      SparsePixelType    pixelType = static_cast<SparsePixelType> (0);
      EUTelVirtualCluster * cluster;

      if ( type == kEUTelFFClusterImpl ) {
        cluster = new EUTelFFClusterImpl ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
      } else if (type == kEUTelBrickedClusterImpl) { //!HACK TAKI
        cluster = new EUTelBrickedClusterImpl ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
      } else if ( type == kEUTelDFFClusterImpl ) {
        cluster = new EUTelDFFClusterImpl ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
      } else if ( type == kEUTelSparseClusterImpl ) {
        // knowing that is a sparse cluster is not enough we need also
        // to know also the sparse pixel type. This information is
        // available in the "original_zsdata" collection. Let's get it!

        //TODO: (Phillip Hamnett) - This is inefficient, must look into a better way of getting the sparse cluster type
        LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
        TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
        CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
        pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));
        if ( pixelType == kEUTelGenericSparsePixel ) {
          cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel > ( static_cast<TrackerDataImpl*> ( pulse->getTrackerData() ) );
        } else {
          streamlog_out ( ERROR4 ) <<  "Unknown sparse cluster type. Sorry for quitting" << endl;
          throw UnknownDataTypeException("Cluster type unknown");
        }
      } else {
        streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
        throw UnknownDataTypeException("Cluster type unknown");
      }

      int detectorID = cluster->getDetectorID();
      //if this is an excluded sensor go to the next element
      bool foundexcludedsensor = false;
      for(size_t i = 0; i < _ExcludedPlanes.size(); ++i)
        {
          if(_ExcludedPlanes[i] == detectorID)
            {
              foundexcludedsensor = true;
	      break;
            }
        }
      if(foundexcludedsensor)
        continue;
      // increment of one unit the event counter for this plane
      eventCounterVec[ _ancillaryIndexMap[ detectorID] ]++;

      // plot the cluster total charge
      (dynamic_cast<AIDA::IHistogram1D*> (_clusterSignalHistos[detectorID]))->fill(cluster->getTotalCharge());
      streamlog_out ( DEBUG4 ) <<  " cluster charge " << cluster->getTotalCharge() << endl;
      // get the cluster size in X and Y separately and plot it:
      int xSize, ySize;
      cluster->getClusterSize(xSize,ySize);
      (dynamic_cast<AIDA::IHistogram1D*> (_clusterSizeXHistos[detectorID]))->fill(xSize);
      (dynamic_cast<AIDA::IHistogram1D*> (_clusterSizeYHistos[detectorID]))->fill(ySize);

      // plot the seed charge
      (dynamic_cast<AIDA::IHistogram1D*> (_seedSignalHistos[detectorID]))->fill(cluster->getSeedCharge());

      vector<float > charges = cluster->getClusterCharge(_clusterSpectraNVector);
      for ( unsigned int i = 0; i < charges.size() ; i++ ) {
        (dynamic_cast<AIDA::IHistogram1D*> (_clusterSignal_NHistos[_clusterSpectraNVector[i]][detectorID]))
          ->fill(charges[i]);
      }

      vector<int >::iterator iter = _clusterSpectraNxNVector.begin();
      while ( iter != _clusterSpectraNxNVector.end() ) {
        (dynamic_cast<AIDA::IHistogram1D*> (_clusterSignal_NxNHistos[*iter][detectorID]))
          ->fill(cluster->getClusterCharge((*iter), (*iter)));
        ++iter;
      }

      int xSeed, ySeed;
      cluster->getCenterCoord(xSeed, ySeed);
      (dynamic_cast<AIDA::IHistogram2D*> (_hitMapHistos[detectorID]))
        ->fill(static_cast<double >(xSeed), static_cast<double >(ySeed), 1.);


      // fill the noise related histograms

      // get the noise TrackerDataImpl corresponding to the detector
      // under analysis and the status matrix as well
      TrackerDataImpl    * noiseMatrix  = dynamic_cast<TrackerDataImpl *>    (noiseCollectionVec->getElementAt(  _ancillaryIndexMap[ detectorID ]) );
      TrackerRawDataImpl * statusMatrix = dynamic_cast<TrackerRawDataImpl *> (statusCollectionVec->getElementAt( _ancillaryIndexMap[ detectorID ]) );

      // prepare also a MatrixDecoder for this matrix
      EUTelMatrixDecoder noiseMatrixDecoder(noiseDecoder, noiseMatrix);
      int minX, minY, maxX, maxY;
      minX = 0;
      minY = 0;

      // this sensorID can be either a reference plane or a DUT, do it
      // differently...
      if ( _layerIndexMap.find( detectorID ) != _layerIndexMap.end() ){
        // this is a reference plane
        maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ detectorID ] ) - 1;
        maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ detectorID ] ) - 1;
      } else if ( _dutLayerIndexMap.find( detectorID ) != _dutLayerIndexMap.end() ) {
        // ok it is a DUT plane
        maxX = _siPlanesLayerLayout->getDUTSensitiveNpixelX() - 1;
        maxY = _siPlanesLayerLayout->getDUTSensitiveNpixelY() - 1;
      } else {
        // this is not a reference plane neither a DUT... what's that?
//        throw  InvalidGeometryException ("Unknown sensorID " + to_string( detectorID ));
        streamlog_out( ERROR5 ) << "Unknown sensorID " << detectorID << ", perhaps your GEAR file is incomplete." << endl;
        continue;
      }

      vector<float > noiseValues;
      //! not sure if this is 100% correct for bricked clusters here
      if ( type == kEUTelFFClusterImpl || type == kEUTelBrickedClusterImpl ) {
        for ( int yPixel = ySeed - ( _ffYClusterSize / 2 ); yPixel <= ySeed + ( _ffYClusterSize / 2 ); yPixel++ ) {
          for ( int xPixel = xSeed - ( _ffXClusterSize / 2 ); xPixel <= xSeed + ( _ffXClusterSize / 2 ); xPixel++ ) {
            // always check we are still within the sensor!!!
            if ( ( xPixel >= minX )  &&  ( xPixel <= maxX ) &&
                 ( yPixel >= minY )  &&  ( yPixel <= maxY ) ) {
              int index = noiseMatrixDecoder.getIndexFromXY(xPixel, yPixel);
              // the corresponding position in the status matrix has to be HITPIXEL
              bool isHit      = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::HITPIXEL );
              bool isBad      = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::BADPIXEL );
              bool isMissing  = ( statusMatrix->getADCValues()[index] == EUTELESCOPE::MISSINGPIXEL );
              if ( !isMissing && !isBad && isHit ) {
                noiseValues.push_back( noiseMatrix->getChargeValues()[index] );
              } else {
                noiseValues.push_back( 0. );
              }
            } else {
              noiseValues.push_back( 0. );
            }
          }
        }
      } else if ( type == kEUTelSparseClusterImpl ) {
        if ( pixelType == kEUTelGenericSparsePixel ) {
          auto_ptr<EUTelGenericSparsePixel> pixel ( new EUTelGenericSparsePixel );
          // this recasting is due to have access to sparse cluster
          // specific methods.
          EUTelSparseClusterImpl<EUTelGenericSparsePixel>* recasted =
            dynamic_cast<EUTelSparseClusterImpl<EUTelGenericSparsePixel>* > (cluster);
          for ( unsigned int iPixel = 0 ; iPixel < recasted->size() ; iPixel++ ) {
            recasted->getSparsePixelAt( iPixel , pixel.get() );
            int index = noiseMatrixDecoder.getIndexFromXY( pixel->getXCoord(), pixel->getYCoord() );
            noiseValues.push_back( noiseMatrix->getChargeValues() [ index ] );

          }
        }
      }

      bool fillSNRSwitch = true;
      if(type != kEUTelDFFClusterImpl)
        try {
          cluster->setNoiseValues(noiseValues);
        } catch ( IncompatibleDataSetException& e ) {
          streamlog_out ( WARNING2 ) << e.what() << endl
                                     << "Continuing without filling the noise histograms" << endl;
          fillSNRSwitch = false;
        }

      if(type == kEUTelDFFClusterImpl)
        fillSNRSwitch = false;

      if ( fillSNRSwitch ) {

        AIDA::IHistogram1D * histo;
        histo = dynamic_cast<AIDA::IHistogram1D* > ( _clusterNoiseHistos[detectorID] );
        if ( histo ) {
          histo->fill( cluster->getClusterNoise() );
        }

        histo = dynamic_cast<AIDA::IHistogram1D* > ( _clusterSNRHistos[detectorID] );
        if ( histo ) {
          histo->fill( cluster->getClusterSNR() );
        }

        histo = dynamic_cast<AIDA::IHistogram1D * > ( _seedSNRHistos[detectorID] );
        if ( histo ) {
          histo->fill( cluster->getSeedSNR() );
        }

        AIDA::IHistogram2D *histo2d = dynamic_cast<AIDA::IHistogram2D * > ( _cluster_vs_seedSNRHistos[detectorID] );
        if ( histo2d ) {
          histo2d->fill( cluster->getSeedSNR(), cluster->getClusterSNR() );
        }

        vector<int >::iterator iter = _clusterSpectraNxNVector.begin();
        while ( iter != _clusterSpectraNxNVector.end() ) {
          histo = dynamic_cast<AIDA::IHistogram1D*> ( _clusterSNR_NxNHistos[*iter][detectorID] ) ;
          if ( histo ) {
            histo->fill(cluster->getClusterSNR( (*iter), (*iter) ));
          }
          ++iter;
        }

        vector<float > snrs = cluster->getClusterSNR(_clusterSpectraNVector);
        for ( unsigned int i = 0; i < snrs.size() ; i++ ) {
          histo = dynamic_cast<AIDA::IHistogram1D * > ( _clusterSNR_NHistos[_clusterSpectraNVector[i]][detectorID] ) ;
          if ( histo ) {
            histo->fill( snrs[i] );
          }
        }
      }

      //delete cluster;
    }

    // fill the event multiplicity here
//     string tempHistoName;
//     for ( int iDetector = 0; iDetector < _noOfDetector; iDetector++ ) {
//       AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D *> ( _eventMultiplicityHistos[_orderedSensorIDVec.at( iDetector)] );
//       if ( histo ) {
//         histo->fill( eventCounterVec[iDetector] );
//       }
//     }
  } catch (lcio::DataNotAvailableException& e) {
    return;
  }

}

#endif

#ifdef MARLIN_USE_AIDA
void EUTelClusteringProcessor::bookHistos() {

  // histograms are grouped in loops and detectors
  streamlog_out ( DEBUG5 )  << "Booking histograms " << endl;
  auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
  EUTelHistogramInfo    * histoInfo;
  bool                    isHistoManagerAvailable;

  try {
    isHistoManagerAvailable = histoMgr->init();
  } catch ( ios::failure& e) {
    streamlog_out ( WARNING2 ) << "I/O problem with " << _histoInfoFileName << "\n"
                               << "Continuing without histogram manager"  << endl;
    isHistoManagerAvailable = false;
  } catch ( ParseException& e ) {
    streamlog_out ( WARNING2 ) << e.what() << "\n"
                               << "Continuing without histogram manager" << endl;
    isHistoManagerAvailable = false;
  }

  // define our histogram names
  std::string _clusterSignalHistoName      = "clusterSignal";
  std::string _clusterSizeXHistoName       = "clusterSizeX";
  std::string _clusterSizeYHistoName       = "clusterSizeY";
  std::string _seedSignalHistoName         = "seedSignal";
  std::string _hitMapHistoName             = "hitMap";
  std::string _seedSNRHistoName            = "seedSNR";
  std::string _clusterNoiseHistoName       = "clusterNoise";
  std::string _clusterSNRHistoName         = "clusterSNR";
  std::string _cluster_vs_seedSNRHistoName = "cluster_vs_seedSNR";
  std::string _eventMultiplicityHistoName  = "eventMultiplicity";


  string tempHistoName;
  string basePath;
  for (size_t iDetector = 0; iDetector < _sensorIDVec.size(); iDetector++) {

    int sensorID = _sensorIDVec.at( iDetector );

    // the min and max information are taken from GEAR.
    int minX, minY, maxX, maxY;
    minX = 0;
    minY = 0;

    // this sensorID can be either a reference plane or a DUT, do it
    // differently...
    if ( _layerIndexMap.find( sensorID ) != _layerIndexMap.end() ){
      // this is a reference plane
      maxX = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[ sensorID ] ) - 1;
      maxY = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[ sensorID ] ) - 1;
    } else if ( _dutLayerIndexMap.find( sensorID )  != _dutLayerIndexMap.end() ) {
      // ok it is a DUT plane
      maxX = _siPlanesLayerLayout->getDUTSensitiveNpixelX() - 1;
      maxY = _siPlanesLayerLayout->getDUTSensitiveNpixelY() - 1;
    } else {
      // this is not a reference plane neither a DUT... what's that?
//      throw  InvalidGeometryException ("Unknown sensorID " + to_string( sensorID ));
      streamlog_out( ERROR5 ) << "Unknown sensorID " << sensorID << ", perhaps your GEAR file is incomplete." << endl;
      continue;
    }

    basePath = "detector_" + to_string( sensorID );
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath.append("/");

    int    clusterNBin  = 1000;
    double clusterMin   = 0.;
    double clusterMax   = 1000.;
    string clusterTitle = "Cluster spectrum with all pixels";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _clusterSignalHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        clusterNBin = histoInfo->_xBin;
        clusterMin  = histoInfo->_xMin;
        clusterMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) clusterTitle = histoInfo->_title;
      }
    }

    int    clusterXNBin  = 1000;
    double clusterXMin   = 0.;
    double clusterXMax   = 1000.;
    string clusterXTitle = "Cluster X spectrum with all pixels";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _clusterSizeXHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        clusterXNBin = histoInfo->_xBin;
        clusterXMin  = histoInfo->_xMin;
        clusterXMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) clusterXTitle = histoInfo->_title;
      }
    }


    int    clusterYNBin  = 1000;
    double clusterYMin   = 0.;
    double clusterYMax   = 1000.;
    string clusterYTitle = "Cluster Y spectrum with all pixels";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _clusterSizeYHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        clusterYNBin = histoInfo->_xBin;
        clusterYMin  = histoInfo->_xMin;
        clusterYMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) clusterYTitle = histoInfo->_title;
      }
    }

    int    clusterSNRNBin  = 300;
    double clusterSNRMin   = 0.;
    double clusterSNRMax   = 200;
    string clusterSNRTitle = "Cluster SNR";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _clusterSNRHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        clusterSNRNBin = histoInfo->_xBin;
        clusterSNRMin  = histoInfo->_xMin;
        clusterSNRMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) clusterSNRTitle = histoInfo->_title;
      }
    }

    int    cluster_vs_seedSNRNBin_X  = 30;
    double cluster_vs_seedSNRMin_X   = 0.;
    double cluster_vs_seedSNRMax_X   = 100;
    int    cluster_vs_seedSNRNBin_Y  = 30;
    double cluster_vs_seedSNRMin_Y   = 0.;
    double cluster_vs_seedSNRMax_Y   = 100;

    string cluster_vs_seedSNRTitle = "Cluster vs Seed SNR";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _cluster_vs_seedSNRHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        cluster_vs_seedSNRNBin_X = histoInfo->_xBin;
        cluster_vs_seedSNRMin_X  = histoInfo->_xMin;
        cluster_vs_seedSNRMax_X = histoInfo->_xMax;
        cluster_vs_seedSNRNBin_Y = histoInfo->_yBin;
        cluster_vs_seedSNRMin_Y  = histoInfo->_yMin;
        cluster_vs_seedSNRMax_Y = histoInfo->_yMax;

        if ( histoInfo->_title != "" ) cluster_vs_seedSNRTitle = histoInfo->_title;
      }
    }


    // cluster signal
    tempHistoName = _clusterSignalHistoName + "_d" + to_string( sensorID );
    _clusterSignalHistos.insert(make_pair(sensorID, 
					  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
												    clusterNBin,clusterMin,clusterMax)));
    _clusterSignalHistos[sensorID]->setTitle(clusterTitle.c_str());

    // cluster signal along X
    tempHistoName = _clusterSizeXHistoName + "_d" + to_string( sensorID );
    _clusterSizeXHistos.insert(make_pair(sensorID, 
					AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
												  clusterXNBin,clusterXMin,clusterXMax)));
    _clusterSizeXHistos[sensorID]->setTitle(clusterXTitle.c_str());

    // cluster signal along Y
    tempHistoName = _clusterSizeYHistoName + "_d" + to_string( sensorID );
    _clusterSizeYHistos.insert(make_pair(sensorID, 
					AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
												  clusterYNBin,clusterYMin,clusterYMax)));
    _clusterSizeYHistos[sensorID]->setTitle(clusterYTitle.c_str());



    // cluster SNR
    tempHistoName = _clusterSNRHistoName + "_d" + to_string( sensorID );
    _clusterSNRHistos.insert( make_pair(sensorID, 
			      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
											clusterSNRNBin, clusterSNRMin, clusterSNRMax)));
    _clusterSNRHistos[sensorID]->setTitle(clusterSNRTitle.c_str());

    
   // cluster vs seed SNR
    tempHistoName = _cluster_vs_seedSNRHistoName + "_d" + to_string( sensorID );
    AIDA::IHistogram2D * cluster_vs_seedSNRHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                cluster_vs_seedSNRNBin_X,cluster_vs_seedSNRMin_X,cluster_vs_seedSNRMax_X,
                                                                cluster_vs_seedSNRNBin_Y,cluster_vs_seedSNRMin_Y,cluster_vs_seedSNRMax_Y);
    _cluster_vs_seedSNRHistos.insert( make_pair(sensorID, cluster_vs_seedSNRHisto) ) ;
    cluster_vs_seedSNRHisto->setTitle(cluster_vs_seedSNRTitle);

    
    vector<int >::iterator iter = _clusterSpectraNVector.begin();
    while ( iter != _clusterSpectraNVector.end() ) {
      // create 'outer' map if not existing for this sprectrum 'N'
      _clusterSignal_NHistos.insert( make_pair( *iter, std::map<int,AIDA::IBaseHistogram*>()));
      _clusterSNR_NHistos.insert( make_pair( *iter, std::map<int,AIDA::IBaseHistogram*>()));
      
      // this is for the signal
      tempHistoName = _clusterSignalHistoName + to_string( *iter ) + "_d" + to_string( sensorID );
      AIDA::IHistogram1D * clusterSignalNHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterNBin, clusterMin, clusterMax);
      string tempTitle = "Cluster spectrum with the " + to_string( *iter ) + " most significant pixels ";
      clusterSignalNHisto->setTitle(tempTitle.c_str());
      _clusterSignal_NHistos.at(*iter).insert(make_pair(sensorID, clusterSignalNHisto) );


      // this is for the SNR
      tempHistoName = _clusterSNRHistoName + to_string( *iter ) + "_d" + to_string( sensorID );
      AIDA::IHistogram1D * clusterSNRNHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterSNRNBin, clusterSNRMin, clusterSNRMax);
      tempTitle = "Cluster SNR with the " + to_string(*iter ) + " most significant pixels";
      clusterSNRNHisto->setTitle(tempTitle.c_str());
      _clusterSNR_NHistos.at(*iter).insert( make_pair( sensorID, clusterSNRNHisto ) );

      ++iter;
    } // while _clusterSpectraNVector

    iter = _clusterSpectraNxNVector.begin();
    while ( iter != _clusterSpectraNxNVector.end() ) {
      // create 'outer' map if not existing for this sprectrum 'NxN'
      _clusterSignal_NxNHistos.insert( make_pair( *iter, std::map<int,AIDA::IBaseHistogram*>()));
      _clusterSNR_NxNHistos.insert( make_pair( *iter, std::map<int,AIDA::IBaseHistogram*>()));

      // first the signal
      tempHistoName = _clusterSignalHistoName + to_string( *iter ) + "x"
        + to_string( *iter ) + "_d" + to_string( sensorID );
      AIDA::IHistogram1D * clusterSignalNxNHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterNBin, clusterMin, clusterMax);
      string tempTitle = "Cluster spectrum with " + to_string( *iter ) + " by " +  to_string( *iter ) + " pixels ";
      clusterSignalNxNHisto->setTitle(tempTitle.c_str());
      _clusterSignal_NxNHistos.at(*iter).insert(make_pair(sensorID, clusterSignalNxNHisto) );
      
      // then the SNR
      tempHistoName = _clusterSNRHistoName + to_string( *iter ) + "x"
        + to_string( *iter ) + "_d" + to_string( sensorID );
      AIDA::IHistogram1D * clusterSNRNxNHisto =
        AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                  clusterSNRNBin, clusterSNRMin, clusterSNRMax);
      tempTitle = "SNR with " + to_string( *iter ) + " by " + to_string( *iter ) + " pixels ";
      clusterSNRNxNHisto->setTitle(tempTitle.c_str());
      _clusterSNR_NxNHistos.at(*iter).insert(make_pair(sensorID, clusterSNRNxNHisto) );

      ++iter;
    } // while _clusterSpectraNxNVector

    tempHistoName = _seedSignalHistoName + "_d" + to_string( sensorID );
    int    seedNBin  = 500;
    double seedMin   = 0.;
    double seedMax   = 500.;
    string seedTitle = "Seed pixel spectrum";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _seedSignalHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        seedNBin = histoInfo->_xBin;
        seedMin  = histoInfo->_xMin;
        seedMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) seedTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * seedSignalHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                seedNBin, seedMin, seedMax);
    _seedSignalHistos.insert(make_pair(sensorID, seedSignalHisto));
    seedSignalHisto->setTitle(seedTitle.c_str());

    tempHistoName = _seedSNRHistoName + "_d" + to_string( sensorID ) ;
    int    seedSNRNBin  =  300;
    double seedSNRMin   =    0.;
    double seedSNRMax   =  200.;
    string seedSNRTitle = "Seed SNR";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _seedSNRHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        seedSNRNBin = histoInfo->_xBin;
        seedSNRMin  = histoInfo->_xMin;
        seedSNRMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) seedSNRTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * seedSNRHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                seedSNRNBin, seedSNRMin, seedSNRMax);
    _seedSNRHistos.insert( make_pair(sensorID, seedSNRHisto ));
    seedSNRHisto->setTitle(seedSNRTitle.c_str());

    tempHistoName = _clusterNoiseHistoName + "_d" + to_string( sensorID );
    int    clusterNoiseNBin  =  300;
    double clusterNoiseMin   =    0.;
    double clusterNoiseMax   =  200.;
    string clusterNoiseTitle = "Cluster noise";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo( _clusterNoiseHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        clusterNoiseNBin = histoInfo->_xBin;
        clusterNoiseMin  = histoInfo->_xMin;
        clusterNoiseMax  = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) clusterNoiseTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * clusterNoiseHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                clusterNoiseNBin, clusterNoiseMin, clusterNoiseMax);
    _clusterNoiseHistos.insert( make_pair(sensorID, clusterNoiseHisto ));
    clusterNoiseHisto->setTitle(clusterNoiseTitle.c_str());

    tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
    int     xBin = maxX - minX + 1;
    double  xMin = static_cast<double >( minX ) - 0.5 ;
    double  xMax = static_cast<double >( maxX ) + 0.5;
    int     yBin = maxY - minY + 1;
    double  yMin = static_cast<double >( minY ) - 0.5;
    double  yMax = static_cast<double >( maxY ) + 0.5;
    AIDA::IHistogram2D * hitMapHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                xBin, xMin, xMax,yBin, yMin, yMax);
    _hitMapHistos.insert(make_pair(sensorID, hitMapHisto));
    hitMapHisto->setTitle("Hit map");

    tempHistoName = _eventMultiplicityHistoName + "_d" + to_string( sensorID );
    int     eventMultiNBin  = 60;
    double  eventMultiMin   =  0.;
    double  eventMultiMax   = 60.;
    string  eventMultiTitle = "Event multiplicity";
    if ( isHistoManagerAvailable ) {
      histoInfo = histoMgr->getHistogramInfo(  _eventMultiplicityHistoName );
      if ( histoInfo ) {
        streamlog_out ( DEBUG2 ) << (* histoInfo ) << endl;
        eventMultiNBin  = histoInfo->_xBin;
        eventMultiMin   = histoInfo->_xMin;
        eventMultiMax   = histoInfo->_xMax;
        if ( histoInfo->_title != "" ) eventMultiTitle = histoInfo->_title;
      }
    }
    AIDA::IHistogram1D * eventMultiHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                eventMultiNBin, eventMultiMin, eventMultiMax);
    _eventMultiplicityHistos.insert( make_pair(sensorID, eventMultiHisto) );
    eventMultiHisto->setTitle( eventMultiTitle.c_str() );

  }
  streamlog_out ( DEBUG5 )  << "end of Booking histograms " << endl;
}

#endif
#endif // USE_GEAR
