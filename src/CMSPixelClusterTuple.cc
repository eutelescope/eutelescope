
// Author: U. Grundler, National Taiwan University

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined(USE_GEAR) && ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

// eutelescope inlcudes
#include "CMSPixelClusterTuple.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseData2Impl.h"
#include "EUTelSparseCluster2Impl.h"

// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/ITreeFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/ITupleFactory.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>


#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static members mainly used to name histograms
std::string CMSPixelClusterTuple::_ClusterTupleName  = "CMSCluster";


CMSPixelClusterTuple::CMSPixelClusterTuple() : Processor("CMSPixelClusterTuple") {

  // modify processor description
  _description = "Prepare n-tuple with cluster results" ;


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACKERPULSE,
                           "InputCollectionName" ,
                           "Name of the input cluster collection"  ,
                           _clusterColName ,
                           std::string("cluster_pixel") ) ;

  // other processor parameters:


  registerProcessorParameter ("DebugEventCount",
                              "Print out every DebugEnevtCount event",
                              _debugCount,  static_cast < int > (100));

  registerProcessorParameter ("MissingValue",
                              "Value used for missing measurements",
                              _missingValue,  static_cast < double > (-100.));

}


void CMSPixelClusterTuple::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;


  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR5> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  // Read geometry information from GEAR

  message<MESSAGE5> ( log() << "Reading telescope geometry description from GEAR ") ;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));


// Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

// Read position in Z (for sorting)

  _planeSort = new int[_nTelPlanes];
  _planePosition   = new double[_nTelPlanes];

  for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
    {
      _planePosition[ipl]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
      _planeSort[ipl]=ipl;
    }

  // Binary sorting

  bool sorted;
  do{
    sorted=false;
    for(int iz=0; iz<_nTelPlanes-1 ; iz++)
      if(_planePosition[iz]>_planePosition[iz+1])
        {
          double _posZ = _planePosition[iz];
          _planePosition[iz] = _planePosition[iz+1];
          _planePosition[iz+1] = _posZ;

          int _idZ = _planeSort[iz];
          _planeSort[iz] = _planeSort[iz+1];
          _planeSort[iz+1] = _idZ;

          sorted=true;
        }

  }while(sorted);

// Book local geometry arrays

  _planeID         = new int[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];

// Fill remaining layer parameters

  for(int iz=0; iz < _nTelPlanes ; iz++)
    {
      int ipl=_planeSort[iz];

      double resolution;

      _planeID[iz]=_siPlanesLayerLayout->getID(ipl);
      resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);

      _isActive[iz] = (resolution > 0);

    }


  // Print out geometry information

  message<MESSAGE5> ( log() << "Telescope configuration with " << _nTelPlanes << " planes" );


  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      stringstream ss ;
      if(_isActive[ipl])
         ss << "Active  plane" ;
      else
         ss << "Passive plane" ;

      ss << "  ID = " << _planeID[ipl]
         << "  at Z [mm] = " << _planePosition[ipl];

      message<MESSAGE5> ( log() << ss.str() );
    }


// Book histograms
  bookHistos();


}

void CMSPixelClusterTuple::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  _runNr = runHeader->getRunNumber();

  message<MESSAGE5> ( log() << "Processing run header " << _nRun
                     << ", run nr " << _runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  message<MESSAGE5> ( log() << detectorName << " : " << detectorDescription ) ;

  int nDet = subDets->size();

  if(nDet)message<MESSAGE5> ( log() << nDet << " subdetectors defined :" );
  stringstream ss;
  for(int idet=0;idet<nDet;idet++)  message<MESSAGE5> (log()  << idet+1 << " : " << subDets->at(idet) );


}

void CMSPixelClusterTuple::processEvent( LCEvent * event ) {

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG5> ( "EORE found: nothing else to do." );
    return;
  }

  bool debug = ( _debugCount>0 && _nEvt%_debugCount == 0);

  _nEvt ++ ;
  _evtNr = event->getEventNumber();

  if(debug)message<DEBUG5> ( log() << "Processing record " << _nEvt << " == event " << _evtNr );

  LCCollectionVec* col;
  try {
     col = static_cast<LCCollectionVec*> (event->getCollection( _clusterColName ));
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection " << _clusterColName << "from event " << event->getEventNumber() << " in run " << event->getRunNumber() << endl;
    return;
  }

  CellIDDecoder<TrackerPulseImpl> cellDecoder( col );

  // Loop over clusters in input collection
  int _noOfDetector = _nTelPlanes;
  vector<unsigned short> eventCounterVec( _noOfDetector, 0 );
  vector<unsigned short> noOfClusters( _noOfDetector, 0);

  int nCluster = col->getNumberOfElements()  ;

  if(debug)message<DEBUG5> ( log() << "Total of " << nCluster << " clusters in input collection " );

  // _ClusterTuple->fill(0,_runNr);
  // _ClusterTuple->fill(1,_evtNr);
  // _ClusterTuple->fill(2,nCluster);

  // AIDA::ITuple *aCluster = _ClusterTuple->getTuple(3);
  for(int iclu=0; iclu< nCluster ; iclu++) {//loop over clusters
     TrackerPulseImpl * pulse = dynamic_cast<TrackerPulseImpl*>( col->getElementAt(iclu) );
     TrackerDataImpl * cluster = dynamic_cast<TrackerDataImpl*>( pulse->getTrackerData() );
     
     int sensorID = ( static_cast<int> ( cellDecoder(pulse)["sensorID"] ));
     if(debug)message<DEBUG5> ( log() << "SensorID " << sensorID );
     ClusterType type = static_cast<ClusterType> ( static_cast<int> (cellDecoder( pulse )["type"]) );
     eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > *pixelCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel >(cluster);	

     int size = pixelCluster->size();
     eventCounterVec[ sensorID ] += size;
     noOfClusters[ sensorID ]++;
      
     int xSize, ySize;
     pixelCluster->getClusterSize(xSize,ySize);

     if(debug)message<DEBUG5> ( log() << "Cluster size (x,y) " << size << " (" << xSize << "," << ySize << ")" );

     // Fill n-tuple
     // aCluster->fill(0,sensorID);

     // aCluster->fill(1,size);
     // aCluster->fill(2,xSize);
     // aCluster->fill(3,ySize);

     int xSeed, ySeed;
     pixelCluster->getCenterCoord(xSeed, ySeed);
     // aCluster->fill(4,xSeed);
     // aCluster->fill(5,ySeed);
     
     // aCluster->fill(6,pixelCluster->getTotalCharge());

     // AIDA::ITuple *aPixel = aCluster->getTuple(7);
     for(int iPix=0; iPix<size; iPix++) {
        EUTelSimpleSparsePixel Pixel;
        pixelCluster->getSparsePixelAt(iPix, &Pixel);

        int icol = 0;
        _ClusterTuple->fill(icol++,_nEvt);
        _ClusterTuple->fill(icol++,_runNr);
        _ClusterTuple->fill(icol++,_evtNr);

        _ClusterTuple->fill(icol++,iclu);
        _ClusterTuple->fill(icol++,sensorID);

        _ClusterTuple->fill(icol++,size);
        _ClusterTuple->fill(icol++,xSize);
        _ClusterTuple->fill(icol++,ySize);

        _ClusterTuple->fill(icol++,xSeed);
        _ClusterTuple->fill(icol++,ySeed);
        _ClusterTuple->fill(icol++,pixelCluster->getTotalCharge());

        _ClusterTuple->fill(icol++,iPix);
        _ClusterTuple->fill(icol++,Pixel.getXCoord());
        _ClusterTuple->fill(icol++,Pixel.getYCoord());
        _ClusterTuple->fill(icol++,Pixel.getSignal());

        _ClusterTuple->addRow();
        // aPixel->fill(0,Pixel.getXCoord());
        // aPixel->fill(1,Pixel.getYCoord());
        // aPixel->fill(2,Pixel.getSignal());
        // aPixel->addRow();
     }
     
     // aCluster->addRow();

  }//loop over clusters
  // _ClusterTuple->addRow();

  return;
}



void CMSPixelClusterTuple::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void CMSPixelClusterTuple::end(){

  //   std::cout << "CMSPixelClusterTuple::end()  " << name()
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;


  message<MESSAGE5> ( log() << "N-tuple with "
                     << _ClusterTuple->rows() << " rows created" );


}



void CMSPixelClusterTuple::bookHistos()
{


  message<MESSAGE5> ( log() << "Booking fit n-tuple \n" );

  std::vector<std::string> _columnNames;
  std::vector<std::string> _columnType;

  _columnNames.push_back("Event");
  _columnType.push_back("int");

  _columnNames.push_back("RunNr");
  _columnType.push_back("int");

  _columnNames.push_back("EvtNr");
  _columnType.push_back("int");

  _columnNames.push_back("clusterID");
  _columnType.push_back("int");
  _columnNames.push_back("sensorID");
  _columnType.push_back("int");
  _columnNames.push_back("clusterSize");
  _columnType.push_back("int");
  _columnNames.push_back("clusterSizeX");
  _columnType.push_back("int");
  _columnNames.push_back("clusterSizeY");
  _columnType.push_back("int");
  _columnNames.push_back("cPosX");
  _columnType.push_back("int");
  _columnNames.push_back("cPosY");
  _columnType.push_back("int");
  _columnNames.push_back("clusterCharge");
  _columnType.push_back("float");
  _columnNames.push_back("pixelID");
  _columnType.push_back("int");
  _columnNames.push_back("pPosX");
  _columnType.push_back("short");
  _columnNames.push_back("pPosY");
  _columnType.push_back("short");
  _columnNames.push_back("pixelCharge");
  _columnType.push_back("float");


  _ClusterTuple=AIDAProcessor::tupleFactory(this)->create(_ClusterTupleName, _ClusterTupleName, _columnNames, _columnType, "");

  //string columnString = "int RunNr = -1, EvtNr = -1; ITuple cluster = {int sensorID=-1, size=-1, widthX=-1, widthY=-1; double posX=-1., posY=-1, charge=-1.; ITuple pixel = {int pixelX, pixelY; double pixelCharge=-1.}";
  //_ClusterTuple= AIDAProcessor::tupleFactory(this)->create(_ClusterTupleName, _ClusterTupleName, columnString, "");


  message<DEBUG5> ( log() << "Booking completed \n\n");

  return;
}

#endif // GEAR && AIDA
