// eutelescope inlcudes
#include "EUTelAPIXTbTrackTuple.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelSparseDataImpl.h"

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
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

#include <string>
#include <vector>
#include <map>
#include <stdlib.h>

#include <UTIL/CellIDEncoder.h>

//TbTrack include
#include <algorithm>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <TVector3.h>

using namespace std;
using namespace marlin ;
using namespace eutelescope;


EUTelAPIXTbTrackTuple::EUTelAPIXTbTrackTuple() : Processor("EUTelAPIXTbTrackTuple") {
  // modify processor description
  _description = "Prepare tbtrack style n-tuple with track fit results" ;
	
  registerInputCollection( LCIO::TRACK,
			   "InputCollectionName",
			   "Name of the input Track collection"  ,
			   _inputTrackColName ,
			   std::string("testfittracks") ) ;

  registerProcessorParameter("InputDutPulseCollectionName",
		    "Name of the input PulseCollection to provide infos about clusters and raw-pixels",
		    _inputDutPulseCollectionName,
		    std::string("cluster_apix"));

  registerProcessorParameter("InputTelPulseCollectionName",
		    "Name of the input PulseCollection to provide infos about clusters and raw-pixels",
		    _inputTelPulseCollectionName,
		    std::string("cluster_m26"));
		    
  registerInputCollection( LCIO::TRACKERHIT,
		"InputTrackerHitCollectionName" ,
		"Name of the plane-wide hit-data hit collection"  ,
		_inputTrackerHitColName ,
		"alignedHit" ) ;
  
  
  registerProcessorParameter ("TelZsColName",
		     "Tel zero surpressed data colection name",
 		     _telZsColName, string ("zsdata_m26"));
  registerProcessorParameter ("DutZsColName",
		     "DUT zero surpressed data colection name",
 		     _dutZsColName, string ("zsdata_apix"));

  registerOptionalParameter("AlignmentCollectionNames", "Names of alignment collections, should be in same order as application", _alignColNames, std::vector<std::string>());

  
  registerOptionalParameter("ReadInClusterBased","Do not readin collection by collection but with the correct mapping between clusters, zsdata and hits", _clusterBased, true);

  // other processor parameters:
  registerProcessorParameter ("OutputPath",
			      "Path/File where root-file should be stored",
			      _path2file,  static_cast < std::string > ("NTuple.root"));
  registerProcessorParameter ("DoScales",
			      "If true assume alignment corrections uses scales, if false assume full rotations",
			      _doScales, true);
}


void EUTelAPIXTbTrackTuple::init() {
  // usually a good idea to
  printParameters();

  _nRun = 0 ;
  _nEvt = 0 ;
  _foundAllign = false;

  message<MESSAGE> ( log() << "Initializing " );
	
  if ( Global::GEAR == NULL ) {
    message<ERROR> ( "GearMgr is not available. Bailing." );
    exit(-1);
  }
	
  message<MESSAGE> ( log() << "Reading telescope geometry description from GEAR ") ;
  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  //z pos of hit -> index
  

  //Mark excluded
  //Mark 
  // Prepare TTree/TFiles
  prepareTree();
  invertGear();
  streamlog_out(MESSAGE3) << "HitCollection is " << _inputTrackerHitColName << endl;
  message<MESSAGE> ( log() << "End of Init" );
}

void EUTelAPIXTbTrackTuple::processRunHeader( LCRunHeader* runHeader) {
  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );
  _nRun++ ;
	
  // Decode and print out Run Header information - just a check
  _runNr = runHeader->getRunNumber();
  message<MESSAGE> ( log() << "Processing run header " << _nRun
		     << ", run nr " << _runNr );
}

void EUTelAPIXTbTrackTuple::processEvent( LCEvent * event ) {
  _nEvt ++ ;
  _evtNr = event->getEventNumber();
  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG> ( "EORE found: nothing else to do." );
    return;
  }

  if( not _foundAllign ) { readAlignment(event); } //Sets _foundAlign to true if found.
  if( not _foundAllign ) {
    streamlog_out  ( ERROR ) << "Have not found the needed alignment collections, will skip this event ( " 
			     << event->getEventNumber() << " )." << endl; 
    return;
  }

  //Clear all event info containers
  clear();

  if( readHits( _inputTrackerHitColName, event ) != 0) { return; }
  /* --- Dump zs info, skip if event if collection not found --- */
  //if( readZsHits( _telZsColName ,  event) != 0) { return; }
   if (_clusterBased) {
   	if( readZsHits( _dutZsColName , event) != 0) { return; }
   } else {
  	if( readZsHitsFromClusters( _dutZsColName , event) != 0) { return; }
   }
  /* --- Dump fitted track params --- */
  if( readTracks(event) != 0) { return; }
  /* ---- Check cluster collections ---*/
  if( readClusters( _inputDutPulseCollectionName , event) != 0 ) { return; } 
  //if( readClusters( _inputTelPulseCollectionName , event) != 0 ) { return; } 
 
  if (_clusterBased) setClusterIdInHits();

  /* Filling tree */
  _zstree->Fill();
  _eutracks->Fill();
  _clutree->Fill();
   if (_clusterBased) _euhits->Fill();
}

void EUTelAPIXTbTrackTuple::end(){
  message<MESSAGE> ( log() << "N-tuple with " << _zstree->GetEntries() << " entries written to" << _path2file.c_str() <<",");
  _file->Write();
}

void EUTelAPIXTbTrackTuple::setClusterIdInHits() {
  size_t endsize = _hitSensorId->size();
  if (_hitClusterId->size() == 0) {
    std::map<IMPL::TrackerDataImpl*, int>::iterator it;
    
    for (size_t ihit =0; ihit < _hitPointerToCluster->size(); ihit++) {
      it = (*_cluPointer).find(_hitPointerToCluster->at(ihit));
      if (it == (*_cluPointer).end()) { 
	//cout << "not in there" << endl; 
	_hitClusterId->push_back(-1);
      } else {
	int cluIndex = (*it).second;
	_hitClusterId->push_back(_cluClusterId->at(cluIndex));
      }
    }
  } else {
    streamlog_out(WARNING) << "Method setClusterIdInHits() may only be called once!" << endl;	
  }
  if (endsize != _hitClusterId->size() ) cout << "ClusterIndex variable is not correct!" << endl;
}

int EUTelAPIXTbTrackTuple::readHits( std::string hitColName, LCEvent* event ) {
  LCCollection *hitCollection = NULL;
  int nHit = -1;
  try {
    hitCollection = event->getCollection( hitColName ); 
  } catch (lcio::DataNotAvailableException& e) {
    //message<MESSAGE> ( log() << "Cluster collection " << colName <<  " not available on event " << event->getEventNumber() << "." );
    return(1);
  }
  nHit = hitCollection->getNumberOfElements();
  _nHits = nHit;
  for(int ihit=0; ihit< nHit ; ihit++) {
    TrackerHit * meshit = dynamic_cast<TrackerHitImpl*>( hitCollection->getElementAt(ihit) ) ;
    const double * pos = meshit->getPosition();	
    LCObjectVec clusterVec = (meshit->getRawHits());
    TrackerDataImpl * trackerData = dynamic_cast < TrackerDataImpl * > ( clusterVec[0] );
    
    int iden = -1;
    for(int plane = 0; plane < _siPlanesLayerLayout->getNLayers(); plane++){
      if( fabs( pos[2] - _siPlanesLayerLayout->getSensitivePositionZ( plane ) ) < 2.0){
	iden = _siPlanesLayerLayout->getID( plane );
      }
    }
    if(iden == -1){ continue; }
    double nominalZpos(0.0);
    for(int plane = 0; plane < _siPlanesLayerLayout->getNLayers(); plane++){
      if( fabs( pos[2] - _siPlanesLayerLayout->getSensitivePositionZ( plane ) ) < 2.0){
	iden = _siPlanesLayerLayout->getID( plane );
	nominalZpos = _siPlanesLayerLayout->getSensitivePositionZ(plane)
	  + 0.5 * _siPlanesLayerLayout->getSensitiveThickness(plane);
      }
    }

    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    reverseAlign(x,y,z,iden, nominalZpos);
    
    _hitXPos->push_back(x);
    _hitYPos->push_back(y);
    _hitZPos->push_back(z);
    _hitSensorId->push_back(iden);
    _hitPointerToCluster->push_back(trackerData);
  }
  
  return(0);
}


int EUTelAPIXTbTrackTuple::readClusters( std::string colName, LCEvent* event ){
  LCCollection* clusterCollectionVec = NULL;
  try {
    clusterCollectionVec = event->getCollection( colName );
  } catch (lcio::DataNotAvailableException& e) {
    //message<MESSAGE> ( log() << "Cluster collection " << colName <<  " not available on event " << event->getEventNumber() << "." );
    return(1);
  }
  CellIDDecoder<TrackerPulseImpl> cellDecoder(clusterCollectionVec);

  for(int cluster = 0 ; cluster < clusterCollectionVec->getNumberOfElements(); cluster++){
    TrackerPulseImpl* pulseFrame   = dynamic_cast<TrackerPulseImpl*> ( clusterCollectionVec->getElementAt(cluster) );
    TrackerDataImpl* clusterFrame = dynamic_cast<TrackerDataImpl*> ( pulseFrame->getTrackerData() );
    int              sensorID     = ( static_cast<int> ( cellDecoder(pulseFrame)["sensorID"] ));
    ClusterType       type         = static_cast<ClusterType> ( static_cast<int> (cellDecoder( pulseFrame )["type"]) );

    int size(0), sizeX(0), sizeY(0);
    int posX(0), posY(0);
    int charge(0), clusterID(0);
    if(type == kEUTelAPIXClusterImpl){
      eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);	
      size = apixCluster->size();
      charge = apixCluster->getTotalCharge();
      apixCluster->getClusterSize(sizeX, sizeY);
      apixCluster->getCenterCoord(posX, posY);
      clusterID = apixCluster->getClusterID();
      delete apixCluster;
    } else if (  type == kEUTelSparseClusterImpl or type == kEUTelM26ClusterImpl or type == kEUTelDFFClusterImpl){
      eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > * telCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel >(clusterFrame);
      size = telCluster->size();
      charge = telCluster->getTotalCharge();
      telCluster->getClusterSize(sizeX, sizeY);
      telCluster->getCenterCoord(posX, posY);
      clusterID = telCluster->getClusterID();
      delete telCluster;
    } else {
      message<WARNING> ( log() << "Unknown cluster type: " << type );
      sensorID = -1;
    }
    if(sensorID > 0){
      _clusize->push_back( size );
      _clusizeX->push_back( sizeX );
      _clusizeY->push_back( sizeY );
      _cluposX->push_back( posX );
      _cluposY->push_back( posY );
      _clucharge->push_back( charge );
      _cluSensorId->push_back( sensorID );
      _cluClusterId->push_back( clusterID );
      (*_cluPointer)[clusterFrame] = _cluClusterId->size() -1; // To enable a mapping between Link and the index in the map (size -1) is the last element of the vector
    }
  }
  return(0);
}

int EUTelAPIXTbTrackTuple::readTracks(LCEvent* event){
  LCCollection* trackCol=NULL;
  try {
    trackCol = event->getCollection( _inputTrackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    //message<MESSAGE> ( log() << "Track collection " << _inputColName << " not available on event " << event->getEventNumber() << "." );
    return(1);
  }

  int nTrackParams=0;
  for(int itrack=0; itrack< trackCol->getNumberOfElements(); itrack++) {
    lcio::Track* fittrack = dynamic_cast<lcio::Track*>( trackCol->getElementAt(itrack) ) ;
    std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();
    double chi2 = fittrack->getChi2();
    double ndof = fittrack->getNdf();
    double dxdz = fittrack->getOmega();
    double dydz = fittrack->getPhi();

    for(unsigned int ihit=0; ihit< trackhits.size(); ihit++) {
      TrackerHit* fittedHit = trackhits.at(ihit);
      const double * pos = fittedHit->getPosition();
      //Type >= 32 for fitted hits
      if( fittedHit->getType() < 32) { continue; }
      int iden = -1;
      double nominalZpos(0.0);
      for(int plane = 0; plane < _siPlanesLayerLayout->getNLayers(); plane++){
	if( fabs( pos[2] - _siPlanesLayerLayout->getSensitivePositionZ( plane ) ) < 2.0){
	  iden = _siPlanesLayerLayout->getID( plane );
	  nominalZpos = _siPlanesLayerLayout->getSensitivePositionZ(plane)
	    + 0.5 * _siPlanesLayerLayout->getSensitiveThickness(plane);
	}
      }
      if(iden == -1){ continue; }
      nTrackParams++;
      double x = pos[0];
      double y = pos[1];
      double z = pos[2];
      reverseAlign(x,y,z,iden, nominalZpos);
      //eutrack tree
      _xPos->push_back(x);
      _yPos->push_back(y);
      _dxdz->push_back(dxdz);
      _dydz->push_back(dydz);
      _trackIden->push_back(iden);
      _trackNum->push_back(itrack);
      _chi2->push_back(chi2);
      _ndof->push_back(ndof);
    }
  }
  _nTrackParams = nTrackParams;
  return(0);
}

int EUTelAPIXTbTrackTuple::readZsHitsFromClusters( std::string colName, LCEvent* event){
  LCCollectionVec* zsInputCollectionVec= NULL;
  try{
    zsInputCollectionVec  = dynamic_cast < LCCollectionVec* > (event->getCollection( colName ));
  } catch (DataNotAvailableException& e){
    //streamlog_out ( ERROR ) << "Could not find collection " << colName << " on event " << event->getEventNumber() << "." << endl;
    return(1);
  }
  CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
  for ( unsigned int cluster = 0 ; cluster < zsInputCollectionVec->size(); cluster++ ) {
    TrackerDataImpl* zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( cluster ) );
    int sensorID = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
    int clusterID = static_cast<int> ( cellDecoder (zsData) ["clusterID"] );
    
    _nPixHits++;
    if (sensorID > 9  ) {
      auto_ptr<EUTelSparseDataImpl<EUTelAPIXSparsePixel> > apixData(new EUTelSparseDataImpl<EUTelAPIXSparsePixel> ( zsData ));
      EUTelAPIXSparsePixel apixPixel;
      for ( unsigned int iHit = 0; iHit < apixData->size(); iHit++ ) {
		apixData->getSparsePixelAt( iHit, &apixPixel);
		p_iden->push_back( sensorID );
		p_chip->push_back( apixPixel.getChip() );
		p_row->push_back( apixPixel.getYCoord() );
		p_col->push_back( apixPixel.getXCoord() );
		p_tot->push_back( apixPixel.getSignal() );
		p_lv1->push_back( apixPixel.getTime() );
		p_clusterId->push_back(clusterID);
      }
		
   
    } else if ( sensorID > 0 && sensorID < 10 ) {
      auto_ptr<EUTelSparseDataImpl<EUTelSimpleSparsePixel> > telData(new EUTelSparseDataImpl<EUTelSimpleSparsePixel> ( zsData ));
      EUTelSimpleSparsePixel telPixel;
      for ( unsigned int iHit = 0; iHit < telData->size(); iHit++ ) {
	telData->getSparsePixelAt( iHit, &telPixel);
	p_iden->push_back( sensorID );
	p_chip->push_back( 0 );
	p_row->push_back( telPixel.getYCoord() );
	p_col->push_back( telPixel.getXCoord() );
	p_tot->push_back( telPixel.getSignal() );
	p_lv1->push_back( 0 );
	p_clusterId->push_back(clusterID);
      }
    }
  }
  if (p_clusterId->size() != p_row->size() ) { streamlog_out(ERROR) << "ClusterID Vector size is not consistent!(" << p_clusterId->size() << " and " << p_row->size() << endl; }
    return(0);
}

int EUTelAPIXTbTrackTuple::readZsHits( std::string colName, LCEvent* event){
  LCCollectionVec* zsInputCollectionVec= NULL;
  try{
    zsInputCollectionVec  = dynamic_cast < LCCollectionVec* > (event->getCollection( colName ));
  } catch (DataNotAvailableException& e){
    //streamlog_out ( ERROR ) << "Could not find collection " << colName << " on event " << event->getEventNumber() << "." << endl;
    return(1);
  }
  CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
  for ( unsigned int plane = 0 ; plane < zsInputCollectionVec->size(); plane++ ) {
    TrackerDataImpl* zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( plane ) );
    SparsePixelType type = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
    int sensorID = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
    if (type == kEUTelAPIXSparsePixel  ) {
      auto_ptr<EUTelSparseDataImpl<EUTelAPIXSparsePixel> > apixData(new EUTelSparseDataImpl<EUTelAPIXSparsePixel> ( zsData ));
      EUTelAPIXSparsePixel apixPixel;
      for ( unsigned int iHit = 0; iHit < apixData->size(); iHit++ ) {
	apixData->getSparsePixelAt( iHit, &apixPixel);
	_nPixHits++;
	p_iden->push_back( sensorID );
	p_chip->push_back( apixPixel.getChip() );
	p_row->push_back( apixPixel.getYCoord() );
	p_col->push_back( apixPixel.getXCoord() );
	p_tot->push_back( apixPixel.getSignal() );
	p_lv1->push_back( apixPixel.getTime() );
      }
    } else if ( type == kEUTelSimpleSparsePixel ) {
      auto_ptr<EUTelSparseDataImpl<EUTelSimpleSparsePixel> > telData(new EUTelSparseDataImpl<EUTelSimpleSparsePixel> ( zsData ));
      EUTelSimpleSparsePixel telPixel;
      for ( unsigned int iHit = 0; iHit < telData->size(); iHit++ ) {
	telData->getSparsePixelAt( iHit, &telPixel);
	p_iden->push_back( sensorID );
	p_chip->push_back( 0 );
	p_row->push_back( telPixel.getYCoord() );
	p_col->push_back( telPixel.getXCoord() );
	p_tot->push_back( telPixel.getSignal() );
	p_lv1->push_back( 0 );
      }
    }
  }
  return(0);
}

gsl_matrix* EUTelAPIXTbTrackTuple::invertLU(int dim, gsl_matrix* matrix){
  gsl_permutation* perm = gsl_permutation_alloc(dim);
  int s = 0;
  gsl_matrix * inverse = gsl_matrix_alloc(dim,dim);
  gsl_linalg_LU_decomp(matrix, perm, &s);
  gsl_linalg_LU_invert(matrix, perm, inverse);
  gsl_permutation_free(perm);
  //gsl_matrix_fprintf(stdout,inverse,"%f");
  return (inverse);
}

void EUTelAPIXTbTrackTuple::invertAlignment(EUTelAlignmentConstant * alignment){
  int iden = alignment->getSensorID();
  
  //Rotations
  gsl_matrix * alignM = gsl_matrix_alloc(3,3);
  //Fill alignment matrix
  if(_doScales){
    gsl_matrix_set( alignM, 0, 0, (1.0 + alignment->getAlpha() ));
    gsl_matrix_set( alignM, 1, 1, (1.0 + alignment->getBeta()  ));
    gsl_matrix_set( alignM, 2, 2, 1.0 );
    gsl_matrix_set( alignM, 0, 1, alignment->getGamma() );
    gsl_matrix_set( alignM, 1, 0, -1 * alignment->getGamma() );
    gsl_matrix_set( alignM, 0, 2, 0.0 );
    gsl_matrix_set( alignM, 2, 0, 0.0 );
    gsl_matrix_set( alignM, 1, 2, 0.0 );
    gsl_matrix_set( alignM, 2, 1, 0.0 );
  } else {
    gsl_matrix_set( alignM, 0, 0, 1.0);
    gsl_matrix_set( alignM, 1, 1, 1.0);
    gsl_matrix_set( alignM, 2, 2, 1.0 );
    gsl_matrix_set( alignM, 0, 1, alignment->getGamma() );
    gsl_matrix_set( alignM, 1, 0, -1 * alignment->getGamma());
    gsl_matrix_set( alignM, 0, 2, alignment->getBeta() );
    gsl_matrix_set( alignM, 2, 0, -1.0 * alignment->getBeta());
    gsl_matrix_set( alignM, 1, 2, alignment->getAlpha() );
    gsl_matrix_set( alignM, 2, 1, -1.0 * alignment->getAlpha());
  }
  message<MESSAGE> ( log() << "Inverting alignment matrix for iden" << iden  ) ;
  gsl_matrix * inverse = invertLU(3, alignM);
  gsl_matrix_free(alignM);
  _alignRot[iden].push_back(inverse);

  //Shifts
  std::vector<double> shifts;
  shifts.push_back( alignment->getXOffset()  );
  shifts.push_back( alignment->getYOffset() );
  shifts.push_back( alignment->getZOffset() );
  _alignShift[iden].push_back(shifts);
  streamlog_out(MESSAGE) << "Iden: " << iden << endl
			 << "X-shift: "<< alignment->getXOffset() << endl
			 << "Y-shift: "<< alignment->getYOffset() << endl;
}

void EUTelAPIXTbTrackTuple::invertGear(){
  for ( int layerIndex = 0 ; layerIndex < _siPlanesParameters->getSiPlanesNumber() ; ++layerIndex ) {
    int iden = _siPlanesLayerLayout->getID( layerIndex );
    //Rotations
    gsl_matrix * gearM = gsl_matrix_alloc(2,2);
    gsl_matrix_set( gearM, 0, 0, _siPlanesLayerLayout->getSensitiveRotation1(layerIndex));
    gsl_matrix_set( gearM, 0, 1, _siPlanesLayerLayout->getSensitiveRotation2(layerIndex));
    gsl_matrix_set( gearM, 1, 0, _siPlanesLayerLayout->getSensitiveRotation3(layerIndex));
    gsl_matrix_set( gearM, 1, 1, _siPlanesLayerLayout->getSensitiveRotation4(layerIndex));

    //Ugly as sin, but following the suggestions of EUTelHitMaker
    double xSign(0.0), ySign(0.0);
    if( (gsl_matrix_get(gearM, 0 ,0) < -0.7) or (gsl_matrix_get(gearM,0,1) < -0.7)){
      xSign = -1.0;
    } else if( (gsl_matrix_get(gearM, 0 ,0) > 0.7) or (gsl_matrix_get(gearM,0,1) > 0.7)){
      xSign = 1.0;
    }
    if( (gsl_matrix_get(gearM, 1 ,0) < -0.7) or (gsl_matrix_get(gearM,1,1) < -0.7)){
      ySign = -1.0;
    } else if( (gsl_matrix_get(gearM, 1 ,0) > 0.7) or (gsl_matrix_get(gearM,1,1) > 0.7)){
      ySign = 1.0;
    }
  
    message<MESSAGE> ( log() << "Inverting gear matrix for iden" << iden  ) ;
    gsl_matrix * inverse = invertLU(2, gearM);
    gsl_matrix_free(gearM);
    _gearRot[iden] = inverse;
    
    //shifts
    vector<double> shifts;
    shifts.push_back(-1 * ( _siPlanesLayerLayout->getSensitivePositionX(layerIndex) 
			    - xSign * 0.5 * _siPlanesLayerLayout->getSensitiveSizeX(layerIndex)));
    shifts.push_back(-1 * ( _siPlanesLayerLayout->getSensitivePositionY(layerIndex)
			    - ySign * 0.5 * _siPlanesLayerLayout->getSensitiveSizeY(layerIndex)));
    shifts.push_back(-1 * _siPlanesLayerLayout->getSensitivePositionZ(layerIndex));
    _gearShift[iden] = shifts;

    //Pitches
    _gearPitch[iden] = make_pair( _siPlanesLayerLayout->getSensitivePitchX(layerIndex),
				  _siPlanesLayerLayout->getSensitivePitchY(layerIndex));

    std::vector<double> gearRotations(3, 0.0);
    //Convert angles to radians
    double conv = 3.1415926/180.0;
    gearRotations.at(0) = conv * _siPlanesLayerLayout->getLayerRotationXY(layerIndex);
    gearRotations.at(1) = conv * _siPlanesLayerLayout->getLayerRotationZX(layerIndex);
    gearRotations.at(2) = conv * _siPlanesLayerLayout->getLayerRotationZY(layerIndex);
    cout << "Plane iden: " << iden << endl;
    cout << "gearRotationsXY = " << gearRotations.at(0) << endl
	 << "gearRotationsZX = " << gearRotations.at(1) << endl
	 << "gearRotationsZY = " << gearRotations.at(2) << endl;
    _gearEulerRot[iden] = gearRotations;
  }
}

void EUTelAPIXTbTrackTuple::clear(){
  /* Clear zsdata */
  p_col->clear();
  p_row->clear();
  p_tot->clear();
  p_iden->clear();
  p_lv1->clear();
  p_chip->clear();
  p_clusterId->clear();
  _nPixHits = 0;
  /* Clear hittrack */
  _xPos->clear();
  _yPos->clear();
  _dxdz->clear();
  _dydz->clear();
  _trackNum->clear();
  _trackIden->clear();
  _chi2->clear();
  _ndof->clear();
  //Clear cluster
  _clusize->clear();
  _clusizeX->clear();
  _clusizeY->clear();
  _cluposX->clear();	 
  _cluposY->clear();	 
  _clucharge->clear();
  _cluSensorId->clear();
  _cluClusterId->clear();
  _cluPointer->clear(); // Clear the map as well!
  //Clear hits
  _hitXPos->clear();
  _hitYPos->clear();
  _hitZPos->clear();
  _hitClusterId->clear();
  _hitSensorId->clear();
  _hitPointerToCluster->clear();
}

void EUTelAPIXTbTrackTuple::reverseAlign(double& x, double& y, double &z, int iden, double nomZpos){
  // Apply alignment translations 
  double xTemp(0.0),yTemp(0.0), zTemp(0.0);
  for( size_t trans = 0; trans < _alignShift[iden].size(); trans++){
    //if( _alignShift.find(iden) != _alignShift.end() ){
    x += _alignShift[iden].at(trans).at(0);
    y += _alignShift[iden].at(trans).at(1);
    z += _alignShift[iden].at(trans).at(2);
    double zShift(z);
    if(not _doScales){ zShift = z - nomZpos;}
    //Apply alignment rotations
    gsl_matrix* m = _alignRot[iden].at(trans);
    xTemp = x * gsl_matrix_get(m,0,0) + y * gsl_matrix_get(m,0,1) + zShift * gsl_matrix_get(m,0,2); 
    yTemp = x * gsl_matrix_get(m,1,0) + y * gsl_matrix_get(m,1,1) + zShift * gsl_matrix_get(m,1,2); 
    zTemp = x * gsl_matrix_get(m,2,0) + y * gsl_matrix_get(m,2,1) + z * gsl_matrix_get(m,2,2); 
    x = xTemp; y = yTemp; z = zTemp;
    //Do not need z's from here, I guess
  }

  //Apply GEAR Euler translations
  zTemp = z - nomZpos;
  TVector3 RotatedSensorHit( x, y, zTemp);
  std::vector<double>& rots = _gearEulerRot[iden];
  if( TMath::Abs(rots.at(0)) > 1e-6 ){
    RotatedSensorHit.RotateZ( -1.0 * rots.at(0) ); // in XY
  }
  if( TMath::Abs(rots.at(1)) > 1e-6 ){
    RotatedSensorHit.RotateY( -1.0 * rots.at(1) ); // in ZX 
  }
  if( TMath::Abs(rots.at(2)) > 1e-6 ){
    RotatedSensorHit.RotateX( -1.0 * rots.at(2) ); // in ZY
  }
  x = RotatedSensorHit.X();
  y = RotatedSensorHit.Y();
  
  //Apply gear local trans
  if(_gearShift.find(iden) != _gearShift.end()){
    x += _gearShift[iden].at(0);
    y += _gearShift[iden].at(1);
  }  //Apply gear rot
  if(_gearRot.find(iden) != _gearRot.end()){
    gsl_matrix* m = _gearRot[iden];
    xTemp = x * gsl_matrix_get(m,0,0) + y * gsl_matrix_get(m,0,1);
    yTemp = x * gsl_matrix_get(m,1,0) + y * gsl_matrix_get(m,1,1);
    x = xTemp; y = yTemp;
  }
  // Change reference from sensor corner to center of pixel 0,0
  x -= 0.5 * _gearPitch[iden].first;
  y -= 0.5 * _gearPitch[iden].second;
}

void EUTelAPIXTbTrackTuple::readAlignment(LCEvent* event){
  //Check for alignment collections, if found invert and store, if not keep _foundAlign false so it will look in next event.
  _foundAllign = true;

  LCCollectionVec * alignmentCollectionVec;
  cout << "Trying to read " << _alignColNames.size() << " collections" << endl;
  for(size_t ii = 0; ii < _alignColNames.size(); ii++){
    size_t index = _alignColNames.size() - ii -1;
    cout << "Trying to read alignment collection " << _alignColNames.at(index) << endl;
    try{
      alignmentCollectionVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _alignColNames.at(index)));
      for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) {
	cout << "Inverting plane " << iPos << endl;
	invertAlignment( static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) ) );
      }
    } catch( DataNotAvailableException& e) {
      streamlog_out ( ERROR ) << "Could not find  alignment collections on event " << event->getEventNumber() << "." << endl;
      _foundAllign = false;
    }
  }
}

void EUTelAPIXTbTrackTuple::prepareTree(){
  _file = new TFile(_path2file.c_str(),"RECREATE");
  cout << "Writing to: " << _path2file.c_str() << endl;
  //Old school tbtrack tree

  _xPos = new vector<double>();	     
  _yPos = new vector<double>();	   
  _dxdz = new vector<double>();	   
  _dydz = new vector<double>();	   
  _trackIden  = new vector<int>();
  _trackNum = new vector<int>();
  _chi2 = new vector<double>();	   
  _ndof = new vector<double>();    
  
  p_col = new vector<int>();
  p_row = new vector<int>();
  p_tot = new vector<int>();
  p_iden = new vector<int>();
  p_lv1 = new vector<int>();
  p_chip = new vector<int>();
  p_clusterId = new vector<int>();
  
  
  _clusize = new vector<int>();
  _clusizeX = new vector<int>();
  _clusizeY = new vector<int>();
  _cluposX = new vector<int>();
  _cluposY = new vector<int>();
  _clucharge = new vector<int>();
  _cluSensorId = new vector<int>();
  _cluClusterId = new vector<int>();
  _cluPointer = new std::map<IMPL::TrackerDataImpl*, int>();
  
  _hitXPos = new vector<double>();
  _hitYPos = new vector<double>();
  _hitZPos = new vector<double>();
  _hitClusterId  = new vector<int>();
  _hitSensorId  = new vector<int>();
  _hitPointerToCluster = new vector<IMPL::TrackerDataImpl*>();
  
  _euhits = new TTree("euhits","euhits");
  _euhits->Branch("nHits", &_nHits);
  _euhits->Branch("xPos", &_hitXPos);
  _euhits->Branch("yPos", &_hitYPos);
  _euhits->Branch("zPos", &_hitZPos);
  _euhits->Branch("clusterId", &_hitClusterId);
  _euhits->Branch("sensorId", &_hitSensorId);
      
  _zstree = new TTree("zspix", "zspix");
  _zstree->Branch("nPixHits", &_nPixHits);
  _zstree->Branch("euEvt",    &_nEvt);
  _zstree->Branch("col",      &p_col);
  _zstree->Branch("row",      &p_row);
  _zstree->Branch("tot",      &p_tot);
  _zstree->Branch("lv1",      &p_lv1);
  _zstree->Branch("iden",     &p_iden);
  _zstree->Branch("chip",     &p_chip);
  _zstree->Branch("clusterId",&p_clusterId);
  //Tree for storing all track param info
  _eutracks = new TTree("eutracks", "eutracks");
  _eutracks->Branch("nTrackParams", &_nTrackParams);
  _eutracks->Branch("euEvt", &_nEvt);
  _eutracks->Branch("xPos", &_xPos);
  _eutracks->Branch("yPos", &_yPos);
  _eutracks->Branch("dxdz", &_dxdz);
  _eutracks->Branch("dydz", &_dydz);
  _eutracks->Branch("trackNum", &_trackNum);
  _eutracks->Branch("iden", &_trackIden);
  _eutracks->Branch("chi2", &_chi2);
  _eutracks->Branch("ndof", &_ndof);
  //TTree for cluster info
  _clutree = new TTree("euclusters", "euclusters");
  _clutree->Branch("euEvt", &_nEvt);
  _clutree->Branch("size", &_clusize);
  _clutree->Branch("sizeX", &_clusizeX);
  _clutree->Branch("sizeY", &_clusizeY);
  _clutree->Branch("posX", &_cluposX);
  _clutree->Branch("posY", &_cluposY);
  _clutree->Branch("charge", &_clucharge);
  _clutree->Branch("iden", &_cluSensorId);
  _clutree->Branch("ID", &_cluClusterId);
  
  _euhits->AddFriend(_clutree);
  _euhits->AddFriend(_zstree);
}
