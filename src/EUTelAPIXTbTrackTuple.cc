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
  
  
  registerProcessorParameter ("TelZsColName",
		     "Tel zero surpressed data colection name",
		     _telZsColName, string ("zsdata_m26"));
  registerProcessorParameter ("DutZsColName",
		     "DUT zero surpressed data colection name",
		     _dutZsColName, string ("zsdata_apix"));

  registerInputCollection (LCIO::LCGENERICOBJECT, "DutAlignmentConstantName",
                           "Alignment constant from the condition file",
                           _dutAlignmentCollectionName, string ("alignment"));
  registerInputCollection (LCIO::LCGENERICOBJECT, "TelAlignmentConstantName",
                           "Alignment constant from the condition file",
                           _telAlignmentCollectionName, string ("alignment"));

  // other processor parameters:
  registerProcessorParameter ("OutputPath",
			      "Path/File where root-file should be stored",
			      _path2file,  static_cast < std::string > ("NTuple.root"));
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
    
    if ( euEvent->getEventType() == kEORE ) 
    {
        message<DEBUG> ( "EORE found: nothing else to do." );
        return;
    }

    if( not _foundAllign ) 
    {
        readAlignment(event);  
        //Sets _foundAlign to true if found.
    }

    // now check once again if the alignment collection has been found
    if( not _foundAllign ) 
    {
        readAlignment(event);  //Sets _foundAlign to true if found.
        
        streamlog_out  ( ERROR ) << "Have not found the needed alignment collections, will skip this event ( " 
            << event->getEventNumber() << " )." << endl; 
        return;
    }

    //Clear all event info containers
    clear();

    /* --- Dump zs info, skip if event if collection not found --- */
    if( readZsHits( _telZsColName ,  event) != 0) 
    {
        streamlog_out  ( ERROR ) << "Have not found the needed telZsColName  " << _telZsColName << " collections, will skip this event ( " 
            << event->getEventNumber() << " )." << endl; 
        return; 
    }
  
    if( readZsHits( _dutZsColName , event) != 0) 
    {
        streamlog_out  ( ERROR ) << "Have not found the needed dutZsColName " << _dutZsColName << " collections, will skip this event ( " 
            << event->getEventNumber() << " )." << endl; 
        return;         
    }
    
    /* --- Dump fitted track params --- */
    if( readTracks(event) != 0) 
    {
        streamlog_out  ( ERROR ) << "Have not found any Tracks, will skip this event ( " 
            << event->getEventNumber() << " )." << endl; 
        return;
    }
    
    /* ---- Check cluster collections ---*/
    if( readClusters( _inputDutPulseCollectionName , event) != 0 ) 
    {
        streamlog_out  ( ERROR ) << "Have not found the needed dutPulseCollectionName " << _inputDutPulseCollectionName << " collections, will skip this event ( " 
            << event->getEventNumber() << " )." << endl; 
        return;
    } 
  
    if( readClusters( _inputTelPulseCollectionName , event) != 0 ) 
    {
        streamlog_out  ( ERROR ) << "Have not found the needed telPulseCollectionName " << _inputTelPulseCollectionName << " collections, will skip this event ( " 
            << event->getEventNumber() << " )." << endl; 
        return; 
    } 
  
    /* Filling tree */
    _zstree->Fill();
    _eutracks->Fill();
    _clutree->Fill();
    
}

void EUTelAPIXTbTrackTuple::end(){
  message<MESSAGE> ( log() << "N-tuple with " << _zstree->GetEntries() << " entries written to" << _path2file.c_str() <<",");
  _file->Write();
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
    int charge(0);
    if(type == kEUTelAPIXClusterImpl){
      eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);	
      size = apixCluster->size();
      charge = apixCluster->getTotalCharge();
      apixCluster->getClusterSize(sizeX, sizeY);
      apixCluster->getCenterCoord(posX, posY);
      delete apixCluster;
    } else if ( type == kEUTelSparseClusterImpl or type == kEUTelDFFClusterImpl){
      eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > * telCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel >(clusterFrame);
      size = telCluster->size();
      charge = telCluster->getTotalCharge();
      telCluster->getClusterSize(sizeX, sizeY);
      telCluster->getCenterCoord(posX, posY);
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
      _cluid->push_back( sensorID );
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
      for(int plane = 0; plane < _siPlanesLayerLayout->getNLayers(); plane++){
	if( fabs( pos[2] - _siPlanesLayerLayout->getSensitivePositionZ( plane ) ) < 2.0){
	  iden = _siPlanesLayerLayout->getID( plane );
	}
      }
      if(iden == -1){ continue; }
      nTrackParams++;
      double x = pos[0];
      double y = pos[1];
      double z = pos[2];
      reverseAlign(x,y,z,iden);
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
  gsl_matrix_fprintf(stdout,inverse,"%f");
  return (inverse);
}

void EUTelAPIXTbTrackTuple::invertAlignment(EUTelAlignmentConstant * alignment){
  int iden = alignment->getSensorID();
  
  //Rotations
  gsl_matrix * alignM = gsl_matrix_alloc(3,3);
  //Fill alignment matrix
  gsl_matrix_set( alignM, 0, 0, (1.0 + alignment->getAlpha() ));
  gsl_matrix_set( alignM, 1, 1, (1.0 + alignment->getBeta()  ));
  gsl_matrix_set( alignM, 2, 2, 1.0 );
  gsl_matrix_set( alignM, 0, 1, alignment->getGamma() );
  gsl_matrix_set( alignM, 1, 0, -1 * alignment->getGamma() );
  gsl_matrix_set( alignM, 0, 2, 0.0 );
  gsl_matrix_set( alignM, 2, 0, 0.0 );
  gsl_matrix_set( alignM, 1, 2, 0.0 );
  gsl_matrix_set( alignM, 2, 1, 0.0 );
  
  message<MESSAGE> ( log() << "Inverting alignment matrix for iden" << iden  ) ;
  gsl_matrix * inverse = invertLU(3, alignM);
  gsl_matrix_free(alignM);
  _alignRot[iden] = inverse;

  //Shifts
  std::vector<double> shifts;
  shifts.push_back( alignment->getXOffset()  );
  shifts.push_back( alignment->getYOffset() );
  shifts.push_back( alignment->getZOffset() );
  _alignShift[iden] = shifts;
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
  _cluid->clear();
}

void EUTelAPIXTbTrackTuple::reverseAlign(double& x, double& y, double &z, int iden){
  // Apply alignment translations 
  
  double xTemp(0.0),yTemp(0.0);
  if( _alignShift.find(iden) != _alignShift.end() ){
    x += _alignShift[iden].at(0);
    y += _alignShift[iden].at(1);
    z += _alignShift[iden].at(2);
  }  //Apply alignment rotations
  if( _alignRot.find(iden) != _alignRot.end() ){
    gsl_matrix* m = _alignRot[iden];
    xTemp = x * gsl_matrix_get(m,0,0) + y * gsl_matrix_get(m,0,1) + z * gsl_matrix_get(m,0,2); 
    yTemp = x * gsl_matrix_get(m,1,0) + y * gsl_matrix_get(m,1,1) + z * gsl_matrix_get(m,1,2); 
    x = xTemp; y = yTemp;
    //Do not need z's from here, I guess
  }
  //This seems to me to be a strange way of doing things, but following the method of EUTelHitMaker in reverse:
  //Apply gear trans
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

  LCCollectionVec * dutAlignmentCollectionVec = NULL; 
  LCCollectionVec * telAlignmentCollectionVec = NULL; 
  try{
    dutAlignmentCollectionVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _dutAlignmentCollectionName));
    telAlignmentCollectionVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _telAlignmentCollectionName));
  } catch( DataNotAvailableException& e) {
    //streamlog_out ( ERROR ) << "Could not find  alignment collections on event " << event->getEventNumber() << "." << endl;
    _foundAllign = false;
  }
  if(not _foundAllign){ return; }
  streamlog_out ( MESSAGE ) << "Found alignment collections on event " << event->getEventNumber() << "."  << endl;
  //The order of these two matters!
  //Read DUT alignment first, Unity rotations of the fixed telescopes are overwritten when tel alignment is read
  for ( size_t iPos = 0; iPos < dutAlignmentCollectionVec->size(); ++iPos ) {
    invertAlignment( static_cast< EUTelAlignmentConstant * > ( dutAlignmentCollectionVec->getElementAt( iPos ) ) );
  }
  for ( size_t iPos = 0; iPos < telAlignmentCollectionVec->size(); ++iPos ) {
    invertAlignment( static_cast< EUTelAlignmentConstant * > ( telAlignmentCollectionVec->getElementAt( iPos ) ) );
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
  
  _clusize = new vector<int>();
  _clusizeX = new vector<int>();
  _clusizeY = new vector<int>();
  _cluposX = new vector<int>();
  _cluposY = new vector<int>();
  _clucharge = new vector<int>();
  _cluid = new vector<int>();
  
  _zstree = new TTree("zspix", "zspix");
  _zstree->Branch("nPixHits", &_nPixHits);
  _zstree->Branch("euEvt",    &_nEvt);
  _zstree->Branch("col",      &p_col);
  _zstree->Branch("row",      &p_row);
  _zstree->Branch("tot",      &p_tot);
  _zstree->Branch("lv1",      &p_lv1);
  _zstree->Branch("iden",     &p_iden);
  _zstree->Branch("chip",     &p_chip);
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
  _clutree->Branch("iden", &_cluid);
}
