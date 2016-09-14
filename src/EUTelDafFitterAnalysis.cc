// based on EUTelDafFitter.cc

// eutelescope includes ".h"
#include "EUTelDafFitterAnalysis.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"
#include "EUTelReferenceHit.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <UTIL/CellIDEncoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/TrackerPulse.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

// ROOT includes
#include "TVector3.h"

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelDafFitter::EUTelDafFitter () : EUTelDafBase("EUTelDafFitterAnalysis"){
    //Child spesific params and description
  dafParams();
}

void EUTelDafFitter::dafParams(){
  _description = "This processor preforms track reconstruction and also outputs residuals of the track with respect to true hits. The tracks are as final track fit for analysis.";

  //Tracker system options
  registerOptionalParameter("AddToLCIO", "Should plots be made and filled?", _addToLCIO, static_cast<bool>(true));
  registerOptionalParameter("FitDuts","Set this to true if you want DUTs to be included in the track fit", _fitDuts, static_cast<bool>(false)); 
  //Track fitter options
	registerInputCollection(LCIO::TRACKERHIT, "TrueHitCollectionName", "Input of True Hit data", _trueHitCollectionName, string("true_hits"));
  registerOutputCollection(LCIO::TRACK,"TrackCollectionName", "Collection name for fitted tracks", _trackCollectionName, string("tracks"));
}

void EUTelDafFitter::dafInit() {
  if(_fitDuts){
    for( size_t ii = 0; ii< _system.planes.size(); ii++){
      if( find(_dutPlanes.begin(), _dutPlanes.end(), _system.planes.at(ii).getSensorID()) != _dutPlanes.end()){ 
	_system.planes.at(ii).include(); 
      }
    }
  }

	bookDiffHistos();
}

void EUTelDafFitter::dafEvent (LCEvent * event) {

	//read in the true hits
	try {

		_trueHitCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_trueHitCollectionName));
		streamlog_out(DEBUG4) << "_trueHitCollectionVec: " << _trueHitCollectionName.c_str() << " found\n";
	}
	catch (lcio::DataNotAvailableException& e) {

		streamlog_out(DEBUG4) << "_trueHitCollectionVec: " << _trueHitCollectionName.c_str() << " not found in event " << event->getEventNumber() << std::endl;
	}

  //Prepare track collection
  if(_addToLCIO){
    // Define output track and hit collections
    _fittrackvec = new LCCollectionVec(LCIO::TRACK);
    _fitpointvec = new LCCollectionVec(LCIO::TRACKERHIT);
    // Set flag for storing track hits in track collection
    LCFlagImpl flag(_fittrackvec->getFlag());
    flag.setBit( LCIO::TRBIT_HITS );
    _fittrackvec->setFlag(flag.getFlag());
  }
  
  //Check found tracks
  for(size_t ii = 0; ii < _system.getNtracks(); ii++ ){
    //run track fitte
    _nCandidates++;
    //Prepare track for DAF fit
    _system.fitPlanesInfoDaf(_system.tracks.at(ii));
    //Check resids, intime, angles
    if(not checkTrack( _system.tracks.at(ii))) { continue;};
    int inTimeHits = checkInTime(_system.tracks.at(ii));
    if( inTimeHits < _nDutHits) { continue;}
 
    //Fill plots
    if(_histogramSwitch){ 
      fillPlots( _system.tracks.at(ii) ); 
      fillDetailPlots( _system.tracks.at(ii) ); 
    }

		fillDiffHistos(_system.tracks.at(ii), event);

    //Dump to LCIO
    if( _addToLCIO) { addToLCIO(_system.tracks.at(ii), _fitpointvec); }
    _nTracks++;
  }

 //Add track collection
  if(_addToLCIO)
  { 
    event->addCollection(_fittrackvec,_trackCollectionName); 
    std::string  sfitpoints = "" ;  

    
    for(int i = 0; i<2000; i++) //TODO (Phillip Hamnett) Why is this hard coded to 1000? Ric changed to 2000.
    {
      sfitpoints = "fitpoints" + i;
      try
      {
       dynamic_cast < LCCollectionVec * > ( event->getCollection( sfitpoints ) )  ;
      }
      catch(...)
      {
        break;
      } 
    }
    event->addCollection(_fitpointvec, sfitpoints );
  }
}

void EUTelDafFitter::addToLCIO(daffitter::TrackCandidate<float,4>& track, LCCollectionVec *lcvec){
  TrackImpl * fittrack = new TrackImpl();
  // Impact parameters are useless and set to 0
  fittrack->setD0(0.);        // impact paramter of the track in (r-phi)
  fittrack->setZ0(0.);        // impact paramter of the track in (r-z)
  fittrack->setTanLambda(0.); // dip angle of the track at reference point

  daffitter::TrackEstimate<float,4>& est = track.estimates.at(0);

  //No good way of storing the track angles, so
  fittrack->setOmega( est.getXdz()); // Storing dxdz as Omega
  fittrack->setPhi( est.getYdz() );   // Storing dx/dy as phi

  fittrack->setChi2(track.chi2);
  fittrack->setNdf(int( round(track.ndof)) );
  // prepare an encoder for the hit collection to store properties
  CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING, lcvec);

  float refpoint[3];
  
  for(size_t plane = 0; plane < _system.planes.size(); plane++){
    daffitter::FitPlane<float>& pl = _system.planes.at(plane);
    daffitter::TrackEstimate<float,4>& estim = track.estimates.at( plane );
    TrackerHitImpl * fitpoint = new TrackerHitImpl();
    // encode and store sensorID
    int sensorID =  _system.planes.at(plane).getSensorID();
    idHitEncoder["sensorID"] = sensorID;
    // set the local/global bit flag property AND the FittedHit property for the hit
    idHitEncoder["properties"] = kHitInGlobalCoord | kFittedHit;
    double pos[3];
    pos[0]= estim.getX() / 1000.0;
    pos[1]= estim.getY() / 1000.0;
    pos[2]= pl.getMeasZ() / 1000.0;

    // overload z coordinate calculation -> important for proper sensor Identification by the hit coordinates based onthe refhit collection
    // if( fabs(pos[2] - getZfromRefHit(plane, sensorID, pos)) > 0.0002 ){
    //   streamlog_out(WARNING) << "Fitted measurement is not in the plane! SensorID " << idHitEncoder["sensorID"] << std::endl;
    //   pos[2] = getZfromRefHit(plane, sensorID, pos);    
    // }
    
    fitpoint->setPosition(pos);
    // Covariance matrix of the fitted position
    // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).
    float cov[TRKHITNCOVMATRIX];
    cov[0]= estim.cov(0,0);
    cov[1]= estim.cov(0,1);
    cov[2]= estim.cov(1,1);
    //Error of z position of fit is unclear to me, this would be a systematic alignment
    //error. Set to 0 along with all covariances.
    cov[3]=cov[4]=cov[5]=0.;
    fitpoint->setCovMatrix(cov);
    // store values
    idHitEncoder.setCellID( fitpoint );
    _fitpointvec->push_back(fitpoint);
    fittrack->addHit(fitpoint);

    streamlog_out(DEBUG3) << " hit : sensorID " << idHitEncoder["sensorID"] << " properties: " << idHitEncoder["properties"]  << std::endl;

    if(plane == 0){
      refpoint[0] = pos[0];
      refpoint[1] = pos[1];
      refpoint[2] = pos[2];
    }
    //Also add measurement point
    for(size_t mm = 0; mm < pl.meas.size(); mm++){
      if( track.weights.at(plane)(mm) < 0.5f){ continue; }
      if( pl.isExcluded()) { continue; }
      TrackerHitImpl* meashit = static_cast<TrackerHitImpl*> ( _hitCollection->getElementAt( pl.meas.at(mm).getIden()) );
      fittrack->addHit(meashit);
    }
  }
  fittrack->setReferencePoint(refpoint);
  _fittrackvec->addElement(fittrack);
}

double EUTelDafFitter::getZfromRefHit(int plane, int sensorID, double *pos){
         
  TVector3 lpoint( pos[0], pos[1], pos[2] );
  TVector3 lvector( 0., 0., 1. );
  TVector3 hitInPlane;
  TVector3 norm2Plane;
 
  //Name is misleading, is actually true if refHit is NOT set 
  if( ReferenceHitVecIsSet() ){
    hitInPlane.SetXYZ( geo::gGeometry().siPlaneXPosition(sensorID), geo::gGeometry().siPlaneYPosition(sensorID), geo::gGeometry().siPlaneZPosition(sensorID) );
    norm2Plane = geo::gGeometry().siPlaneNormal(sensorID);
  } else {
    EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(plane) ) ;
    hitInPlane.SetXYZ( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
    norm2Plane.SetXYZ( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
  } 
  
  TVector3 point( 1.,1.,1. );
          
  double linecoord_numenator   = norm2Plane.Dot(hitInPlane-lpoint);
  double linecoord_denumenator = norm2Plane.Dot(lvector);
  point = (linecoord_numenator/linecoord_denumenator)*lvector + lpoint;

  return point(2);
}

int EUTelProcessorTrueHitAnalysis::findPairIndex(float x, float y, int planeID) {

	CellIDDecoder<TrackerHitImpl> trueHitDecoder(_trueHitCollectionVec);

	if (_trueHitCollectionVec.size() == 1) return 0;
	else {

		int min_dist_index = -1;
		double min_distance = 0;
		for (int i = 0; i < _trueHitCollectionVec->getNumberOfElements(); i++) {

			if (static_cast<int>(trueHitDecoder(trueHit)["sensorID"]) == planeID) {

				if (min_dist_index == -1) {

					min_distance = sqrt(pow(x - _trueHitCollectionVec->getElementAt(i)->getPosition()[0], 2) + pow(y - _trueHitCollectionVec->getElementAt(i)->getPosition()[1], 2));
					min_dist_index = i;
				}
				else {

					double current_distance = sqrt(pow(x - _trueHitCollectionVec->getElementAt(i)->getPosition()[0], 2) + pow(y - _trueHitCollectionVec->getElementAt(i)->getPosition()[1], 2));

					if (current_distance < min_distance) {

						min_dist_index = i;
						min_distance = current_distance;
					}
				}
			}
		}

		return min_dist_index;
	}
}

void EUTelDafFitterAnalysis::fillDiffHistos(daffitter::TrackCandidate<float,4>& track, LCEvent* event) {

	/*_aidaHistoMap["chi2"]->fill( track.chi2);
	_aidaHistoMap["logchi2"]->fill( std::log10(track.chi2));
	_aidaHistoMap["ndof"]->fill( track.ndof);
	_aidaHistoMap["chi2overndof"]->fill( track.chi2 / track.ndof);*/

	//Fill plots per plane
	for( size_t i = 0; i < _system.planes.size() ; i++){

		daffitter::FitPlane<float>& plane = _system.planes.at(i);
		/*char iden[4];
		sprintf(iden, "%d", plane.getSensorID());*/
		string bname = static_cast< string >("pl") + to_string(plane.getSensorID()) + "_trueHit_";
		//Plot resids, angles for all hits with > 50% includion in track.
		//This should be one measurement per track

		daffitter::TrackEstimate<float,4>& estim = track.estimates.at(i);
		int pairIndex = findPairIndex(estim.getX(), estim.getY(), i);

		if (pairIndex == -1) {

			streamlog_out(WARNING2) << "found an unpaired track hit estimate at event " << event->getEventNumber() << " on plane " << plane.getSensorID() << std::endl;
			continue;
		}

		_aidaHistoMap[bname+"residualX"]->fill((estim.getX()-_trueHitCollectionVec->getElementAt(pairIndex)->getPosition()[0])*1e-3);
		_aidaHistoMap[bname+"residualY"]->fill((estim.getY()-_trueHitCollectionVec->getElementAt(pairIndex)->getPosition()[1])*1e-3);

		_aidaHistoMapProf1D[bname+"residualdXvsX"]->fill(estim.getX(), estim.getX()-_trueHitCollectionVec->getElementAt(pairIndex)->getPosition()[0]);
		_aidaHistoMapProf1D[bname+"residualdYvsX"]->fill(estim.getX(), estim.getY()-_trueHitCollectionVec->getElementAt(pairIndex)->getPosition()[1]);
		_aidaHistoMapProf1D[bname+"residualdXvsY"]->fill(estim.getY(), estim.getX()-_trueHitCollectionVec->getElementAt(pairIndex)->getPosition()[0]);
		_aidaHistoMapProf1D[bname+"residualdYvsY"]->fill(estim.getY(), estim.getY()-_trueHitCollectionVec->getElementAt(pairIndex)->getPosition()[1]);

		//finished with true hit, so remove it
		_trueHitCollectionVec->removeElementAt(pairIndex);

		/*for(size_t w = 0; w < plane.meas.size(); w++){

			if( track.weights.at(ii)(w) < 0.5f ) {  continue; }

			daffitter::Measurement<float>& meas = plane.meas.at(w);
			//Resids 
			_aidaHistoMap[bname + "residualX"]->fill( (estim.getX() - meas.getX())*1e-3 );
			_aidaHistoMap[bname + "residualY"]->fill( (estim.getY() - meas.getY())*1e-3 );*/

			//Resids 
			/*_aidaHistoMapProf1D[bname + "residualdXvsX"]->fill(estim.getX(), estim.getX() - meas.getX() );
			_aidaHistoMapProf1D[bname + "residualdYvsX"]->fill(estim.getX(), estim.getY() - meas.getY() );
			_aidaHistoMapProf1D[bname + "residualdXvsY"]->fill(estim.getY(), estim.getX() - meas.getX() );
			_aidaHistoMapProf1D[bname + "residualdYvsY"]->fill(estim.getY(), estim.getY() - meas.getY() );
			_aidaHistoMapProf1D[bname + "residualdZvsX"]->fill(estim.getX(), plane.getMeasZ() - meas.getZ()  );
			_aidaHistoMapProf1D[bname + "residualdZvsY"]->fill(estim.getY(), plane.getMeasZ() - meas.getZ()  );
			_aidaHistoMap2D[bname + "residualmeasZvsmeasX"]->fill(  meas.getZ()/1000., meas.getX()  );
			_aidaHistoMap2D[bname + "residualmeasZvsmeasY"]->fill(  meas.getZ()/1000., meas.getY()  );
			_aidaHistoMap2D[bname + "residualfitZvsmeasX"]->fill( plane.getMeasZ()/1000., meas.getX() );
			_aidaHistoMap2D[bname + "residualfitZvsmeasY"]->fill( plane.getMeasZ()/1000., meas.getY() );
 
			_aidaHistoMap2D[ "AllResidmeasZvsmeasX"]->fill(  meas.getZ()/1000., meas.getX()  );
			_aidaHistoMap2D[ "AllResidmeasZvsmeasY"]->fill(  meas.getZ()/1000., meas.getY()  );
			_aidaHistoMap2D[ "AllResidfitZvsmeasX"]->fill( plane.getMeasZ()/1000., meas.getX() );
			_aidaHistoMap2D[ "AllResidfitZvsmeasY"]->fill( plane.getMeasZ()/1000., meas.getY() );
			//Angles
			_aidaHistoMap[bname + "dxdz"]->fill( estim.getXdz() );
			_aidaHistoMap[bname + "dydz"]->fill( estim.getYdz() );

			if( ii != 4) { continue; }

			_aidaZvHitX->fill(estim.getX(), meas.getZ() - plane.getZpos());
			_aidaZvFitX->fill(estim.getX(), (plane.getMeasZ() - plane.getZpos()) - (meas.getZ() - plane.getZpos()));
			_aidaZvHitY->fill(estim.getY(), meas.getZ() - plane.getZpos());
			_aidaZvFitY->fill(estim.getY(), (plane.getMeasZ() - plane.getZpos()) - (meas.getZ() - plane.getZpos()));*/
		}
	}
}

void EUTelDafFitterAnalysis::bookDiffHistos() {

	streamlog_out(DEBUG5) << "Booking histograms\n";

	string basePath = "TrueHitResiduals";
	AIDAProcessor::tree(this)->mkdir(basePath.c_str());
	basePath.append("/");

	double residminX = -0.3;
	double residmaxX =  0.3;

	for (size_t i = 0; i < _system.planes.size(); i++) {

		string bname = static_cast< string >("pl") + to_string(plane.getSensorID()) + "_trueHit_";

		_aidaHistoMap[bname+"residualX"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D(basePath+bname+"residualX",600, residminX, residmaxX);
		_aidaHistoMap[bname+"residualY"] =  AIDAProcessor::histogramFactory(this)->createHistogram1D(basePath+bname+"residualY",600, residminX, residmaxX);
		//Resids 2D // profiles
		_aidaHistoMapProf1D[bname+"residualdXvsX"] = AIDAProcessor::histogramFactory(this)->createProfile1D(basePath+bname+"dXvsX", 200, -10000., 10000.,   residminX, residmaxX );
		_aidaHistoMapProf1D[bname+"residualdYvsX"] = AIDAProcessor::histogramFactory(this)->createProfile1D(basePath+bname+"dXvsY", 200, -10000., 10000.,   residminX, residmaxX );
		_aidaHistoMapProf1D[bname+"residualdXvsY"] = AIDAProcessor::histogramFactory(this)->createProfile1D(basePath+bname+"dYvsX", 200, -10000., 10000.,   residminX, residmaxX );
		_aidaHistoMapProf1D[bname+"residualdYvsY"]= AIDAProcessor::histogramFactory(this)->createProfile1D(basePath+bname+"dYvsY", 200, -10000., 10000.,   residminX, residmaxX );
	}
}

void EUTelDafFitter::dafEnd() {
}
