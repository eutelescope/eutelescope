/* 
 * File:   EUTelGBLFitter.cc
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Created on January 25, 2013, 2:53 PM
 */

#ifdef USE_GBL

// its own header 
#include "EUTelGBLFitter.h"

// eutelescope includes ".h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelUtilityRungeKutta.h"
#include "EUTELESCOPE.h"

// marlin util includes
#include "mille/Mille.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// GBL
#include "include/GblTrajectory.h"
#include "include/GblPoint.h"
#include "include/GblData.h"
#include "include/BorderedBandMatrix.h"
#include "include/MilleBinary.h"
#include "include/VMatrix.h"


// system includes <>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <iterator>
#include <algorithm>

namespace eutelescope {

    EUTelGBLFitter::EUTelGBLFitter() : EUTelTrackFitter("GBLTrackFitter"),
    _trackCandidates(),
    _gblTrackCandidates(),
    _parPropJac(5, 5),
    _beamQ(-1),
    _eBeam(4.),
    _hitId2GblPointLabel(),
    _hitId2GblPointLabelMille(),
		_PointToState(),
		_binaryname(),
    _alignmentMode( Utility::noAlignment ), // default - no Alignment 
    _mEstimatorType(),
    _mille(0),
    _parameterIdXShiftsMap(),
    _parameterIdYShiftsMap(),
    _parameterIdZShiftsMap(),
    _parameterIdXRotationsMap(),
    _parameterIdYRotationsMap(),
    _parameterIdZRotationsMap(),
    _excludeFromFit(),
    _chi2cut(1000.),
    _eomIntegrator( new EUTelUtilityRungeKutta() ),
    _eomODE( 0 )
    {
       streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter default constructor FitTrackVec: " << getFitTrackVec() << std::endl;
                // Initialise ODE integrators for eom and jacobian
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
       // reset parent members: 
       streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter constructor over" << std::endl;
    }

    EUTelGBLFitter::EUTelGBLFitter(std::string name) : EUTelTrackFitter(name),
    _trackCandidates(),
    _gblTrackCandidates(),
    _parPropJac(5, 5),
    _beamQ(-1),
    _eBeam(-1.),
    _hitId2GblPointLabel(),
    _hitId2GblPointLabelMille(),
		_PointToState(),
		_binaryname(),
    _alignmentMode( Utility::noAlignment ),
    _mEstimatorType(),
    _mille(0),
    _parameterIdXShiftsMap(),
    _parameterIdYShiftsMap(),
    _parameterIdZShiftsMap(),
    _parameterIdXRotationsMap(),
    _parameterIdYRotationsMap(),
    _parameterIdZRotationsMap(),
    _excludeFromFit(),
    _chi2cut(1000.),
    _eomIntegrator( new EUTelUtilityRungeKutta() ),
    _eomODE(  0 )
    {
        streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter " << name << " constructor FitTrackVec: " << getFitTrackVec() << std::endl;
               // Initialise ODE integrators for eom and jacobian
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
        streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter constructor over" << std::endl;
   }

    EUTelGBLFitter::~EUTelGBLFitter() {
    }
    
    void EUTelGBLFitter::setParamterIdPlaneVec( const std::vector<int>& vector)
    {
      _parameterIdPlaneVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdXResolutionVec( const std::vector<float>& vector)
    {
      _parameterIdXResolutionVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdYResolutionVec( const std::vector<float>& vector)
    {
      _parameterIdYResolutionVec = vector;
    }
       

    void EUTelGBLFitter::setParamterIdXRotationsMap( const std::map<int, int>& map ) {
        _parameterIdXRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdYRotationsMap( const std::map<int, int>& map ) {
        _parameterIdYRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdZRotationsMap( const std::map<int, int>& map ) {
        _parameterIdZRotationsMap = map;
    }

    void EUTelGBLFitter::setParamterIdZShiftsMap( const std::map<int, int>& map ) {
        _parameterIdZShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdYShiftsMap( const std::map<int, int>& map ) {
        _parameterIdYShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdXShiftsMap( const std::map<int, int>& map ) {
        _parameterIdXShiftsMap = map;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdXRotationsMap() const {
        return _parameterIdXRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdYRotationsMap() const {
        return _parameterIdYRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdZRotationsMap() const {
        return _parameterIdZRotationsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdZShiftsMap() const {
        return _parameterIdZShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdYShiftsMap() const {
        return _parameterIdYShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdXShiftsMap() const {
        return _parameterIdXShiftsMap;
    }

    void EUTelGBLFitter::setMEstimatorType( const std::string& mEstimatorType ) {
        std::string mEstimatorTypeLowerCase = mEstimatorType;
        std::transform( mEstimatorType.begin(), mEstimatorType.end(), mEstimatorTypeLowerCase.begin(), ::tolower);
        
        if ( mEstimatorType.size() != 1 ) {
            streamlog_out( WARNING1 ) << "More than one character supplied as M-estimator option" << std::endl;
            streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
            return;
        }
        
        if ( mEstimatorType.compare("t") == 0 ||
             mEstimatorType.compare("h") == 0 ||
             mEstimatorType.compare("c") == 0   ) this->_mEstimatorType = _mEstimatorType;
        else {
            streamlog_out( WARNING1 ) << "M-estimator option " << mEstimatorType << " was not recognized" << std::endl;
            streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
        }
    }

    std::string EUTelGBLFitter::getMEstimatorType( ) const {
        return _mEstimatorType;
    }

    const std::map<long, int>& EUTelGBLFitter::getHitId2GblPointLabel( ) const {
        return _hitId2GblPointLabel;
    }

    void EUTelGBLFitter::setExcludeFromFitPlanes( const std::vector<int>& excludedPlanes ) {
        this->_excludeFromFit = excludedPlanes;
    }

    std::vector<int> EUTelGBLFitter::getExcludeFromFitPlanes( ) const {
        return _excludeFromFit;
    }
     
    void EUTelGBLFitter::SetTrackCandidates( vector<const IMPL::TrackImpl*>& trackCandidatesVec) {

        this->_trackCandidatesVec = trackCandidatesVec	;
        return;
    }
       
    void EUTelGBLFitter::SetTrackCandidates( const EVENT::TrackVec& trackCandidates) {

        this->_trackCandidates = trackCandidates;
        return;
    }


    void EUTelGBLFitter::Clear() {

				_points.clear();
				_states.clear();
				_counter_num_pointer = 1;
    }

    void EUTelGBLFitter::pushBackPointandState( std::vector< gbl::GblPoint >* pointListTrack, gbl::GblPoint pointTrack, EUTelTrackStateImpl *state) {
				pointTrack.setLabel(_counter_num_pointer);
        pointListTrack->push_back(pointTrack);
       
        streamlog_out(DEBUG0) << endl << "pushBackPoint: " << pointListTrack->size() <<  std::endl;
        // store point's GBL label for future reference
			 	streamlog_out(DEBUG0) << endl << "This is the state and point " << state <<"," <<&(pointListTrack->back())<<"State hit: "<<state->getHit()<<std::endl;
				//This does not work. Maps I hate!!!!!!!!!!!!
        //_PointToState[*state] = &(pointListTrack->back()); 
				//OutputMap(_PointToState);
				_states.push_back(*state);
				_points.push_back(pointTrack);
				_counter_num_pointer++;
				
}

void EUTelGBLFitter::OutputMap(std::map< EUTelTrackStateImpl, gbl::GblPoint*, compare_points > mapex){
typedef std::map<EUTelTrackStateImpl,gbl::GblPoint*, compare_points >::const_iterator MapIterator;
for (MapIterator iter = mapex.begin(); iter != mapex.end(); iter++)
{
    cout << "Key: " << &(iter->first) << "Values:  "<< iter->second << endl;
   
}

}

void EUTelGBLFitter::UpdateTrackFromGBLTrajectory (gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >* pointList){
int i=0;
typedef std::vector<EUTelTrackStateImpl>::iterator MapIterator ;
for (MapIterator iter = _states.begin(); iter != _states.end(); iter++)
{
			EUTelTrackStateImpl & state = *iter;

			TVectorD corrections(5);
			TMatrixDSym correctionsCov(5,5);
			unsigned int pointNum = _points[i].getLabel();
      traj->getResults(pointNum, corrections, correctionsCov );

     	streamlog_out(DEBUG3) << endl << "State before we have added corrections: " << std::endl;
			state.Print();
			state.setX( state.getX() + corrections[3]);
			state.setY( state.getY() + corrections[4]);
			state.setTx( state.getTx() + corrections[1]);
			state.setTy( state.getTy() + corrections[2]);
			state.setInvP( state.getInvP() + corrections[0]);

			//This will wor for now but the when we tilt sensor no longer. Need to start troring z parameter. Since we need this to transfrom from global to local coordinates.
			float ref[3];
			ref[0] = state.getReferencePoint()[0]; ref[1] = state.getReferencePoint()[1];			ref[2] = 0;
			state.setReferencePoint(ref);
			
     	streamlog_out(DEBUG3) << endl << "State after we have added corrections: " << std::endl;
			state.Print();

		
i++;
	}//END OF LOOP OVER POINTS

}

//This used after trackfit will fill a map between (sensor ID and residual), (sensor ID and residual error).
void EUTelGBLFitter::getResidualOfTrackandHits(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >* pointList, map< int, map< float, float > > &  SensorResidualError){

int i=0;
typedef std::vector<EUTelTrackStateImpl>::iterator MapIterator ;
for (MapIterator iter = _states.begin(); iter != _states.end(); iter++)
{
			EUTelTrackStateImpl & state = *iter;
		if(state.getHit() != NULL){
			streamlog_out(DEBUG0) << endl << "There is a hit on the state. Hit pointer: "<< state.getHit()<<" Find update Residuals!" << std::endl;
  	unsigned int numData; //Not sure what this is used for??????
		TVectorD aResiduals(2);
		TVectorD aMeasErrors(2);
		TVectorD aResErrors(2,2);
		TVectorD aDownWeights(2); 
		unsigned int pointNum = _points[i].getLabel();
		traj->getMeasResults(pointNum, numData, aResiduals, aMeasErrors, aResErrors, aDownWeights);
		//Store the x and y component of residual. Need to change the naem of the container
		map<float, float> res_err; //This is create on the stack but will pass thisa by value to the new map so it 
		res_err.insert(make_pair(aResiduals[0],aResiduals[1]));
		SensorResidualError.insert(make_pair(state.getLocation(), res_err));		
		
		}
		else{
			streamlog_out(DEBUG0) << "The hit is NULL. State pointer: "<<&state<< " Hit pointer " << state.getHit() <<" So will not get residual" << std::endl;
		}
i++;
	}

} 

void EUTelGBLFitter::CreateTrajectoryandFit(std::vector< gbl::GblPoint >* pointList,  gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr){
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::CreateTrajectoryandFit -- BEGIN " << endl;

	double loss = 0.;

		streamlog_out ( DEBUG0 ) << "This is the trajectory we are just about to fit: " << endl;
	  streamlog_message( DEBUG0, traj->printTrajectory(10);, std::endl; );
		streamlog_out ( DEBUG0 ) << "This is the points in that trajectory " << endl;
	  streamlog_message( DEBUG0, traj->printPoints(10);, std::endl; );

	if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( *chi2, *ndf, loss, _mEstimatorType );
  else ierr = traj->fit( *chi2, *ndf, loss );

	if( ierr != 0 ){
		streamlog_out(MESSAGE0) << "Fit failed! Track error: "<< ierr << " and chi2: " << *chi2 << std::endl;
	}
	else{
  streamlog_out(MESSAGE0) << "Fit Successful! Track error; "<< ierr << " and chi2: " << *chi2 << std::endl;
	}
	

	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::CreateTrajectoryandFit -- END " << endl;
}







// convert input TrackCandidates and TrackStates into a GBL Trajectory
void EUTelGBLFitter::FillInformationToGBLPointObject(EUTelTrackImpl* EUtrack, std::vector< gbl::GblPoint >* pointList){
	// sanity check. Mustn't happen in principle. That the number of hits is greater than the number of hits
  if ( EUtrack->getTrackerHits().size() > geo::gGeometry().nPlanes() ){
  	streamlog_out(ERROR) << "Sanity check. This should not happen in principle. Number of hits is greater then number of planes" << std::endl;
   	return;
  }
  //Create the jacobian
  TMatrixD jacPointToPoint(5, 5);
  jacPointToPoint.UnitMatrix();
 	////////////////////////////////////////////////////////////////////////////////////////////////// loop through all states.
  for(int i=0;i < EUtrack->getTrackStates().size(); i++){		
		/////////////////////////////////////////////////////////////////////////////////////////////BEGIN to create GBL point 
  	streamlog_out(DEBUG3) << "The first GBL point is made from this jacobian:" << std::endl;
		TMatrixD output(5,5);
	  streamlog_message( DEBUG0, jacPointToPoint.Print();, std::endl; );
		changejacobainGBL(jacPointToPoint, output);
		gbl::GblPoint point(output);
  		EUTelTrackStateImpl* state = EUtrack->getTrackStates().at(i); //get the state for this track. Static cast from EVENT::TrackState to derived class IMPL::TrackStateImpl.
  		streamlog_out(DEBUG3) << "This is the track state being used in creation of GBL points" << std::endl;
			state->Print();

		//Need to find hit that this state may be associated with. Note this is a problem for two reasons. Not all states have a hit. Furthermore we can not associate a hit with a state with the current LCIO format. This must be fixed
		EVENT::TrackerHit* hit = state->getHit();
		//FindHitIfThereIsOne(EUtrack, hit, state); //This will point the hit to the correct hit object associated with this state. If non exists then point it will remain pointed to NULL. Not needed state holds hit
		double fitPointLocal[] = {0.,0.,0.};
  	fitPointLocal [0] = state->getReferencePoint()[0] ;
  	fitPointLocal [1] = state->getReferencePoint()[1] ;
  	fitPointLocal [2] = state->getReferencePoint()[2] ;
		//addSiPlaneScattererGBL(point, state->getLocation()); //This we still functions still assumes silicon is the thin scatterer. This can be easily changed when we have the correct gear file. However we will always assume that states will come with scattering information. To take into account material between states this will be dealt with latter. 
		double fitPointGlobal[3];
		geo::gGeometry().local2Master( state->getLocation(), fitPointLocal, fitPointGlobal);
		state->setZParameter(fitPointGlobal[2]); //This is needed to calculate jacobian for some reason. Need to check this.
		streamlog_out(DEBUG3) << "This is the global position of the track state. Should be the same x,y as above: " <<fitPointGlobal[0]<<","<<fitPointGlobal[1]<<","<<fitPointGlobal[2]<< std::endl;	
		if(hit != NULL){
			SetHitCovMatrixFromFitterGBL(state);
			double cov[4] ;
			state->getTrackStateHitCov(cov); //This part should not be done in the way it has. MUST FIX! Hit cov should be part of hits own class. Must learn more about LCIO data format
			
			addMeasurementGBL(point, hit->getPosition(),  fitPointLocal, cov, state->getH());
		streamlog_out(DEBUG3) << "Just before adding state hit pointer is " << state->getHit() <<std::endl; 		
			pushBackPointandState(pointList, point, state);

		}else{
			streamlog_out(DEBUG3) << "The state has no hit so just just at point with no measurement" <<std::endl;
			pushBackPointandState(pointList, point, state);
		}

		////////////////////////////////////////////////////////////////////////////////START TO CREATE SCATTERS BETWEEN PLANES
		if(i != (EUtrack->getTrackStates().size()-1)){
  			EUTelTrackStateImpl* state_next = EUtrack->getTrackStates().at(i+1); //get the next tracks state to determine dz between the two states
  		streamlog_out(DEBUG3) << "This is the track state that is one state ahead" << std::endl;
			state->Print();
			double fitPointLocal_next[] = {0.,0.,0.};  //Need this since we dont save the z parameter as state variable
			fitPointLocal_next [0] = state_next->getReferencePoint()[0] ;
  			fitPointLocal_next [1] = state_next->getReferencePoint()[1] ;
	  		fitPointLocal_next [2] = state_next->getReferencePoint()[2] ;

			double fitPointGlobal_next[3];
			geo::gGeometry().local2Master( state_next->getLocation(), fitPointLocal_next, fitPointGlobal_next );
			state_next->setZParameter(fitPointGlobal_next[2]); //This is needed to calculate jacobian for some reason. Need to check this. This is not really needed for state_next
		streamlog_out(DEBUG3) << "This is the global position of the track state. Should be the same x,y as above: " <<fitPointGlobal_next[0]<<","<<fitPointGlobal_next[1]<<","<<fitPointGlobal_next[2]<< std::endl;	
			float rad = geo::gGeometry().findRadLengthIntegral( fitPointGlobal, fitPointGlobal_next, true ); //We need to skip the volumes that contain the hits since this has already been counted. Must check this functions as expect????
			streamlog_out(DEBUG3) << "This is the radiation length between the two points  " << fitPointGlobal[0]<<","<<fitPointGlobal[1]<<","<<fitPointGlobal[2]<<" and  " <<fitPointGlobal_next[0]<<","<<fitPointGlobal_next[1]<<","<<fitPointGlobal_next[2] <<"Radition length:  "<< rad <<std::endl;

			///////////////////////////////////////////////////////////////////////////////////////////////////////BEGIN THE FIRST SCATTERING PLANE
			//These distances are from the last state plane. There are where the next scatterer should be
			//float distance1 = (fitPointGlobal_next[2] + fitPointGlobal[2])/2 - (fitPointGlobal_next[2] - fitPointGlobal[2])/sqrt(12); 
			//The original plane should always be 0 so the above expression does not work
			float distance1 = (fitPointGlobal_next[2] - fitPointGlobal[2])/2 - (fitPointGlobal_next[2] - fitPointGlobal[2])/sqrt(12); 
			streamlog_out(DEBUG3) << "This is the distance to the first scatterer: " << distance1 <<std::endl;	
			//Note the distance is used along the track and not from the scattering plane. How should this be dealt with?
			TMatrixD jacobianScat1(5,5); jacobianScat1 = state->getPropagationJacobianF(distance1);
			TVectorD stateVec = state->getTrackStateVec();
			TVectorD stateVecOnScat1 = jacobianScat1 * stateVec; 
  		streamlog_out(DEBUG3) << "The first scattering point point is made from this jacobian:" << std::endl;
	  	streamlog_message( DEBUG0, jacobianScat1.Print();, std::endl; );
			changejacobainGBL(jacobianScat1, output);
			gbl::GblPoint pointScat1(output);
			TVectorD scat(2);
			scat[0] = 0.0; //This should always be 0 right? If not then it should be given as a parameter
			scat[1] = 0.0; 

 			const double scatvariance  = Utility::getThetaRMSHighland(GetBeamEnergy(), rad/2);
			streamlog_out(DEBUG3) << "Scattering mean angle and variance. Average angle: " << scat[0] <<"," <<scat[1] << ".Variance " <<scatvariance<<std::endl;	
			TVectorD scatPrecSensor(2);
 			scatPrecSensor[0] = 1.0 / (scatvariance * scatvariance );
 			scatPrecSensor[1] = 1.0 / (scatvariance * scatvariance );

  		pointScat1.addScatterer(scat, scatPrecSensor);
		//	pushBackPointandState(pointList, pointScat1, NULL);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////END THE FIRST SCATTERING PLANE
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////BEGIN THE SECOND SCATTERING PLANE
			//float distance2 = (fitPointGlobal_next[2] + fitPointGlobal[2])/2 + (fitPointGlobal_next[2] - fitPointGlobal[2])/sqrt(12);
			//The original plane should always be 0 so the above expression does not work
			float distance2 = (fitPointGlobal_next[2] - fitPointGlobal[2])/2 + (fitPointGlobal_next[2] - fitPointGlobal[2])/sqrt(12) - distance1; //Since we want the distance between the scatterers
			double p  = 1. / (stateVecOnScat1[4] * -1);
  		double px = p*stateVecOnScat1[2] / sqrt( 1. + stateVecOnScat1[2]*stateVecOnScat1[2] + stateVecOnScat1[3]*stateVecOnScat1[3] );
  		double py = p*stateVecOnScat1[3] / sqrt( 1. + stateVecOnScat1[2]*stateVecOnScat1[2] + stateVecOnScat1[3]*stateVecOnScat1[3] );
  		double pz = p    / sqrt( 1. + stateVecOnScat1[2]*stateVecOnScat1[2] + stateVecOnScat1[3]*stateVecOnScat1[3] );
			TMatrixD jacobianScat2(5,5); jacobianScat2 = geo::gGeometry().getPropagationJacobianF( stateVecOnScat1[0], stateVecOnScat1[1], 0, px, py, pz, -1, distance2 ); //Set z to zero since the magnetic field is assumed homogenious and this is all that z is used for
			TVectorD stateVecOnScat2 = jacobianScat2 * stateVecOnScat1; 
  		streamlog_out(DEBUG3) << "The second scattering point point is made from this jacobian:" << std::endl;
	  	streamlog_message( DEBUG0, jacobianScat2.Print();, std::endl; );
				changejacobainGBL(jacobianScat2, output);
			gbl::GblPoint pointScat2(output);


  		pointScat2.addScatterer(scat, scatPrecSensor);
		//	pushBackPointandState(pointList, pointScat2, NULL);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////END OF SECOND SCATTERING PLANE

			p  = 1. / (stateVecOnScat2[4] * -1);
  		px = p*stateVecOnScat2[2] / sqrt( 1. + stateVecOnScat2[2]*stateVecOnScat2[2] + stateVecOnScat2[3]*stateVecOnScat2[3] );
  		py = p*stateVecOnScat2[3] / sqrt( 1. + stateVecOnScat2[2]*stateVecOnScat2[2] + stateVecOnScat2[3]*stateVecOnScat2[3] );
  		pz = p    / sqrt( 1. + stateVecOnScat2[2]*stateVecOnScat2[2] + stateVecOnScat2[3]*stateVecOnScat2[3] );
			//jacPointToPoint =  geo::gGeometry().getPropagationJacobianF( stateVecOnScat2[0], stateVecOnScat2[1], 0, px, py, pz, -1, fitPointGlobal_next[2] - distance2 );

//This part is just to check without scatterers
float distance = (fitPointGlobal_next[2] -  fitPointGlobal[2]);
jacPointToPoint = state->getPropagationJacobianF(distance);



		}  
		/////////////////////////////////////////////////////////////////////////////////////////END OF CREATE SCATTERERS BETWEEN PLANES
	
	}//END OF LOOP THROUGH ALL PLANES
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
}

//This matrix is a hack. It changes the Jacobian in geometry to what is needed as input to GBL. 
void EUTelGBLFitter::changejacobainGBL(TMatrixD & jacobianF, TMatrixD & _parPropJac){
	_parPropJac.UnitMatrix();

_parPropJac[1][0] =	jacobianF[2][4]; _parPropJac[1][2] = 	jacobianF[2][3];	
_parPropJac[2][0] = jacobianF[3][4]; _parPropJac[2][1] = jacobianF[3][2];	
_parPropJac[3][0] = jacobianF[0][4]; _parPropJac[3][1] = jacobianF[0][2];	_parPropJac[3][2] = jacobianF[0][3];	
_parPropJac[4][0] = jacobianF[1][4]; _parPropJac[4][1] = jacobianF[1][2];	_parPropJac[4][2] = jacobianF[1][3];	


} 

void EUTelGBLFitter::addSiPlaneScattererGBL(gbl::GblPoint& point, int iPlane) {

	streamlog_out(MESSAGE1) << " addSiPlaneScattererGBL ------------- BEGIN --------------  " << std::endl;


	TVectorD scatPrecSensor(2);
	TVectorD scat(2); 
	scat[0] = 0.0; scat[1]=0.0; //This should always be 0 right? If not then it should be given as a parameter
	const double radlenSi           = geo::gGeometry().siPlaneRadLength(iPlane);
	const double thicknessSi        = geo::gGeometry().siPlaneZSize(iPlane);

	const double X0Si  = thicknessSi / radlenSi; // Si 
               
	
	const double tetSi  = Utility::getThetaRMSHighland(GetBeamEnergy(), X0Si);

        scatPrecSensor[0] = 1.0 / (tetSi * tetSi );
        scatPrecSensor[1] = 1.0 / (tetSi * tetSi );

        point.addScatterer(scat, scatPrecSensor);

	streamlog_out(MESSAGE1) << " addSiPlaneScattererGBL  ------------- END ----------------- " << std::endl;
}



//This will add the measurement of the hit and predicted position. Using the covariant matrix of the hit. NOT! the residual.
void EUTelGBLFitter::addMeasurementGBL(gbl::GblPoint& point, const double *hitPos, const double *statePos, double hitCov[4], TMatrixD HMatrix){
     
	streamlog_out(MESSAGE1) << " addMeasurementsGBL ------------- BEGIN --------------- " << std::endl;

 	TVectorD meas(2);
	meas[0] = hitPos[0] - statePos[0];
        meas[1] = hitPos[1] - statePos[1];

	TVectorD measPrec(5); //Precision matrix is symmetric. The vector option that was here was silly since there could be correlation between variance and x/y. However for now leave it.
        measPrec[0] = 1. / hitCov[0];	// cov(x,x)
       // measPrec[0][1] = 1. / hitCov[1];	// cov(x,x)
       // measPrec[1][0] = 1. / hitCov[2];	// cov(x,x)
        measPrec[1] = 1. / hitCov[3];	// cov(y,y)
	//measPrec[0][1] = 1. / hitCov[1];  //cov(x,y)
	streamlog_out(DEBUG4) << "This is what we add to the measured point:" << std::endl;
	streamlog_out(DEBUG4) << "Residuals and covariant matrix for the hit:" << std::endl;
        streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] <<"," << std::endl;
        streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20)  <<"," << measPrec[1] << std::endl;
	streamlog_out(DEBUG4) << "This H matrix:" << std::endl;
	streamlog_message( DEBUG0, HMatrix.Print();, std::endl; );

//The gbl library creates 5 measurement vector and 5x5 propagation matrix automatically. 
//So need to alter how HMatrix is added to this. 
TMatrixD proM2l(2, 2);
proM2l.UnitMatrix();
proM2l[0][0]=HMatrix[0][0]; proM2l[0][1]=HMatrix[0][1];
proM2l[1][0]=HMatrix[1][0]; proM2l[1][1]=HMatrix[1][1];

        point.addMeasurement(proM2l, meas, measPrec);

	streamlog_out(MESSAGE1) << " addMeasurementsGBL ------------- END ----------------- " << std::endl;
}

//This is use in alignment.
void EUTelGBLFitter::CreateAlignmentToMeasurementJacobian(std::vector< gbl::GblPoint >& pointList ){

int i=0;
typedef std::vector<EUTelTrackStateImpl>::iterator MapIterator ;
for (MapIterator iter = _states.begin(); iter != _states.end(); iter++)
{
			EUTelTrackStateImpl & state = *iter;
		
		if(state.getHit() != NULL){
		_MilleInterface->CreateAlignmentToMeasurementJacobian(&state); //Get the jacobain for that state and sensor
		_MilleInterface->CreateGlobalLabels(&state);  //Gets the coorect labels for that sensor
		TMatrixD &  Jac = _MilleInterface->getAlignmentJacobian();
		std::vector<int> labels =  _MilleInterface->getGlobalParameters();
		
			
		pointList[i].addGlobals(labels, Jac);

		

		}//END OF IF STATEMENT IF THE THERE IS A STATE WITH HIT.  REMOVE THIS FOR NOW
		else{		
			streamlog_out(DEBUG0)<<"There was no hit so will not create jacobian and global variables for state"<<std::endl;
		}
i++;
	}


}


//This should be in hit object
void EUTelGBLFitter::SetHitCovMatrixFromFitterGBL(EUTelTrackStateImpl *state){

	double hitcov[4];

	int izPlane = state->getLocation();
	if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 ){
		if( state->getHit() != NULL ){
  		hitcov[0] = _parameterIdXResolutionVec[izPlane];
    	hitcov[3] = _parameterIdYResolutionVec[izPlane];

   		hitcov[0] *= hitcov[0]; // squared !
   		hitcov[2] *= hitcov[2]; // squared !
  	}
	}
	else{
		hitcov[0]=1;
		hitcov[1]=0;
		hitcov[2]=0;
		hitcov[3]=1;
	}

state->setTrackStateHitCov(hitcov);

}

} // namespace eutelescope


#endif
