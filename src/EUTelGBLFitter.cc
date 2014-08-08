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


    void EUTelGBLFitter::resetPerTrack() {
			_counter_num_pointer=1;
			_measurementStatesInOrder.clear();
    }
//Not all points are states so we need a way to do:
//state->label->point. This is what this function does.  
void EUTelGBLFitter::setPointVec( std::vector< gbl::GblPoint >& pointList, gbl::GblPoint& point) {
	point.setLabel(_counter_num_pointer);
	_counter_num_pointer++;
	pointList.push_back(point);
	streamlog_out(DEBUG0) << endl << "pushBackPoint size: " << pointList.size() <<  std::endl;
}

void EUTelGBLFitter::OutputMap(std::map< EUTelTrackStateImpl, gbl::GblPoint*, compare_points > mapex){
typedef std::map<EUTelTrackStateImpl,gbl::GblPoint*, compare_points >::const_iterator MapIterator;
	for (MapIterator iter = mapex.begin(); iter != mapex.end(); iter++)
	{
    cout << "Key: " << &(iter->first) << "Values:  "<< iter->second << endl;
   
	}
}
//We use the labels for trajectory and equate them to the pointList index.
//We create them when creating the points. However trajectory will overwrite these.
void EUTelGBLFitter::setListStateAndLabelAfterTrajectory(std::vector< gbl::GblPoint >& pointList, gbl::GblTrajectory* traj){
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setListStateAndLabelAfterTrajectory-- BEGIN " << endl;
	_vectorOfPairsMeasurementStatesAndLabels.clear();
	streamlog_out(DEBUG5)<<"The number of measurement states is: "<< _measurementStatesInOrder.size()<<std::endl;
	int	counter=0;
	std::vector<unsigned int> labels;
	traj->getLabels(labels);
	if((_counter_num_pointer-1) != labels.size() and (_counter_num_pointer-1) != labels.at((labels.size()-1))){
		throw(lcio::Exception(Utility::outputColourString("The labels given by trajectory and what we labeled does not seem to watch .", "RED")));
	}
	for(int i=0; i<pointList.size();++i){
		if(pointList.at(i).hasMeasurement()>0){//Here we assume that traj has ordered the pointList the same.
			streamlog_out(DEBUG0)<<"Measurement found! Pair (state,label)  ("<<  &(_measurementStatesInOrder.at(counter))<<","<<labels.at(i)<<")"<<std::endl; 
			_vectorOfPairsMeasurementStatesAndLabels.push_back(make_pair(_measurementStatesInOrder.at(counter), labels.at(i)));	
			counter++;
		}
		if(counter == (_measurementStatesInOrder.size()+1)){//Add one since you will make an extra loop
			throw(lcio::Exception(Utility::outputColourString("The counter is larger than the number of states saved as measurements.", "RED")));
		}
	}
	if(counter !=  (_measurementStatesInOrder.size())){//Since we make an extra loop and counter++ again
		throw(lcio::Exception(Utility::outputColourString("We did not add all the states .", "RED")));
	}
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setListStateAndLabelAfterTrajectory-- END " << endl;
}
//This code must be repeated since we need to create the ling between states and labels before trajectory in alignment
void EUTelGBLFitter::setListStateAndLabelBeforeTrajectory(std::vector< gbl::GblPoint >& pointList){
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setListStateAndLabelBeforeTrajectory-- BEGIN " << endl;
	_vectorOfPairsMeasurementStatesAndLabels.clear();
	streamlog_out(DEBUG5)<<"The number of measurement states is: "<< _measurementStatesInOrder.size()<<std::endl;
	int	counter=0;
	for(int i=0; i<pointList.size();++i){
		if(pointList.at(i).hasMeasurement()>0){//Here we assume that traj has ordered the pointList the same.
			streamlog_out(DEBUG0)<<"Measurement found! Pair (state,label)  ("<<  &(_measurementStatesInOrder.at(counter))<<","<<pointList.at(i).getLabel()<<")"<<std::endl; 
			_vectorOfPairsMeasurementStatesAndLabels.push_back(make_pair(_measurementStatesInOrder.at(counter), pointList.at(i).getLabel()));	
			counter++;
		}
		if(counter == (_measurementStatesInOrder.size()+1)){//Add one since you will make an extra loop
			throw(lcio::Exception(Utility::outputColourString("The counter is larger than the number of states saved as measurements.", "RED")));
		}
	}
	if(counter !=  (_measurementStatesInOrder.size())){//Since we make an extra loop and counter++ again
		throw(lcio::Exception(Utility::outputColourString("We did not add all the states .", "RED")));
	}
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setListStateAndLabelBeforeTrajectory-- END " << endl;
}

void EUTelGBLFitter::UpdateTrackFromGBLTrajectory (gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >& pointList,EUTelTrack &track){
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::UpdateTrackFromGBLTrajectory-- BEGIN " << endl;
	for(int i=0;i < track.getStatesPointers().size(); i++){//We get the pointers no since we want to change the track state contents		
		EUTelState* state = track.getStatesPointers().at(i);
		TVectorD corrections(5);
		TMatrixDSym correctionsCov(5,5);
		for(int j=0 ; j< _vectorOfPairsMeasurementStatesAndLabels.size();++j){
			if(_vectorOfPairsMeasurementStatesAndLabels.at(j).first == *state){
				streamlog_out(DEBUG0)<<"The loop number for states with measurements is: " << j << ". The label is: " << _vectorOfPairsMeasurementStatesAndLabels.at(j).second <<endl; 
				if(getLabelToPoint(pointList,_vectorOfPairsMeasurementStatesAndLabels.at(j).second).hasMeasurement() == 0){//TO DO: Some states will not have hits in the future so should remove. Leave for now to test
					throw(lcio::Exception(Utility::outputColourString("This point does not contain a measurements. Labeling of the state must be wrong ", "RED")));
				} 
				streamlog_out(DEBUG0)<<"To update track we use label: "<<_vectorOfPairsMeasurementStatesAndLabels.at(j).second<<std::endl; 
				traj->getResults(_vectorOfPairsMeasurementStatesAndLabels.at(j).second, corrections, correctionsCov );
				streamlog_out(DEBUG3) << endl << "State before we have added corrections: " << std::endl;
				state->print();
				TVectorD newStateVec(5);
				newStateVec[0] = state->getOmega() + corrections[0];
				newStateVec[1] = state->getDirectionXY()+corrections[1];
				newStateVec[2] = state->getDirectionYZ()+corrections[2];
				newStateVec[3] = state->getPosition()[0]+corrections[3];
				newStateVec[4] = state->getPosition()[1]+corrections[4]; 
				state->setTrackStateVecPlusZParameter(newStateVec,state->getPosition()[2]);//TO DO:We set the z position to be the same. However if the sensor is tilted this should change. How should we deal with this
				streamlog_out(DEBUG3) << endl << "State after we have added corrections: " << std::endl;
				state->print();
				break;
			}
		}//END of loop of all states with hits	
	}//END of loop of all states
}

//This used after trackfit will fill a map between (sensor ID and residualx/y).
void EUTelGBLFitter::getResidualOfTrackandHits(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint > pointList,EUTelTrack& track, map< int, map< float, float > > &  SensorResidual){
	for(int j=0 ; j< _vectorOfPairsMeasurementStatesAndLabels.size();j++){
		EUTelState state = _vectorOfPairsMeasurementStatesAndLabels.at(j).first;
		if(getLabelToPoint(pointList,_vectorOfPairsMeasurementStatesAndLabels.at(j).second).hasMeasurement() == 0){
			throw(lcio::Exception(Utility::outputColourString("This point does not contain a measurements. Labeling of the state must be wrong ", "RED")));
		} 
		streamlog_out(DEBUG0) << endl << "There is a hit on the state. Hit pointer: "<< state.getTrackerHits()[0]<<" Find updated Residuals!" << std::endl;
		unsigned int numData; //Not sure what this is used for??????
		TVectorD aResiduals(2);
		TVectorD aMeasErrors(2);
		TVectorD aResErrors(2);
		TVectorD aDownWeights(2); 
		streamlog_out(DEBUG0)<<"To get residual of states we use label: "<<_vectorOfPairsMeasurementStatesAndLabels.at(j).second<<std::endl; 
		traj->getMeasResults(_vectorOfPairsMeasurementStatesAndLabels.at(j).second, numData, aResiduals, aMeasErrors, aResErrors, aDownWeights);
		streamlog_out(DEBUG0) <<"State location: "<<state.getLocation()<<" The residual x " <<aResiduals[0]<<" The residual y " <<aResiduals[1]<<endl;
		map<float, float> res; //This is create on the stack but will pass thisa by value to the new map so it 
		res.insert(make_pair(aResiduals[0],aResiduals[1]));
		SensorResidual.insert(make_pair(state.getLocation(), res));		
	}
}

void EUTelGBLFitter::computeTrajectoryAndFit(std::vector< gbl::GblPoint >& pointList,  gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr){
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::computeTrajectoryAndFit-- BEGIN " << endl;
	double loss = 0.;
	streamlog_out ( DEBUG0 ) << "This is the trajectory we are just about to fit: " << endl;
	streamlog_message( DEBUG0, traj->printTrajectory(10);, std::endl; );
	streamlog_out ( DEBUG0 ) << "This is the points in that trajectory " << endl;
	streamlog_message( DEBUG0, traj->printPoints(10);, std::endl; );


	if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( *chi2, *ndf, loss, _mEstimatorType );
  else ierr = traj->fit( *chi2, *ndf, loss );

	if( ierr != 0 ){
		streamlog_out(MESSAGE0) << Utility::outputColourString("Fit failed!","YELLOW")<<" Track error: "<< ierr << " and chi2: " << *chi2 << std::endl;
	}
	else{
  streamlog_out(MESSAGE0) << Utility::outputColourString("Fit Successful!","GREEN")<<" Track error; "<< ierr << " and chi2: " << *chi2 << std::endl;
	}
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::computeTrajectoryAndFit -- END " << endl;
}

void EUTelGBLFitter::testTrack(EUTelTrack& track){
	streamlog_out(DEBUG4)<<"EUTelGBLFitter::testTrack------------------------------------BEGIN"<<endl;
	if(track.getStates().size() == 0 ){
		throw(lcio::Exception(Utility::outputColourString("The number of states is zero.", "RED")));
	}
///Note we do not use excluded planes here. This should be dealt with in pattern recognition.
  if (track.getNumberOfHitsOnTrack() > geo::gGeometry().nPlanes() ){
		throw(lcio::Exception(Utility::outputColourString("The number of hits on the track is greater than the number of planes.", "RED"))); 	
  }
	streamlog_out(DEBUG5)<< Utility::outputColourString("Input track passed tests!", "GREEN")<<std::endl;
	streamlog_out(DEBUG4)<<"EUTelGBLFitter::testTrack------------------------------------END"<<endl;

} 
// convert input TrackCandidates and TrackStates into a GBL Trajectory
void EUTelGBLFitter::setInformationForGBLPointList(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList){
	streamlog_out(DEBUG4)<<"EUTelGBLFitter::setInformationForGBLPointList-------------------------------------BEGIN"<<endl;
  TMatrixD jacPointToPoint(5, 5);
  jacPointToPoint.UnitMatrix();
  for(int i=0;i < track.getStates().size(); i++){		
		streamlog_out(DEBUG3) << "The jacobian to get to this state jacobian on state number: " << i<<" Out of a total of states "<<track.getStates().size() << std::endl;
	  streamlog_message( DEBUG0, jacPointToPoint.Print();, std::endl; );
		gbl::GblPoint point(jacPointToPoint);
		EUTelState state = track.getStates().at(i);
		EUTelState nextState;
		if(i != (track.getStates().size()-1)){
			nextState = track.getStates().at(i+1);
		}
		setScattererGBL(point,state.getLocation());//Every sensor will have scattering due to itself. 
		if(state.getTrackerHits().size() == 0 ){
			streamlog_out(DEBUG3)  << Utility::outputColourString("This state does not have a hit. ", "YELLOW")<<std::endl;
			setPointVec(pointList, point);//This creates the vector of points and keeps a link between states and the points they created
		}else{
			double localPositionForState[3];
			//TO DO:Fix master2Localtwo so that is is the only one. This is currently a hack
			const double referencePoint[]	= {state.getReferencePoint()[0], state.getReferencePoint()[1],state.getReferencePoint()[2]};//Need this since geometry works with const doubles not floats 
			geo::gGeometry().master2Localtwo( state.getLocation(), referencePoint, localPositionForState );
			setMeasurementCov(state);
			double cov[4] ;
			state.getCombinedHitAndStateCovMatrixInLocalFrame(cov);
			setMeasurementGBL(point, state.getTrackerHits()[0]->getPosition(),  localPositionForState,  cov, state.getProjectionMatrix());
			_measurementStatesInOrder.push_back(state);//This is list of measurements states in the correct order. This is used later to associate states with point labels
			setPointVec(pointList, point);
		}//End of else statement if there is a hit.

		if(i != (track.getStates().size()-1)){//We do not produce scatterers after the last plane
			double stateReferencePoint[3];
			stateReferencePoint[0]=state.getPosition()[0];//Position is in mm
			stateReferencePoint[1]=state.getPosition()[1];
			stateReferencePoint[2]=state.getPosition()[2];
			double nextStateReferencePoint[3];
			nextStateReferencePoint[0]=nextState.getPosition()[0];
			nextStateReferencePoint[1]=nextState.getPosition()[1];
			nextStateReferencePoint[2]=nextState.getPosition()[2];
			testDistanceBetweenPoints(stateReferencePoint,nextStateReferencePoint);
			float rad = geo::gGeometry().findRadLengthIntegral(stateReferencePoint,nextStateReferencePoint, false );//TO DO: This adds the radiation length of the plane again. If you chose true the some times it returns 0 
			findScattersZPositionBetweenTwoStates(state, nextState); 
			jacPointToPoint=findScattersJacobians(state,nextState);
			setPointListWithNewScatterers(pointList, rad);
		}else{
			streamlog_out(DEBUG3)<<"We have reached the last plane"<<std::endl;
		}
	} 
	streamlog_out(DEBUG4)<<"EUTelGBLFitter::setInformationForGBLPointList-------------------------------------END"<<endl;

}
void EUTelGBLFitter::testDistanceBetweenPoints(double* position1,double* position2){
	TVectorD displacement(3);
	displacement(0) = position1[0] - position2[0];
	displacement(1) = position1[1] - position2[1];
	displacement(2) = position1[2] - position2[2];
	float distance= sqrt(displacement.Norm2Sqr());
	streamlog_out(DEBUG0)<<"The distance between the points to calculate radiation length is: "<<distance<<std::endl;
	if(distance<1){
		throw(lcio::Exception(Utility::outputColourString("The distance between the states is too small.", "RED"))); 	
	}
	if(distance>1000){//The planes should not be over a 1m in width
		throw(lcio::Exception(Utility::outputColourString("The distance between the planes is over 1m.", "RED"))); 	
	}
}

void EUTelGBLFitter::setPointListWithNewScatterers(std::vector< gbl::GblPoint >& pointList, float rad ){
	if(_scattererJacobians.size() != _scattererPositions.size()){
		throw(lcio::Exception(Utility::outputColourString("The size of the scattering positions and jacobians is different.", "RED"))); 	
	}
	if(rad == 0){
		throw(lcio::Exception(Utility::outputColourString("The provides radiation length to create scatterers it zero.", "RED"))); 	
	}
	for(int i = 0 ;i < _scattererJacobians.size()-1;++i){//The last jacobain is used to get to the plane! So only loop over to (_scatterJacobians-1)
		gbl::GblPoint point(_scattererJacobians[i]);
		point.setLabel(_counter_num_pointer);
		_counter_num_pointer++;
		setScattererGBL(point,rad/2);//TO DO:Simply dividing by 2 will work when there is only two planes
		pointList.push_back(point);
	}
}
//We want to create a jacobain from (Plane1 -> scatterer1) then (scatterer1->scatterer2) then (scatter2->plane2). We return the last jacobain
TMatrixD EUTelGBLFitter::findScattersJacobians(EUTelState state, EUTelState nextState){
	_scattererJacobians.clear();
	EUTelState loopState = state;
	for(int i=0;i<_scattererPositions.size();i++){
		_scattererJacobians.push_back(loopState.computePropagationJacobianFromStateToThisZLocation(_scattererPositions[i]));
		TVectorD nextStateVec = _scattererJacobians.at(i) * loopState.getTrackStateVec(); 
		streamlog_out(DEBUG5)<<"The new state vector on the scatter is: " <<std::endl;
		//nextStateVec.Print();
		streamlog_out(DEBUG5)<<"The z position of this scatter is: "<< _scattererPositions[i] <<std::endl;
		loopState.setTrackStateVecPlusZParameter(nextStateVec, _scattererPositions[i]);
		streamlog_out(DEBUG5)<<"The z position of loop state: "<< loopState.getPosition()[2]<<std::endl;

	}
	if(_scattererJacobians.size() != 3){
		throw(lcio::Exception(Utility::outputColourString("There is not 3 jacobians produced by scatterers!.", "RED"))); 	
	}
	return _scattererJacobians.back();//return the last jacobian so the next state can use this
}
//The distance from the first state to the next scatterer and then from that scatterer to the next all the way to the next state. 
//TO DO: This uses the optimum positions as described by Claus.  However for non homogeneous material distribution this might not be the case.
void EUTelGBLFitter::findScattersZPositionBetweenTwoStates(EUTelState& state, EUTelState& nextState){
	_scattererPositions.clear();	
	float distance1 = (nextState.getReferencePoint()[2] + state.getReferencePoint()[2])/2 - (nextState.getReferencePoint()[2] - state.getReferencePoint()[2])/sqrt(12); 
	_scattererPositions.push_back(distance1);//Z position of 1st scatter	
	float distance2 = (nextState.getReferencePoint()[2] + state.getReferencePoint()[2])/2 + (nextState.getReferencePoint()[2] - state.getReferencePoint()[2])/sqrt(12);//Z position of 2nd scatter 
	_scattererPositions.push_back(distance2);
	_scattererPositions.push_back(nextState.getReferencePoint()[2]); 
}
//Note that we take the planes themselfs at scatters and also add scatterers to simulate the medium inbetween. 
void EUTelGBLFitter::setScattererGBL(gbl::GblPoint& point, int iPlane) {
	streamlog_out(DEBUG1) << " setScattererGBL ------------- BEGIN --------------  " << std::endl;
	TVectorD scatPrecSensor(2);
	TVectorD scat(2); 
	scat[0] = 0.0; scat[1]=0.0;//TO DO: This will depend on the direction of the beam.  
	//TO DO:Should make the reading for any planes not just si. 
	const double radiationLength  = geo::gGeometry().siPlaneRadLength(iPlane);//This is in mm.
	const double thickness        = geo::gGeometry().siPlaneZSize(iPlane);//TO DO: Need to get the correct thickness
	streamlog_out(DEBUG5)<<"The radiation length for the sensor is: " << radiationLength <<std::endl; 
	if( radiationLength== 0){
		throw(lcio::Exception(Utility::outputColourString("Radiation length is zero.", "RED"))); 	
	}
	streamlog_out(DEBUG5)<<"The thickness for the sensor is: " << thickness <<std::endl; 
	if(thickness == 0 ){ 
		throw(lcio::Exception(Utility::outputColourString("thickness is zero", "RED"))); 	
	}
	const double x0  = thickness / radiationLength; 
	streamlog_out(DEBUG5)<<"The X0 for the sensor is: " << x0 <<std::endl; 
	if(x0 == 0){
		streamlog_out(MESSAGE0)<<"The thickness: " << thickness << "The radiation length is: "<< radiationLength <<std::endl; 
		throw(lcio::Exception(Utility::outputColourString("X0 is zero for the plane itself.", "RED"))); 	
	}
	if(getBeamEnergy() == 0 ){ 
		throw(lcio::Exception(Utility::outputColourString("Beam energy is zero", "RED"))); 	
	}

	const double scatteringVariance  = Utility::getThetaRMSHighland(getBeamEnergy(), x0);

	scatPrecSensor[0] = 1.0 / (scatteringVariance *scatteringVariance);
	scatPrecSensor[1] = 1.0 / (scatteringVariance * scatteringVariance);

	point.addScatterer(scat, scatPrecSensor);

	streamlog_out(DEBUG1) << "  setScattererGBL  ------------- END ----------------- " << std::endl;
}
//This is used when the we know the radiation length already
void EUTelGBLFitter::setScattererGBL(gbl::GblPoint& point, float x0 ) {
	streamlog_out(MESSAGE1) << " setScattererGBL ------------- BEGIN --------------  " << std::endl;
	TVectorD scatPrecSensor(2);
	TVectorD scat(2); 
	scat[0] = 0.0; scat[1]=0.0;//TO DO: This will depend on the direction of the beam.  
	//TO DO:Should make the reading for any planes not just si. 
	streamlog_out(DEBUG5)<<"The X0 for the sensor is: " << x0 <<std::endl; 
	if(x0 == 0){
		throw(lcio::Exception(Utility::outputColourString("X0 is zero for the medium inbetween the planes .", "RED"))); 	
	}
	if(getBeamEnergy() == 0 ){ 
		throw(lcio::Exception(Utility::outputColourString("Beam energy is zero", "RED"))); 	
	}
	const double scatteringVariance  = Utility::getThetaRMSHighland(getBeamEnergy(), x0);
	scatPrecSensor[0] = 1.0 / (scatteringVariance *scatteringVariance);
	scatPrecSensor[1] = 1.0 / (scatteringVariance * scatteringVariance);

	point.addScatterer(scat, scatPrecSensor);

	streamlog_out(MESSAGE1) << "  setScattererGBL  ------------- END ----------------- " << std::endl;
}
//This will add measurement information to the GBL point
//Note that if we have a strip sensor then y will be ignored using projection matrix.
void EUTelGBLFitter::setMeasurementGBL(gbl::GblPoint& point, const double *hitPos,  double statePos[3], double combinedCov[4], TMatrixD projection){
	streamlog_out(DEBUG1) << " setMeasurementGBL ------------- BEGIN --------------- " << std::endl;
 	TVectorD meas(5);//Remember we need to pass the same 5 since gbl expects this due to jacobian
	meas.Zero();
	meas[3] = hitPos[0] - statePos[0];
	meas[4] = hitPos[1] - statePos[1];
	TVectorD measPrec(5); 
	//TO DO: Can make the variance vector a matrix to explore relationships between x/y and variance.  
	measPrec[3] = 1. / combinedCov[0];	// cov(x,x)
	measPrec[4] = 1. / combinedCov[3];	// cov(y,y)
	streamlog_out(DEBUG4) << "This is what we add to the measured point:" << std::endl;
	streamlog_out(DEBUG4) << "Residuals and covariant matrix for the hit:" << std::endl;
	streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] <<"," << std::endl;
	streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20)  <<"," << measPrec[1] << std::endl;
	streamlog_out(DEBUG4) << "This H matrix:" << std::endl;
	streamlog_message( DEBUG0, projection.Print();, std::endl; );
	//The gbl library creates 5 measurement vector and 5x5 propagation matrix automatically. If  
	point.addMeasurement(projection, meas, measPrec, 0);

	streamlog_out(DEBUG1) << " setMeasurementsGBL ------------- END ----------------- " << std::endl;
}

//This function calculates the alignment jacobian and labels and attaches this to the point
void EUTelGBLFitter::setAlignmentToMeasurementJacobian(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList ){
	for (int i = 0;i<_vectorOfPairsMeasurementStatesAndLabels.size() ;++i ){
		if(_vectorOfPairsMeasurementStatesAndLabels.at(i).first.getTrackerHits().empty()){
			throw(lcio::Exception(Utility::outputColourString("One of the points on the list of measurements states has no hit.", "RED")));
		}
		for(int j = 0; j < pointList.size(); ++j){
			if( _vectorOfPairsMeasurementStatesAndLabels.at(i).second == pointList.at(j).getLabel()){
				EUTelState state = _vectorOfPairsMeasurementStatesAndLabels.at(i).first;
				streamlog_out(DEBUG0)<<"This point has "<<pointList.at(j).hasMeasurement()<<" measurements."<<std::endl;
				if(getLabelToPoint(pointList,_vectorOfPairsMeasurementStatesAndLabels.at(i).second).hasMeasurement() == 0){
					throw(lcio::Exception(Utility::outputColourString("This point does not contain a measurements. Labeling of the state must be wrong ", "RED")));
				} 
				_MilleInterface->computeAlignmentToMeasurementJacobian(state);//We calculate the jacobian. 
				_MilleInterface->setGlobalLabels(state); //Get the correct label for the sensors x,y,z shift and rotations. Depending on alignment mode and sensor the plane is on 
				TMatrixD&  alignmentJacobian = _MilleInterface->getAlignmentJacobian();//Return what was calculated by computeAlignmentToMeasurementJacobian
				std::vector<int> labels =  _MilleInterface->getGlobalParameters();//Return what was set by setGlobalLabels
				streamlog_out(DEBUG0)<<"This is the alignment jacobian we are about to add to the point at location "<<state.getLocation()<<std::endl;
				streamlog_message( DEBUG0, alignmentJacobian.Print();, std::endl; );
				streamlog_out(DEBUG0)<<"Here are the labels for these alignment parameters for the state at location: "<<state.getLocation()<<std::endl;
				for(int k=0; k<labels.size();++k){
					streamlog_out(DEBUG0)<<"Label: "<<labels.at(k) <<std::endl;
				}
				pointList.at(j).addGlobals(labels, alignmentJacobian);
				streamlog_out(DEBUG0)<<"The number of global parameters for this point are "<<pointList[i].getNumGlobals()<<std::endl;
				streamlog_out(DEBUG0)<<"The alignment matrix from GBL looks like: "<<std::endl;
				streamlog_message( DEBUG0, pointList.at(i).getGlobalDerivatives().Print() ;, std::endl; );
				streamlog_out(DEBUG0)<<"The state associated with this alignment jacobian:  "<<std::endl;
				streamlog_message( DEBUG0, state.print() ;, std::endl; );
				break;
			}
			streamlog_out(DEBUG0)<<"This point has "<<pointList.at(j).hasMeasurement()<<" measurements. It has label "<<pointList.at(j).getLabel() <<std::endl;
			if(j == (pointList.size()-1)){
				throw(lcio::Exception(Utility::outputColourString("There is no point associated with this state.", "RED")));
			}
		}
	}
}
gbl::GblPoint EUTelGBLFitter::getLabelToPoint(std::vector<gbl::GblPoint> & pointList, int label){
	for(int i = 0; i< pointList.size();++i){
		streamlog_out(DEBUG0)<<"This point has label : "<<pointList.at(i).getLabel()<<". The label number: "<<label<<std::endl;
		if(pointList.at(i).getLabel() == label){
			return pointList.at(i);
		}
		if(i == (pointList.size()-1)){ //If we reach the end of the loop and still no match then this label must belong to no point
			streamlog_out(MESSAGE5)<<"The size of pointList: "<<pointList.size()<<". The label number: "<<label<<std::endl;
			throw(lcio::Exception(Utility::outputColourString("There is no point with this label", "RED")));
		}
	}
}
//TO DO: This is the iterative alignment part that varies the precision of the measurement to allow the GBL tracks to be fitted. The precision matrix itself is an input by us. In the long run can we not do proper error calculations? 
//TO DO: This does not have to be done for each state since it only varies per plane. So could input this data another way. However in the long run this may not be favourable since we could have correct error analysis eventually. 
void EUTelGBLFitter::setMeasurementCov(EUTelState& state){
	double hitcov[]= {0,0,0,0};
	int izPlane = state.getLocation();
	if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 ){
		hitcov[0] = _parameterIdXResolutionVec[izPlane];
		hitcov[3] = _parameterIdYResolutionVec[izPlane];

		hitcov[0] *= hitcov[0]; 
		hitcov[3] *= hitcov[3];
	}
	else{//Default in case you have not specified a variance  
		throw(lcio::Exception(Utility::outputColourString("There is no measurement variance specified.", "RED"))); 	
	}
	state.setCombinedHitAndStateCovMatrixInLocalFrame(hitcov);
}
} // namespace eutelescope


#endif
