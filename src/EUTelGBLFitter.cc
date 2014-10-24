/* 
 * File:   EUTelGBLFitter.cc
 * Contact: alexander.morton975@gmail.com
 * 
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

	EUTelGBLFitter::EUTelGBLFitter() :
	_beamQ(-1),
	_eBeam(4.),
	_mEstimatorType(),
	_mille(0),
	_alignmentMode(0),
	_counter_num_pointer(1),
	_parameterIdXShiftsMap(),
	_parameterIdYShiftsMap(),
	_parameterIdZShiftsMap(),
	_parameterIdXRotationsMap(),
	_parameterIdYRotationsMap(),
	_parameterIdZRotationsMap()
	{}

	EUTelGBLFitter::~EUTelGBLFitter() {
	}
	//THIS IS THE SETTERS. Not simple get function but the action of all these functions it in the end to set a member variable to something.	
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
	//Note that we take the planes themselfs at scatters and also add scatterers to simulate the medium inbetween. 
	void EUTelGBLFitter::setScattererGBL(gbl::GblPoint& point, EUTelState & state ) {
		streamlog_out(DEBUG1) << " setScattererGBL ------------- BEGIN --------------  " << std::endl;
		TVectorD scat(2); 
		scat[0] = 0.0; scat[1]=0.0;//TO DO: This will depend on if the initial track guess has kinks. Only important if we reiterate GBL track in high radiation length enviroments  
		TMatrixDSym precisionMatrix =  state.getScatteringVarianceInLocalFrame();
		streamlog_out(MESSAGE1) << "The precision matrix being used for the sensor  "<<state.getLocation()<<":" << std::endl;
		streamlog_message( DEBUG0, precisionMatrix.Print();, std::endl; );
		point.addScatterer(scat, precisionMatrix);
		streamlog_out(DEBUG1) << "  setScattererGBL  ------------- END ----------------- " << std::endl;
		}
		//This is used when the we know the radiation length already
		void EUTelGBLFitter::setScattererGBL(gbl::GblPoint& point,EUTelState & state, float  percentageRadiationLength) {
		streamlog_out(MESSAGE1) << " setScattererGBL ------------- BEGIN --------------  " << std::endl;
		TVectorD scat(2); 
		scat[0] = 0.0; scat[1]=0.0;  
		TMatrixDSym precisionMatrix =  state.getScatteringVarianceInLocalFrame(percentageRadiationLength);
		streamlog_out(MESSAGE1) << "The precision matrix being used for the scatter:  " << std::endl;
		streamlog_message( DEBUG0, precisionMatrix.Print();, std::endl; );
		point.addScatterer(scat, precisionMatrix);
		streamlog_out(MESSAGE1) << "  setScattererGBL  ------------- END ----------------- " << std::endl;
	}
	//This will add measurement information to the GBL point
	//Note that if we have a strip sensor then y will be ignored using projection matrix.
	void EUTelGBLFitter::setMeasurementGBL(gbl::GblPoint& point, const double *hitPos,  double statePos[3], double combinedCov[4], TMatrixD projection){
		streamlog_out(DEBUG1) << " setMeasurementGBL ------------- BEGIN --------------- " << std::endl;
		TVectorD meas(2);//Remember we need to pass the same 5 since gbl expects this due to jacobian
		meas.Zero();
		meas[0] = hitPos[0] - statePos[0];//This is how residuals are shown as output after fit (measurement-prediction)
		meas[1] = hitPos[1] - statePos[1];
		TVectorD measPrec(2); 
		//TO DO: Can make the variance vector a matrix to explore relationships between x/y and variance.  
		measPrec[0] = 1. / combinedCov[0];	// cov(x,x)
		measPrec[1] = 1. / combinedCov[3];	// cov(y,y)
		streamlog_out(DEBUG4) << "This is what we add to the measured point:" << std::endl;
		streamlog_out(DEBUG4) << "Residuals and covariant matrix for the hit:" << std::endl;
		streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] <<"," << std::endl;
		streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20)  <<"," << measPrec[1] << std::endl;
		streamlog_out(DEBUG4) << "This H matrix:" << std::endl;
		streamlog_message( DEBUG0, projection.Print();, std::endl; );
		//The gbl library creates 5 measurement vector and 5x5 propagation matrix automatically. If  
		point.addMeasurement(projection, meas, measPrec, 0);//The last zero is the minimum precision before this is set to 0. TO DO:Remove this magic number
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
					streamlog_out(DEBUG0)<<"Point needs global parameters. It has  "<<pointList.at(j).hasMeasurement()<<" measurements."<<std::endl;
					if(getLabelToPoint(pointList,_vectorOfPairsMeasurementStatesAndLabels.at(i).second).hasMeasurement() == 0){
						throw(lcio::Exception(Utility::outputColourString("This point does not contain a measurements. Labeling of the state must be wrong ", "RED")));
					} 
					_MilleInterface->computeAlignmentToMeasurementJacobian(state);//We calculate the jacobian. 
					_MilleInterface->setGlobalLabels(state); //Get the correct label for the sensors x,y,z shift and rotations. Depending on alignment mode and sensor the plane is on 
					TMatrixD&  alignmentJacobian = _MilleInterface->getAlignmentJacobian();//Return what was calculated by computeAlignmentToMeasurementJacobian
					std::vector<int> labels =  _MilleInterface->getGlobalParameters();//Return what was set by setGlobalLabels
					streamlog_out(DEBUG0)<<"The state associated with this alignment jacobian:  "<<std::endl;
					streamlog_message( DEBUG0, state.print() ;, std::endl; );
					if(!state.getIsThereAHit()){
						throw(lcio::Exception(Utility::outputColourString("You are just about to use a state with no hit to determine alignment jacobian.", "RED")));
					}
					streamlog_out(DEBUG0)<<"This is the alignment jacobian we are about to add to the point at location "<<state.getLocation()<<std::endl;
					streamlog_message( DEBUG0, alignmentJacobian.Print();, std::endl; );
					streamlog_out(DEBUG0)<<"Here are the labels for these alignment parameters for the state at location: "<<state.getLocation()<<std::endl;
					for(int k=0; k<labels.size();++k){
						streamlog_out(DEBUG0)<<"Label just before adding: "<<labels.at(k) <<std::endl;
					}
					pointList.at(j).addGlobals(labels, alignmentJacobian);
					streamlog_out(DEBUG0)<<"The number of global parameters for this point are "<<pointList[i].getNumGlobals()<<std::endl;
					for(int k=0; k<pointList.at(j).getGlobalLabels().size();++k){
						streamlog_out(DEBUG0)<<"Label just after adding: "<<pointList.at(j).getGlobalLabels().at(k) <<std::endl;
					}
					streamlog_out(DEBUG0)<<"The alignment matrix after adding to point: "<<std::endl;
					streamlog_message( DEBUG0, pointList.at(j).getGlobalDerivatives().Print() ;, std::endl; );
					break;
				}
				//streamlog_out(DEBUG0)<<"This point has "<<pointList.at(j).hasMeasurement()<<" measurements. It has label "<<pointList.at(j).getLabel() <<std::endl;
				if(j == (pointList.size()-1)){
					throw(lcio::Exception(Utility::outputColourString("There is no point associated with this state.", "RED")));
				}
			}
		}
	}
	//This creates the scatters point to simulate the passage through air. 
	void EUTelGBLFitter::setPointListWithNewScatterers(std::vector< gbl::GblPoint >& pointList,EUTelState & state, float  percentageRadiationLength){
		if(_scattererJacobians.size() != _scattererPositions.size()){
			throw(lcio::Exception(Utility::outputColourString("The size of the scattering positions and jacobians is different.", "RED"))); 	
		}
		for(int i = 0 ;i < _scattererJacobians.size()-1;++i){//The last jacobain is used to get to the plane! So only loop over to (_scatterJacobians-1)
			gbl::GblPoint point(_scattererJacobians[i]);
			point.setLabel(_counter_num_pointer);
			_counter_num_pointer++;
			setScattererGBL(point,state, percentageRadiationLength/2);//TO DO:Simply dividing by 2 will work when there is only two planes
			pointList.push_back(point);
		}
	}
	//This set the estimate resolution for each plane in the X direction.
	void EUTelGBLFitter::setParamterIdXResolutionVec( const std::vector<float>& vector)
	{
		//We have a similar check after this to see that number of planes and elements in resolution vector are the same. We need this here since if they are different then it will just give an exception from the vector tryign to access a element that does not exist.
		if ( geo::gGeometry().sensorZOrdertoIDs().size() != vector.size() ){
			streamlog_out( ERROR5 ) << "The number of planes: "<< geo::gGeometry().sensorZOrdertoIDs().size()<<"  The size of input resolution vector: "<<vector.size()  << std::endl;
			throw(lcio::Exception(Utility::outputColourString("The size of the resolution vector and the total number of planes is different for x axis.", "RED")));
		}

		for(int i=0; i < geo::gGeometry().sensorZOrdertoIDs().size(); ++i){
			_parameterIdXResolutionVec[ geo::gGeometry().sensorZOrderToID(i)] = vector.at(i);
		}
	}
	//This sets the estimated resolution for each plane in the Y direction.
	void EUTelGBLFitter::setParamterIdYResolutionVec( const std::vector<float>& vector)
	{
		if ( geo::gGeometry().sensorZOrdertoIDs().size() != vector.size() ){
			streamlog_out( ERROR5 ) << "The number of planes: "<< geo::gGeometry().sensorZOrdertoIDs().size()<<"  The size of input resolution vector: "<<vector.size()  << std::endl;
			throw(lcio::Exception(Utility::outputColourString("The size of the resolution vector and the total number of planes is different for y axis.", "RED")));
		}
		for(int i=0; i < geo::gGeometry().sensorZOrdertoIDs().size(); ++i){
			_parameterIdYResolutionVec[ geo::gGeometry().sensorZOrderToID(i)] = vector.at(i);
		}
	}

	//This is used to deal with downweighting of ouliers. You must provide as input t,h or c. This specifies the function that will be used to do the downweigting.
	//TO DO: Check that this function works as expected since it has never been used. 
	void EUTelGBLFitter::setMEstimatorType( const std::string& mEstimatorType ) {
		std::string mEstimatorTypeLowerCase = mEstimatorType;
		std::transform( mEstimatorType.begin(), mEstimatorType.end(), mEstimatorTypeLowerCase.begin(), ::tolower);//Make the character lower case
		if ( mEstimatorType.size() != 1 ) {
			streamlog_out( WARNING1 ) << "More than one character supplied as M-estimator option" << std::endl;
			streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
			return;
		}
		//Compare character to the ones that are accepted and if one is the same then set out member variable equal to it.
		if ( mEstimatorTypeLowerCase.compare("t") == 0 ||
			mEstimatorTypeLowerCase.compare("h") == 0 ||
			mEstimatorTypeLowerCase.compare("c") == 0   ) this->_mEstimatorType = _mEstimatorType;
		else {
			streamlog_out( WARNING1 ) << "M-estimator option " << mEstimatorType << " was not recognized" << std::endl;
			streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
		}
	}
	//Not all points are states so we need a way to do:
	//state->label->point. This is what this function does.  
	void EUTelGBLFitter::setPointVec( std::vector< gbl::GblPoint >& pointList, gbl::GblPoint& point) {
	point.setLabel(_counter_num_pointer);
	_counter_num_pointer++;
	pointList.push_back(point);
	streamlog_out(DEBUG0) << endl << "pushBackPoint size: " << pointList.size() <<  std::endl;
	}
	//We use the labels for trajectory and equate them to the pointList index.
	//We create them when creating the points. However trajectory will overwrite these. However what we write should be the same as trajectory.
	//Note this is used by track fitting so we need to associate ANY state to the correct GBL point.
	//It is important to note the different between  setPairAnyStateAndPointLabelVec and setPairMeasurementStateAndPointLabelVec. One is for any state even if it has no hit and the other is used only for states that have a hit.
	void EUTelGBLFitter::setPairAnyStateAndPointLabelVec(std::vector< gbl::GblPoint >& pointList, gbl::GblTrajectory* traj){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairAnyStateAndPointLabelVec-- BEGIN " << endl;
		_vectorOfPairsStatesAndLabels.clear();
		streamlog_out(DEBUG5)<<"The number of states is: "<< _statesInOrder.size()<<std::endl;
		std::vector<unsigned int> labels;
		traj->getLabels(labels);
		for(int i=0; i<_statesInOrder.size();++i){
			int threeiPlus1 = 3*i; //This is done since we always have two scatters between states
			streamlog_out(DEBUG0)<<"Pair (state,label)  ("<<  &(_statesInOrder.at(i))<<","<<labels.at(threeiPlus1)<<")"<<std::endl; 
			_vectorOfPairsStatesAndLabels.push_back(make_pair(_statesInOrder.at(i), labels.at(threeiPlus1)));	
		}
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairAnyStateAndPointLabelVec- END " << endl;
	}
	//This code must be repeated since we need to create the ling between states and labels before trajectory in alignment
	//Note here we create the link between MEASUREMENT states and label. 
	void EUTelGBLFitter::setPairMeasurementStateAndPointLabelVec(std::vector< gbl::GblPoint >& pointList){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairMeasurementStateAndPointLabelVec-- BEGIN " << endl;
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
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setPairMeasurementStateAndPointLabelVec------------ END " << endl;
	}
	//As a track passes through a scatterer it will be kinked. The initial guessed trajectory has to provide GBL this information from pattern recognition. These come effectively from the states at each plane and can be calculated from these. However we store these number in the lcio file since the calculation is rather arduous
	void EUTelGBLFitter::setKinkInformationToTrack(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >& pointList,EUTelTrack &track){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setKinkInformationToTrack-- BEGIN " << endl;
		for(int i=0;i < track.getStatesPointers().size(); i++){//We get the pointers no since we want to change the track state contents		
			EUTelState* state = track.getStatesPointers().at(i);
			TVectorD corrections(5);
			TMatrixDSym correctionsCov(5,5);
			for(int j=0 ; j< _vectorOfPairsStatesAndLabels.size();++j){
				if(_vectorOfPairsStatesAndLabels.at(j).first == *state){
					streamlog_out(DEBUG0)<<"The loop number for states with measurements for kink update is: " << j << ". The label is: " << _vectorOfPairsStatesAndLabels.at(j).second <<endl; 
					if(getLabelToPoint(pointList,_vectorOfPairsStatesAndLabels.at(j).second).hasMeasurement() == 0){//TO DO: Some states will not have hits in the future so should remove. Leave for now to test
						throw(lcio::Exception(Utility::outputColourString("This point does not contain a measurements. So can not add kink information. Labeling of the state must be wrong ", "RED")));
					} 
					streamlog_out(DEBUG0)<<"To update track kink we use label: "<<_vectorOfPairsStatesAndLabels.at(j).second<<std::endl; 
					unsigned int numData; //Not sure what this is used for??????
					TVectorD aResidualsKink(2);//Measurement - Prediction
					TVectorD aMeasErrorsKink(2);
					TVectorD aResErrorsKink(2);
					TVectorD aDownWeightsKink(2); 
					traj->getMeasResults(_vectorOfPairsMeasurementStatesAndLabels.at(j).second, numData, aResidualsKink, aMeasErrorsKink, aResErrorsKink, aDownWeightsKink);
					streamlog_out(DEBUG3) << endl << "State before we have added corrections: " << std::endl;
					state->print();
					TVectorD kinks = state->getKinks();	
					TVectorD updateKinks(2);
					updateKinks(0) = kinks(0);
					updateKinks(1) = kinks(1);
					state->setKinks(updateKinks);
					streamlog_out(DEBUG3) << endl << "State after we have added corrections: " << std::endl;
					state->print();
					break;
				}
			}//END of loop of all states with hits	
		}//END of loop of all states
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::setKinkInformationToTrack-- END " << endl;
	}
	// Convert input TrackCandidates and TrackStates into a GBL Trajectory
	// This is done using the geometry setup, the scattering and the hits + predicted states.
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
			if(i != (track.getStates().size()-1)){//Since we don't want to propagate from the last state.
				nextState = track.getStates().at(i+1);
			}
			setScattererGBL(point,state);//Every sensor will have scattering due to itself. 
			_statesInOrder.push_back(state);//This is list of measurements states in the correct order. This is used later to associate ANY states with point labels
			if(state.getTrackerHits().size() == 0 ){
				streamlog_out(DEBUG3)  << Utility::outputColourString("This state does not have a hit. ", "YELLOW")<<std::endl;
				setPointVec(pointList, point);//This creates the vector of points and keeps a link between states and the points they created
			}else{
				double localPositionForState[]	= {state.getPosition()[0], state.getPosition()[1],state.getPosition()[2]};//Need this since geometry works with const doubles not floats 
				setMeasurementCov(state);
				double cov[4] ;
				state.getCombinedHitAndStateCovMatrixInLocalFrame(cov);
				setMeasurementGBL(point, state.getTrackerHits()[0]->getPosition(),  localPositionForState,  cov, state.getProjectionMatrix());
				_measurementStatesInOrder.push_back(state);//This is list of measurements states in the correct order. This is used later to associate MEASUREMENT states with point labels in alignment
				setPointVec(pointList, point);
			}//End of else statement if there is a hit.

			if(i != (track.getStates().size()-1)){//We do not produce scatterers after the last plane
				//Note here to determine the scattering we use a straight line approximation between the two points the particle will travel through. However to determine were to place the scatterer we use the exact arc length. We do this since to change the TGeo radiation length would be a lot of work for a very small change. 
				const double stateReferencePoint[] = {state.getPosition()[0], state.getPosition()[1],state.getPosition()[2]};
				double globalPosSensor1[3];
				geo::gGeometry().local2Master(state.getLocation(),stateReferencePoint , globalPosSensor1 );
				const double nextStateReferencePoint[] = {nextState.getPosition()[0], nextState.getPosition()[1],nextState.getPosition()[2]};
				double globalPosSensor2[3];
				geo::gGeometry().local2Master(nextState.getLocation(),nextStateReferencePoint , globalPosSensor2 );
				testDistanceBetweenPoints(globalPosSensor1,globalPosSensor2);
				float percentageRadiationLength  = geo::gGeometry().findRadLengthIntegral(globalPosSensor1,globalPosSensor2, false );//TO DO: This adds the radiation length of the plane again. If you chose true the some times it returns 0. The could be the reason for the slightly small residuals. 
				if(percentageRadiationLength == 0){
					streamlog_out(MESSAGE9)<<"The positions between the scatterers are: "<<endl;
					streamlog_out(MESSAGE9)<<"Start: "<<globalPosSensor1[0]<<" "<<globalPosSensor1[1]<<" "<<globalPosSensor1[2]<<endl;
					streamlog_out(MESSAGE9)<<"End: "<<globalPosSensor2[0]<<" "<<globalPosSensor2[1]<<" "<<globalPosSensor2[2]<<endl;
					throw(lcio::Exception(Utility::outputColourString("The provides radiation length to create scatterers it zero.", "RED"))); 	
				} 
				findScattersZPositionBetweenTwoStates(state);//We use the exact arc length between the two states to place the scatterers. 
				jacPointToPoint=findScattersJacobians(state,nextState);
				setPointListWithNewScatterers(pointList,state, percentageRadiationLength);//We assume that on all scattering planes the incidence angle is the same as on the last measurement state. Not a terrible approximation and will be corrected by GBL anyway.
			}else{
				streamlog_out(DEBUG3)<<"We have reached the last plane"<<std::endl;
			}
		} 
		streamlog_out(DEBUG4)<<"EUTelGBLFitter::setInformationForGBLPointList-------------------------------------END"<<endl;

	}
	//THIS IS THE GETTERS
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
	//This used after trackfit will fill a map between (sensor ID and residualx/y). 
	void EUTelGBLFitter::getResidualOfTrackandHits(gbl::GblTrajectory* traj, std::vector< gbl::GblPoint > pointList,EUTelTrack& track, map< int, map< float, float > > &  SensorResidual, map< int, map< float, float > >& sensorResidualError ){
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
			map<float, float> resError; //This is create on the stack but will pass thisa by value to the new map so it 
			resError.insert(make_pair(aResErrors[0],aResErrors[1]));
			sensorResidualError.insert(make_pair(state.getLocation(), resError));	
		}
	}
	std::string EUTelGBLFitter::getMEstimatorType( ) const {
			return _mEstimatorType;
	}
	//COMPUTE
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
	//TEST
	void EUTelGBLFitter::testUserInput(){
		if(_parameterIdXResolutionVec.size() != _parameterIdYResolutionVec.size()){
				throw(lcio::Exception(Utility::outputColourString("The vector for resolutions for X and Y are different sizes.", "RED")));
		}
		if(_parameterIdXResolutionVec.size() != geo::gGeometry().nPlanes() ){
				throw(lcio::Exception(Utility::outputColourString("The total number of planes and the resolution of the planes vector are different sizes.", "RED")));
		}
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

	//OTHER FUNCTIONS
	//We want to create a jacobain from (Plane1 -> scatterer1) then (scatterer1->scatterer2) then (scatter2->plane2). We return the last jacobain
	TMatrixD EUTelGBLFitter::findScattersJacobians(EUTelState state, EUTelState nextState){
		_scattererJacobians.clear();
		TVector3 position = state.getPositionGlobal();
		TVector3 momentum = state.computeCartesianMomentum();
		TVector3 newMomentum;
		int location = state.getLocation();
		int locationEnd = -999; //This will create a scatterer with normal in beam direction
		const gear::BField&   Bfield = geo::gGeometry().getMagneticFiled();
		gear::Vector3D vectorGlobal(position[0],position[1],position[1]);//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
		const double Bx = (Bfield.at( vectorGlobal ).x());//We times bu 0.3 due to units of other variables. See paper. Must be Tesla
		const double By = (Bfield.at( vectorGlobal ).y());
		const double Bz = (Bfield.at( vectorGlobal ).z());
		TVector3 B;
		B[0]=Bx; B[1]=By; B[2]=Bz;
		for(int i=0;i<_scattererPositions.size();i++){
			newMomentum = geo::gGeometry().getXYZMomentumfromArcLength(momentum, position,state.getBeamCharge(), _scattererPositions[i] );
			TMatrixD curvilinearJacobian = geo::gGeometry().getPropagationJacobianCurvilinear(_scattererPositions[i], state.getOmega(), momentum.Unit(),newMomentum.Unit());
			streamlog_out(DEBUG0)<<"This is the curvilinear jacobian at sensor : " << location << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, curvilinearJacobian.Print();, std::endl; );
			TMatrixD localToCurvilinearJacobianStart =  geo::gGeometry().getLocalToCurvilinearTransformMatrix(momentum, location ,state.getBeamCharge() );
			streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor : " << location << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, localToCurvilinearJacobianStart.Print();, std::endl; );
			TMatrixD localToCurvilinearJacobianEnd =  geo::gGeometry().getLocalToCurvilinearTransformMatrix(newMomentum,locationEnd ,state.getBeamCharge() );
			streamlog_out(DEBUG0)<<"This is the local to curvilinear jacobian at sensor : " << locationEnd << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, localToCurvilinearJacobianEnd.Print();, std::endl; );
			TMatrixD curvilinearToLocalJacobianEnd = localToCurvilinearJacobianEnd.Invert();
			streamlog_out(DEBUG0)<<"This is the curvilinear to local jacobian at sensor : " << locationEnd << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, curvilinearToLocalJacobianEnd.Print();, std::endl; );
			TMatrixD localToNextLocalJacobian = curvilinearToLocalJacobianEnd*curvilinearJacobian*localToCurvilinearJacobianStart;
			streamlog_out(DEBUG0)<<"This is the full jacobian : " << locationEnd << " or scatter: "<< i << std::endl; 
			streamlog_message( DEBUG0, localToNextLocalJacobian.Print();, std::endl; );
			_scattererJacobians.push_back(localToNextLocalJacobian);//To DO if scatter then plane is always parallel to z axis
			momentum[0]=newMomentum[0]; momentum[1]=newMomentum[1];	momentum[2]=newMomentum[2];
			location = -999;//location will always be a scatter after first loop 
			if(i == (_scattererPositions.size()-2)){//On the last loop we want to create the jacobain to the next plane
				locationEnd = nextState.getLocation();
			}
		}
		if(_scattererJacobians.size() != 3){
			throw(lcio::Exception(Utility::outputColourString("There are not 3 jacobians produced by scatterers!.", "RED"))); 	
		}
		return _scattererJacobians.back();//return the last jacobian so the next state can use this
	}
	//The distance from the first state to the next scatterer and then from that scatterer to the next all the way to the next state. 
	//TO DO: This uses the optimum positions as described by Claus.  However for non homogeneous material distribution this might not be the case.
	void EUTelGBLFitter::findScattersZPositionBetweenTwoStates(EUTelState& state){
		streamlog_out(DEBUG1) << "  findScattersZPositionBetweenTwoStates------------- BEGIN --------------  " << std::endl;
		_scattererPositions.clear();	
		float arcLength = state.getArcLengthToNextState();
		streamlog_out(DEBUG1) << "The arc length to the next state is: " << arcLength << std::endl;
		float distance1 =arcLength/2 - arcLength/sqrt(12); 
		_scattererPositions.push_back(distance1);//Z position of 1st scatter	
		float distance2 = arcLength/2 + arcLength/sqrt(12);//Z position of 2nd scatter 
		_scattererPositions.push_back(distance2);
		_scattererPositions.push_back(arcLength); 
		streamlog_out(DEBUG1) << "  findScattersZPositionBetweenTwoStates------------- END --------------  " << std::endl;
	}
	void EUTelGBLFitter::resetPerTrack() {
		_counter_num_pointer=1;
		_measurementStatesInOrder.clear();
		_statesInOrder.clear();

	}
	//This function will take the estimate track from pattern recognition and add a correction to it. This estimated track + correction is you final GBL track.
	void EUTelGBLFitter::updateTrackFromGBLTrajectory (gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >& pointList,EUTelTrack &track, map<int, vector<double> > &  mapSensorIDToCorrectionVec){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::UpdateTrackFromGBLTrajectory-- BEGIN " << endl;
		for(int i=0;i < track.getStatesPointers().size(); i++){//We get the pointers no since we want to change the track state contents		
			EUTelState* state = track.getStatesPointers().at(i);
			TVectorD corrections(5);
			TMatrixDSym correctionsCov(5,5);
			for(int j=0 ; j< _vectorOfPairsStatesAndLabels.size();++j){
				if(_vectorOfPairsStatesAndLabels.at(j).first == *state){
					streamlog_out(DEBUG0)<<"The loop number for states with measurements is: " << j << ". The label is: " << _vectorOfPairsStatesAndLabels.at(j).second <<endl; 
		//				if(getLabelToPoint(pointList,_vectorOfPairsStatesAndLabels.at(j).second).hasMeasurement() == 0){//TO DO: Some states will not have hits in the future so should remove. Leave for now to test
		//					throw(lcio::Exception(Utility::outputColourString("This point does not contain a measurements. Labeling of the state must be wrong ", "RED")));
		//				} 
					streamlog_out(DEBUG0)<<"To update track we use label: "<<_vectorOfPairsStatesAndLabels.at(j).second<<std::endl; 
					traj->getResults(_vectorOfPairsStatesAndLabels.at(j).second, corrections, correctionsCov );
					streamlog_out(DEBUG3) << endl << "State before we have added corrections: " << std::endl;
					state->print();
					TVectorD newStateVec(5);
					newStateVec[0] = state->getOmega() + corrections[0];
					newStateVec[1] = state->getIntersectionLocalXZ()+corrections[1];
					newStateVec[2] = state->getIntersectionLocalYZ()+corrections[2];
					newStateVec[3] = state->getPosition()[0]+corrections[3];
					newStateVec[4] = state->getPosition()[1]+corrections[4]; 
					//Here we collect the total correction for every track and state. This we use to work out an average to make sure that on each iteration the track corrections are reducing
					_omegaCorrections = _omegaCorrections + corrections[0];	
					_intersectionLocalXZCorrections= _intersectionLocalXZCorrections+ corrections[1];	
					_intersectionLocalYZCorrections = _intersectionLocalYZCorrections + corrections[2];	
					_localPosXCorrections = _localPosXCorrections + corrections[3];
					_localPosYCorrections = _localPosYCorrections + corrections[4];

					std::vector<double> correctionVec;
					correctionVec.push_back(corrections[0]);	
					correctionVec.push_back(corrections[1]);	
					correctionVec.push_back(corrections[2]);	
					correctionVec.push_back(corrections[3]);	
					correctionVec.push_back(corrections[4]);	

					mapSensorIDToCorrectionVec[state->getLocation()] = correctionVec;
					state->setStateVec(newStateVec);
					streamlog_out(DEBUG3) << endl << "State after we have added corrections: " << std::endl;
					state->print();
					break;
				}
			}//END of loop of all states with hits	
		}//END of loop of all states
	}

} // namespace eutelescope


#endif
