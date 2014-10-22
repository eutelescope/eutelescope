//This class contains all the pattern recognition functions. These are used within EUTelProcessorPatternRecogntion. They use basic input about the particles flight. So the particle charge, the magnetic field strength and beam energy. This combined with seeds states on the first plane (Can be others but the first by default) determined what set of hits come from a single track. 
#include "EUTelPatternRecognition.h"
namespace eutelescope {
    
	EUTelPatternRecognition::EUTelPatternRecognition() :  
	_allowedMissingHits(0),
	_AllowedSharedHitsOnTrackCandidate(0),
	_totalNumberOfHits(0),
	_totalNumberOfSharedHits(0),
	_numberOfTracksTotal(0),
	_numberOfTracksAfterHitCut(0),
	_numberOfTracksAfterPruneCut(0),
	_firstExecution(true),
	_beamE(-1.),
	_beamQ(-1.)
	{}

	EUTelPatternRecognition::~EUTelPatternRecognition() { 
	}
 
std::vector<EUTelTrack>& EUTelPatternRecognition::getTracks(){
	return	_finalTracks; 
}
void EUTelPatternRecognition::clearEveryRun(){
	int _totalNumberOfHits=0;
	int _totalNumberOfSharedHits=0;
}


void EUTelPatternRecognition::testTrackCandidates(){
	for(int i=0; i < _tracks.size(); ++i){
		int idBefore=-999;
		for(int j = 0; j<_tracks.at(i).getTracks().size();++j){
			if(_tracks.at(i).getTracks().at(j)->getTrackerHits().size() > 1){
				throw(lcio::Exception( Utility::outputColourString("The number of hits for each state is greater than 1","RED")));
			}
			if(!_tracks.at(i).getTracks().at(j)->getTrackerHits().empty()){//Since some states will have not hits
				const EVENT::TrackerHit* hit =  _tracks.at(i).getTracks().at(j)->getTrackerHits().at(0);
				int id = hit->id();
				if(j>1 and id == idBefore){
					streamlog_out(MESSAGE5) << "The IDs of the hits are: " <<id <<" and before  "<<idBefore<<std::endl; 
					//note z0 stores the position of the state.
					streamlog_out(MESSAGE5) << "The state locations are " << _tracks.at(i).getTracks().at(j)->getZ0() <<" and  "<<_tracks.at(i).getTracks().at(j-1)->getZ0()<<std::endl; 
					throw(lcio::Exception( Utility::outputColourString("Some states have the same hits. ","RED")));
				}
				idBefore=id;
			}
		}
	}
}
//This is the work horse of the class. Using seeds it propagates the track forward using equations of motion. This can be with or without magnetic field.
void EUTelPatternRecognition::propagateForwardFromSeedState( EUTelState& stateInput, EUTelTrack & track    ){
	EUTelState *state = &stateInput;//Make it a pointer so we can change this to newState after.
	streamlog_out ( DEBUG1 ) << "EUTelPatternRecognition::propagateForwardFromSeedState-----BEGIN "<< endl;
	//TO DO: To delete this in smart way. 
	EUTelState *firstState = new EUTelState(stateInput);//Need to create an initial state that will not be deleted outside this scope  
	streamlog_out(DEBUG2) << "This is the memory location of the state: "<< firstState << std::endl;
	track.addTrack(static_cast<EVENT::Track*>(firstState));//Note we do not have to create new since this object State is saved in class member scope
	//Here we loop through all the planes not excluded. We begin at the seed which might not be the first. Then we stop before the last plane, since we do not want to propagate anymore
	bool firstLoop =true;//TO DO:: This works but is very stupid. Must fix
	for(int i = geo::gGeometry().sensorIDToZOrderWithoutExcludedPlanes()[state->getLocation()]; i < (geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()-1); ++i){
	
		float globalIntersection[3];
		TVector3 momentumAtIntersection;
		float arcLength;
		int newSensorID = state->findIntersectionWithCertainID(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i+1], globalIntersection, momentumAtIntersection, arcLength);
		if(arcLength == 0 or arcLength < 0 ){ 
			throw(lcio::Exception( Utility::outputColourString("The arc length is less than or equal to zero. ","RED"))); 
		}
		if(firstLoop){
			firstState->setArcLengthToNextState(arcLength); 
			firstLoop =false;
		}else{
			state->setArcLengthToNextState(arcLength);
		}
		//cout<<"HERE1: "<<state->getPosition()[0]<<","<<state->getPosition()[1]<<","<<state->getPosition()[2]<<","<<state->getLocation()<<endl;
		int sensorIntersection = geo::gGeometry( ).getSensorID(globalIntersection);
		if(newSensorID < 0 or sensorIntersection < 0 ){
			streamlog_out ( MESSAGE5 ) << "Intersection point on infinite plane: " <<  globalIntersection[0]<<" , "<<globalIntersection[1] <<" , "<<globalIntersection[2]<<std::endl;
			streamlog_out ( MESSAGE5 ) << "Momentum on next plane: " <<  momentumAtIntersection[0]<<" , "<<momentumAtIntersection[1] <<" , "<<momentumAtIntersection[2]<<std::endl;
			streamlog_out(MESSAGE5) <<" From ID= " <<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i]<< " to " <<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i+1]  <<std::endl;
			streamlog_out(MESSAGE5)<<"Was there intersection on plane: "<<newSensorID<<" Was there intersection in sensitive area: "<< sensorIntersection <<std::endl;
			streamlog_out(MESSAGE5) << Utility::outputColourString("No intersection found moving on plane. Move to next plane and look again.","YELLOW")<<std::endl; 
			streamlog_out(MESSAGE5)<<"This is for event number " <<getEventNumber()<<std::endl;
			continue;//So if there is no intersection look on the next plane. Important since two planes could be at the same z position
		}
		streamlog_out ( DEBUG0 ) << "Intersection point on infinite plane: " <<  globalIntersection[0]<<" , "<<globalIntersection[1] <<" , "<<globalIntersection[2]<<std::endl;
		streamlog_out ( DEBUG0 ) << "Momentum on next plane: " <<  momentumAtIntersection[0]<<" , "<<momentumAtIntersection[1] <<" , "<<momentumAtIntersection[2]<<std::endl;

		//So we have intersection lets create a new state
		EUTelState *newState = new EUTelState();//Need to create this since we save the pointer and we would be out of scope when we leave this function. Destroying this object. 
		newState->setDimensionSize(_planeDimensions[newSensorID]);//We set this since we need this information for later processors
		newState->setBeamCharge(_beamQ);
		newState->setLocation(newSensorID);
		newState->setPositionGlobal(globalIntersection);
		newState->setLocalXZAndYZIntersectionAndCurvatureUsingGlobalMomentum(momentumAtIntersection);
		if(_mapHitsVecPerPlane[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i+1]].size() == 0){
			streamlog_out(DEBUG5) << Utility::outputColourString("There are no hits on the plane with this state. Add state to track as it is and move on ","YELLOW");
			track.addTrack(static_cast<EVENT::Track*>(newState));//Need to return this to LCIO object. Loss functionality but retain information 
			state = newState;
			continue;
		}
		EVENT::TrackerHit* closestHit = const_cast< EVENT::TrackerHit* > ( findClosestHit( *newState )); //This will look for the closest hit but not if it is within the excepted range		
		const double* hitPosition = closestHit->getPosition();
		const double distance = sqrt(computeResidual(*newState, closestHit ).Norm2Sqr()); //Determine the residual of state and hit  //Distance is in mm. Norm2Sqr does not square root for some reason
		const double DCA = getXYPredictionPrecision( *newState ); //This does nothing but return a number specified by user. In the future this should use convariance matrix information TO DO: FIX
		streamlog_out ( DEBUG1 ) <<"At plane: "<<newState->getLocation() << ". Distance between state and hit: "<< distance <<" Must be less than: "<<DCA<< endl;
		streamlog_out(DEBUG0) <<"Closest hit position: " << closestHit->getPosition()[0]<<" "<< closestHit->getPosition()[1]<<"  "<< closestHit->getPosition()[2]<<endl;
		if ( distance > DCA ) {
			streamlog_out ( DEBUG1 ) << "Closest hit is outside of search window." << std::endl;
			track.addTrack(static_cast<EVENT::Track*>(newState));//Need to return this to LCIO object. Loss functionality but retain information 
			state = newState;
			continue;
		}	
		if(closestHit == NULL){
			throw(lcio::Exception( Utility::outputColourString("The closest hit you are trying to add is NULL. This can not be correct","RED")));
		}
		streamlog_out ( DEBUG1 ) << "Found a hit with memory address: " << closestHit<<" and ID of " <<closestHit->id() <<" At a Distance: "<< distance<<" from state." << endl;
		newState->addHit(closestHit);
		_totalNumberOfHits++;//This is used for test of the processor later.   
		streamlog_out(DEBUG2) << "This is the memory location of the state: "<< newState << std::endl;

		track.addTrack(static_cast<EVENT::Track*>(newState));//Need to return this to LCIO object. Loss functionality but retain information 
		track.print();
		streamlog_out ( DEBUG1 ) << "The number of hits on the track now is "<< track.getNumberOfHitsOnTrack()<< endl;
		//cout<<"HERE3: "<<newState->getPosition()[0]<<","<<newState->getPosition()[1]<<","<<newState->getPosition()[2]<<","<<newState->getLocation()<<endl;
		//cout<<"HERE5: "<<state->getPosition()[0]<<","<<state->getPosition()[1]<<","<<state->getPosition()[2]<<","<<state->getLocation()<<endl;
		state = newState;
		//cout<<"HERE4: "<<state->getPosition()[0]<<","<<state->getPosition()[1]<<","<<state->getPosition()[2]<<","<<state->getLocation()<<endl;
		streamlog_out ( DEBUG1 ) << "End of loop "<< endl;

	}
	streamlog_out ( DEBUG1 ) << "EUTelPatternRecognition::propagateForwardFromSeedState-----END "<< endl;

}	
void EUTelPatternRecognition::printTrackCandidates(){
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::printTrackCandidates----BEGIN "<< endl;
	for(int i = 0; i < _tracks.size();++i){
		streamlog_out(DEBUG5)<<"Track number "<<i<<" Out of " <<_tracks.size()<<std::endl; 
		_tracks.at(i).print();
	}
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::printTrackCandidates----END "<< endl;
}
void EUTelPatternRecognition::clearFinalTracks(){
	if(!_finalTracks.empty()){
		_finalTracks.clear();
	}
}
//It is important to note what the output of this tells you. If your get that 25% of tracks passed pruning this does not mean that only 25% of events had a track. It simply means that out of all tracks that began from seed hit. Only 25% made a full track in the end.
void EUTelPatternRecognition::testTrackQuality(){
	_numberOfTracksTotal = _numberOfTracksTotal + _tracks.size();
	_numberOfTracksAfterHitCut = _numberOfTracksAfterHitCut + _tracksAfterEnoughHitsCut.size();
	_numberOfTracksAfterPruneCut =_numberOfTracksAfterPruneCut + 	_finalTracks.size();

	if(_numberOfTracksTotal % 1000 == 0){
		streamlog_out(MESSAGE5)<<"///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
		streamlog_out(MESSAGE5)<<" Total Tracks "<<" Pass Hit Cut " <<"                           Pass Prune Cut " <<std::endl;
		streamlog_out(MESSAGE5)<<"    "<<_numberOfTracksTotal<<"    "<<"    "<<_numberOfTracksAfterHitCut<<"                                                     "<<_numberOfTracksAfterPruneCut<<std::endl;
		float percentAfterHitCut = (static_cast<float>(_numberOfTracksAfterHitCut)/static_cast<float>(_numberOfTracksTotal))*100;
		float percentAfterPruneCut = (static_cast<float>(_numberOfTracksAfterPruneCut)/static_cast<float>(_numberOfTracksTotal))*100;
		streamlog_out(MESSAGE5)<<"               "<<percentAfterHitCut<<"                                                      "<<percentAfterPruneCut<<std::endl;
		float averageNumberOfHitsOnTrack = static_cast<float>(_totalNumberOfHits)/static_cast<float>(_numberOfTracksTotal);
		float averageNumberOfSharedHitsOnTrack = static_cast<float>(_totalNumberOfSharedHits)/static_cast<float>( _numberOfTracksTotal);
		streamlog_out(MESSAGE5)<<"The average number of hits on a track: "<< averageNumberOfHitsOnTrack<<std::endl;
		streamlog_out(MESSAGE5)<<"The average number of shared hits on a track with other tracks: "<< averageNumberOfSharedHitsOnTrack<<std::endl;//Note this can change with cust since we remove tracks before we have counted all similar hits. To see all similar hits make cut very large and then run
		if(percentAfterHitCut < 0.1){
			streamlog_out(MESSAGE5)<< Utility::outputColourString("The percentage of track making the hit cut is very low at the moment ","YELLOW")<<std::endl;
		}
		if(percentAfterPruneCut < 0.1){
			streamlog_out(MESSAGE5)<< Utility::outputColourString("The percentage of track making the prune cut is very low at the moment ","YELLOW")<<std::endl;
		}	
	}
}
//I have observed a difference of just over 6 microns between the two methods over 20 cm distance. So we will set the maximum distance to 10 microns difference
void EUTelPatternRecognition::testPositionEstimation(float position1[], float position2[]){
	streamlog_out(DEBUG0) << " The distance between the jacobain methods and simple equations of motion are for the x direction " <<  position1[0] -position2[0] << std::endl;
	if(fabs(position1[0] -position2[0]) > 0.01){ //in mm
		throw(lcio::Exception( Utility::outputColourString("The positions between the two methods is different in the x direction. ","RED"))); 
	}
	streamlog_out(DEBUG0) << " The distance between the jacobain methods and simple equations of motion are for the y direction " <<  position1[1] -position2[1] << std::endl;
	if(fabs(position1[1] -position2[1]) > 0.01){ //in mm
		throw(lcio::Exception( Utility::outputColourString("The positions between the two methods is different in the y direction. ","RED"))); 
	}
	streamlog_out(DEBUG0) << " The distance between the jacobain methods and simple equations of motion are for the z  direction " <<  position1[2] -position2[2] << std::endl;
	if(fabs(position1[2] -position2[2]) > 0.01){ //in mm
		throw(lcio::Exception( Utility::outputColourString("The positions between the two methods is different in the z direction. ","RED"))); 
	} 
}

     // Perform track pruning this removes tracks that have the same hits used to create the track on some planes
		//TO DO: consider a better approach. Since we could remove tracks that have a better residual to hits just because that have came first. For example was will always add the final track in the list. Could this introduce bias? 
		//TO DO:We compare all states to all other states on a track. We don't need to do that since hits on differents planes must be different.
void EUTelPatternRecognition::findTrackCandidatesWithSameHitsAndRemove(){
	streamlog_out(MESSAGE1) << "EUTelPatternRecognition::findTrackCandidatesWithSameHitsAndRemove----BEGIN" << std::endl;
	for(int i =0; i < _tracksAfterEnoughHitsCut.size();++i){//LOOP through all tracks 
		streamlog_out(DEBUG1) <<  "Loop at track number: " <<  i <<". Must loop over " << _tracksAfterEnoughHitsCut.size()<<" tracks in total."   << std::endl;
		std::vector<EUTelState> iStates = _tracksAfterEnoughHitsCut.at(i).getStates();
		//Now loop through all tracks one ahead of the original track itTrk. This is done since we want to compare all the track to each other to if they have similar hits     
		for(int j =i+1; j < _tracksAfterEnoughHitsCut.size();++j){ //LOOP over all track again.
	//		cout<<"Increase track to: "<<j <<endl;
			int hitscount=0;
			std::vector<EUTelState> jStates = _tracksAfterEnoughHitsCut[j].getStates();
			for(int i=0;i<iStates.size();i++){
		//		cout<<"I top states "<< i<<endl;
			
				EVENT::TrackerHit* ihit;
				if(!iStates[i].getTrackerHits().empty()){//Need since we could have tracks that have a state but no hits here.
					ihit = iStates[i].getTrackerHits()[0];
				}else{
					continue;
				}
				int ic = ihit->id();
				for(int j=0;j<jStates.size();j++){
		//		cout<<"I bottom states "<< j<<endl;

					EVENT::TrackerHit* jhit;
					if(!jStates[j].getTrackerHits().empty()){//Need since we could have tracks that have a state but no hits here.
						jhit = jStates[j].getTrackerHits()[0];
					}else{
						continue;
					}
					int jc = jhit->id();
					if(ic == jc ){
				//		cout<<"This is ic and jc " <<ic<<","<<jc<<endl;
						_totalNumberOfSharedHits++;
						hitscount++; 
						streamlog_out(MESSAGE1) <<  "Hit number on track you are comparing all other to :" << i << ". Hit ID: " << ic << ". Hit number of comparison: " << j << ". Hit ID of this comparison : " << jc << ". Number of common hits: " << hitscount << std::endl; 
					}
				}

			} 
			//If for any track we are comparing to the number of similar hits is to high then we move to the next track and do not add this track to the new list of tracks
	//		cout<<"Number of hits: " <<hitscount<<" allowed: " <<_AllowedSharedHitsOnTrackCandidate<<endl;
			if(hitscount > _AllowedSharedHitsOnTrackCandidate) {   
				streamlog_out(DEBUG1)<<"Tracks has too many similar hits remove "<<std::endl;
				break;
			}
			if(j == (_tracksAfterEnoughHitsCut.size()-1)){//If we have loop through all and not breaked then track must be good.
	//			cout<<"Enter the addition of track"<<endl;
				_finalTracks.push_back(_tracksAfterEnoughHitsCut[i]);
				streamlog_out(DEBUG1)<<"Track made prune tracks cut"<<std::endl;
			}
		}
		//We need to add the last track here since the inner loop j+1 will never be entered. We always add the last track since if it has similar hits to past tracks then those tracks have been removed.
		if(i == (_tracksAfterEnoughHitsCut.size()-1)){//If we have loop through all and not breaked then track must be good.
//		cout<<"Enter the last addition of track"<<endl;
		_finalTracks.push_back(_tracksAfterEnoughHitsCut[i]);
	//	streamlog_out(DEBUG1)<<"Track made prune tracks cut"<<std::endl;
		}

	}
	streamlog_out(MESSAGE1) << "------------------------------EUTelPatternRecognition::findTrackCandidatesWithSameHitsAndRemove()---------------------------------END" << std::endl;
}
//This function is used to check that the input data is as expect. If not then we end the processor by throwing a exception.
////TO DO: should check if seed planes are also excluded	
void EUTelPatternRecognition::testUserInput() {
	streamlog_out(DEBUG2) << "EUTelPatternRecognition::testUserInput()" << std::endl;
	if ( _beamE < 1.E-6 ) {
		throw(lcio::Exception( Utility::outputColourString("Beam direction was set incorrectly","RED"))); 
	}
	else{
	 streamlog_out(DEBUG0) << Utility::outputColourString("Beam energy is reasonable", "GREEN") << std::endl;
	}
	if(_createSeedsFromPlanes.size() == 0){
		throw(lcio::Exception( Utility::outputColourString("The number of planes to make seeds from is 0. We need at least one plane", "RED")));
	}

	if(_createSeedsFromPlanes.size() >= geo::gGeometry().sensorIDstoZOrder().size()){
		throw(lcio::Exception( Utility::outputColourString("You have specified all planes or more than that. This is too many planes for seeds", "RED")));
	}else{
		streamlog_out(DEBUG0) << Utility::outputColourString("The number of planes to make seeds from is good " +  to_string(_createSeedsFromPlanes.size()),"GREEN") << std::endl;
	}
	testPlaneDimensions();
}	
void  EUTelPatternRecognition::testHitsVec(){
	if(_allHitsVec.size() == 0){
		throw(lcio::Exception( "The hit input is wrong!"));
	}

}

/** 
*   Using a list of planes use the hits on those planes to create the initial seeds. 
*   These seeds are used then to propagate to the next plane.
* 
*/
//TO DO:This will not work with a strip sensor creating a seed. Since the strip sensor creates a line of possible hit positions not just a single point. Therefore you need to loop through  all points on the line and determine its trajectory and then see if the intersection is close enough to the hit on the next plane to add that hit on the track. This is a lot of work for something which is not that important.   
//There is no way round this since the strip give a range of possible positions in 2D/3D space.      //To allow strip sensor to be seed would involve major change to how this is done. For now just don't do it  

void EUTelPatternRecognition::initialiseSeeds() {
	streamlog_out(DEBUG2) << "EUTelPatternRecognition::initialiseSeeds()" << std::endl;
	_mapSensorIDToSeedStatesVec.clear();
	for( int iplane = 0; iplane < _createSeedsFromPlanes.size() ; iplane++) {

		streamlog_out(DEBUG1) << "We are using plane: " <<  _createSeedsFromPlanes[iplane] << " to create seeds" << std::endl;

		EVENT::TrackerHitVec& hitFirstLayer = _mapHitsVecPerPlane[_createSeedsFromPlanes[iplane]];
		streamlog_out(DEBUG1) << "N hits on sensor : " << hitFirstLayer.size() << std::endl;
		if(hitFirstLayer.size()== 0){
			continue;
		}
		std::vector<EUTelState> stateVec;
		EVENT::TrackerHitVec::iterator itHit;
		for ( itHit = hitFirstLayer.begin(); itHit != hitFirstLayer.end(); ++itHit ) {
			EUTelState state;//Here we create a track state. This is a point on a track that we can give a position,momentum and direction. We combine these to create a track. 
			double temp[] = { (*itHit)->getPosition()[0], (*itHit)->getPosition()[1], (*itHit)->getPosition()[2] };
			float posLocal[] = { static_cast<float>(temp[0]), static_cast<float>(temp[1]), static_cast<float>(temp[2]) };
			if( posLocal[2] != 0){//TO DO:This should be in test hits
				streamlog_out(MESSAGE5) << "The local position of the hit is " << posLocal[0]<<","<< posLocal[1]<<","<< posLocal[2]<<","<<std::endl;
				throw(lcio::Exception(Utility::outputColourString("The position of this local coordinate in the z direction is non 0", "RED"))); 	
			}
			state.setLocation(_createSeedsFromPlanes[iplane]);//This stores the location as a float in Z0 since no location for track LCIO. This is preferable to problems with storing hits.  
			state.setPositionLocal(posLocal); //This will automatically take cartesian coordinate system and save it to reference point. //This is different from most LCIO applications since each track will have own reference point. 		
			float xzDirection;
			float yzDirection;
			state.setBeamCharge(_beamQ);//this is set for each state. to do: is there a more efficient way of doing this since we only need this stored once?
			TVector3 momentum = computeInitialMomentumGlobal(); 
			state.setLocalXZAndYZIntersectionAndCurvatureUsingGlobalMomentum(momentum); 
			state.addHit(*itHit);
			_totalNumberOfHits++;//This is used for test of the processor later.   
			state.setDimensionSize(_planeDimensions[state.getLocation()]);
			stateVec.push_back(state);
		}
		_mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[iplane]] = stateVec; 
	}
	streamlog_out(DEBUG2) << "--------------------------------EUTelPatternRecognition::initialiseSeeds()---------------------------" << std::endl;
}
TVector3 EUTelPatternRecognition::computeInitialMomentumGlobal(){
	//We assume that the arc length is the displacement in the z direction. The assumption should be a valid one in most cases
	TVector3 position(0,0,0);//The position we start from does not matter since the magnetic field is homogeneous.
	TVector3 momentum(0,0,_beamE);//Assume the beam starts in a straight line
	float arcLength= geo::gGeometry().getInitialDisplacementToFirstPlane();
	TVector3 momentumEnd = geo::gGeometry().getXYZMomentumfromArcLength(momentum, position, _beamQ, arcLength);
	streamlog_out(DEBUG2) << "Momentum on the first sensor: px,py,pz "<<momentumEnd[0] <<","<<momentumEnd[1]<<","<<momentumEnd[2]<<","<<"At an arc length of "<<arcLength<<std::endl;
	return momentumEnd;
}
void EUTelPatternRecognition::testInitialSeeds(){
	//Even if there are no hits then we can compare these. Since it will contain stateVec with size 0.      
//	if(_mapSensorIDToSeedStatesVec.size() != _createSeedsFromPlanes.size()){TO DO:Fix strange error with some event gives no seed planes
	//	streamlog_out(MESSAGE5) <<"The size of sensors with seeds: " << _mapSensorIDToSeedStatesVec.size() <<" The number of planes you are suppose to create seeds from: "<<_createSeedsFromPlanes.size()<< " For event number: " <<getEventNumber()<<std::endl;
//		throw(lcio::Exception(Utility::outputColourString("The size of intial state seeds planes and the number to use at the start are different", "RED"))); 	
//	}
	for(int i = 0 ; i < _mapSensorIDToSeedStatesVec.size(); ++i){
		std::vector<EUTelState> StatesVec =  _mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[i]]; 	
		for(int j = 0 ; j < StatesVec.size() ; ++j){
			if(StatesVec[j].getTrackerHits()[0] == NULL ){
				throw(lcio::Exception(Utility::outputColourString("The hit is NULL. All seeds must have hits.", "RED"))); 	
			}
		}
	}
}
//Loop through each plane that contains seeds and then each seed. From that seed you then create a track.
void EUTelPatternRecognition::findTrackCandidates() {
	streamlog_out(MESSAGE1) << "EUTelPatternRecognition::findTrackCandidates()" << std::endl;
	clearTrackAndTrackStates(); //Clear all past track information
	for(int i = 0 ; i < _mapSensorIDToSeedStatesVec.size(); ++i){
		std::vector<EUTelState> statesVec =  _mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[i]]; 	
		if(statesVec.size() == 0){
			streamlog_out(MESSAGE5) << Utility::outputColourString("The size of state Vector seeds is zero. try next seed plane","YELLOW")<<std::endl; 
			continue;
		}
		for(int j = 0 ; j < statesVec.size() ; ++j){
			EUTelTrack track;
			propagateForwardFromSeedState(statesVec[j], track);
			streamlog_out ( DEBUG1 ) << "Before adding track to vector "<< endl;
			_tracks.push_back(track);//Here we create a long list of possible tracks
			streamlog_out ( DEBUG1 ) << "After adding track to vector "<< endl;


		}
	}
	streamlog_out(MESSAGE1) << "EUTelPatternRecognition::findTrackCandidates()------END" << std::endl;
}

//We need to delete the states AND the track. Since we have allocated memory to store the state.
void EUTelPatternRecognition::clearTrackAndTrackStates(){
	for(int i=0; i < _tracks.size(); ++i){
		for(int j=0; j< _tracks.at(i).getTracks().size(); ++j){   
		//	delete [] _tracks.at(i).getTracks().at(j); TO DO: Delete memory properly 
		}
	}
	_tracks.clear();
}

void EUTelPatternRecognition::findTracksWithEnoughHits(){
	streamlog_out(DEBUG1) << "EUTelPatternRecognition::findTracksWithEnoughHits()------BEGIN" << std::endl;
	_tracksAfterEnoughHitsCut.clear();
	if(_tracks.size() == 0 ){
		streamlog_out(MESSAGE5) <<"This is event: " <<getEventNumber()<<std::endl;   
		streamlog_out(MESSAGE5) << Utility::outputColourString("The number of tracks for this event is zero ","YELLOW")<<std::endl; 
	}
	for(int i = 0 ; i<_tracks.size(); ++i){
		EUTelTrack& track = _tracks[i];
		streamlog_out ( DEBUG2 ) << "Number of hits on the track: " <<track.getNumberOfHitsOnTrack()<<" Number needed: " <<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - _allowedMissingHits << std::endl;
		if(track.getNumberOfHitsOnTrack()>= geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - _allowedMissingHits){
			streamlog_out(DEBUG5) << "There are enough hits. So attach this track makes the cut!"<<endl;
			_tracksAfterEnoughHitsCut.push_back(track);
		}
	}
	streamlog_out(DEBUG1) << "EUTelPatternRecognition::findTracksWithEnoughHits()-----END" << std::endl;

}
//This creates map between plane ID and hits on that plane. 
//We also order the map correcly with geometry.
void EUTelPatternRecognition::setHitsVecPerPlane(){
	streamlog_out(DEBUG0) <<"EUTelPatternRecognition::setHitsVecPerPlane()----------------------------BEGIN" <<std::endl;
	_mapHitsVecPerPlane.clear();
	int numberOfPlanes = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size();//Note should not make this a class data member since we call this by reference in geometry so defacto it is already accessed directly each time. By this I mean we do not create a new copy each time we call the function. 
	streamlog_out(DEBUG0) <<"The number of planes that we will add hits to: "<< numberOfPlanes  <<std::endl;
	if(numberOfPlanes == 0){
		throw(lcio::Exception( "The number of planes is 0 to get hits from."));
	}
	if( _allHitsVec.size()== 0){
		throw(lcio::Exception( "The number of hits is zero."));
	}
	for(int i=0 ; i<numberOfPlanes;++i){
		EVENT::TrackerHitVec tempHitsVecPlaneX; 
		for(int j=0 ; j<_allHitsVec.size();++j){
			if(Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*>(_allHitsVec[j]) ) ==  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)){
				tempHitsVecPlaneX.push_back(_allHitsVec.at(j));
			}
		}		
	_mapHitsVecPerPlane[ geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)] = 	tempHitsVecPlaneX;
	}	
	streamlog_out(DEBUG0) <<"EUTelPatternRecognition::setHitsVecPerPlane()----------------------------END" <<std::endl;
}
//Note loop through all planes. Even the excluded. This is easier since you don't have to change this input each time then.
//TO DO: sensors z position as used here only works if there is sufficient difference between planes. In the order of 1mm in gear file. This is too large.
void EUTelPatternRecognition::setPlaneDimensionsVec(EVENT::IntVec planeDimensions){
	streamlog_out(DEBUG0) <<"EUTelPatternRecognition::setPlaneDimensionsVec()----------------------------BEGIN" <<std::endl;
	if(planeDimensions.size() != geo::gGeometry().sensorZOrdertoIDs().size()){
		streamlog_out(ERROR) << "The size of planesDimensions input is: "<< planeDimensions.size()<<" The size of sensorZOrdertoIDs is: " << geo::gGeometry().sensorZOrdertoIDs().size()<< std::endl;
		throw(lcio::Exception( Utility::outputColourString("The input dimension vector not the same as the number of planes!","RED")));
	}
	_planeDimensions.clear();
	for(int i=0; i<geo::gGeometry().sensorZOrdertoIDs().size(); ++i){//Loop through each plane
		int planeID = geo::gGeometry().sensorZOrdertoIDs().at(i);
		if (_planeDimensions.find(planeID) == _planeDimensions.end()){//This is to check that we don't try to map the same sensor to two different plane dimensions.
			_planeDimensions[planeID] = planeDimensions.at(i);
		}else{
			streamlog_out(ERROR5) <<"The z position is : "<< i <<" The plane ID is: " << planeID <<std::endl;
			throw(lcio::Exception( Utility::outputColourString("You are trying to map the same sensor ID to two different plane dimensions. There is something wrong with you gear file input. Make sure there is some distance between your planes in the gear file!","RED")));
		}
	}//END of loop over planes
	streamlog_out(DEBUG0) <<"EUTelPatternRecognition::setPlaneDimensionsVec()----------------------------END" <<std::endl;
}	    


void EUTelPatternRecognition::testHitsVecPerPlane(){
	streamlog_out(DEBUG4) <<"EUTelPatternRecognition::testHitsVecPerPlane----------------------------BEGIN" <<std::endl;
	if(_mapHitsVecPerPlane.size() !=  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()){
		streamlog_out(ERROR0) <<"The size of the planes with hits " << _mapHitsVecPerPlane.size() <<"  Sensors from Geometry with no excluded planes: "<<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()<<std::endl;
		throw(lcio::Exception(Utility::outputColourString("The number of planes that could contain hits and the number of planes is different. Problem could be the gear file has to planes at the same z position.", "RED"))); 	

	}
	for(int i=0 ;i<_mapHitsVecPerPlane.size();++i){
		int planeID = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i);
		if(_mapHitsVecPerPlane.at(planeID).size() <= 0){
			streamlog_out(DEBUG0) << Utility::outputColourString("One plane has no hits at all. Is this correct?","YELLOW") << std::endl;
		}
	}
	streamlog_out(DEBUG4) <<"EUTelPatternRecognition::testHitsVecPerPlane----------------------------END" <<std::endl;

}
//This function makes sure that the input for the dimension size is correct.
//NOTE:This input is alway for ALL planes and not just no excluded ones. This means we don't have to change this if we want to exclude planes. 
void EUTelPatternRecognition::testPlaneDimensions(){
	streamlog_out(DEBUG4) <<"EUTelPatternRecognition::testPlaneDimensions----------------------------BEGIN" <<std::endl;
	if(_planeDimensions.size() != geo::gGeometry().sensorZOrdertoIDs().size()){
		streamlog_out(ERROR5) << "The size of _planesDimensions is: "<< _planeDimensions.size()<<" The size of sensorZOrdertoIDs is: " << geo::gGeometry().sensorZOrdertoIDs().size()<< std::endl;
		throw(lcio::Exception( Utility::outputColourString("The output dimension vector is not the same size as the number of planes. Could be something to do with the gear file. Make sure you have some distances between you planes!","RED")));
	}
	for(int i=0;i<_planeDimensions.size();++i){
		if(_planeDimensions.at(geo::gGeometry().sensorZOrdertoIDs().at(i))>2 or _planeDimensions.at(geo::gGeometry().sensorZOrdertoIDs().at(i))<=0){
			throw(lcio::Exception( "The number of dimension for one of your planes is greater than 2 or less than 0. If this is not a mistake collect you nobel prize now!"));
		}
	}
	streamlog_out(DEBUG4) <<"EUTelPatternRecognition::testPlaneDimensions----------------------------END" <<std::endl;
}

    /** Find the hit closest to the intersection of a track with given sensor
     * 
     * @param ts track state
     * @return hit closest to the intersection of the track with the sensor plane
     * 
	 */
const EVENT::TrackerHit* EUTelPatternRecognition::findClosestHit(EUTelState & state ) {
	streamlog_out(DEBUG2) << "EUTelPatternRecognition::findClosestHit()" << std::endl;
	EVENT::TrackerHitVec& hitInPlane = _mapHitsVecPerPlane[state.getLocation()];
	double maxDistance = std::numeric_limits<double>::max();
	EVENT::TrackerHitVec::const_iterator itClosestHit;
	EVENT::TrackerHitVec::const_iterator itHit;
	streamlog_out(DEBUG0) << "N hits in plane " << state.getLocation() << ": " << hitInPlane.size() << std::endl;
	for ( itHit = hitInPlane.begin(); itHit != hitInPlane.end(); ++itHit ) {
		const double distance = computeResidual( state, *itHit ).Norm2Sqr();//This distance could be 2D or 1D depending on if you have a strip or pixel sensor.
		streamlog_out(DEBUG0) << "Distance^2 between hit and track intersection: " << distance << std::endl;
		if ( distance < maxDistance ) {
			itClosestHit = itHit;
			maxDistance = distance;
		}
	}
	streamlog_out(DEBUG0) << "Minimal distance^2 between hit and track intersection: " << maxDistance << std::endl;
	streamlog_out(DEBUG2) << "----------------------EUTelPatternRecognition::findClosestHit()------------------------" << std::endl;

	return *itClosestHit;
}
//This is not very useful at the moment since the covariant matrix for the hit is guess work at the moment.   
double EUTelPatternRecognition::getXYPredictionPrecision(EUTelState& ts ) const {
	streamlog_out(DEBUG2) << "EUTelPatternRecognition::getXYPredictionPrecision()---BEGIN" << std::endl;
	// TO DO: Need proper error analysis to calculate this rather than providing an answer.  
	double xyPrec = getWindowSize();   //sqrt( Ckkm1[0][0]*Ckkm1[0][0] + Ckkm1[1][1]*Ckkm1[1][1] );
	streamlog_out(DEBUG0) << "Minimal combined UV resolution : " << xyPrec << std::endl;
	streamlog_out(DEBUG2) << "----------------------EUTelPatternRecognition::getXYPredictionPrecision()------------------------END" << std::endl;
	return xyPrec;
}

        
/** Calculate residual vector between given track and hit
* 
* @param ts track state
* @param hit hit
* @return 
*/
TVectorD EUTelPatternRecognition::computeResidual(EUTelState& state , const EVENT::TrackerHit* hit ) const {
	streamlog_out( DEBUG2 ) << "EUTelPatternRecognition::computeResidual()---BEGIN" << std::endl;
	const double* hitPosition = hit->getPosition();//In local coordinates

	streamlog_out( DEBUG3 ) << "Hit (id=" << hit->id() << ") local(u,v) coordinates of hit: (" << hitPosition[0] << "," << hitPosition[1] <<","<<hitPosition[2] << ")" << std::endl;

	double localPosition [3];
	localPosition[0]=state.getPosition()[0];	localPosition[1]=state.getPosition()[1];	localPosition[2]=state.getPosition()[2];
	streamlog_out( DEBUG3 ) << "	Prediction for hit (id=" << hit->id() << ") local(u,v) coordinates of state: ("  << localPosition[0] << "," << localPosition[1] <<","<<localPosition[2] << ")" << std::endl;

	TVectorD residual(2);
	residual[0] = 0 ; residual[1] = 0;
	streamlog_out( DEBUG2 ) << "The size of this state dimension: "<< _planeDimensions.at(state.getLocation())<< std::endl;
	//This loop is used in case we have a strip sensor. So we should use only 1 dimension of information
	for(int i=0; i<_planeDimensions.at(state.getLocation());++i){
		residual[i] = hitPosition[i] - localPosition[i];
	}	
	
	if ( streamlog_level(DEBUG2) ){
			streamlog_out( DEBUG2 ) << "	Residual vector residual: (" <<residual[0] <<","<<residual[1]<<")"<< std::endl;
			residual.Print();
	}
	
	streamlog_out( DEBUG2 ) << "----------------------------------EUTelPatternRecognition::computeResidual()------------------------------------END" << std::endl;
	
	return residual;
}

	
void EUTelPatternRecognition::findHitsOrderVec(LCCollection* lcCollection,EVENT::TrackerHitVec& hitsOrderVec) {
	for (int iHit = 0; iHit < lcCollection->getNumberOfElements(); iHit++) {
		TrackerHitImpl * hit = static_cast<TrackerHitImpl*> (lcCollection->getElementAt(iHit));

		//TO DO:The determination of HitID is inefficient. Since you create a decoder for each event.
		const int localSensorID = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (hit) );
		if ( localSensorID >= 0 ) hitsOrderVec.push_back( hit );

	} // end loop over all hits in lcCollection

}
void EUTelPatternRecognition::printHits(){
	streamlog_out(MESSAGE0) << "EUTelPatternRecognition::prinitHit: BEGIN ==============" << std::endl;
	EVENT::TrackerHitVec::const_iterator itHit;
	for ( itHit = _allHitsVec.begin() ; itHit != _allHitsVec.end(); ++itHit ) {
		const double* uvpos = (*itHit)->getPosition();
		const int sensorID = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHit) );
		streamlog_out(MESSAGE0) << "Hit (id=" << setw(3) << sensorID << ") local(u,v) coordinates: ("<< setw(7) << setprecision(4) << uvpos[0] << "," << setw(7) << setprecision(4) << uvpos[1] << ")" << std::endl;
	}
streamlog_out(MESSAGE0) << "EUTelPatternRecognition::printHits: END ==============" << std::endl;
}
    
} // namespace eutelescope
