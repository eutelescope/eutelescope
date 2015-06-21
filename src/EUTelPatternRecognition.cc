//This class contains all the pattern recognition functions. These are used within EUTelProcessorPatternRecogntion. 
//They use basic input about the particles flight. So the particle charge, the magnetic field strength and beam energy. 
//This combined with seeds states on the first plane (Can be others but the first by default) determined what set 
//of hits come from a single track. 
#include "EUTelPatternRecognition.h"
#include "EUTelNav.h"

namespace eutelescope {

EUTelPatternRecognition::EUTelPatternRecognition():  
_totalNumberOfHits(0),
_totalNumberOfSharedHits(0),
_firstExecution(true),
_numberOfTracksTotal(0),
_numberOfTracksAfterHitCut(0),
_numberOfTracksAfterPruneCut(0),
_allowedMissingHits(0),
_AllowedSharedHitsOnTrackCandidate(0),
_beamE(-1.),
_beamQ(-1.)
{}
EUTelPatternRecognition::~EUTelPatternRecognition()  
{}


std::vector<EUTelTrack>& EUTelPatternRecognition::getTracks()
{
	return	_finalTracks; 
}

//void EUTelPatternRecognition::testTrackCandidates(){
//	for(size_t i=0; i < _tracks.size(); ++i){
//		int idBefore=-999;
//		for(size_t j = 0; j<_tracks.at(i).getStates().size();++j){
//			if(_tracks.at(i).getStates().at(j)->getIsThereAHit()){//Since some states will have not hits
//				EUTelHit hit =  _tracks.at(i).getStates().at(j)->getHit();
//				int id = hit->id();
//				if(j>1 and id == idBefore){
//					streamlog_out(MESSAGE5) << "The IDs of the hits are: " <<id <<" and before  "<<idBefore<<std::endl; 
//					//note z0 stores the position of the state.
//					streamlog_out(MESSAGE5) << "The state locations are " << _tracks.at(i).getTracks().at(j)->getZ0() <<" and  "<<_tracks.at(i).getTracks().at(j-1)->getZ0()<<std::endl; 
//					throw(lcio::Exception( "Some states have the same hits. "));
//				}
//				idBefore=id;
//			}
//		}
//	}
//}

//This is the work horse of the class. Using seeds it propagates the track forward using equations of motion. This can be with or without magnetic field.
void EUTelPatternRecognition::propagateForwardFromSeedState(EUTelState& stateInput, EUTelTrack& track)
{
    streamlog_out ( DEBUG1 ) << "Initial Seed: "<< std::endl;
    stateInput.print();

	std::map<const int,double>  mapSensor;
	std::map<const int ,double>  mapAir;
	double rad =	stateInput.computeRadLengthsToEnd(mapSensor, mapAir);
	EUTelState state = stateInput;
	
	if(rad == 0 ){ //If the estimated radiation length is 0 then we do not use the track.
		return;
	}
	//Here we loop through all the planes not excluded. We begin at the seed which might not be the first. Then we stop before the last plane, since we do not want to propagate anymore
	bool firstLoop =true;//This is needed so we get the arclength to the next state on the first. Completing the state and adding.
    bool calcDirection =true;// This is used to have an estimate of the direction of the particle using the first two hits associated together. 
    EUTelState newState; //Must exist after exiting loop. 
	for(size_t i = geo::gGeometry().sensorIDToZOrderWithoutExcludedPlanes().at(state.getLocation()); i < (geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()); ++i){
        //Loop one more than the last plane to add the last plane on the next loop then end 
        if(i == geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()-1){
            //Do not add arclength to the last plane 
            track.setState(newState); 
            newState.clear();
            streamlog_out ( DEBUG1 ) << "ADD STATES (Add scattering later)  Event number: "<<getEventNumber() << std::endl;
            streamlog_out ( DEBUG1 ) << "Finished track!"<< std::endl;
            track.print();
            break;
        }

	
		float globalIntersection[3];
		TVector3 momentumAtIntersection;
		float arcLength;
		int newSensorID = 0;
		bool foundNextIntersection = state.findIntersectionWithCertainID(	geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i+1), 
											globalIntersection, momentumAtIntersection, arcLength, newSensorID);

		if(!foundNextIntersection)
		{
			streamlog_out(DEBUG5) 	<< "INTERSECTION NOT FOUND! Intersection point on infinite plane: " 
						<<  globalIntersection[0] << ", " <<globalIntersection[1] << ", " << globalIntersection[2] << std::endl
						<< "Momentum on next plane: " 
						<<  momentumAtIntersection[0] << ", " << momentumAtIntersection[1] << ", " << momentumAtIntersection[2] << std::endl
						<< "From ID: " << geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i) << " to " 
						<<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i+1) << std::endl
						<< "This is for event number: " << getEventNumber() << std::endl;
			//So if there is no intersection look on the next plane.
			//Important since two planes could be at the same z position
			continue;
		}

		streamlog_out(DEBUG5) 	<< "INTERSECTION FOUND! From ID: " << geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)
					<< " to " << geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i+1) << std::endl
					<< "Intersection point on infinite plane: " 
					<<  globalIntersection[0] << ", " << globalIntersection[1] << ", " << globalIntersection[2] << std::endl
					<< "Momentum on next plane: " 
					<<  momentumAtIntersection[0]<< ", "<<momentumAtIntersection[1] << ", " << momentumAtIntersection[2] << std::endl;
        streamlog_out ( DEBUG1 ) << "ADD STATES (Add scattering later)  Event number: "<<getEventNumber() << std::endl;
		if(firstLoop)//This is where we add the state with hit collected on first entry to the loop.
		{
			state.setArcLengthToNextState(arcLength); 
            track.setState(state);//All we need to add to the first state is the arclength to the next plane.
            streamlog_out ( DEBUG1 ) << "Begin by adding first state... "<< std::endl;
            track.print();
			firstLoop =false;
		}
		else
		{
			newState.setArcLengthToNextState(arcLength);//This variable must be calculated one loop after all other variables.
			track.setState(newState); 
            newState.clear();
            streamlog_out ( DEBUG1 ) << "Update track to..."<< std::endl;
            track.print();
		}
		//So we have intersection lets create a new state
		newState.setDimensionSize(_planeDimensions[newSensorID]);//We set this since we need this information for later processors
		newState.setLocation(newSensorID);
		newState.setPositionGlobal(globalIntersection);
		newState.setLocalDirGlobalDir(momentumAtIntersection);
        state = newState;//Set state here ready to propagate. It does not need hit information to do this.

		if(_mapHitsVecPerPlane[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i+1)].empty()){
			streamlog_out(DEBUG5) << "There are no hits on the plane with this state. Add state to track as it is and move on." << std::endl;
//			track.setState(newState); 
//			state = newState;
			continue;
		}
		EVENT::TrackerHit* closestHit = const_cast<EVENT::TrackerHit*>( findClosestHit(newState)); //This will look for the closest hit but not if it is within the excepted range		
		double distance;
		if(newState.getDimensionSize() == 2)
		{
			//This distance could be 2D or 1D depending on if you have a strip or pixel sensor. 
			//Norm2Sqr does not square root for some reason.
			distance = sqrt(computeResidual( newState, closestHit ).Norm2Sqr());
		}
		else if(newState.getDimensionSize() == 1)
		{
			//If strip sensor then use only displacement along strips, which should be x axis.
			distance = computeResidual( newState, closestHit )[0];
		}
		else
		{
			throw(lcio::Exception( "The closest hit is not on a pixel or strip sensor. Since the dimensionality if less than 1 or greater than 2."));
		}

		const double DCA = getXYPredictionPrecision(newState);

		streamlog_out ( DEBUG1 ) <<"At plane: "<<newState.getLocation() << ". Distance between state and hit: "<< distance <<" Must be less than: "<<DCA<< std::endl;
		streamlog_out(DEBUG0) <<"Closest hit position: " << closestHit->getPosition()[0]<<" "<< closestHit->getPosition()[1]<<"  "<< closestHit->getPosition()[2]<<std::endl;

		if ( distance > DCA ) {
			streamlog_out ( DEBUG1 ) << "Closest hit is outside of search window." << std::endl;
//			track.setState(newState); 
//			state = newState;
			continue;
		}	
		if(closestHit == NULL){
			throw(lcio::Exception( "The closest hit you are trying to add is NULL. This can not be correct"));
		}

		streamlog_out ( DEBUG1 ) << "Found a hit with memory address: " << closestHit<<" and ID of " <<closestHit->id() <<" At a Distance: "<< distance<<" from state." << std::endl;
		newState.setHit(closestHit);
        if(calcDirection){//Here we update the direction of the track from now on. This will be correct when we have the full track.
            TVector3 momGlobal = getGlobalMomBetweenStates(stateInput,newState);
          //  momGlobal.Print();
            state.setLocalDirGlobalDir(momGlobal);
            calcDirection=false;
        }
		_totalNumberOfHits++;//This is used for test of the processor later.   
//		streamlog_out(DEBUG2) << "This is the memory location of the state: "<< newState << std::endl;

//		track.setState(newState);//Need to return this to LCIO object. Loss functionality but retain information 
//		track.print();
		streamlog_out ( DEBUG1 ) << "The number of hits on the track now is "<< track.getNumberOfHitsOnTrack()<< std::endl;
//		state = newState;
		streamlog_out ( DEBUG1 ) << "End of loop "<< std::endl;

	}
	//NOW WE ASSOCIATE THE STATES TO A SCATTERING LENGTH.

	setRadLengths(track, mapSensor, mapAir, rad);
    streamlog_out ( DEBUG1 ) << "ADD SCATTERING TO TRACKS: "<< std::endl;
    track.print();

}
//setRadLengths: This will determine the variance fraction each scatterer will get. Note this comes in two parts. The first is the plane and the next scattering from the air.    
void EUTelPatternRecognition::setRadLengths(EUTelTrack & track,	std::map<const int,double>  mapSensor, std::map<const int ,double>  mapAir, double rad ){
	//THE FINAL WEIGHT WE HAVE WILL BE A FRACTION PERCENTAGE OF THE TOTAL RADIATION LENGTH
	std::vector<EUTelState>& states = track.getStates();
	const double var  = pow( Utility::getThetaRMSHighland(track.getBeamEnergy(), rad) , 2);
	for(size_t i =0; i < track.getStates().size();++i){ //LOOP over all track again.
		streamlog_out(DEBUG0)<< std::scientific << " Values placed in variance using Highland formula corrected. (SENSOR) : " << (mapSensor[states.at(i).getLocation()]/rad)*var << "  (AIR)  " << (mapAir[states.at(i).getLocation()]/rad)*var <<std::endl;
		states.at(i).setRadFrac((mapSensor[states.at(i).getLocation()]/rad)*var,(mapAir[states.at(i).getLocation()]/rad)*var);//We input the fraction percentage.
	}
	//NOW DETERMINE THE VARIANCE DUE TO THE RADIATION LENGTH. THIS IN THE END WILL BE DIVIDED AMOUNG THE SCATTERERS.
	track.setTotalVariance(var);
}


void EUTelPatternRecognition::printTrackCandidates(){
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::printTrackCandidates----BEGIN "<< std::endl;
	for(size_t i = 0; i < _tracks.size();++i){
		streamlog_out(DEBUG5)<<"Track number "<<i<<" Out of " <<_tracks.size()<<std::endl; 
		_tracks.at(i).print();
	}
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::printTrackCandidates----END "<< std::endl;
}
void EUTelPatternRecognition::clearFinalTracks(){
	_finalTracks.clear();
}

//It is important to note what the output of this tells you. If your get that 25% of tracks passed pruning,
//this does not mean that only 25% of events had a track. It simply means that out of all tracks that began 
//from seed hit, only 25% made a full track in the end.
void EUTelPatternRecognition::testTrackQuality()
{
	_numberOfTracksTotal = _numberOfTracksTotal + _tracks.size();
	_numberOfTracksAfterHitCut = _numberOfTracksAfterHitCut + _tracksAfterEnoughHitsCut.size();
	_numberOfTracksAfterPruneCut =_numberOfTracksAfterPruneCut + 	_finalTracks.size();

	if(_numberOfTracksTotal % 1000 == 0)
	{
		float percentAfterHitCut = (static_cast<float>(_numberOfTracksAfterHitCut)/static_cast<float>(_numberOfTracksTotal))*100;
		float percentAfterPruneCut = (static_cast<float>(_numberOfTracksAfterPruneCut)/static_cast<float>(_numberOfTracksTotal))*100;
		float averageNumberOfHitsOnTrack = static_cast<float>(_totalNumberOfHits)/static_cast<float>(_numberOfTracksTotal);
		float averageNumberOfSharedHitsOnTrack = static_cast<float>(_totalNumberOfSharedHits)/static_cast<float>( _numberOfTracksTotal);
		streamlog_out(MESSAGE5) << "//////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl
				<< "Total Tracks: " << _numberOfTracksTotal << " Pass Hit Cut: " <<  _numberOfTracksAfterHitCut << " Pass Prune Cut: " << _numberOfTracksAfterPruneCut << std::endl
				<< "Percentage after Hit Cut: " << percentAfterHitCut << " Percentage after Prune Cut: " << percentAfterPruneCut << std::endl
				<< "The average number of hits on a track: " << averageNumberOfHitsOnTrack << std::endl
				<< "The average number of shared hits on a track with other tracks: " << averageNumberOfSharedHitsOnTrack << std::endl;
				//(Note this can change with cust since we remove tracks before we have counted all similar hits. To see all similar hits make cut very large and then run)
	
		if(percentAfterHitCut < 0.1) streamlog_out(DEBUG5)<< "The percentage of track making the hit cut is very low at the moment "<<std::endl;
		if(percentAfterPruneCut < 0.1)streamlog_out(DEBUG5)<< "The percentage of track making the prune cut is very low at the moment "<<std::endl;
	}
}

//I have observed a difference of just over 6 microns between the two methods over 20 cm distance. So we will set the maximum distance to 10 microns difference
void EUTelPatternRecognition::testPositionEstimation(float position1[], float position2[]){
	streamlog_out(DEBUG0) << " The distance between the jacobain methods and simple equations of motion are for the x direction " <<  position1[0] -position2[0] << std::endl;
	if(fabs(position1[0] -position2[0]) > 0.01){ //in mm
		throw(lcio::Exception( "The positions between the two methods is different in the x direction. ")); 
	}
	streamlog_out(DEBUG0) << " The distance between the jacobain methods and simple equations of motion are for the y direction " <<  position1[1] -position2[1] << std::endl;
	if(fabs(position1[1] -position2[1]) > 0.01){ //in mm
		throw(lcio::Exception( "The positions between the two methods is different in the y direction. ")); 
	}
	streamlog_out(DEBUG0) << " The distance between the jacobain methods and simple equations of motion are for the z  direction " <<  position1[2] -position2[2] << std::endl;
	if(fabs(position1[2] -position2[2]) > 0.01){ //in mm
		throw(lcio::Exception( "The positions between the two methods is different in the z direction. ")); 
	} 
}

     // Perform track pruning this removes tracks that have the same hits used to create the track on some planes
		//TO DO: consider a better approach. Since we could remove tracks that have a better residual to hits just because that have came first. For example was will always add the final track in the list. Could this introduce bias? 
		//TO DO:We compare all states to all other states on a track. We don't need to do that since hits on differents planes must be different.
void EUTelPatternRecognition::findTrackCandidatesWithSameHitsAndRemove(){
	streamlog_out(MESSAGE1) << "EUTelPatternRecognition::findTrackCandidatesWithSameHitsAndRemove----BEGIN" << std::endl;
	for(size_t i =0; i < _tracksAfterEnoughHitsCut.size();++i){//LOOP through all tracks 
		streamlog_out(DEBUG1) <<  "Loop at track number: " <<  i <<". Must loop over " << _tracksAfterEnoughHitsCut.size()<<" tracks in total."   << std::endl;
		std::vector<EUTelState> iStates = _tracksAfterEnoughHitsCut.at(i).getStates();
		//Now loop through all tracks one ahead of the original track itTrk. This is done since we want to compare all the track to each other to if they have similar hits     
		for(size_t j =i+1; j < _tracksAfterEnoughHitsCut.size();++j){ //LOOP over all track again.
			int hitscount=0;
			std::vector<EUTelState> jStates = _tracksAfterEnoughHitsCut[j].getStates();
			for(size_t k=0;k<iStates.size();k++)
			{
					EUTelHit ihit;
					//Need since we could have tracks that have a state but no hits here.
					if(iStates.at(k).getStateHasHit())
					{
							ihit = iStates[k].getHit();
					}
					else
					{
							continue;
					}
					int ic = ihit.getID();
					
					for(size_t l=0;l<jStates.size();l++)
					{
							EUTelHit jhit;
							//Need since we could have tracks that have a state but no hits here.
							if(jStates.at(l).getStateHasHit())
							{
									jhit = jStates.at(l).getHit();
							}
							else
							{
									continue;
							}
							int jc = jhit.getID();
							if(ic == jc )
							{
									_totalNumberOfSharedHits++;
									hitscount++; 
									streamlog_out(MESSAGE1) <<  "Hit number on track you are comparing all other to: " << i << ". Hit ID: " 
															<< ic << ". Hit number of comparison: " << j << ". Hit ID of this comparison : " << jc 
															<< ". Number of common hits: " << hitscount << std::endl; 
							}
					}

			} 
			//If for any track we are comparing to the number of similar hits is to high then we move to the next track and do not add this track to the new list of tracks
			if(hitscount > _AllowedSharedHitsOnTrackCandidate) {   
				streamlog_out(DEBUG1)<<"Tracks has too many similar hits remove "<<std::endl;
				break;
			}
			if(j == (_tracksAfterEnoughHitsCut.size()-1)){//If we have loop through all and not breaked then track must be good.
				_finalTracks.push_back(_tracksAfterEnoughHitsCut[i]);
				streamlog_out(DEBUG1)<<"Track made prune tracks cut"<<std::endl;
			}
		}
		//We need to add the last track here since the inner loop j+1 will never be entered. We always add the last track since if it has similar hits to past tracks then those tracks have been removed.
		if(i == (_tracksAfterEnoughHitsCut.size()-1)){//If we have loop through all and not breaked then track must be good.
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
		throw(lcio::Exception( "Beam direction was set incorrectly")); 
	}
	else{
	 streamlog_out(DEBUG0) << "Beam energy is reasonable" << std::endl;
	}
	if(_createSeedsFromPlanes.size() == 0){
		throw(lcio::Exception( "The number of planes to make seeds from is 0. We need at least one plane"));
	}

	if(_createSeedsFromPlanes.size() >= geo::gGeometry().sensorIDstoZOrder().size()){
		throw(lcio::Exception( "You have specified all planes or more than that. This is too many planes for seeds"));
	}else{
		streamlog_out(DEBUG0) << "The number of planes to make seeds from is good " +  to_string(_createSeedsFromPlanes.size()) << std::endl;
	}
	testPlaneDimensions();
}	

/** 
*   Using a list of planes use the hits on those planes to create the initial seeds. 
*   These seeds are used then to propagate to the next plane.
* 
*/
//TO DO:This will not work with a strip sensor creating a seed. Since the strip sensor creates a line of possible hit positions not just a single point. Therefore you need to loop through  all points on the line and determine its trajectory and then see if the intersection is close enough to the hit on the next plane to add that hit on the track. This is a lot of work for something which is not that important.   
//There is no way round this since the strip give a range of possible positions in 2D/3D space.      //To allow strip sensor to be seed would involve major change to how this is done. For now just don't do it  

void EUTelPatternRecognition::initialiseSeeds()
{
	_mapSensorIDToSeedStatesVec.clear();
	for( size_t iplane = 0; iplane < _createSeedsFromPlanes.size(); iplane++) 
	{
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
				throw(lcio::Exception("The position of this local coordinate in the z direction is non 0")); 	
			}
			state.setLocation(_createSeedsFromPlanes[iplane]);  
			state.setPositionLocal(posLocal);  		
			TVector3 momentum = computeInitialMomentumGlobal(); 
			state.setLocalDirGlobalDir(momentum); 
			state.setHit(*itHit);
			_totalNumberOfHits++;//This is used for test of the processor later.   
			state.setDimensionSize(_planeDimensions[state.getLocation()]);
			stateVec.push_back(state);

		}
        streamlog_out(DEBUG0) << "States to create tracks from: "<< std::endl;       
        for(size_t i = 0; i<stateVec.size(); i++){
            stateVec.at(i).print();
        }
		_mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[iplane]] = stateVec; 
	}
}
TVector3 EUTelPatternRecognition::computeInitialMomentumGlobal(){
	//We assume that the arc length is the displacement in the z direction. The assumption should be a valid one in most cases
	TVector3 position(0,0,0);//The position we start from does not matter since the magnetic field is homogeneous.
	TVector3 momentum(0,0,_beamE);//Assume the beam starts in a straight line
	float arcLength= geo::gGeometry().getInitialDisplacementToFirstPlane();
	TVector3 momentumEnd = EUTelNav::getMomentumfromArcLength(momentum, _beamQ, arcLength);
	streamlog_out(DEBUG2) << "Momentum on the first sensor: px,py,pz "<<momentumEnd[0] <<","<<momentumEnd[1]<<","<<momentumEnd[2]<<","<<"At an arc length of "<<arcLength<<std::endl;
	return momentumEnd;
}
void EUTelPatternRecognition::testInitialSeeds(){
	//Even if there are no hits then we can compare these. Since it will contain stateVec with size 0.      
//	if(_mapSensorIDToSeedStatesVec.size() != _createSeedsFromPlanes.size()){TO DO:Fix strange error with some event gives no seed planes
	//	streamlog_out(MESSAGE5) <<"The size of sensors with seeds: " << _mapSensorIDToSeedStatesVec.size() <<" The number of planes you are suppose to create seeds from: "<<_createSeedsFromPlanes.size()<< " For event number: " <<getEventNumber()<<std::endl;
//		throw(lcio::Exception("The size of intial state seeds planes and the number to use at the start are different")); 	
//	}
	for(size_t i = 0 ; i < _mapSensorIDToSeedStatesVec.size(); ++i){
		std::vector<EUTelState> StatesVec =  _mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[i]]; 	
		for(size_t j = 0 ; j < StatesVec.size() ; ++j){
			if(!StatesVec[j].getStateHasHit()){
				throw(lcio::Exception("The hit is not on first state. All seeds must have hits.")); 	
			}
		}
	}
}
//Loop through each plane that contains seeds and then each seed. From that seed you then create a track.
void EUTelPatternRecognition::findTrackCandidates() {
	streamlog_out(MESSAGE1) << "EUTelPatternRecognition::findTrackCandidates()" << std::endl;
	clearTrackAndTrackStates(); //Clear all past track information
	for(size_t i = 0 ; i < _mapSensorIDToSeedStatesVec.size(); ++i){
		std::vector<EUTelState> statesVec =  _mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[i]]; 	
		if(statesVec.size() == 0){
			streamlog_out(MESSAGE5) << "The size of state Vector seeds is zero. try next seed plane"<<std::endl; 
			continue;
		}
		for(size_t j = 0 ; j < statesVec.size() ; ++j){
			EUTelTrack track;
			propagateForwardFromSeedState(statesVec[j], track);
//			streamlog_out ( DEBUG1 ) << "Before adding track to vector "<< std::endl;
			_tracks.push_back(track);//Here we create a long list of possible tracks
//			streamlog_out ( DEBUG1 ) << "After adding track to vector "<< std::endl;


		}
	}
	streamlog_out(MESSAGE1) << "EUTelPatternRecognition::findTrackCandidates()------END" << std::endl;
}

//We need to delete the states AND the track. Since we have allocated memory to store the state.
void EUTelPatternRecognition::clearTrackAndTrackStates(){
	for(size_t i=0; i < _tracks.size(); ++i){
		for(size_t j=0; j< _tracks.at(i).getStates().size(); ++j){   
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
		streamlog_out(MESSAGE5) << "The number of tracks for this event is zero "<<std::endl; 
	}
	for(size_t i = 0 ; i<_tracks.size(); ++i){
		EUTelTrack& track = _tracks[i];
		streamlog_out ( DEBUG2 ) << "Number of hits on the track: " <<track.getNumberOfHitsOnTrack()<<" Number needed: " <<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - _allowedMissingHits << std::endl;
		if(track.getNumberOfHitsOnTrack() >= (int)geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - _allowedMissingHits){
			streamlog_out(DEBUG5) << "There are enough hits. So attach this track makes the cut!"<<std::endl;
			_tracksAfterEnoughHitsCut.push_back(track);
		}
	}
	streamlog_out(DEBUG1) << "EUTelPatternRecognition::findTracksWithEnoughHits()-----END" << std::endl;

}
//This creates map between plane ID and hits on that plane. 
//We also order the map correcly with geometry.
void EUTelPatternRecognition::setHitsVecPerPlane()
{
	_mapHitsVecPerPlane.clear();
	int numberOfPlanes = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size();
	
	if(numberOfPlanes == 0)
	{
		throw(lcio::Exception( "The number of planes is 0 to get hits from."));
	}
	if( _allHitsVec.size()== 0)
	{
		throw(lcio::Exception( "The number of hits is zero."));
	}

	for(int i=0 ; i<numberOfPlanes;++i)
	{
		EVENT::TrackerHitVec tempHitsVecPlaneX; 
		for(size_t j=0 ; j<_allHitsVec.size();++j)
		{
			if(Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*>(_allHitsVec[j]) ) ==  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i))
			{
				tempHitsVecPlaneX.push_back(_allHitsVec.at(j));
			}
		}
		_mapHitsVecPerPlane[ geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)] = 	tempHitsVecPlaneX;
	}	
}

std::vector<EUTelTrack> EUTelPatternRecognition::getSeedTracks(){
    std::vector<EUTelTrack> seededTracks;
    std::vector<EUTelTrack> tracksOriginal = _finalTracks; //Must make copy here so we do not change _finalTracks
    for(size_t i=0 ; i <tracksOriginal.size(); i++ ){
//        streamlog_out(DEBUG1) <<"Track before seed running:  "  <<std::endl;
//        tracksOriginal.at(i).print();
        bool found=true;
        EUTelTrack trackOut;
        found = seedTrackOuterHits(tracksOriginal.at(i), trackOut );
        if(found){
            streamlog_out(DEBUG1) <<"Track before seed:  "  <<std::endl;
            _finalTracks.at(i).print();
            streamlog_out(DEBUG1) <<"Track after seed:  "  <<std::endl;
            trackOut.print();
            seededTracks.push_back(trackOut);
        }

    }
    return seededTracks;
}
bool EUTelPatternRecognition::seedTrackOuterHits(EUTelTrack track,EUTelTrack & trackOut){
    //Deterimine last state with hit//
    int lastStateWithHit=0;
    for(size_t i=0 ; i < track.getStates().size() ; i++){
      if(track.getStates().at(i).getStateHasHit()){
          lastStateWithHit=i;
      }
    }

    EUTelState firstState = track.getStates().at(0);
    EUTelState lastState = track.getStates().at(lastStateWithHit);

    TVector3 momGlobal = getGlobalMomBetweenStates(firstState, lastState);
    for(size_t i=0 ; i < track.getStates().size() ; i++){
        track.getStates().at(i).setLocalDirGlobalDir(momGlobal);
    }
    for(size_t i=0 ; i < (track.getStates().size()-1) ; i++){
        float intersectionPoint[3];
        TVector3 momentumAtIntersection;
        float arcLength;
        int holder; //This is used to return the plane which is found.
        bool found =track.getStates().at(i).findIntersectionWithCertainID(track.getStates().at(i+1).getLocation(), intersectionPoint, momentumAtIntersection,arcLength,holder );
        if(!found){
            return false; //DO NOT NEED TO RETURN TRUE. TRUE BY DEFAULT
        }
        track.getStates().at(i+1).setPositionGlobal(intersectionPoint);
    }
    trackOut = EUTelTrack(track);

}
TVector3 EUTelPatternRecognition::getGlobalMomBetweenStates(EUTelState firstState, EUTelState lastState){
    const double * firstPos = firstState.getHit().getPosition();
    double firstPosGlobal[3];
	geo::gGeometry().local2Master(firstState.getLocation() ,firstPos,firstPosGlobal);
    const double * lastPos = lastState.getHit().getPosition();
    double lastPosGlobal[3];
	geo::gGeometry().local2Master(lastState.getLocation() ,lastPos,lastPosGlobal);

    float incidenceXZ = (lastPosGlobal[0] - firstPosGlobal[0])/(lastPosGlobal[2] - firstPosGlobal[2]);
    float incidenceYZ = (lastPosGlobal[1] - firstPosGlobal[1])/(lastPosGlobal[2] - firstPosGlobal[2]);

    TVector3 momGlobal;
    momGlobal[0] = (firstState.getDirLocal().Mag())*incidenceXZ;  
    momGlobal[1] = (firstState.getDirLocal().Mag())*incidenceYZ;  
    momGlobal[2] = sqrt(pow(firstState.getDirLocal().Mag(),2) - pow(firstState.getDirLocalX(),2) - pow(firstState.getDirLocalY(),2));
    return momGlobal;
}


//Note loop through all planes. Even the excluded. This is easier since you don't have to change this input each time then.
//TO DO: sensors z position as used here only works if there is sufficient difference between planes. In the order of 1mm in gear file. This is too large.
void EUTelPatternRecognition::setPlaneDimensionsVec(EVENT::IntVec planeDimensions){
	if(planeDimensions.size() != geo::gGeometry().sensorZOrdertoIDs().size()){
		streamlog_out(ERROR) << "The size of planesDimensions input is: "<< planeDimensions.size()<<" The size of sensorZOrdertoIDs is: " << geo::gGeometry().sensorZOrdertoIDs().size()<< std::endl;
		throw(lcio::Exception( "The input dimension vector not the same as the number of planes!"));
	}
	_planeDimensions.clear();
	for(size_t i=0; i<geo::gGeometry().sensorZOrdertoIDs().size(); ++i){//Loop through each plane
		int planeID = geo::gGeometry().sensorZOrdertoIDs().at(i);
		if (_planeDimensions.find(planeID) == _planeDimensions.end()){//This is to check that we don't try to map the same sensor to two different plane dimensions.
			_planeDimensions[planeID] = planeDimensions.at(i);
		}else{
			streamlog_out(ERROR5) <<"The z position is : "<< i <<" The plane ID is: " << planeID <<std::endl;
			throw(lcio::Exception( "You are trying to map the same sensor ID to two different plane dimensions. There is something wrong with you gear file input. Make sure there is some distance between your planes in the gear file!"));
		}
	}//END of loop over planes
}	    


void EUTelPatternRecognition::testHitsVecPerPlane(){
	if(_mapHitsVecPerPlane.size() !=  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()){
		streamlog_out(ERROR0) <<"The size of the planes with hits " << _mapHitsVecPerPlane.size() <<"  Sensors from Geometry with no excluded planes: "<<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()<<std::endl;
		throw(lcio::Exception("The number of planes that could contain hits and the number of planes is different. Problem could be the gear file has to planes at the same z position.")); 	

	}
	for(size_t i=0 ;i<_mapHitsVecPerPlane.size();++i){
		int planeID = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i);
		if(_mapHitsVecPerPlane.at(planeID).size() <= 0){
			streamlog_out(DEBUG0) << "One plane has no hits at all. Is this correct?" << std::endl;
		}
	}
}

//This function makes sure that the input for the dimension size is correct.
//NOTE:This input is alway for ALL planes and not just no excluded ones. This means we don't have to change this if we want to exclude planes. 
void EUTelPatternRecognition::testPlaneDimensions(){
	if(_planeDimensions.size() != geo::gGeometry().sensorZOrdertoIDs().size()){
		streamlog_out(ERROR5) << "The size of _planesDimensions is: "<< _planeDimensions.size()<<" The size of sensorZOrdertoIDs is: " << geo::gGeometry().sensorZOrdertoIDs().size()<< std::endl;
		throw(lcio::Exception( "The output dimension vector is not the same size as the number of planes. Could be something to do with the gear file. Make sure you have some distances between you planes!"));
	}
	for(size_t i=0;i<_planeDimensions.size();++i){
		if(_planeDimensions.at(geo::gGeometry().sensorZOrdertoIDs().at(i))>2 or _planeDimensions.at(geo::gGeometry().sensorZOrdertoIDs().at(i))<=0){
			throw(lcio::Exception( "The number of dimension for one of your planes is greater than 2 or less than 0. If this is not a mistake collect you nobel prize now!"));
		}
	}
}

    /** Find the hit closest to the intersection of a track with given sensor
     * 
     * @param ts track state
     * @return hit closest to the intersection of the track with the sensor plane
     * 
	 */
const EVENT::TrackerHit* EUTelPatternRecognition::findClosestHit(EUTelState& state)
{
	EVENT::TrackerHitVec& hitInPlane = _mapHitsVecPerPlane[state.getLocation()];
	double maxDistance = std::numeric_limits<double>::max();
	EVENT::TrackerHitVec::const_iterator itClosestHit;
	EVENT::TrackerHitVec::const_iterator itHit;

	streamlog_out(DEBUG0) << "N hits in plane " << state.getLocation() << ": " << hitInPlane.size() << std::endl;

	for ( itHit = hitInPlane.begin(); itHit != hitInPlane.end(); ++itHit ) {
		double distance;
		if(state.getDimensionSize() == 2){
			distance = sqrt(computeResidual( state, *itHit ).Norm2Sqr());//This distance could be 2D or 1D depending on if you have a strip or pixel sensor. Norm2Sqr does not square toot for some reason.
		}else if(state.getDimensionSize() == 1){
			distance = computeResidual( state, *itHit )[0];//If strip sensor then use only displacement along strips. Which should be x axis.
		}else{
			throw(lcio::Exception( "When finding the closest hit to predicted state we find a hit which is not a strip or pixel sensor. Since the dimensionality if less than 1 or greater than 2."));
		}
		streamlog_out(DEBUG0) << "Distance^2 between hit and track intersection: " << distance << std::endl;
		if ( distance < maxDistance ) {
			itClosestHit = itHit;
			maxDistance = distance;
		}
	}
	streamlog_out(DEBUG0) << "Minimal distance^2 between hit and track intersection: " << maxDistance << std::endl;

	return *itClosestHit;
}

//TODO: Need proper error analysis to calculate this rather than providing an answer.  
//This is not very useful at the moment since the covariant matrix for the hit is guess work at the moment.   
double EUTelPatternRecognition::getXYPredictionPrecision(EUTelState& /*ts*/) const 
{
	//sqrt( Ckkm1[0][0]*Ckkm1[0][0] + Ckkm1[1][1]*Ckkm1[1][1] );
	return getWindowSize();   
}

        
/** Calculate residual vector between given track and hit
* 
* @param ts track state
* @param hit hit
* @return 
*/
TVectorD EUTelPatternRecognition::computeResidual(EUTelState& state, const EVENT::TrackerHit* hit) const
{
	const double* hitPosition = hit->getPosition();//In local coordinates

	double localPosition [3];
	localPosition[0] = state.getPosition()[0];
	localPosition[1] = state.getPosition()[1];
	localPosition[2] = state.getPosition()[2];

	streamlog_out(DEBUG3)	<< "Hit (id=" << hit->id() << ") local(u,v) coordinates of hit: (" 
				<< hitPosition[0] << ", " << hitPosition[1] << ", " << hitPosition[2] << ")" << std::endl
				<< "Prediction for hit (id=" << hit->id() << ") local(u,v) coordinates of state: ("  
				<< localPosition[0] << ", " << localPosition[1] << ", " <<localPosition[2] << ")" << std::endl;

	TVectorD residual(2);
	residual[0] = 0; 
	residual[1] = 0;
	
	streamlog_out(DEBUG2) << "The size of this state dimension: "<< _planeDimensions.at(state.getLocation())<< std::endl;
	
	//This loop is used in case we have a strip sensor. So we should use only 1 dimension of information
	for(int i=0; i<_planeDimensions.at(state.getLocation()); ++i)
	{
		residual[i] = hitPosition[i] - localPosition[i];
	}	
	
	if (streamlog_level(DEBUG2))
	{
			streamlog_out(DEBUG2) << "Residual vector residual: (" << residual[0] << ", " << residual[1] << ")" << std::endl;
			residual.Print();
	}
	return residual;
}

void EUTelPatternRecognition::printHits()
{
	streamlog_out(MESSAGE0) << "EUTelPatternRecognition::prinitHit: BEGIN ==============" << std::endl;
	EVENT::TrackerHitVec::const_iterator itHit;
	for(itHit = _allHitsVec.begin(); itHit != _allHitsVec.end(); ++itHit )
	{
		const double* uvpos = (*itHit)->getPosition();
		const int sensorID = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHit) );
		streamlog_out(MESSAGE0)	<< "Hit (id=" << std::setw(3) << sensorID << ") local(u,v) coordinates: (" 
					<< std::setw(7) << std::setprecision(4) << uvpos[0] << ", " << std::setw(7) 
					<< std::setprecision(4) << uvpos[1] << ")" << std::endl;
	}
	streamlog_out(MESSAGE0) << "EUTelPatternRecognition::printHits: END ==============" << std::endl;
}

} // namespace eutelescope
