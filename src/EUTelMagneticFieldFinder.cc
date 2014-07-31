// Created on July 2, 2013, 12:59 PM
// This code should take always a trackerHit object. Since it will never need to use a derived hit object with other functions since it only needs basic information about the hit.
#include "EUTelMagneticFieldFinder.h"


namespace eutelescope {

    namespace {

        struct hasSensorID
        {
            int _requiredSensorID;

            bool operator( )( MeasurementLayer* layer ) {
                return ( layer->sensorID( ) == _requiredSensorID );
            }
        } ;
    }
    
    EUTelKalmanFilter::EUTelKalmanFilter() :  
            _tracksCartesian(), 
            _trackStates(), 
            _allHits(),
            _allMeasurements(),
            _isHitsOK(false),
            _allowedMissingHits(0),
						_AllowedSharedHitsOnTrackCandidate(0),
            _maxTrackCandidates(0),
            _beamE(-1.),
            _beamQ(-1.),
            _beamEnergyUncertainty(0.),
            _beamAngularSpread(2,-1.),
	    _jacobianF(5,5),
            _trkParamCovCkkm1(5,5),
            _processNoiseQ(5,5),
            _residualCovR(2,2),
            _eomIntegrator( new EUTelUtilityRungeKutta() ),
            _jacobianIntegrator( new EUTelUtilityRungeKutta() ),
            _eomODE( 0 ),
            _jacobianODE( 0 ),
            _planesForPR(1) 
            {
                // Initialise ODE integrators for eom and jacobian       
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }

                {
                    _jacobianIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
            }

    EUTelKalmanFilter::EUTelKalmanFilter( std::string name ) : 
            _tracksCartesian(),
            _trackStates(), 
            _allHits(),
            _allMeasurements(),
            _isHitsOK(false),
            _allowedMissingHits(0),
            _maxTrackCandidates(0),
            _beamE(-1.),
            _beamQ(-1.),
            _beamEnergyUncertainty(0.),
            _beamAngularSpread(2,-1.),
	    _jacobianF(5,5),
            _trkParamCovCkkm1(5,5),
            _processNoiseQ(5,5),
            _residualCovR(2,2),
            _eomIntegrator( new EUTelUtilityRungeKutta() ),
            _jacobianIntegrator( new EUTelUtilityRungeKutta() ),
            _eomODE( 0 ),
            _jacobianODE( 0 ),
            _planesForPR( 1 ) 
            {
                // Initialise ODE integrators for eom and jacobian       
                {
                    _eomODE = new eom::EOMODE(5);
                    // _eomIntegrator integrator becomes the owner of _eomODE and ButcherTableauDormandPrince
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }

                {
                    _jacobianIntegrator->setRhs( _jacobianODE );
                    _jacobianIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
            }

        EUTelKalmanFilter::~EUTelKalmanFilter() { 
        delete _eomIntegrator;
        delete _jacobianIntegrator;
    }
 
    /** */
    int EUTelKalmanFilter::findNextPlaneEntrance(  EUTelTrackStateImpl* ts , int sensorID  ){

// how to make sure that ts is in localframe system??

      const double x0 = ts->getReferencePoint()[0];
      const double y0 = ts->getReferencePoint()[1];
      const double z0 = ts->getReferencePoint()[2];
      int tsPlaneID   = ts->getLocation();
	

      // get the very first sensor : why?
      // remember last hit-point from the track candidate below
      double lpoint[3]   = { x0, y0, z0}; // 
      double gpoint[3]   = {0.,0.,0.};  // initialise output point-vector global frame (World)

      geo::gGeometry().local2Master( tsPlaneID, lpoint, gpoint );
      streamlog_out ( DEBUG2 ) << "findNextPlaneEntrance: " << tsPlaneID << " " <<  lpoint[0] << " " <<  lpoint[1] << " " << lpoint[2] <<  " :"
                                << gpoint[0] << " "    << gpoint[1] << " "   << gpoint[2] << " " << endl; 

      gpoint[2] += -100.; 
//      float start[3]     = {static_cast<float> (gpoint[0]),static_cast<float> (gpoint[1]),static_cast<float> (gpoint[2]) };  

      double dir[3]      = {0.,0.,1.};  // as good as any other value along z axis.

      float nextPoint[3] = {0.,0.,0.};  // initialise output point-vector local frame (sensor)

            int newSensorID   = geo::gGeometry( ).findNextPlaneEntrance( gpoint,  dir, sensorID, nextPoint ) ;

// get to local frame for sensor newSensorID:
            double dGpoint[3] = {static_cast<double> (nextPoint[0]),static_cast<double> (nextPoint[1]),static_cast<double> (nextPoint[2]) };
            double dLpoint[3] = {0.,0.,0.};  // initialise output point-vector local frame (sensor)
            geo::gGeometry().master2Local( dGpoint, dLpoint);

            if( newSensorID < 0 ) 
            {
              streamlog_out ( DEBUG4 ) << "no Entrance: " <<  lpoint[0] << " " <<  lpoint[1] << " " << lpoint[2] << " err:"<< newSensorID << endl;
            } else {     
              streamlog_out ( DEBUG2 ) << "identified NextPlane Entrance: " <<  lpoint[0] << " " <<  lpoint[1] << " " << lpoint[2] <<  " at : "<< newSensorID << endl;
              const float opoint[3] = { lpoint[0], lpoint[1], lpoint[2] };
              ts->setReferencePoint( opoint );  
              ts->setLocation( newSensorID );
           }

      return newSensorID;
    }

    /** */


std::vector<EUTelTrack>& EUTelKalmanFilter::getTracks(){
	return	_trackAfterSameHitsCut; 
}
void EUTelKalmanFilter::clearEveryRun(){
	int _totalNumberOfHits=0;
	int _totalNumberOfSharedHits=0;
}

//This is the work horse of the class. Using seeds it propagates the track forward using equations of motion. This can be with or without magnetic field.
void EUTelKalmanFilter::propagateForwardFromSeedState( EUTelTrack& state, EUTelTrack & track    ){
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::propagateForwardFromSeedState-----BEGIN "<< endl;
	//TO DO: To delete this in smart way. 
	EUTelTrack *firstState = new EUTelTrack(state);//Need to create an initial state that will not be deleted outside this scope  
	track.addTrack(static_cast<EVENT::Track*>(firstState));//Note we do not have to create new since this object State is saved in class member scope
	//Here we loop through all the planes not excluded. We begin at the seed which might not be the first. Then we stop before the last plane, since we do not want to propagate anymore
	for(int i = geo::gGeometry().sensorIDToZOrderWithoutExcludedPlanes()[state.getLocation()]; i < (geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()-1); ++i){
	
		float globalIntersection[3];
		int newSensorID = state.findIntersectionWithCertainID(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i+1], globalIntersection);
		int sensorIntersection = geo::gGeometry( ).getSensorID(globalIntersection);//This is needed with the jacobian since the jacobian needs the z distance to the next point to do the calculation
		if(newSensorID < 0 ){ //TO DO Fix geo: or sensorIntersection < 0 ){
			streamlog_out ( MESSAGE5 ) << "Intersection point on infinite plane: " <<  globalIntersection[0]<<" , "<<globalIntersection[1] <<" , "<<globalIntersection[2]<<std::endl;
			streamlog_out(MESSAGE5) <<" From ID= " <<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i]<< " to " <<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i+1]  <<std::endl;
			streamlog_out(MESSAGE5)<<"Was there intersection on plane: "<<newSensorID<<" Was there intersection in sensitive area: "<< sensorIntersection <<" Both number -999 if no intersection" << std::endl;
			streamlog_out(MESSAGE5) << Utility::outputColourString("No intersection found moving on plane. Move to next plane and look again.","YELLOW")<<std::endl; 
			continue;//So if there is no intersection look on the next plane. Important since two planes could be at the same z position
		}
		//So we have intersection lets create a new state
		EUTelTrack *newState = new EUTelTrack();//Need to create this since we save the pointer and we would be out of scope when we leave this function. Destroying this object. 
		newState->setDimensionSize(_planeDimensions[newSensorID]);//We set this since we need this information for later processors
		newState->setBeamEnergy(_beamE);
		newState->setBeamCharge(_beamQ);
		newState->initialiseCurvature();
		newState->setLocation(newSensorID);
		TMatrixD jacobian(5,5);
		jacobian.Zero();
		jacobian = newState->computePropagationJacobianFromStateToThisZLocation(globalIntersection[2]);
		findNextStateUsingJacobian(state,jacobian,globalIntersection[2], *newState);
		testPositionEstimation(const_cast<float*>(newState->getReferencePoint()),globalIntersection);//Here we compare the estimation from the simple equation of motion and jacobian.
		if(_mapHitsVecPerPlane[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i+1]].size() == 0){
			streamlog_out(DEBUG5) << Utility::outputColourString("There are no hits on the plane with this state. Add state to track as it is and move on ","YELLOW");
			track.addTrack(static_cast<EVENT::Track*>(newState));//Need to return this to LCIO object. Loss functionality but retain information 
			state = *newState;
			continue;
		}
		EVENT::TrackerHit* closestHit = const_cast< EVENT::TrackerHit* > ( findClosestHit( *newState )); //This will look for the closest hit but not if it is within the excepted range		
		const double* hitPosition = closestHit->getPosition();
		const double distance = std::sqrt(computeResidual(*newState, closestHit ).Norm2Sqr()); //Determine the residual to it. //Distance is in mm.
		const double DCA = getXYPredictionPrecision( *newState ); //This does nothing but return a number specified by user. In the future this should use convariance matrix information TO DO: FIX
		streamlog_out ( DEBUG1 ) <<"At plane: "<<newState->getLocation() << ". Distance between state and hit: "<< distance <<" Must be less than: "<<DCA<< endl;
		if ( distance > DCA ) {
			streamlog_out ( DEBUG1 ) << "Closest hit is outside of search window." << std::endl;
			streamlog_out ( DEBUG1 ) << "Skipping current plane. Covariant Matrix and position already updated to this point. Add to track " << std::endl;
			track.addTrack(static_cast<EVENT::Track*>(newState));//Need to return this to LCIO object. Loss functionality but retain information 
			state = *newState;
			continue;
		}	
		if(closestHit == NULL){
			throw(lcio::Exception( Utility::outputColourString("The closest hit you are trying to add is NULL. This can not be correct","RED")));
		}
		streamlog_out ( DEBUG1 ) << "Found a hit! Distance: "<< distance << endl;
		newState->addHit(closestHit);
		_totalNumberOfHits++;//This is used for test of the processor later.   
		streamlog_out ( DEBUG1 ) << "Add state to track. "<< endl;
		track.addTrack(static_cast<EVENT::Track*>(newState));//Need to return this to LCIO object. Loss functionality but retain information 
		streamlog_out ( DEBUG1 ) << "Make new state just state. "<< endl;
		state = *newState;
		streamlog_out ( DEBUG1 ) << "End of loop "<< endl;

	}//TO DO: Must add the kalman filter back into the code. 
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::propagateForwardFromSeedState-----END "<< endl;

}	
		/*	
	           	streamlog_out ( DEBUG1 ) << "NextPlane " << newSensorID << " " << uvpos[0] << " " << uvpos[1] << " " << uvpos[2] << " resid:" << distance << " ResidualCut: " << DCA << endl;
           				else{
            			streamlog_out (DEBUG1 ) << "NextPlane MATCHED. Position of Hit (Local) " <<  uvpos[0] << " " <<  uvpos[1] << " " << uvpos[2] 
                                                        <<" Position of state (Global) " << state_new->getX() 
                                                        << "," << state_new->getY()
                                                        <<","<< state_new->getZParameter() 
                                                        <<" Distance between them: " << distance 
                                                        << " Sensor ID:" << state_new->getLocation() 
                                                        << " Seed we began at: " << (*itTrk) << endl;
	              		streamlog_out ( DEBUG1 ) << "Will now alter Cov matrix and state variables using hit information " << std::endl;
				TMatrixD HMatrix = state_new->getH(); //We need to be able to move from the measurement to the state space
				TMatrixD GainMatrix = updateGainK( state_new, closestHit ); //This is a matrix that tells you how much the state should be changed with the information from the hit
				UpdateStateUsingHitInformation( state_new ,closestHit, jacobian, GainMatrix, HMatrix); //Update the state on the track
				streamlog_out (DEBUG0) << "Here is the state variable after update with hit and propagation (state_new): " <<std::endl;
				state_new->Print();
				UpdateTrackUsingHitInformation( state_new , closestHit, (*itTrk), jacobian, GainMatrix,HMatrix); //Update the track itself.									
            		}

		}
                if( streamlog_level(DEBUG1) ) state_new->Print();
  		(*itTrk)->addTrackState( new EUTelTrackStateImpl(*state_new) ); //New memory allocation here should be deleted by LCIO memory management. I think...
		streamlog_out ( DEBUG1 ) << "Memory location of initial state after all allocation (Should be the same): " << state << endl; 
		streamlog_out ( DEBUG1 ) << "Memory location of state_new: " << state_new << endl; 
						
	 	state = state_new; //MUST DOUBLE CHECK THIS MEMORY MANAGEMENT DOES NOT LOOK RIGHT TO ME
		streamlog_out ( DEBUG1 ) << "Memory location of initial state after it was made equal to state_new: " << state << endl; 		

          }
	  ++iter;	
       
 
    }
*/
void EUTelKalmanFilter::testTrackQuality(){
	_numberOfTracksTotal = _numberOfTracksTotal + _tracks.size();
	_numberOfTracksAfterHitCut = _numberOfTracksAfterHitCut + _tracksAfterEnoughHitsCut.size();
	_numberOfTracksAfterPruneCut = _numberOfTracksAfterPruneCut +	_trackAfterSameHitsCut.size();

	if(_numberOfTracksTotal % 1000 == 0){
		streamlog_out(MESSAGE5)<<"///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
		streamlog_out(MESSAGE5)<<" Total Tracks "<<" Pass Hit Cut " <<"                           Pass Prune Cut " <<std::endl;
		streamlog_out(MESSAGE5)<<"    "<<_numberOfTracksTotal<<"    "<<"    "<<_numberOfTracksAfterHitCut<<"                                                     "<<_numberOfTracksAfterPruneCut<<std::endl;
		float percentAfterHitCut = (_numberOfTracksAfterHitCut/_numberOfTracksTotal)*100;
		float percentAfterPruneCut = (_numberOfTracksAfterPruneCut/_numberOfTracksTotal)*100;
		streamlog_out(MESSAGE5)<<"               "<<percentAfterHitCut<<"                                                      "<<percentAfterPruneCut<<std::endl;
		float averageNumberOfHitsOnTrack = _totalNumberOfHits/_numberOfTracksTotal;
		float averageNumberOfSharedHitsOnTrack = _totalNumberOfSharedHits/ _numberOfTracksTotal;
		streamlog_out(MESSAGE5)<<"The average number of hits on a track: "<< averageNumberOfHitsOnTrack<<std::endl;
		streamlog_out(MESSAGE5)<<"The average number of shared hits on a track with other tracks: "<< averageNumberOfSharedHitsOnTrack<<std::endl;
		if(percentAfterHitCut < 0.1){
			streamlog_out(MESSAGE5)<< Utility::outputColourString("The percentage of track making the hit cut is very low at the moment ","YELLOW")<<std::endl;
		}
		if(percentAfterPruneCut < 0.1){
			streamlog_out(MESSAGE5)<< Utility::outputColourString("The percentage of track making the prune cut is very low at the moment ","YELLOW")<<std::endl;
		}	
	}
}
void EUTelKalmanFilter::testPositionEstimation(float position1[], float position2[]){
	if(position1[0] -position2[0] > 0.004){ //in mm
		throw(lcio::Exception( Utility::outputColourString("The positions between the two methods is different in the x direction. ","RED"))); 
	}
if(position1[1] -position2[1] > 0.004){ //in mm
		throw(lcio::Exception( Utility::outputColourString("The positions between the two methods is different in the y direction. ","RED"))); 
	}
 if(position1[2] -position2[2] > 0.004){ //in mm
		throw(lcio::Exception( Utility::outputColourString("The positions between the two methods is different in the z direction. ","RED"))); 
	} 
}

    void EUTelKalmanFilter::Print( std::string Name, std::vector< EUTelTrackImpl*> & _collection) 
    {
       int itrk = 0; 
       int size_itTrk  = _collection.size();

       std::vector< EUTelTrackImpl* >::iterator itTrk;
 
       for ( itTrk = _collection.begin(); itTrk != _collection.end(); itTrk++, itrk++ ) 
       {
          if( (*itTrk) == 0 ) 
          {         
             streamlog_out(WARNING1) << " Track vector: " << Name.c_str() << " contains no information at track = " << itrk << "." << std::endl;
             continue; 
          }
          const EVENT::TrackerHitVec& ihits = (*itTrk)->getTrackerHits();
          int nhits =  ihits.size( ) ;
          int expec =  geo::gGeometry( ).nPlanes( ) - _allowedMissingHits;
          streamlog_out(DEBUG5) << " Track vector " << Name.c_str() << " at track number " <<  itrk << " of size " << size_itTrk << " with " << nhits << " and expecting at least " << expec << " hits " << std::endl;
          (*itTrk)->Print();
       }
    }
 
//    void EUTelKalmanFilter::Prune( std::vector< EUTelTrackImpl*> &_collection, std::vector< EUTelTrackImpl*> &_collection_to_delete ) 
//    {
    // Perform track pruning this removes tracks that have the same hits used to create the track on some planes
void EUTelKalmanFilter::findTrackCandidatesWithSameHitsAndRemove() {

	streamlog_out(MESSAGE1) << "EUTelKalmanFilter::findTrackCandidatesWithSameHitsAndRemove" << std::endl;
	for(int i =0; i < _tracksAfterEnoughHitsCut.size();++i){ 
		streamlog_out(DEBUG1) <<  "Loop at track number:" <<  i << std::endl;
		const EVENT::TrackerHitVec& ihits = _tracksAfterEnoughHitsCut[i].getTrackerHits(); //get the hits contained within this track object
 
		//Now loop through all tracks one ahead of the original track itTrk. This is done since we want to compare all the track to each other to if they have similar hits     
		for(int j =0; j < _tracksAfterEnoughHitsCut.size();++j){ 
			int hitscount=0;
			const EVENT::TrackerHitVec& jhits = _tracksAfterEnoughHitsCut[j].getTrackerHits();
			for(int i=0;i<ihits.size();i++){
				EVENT::TrackerHit* ihit = ihits[i];
				int ic = ihit->id();
				for(int j=0;j<jhits.size();j++){
					EVENT::TrackerHit* jhit = jhits[j];
					int jc = jhit->id();
					if(ic == jc ){
						_totalNumberOfSharedHits++;
						hitscount++; 
						streamlog_out(MESSAGE1) <<  "Hit number on track you are comparing all other to :" << i << ". Hit ID: " << ic << ". Hit number of comparison: " << j << ". Hit ID of this comparison : " << jc << ". Number of common hits: " << hitscount << std::endl; 
					}
				}
			} 
			if(hitscount > _AllowedSharedHitsOnTrackCandidate) {   
				_trackAfterSameHitsCut.push_back(_tracksAfterEnoughHitsCut[i]);
			}
			continue;
		}
	}

	streamlog_out(MESSAGE1) << "------------------------------EUTelKalmanFilter::findTrackCandidatesWithSameHitsAndRemove()---------------------------------" << std::endl;
}

//This function can be expanded to check other user input that may be needed. 
//Furthermore you can put other functions that cause _userInputGood to fail.        
void EUTelKalmanFilter::testUserInput() {
	streamlog_out(DEBUG2) << "EUTelKalmanFilter::testUserInput()" << std::endl;

	// Check the validity of supplied beam energy 
	if ( _beamE < 1.E-6 ) {
		throw(lcio::Exception( Utility::outputColourString("Beam direction was set incorrectly","RED"))); 
	}
	else{
	 streamlog_out(DEBUG0) << Utility::outputColourString("Beam energy is reasonable", "GREEN") << std::endl;
	}

	if ( _beamEnergyUncertainty < 0 ) {
		throw(lcio::Exception( Utility::outputColourString("Beam uncertainty is negative. Check supplied values","RED"))); 
	}else{
	 streamlog_out(DEBUG0) << Utility::outputColourString("Beam Energy uncertainty is reasonable","GREEN") << std::endl;
	}

	if(_beamAngularSpread.size() == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam spread size is zero.","RED"))); 
	}
	if(_beamAngularSpread[0] <= 0  or  _beamAngularSpread[1] <= 0){
		throw(lcio::Exception( Utility::outputColourString("The beam spread is zero.","RED"))); 
	}else{
		streamlog_out(DEBUG0) << Utility::outputColourString("The size of beam spread is" +  to_string(_beamAngularSpread.size()),"GREEN") << std::endl;
	}
	if(_createSeedsFromPlanes.size() == 0){
		throw(lcio::Exception( Utility::outputColourString("The number of planes to make seeds from is 0. We need at least one plane", "RED")));
	}

	if(_createSeedsFromPlanes.size() >= geo::gGeometry().sensorIDstoZOrder().size()){
		throw(lcio::Exception( Utility::outputColourString("You have specified all planes or more than that. This is too many planes for seeds", "RED")));
	}else{
		streamlog_out(DEBUG0) << Utility::outputColourString("The number of planes to make seeds from is good " +  to_string(_createSeedsFromPlanes.size()),"GREEN") << std::endl;
	}
	if(_excludePlanes.size() >= geo::gGeometry().sensorIDstoZOrder().size()){//TO DO: should check if seed planes are also excluded	
	 streamlog_out(DEBUG0) << Utility::outputColourString("Beam energy is reasonable", "GREEN") << std::endl;
	}

	if ( _beamEnergyUncertainty < 0 ) {
		throw(lcio::Exception( Utility::outputColourString("Beam uncertainty is negative. Check supplied values","RED"))); 
	}else{
	 streamlog_out(DEBUG0) << Utility::outputColourString("Beam Energy uncertainty is reasonable","GREEN") << std::endl;
	}

	if(_beamAngularSpread.size() == 0){
		throw(lcio::Exception( Utility::outputColourString("The beam spread size is zero.","RED"))); 
	}
	if(_beamAngularSpread[0] <= 0  or  _beamAngularSpread[1] <= 0){
		throw(lcio::Exception( Utility::outputColourString("The beam spread is zero.","RED"))); 
	}else{
		streamlog_out(DEBUG0) << Utility::outputColourString("The size of beam spread is" +  to_string(_beamAngularSpread.size()),"GREEN") << std::endl;
	}
	if(_createSeedsFromPlanes.size() == 0){
		throw(lcio::Exception( Utility::outputColourString("The number of planes to make seeds from is 0. We need at least one plane", "RED")));
	}

	if(_createSeedsFromPlanes.size() >= geo::gGeometry().sensorIDstoZOrder().size()){
		throw(lcio::Exception( Utility::outputColourString("You have specified all planes or more than that. This is too many planes for seeds", "RED")));
	}else{
		streamlog_out(DEBUG0) << Utility::outputColourString("The number of planes to make seeds from is good " +  to_string(_createSeedsFromPlanes.size()),"GREEN") << std::endl;
	}
	if(_excludePlanes.size() >= geo::gGeometry().sensorIDstoZOrder().size()){//TO DO: should check if seed planes are also excluded	
		throw(lcio::Exception( Utility::outputColourString("The number of excluded planes is too large. We can not fit tracks with magic.", "RED")));
	}
	if(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() != (geo::gGeometry().sensorIDstoZOrder().size()-_excludePlanes.size())){
		throw(lcio::Exception( Utility::outputColourString("The number of Planes-Excluded is not correct. This could be a problem with geometry", "RED")));
	}else{
		streamlog_out(DEBUG0) << Utility::outputColourString("The number of excluded planes is" +  to_string(_excludePlanes.size()),"GREEN") << std::endl;
	}
}	
void  EUTelKalmanFilter::testHitsVec(){
	if(_allHitsVec.size() == 0){
		throw(lcio::Exception( "The hit input is wrong!"));
	}

}

		
		/**
		 *
		 * Prune seed track states necessary to
		 *
     * start Kalman filter
     * 
     * 
     */
    void EUTelKalmanFilter::pruneSeeds() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::pruneSeeds()" << std::endl;

        if ( _allMeasurements.empty() ) {
            streamlog_out(WARNING1) << "Can't initialise track seeds for the finder. No hits in this event." << std::endl;
            return;
        }


	// Start Kalman filter
	std::vector< EUTelTrackImpl* >::iterator itTrk;
	for ( itTrk = _tracksCartesian.begin(); itTrk != _tracksCartesian.end(); itTrk++ ) {
            bool isDuplTrack = false; 
             
            std::vector< EUTelTrackImpl* >::iterator itTrk_in;
	    for ( itTrk_in = _tracksCartesian.begin(); itTrk_in != _tracksCartesian.end(); itTrk_in++ ) {
              bool isDuplTrack = false; 

            }
            
        }
 
    }
   
/** 
*   Using a list of planes use the hits on those planes to create the initial seeds. 
*   These seeds are used then to propagate to the next plane.
* 
*/
//TO DO:This will not work with a strip sensor creating a seed. Since the strip sensor creates a line of possible hit positions not just a single point. Therefore you need to loop through  all points on the line and determine its trajectory and then see if the intersection is close enough to the hit on the next plane to add that hit on the track. This is a lot of work for something which is not that important.   
void EUTelKalmanFilter::initialiseSeeds() {
	streamlog_out(DEBUG2) << "EUTelKalmanFilter::initialiseSeeds()" << std::endl;
	_mapSensorIDToSeedStatesVec.clear();
	for( int iplane = 0; iplane < _createSeedsFromPlanes.size() ; iplane++) {

		streamlog_out(DEBUG1) << "We are using plane: " <<  _createSeedsFromPlanes[iplane] << " to create seeds" << std::endl;

		EVENT::TrackerHitVec& hitFirstLayer = _mapHitsVecPerPlane[_createSeedsFromPlanes[iplane]];
		streamlog_out(DEBUG1) << "N hits on sensor : " << hitFirstLayer.size() << std::endl;
		if(hitFirstLayer.size()== 0){
			continue;
		}
		std::vector<EUTelTrack> stateVec;
		EVENT::TrackerHitVec::iterator itHit;
		for ( itHit = hitFirstLayer.begin(); itHit != hitFirstLayer.end(); ++itHit ) {
			EUTelTrack state;//Here we create a track state. This is a point on a track that we can give a position,momentum and direction. We combine these to create a track. 
			double posLocal[] =  { (*itHit)->getPosition()[0], (*itHit)->getPosition()[1], (*itHit)->getPosition()[2] };
			if(posLocal[2] != 0){//TO DO:This should be in test hits
				streamlog_out(MESSAGE5) << "The local position of the hit is " << posLocal[0]<<","<< posLocal[1]<<","<< posLocal[2]<<","<<std::endl;
				throw(lcio::Exception(Utility::outputColourString("The position of this local coordinate in the z direction is non 0", "RED"))); 	
			}
			double temp[] = {0.,0.,0.};
			geo::gGeometry().local2Master(_createSeedsFromPlanes[iplane] , posLocal, temp);//IMPORTANT:For strip sensors this will make the hit strip look like a pixel at (Xstriplocal,somevalue,somevalue).
			//There is no way round this since the strip give a range of possible positions in 2D/3D space.      //To allow strip sensor to be seed would involve major change to how this is done. For now just don't do it  
			float posGlobal[] = { static_cast<float>(temp[0]), static_cast<float>(temp[1]), static_cast<float>(temp[2]) };
			streamlog_out ( DEBUG2 ) << "pick next hit " << _createSeedsFromPlanes[iplane]<< " " <<  posLocal[0] << " " << posLocal[1] << " " << posLocal[2] <<  " : "<< posGlobal[0] << " "    << posGlobal[1] << " "   << posGlobal[2] << " " << endl;
			state.setLocation(_createSeedsFromPlanes[iplane]);//This stores the location as a float in Z0 since no location for track LCIO. This is preferable to problems with storing hits.  
			state.setPosition(posGlobal); //This will automatically take cartesian coordinate system and save it to reference point. //This is different from most LCIO applications since each track will have own reference point. 		
			state.setDirectionXY(0.);  // seed with 0 at first hit. This one is stored as tan (x/y) rather than angle as in LCIO format TO DO:This really should be give as a beam direction.  
			state.setDirectionYZ(0.);  // seed with 0 at first hit.
			state.setBeamCharge(_beamQ);//this is set for each state. to do: is there a more efficient way of doing this since we only need this stored once?
			state.setBeamEnergy(_beamE);//this is saved in gev. 
			state.initialiseCurvature(); //this will perform the calculation _beamq/_beame ad place in invp
			//TO DO We just fill the cov matrix with junk. We don't use it at the moment 
	/*		float trkCov[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};//This is the covariance matrix of the track state: (1/p, Tx,Ty,x,y)X(1/p, Tx,Ty,x,y) Z is not included since this is a parameter.This is also a symmetric matrix so only need 15 entries. 
			trkCov[0] = ( _beamEnergyUncertainty * _beamE ) * ( _beamEnergyUncertainty * _beamE );               //cov(q/p,x)=0, cov(q/p,y)=0, cov(q/p,tx)=0, cov(q/p,ty)=0, cov(q/p,q/p)

			trkCov[2] = _beamAngularSpread[0] * _beamAngularSpread[0];          //cov(tx,x)=0, cov(tx,y)=0, cov(tx,tx)
			trkCov[5] = _beamAngularSpread[1] * _beamAngularSpread[1];          //cov(ty,x)=0, cov(ty,y)=0, cov(ty,tx)=0, cov(ty,ty)

			//TO DO: Need to get proper convariance matrix for hits. Since we make the hits error so large the Kalman filter will not work as intended. This would also make life easier to deal with strip sensors since we could make the variance very large for the dimension we do not know and use the same code as pixel sensor.  
			const EVENT::FloatVec uvcov = (*itHit)->getCovMatrix();
			trkCov[13] = uvcov[0]*1.E4;                               //cov(x,x)
			trkCov[14] = uvcov[1]*1.E4;   trkCov[15] = uvcov[2]*1.E4;       //cov(y,x), cov(y,y)
			state.setCovMatrix(trkCov);*/
			state.addHit(*itHit);
			state.setDimensionSize(_planeDimensions[state.getLocation()]);
			stateVec.push_back(state);
		}
		_mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[iplane]] = stateVec; 
	//	streamlog_out(MESSAGE5) << "The size of stateVec: "<<stateVec.size()<<"The size of the map " <<_mapSensorIDToSeedStatesVec.size() << std::endl;
	}
	streamlog_out(DEBUG2) << "--------------------------------EUTelKalmanFilter::initialiseSeeds()---------------------------" << std::endl;
}
void EUTelKalmanFilter::testInitialSeeds(){
	//Even if there are no hits then we can compare these. Since it will contain stateVec with size 0.      
//	if(_mapSensorIDToSeedStatesVec.size() != _createSeedsFromPlanes.size()){TO DO:Fix strange error with some event gives no seed planes
	//	streamlog_out(MESSAGE5) <<"The size of sensors with seeds: " << _mapSensorIDToSeedStatesVec.size() <<" The number of planes you are suppose to create seeds from: "<<_createSeedsFromPlanes.size()<< " For event number: " <<getEventNumber()<<std::endl;
//		throw(lcio::Exception(Utility::outputColourString("The size of intial state seeds planes and the number to use at the start are different", "RED"))); 	
//	}
	for(int i = 0 ; i < _mapSensorIDToSeedStatesVec.size(); ++i){
		std::vector<EUTelTrack> StatesVec =  _mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[i]]; 	
		for(int j = 0 ; j < StatesVec.size() ; ++j){
			if(StatesVec[j].getTrackerHits()[0] == NULL ){
				throw(lcio::Exception(Utility::outputColourString("The hit is NULL. All seeds must have hits.", "RED"))); 	
			}
		}
	}
}
//Loop through each plane that contains seeds and then each seed. From that seed you then create a track.
void EUTelKalmanFilter::findTrackCandidates() {
	streamlog_out(MESSAGE1) << "EUTelKalmanFilter::findTrackCandidates()" << std::endl;
	clearTrackAndTrackStates(); //Clear all past track information
	for(int i = 0 ; i < _mapSensorIDToSeedStatesVec.size(); ++i){
		std::vector<EUTelTrack> statesVec =  _mapSensorIDToSeedStatesVec[_createSeedsFromPlanes[i]]; 	
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
	streamlog_out(MESSAGE1) << "EUTelKalmanFilter::findTrackCandidates()------END" << std::endl;
}
//We need to delete the states AND the track. Since we have allocated memory to store the state.
void EUTelKalmanFilter::clearTrackAndTrackStates(){
	for(int i=0; i < _tracks.size(); ++i){
		for(int j=0; j< _tracks.at(i).getTracks().size(); ++j){   
		//	delete [] _tracks.at(i).getTracks().at(j); TO DO: Delete memory properly 
		}
	}
	_tracks.clear();
}

void EUTelKalmanFilter::findTracksWithEnoughHits(){
	streamlog_out(DEBUG1) << "EUTelKalmanFilter::findTracksWithEnoughHits()------BEGIN" << std::endl;
	_tracksAfterEnoughHitsCut.clear();
	if(_tracks.size() == 0 ){
		streamlog_out(MESSAGE5) <<"This is event: " <<getEventNumber()<<std::endl;   
		streamlog_out(MESSAGE5) << Utility::outputColourString("The number of tracks for this event is zero ","YELLOW")<<std::endl; 
	}
	for(int i = 0 ; i<_tracks.size(); ++i){
		EUTelTrack& track = _tracks[i];
		streamlog_out ( DEBUG2 ) << "Number of hits on the track: " <<track.getNumberOfHitsOnTrack()<<" Number needed: " <<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - _allowedMissingHits << std::endl;
		if(track.getNumberOfHitsOnTrack()> geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - _allowedMissingHits){
			streamlog_out(DEBUG5) << "There are enough hits. So attach this track makes the cut!"<<endl;
			_tracksAfterEnoughHitsCut.push_back(track);
		}
	}
	streamlog_out(DEBUG1) << "EUTelKalmanFilter::findTracksWithEnoughHits()-----END" << std::endl;

}
/*
//Check the number of hit on the track after propagation and collecting hits is over the minimum
if ( isGoodTrack && ( *itTrk )->getHitsOnTrack().size( ) < geo::gGeometry( ).nPlanes( ) - _allowedMissingHits ) {
	streamlog_out ( DEBUG5 ) << "Track candidate has to many missing hits." << std::endl;
	streamlog_out ( DEBUG5 ) << "Removing this track candidate from further consideration." << std::endl;
(*itTrk)->Print();
	delete (*itTrk);  //Is this really nee:ded?
	isGoodTrack = false;
	itTrk = _tracksCartesian.erase( itTrk ); itTrk--; 
}

if ( isGoodTrack ) {
//               state->setLocation( EUTelTrackStateImpl::AtLastHit );
	int nstates = (*itTrk)->getTrackStates().size();

	streamlog_out(DEBUG5) << "'Tracks' looped through. I.e initial hit seed as a track. (Ignore states with no hits): "  << local_itTrk << ". Number of seeds: " << size_itTrk << ". At seed: " << *itTrk << " after propagation with: " << ( *itTrk )->getHitsOnTrack().size( ) << " hits collected on track.  Expecting at least " << geo::gGeometry( ).nPlanes( ) - _allowedMissingHits << " The states. I.e Including planes with no hits: "<< nstates << std::endl << std::endl << std::endl << std::endl;

streamlog_out ( DEBUG5 ) << "Successful track state after propagation: " << endl; 
(*itTrk)->Print();
}

local_itTrk++;
}

streamlog_out(MESSAGE0) << "Finished looping through all seeds! Looped at total of : " << local_itTrk << " after seeds with no state are deduced. Total seeds were: " << size_itTrk << endl;

streamlog_out(MESSAGE1) << "------------------------------EUTelKalmanFilter::findTrackCandidates()---------------------------------" << std::endl;
*/
//This creates map between plane ID and hits on that plane. 
//We also order the map correcly with geometry.
void EUTelKalmanFilter::setHitsVecPerPlane(){
	_mapHitsVecPerPlane.clear();
	int numberOfPlanes = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size();//Note should not make this a class data member since we call this by reference in geometry so defacto it is already accessed directly each time. By this I mean we do not create a new copy each time we call the function. 
	if(numberOfPlanes == 0){
		throw(lcio::Exception( "The number of planes is 0 to get hits from."));
	}
	if( _allHitsVec.size()== 0){
		throw(lcio::Exception( "The number of hits is zero."));
	}
	for(int i=0 ; i<numberOfPlanes;++i){
		EVENT::TrackerHitVec tempHitsVecPlaneX; 
		for(int j=0 ; j<_allHitsVec.size();++j){
			if(Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*>(_allHitsVec[j]) ) ==  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i]){
				tempHitsVecPlaneX.push_back(_allHitsVec[j]);
			}
		}		
	_mapHitsVecPerPlane[ geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i]] = 	tempHitsVecPlaneX;
	}	
}
//This member function will loop through each position array and determine if some entries are all one number (Usually 0). This is done so we can determine the dimension of the hit.
//We create a vector that excluded planes that are notoncluded
//TO DO: This way of calculating the dimesions is silly since y ou way have only 1 hit in the first event on any plane. Need to fix this. 
void EUTelKalmanFilter::setPlaneDimensionsVec(){
	_planeDimensions.clear();
	const double* positionBefore=NULL;//This is a pointer to a const not a const pointer
	const double* position=NULL;
	int numberOfPlanes = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size();
	for(int i=0; i<numberOfPlanes; ++i){//Loop through each plane
		int planeID = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i];
		int numberOfDimensions=0;
		bool sensorIsAtLeast1D=false;
		bool sensorIsAtLeast2D=false;
		bool sensorIsAtLeast3D=false;
		for(int j = 0 ; j< _mapHitsVecPerPlane.at(i).size();++j){//Loop through all hits on each plane
			streamlog_out(DEBUG0) << "Entering loop over hits. "<< std::endl;
			position = (_mapHitsVecPerPlane.at(planeID)).at(j)->getPosition();
			if(position == NULL){
				throw(lcio::Exception( "Must exit since one position vector within hit is NULL!"));
			}
			streamlog_out(DEBUG0) << "Here is the information for this hit about to be checked:"<< std::endl;
			streamlog_out(DEBUG0) << "Position: "<<position[0]<<","<<position[1]<<","<< position[2]<< " Plane: "<< j <<std::endl;
			if(j !=0){
				if(positionBefore == NULL){
				throw(lcio::Exception( "Must exit since position before vector within hit is NULL!"));
				}
				streamlog_out(DEBUG0) << "Position: "<<positionBefore[0]<<","<<positionBefore[1]<<","<< positionBefore[2]<< " Plane: "<< j <<std::endl;
				streamlog_out(DEBUG0) << "Position 0 check enter"<< std::endl;
				if(position[0] != positionBefore[0]){//The last condition j != 0 is to ensure we do not compare the last sensor with the next. 
				sensorIsAtLeast1D=true;
				}			
				streamlog_out(DEBUG0) << "Position 1 check enter"<< std::endl;
				if(position[1] != positionBefore[1]){//The last condition j != 0 is to ensure we do not compare the last sensor with the next. 
				sensorIsAtLeast2D=true;
				}			
				streamlog_out(DEBUG0) << "Position 2 check enter"<< std::endl;
				if(position[2] != positionBefore[2]){//The last condition j != 0 is to ensure we do not compare the last sensor with the next. 
				sensorIsAtLeast3D=true;
				}			
			}//END of if j !=0;
				positionBefore= position;
		}	//END of loop over hits
		if(sensorIsAtLeast1D){
			numberOfDimensions++;
		}
		if(sensorIsAtLeast2D){
			numberOfDimensions++;
		}
		if(sensorIsAtLeast3D){
			numberOfDimensions++;
		}
		streamlog_out(DEBUG2) << "The size of Dimensions : "<< numberOfDimensions << std::endl;
		_planeDimensions[planeID] = 2;
	}//END of loop over planes
}	    

void EUTelKalmanFilter::testHitsVecPerPlane(){
	if(_mapHitsVecPerPlane.size() !=  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()){
		streamlog_out(ERROR0) <<"The size of the planes with hits " << _mapHitsVecPerPlane.size() <<"  Sensors from Geometry with no excluded planes: "<<  geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()<<std::endl;
		throw(lcio::Exception(Utility::outputColourString("The number of planes with hits and the number of planes is different", "RED"))); 	

	}
	for(int i=0 ;i<_mapHitsVecPerPlane.size();++i){
		if(_mapHitsVecPerPlane.at(i).size() <= 0){
			streamlog_out(DEBUG0) << Utility::outputColourString("One plane has no hits at all. Is this correct?","YELLOW") << std::endl;
		}
	}
}

void EUTelKalmanFilter::testPlaneDimensions(){
	if(_planeDimensions.size() != geo::gGeometry().sensorIDstoZOrder().size()){
		streamlog_out(ERROR5) << "The size of planesDimensions is: "<< _planeDimensions.size()<<" The size of sensorIDtoOrderZ is: " << geo::gGeometry().sensorIDstoZOrder().size()<< std::endl;
		throw(lcio::Exception( Utility::outputColourString("The size of your dimesion vector is not the same as the number of planes. Something must be wrong!","RED")));
	}
	for(int i=0;i<_planeDimensions.size();++i){
		if(_planeDimensions.at(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i])>3 or _planeDimensions.at(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i])<=0){
			throw(lcio::Exception( "The number of dimension for one of your planes is greater than 3 or less than 0. If this is not a mistake collect you nobel prize now!"));
		}
	}
	for(int i=0;i<_planeDimensions.size();++i){
		if(_planeDimensions.at(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes()[i]) != 2){
			throw(lcio::Exception( "The number of dimensions is not two. If you have a strip sensor remove this test!!!!!!!!   "));
		}
	}

}
void EUTelKalmanFilter::onlyRunOnce(){
	if(_firstExecution){
		setPlaneDimensionsVec();
	}
	_firstExecution=false;
}


    /** Find the hit closest to the intersection of a track with given sensor
     * 
     * @param ts track state
     * @return hit closest to the intersection of the track with the sensor plane
     * 
	 */
const EVENT::TrackerHit* EUTelKalmanFilter::findClosestHit(EUTelTrack & state ) {
	streamlog_out(DEBUG2) << "EUTelKalmanFilter::findClosestHit()" << std::endl;
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
	streamlog_out(DEBUG2) << "----------------------EUTelKalmanFilter::findClosestHit()------------------------" << std::endl;

	return *itClosestHit;
}
//This is not very useful at the moment since the covariant matrix for the hit is guess work at the moment.   
    double EUTelKalmanFilter::getXYPredictionPrecision(EUTelTrack& ts ) const {
      streamlog_out(DEBUG2) << "EUTelKalmanFilter::getXYPredictionPrecision()---BEGIN" << std::endl;

     // TO DO: Need proper error analysis to calculate this rather than providing an answer.  
     // TMatrixDSym Ckkm1(5); 
		//	Ckkm1 = ts.getTrackStateCov();
      double xyPrec = getWindowSize();   //sqrt( Ckkm1[0][0]*Ckkm1[0][0] + Ckkm1[1][1]*Ckkm1[1][1] );
      
      streamlog_out(DEBUG0) << "Minimal combined UV resolution : " << xyPrec << std::endl;
      streamlog_out(DEBUG2) << "----------------------EUTelKalmanFilter::getXYPredictionPrecision()------------------------END" << std::endl;

      return xyPrec;
    }

        
    /**
     * Propagate track state k-1 -> k
     * @param ts track state to update
     */
    void EUTelKalmanFilter::propagateTrackState( EUTelTrack &ts ) {
       /* streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::propagateTrackState()" << std::endl;
        TVectorD xkm1 = ts->getTrackStateVec();
        TVectorD xkkm1 = _jacobianF * xkm1;
        
        streamlog_message( DEBUG0, xkkm1.Print();, std::endl; );
        
        ts->setX( xkkm1[0] );
        ts->setY( xkkm1[1] );
        ts->setTx( xkkm1[2] );
        ts->setTy( xkkm1[3] );
        ts->setInvP( xkkm1[4] );
        
        streamlog_out( DEBUG2 ) << "-----------------------------------EUTelKalmanFilter::propagateTrackState()----------------------------------" << std::endl;
  */  }

	void EUTelKalmanFilter::findNextStateUsingJacobian(EUTelTrack& state,TMatrixD & jacobian,float zPosition, EUTelTrack& newState){
	streamlog_out( DEBUG5 ) << "EUTelKalmanFilter::findNextStateUsingJacobianFinder()" << std::endl;

	streamlog_out( DEBUG5 ) << "The jacobian that will transfrom the state and covariance matrix" << std::endl;
	streamlog_message( DEBUG4, jacobian.Print();, std::endl; );

	TVectorD stateVec = state.getTrackStateVec();
	streamlog_out( DEBUG4 ) << "Before transformation state variables" << std::endl;
	streamlog_message( DEBUG4, stateVec.Print();, std::endl; );
	TVectorD stateVecNew = jacobian * stateVec; 
	streamlog_out( DEBUG4 ) << "After transformation state variables" << std::endl;
	streamlog_message( DEBUG4, stateVecNew.Print();, std::endl; );
	float referencePoint[] = {stateVecNew[3],stateVecNew[4],zPosition};//Z position is not a state parameter but should add here for use later. 
	newState.setPosition(referencePoint);		
///////////////////////////////////////////////////////////Here we update the Covariant matrix
	TMatrixDSym TrackCovState = state.getTrackStateCov();
	streamlog_out( DEBUG4 ) << "Before transformation covMatrix" << std::endl;
	streamlog_message( DEBUG4, TrackCovState.Print();, std::endl; );
	TMatrixDSym TrackCovStateNew  = TrackCovState.Similarity( jacobian );
	streamlog_out( DEBUG4 ) << "After transformation covMatrix" << std::endl;
	streamlog_message( DEBUG4, TrackCovStateNew.Print();, std::endl; ); 

		float trkCov[15] = { static_cast<float>(TrackCovStateNew[0][0]), static_cast<float>(TrackCovStateNew[1][0]), static_cast<float>(TrackCovStateNew[1][1]), 
												 static_cast<float>(TrackCovStateNew[2][0]), static_cast<float>(TrackCovStateNew[2][1]), static_cast<float>(TrackCovStateNew[2][2]),
												 static_cast<float>(TrackCovStateNew[3][0]), static_cast<float>(TrackCovStateNew[3][1]), static_cast<float>(TrackCovStateNew[3][2]),
												 static_cast<float>(TrackCovStateNew[3][3]), static_cast<float>(TrackCovStateNew[4][0]), static_cast<float>(TrackCovStateNew[4][1]),
												 static_cast<float>(TrackCovStateNew[4][2]), static_cast<float>(TrackCovStateNew[4][3]), static_cast<float>(TrackCovStateNew[4][4]) };

	newState.setCovMatrix( trkCov );    
	
	
	streamlog_out( DEBUG5 ) << "-----------------------------------EUTelKalmanFilter::findNextStateUsingJacobianFinder()----------------------------------" << std::endl;
}
				    
    /** Retrieve hit covariance matrix from hit object. Useful for matrix operations
     * 
     * @param hit
     * @return hit covariance matrix
     */
    TMatrixDSym EUTelKalmanFilter::getHitCov( const EVENT::TrackerHit* hit ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getHitCov()" << std::endl;
        const EVENT::FloatVec uvcov = hit->getCovMatrix();
        TMatrixDSym V(2);
        V[0][0] = uvcov[0];                             //cov(x,x)
        V[1][0] = uvcov[1];   V[1][1] = uvcov[2];       //cov(y,x), cov(y,y)
        
        if ( streamlog_level(DEBUG0) ){
            streamlog_out( DEBUG0 ) << "Hit covariance matrix:" << std::endl;
            V.Print();
        }
        
        streamlog_out( DEBUG2 ) << "--------------------------------------EUTelKalmanFilter::getHitCov()-----------------------------------------" << std::endl;
        
        return V;
    }
	
/** Calculate residual vector between given track and hit
* 
* @param ts track state
* @param hit hit
* @return 
*/
TVectorD EUTelKalmanFilter::computeResidual(EUTelTrack& state , const EVENT::TrackerHit* hit ) const {
	streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::computeResidual()" << std::endl;
	const double* hitPosition = hit->getPosition();//In local coordinates

	streamlog_out( DEBUG3 ) << "Hit (id=" << hit->id() << ") local(u,v) coordinates of hit: (" << hitPosition[0] << "," << hitPosition[1] <<","<<hitPosition[2] << ")" << std::endl;

	double localPosition [3];
	//TO DO:Fix master2Localtwo so that is is the only one. This is currently a hack
	const double referencePoint[]	= {state.getReferencePoint()[0], state.getReferencePoint()[1],state.getReferencePoint()[2]};//Need this since geometry works with const doubles not floats 
	geo::gGeometry().master2Localtwo( state.getLocation(), referencePoint, localPosition );

	streamlog_out( DEBUG3 ) << "	Prediction for hit (id=" << hit->id() << ") local(u,v) coordinates of state: ("  << localPosition[0] << "," << localPosition[1] <<","<<localPosition[2] << ")" << std::endl;

	TVectorD residual(2);
	residual[0] = 0 ; residual[1] = 0;
	//This loop is used in case we have a strip sensor. So we should use only 1 dimension of information
	for(int i=0; i<_planeDimensions.at(state.getLocation());++i){
		residual[i] = hitPosition[i] - localPosition[i];
	}	
	
	if ( streamlog_level(DEBUG2) ){
			streamlog_out( DEBUG2 ) << "	Residual vector residual: (" <<residual[0] <<","<<residual[1]<<")"<< std::endl;
			residual.Print();
	}
	
	streamlog_out( DEBUG2 ) << "----------------------------------EUTelKalmanFilter::computeResidual()------------------------------------" << std::endl;
	
	return residual;
}

    /** Retrieve residuals covariance
     * 
     * @param ts track state
     * @param hit
     * @return 
     */
    TMatrixDSym EUTelKalmanFilter::getResidualCov( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getResidualCov()" << std::endl;
        TMatrixD Hk = ts->getH();
        
        _processNoiseQ.Zero();
        TMatrixDSym Ckm1 = ts->getTrackStateCov( ); // 5x5
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF ); // 5x5        //Ckkm1 += _processNoiseQ;       
        TMatrixDSym Rkkm1 = Ckkm1.Similarity(Hk); // 5x5
        Rkkm1 += getHitCov(hit); // hitCov is 2x2 ??
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Residual covariance matrix:" << std::endl;
            Rkkm1.Print();
        }
        
        streamlog_out( DEBUG2 ) << "-----------------------------------EUTelKalmanFilter::getResidualCov()------------------------------------" << std::endl;
        
        return Rkkm1;
    }
    
    /** Retrieve Kalman gain matrix. This is a 5x2 matrix so takes a state 5x1 as an argument.
     * 
     * @param ts track state
     * @param hit
     * 
     * @return Gain matrix K
     */
    TMatrixD EUTelKalmanFilter::updateGainK( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::updateGainK()" << std::endl;

	TMatrixD gainK(5,2);
        _processNoiseQ.Zero(); //Not sure the need for this?
        TMatrixDSym Ckm1 = ts->getTrackStateCov( );
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF ); //Transform covariant matrix from one position to the next      
        TMatrixD Ht(5,2);     Ht = Ht.Transpose( ts->getH() );//This matrix transforms from local to global.
        
        gainK = Ckkm1 * Ht * getResidualCov( ts, hit ).Invert();//take residual of hit and track in local frame. Transform to global. Then multiply by global state covariant matrix
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Gain matrix:" << std::endl;
            gainK.Print();
        }
        
        streamlog_out( DEBUG2 ) << "----------------------------------------EUTelKalmanFilter::updateGainK()------------------------------------" << std::endl;
        
        return gainK;
    }
    
    /** Update track state given new hit
     * 
     * @param ts track state
     * @param hit
     */


		void EUTelKalmanFilter::UpdateStateUsingHitInformation(EUTelTrackStateImpl* input,EVENT::TrackerHit* hit, const TMatrixD& jacobian, TMatrixD & KGain, TMatrixD & HMatrix){
       streamlog_out( DEBUG2 ) << "-----------------------EUTelKalmanFilter::UpdateStateUsingHitInformation()-------------------------------START" << std::endl;
				//Get the residual of the hit and the track and the state vector/////////
       // TVectorD residual = computeResidual( input, hit ); //This is just the components of distance in x and y
        TVectorD state = input->getTrackStateVec( );       
				///////////////////////////////////////////////////////////////////////

				///////////////////////////////////////////////////////First the state. Note that if the hit is accurate you want to change the state to that position
       // state += KGain * residual;   //If the hit is very certain KGain will be unity so the new state position will be the hit position. If the hit is very uncertain then KGain is very small and has no effect. Note matrix (5x1)=+(5x2)*(2x1). 
        input -> setX( state[0] );
        input -> setY( state[1] );
        input -> setTx( state[2] );
        input -> setTy( state[3] );
        input -> setInvP( state[4] );
				///////////////////////////////////////////////////////////////////////////
				
		/////////////////////////////////////////////////////////////////Now the Cov Matrix. Note that if the hit is very accurate it will reduce the Cov matrix
		TMatrixDSym CovMatrix = input->getTrackStateCov();
		TMatrixD I(5,5);
		I.UnitMatrix();
		I -= KGain*HMatrix;  //If we have low uncertainty in the hit then I will be 0. So the covariant will be 0 at this state // Matrix    (5x5)=(5x2)*(2x5). Simply transform gain from
        	TMatrixD newCovMatrix = I*CovMatrix;

	        float trkCov[15] = { static_cast<float>(newCovMatrix[0][0]), static_cast<float>(newCovMatrix[1][0]), static_cast<float>(newCovMatrix[1][1]), 
                             static_cast<float>(newCovMatrix[2][0]), static_cast<float>(newCovMatrix[2][1]), static_cast<float>(newCovMatrix[2][2]),
                             static_cast<float>(newCovMatrix[3][0]), static_cast<float>(newCovMatrix[3][1]), static_cast<float>(newCovMatrix[3][2]),
                             static_cast<float>(newCovMatrix[3][3]), static_cast<float>(newCovMatrix[4][0]), static_cast<float>(newCovMatrix[4][1]),
                             static_cast<float>(newCovMatrix[4][2]), static_cast<float>(newCovMatrix[4][3]), static_cast<float>(newCovMatrix[4][4]) };

        input->setCovMatrix( trkCov );
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				input->setHit(hit);
        streamlog_out( DEBUG2 ) << "-----------------EUTelKalmanFilter::UpdateStateUsingHitInformation()-------------------------------END" << std::endl;

	}

	
	void EUTelKalmanFilter::UpdateTrackUsingHitInformation( EUTelTrackStateImpl* input,const EVENT::TrackerHit* hit, EUTelTrackImpl* track, const TMatrixD& jacobian, TMatrixD & KGain, TMatrixD & HMatrix){
				streamlog_out( DEBUG2 ) << "-----------------EUTelKalmanFilter::UpdateTrackUsingHitInformation()-------------------------------BEGIN" << std::endl;			
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Now determine the chi2 of the track.
				TMatrixDSym newCovMatrix = input->getTrackStateCov( );
        TVectorD residual(2);  // residual = computeResidual( input, hit ); //This is just the components of distance in x and y

        TMatrixD HMatrixTranspose(5,2);      HMatrixTranspose.Transpose( HMatrix );
        TMatrixD newCovMatrixMeas(2,2);  newCovMatrixMeas = HMatrix * newCovMatrix * HMatrixTranspose; //This is the new state covariant matrix in measurements space.
				TMatrixD hitCov(2,2);       hitCov = getHitCov( hit );
				hitCov -= newCovMatrixMeas;
				TMatrixD I(2,2);
				I.UnitMatrix();
				I -= HMatrix*KGain;
				TVectorD change_residual(2); change_residual =  I*residual;

				double chi2 = hitCov.Invert().Similarity(change_residual);
        track->setChi2( chi2 );
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				streamlog_out( DEBUG5 ) << "Add hit: " << hit << " to track " << track << std::endl;
				streamlog_out( DEBUG2 ) << "-----------------EUTelKalmanFilter::UpdateTrackUsingHitInformation()-------------------------------END" << std::endl;
	
	}



    double EUTelKalmanFilter::updateTrackState( EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::updateTrackState()" << std::endl;
        TVectorD rkkm1(2);   //   rkkm1 = computeResidual( ts, hit );
        TVectorD xk(5);         xk = ts->getTrackStateVec( );
        
        TMatrixD Kk = updateGainK( ts, hit );
        xk += Kk * rkkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Updated track parameters:" << std::endl;
            xk.Print();
        }
        
        TMatrixD Hk = ts->getH();
        TMatrixD I(5,5);     I.UnitMatrix();
        I -= Kk*Hk;

        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Track parameters projection matrix Hk:" << std::endl;
            Hk.Print();
            streamlog_out( DEBUG0 ) << "Gain matrix Kk:" << std::endl;
            Kk.Print();
        }
        
        _processNoiseQ.Zero();
        TMatrixDSym Ckm1 = ts->getTrackStateCov( );
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF );        //Ckkm1 += _processNoiseQ;
        
        TMatrixD Ck(5,5);
        Ck = I*Ckkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Updated track covariance matrix:" << std::endl;
            Ck.Print();
        }
        
        ts -> setX( xk[0] );
        ts -> setY( xk[1] );
        ts -> setTx( xk[2] );
        ts -> setTy( xk[3] );
        ts -> setInvP( xk[4] );
        
        float trkCov[15] = { static_cast<float>(Ck[0][0]), static_cast<float>(Ck[1][0]), static_cast<float>(Ck[1][1]), 
                             static_cast<float>(Ck[2][0]), static_cast<float>(Ck[2][1]), static_cast<float>(Ck[2][2]),
                             static_cast<float>(Ck[3][0]), static_cast<float>(Ck[3][1]), static_cast<float>(Ck[3][2]),
                             static_cast<float>(Ck[3][3]), static_cast<float>(Ck[4][0]), static_cast<float>(Ck[4][1]),
                             static_cast<float>(Ck[4][2]), static_cast<float>(Ck[4][3]), static_cast<float>(Ck[4][4]) };

        ts->setCovMatrix( trkCov );
        ts->setLocation( EUTelTrackStateImpl::AtOther );

        const float newPos[] = {static_cast<float>(xk[0]), static_cast<float>(xk[1]), static_cast<float>(ts->getReferencePoint()[2])};
        ts->setReferencePoint( newPos );
        
        double chi2 = 0.;
        
        TVectorD rk(2);
        TMatrixD Ir(2,2);        Ir.UnitMatrix();
        Ir -= Hk * Kk;
        rk = Ir * rkkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Residual vector rk:" << std::endl;
            rk.Print();
        }
        
        TMatrixD Rk(2,2);       Rk = getHitCov( hit );
        TMatrixD HkT(5,2);      HkT.Transpose( Hk );
        TMatrixD HkCkHkT(2,2);  HkCkHkT = Hk * Ck * HkT;
        
        Rk -= HkCkHkT;

        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Residual covariance matrix Rk:" << std::endl;
            Rk.Print();
        }
        
        chi2 = Rk.Invert().Similarity(rk);
        
        streamlog_out( DEBUG0 ) << "Chi2 contribution of the hit: " << chi2 << std::endl;
        
        streamlog_out( DEBUG2 ) << "------------------------------------------EUTelKalmanFilter::updateTrackState()----------------------------------" << std::endl;
        
        return chi2;
    }
    
    IMPL::TrackImpl* EUTelKalmanFilter::cartesian2LCIOTrack( EUTelTrackImpl* track ) const {

        IMPL::TrackImpl* LCIOtrack = new IMPL::TrackImpl;


        int nstates =  track->getTrackStates().size();
        for(int i=0;i < nstates; i++) 
        {
//           EUTelTrackStateImpl* nexttrackstate = (track->getTrackStates().at(i)) ;  
          EUTelTrackStateImpl* nexttrackstate = new EUTelTrackStateImpl( *(track->getTrackStates().at(i)) ); // memory leak source ?
          IMPL::TrackStateImpl* implstate     = static_cast <IMPL::TrackStateImpl*> (nexttrackstate );
          LCIOtrack->addTrackState( implstate );
        }

        // Assign hits to LCIO TRACK
        const EVENT::TrackerHitVec& trkcandhits = track->getTrackerHits();
        EVENT::TrackerHitVec::const_iterator itrHit;
        for ( itrHit = trkcandhits.begin(); itrHit != trkcandhits.end(); ++itrHit ) {
             LCIOtrack->addHit( *itrHit );
        }

        return LCIOtrack;

    }
void EUTelKalmanFilter::findHitsOrderVec(LCCollection* lcCollection,EVENT::TrackerHitVec& hitsOrderVec) {
	for (int iHit = 0; iHit < lcCollection->getNumberOfElements(); iHit++) {
		TrackerHitImpl * hit = static_cast<TrackerHitImpl*> (lcCollection->getElementAt(iHit));

		//TO DO:The determination of HitID is inefficient. Since you create a decoder for each event.
		const int localSensorID = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (hit) );
		if ( localSensorID >= 0 ) hitsOrderVec.push_back( hit );

	} // end loop over all hits in lcCollection

}
void EUTelKalmanFilter::printHits(){
	streamlog_out(MESSAGE0) << "EUTelKalmanFilter::prinitHit: BEGIN ==============" << std::endl;
	EVENT::TrackerHitVec::const_iterator itHit;
	for ( itHit = _allHitsVec.begin() ; itHit != _allHitsVec.end(); ++itHit ) {
		const double* uvpos = (*itHit)->getPosition();
		const int sensorID = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHit) );
		streamlog_out(MESSAGE0) << "Hit (id=" << setw(3) << sensorID << ") local(u,v) coordinates: ("<< setw(7) << setprecision(4) << uvpos[0] << "," << setw(7) << setprecision(4) << uvpos[1] << ")" << std::endl;
	}
streamlog_out(MESSAGE0) << "EUTelKalmanFilter::printHits: END ==============" << std::endl;
}




    MeasurementLayer::MeasurementLayer() : _id(-1), _allHits() {}
    
    MeasurementLayer::MeasurementLayer( int id ) : _id(id), _allHits() {}
    
    MeasurementLayer::~MeasurementLayer() { }
    
    void MeasurementLayer::addHit( EVENT::TrackerHit* hit ) {
        _allHits.push_back( hit );
    }    
    
} // namespace eutelescope
