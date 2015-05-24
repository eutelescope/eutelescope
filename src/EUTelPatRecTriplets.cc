#include "EUTelPatRecTriplets.h"
#include "EUTelNav.h"

namespace eutelescope {

EUTelPatRecTriplets::EUTelPatRecTriplets():  
_totalNumberOfHits(0),
_totalNumberOfSharedHits(0),
_firstExecution(true),
_numberOfTracksTotal(0),
_numberOfTracksAfterHitCut(0),
_numberOfTracksAfterPruneCut(0),
_allowedMissingHits(0),
_tripletSlopeCuts(0,0),
_beamE(-1.),
_beamQ(-1.)
{}
EUTelPatRecTriplets::~EUTelPatRecTriplets()  
{}

void EUTelPatRecTriplets::setPlaneDimensionsVec(EVENT::IntVec planeDimensions){
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

//setRadLengths: This will determine the variance fraction each scatterer will get. Note this comes in two parts. The first is the plane and the next scattering from the air.    
void EUTelPatRecTriplets::setRadLengths(EUTelTrack & track,	std::map<const int,double>  mapSensor, std::map<const int ,double>  mapAir, double rad ){
	//THE FINAL WEIGHT WE HAVE WILL BE A FRACTION PERCENTAGE OF THE TOTAL RADIATION LENGTH
	std::vector<EUTelState>& states = track.getStates();
	const double var  = pow( Utility::getThetaRMSHighland(states.at(0).getMomLocal().Mag(), rad) , 2);
	for(size_t i =0; i < track.getStates().size();++i){ //LOOP over all track again.
		streamlog_out(DEBUG0)<< std::scientific << " Values placed in variance using Highland formula corrected. (SENSOR) : " << (mapSensor[states.at(i).getLocation()]/rad)*var << "  (AIR)  " << (mapAir[states.at(i).getLocation()]/rad)*var <<std::endl;
		states.at(i).setRadFrac((mapSensor[states.at(i).getLocation()]/rad)*var,(mapAir[states.at(i).getLocation()]/rad)*var);//We input the fraction percentage.
	}
	//NOW DETERMINE THE VARIANCE DUE TO THE RADIATION LENGTH. THIS IN THE END WILL BE DIVIDED AMOUNG THE SCATTERERS.
	track.setTotalVariance(var);
}


void EUTelPatRecTriplets::printTrackCandidates(){
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::printTrackCandidates----BEGIN "<< std::endl;
	for(size_t i = 0; i < _tracks.size();++i){
		streamlog_out(DEBUG5)<<"Track number "<<i<<" Out of " <<_tracks.size()<<std::endl; 
		_tracks.at(i).print();
	}
	streamlog_out ( DEBUG1 ) << "EUTelKalmanFilter::printTrackCandidates----END "<< std::endl;
}
void EUTelPatRecTriplets::clearFinalTracks(){
	_finalTracks.clear();
}

//It is important to note what the output of this tells you. If your get that 25% of tracks passed pruning,
//this does not mean that only 25% of events had a track. It simply means that out of all tracks that began 
//from seed hit, only 25% made a full track in the end.
void EUTelPatRecTriplets::testTrackQuality()
{
	_numberOfTracksTotal = _numberOfTracksTotal + _tracks.size();

	if(_numberOfTracksTotal % 1000 == 0)
	{
        streamlog_out(MESSAGE5) << "Number of tracks per event: " << static_cast<float>(_numberOfTracksTotal)/static_cast<float>(getEventNumber() )<< std::endl;
	}
}

//I have observed a difference of just over 6 microns between the two methods over 20 cm distance. So we will set the maximum distance to 10 microns difference
void EUTelPatRecTriplets::testPositionEstimation(float position1[], float position2[]){
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

//This function is used to check that the input data is as expect. If not then we end the processor by throwing a exception.
////TO DO: should check if seed planes are also excluded	
void EUTelPatRecTriplets::testUserInput() {
	streamlog_out(DEBUG2) << "EUTelPatRecTriplets::testUserInput()" << std::endl;
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
}	

void EUTelPatRecTriplets::createTriplets()
{
	 streamlog_out(DEBUG0) << "Create triplets..." << std::endl;

    _tripletsVec.clear();
    float omega = -1.0/_beamE;
	const gear::BField& Bfield = geo::gGeometry().getMagneticField();
	gear::Vector3D vectorGlobal(0.1,0.1,0.1);
	const double Bx = (Bfield.at( vectorGlobal ).x());  
	const double By = (Bfield.at( vectorGlobal ).y());
	const double Bz = (Bfield.at( vectorGlobal ).z());
    TVector3 B(Bx, By, Bz );
    TVector3 H = (B.Unit());
    TVector3 curv = H.Cross(TVector3(0,0,1));
    //Note the cross product
    double curvX = 0.0003*curv[0]*omega; 
    double curvY = 0.0003*curv[1]*omega; 
    std::vector<unsigned int> cenID;
    cenID.push_back(1);
    cenID.push_back(4);
    for(unsigned int i=0 ; i< cenID.size(); i++){
        streamlog_out(DEBUG0) << "Centre sensor ID: " << cenID.at(i)  << std::endl;
        EVENT::TrackerHitVec& hitCentre = _mapHitsVecPerPlane[cenID.at(i)];
        EVENT::TrackerHitVec& hitCentreLeft = _mapHitsVecPerPlane[cenID.at(i) - 1];
        EVENT::TrackerHitVec& hitCentreRight = _mapHitsVecPerPlane[cenID.at(i) + 1];

		EVENT::TrackerHitVec::iterator itHit;
		EVENT::TrackerHitVec::iterator itHitLeft;
		EVENT::TrackerHitVec::iterator itHitRight;
		for ( itHitLeft = hitCentreLeft.begin(); itHitLeft != hitCentreLeft.end(); ++itHitLeft ) {
            const int hitLeftLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHitLeft) );
            double hitLeftPos[] = { (*itHitLeft)->getPosition()[0], (*itHitLeft)->getPosition()[1], (*itHitLeft)->getPosition()[2] };
            double hitLeftPosGlobal[3];
            geo::gGeometry().local2Master(hitLeftLoc ,hitLeftPos,hitLeftPosGlobal);
            EUTelState stateLeft;
            stateLeft.setLocation(hitLeftLoc);
            stateLeft.setPositionGlobal(hitLeftPosGlobal);
            stateLeft.setHit(*itHitLeft);
            for ( itHitRight = hitCentreRight.begin(); itHitRight != hitCentreRight.end(); ++itHitRight ) {
                const int hitRightLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHitRight) );
                double hitRightPos[] = { (*itHitRight)->getPosition()[0], (*itHitRight)->getPosition()[1], (*itHitRight)->getPosition()[2] };
                double hitRightPosGlobal[3];
                geo::gGeometry().local2Master(hitRightLoc ,hitRightPos,hitRightPosGlobal);

                EUTelState stateRight;
                stateRight.setLocation(hitRightLoc);
                stateRight.setPositionGlobal(hitRightPosGlobal);
                stateRight.setHit(*itHitRight);
                doublets doublet;
                doublet = getDoublet(hitLeftPosGlobal,hitRightPosGlobal,curvX,curvY);
                //Add State slopes after doublet creation.
                stateLeft.setMomGlobalIncEne(doublet.slope,getBeamMomentum() );
                stateRight.setMomGlobalIncEne(doublet.slope,getBeamMomentum());
                streamlog_out(DEBUG0) << "Doublet delta X: "<< std::abs(doublet.diff.at(0)) << " Cut X: "<<_doubletDistCut.at(0) <<" Delta Y: " << doublet.diff.at(1)<< " Cut Y: " << _doubletDistCut.at(1)  << std::endl;

                if(fabs(doublet.diff.at(0)) >  _doubletDistCut.at(0) or fabs(doublet.diff.at(1)) >  _doubletDistCut.at(1) ){
                    continue;
                }
                streamlog_out(DEBUG0) << "PASS 1!! " << std::endl;

                //Now loop through all hits on plane between two hits which create doublets. 
                for ( itHit = hitCentre.begin(); itHit != hitCentre.end(); ++itHit ) {
                    const int hitLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHit) );
                    double hitPos[] = { (*itHit)->getPosition()[0], (*itHit)->getPosition()[1], (*itHit)->getPosition()[2] };
                    double hitPosGlobal[3];
                    geo::gGeometry().local2Master(hitLoc ,hitPos,hitPosGlobal);
                    EUTelState state;
                    state.setLocation(hitLoc);
                    state.setPositionGlobal(hitPosGlobal);
                    state.setHit(*itHit);
                    state.setMomGlobalIncEne(doublet.slope,getBeamMomentum());
                    float initDis = geo::gGeometry().getInitialDisplacementToFirstPlane();
                    float x1 = hitPos[0] - 0.5*curvX*pow(hitPos[2] - initDis, 2);
                    float y1 = hitPos[1] - 0.5*curvY*pow(hitPos[2] - initDis, 2);
                //    streamlog_out(DEBUG0) << "Centre Doublet/Hit positions: "<<doublet.pos.at(0)<<"/" <<x1<<"  " <<doublet.pos.at(1)<<"/"<<y1 << std::endl;
                    double delX = doublet.pos.at(0) - x1;
                    double delY = doublet.pos.at(1) - y1;
                    streamlog_out(DEBUG0) << "Doublet delta X to centre hit: "<< delX << " Cut X: "<<_doubletCenDistCut.at(0) <<" Distance Y: " << delY<< " Cut Y: " << _doubletCenDistCut.at(1)  << std::endl;

                    if(fabs(delX) >  _doubletCenDistCut.at(0) or fabs(delY) >  _doubletCenDistCut.at(1) ){
                        continue;
                    }
                    streamlog_out(DEBUG0) << "PASS 2!! " << std::endl;

                    //Set distance
//                    TVector3 vectorDist = state.getPositionGlobal()-stateLeft.getPositionGlobal();
                  //  float dist = sqrt(pow(state.getPositionGlobal()[0]-stateLeft.getPositionGlobal()[0],2)+pow(state.getPositionGlobal()[1]-stateLeft.getPositionGlobal()[1],2)+pow(state.getPositionGlobal()[2]-stateLeft.getPositionGlobal()[2],2));
                    float dist = (state.getPositionGlobal() - stateLeft.getPositionGlobal()).Mag();
                    stateLeft.setArcLengthToNextState(dist);
                    dist = (stateRight.getPositionGlobal() - state.getPositionGlobal()).Mag();
                    state.setArcLengthToNextState(dist);
                    streamlog_out(DEBUG0) << "LEFT: " << cenID.at(i)  << std::endl;
                    stateLeft.print();
                    streamlog_out(DEBUG0) << "CENTRE: " << cenID.at(i)  << std::endl;
                    state.print();
                    streamlog_out(DEBUG0) << "RIGHT: " << cenID.at(i)  << std::endl;
                    stateRight.print();
                    triplets triplet;
                    triplet.states.push_back(stateLeft);
                    triplet.states.push_back(state);
                    triplet.states.push_back(stateRight);
                    triplet.pos.push_back(doublet.pos.at(0));
                    triplet.pos.push_back(doublet.pos.at(1));
                    triplet.pos.push_back(doublet.pos.at(2));
                    triplet.cenPlane = cenID.at(i);
                    triplet.slope.push_back(doublet.slope.at(0));
                    triplet.slope.push_back(doublet.slope.at(1));
                    _tripletsVec.push_back(triplet); 
                }
            }
        }
        if(cenID.size()-1 == i and _tripletsVec.size() == 0){
            streamlog_out(DEBUG0) << "FOUND NO TRIPLETS FOR EVENT: " << getEventNumber()  << std::endl;
            break;
        }
    }
}

std::vector<EUTelTrack> EUTelPatRecTriplets::getTracks( ){
    //Adding all the elements needed for track fitting.
    _tracks.clear();
    //Determine hits which form the track measured by telescope. Cuts are optimum at different values for different geometric setups. 
    createTriplets();
    //Link triplets to form the basis of EUTelTrack objects _tracks. Incomplete state, since scattering information must be added for GBL fit. 
    findTrackFromTriplets();
    _tracksWithDUTHit.clear();
    streamlog_out(DEBUG0) << "TRACKS AFTER: setTrackStatesHits(). Event: " << getEventNumber()  << std::endl;
    for(unsigned int i = 0 ; i < _tracks.size();++i){
        _tracks.at(i).print();
    }
    //If we have dut planes then get hits. If not then just pass mimosa tracks.
    unsigned int  dutNum = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - 6;
    if(dutNum != 0){
        getDUTHit();
    }else{
        _tracksWithDUTHit = _tracks;
    }
    streamlog_out(DEBUG0) << "TRACKS AFTER: getDUTHit(). Event: " << getEventNumber()  << std::endl;
    for(unsigned int i = 0 ; i < _tracksWithDUTHit.size();++i){
        _tracksWithDUTHit.at(i).print();
    }
    //Determine the scattering for the particles trajectory.  
    setScattering();
    streamlog_out(DEBUG0) << "TRACKS AFTER: getScattering(). Event: " << getEventNumber()  << std::endl;
    for(unsigned int i = 0 ; i < _tracksWithDUTHit.size();++i){
        _tracksWithDUTHit.at(i).print();
    }
    setStraightLineFit(); //Must create straight line fit since kink angle estimations are not calculated here.
    streamlog_out(DEBUG0) << "TRACKS AFTER: getStraightLineFit(). Event: " << getEventNumber()  << std::endl;
    for(unsigned int i = 0 ; i < _tracksWithDUTHit.size();++i){
        _tracksWithDUTHit.at(i).print();
    }
    return _tracksWithDUTHit;
}
void EUTelPatRecTriplets::setStraightLineFit(){
    std::vector<EUTelTrack>::iterator itTrack;
    for(itTrack = _tracksWithDUTHit.begin(); itTrack != _tracksWithDUTHit.end();itTrack++){
        EUTelState State1 = itTrack->getStates().at(1);
        EUTelState State2 = itTrack->getStates().at(4);
        double localAveX = (State1.getMomLocalX() + State2.getMomLocalX())/2.0;
        double localAveY = (State1.getMomLocalY() + State2.getMomLocalY())/2.0;
        double localAveZ = (State1.getMomLocalZ() + State2.getMomLocalZ())/2.0;
        itTrack->getStates().at(0).setMomLocalX(localAveX);
        itTrack->getStates().at(0).setMomLocalY(localAveY);
        itTrack->getStates().at(0).setMomLocalZ(localAveZ);

        for(unsigned int i =0 ; i < itTrack->getStates().size()-1; i++){
            float intersectionPoint[3];
            TVector3 momentumAtIntersection;
            float arcLength;
            int holder; //This is used to return the plane which is found.
            bool found =itTrack->getStates().at(i).findIntersectionWithCertainID(itTrack->getStates().at(i+1).getLocation(), intersectionPoint, momentumAtIntersection,arcLength,holder );
            itTrack->getStates().at(i).setArcLengthToNextState(arcLength);
            itTrack->getStates().at(i+1).setPositionGlobal(intersectionPoint);
			itTrack->getStates().at(i+1).setLocalMomentumGlobalMomentum(momentumAtIntersection); 


        }
    }

}

void EUTelPatRecTriplets::setScattering(){
    std::vector<EUTelTrack>::iterator itTrack;
    for(itTrack = _tracksWithDUTHit.begin(); itTrack != _tracksWithDUTHit.end();itTrack++){
        std::map<const int,double>  mapSensor;
        std::map<const int ,double>  mapAir;
        double rad =	itTrack->getStates().at(0).computeRadLengthsToEnd(mapSensor, mapAir);
        if(rad == 0){
            throw(std::string("Radiation length is zero."));
        }
        itTrack->print();
        setRadLengths(*itTrack, mapSensor, mapAir, rad);

    }
}

void EUTelPatRecTriplets::findTrackFromTriplets(){
    streamlog_out(DEBUG0) << "Set track states and hits... " << std::endl;

    std::vector<triplets>::iterator itTriplet;
    std::vector<triplets> leftTriplets;
    std::vector<triplets> rightTriplets;

    for(itTriplet = _tripletsVec.begin();itTriplet != _tripletsVec.end();  itTriplet++){
        if(itTriplet->cenPlane == 1 ){
            leftTriplets.push_back(*itTriplet);
        }else if(itTriplet->cenPlane == 4){
            rightTriplets.push_back(*itTriplet);
        }else{
            throw(lcio::Exception( "Triplet are not from the left anre right arms!"  ));
        }
    }
    std::vector<triplets>::iterator itLeftTriplet;
    std::vector<triplets>::iterator itRightTriplet;
    streamlog_out(DEBUG0) << "Use triplets... " << std::endl;
    for(itLeftTriplet = leftTriplets.begin();itLeftTriplet != leftTriplets.end();  itLeftTriplet++){
        for(itRightTriplet = rightTriplets.begin();itRightTriplet != rightTriplets.end();  itRightTriplet++){
                    streamlog_out(DEBUG0) << "Triplet slope delta match Cut: " <<"X delta: " << fabs(itRightTriplet->slope.at(0) - itLeftTriplet->slope.at(0)) <<" Cut: " << _tripletSlopeCuts.at(0) << " Y delta: " <<  fabs(itRightTriplet->slope.at(1) - itLeftTriplet->slope.at(1))<<" Cut: " <<  _tripletSlopeCuts.at(1)  << std::endl;

            if(fabs(itRightTriplet->slope.at(0) - itLeftTriplet->slope.at(0)) > _tripletSlopeCuts.at(0)  or fabs(itRightTriplet->slope.at(1) - itLeftTriplet->slope.at(1)) >_tripletSlopeCuts.at(1)  ){
                continue;
            }
            streamlog_out(DEBUG0) << "PASS 3!! " << std::endl;
            float aveZPosTrip = (itLeftTriplet->pos.at(2)+ itRightTriplet->pos.at(2))/2.0;
            streamlog_out(DEBUG0) << "Triplet Propagation..." <<std::endl;   
            streamlog_out(DEBUG0) << "LEFT" <<std::endl;   
            std::vector<float> posLeftAtZ = getTripPosAtZ(*itLeftTriplet,aveZPosTrip); 
            streamlog_out(DEBUG0) << "RIGHT" <<std::endl;   
            std::vector<float> posRightAtZ = getTripPosAtZ(*itRightTriplet,aveZPosTrip); 
            streamlog_out(DEBUG0) << "Predicted position of triplet1/triplet2 " << posLeftAtZ.at(0) <<"/"<<posRightAtZ.at(0) << "  " << posLeftAtZ.at(1) <<"/"<<posRightAtZ.at(1)<< "  " << posLeftAtZ.at(2) <<"/"<<posRightAtZ.at(2) <<std::endl;
            streamlog_out(DEBUG0) << "Delta between Triplets X: "<< fabs(posLeftAtZ.at(0)- posRightAtZ.at(0)) <<" Cut X: " << _tripletConnectDistCut.at(0) << " Delta Y: " << fabs(posLeftAtZ.at(1)- posRightAtZ.at(1)) << " Cut Y: " << _tripletConnectDistCut.at(1) << std::endl;
            if(fabs(posLeftAtZ.at(0)- posRightAtZ.at(0)) > _tripletConnectDistCut.at(0) or fabs(posLeftAtZ.at(1)- posRightAtZ.at(1)) > _tripletConnectDistCut.at(1)){
                continue;
            }
            //Pass without DUT
            streamlog_out(DEBUG0) << "PASS 4!! " << std::endl;
            float dist = (itRightTriplet->states.at(0).getPositionGlobal() - itLeftTriplet->states.at(2).getPositionGlobal()).Mag();
            itLeftTriplet->states.at(2).setArcLengthToNextState(dist);
            _tracks.push_back(getTrack(*itLeftTriplet,*itRightTriplet));

        }
    }
}
EUTelTrack EUTelPatRecTriplets::printTrack(std::vector<EUTelTrack>& tracks){
    for(unsigned int i = 0; i < _tracks.size(); ++i){
        _tracks.at(i).print();
    }
}
EUTelTrack EUTelPatRecTriplets::getTrack(triplets tripLeft,triplets tripRight){
    std::vector<EUTelState> states = tripLeft.states;
    states.insert( states.end(), tripRight.states.begin(), tripRight.states.end() );
    std::vector<EUTelState>::iterator itState;
    std::vector<float> slope;
//    std::cout<<"Hit position "<<states.at(5).getHit().getPosition()[0]<<std::endl; 
    slope.push_back((states.at(5).getHitPositionGlobal()[0]-states.at(0).getHitPositionGlobal()[0])/(states.at(5).getHitPositionGlobal()[2]-states.at(0).getHitPositionGlobal()[2]));
    slope.push_back((states.at(5).getHitPositionGlobal()[1]-states.at(0).getHitPositionGlobal()[1])/(states.at(5).getHitPositionGlobal()[2]-states.at(0).getHitPositionGlobal()[2]));

    EUTelTrack track;
    for(itState = states.begin();itState != states.end();  itState++){
        float intersectionPoint[3];
        TVector3 momentumAtIntersection;
        float arcLength;
        int holder; //This is used to return the plane which is found.
        itState->setMomGlobalIncEne(slope,getBeamMomentum());
        itState->print();
        if(itState != (states.end()-1)){
            bool found =itState->findIntersectionWithCertainID((itState+1)->getLocation() , intersectionPoint, momentumAtIntersection,arcLength,holder );
            (itState+1)->setPositionGlobal(intersectionPoint);
        }
        (itState)->setArcLengthToNextState(arcLength);

        track.setState(*itState);
    }
    return track;

}
void EUTelPatRecTriplets::getDUTHit(){
    float dist;
    float distBest=10000;
    int hitDUTLocBest;
    double hitDUTPosBest[3];
    double hitDUTPosGlobalBest[3];
    float intersectionPointBest[3];
    float arcLengthBest;
    TVector3 momentumAtIntersectionBest;
    EUTelHit hit;

    std::vector<EUTelTrack>::iterator itTrack;
    for(itTrack = _tracks.begin(); itTrack != _tracks.end();itTrack++){
        unsigned int  dutNum = geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - 6;
        for(unsigned int i = 0 ; i < dutNum ; i++){
            float intersectionPoint[3];
            TVector3 momentumAtIntersection;
            float arcLength;
            int holder; //This is used to return the plane which is found.
            //Propagate always from the plane closest to the DUTs
            itTrack->print();
            bool found =itTrack->getStates().at(2).findIntersectionWithCertainID(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(3+i), intersectionPoint, momentumAtIntersection,arcLength,holder );
            if(!found){
                streamlog_out(DEBUG0) << "Can not find estimated intersection on DUT!"  << std::endl;
                continue;
            }else{
                streamlog_out(DEBUG0) << "INTERSECTION FOUND ON PLANE: "<< holder  << std::endl;
            }
            EVENT::TrackerHitVec& hitDUT = _mapHitsVecPerPlane[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(3+i)];
            EVENT::TrackerHitVec::iterator itHitDUT;
            for ( itHitDUT = hitDUT.begin(); itHitDUT != hitDUT.end(); ++itHitDUT ) {
                streamlog_out(DEBUG0) << "Hit on DUT "  << std::endl;
                const int hitDUTLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHitDUT) );
                double hitDUTPos[] = { (*itHitDUT)->getPosition()[0], (*itHitDUT)->getPosition()[1], (*itHitDUT)->getPosition()[2] };
                double hitDUTPosGlobal[3];
                geo::gGeometry().local2Master(hitDUTLoc ,hitDUTPos,hitDUTPosGlobal);
                if(_planeDimensions[hitDUTLoc] == 2){ 
                    streamlog_out(DEBUG0) << "Triplet DUT Match Cut Pixel: " <<"X delta: " << fabs(intersectionPoint[0]-hitDUTPosGlobal[0]) << " Y delta: " << fabs(intersectionPoint[1]- hitDUTPosGlobal[1]) << std::endl;
                    dist = sqrt(pow(intersectionPoint[0]-hitDUTPosGlobal[0],2)+pow(intersectionPoint[1]-hitDUTPosGlobal[1],2));
                }else if(_planeDimensions[hitDUTLoc] == 1){
                    streamlog_out(DEBUG0) << "Triplet DUT Match Cut Strip: " <<"X delta: " << fabs(intersectionPoint[0]-hitDUTPosGlobal[0]) <<" Cut: " <<  _tripletConnectDistCut.at(0) << std::endl;
                    dist = sqrt(pow(intersectionPoint[0]-hitDUTPosGlobal[0],2));

                }else{
                    throw(lcio::Exception( "This is not a strip or pixel sensor!"));
                }
                //Will enter always on the first loop if reached. Need to find intersection and hit on DUT.
                if(dist < distBest){
                    streamlog_out(DEBUG0) << "Save hit information"  << std::endl;
                    if(itHitDUT == hitDUT.begin()){
                        streamlog_out(DEBUG0) << "First hit used as initial DUT hit" << std::endl;
                    }
                    if(dist < distBest){
                        streamlog_out(DEBUG0) << "Improvement! Use this hit." << std::endl;
                    }
                    streamlog_out(DEBUG0) << "Dist: "<< dist << std::endl;
                    streamlog_out(DEBUG0) << "Arclength: "<< arcLength << std::endl;
                    streamlog_out(DEBUG0) << "Location: "<< hitDUTLoc << std::endl;
                    streamlog_out(DEBUG0) << "Position Global: "<< intersectionPoint[0] << " " << intersectionPoint[1] << "  " <<intersectionPoint[2] <<std::endl;
                    streamlog_out(DEBUG0) << "Momentum Global: "<< momentumAtIntersection[0] << " " << momentumAtIntersection[1] << "  " <<momentumAtIntersection[2] <<std::endl;
                    distBest = dist;
                    arcLengthBest=arcLength;
                    hitDUTLocBest = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHitDUT) );
                    hitDUTPosBest[0] = (*itHitDUT)->getPosition()[0]; hitDUTPosBest[1] = (*itHitDUT)->getPosition()[1];hitDUTPosBest[2] = (*itHitDUT)->getPosition()[2];
                    hitDUTPosGlobalBest[0] =  hitDUTPosGlobal[0];hitDUTPosGlobalBest[1] =  hitDUTPosGlobal[1];hitDUTPosGlobalBest[2] =  hitDUTPosGlobal[2];
                    intersectionPointBest[0] = intersectionPoint[0];intersectionPointBest[1] = intersectionPoint[1];intersectionPointBest[2] = intersectionPoint[2];
                    momentumAtIntersectionBest = momentumAtIntersection;
                    hit = *itHitDUT; 
                }

            }
            if(hitDUT.size() != 0 ){
                streamlog_out(DEBUG0) << "ADD DUT HIT! " << std::endl;

                EUTelState stateDUT;
                stateDUT.setLocation(hitDUTLocBest);
                stateDUT.setArcLengthToNextState(arcLengthBest);
                stateDUT.setPositionGlobal(intersectionPointBest);
                stateDUT.setLocalMomentumGlobalMomentum(momentumAtIntersectionBest);
                stateDUT.setHit(&hit);
                EUTelTrack track = getTrackDUTHit(itTrack,stateDUT);
                //DUT hit will be added at position 3. Propagate forward to next mimosa plane. 
                found =itTrack->getStates().at(3).findIntersectionWithCertainID(itTrack->getStates().at(4).getLocation(), intersectionPoint, momentumAtIntersection,arcLength,holder );
                if(!found){
                    continue;
                    streamlog_out(DEBUG0) << "FAIL AT ADDING TRACK! Can not connect to mimosa plane." << std::endl;
                }
                track.getStates().at(3).setArcLengthToNextState(arcLength);
                _tracksWithDUTHit.push_back(track);
            }
        }
    }
}

EUTelTrack EUTelPatRecTriplets::getTrackDUTHit(std::vector<EUTelTrack>::iterator itTrack, EUTelState stateDUT ){
    EUTelTrack track;
    std::vector<EUTelState>::iterator itState;
    for(itState = itTrack->getStates().begin();itState != itTrack->getStates().end();  itState++){
        if(itState->getLocation() == 3){
            track.setState(stateDUT);
        }
        track.setState(*itState);
    }
    return track;
}

std::vector<float>  EUTelPatRecTriplets::getTripPosAtZ(triplets trip, float posZ ){
	streamlog_out(DEBUG0) << "Slope x/y: " <<trip.slope.at(0) << "  " <<trip.slope.at(1) <<" Position ave: " <<trip.pos.at(0)<<" "<<trip.pos.at(1)<<" "<<trip.pos.at(2) <<std::endl;
    float dz = posZ - trip.pos.at(2);
    float x = trip.pos.at(0) + trip.slope.at(0)*dz;
    float y = trip.pos.at(1) + trip.slope.at(1)*dz;
    std::vector<float> position;
    position.push_back(x);position.push_back(y);position.push_back(posZ);
    return position;
}


EUTelPatRecTriplets::doublets EUTelPatRecTriplets::getDoublet(double hitLeftPos[3], double hitRightPos[3],double curvX,double curvY )
{
    float initDis = geo::gGeometry().getInitialDisplacementToFirstPlane();
    //Remove the curvature as a factor between hits. Therefore only the slope will displace the hit position from plane to plane.
    float x1 = hitLeftPos[0] - 0.5*curvX*pow(hitLeftPos[2] - initDis, 2);
    float y1 = hitLeftPos[1] - 0.5*curvY*pow(hitLeftPos[2] - initDis, 2);
    float x2 = hitRightPos[0] - 0.5*curvX*pow(hitRightPos[2] - initDis, 2);
    float y2 = hitRightPos[1] - 0.5*curvY*pow(hitRightPos[2] - initDis, 2);

    doublets doublet;
    doublet.pos.push_back((x2 + x1)/2.0);
    doublet.pos.push_back((y2 + y1)/2.0);
    doublet.pos.push_back((hitRightPos[2] + hitLeftPos[2])/2.0);

    doublet.diff.push_back(x2 - x1);
    doublet.diff.push_back(y2 - y1);
    doublet.diff.push_back( hitRightPos[2] - hitLeftPos[2]);

    doublet.slope.push_back( doublet.diff.at(0)/doublet.diff.at(2)); 
    doublet.slope.push_back( doublet.diff.at(1)/doublet.diff.at(2)); 

    return doublet;
}

TVector3 EUTelPatRecTriplets::computeInitialMomentumGlobal(){
	//We assume that the arc length is the displacement in the z direction. The assumption should be a valid one in most cases
	TVector3 position(0,0,0);//The position we start from does not matter since the magnetic field is homogeneous.
	TVector3 momentum(0,0,_beamE);//Assume the beam starts in a straight line
	float arcLength= geo::gGeometry().getInitialDisplacementToFirstPlane();
	TVector3 momentumEnd = EUTelNav::getMomentumfromArcLength(momentum, _beamQ, arcLength);
	streamlog_out(DEBUG2) << "Momentum on the first sensor: px,py,pz "<<momentumEnd[0] <<","<<momentumEnd[1]<<","<<momentumEnd[2]<<","<<"At an arc length of "<<arcLength<<std::endl;
	return momentumEnd;
}


//We need to delete the states AND the track. Since we have allocated memory to store the state.
void EUTelPatRecTriplets::clearTrackAndTrackStates(){
    _tracks.clear();
}

//This creates map between plane ID and hits on that plane. 
//We also order the map correcly with geometry.
void EUTelPatRecTriplets::setHitsVecPerPlane()
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

void EUTelPatRecTriplets::testHitsVecPerPlane(){
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
        
void EUTelPatRecTriplets::printHits()
{
	streamlog_out(MESSAGE0) << "EUTelPatRecTriplets::prinitHit: BEGIN ==============" << std::endl;
	EVENT::TrackerHitVec::const_iterator itHit;
	for(itHit = _allHitsVec.begin(); itHit != _allHitsVec.end(); ++itHit )
	{
		const double* uvpos = (*itHit)->getPosition();
		const int sensorID = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHit) );
		streamlog_out(MESSAGE0)	<< "Hit (id=" << std::setw(3) << sensorID << ") local(u,v) coordinates: (" 
					<< std::setw(7) << std::setprecision(4) << uvpos[0] << ", " << std::setw(7) 
					<< std::setprecision(4) << uvpos[1] << ")" << std::endl;
	}
	streamlog_out(MESSAGE0) << "EUTelPatRecTriplets::printHits: END ==============" << std::endl;
}

} // namespace eutelescope
