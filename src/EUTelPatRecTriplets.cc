#include "EUTelPatRecTriplets.h"
namespace eutelescope {

EUTelPatRecTriplets::EUTelPatRecTriplets():  
_totalNumberOfHits(0),
_totalNumberOfSharedHits(0),
_firstExecution(true),
_numberOfTracksTotal(0),
_numberOfTracksTotalWithDUT(0),
_tracksWithoutHit(0),
_numberTripletsLeft(0),
_numberTripletsRight(0),
_allowedMissingHits(0),
_tripletSlopeCuts(0,0),
_beamE(-1.),
_beamQ(-1.)
{}
EUTelPatRecTriplets::~EUTelPatRecTriplets()  
{}

// This function will determine the planes which are included in the final track. 
/// To exclude a plane but still create a track then increase resolution to infinity.
/// The hit will therefore have no pull on the track.
/// I.e what planes will get position and incidence information
/// Radiation length is still included in the scattering calculation
void EUTelPatRecTriplets::setPlaneExclude(IntVec& planeIDs){   
    _senNotExc = planeIDs;
    unsigned  int counter=0;
    for(std::vector<int>::const_iterator itID = geo::gGeometry().sensorIDsVec().begin(); itID != geo::gGeometry().sensorIDsVec().end(); ++itID){
        bool excluded=false;
        for(size_t j =0; j< planeIDs.size(); ++j){
            if(*itID == planeIDs[j]){
                excluded=true;
                break;
            }   
        }   
        if(!excluded){
            _senZOrderToIDWithoutExcPla[counter]= *itID;
            counter++;
        }   
    }   
    //Check if the number of excluded planes set is the same as (total-number of plane IDs inputed that should be excluded)
    if(_senZOrderToIDWithoutExcPla.size() != (geo::gGeometry().sensorIDsVec().size()-planeIDs.size())){
        throw(lcio::Exception( "The number of Planes-Excluded is not correct. This could be a problem with geometry."));
    }else{
        streamlog_out(DEBUG0) <<"The correct number of planes have been excluded" << std::endl;
    }   
}  

void EUTelPatRecTriplets::setPlaneDimensionsVec(EVENT::IntVec& planeDimensions){
	if(planeDimensions.size() != geo::gGeometry().sensorIDsVec().size()){
		streamlog_out(ERROR) << "The size of planesDimensions input is: "<< planeDimensions.size()<<" The size of sensorIDsVec is: " << geo::gGeometry().sensorIDsVec().size()<< std::endl;
		throw(lcio::Exception( "The input dimension vector not the same as the number of planes!"));
	}
	_planeDimensions.clear();
	for(size_t i=0; i<geo::gGeometry().sensorIDsVec().size(); ++i){
		const int planeID = geo::gGeometry().sensorIDsVec().at(i);
		if (_planeDimensions.find(planeID) == _planeDimensions.end()){//This is to check that we don't try to map the same sensor to two different plane dimensions.
			_planeDimensions[planeID] = planeDimensions.at(i);
		}else{
			streamlog_out(ERROR5) <<"The z position is : "<< i <<" The plane ID is: " << planeID <<std::endl;
			throw(lcio::Exception( "You are trying to map the same sensor ID to two different plane dimensions. There is something wrong with you gear file input. Make sure there is some distance between your planes in the gear file!"));
		}
	}
}	    

/*setRadLengths: This will determine the variance fraction each scatterer will get. Note this comes in two parts. The first is the plane and the next scattering from the air.   */
void EUTelPatRecTriplets::setRadLengths(EUTelTrack & track,	std::map<const int,double>&  mapSensor, std::map<const int ,double>&  mapAir, double & rad ){
	//THE FINAL WEIGHT WE HAVE WILL BE A FRACTION PERCENTAGE OF THE TOTAL RADIATION LENGTH
	std::vector<EUTelState>& states = track.getStates();
	const double var  = pow( Utility::getThetaRMSHighland(track.getBeamEnergy(), rad) , 2);
	for(size_t i =0; i < track.getStates().size();++i){ 
		streamlog_out(DEBUG0)<< std::scientific << " Values placed in variance using Highland formula corrected. (SENSOR) : " << (mapSensor[states.at(i).getLocation()]/rad)*var << "  (AIR)  " << (mapAir[states.at(i).getLocation()]/rad)*var <<std::endl;
		states.at(i).setRadFrac((mapSensor[states.at(i).getLocation()]/rad)*var,(mapAir[states.at(i).getLocation()]/rad)*var);//We input the fraction percentage.
	}
	track.setTotalVariance(var);
}



void EUTelPatRecTriplets::testTrackQuality(std::vector<EUTelTrack>&  tracksWithDUTs )
{
	_numberOfTracksTotal = _numberOfTracksTotal + tracksWithDUTs.size();
	_numberOfTracksTotalWithDUT = _numberOfTracksTotalWithDUT + _tracksWithDUTHit.size();

	if(_numberOfTracksTotal % 5000 == 0){
        streamlog_out(MESSAGE5) << "Percentage tracks without DUT hit: " << static_cast<float>(_tracksWithoutHit)/static_cast<float>(_numberOfTracksTotal)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of tracks per event: " << static_cast<float>(_numberOfTracksTotal)/static_cast<float>(getEventNumber() +1)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of left arm triplets per event: " << static_cast<float>(_numberTripletsLeft)/static_cast<float>(getEventNumber() +1)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of right arm triplets per event: " << static_cast<float>(_numberTripletsRight)/static_cast<float>(getEventNumber() +1)<< std::endl;
    }
}
void EUTelPatRecTriplets::testUserInput() {
	streamlog_out(DEBUG2) << "EUTelPatRecTriplets::testUserInput()" << std::endl;
	if ( _beamE < 1.E-6 ) {
		throw(lcio::Exception( "Beam direction was set incorrectly")); 
	}
	else{
	 streamlog_out(DEBUG1) << "Beam energy is reasonable" << std::endl;
	}
}	
///Curvature in x/y is passed to naviagation.
std::vector<double>  EUTelPatRecTriplets::getCurvXY(){
    //Defined the same as saved in track parameters.
    const float omega = -1.0/_beamE;
    TVector3 bFac = getBFac();
    streamlog_out(DEBUG0) << "BFac field unit: " << bFac[0] << "  " << bFac[1] <<"  "<< bFac[2] << "  Omega: " << omega << std::endl;
    //Note the cross product
    const double curvX = bFac[0]*omega; 
    const double curvY = bFac[1]*omega; 
    std::vector<double> curv;
    curv.push_back(curvX); curv.push_back(curvY);
    //streamlog_out(DEBUG0) << "The curvature calculated: " << curv.at(0) << "  " << curv.at(1) << std::endl;

    return curv;

}
///cZxB is passed to navigation.
TVector3  EUTelPatRecTriplets::getBFac(){
    const gear::BField& Bfield = geo::gGeometry().getMagneticField();
    gear::Vector3D vectorGlobal(0.1,0.1,0.1);
    const double Bx = (Bfield.at( vectorGlobal ).x());  
    const double By = (Bfield.at( vectorGlobal ).y());
    const double Bz = (Bfield.at( vectorGlobal ).z());
    TVector3 B(Bx, By, Bz );
    TVector3 H = (B.Unit());
    TVector3 bFac = 0.0003*(TVector3(0,0,1).Cross(H));
    return bFac;
}


std::vector<EUTelPatRecTriplets::triplets> EUTelPatRecTriplets::getTriplets()
{
	 streamlog_out(DEBUG1) << "Create triplets..." << std::endl;
    std::vector<triplets> tripletsVec;
    const double curvX = EUTelNav::_curv.at(0); 
    const double curvY =  EUTelNav::_curv.at(1); 
    std::vector<unsigned int> cenID;
    cenID.push_back(1);
    cenID.push_back(4);
    for(size_t i=0 ; i< cenID.size(); i++){
        streamlog_out(DEBUG1) << "Centre sensor ID: " << cenID.at(i)  << std::endl;
        std::vector<EUTelHit>& hitCentre = _mapHitsVecPerPlane[cenID.at(i)];
        std::vector<EUTelHit>& hitCentreLeft = _mapHitsVecPerPlane[cenID.at(i) - 1];
        std::vector<EUTelHit>& hitCentreRight = _mapHitsVecPerPlane[cenID.at(i) + 1];

		std::vector<EUTelHit>::iterator itHit;
		std::vector<EUTelHit>::iterator itHitLeft;
		std::vector<EUTelHit>::iterator itHitRight;
		for ( itHitLeft = hitCentreLeft.begin(); itHitLeft != hitCentreLeft.end(); ++itHitLeft ) {
            for ( itHitRight = hitCentreRight.begin(); itHitRight != hitCentreRight.end(); ++itHitRight ) {
                doublets doublet;
                bool pass = getDoublet(*itHitLeft,*itHitRight,doublet);
                if(!pass) continue;
                streamlog_out(DEBUG1) << "PASS 1!! " << std::endl;

                //Now loop through all hits on plane between two hits which create doublets. 
                for ( itHit = hitCentre.begin(); itHit != hitCentre.end(); ++itHit ) {
                    triplets triplet
                    bool pass2 = getTriplet(doublet,*itHit, triplet);
                    if(!pass) continue;
                    streamlog_out(DEBUG1) << "PASS 2!! " << std::endl;
                     if(triplet.cenPlane == 1){
                        _numberTripletsLeft++;
                    }else if(triplet.cenPlane == 4){
                        _numberTripletsRight++;
                    }
                    tripletsVec.push_back(triplet); 
                }
            }
        }
        if(cenID.size()-1 == i and _tripletsVec.size() == 0){//If not triplets found on first plane end search.
            streamlog_out(DEBUG1) << "FOUND NO TRIPLETS FOR EVENT: " << getEventNumber()  << std::endl;
            break;
        }
    }
    return tripletsVec;

}
bool EUTelPatRecTriplets::getTriplet(doublets& doublet, EUTelHit& hitCen){

    const float initDis = geo::gGeometry().getInitialDisplacementToFirstPlane();
    const float x1 = hitCen.getPositionGlobal()[0] - 0.5*curvX*pow(hitCen.getPositionGlobal()[2] - initDis, 2);
    const float y1 = hitCen.getPositionGlobal()[1] - 0.5*curvY*pow(hitCen.getPositionGlobal()[2] - initDis, 2);
    std::vector<double> pos = getDoubPosAtZ(doublet, itHit->getPositionGlobal()[2]);
    const double delX = pos.at(0) - x1;
    const double delY = pos.at(1) - y1;
    streamlog_out(DEBUG1) << "Doublet delta X to centre hit: "<< delX << " Cut X: "<<_doubletCenDistCut.at(0) <<" Distance Y: " << delY<< " Cut Y: " << _doubletCenDistCut.at(1)  << std::endl;

    if(fabs(delX) >  _doubletCenDistCut.at(0) or fabs(delY) >  _doubletCenDistCut.at(1) ){
        return false;
    }
    streamlog_out(DEBUG1) << "PASS 2!! " << std::endl;

    triplets triplet;
    triplet.matches = 0;
    triplet.fitID = -999;
    triplet.pos.push_back(doublet.pos.at(0));
    triplet.pos.push_back(doublet.pos.at(1));
    triplet.pos.push_back(doublet.pos.at(2));
    ///Used to find final tracks.
    triplet.cenPlane = hitCen.getLocation();
    triplet.slope.push_back(doublet.slope.at(0));
    triplet.slope.push_back(doublet.slope.at(1));
    triplet.hits =  std::vector<EUTelHit>(doublet.hits.at(0),hitCen,doublet.hits.at(1));
    streamlog_out(DEBUG1) << "triplet created." << std::endl;
    return true;

}

/*This is the function which should be used in processEvent
*/
std::vector<EUTelTrack> EUTelPatRecTriplets::getTracks( ){
    std::vector<EUTelTrack>  tracks;
    bool _getTriplets = true;
    if(_getTriplets){
        EUTelNav::_intBeamE = getBeamMomentum();
        std::vector<EUTelPatRecTriplets::triplets> tripletVec = getTriplets();
        streamlog_out(DEBUG1) << "Total number of triplets found:  " << tripletVec.size()  << std::endl;
        std::map<int,std::vector<EUTelHit> > id_Hits  = getTrackHitsFromTriplets(tripletVec);
        for(std::map<int, vector<EUTelHit> > itID_Hit = idHits.begin(); itID_Hit != id_Hits.end();itID_Hit){
            doublets doub;
            bool pass = getDoublet(*(itID_Hit->second.begin()),*(itID_Hit->second.rbegin()),doub);
            std::vector<EUTelHit> newHits;
            std::vector<int> dut;
            for(size_t i = 0 ; i < _senNotExc.size(); ++i){
                if(_senNotExc.at(i) > 5){
                    dut.push_back(_senNotExc.at(i));
                }
            }
            bool pass =  getDoubHitOnTraj(doub, dut, newHits);
            newHits.insert( newHits.end(), itID_Hit->second.begin(), itID_Hit->second.end() );
            ///Hit order does not matter
            tracks.push_back(getTrack(newHits));
        }
    }else{
        ///Basic pattern recognition used to find DUT hit.
        /// Can also use basic pattern recognition to find more tracks than triplets.
    }
    streamlog_out(DEBUG1) << tracksWithDUTs.size()<<" TRACKS: " << getEventNumber()  << std::endl;
    for(size_t i = 0 ; i < tracks.size();++i){
        tracks.at(i).print();
    }
    streamlog_out( DEBUG1 ) << "tracks.size() = "<<tracks.size()<<std::endl;
    testTrackQuality(tracks );
    return tracksWithDUTs;
}
///MUST USE A MIMIMUM OF FOUR HITS!!!
EUTelTrack EUTelPatRecTriplets::getTrackFourHits(std::vector<EUTelHit> hits ){
	streamlog_out(DEBUG1) << "HITS TO FORM TRACK FROM: " << std::endl;
    std::vector<EUTelHit>::iterator itHit;
    for(itHit = hits.begin(); itHit != hits.end(); itHit++){
        itHit->print(); 
    }
    //Always use mimosa planes to create initial track parameterisation.
    EUTelHit hitArmOne1;
    EUTelHit hitArmOne2;
    EUTelHit hitArmTwo1;
    EUTelHit hitArmTwo2;

    for(itHit = hits.begin(); itHit != hits.end(); itHit++){
        if(itHit->getLocation() == 0){
            hitArmOne1 = *itHit;
        }
        if(itHit->getLocation() == 2){
            hitArmOne2 = *itHit;
        }
        if(itHit->getLocation() == 3){
            hitArmTwo1 = *itHit;
        }
        if(itHit->getLocation() == 5){
            hitArmTwo2 = *itHit;
        }

    }
    double curvCorr;
    //Find correction of curvature through slope change.
    const gear::BField& B = geo::gGeometry().getMagneticField();
    const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();
    if ( Bmag < 1.E-6 ){
        curvCorr = 0;
    }else{
        curvCorr =  EUTelNav::getCorr(hitArmOne1,hitArmOne2,hitArmTwo1,hitArmTwo2);
    }
    //NOW CREATE TRACK CANDIDATE
    std::vector<double> offset;
    std::vector<double> trackSlope; 
    EUTelNav::getTrackAvePara(hitArmOne1, hitArmTwo2, offset, trackSlope);
    EUTelTrack track = getTrack(hits,offset,trackSlope,curvCorr);
    return track;

}



std::map<int,std::vector<EUTelHit> > EUTelPatRecTriplets::getTrackHitsFromTriplets(std::vector<EUTelPatRecTriplets::triplets>& tripletsVec){
    streamlog_out(DEBUG1) << "Set track states and hits... " << std::endl;
    std::map<int,std::vector<EUTelHit> > fitID_Hits;
    std::vector<triplets>::iterator itTriplet;
    std::vector<triplets> leftTriplets;
    std::vector<triplets> rightTriplets;
    unsigned int fitID = 0;
    for(itTriplet = tripletsVec.begin();itTriplet != tripletsVec.end();  itTriplet++){
      streamlog_out(DEBUG0) << "itTriplet->matches = " <<itTriplet->matches<< std::endl;
        if(itTriplet->cenPlane == 1 ){
            leftTriplets.push_back(*itTriplet);
        }else if(itTriplet->cenPlane == 4){
            rightTriplets.push_back(*itTriplet);
        }else{
            throw(lcio::Exception( "Triplet are not from the left and right arms!"  ));
        }
    }
    std::vector<triplets>::iterator itLeftTriplet;
    std::vector<triplets>::iterator itRightTriplet;
    streamlog_out(DEBUG1) << "Use triplets... " << std::endl;
    for(itLeftTriplet = leftTriplets.begin();itLeftTriplet != leftTriplets.end();  itLeftTriplet++){
        for(itRightTriplet = rightTriplets.begin();itRightTriplet != rightTriplets.end();  itRightTriplet++){
	  streamlog_out(DEBUG0) << "itLeftTriplet->matches = " <<itLeftTriplet->matches<< "  itRightTriplet->matches = " <<itRightTriplet->matches<< std::endl;

                    streamlog_out(DEBUG1) << "Triplet slope delta match Cut: " <<"X delta: " << fabs(itRightTriplet->slope.at(0) - itLeftTriplet->slope.at(0)) <<" Cut: " << _tripletSlopeCuts.at(0) << " Y delta: " <<  fabs(itRightTriplet->slope.at(1) - itLeftTriplet->slope.at(1))<<" Cut: " <<  _tripletSlopeCuts.at(1)  << std::endl;

            if(fabs(itRightTriplet->slope.at(0) - itLeftTriplet->slope.at(0)) > _tripletSlopeCuts.at(0)  or fabs(itRightTriplet->slope.at(1) - itLeftTriplet->slope.at(1)) >_tripletSlopeCuts.at(1)  ){
                continue;
            }
            streamlog_out(DEBUG1) << "PASS 3!! " << std::endl;
            const float aveZPosTrip = (itLeftTriplet->pos.at(2)+ itRightTriplet->pos.at(2))/2.0;
            streamlog_out(DEBUG1) << "Triplet Propagation..." <<std::endl;   
            streamlog_out(DEBUG1) << "LEFT" <<std::endl;   
            std::vector<float> posLeftAtZ = getTripPosAtZ(*itLeftTriplet,aveZPosTrip); 
            streamlog_out(DEBUG1) << "RIGHT" <<std::endl;   
            std::vector<float> posRightAtZ = getTripPosAtZ(*itRightTriplet,aveZPosTrip); 
            streamlog_out(DEBUG1) << "Predicted position of triplet1/triplet2 " << posLeftAtZ.at(0) <<"/"<<posRightAtZ.at(0) << "  " << posLeftAtZ.at(1) <<"/"<<posRightAtZ.at(1)<< "  " << posLeftAtZ.at(2) <<"/"<<posRightAtZ.at(2) <<std::endl;
            streamlog_out(DEBUG1) << "Delta between Triplets X: "<< fabs(posLeftAtZ.at(0)- posRightAtZ.at(0)) <<" Cut X: " << _tripletConnectDistCut.at(0) << " Delta Y: " << fabs(posLeftAtZ.at(1)- posRightAtZ.at(1)) << " Cut Y: " << _tripletConnectDistCut.at(1) << std::endl;
            if(fabs(posLeftAtZ.at(0)- posRightAtZ.at(0)) > _tripletConnectDistCut.at(0) or fabs(posLeftAtZ.at(1)- posRightAtZ.at(1)) > _tripletConnectDistCut.at(1)){
                continue;
            }
            //Pass without DUT
            streamlog_out(DEBUG1) << "PASS 4!! " << std::endl;
            const float dist = (itRightTriplet->states.at(0).getPositionGlobal() - itLeftTriplet->states.at(2).getPositionGlobal()).Mag();
            itLeftTriplet->matches = itLeftTriplet->matches + 1;
            itRightTriplet->matches = itRightTriplet->matches + 1;
            itLeftTriplet->fitID = fitID;
            itRightTriplet->fitID = fitID;
            fitID++;
        }
    }
    for(itLeftTriplet = leftTriplets.begin();itLeftTriplet != leftTriplets.end();  itLeftTriplet++){streamlog_out(DEBUG1) << "1fitID = " <<fitID<< std::endl;
        for(itRightTriplet = rightTriplets.begin();itRightTriplet != rightTriplets.end();  itRightTriplet++){streamlog_out(DEBUG1) << "itLeftTriplet->fitID = " <<itLeftTriplet->fitID<< "   itRightTriplet->fitID = " <<itRightTriplet->fitID<< std::endl;
            if(itLeftTriplet->fitID == itRightTriplet->fitID){streamlog_out(DEBUG1) << "itLeftTriplet->matches = " <<itLeftTriplet->matches<< "  itRightTriplet->matches = " <<itRightTriplet->matches<< std::endl;
                if(itLeftTriplet->matches == 1 and itRightTriplet->matches == 1 ){
                    streamlog_out(DEBUG1) << "FOUND TRACK FROM TRIPLETS!" << std::endl;
                    hits.insert( hits.end(), itLeftTriplet->hits.begin(), itLeftTriplet->hits.end() );
                    hits.insert( hits.end(), itRightTriplet->hits.begin(), itRightTriplet->hits.end() );
                    fitID_Hits[itLeftTriplet->fitID] = hits;
                }
            }
        }
    }
    return fitID_Hits;
}
EUTelTrack EUTelPatRecTriplets::getTrack(triplets tripLeft,triplets tripRight){
    std::vector<EUTelHit> hits;
	streamlog_out(DEBUG1) << "Fill using left arm..."  << "    Number of states in left arm: " << tripLeft.states.size() << std::endl;
    for(unsigned int i = 0 ; i < tripLeft.states.size() ; ++i){
        hits.push_back(tripLeft.states.at(i).getHit());
    }
	streamlog_out(DEBUG1) << "Fill using right arm..." << "    Number of states in right arm: " << tripRight.states.size()<< std::endl;
    for(unsigned int i = 0 ; i < tripRight.states.size() ; ++i){
        hits.push_back(tripRight.states.at(i).getHit());
    }
    streamlog_out(DEBUG1) << "Got hits from triplet." << std::endl;

    EUTelTrack track = getTrackFourHits(hits);
    return track;


}
EUTelTrack EUTelPatRecTriplets::getTrack(std::vector<EUTelHit> hits, std::vector<double> offset, std::vector<double> trackSlope,double qOverPCorr){
    EUTelTrack track;
    track.setQOverP(-1.0/getBeamMomentum() + qOverPCorr);
    //Calculate prediction using properties of hits only. 
    // loop around planes
    // loop around hits -> if hit matched to plane record hit, if not record my intersection
    for(unsigned int  j = 0; j < (_senZOrderToIDWithoutExcPla.size()); ++j){
        unsigned int sensorID = _senZOrderToIDWithoutExcPla.at(j);
        streamlog_out(DEBUG1) << "The Z position " << j << " sensor ID: " << sensorID  <<std::endl;
        bool hitOnPlane=false;
        EUTelState state; //Create a state for each plane included in the fit.
        for(unsigned int i = 0; i < hits.size(); ++i){//Check the list of hits to see if we have one on this plane.
            if(hits.at(i).getLocation()==sensorID){
                hitOnPlane=true;
                state.setLocation(hits.at(i).getLocation());
                state.setHit(hits.at(i));
                double dz1 = hits.at(i).getPositionGlobal()[2] - offset.at(2);
                double dz2 = hits.at(i).getPositionGlobal()[2] - offset.at(3); 
                double posX = offset.at(0) + dz1*trackSlope.at(0) + 0.5*dz1*dz2*qOverPCorr*EUTelNav::_bFac[0];
                double posY = offset.at(1) + dz1*trackSlope.at(1) + 0.5*dz1*dz2*qOverPCorr*EUTelNav::_bFac[1];
                double dz = (dz1 + dz2)/2.0;
                std::vector<double> slope;
                slope.push_back(trackSlope.at(0)+dz*qOverPCorr*EUTelNav::_bFac[0]);
                slope.push_back(trackSlope.at(1)+dz*EUTelNav::_bFac[0]);
                float intersectionPoint[3];
                intersectionPoint[0] = posX;  intersectionPoint[1] = posY; intersectionPoint[2] = hits.at(i).getPositionGlobal()[2];
                //intersection might not be inside a volume. 
                state.setPositionGlobal(intersectionPoint);
                state.setDirFromGloSlope(slope);
                track.setState(state);
            }
        }
        if(hitOnPlane){streamlog_out(DEBUG1) <<"this should already be recorded as there is a hit on sensor "<<sensorID<<std::endl;}
        else {
            //if no hit match - find intersection
            //loop over duts - loop over sensorid
            //   //can i just add the global position to the track ?
            //think so, look at triplet creation
            double z_dut=geo::gGeometry().getOffsetVector(sensorID)[2];
            double dz1 = z_dut - offset.at(2);
            double dz2 = z_dut - offset.at(3); 
            double dz = (dz1+dz2)/2.0;//how to rewrite this line?
            std::vector<double> slope;
            slope.push_back(trackSlope.at(0) + dz*qOverPCorr*EUTelNav::_bFac[0]);
            slope.push_back(trackSlope.at(1) + dz*qOverPCorr*EUTelNav::_bFac[1]);
            state.setDirFromGloSlope(slope);
            //   TVector3 hitPosGlo = hit.getPositionGlobal();
            double posX = offset.at(0) + dz1*trackSlope.at(0) + 0.5*dz1*dz2*qOverPCorr*EUTelNav::_bFac[0];
            double posY = offset.at(1) + dz1*trackSlope.at(1) + 0.5*dz1*dz2*qOverPCorr*EUTelNav::_bFac[1];
            float intersectionPoint[3];
            intersectionPoint[0] = posX;  intersectionPoint[1] = posY; intersectionPoint[2] = z_dut;
            //   //intersection might not be inside a volume. 
            streamlog_out(DEBUG1)<<"intersection point on sensorID "<<sensorID<<" = "<<	intersectionPoint[0]<<", "<<intersectionPoint[1]<<", "<<intersectionPoint[2]<<std::endl;
            //add explicit check that it intersects with sensor?   but i want edge effects?   
            //add arc length thingy
            state.setPositionGlobal(intersectionPoint);
            state.setLocation(sensorID);
            track.setState(state);

        }//else
    }//loop about planes, j iterator
    setArcLengths(track);
//    std::map<const int,double>  mapSensor;
 //   std::map<const int ,double>  mapAir;
//    double rad =	track.getStates().at(0).computeRadLengthsToEnd(mapSensor, mapAir);
 //   if(rad == 0){
 //       throw(std::string("Radiation length is zero for mimosa tracks."));
//    }
    track.print();
//    setRadLengths(track, mapSensor, mapAir, rad);
    return track;
}
void EUTelPatRecTriplets::setArcLengths(EUTelTrack & track){
    for(std::vector<EUTelState>::iterator itState = track.getStates().begin(); itState != --track.getStates().end(); itState++){
        TVector3 gPos1 = itState->getPositionGlobal();
        TVector3 gPos2 = (itState+1)->getPositionGlobal();
        TVector3 diff = gPos2 - gPos1;
        double arc = diff.Mag();
        itState->setArcLengthToNextState(arc);
    }
}

std::vector<EUTelTrack> EUTelPatRecTriplets::getDUTHit(std::vector<EUTelTrack> & tracks){
    std::vector<EUTelTrack> tracksWithDUTHit;
    if(tracks.size() != 0 ){
        std::map<int ,std::vector<EUTelHit>> ::iterator itIDHit;
        for(itIDHit = _mapHitsVecPerPlane.begin(); itIDHit != _mapHitsVecPerPlane.end(); ++itIDHit) {//Should find a better way than this loop
            if(itIDHit->first <= 5){
                streamlog_out(DEBUG1) << "Mimosa hit."  << std::endl;
                continue;
            }else{
                std::vector<EUTelHit> hits = itIDHit->second;
                streamlog_out(DEBUG1) << "Not a Mimosa hit.  hits.size() = "  << hits.size() <<", itIDHit->first = "<<itIDHit->first<<std::endl;
                //Make sure we have hits associated with this plane
                if(hits.size() != 0 ){
                    std::vector<EUTelHit>::iterator itHit;
                    EUTelHit hitBest;
                    streamlog_out(DEBUG0) << "New DUT!"  << std::endl;
                    for ( itHit = hits.begin(); itHit != hits.end(); ++itHit ) {
                        streamlog_out(DEBUG1) << "DUT hit."  << std::endl;
                        EUTelHit& hit = *itHit;
                        //set for each new hit to large value to pass first track.
                        double distBest=10000000;
                        //Find closest track.
                        EUTelTrack bestTrack;
                        std::vector<EUTelTrack>::iterator itTrack;
                        for(itTrack = tracks.begin(); itTrack != tracks.end();itTrack++){//Need to allow more hits to be collected here.
                            //Get track information.
                            std::vector<double> offset;
                            std::vector<double> trackSlope; 
                            EUTelNav::getTrackAvePara(itTrack->getStates().at(0).getHit(), itTrack->getStates().at(5).getHit(), offset, trackSlope);
                            //find state information at hit location.
                            double dzToHit = hit.getPositionGlobal()[2] - 0.5*( offset.at(2)+ offset.at(3));
                            std::vector<float> slope;
                            //TO DO: We use the old curv to find the tracks but the corrected in the parameterisation. Should use corrected here too.
                            slope.push_back(trackSlope.at(0) + dzToHit*getCurvXY()[0]);
                            slope.push_back(trackSlope.at(1) + dzToHit*getCurvXY()[1]);
                            double dz1 = hit.getPositionGlobal()[2] - offset.at(2);
                            double dz2 = hit.getPositionGlobal()[2] - offset.at(3); 
                            double posX = offset.at(0) + dz1*trackSlope.at(0) + 0.5*dz1*dz2*getCurvXY()[0];
                            double posY = offset.at(1) + dz1*trackSlope.at(1) + 0.5*dz1*dz2*getCurvXY()[1];
                            double dist=1000000;
                            /// Transform prediction to local frame of DUT. 
                            /// This is done so the prediction and measured hit is compared in the correct axis.
                            /// Since strip sensors will only have information about position in the local DUT x-axis by default. 
                            double locPos [3];
                            TVector3 hitPosGlo = hit.getPositionGlobal();
                            const double referencePoint[]	= {posX,posY ,hit.getPositionGlobal()[2]};
                            geo::gGeometry().master2Local( hit.getLocation(), referencePoint, locPos);
                            if(_planeDimensions[hit.getLocation()] == 2){ 
                                streamlog_out(DEBUG0) << "Triplet DUT Match Cut Pixel: " <<"X delta: " << fabs(posX-hitPosGlo[0]) << " Y delta: " << fabs(posY - hitPosGlo[1]) << std::endl;
                                double dist = sqrt(pow(locPos[0]-hit.getPosition()[0],2)+pow(locPos[1]-hit.getPosition()[1],2));
                            }else if(_planeDimensions[hit.getLocation()] == 1){
                                streamlog_out(DEBUG0) << "Triplet DUT Match Cut Strip: " <<"X delta: " << fabs(posX-hitPosGlo[0]) << std::endl;
                                ///Keep it positive, the distance that is!!!
                                dist = sqrt(pow(locPos[0]-hit.getPosition()[0],2));

                            }else{
                                throw(lcio::Exception( "This is not a strip or pixel sensor!"));
                            }
                            //Will enter always on the first loop if reached. Need to find intersection and hit on DUT.
                            if(dist < distBest){
                                streamlog_out(DEBUG0) << "Save track information"  << std::endl;
                                if(itTrack == tracks.begin()){
                                    streamlog_out(DEBUG0) << "First track used!" << std::endl;
                                }
                                if(dist < distBest){
                                    streamlog_out(DEBUG0) << "Improvement! Use this track." << std::endl;
                                }
                                bestTrack = *itTrack;
                            }
                        }
                        //Create new track from dutHit and mimosa track.
                        std::vector<EUTelHit> dutHits;
                        dutHits.push_back(*itHit);
                        std::vector<EUTelHit> hits;
                        hits = getDUTHitsOrder(bestTrack,dutHits);
                        tracksWithDUTHit.push_back(getTrackFourHits(hits));
                    }
                }else{//If there are no hit on DUT plane.
                    _tracksWithoutHit = _tracksWithoutHit + tracks.size();
                }
            }
        }
    }
    return tracksWithDUTHit;
}
std::vector<EUTelHit> EUTelPatRecTriplets::getDUTHitsOrder(EUTelTrack track, std::vector<EUTelHit> dutHits ){
    std::vector<EUTelHit> finalHits;
    for(unsigned int  i = 0; i < (_senZOrderToIDWithoutExcPla.size()); ++i){
        unsigned int sensorID = _senZOrderToIDWithoutExcPla.at(i);
        streamlog_out(DEBUG0) << "The Z position " << i << " sensor ID: " << sensorID  <<std::endl;
        std::vector<EUTelState>::iterator itState;
        for( itState = track.getStates().begin(); itState != track.getStates().end(); ++itState){
            if(itState->getLocation() == sensorID){
                streamlog_out(DEBUG0) << "Add mimosa hit " <<std::endl; 
                itState->print();
                if(itState->getStateHasHit()){
                    finalHits.push_back(itState->getHit());
                }
                streamlog_out(DEBUG0) << "Added! " <<std::endl; 
                break;
            }

        }
        dutHits.at(0).print();
        std::vector<EUTelHit>::iterator itDUTHit;
        for(itDUTHit = dutHits.begin(); itDUTHit != dutHits.end(); ++itDUTHit){
            itDUTHit->print();
            if(itDUTHit->getLocation() == sensorID){
                streamlog_out(DEBUG0) << "Add DUT hit " <<std::endl; 
                finalHits.push_back(*itDUTHit);
                streamlog_out(DEBUG0) << "Added! " <<std::endl; 
                break;
            }
        }
    }
    return finalHits;
}
std::vector<EUTelHit> EUTelPatRecTriplets::getCorrHitOrder(std::vector<EUTelHit> hits ){
    std::vector<EUTelHit> finalHits;
    for(unsigned int  i = 0; i < (_senZOrderToIDWithoutExcPla.size()); ++i){
        unsigned int sensorID = _senZOrderToIDWithoutExcPla.at(i);
        std::vector<EUTelHit>::iterator itHit;
        for(itHit = hits.begin(); itHit != hits.end(); ++itHit){
            itHit->print();
            if(itHit->getLocation() == sensorID){
                streamlog_out(DEBUG0) << "Add hit " <<std::endl; 
                finalHits.push_back(*itHit);
                streamlog_out(DEBUG0) << "Added! " <<std::endl; 
                break;
            }
        }
    }
    return finalHits;
}



std::vector<float>  EUTelPatRecTriplets::getTripPosAtZ(triplets trip, float posZ ){
	streamlog_out(DEBUG1) << "Slope x/y: " <<trip.slope.at(0) << "  " <<trip.slope.at(1) <<" Position ave: " <<trip.pos.at(0)<<" "<<trip.pos.at(1)<<" "<<trip.pos.at(2) <<std::endl;
    float dz = posZ - trip.pos.at(2);
    float x = trip.pos.at(0) + trip.slope.at(0)*dz;
    float y = trip.pos.at(1) + trip.slope.at(1)*dz;
    std::vector<float> position;
    position.push_back(x);position.push_back(y);position.push_back(posZ);
    return position;
}
bool EUTelPatRecTriplets::getDoubHitOnTraj(const doublets& doub, const std::vector<unsigned int> & sen,std::vector<EUTelHit>& newHits   ){
    for(std::vector::iterator<unsigned int> itID = sen.begin(); itID != sen.end(); ++itID){
        std::vector<EUTelHit> hits =  _mapHitsVecPerPlane.at(*itID);
        float distBest = 10000000;
        EUTelHit hitBest;
        for(std::vector::iterator<EUTelHit> itHit = hits.begin(); itHit != hits.end(); ++itHit){
            double hitPosX = itHit->getPositionGlobal()[0]; 
            double hitPosY = itHit->getPositionGlobal()[1]; 
            double hitPosZ = itHit->getPositionGlobal()[2]; 
            std::vector<float>  pos = getDoubPosAtZ(doub, hitPosZ);/// Could calculate this once. Might be a bit off for tilted sensors.
            double locPos [3];
            TVector3 hitPosGlo = hit.getPositionGlobal();
            float dist = getDistLocal(hit.getLocation(), pos, hitPosGlo);
            if(itHit == hits.begin()){
                hitBest = *itHit;
                distBest = dist;
            }
            if(dist < distBest){
                hitBest = *itHit;
                distBest = dist;
            }
        }
        if(distBest >  _doubletCenDistCut.at(0)){
            continue;
        }
        streamlog_out(DEBUG1) << "PASS Doublet cut!! " << std::endl;
        newHits.push_back(hitBest);
    }
    ///Number of hits cut here.
    return true;
}
float EUTelPatRecTriplets::getDistLocal(int& location, std::vector<float>& pos, TVector3& hitPosGlo){
    const double referencePoint[]	= {pos.at(0),pos.at(1) ,hitPosZ};
    geo::gGeometry().master2Local( hit.getLocation(), referencePoint, locPos);

    if(_planeDimensions[hit.getLocation()] == 2){ 
        streamlog_out(DEBUG0) << "Triplet DUT Match Cut Pixel: " <<"X delta: " << fabs(pos.at(0)-hitPosGlo[0]) << " Y delta: " << fabs(pos.at(1) - hitPosGlo[1]) << std::endl;
        double dist = sqrt(pow(locPos[0]-itHit.getPosition()[0],2)+pow(locPos[1]-itHit.getPosition()[1],2));
    }else if(_planeDimensions[hit.getLocation()] == 1){
        streamlog_out(DEBUG0) << "Triplet DUT Match Cut Strip: " <<"X delta: " << fabs(pos.at(0)-hitPosGlo[0]) << std::endl;
        ///Keep it positive, the distance that is!!!
        dist = sqrt(pow(locPos[0]-hit.getPosition()[0],2));
    }else{
        throw(lcio::Exception( "This is not a strip or pixel sensor!"));
    }
    return dist;
}


std::vector<float>  EUTelPatRecTriplets::getDoubPosAtZ(doublets doub, float posZ){
    float dz = posZ - doub.pos.at(2);
    float x = doub.pos.at(0) + doub.slope.at(0)*dz;
    float y = doub.pos.at(1) + doub.slope.at(1)*dz;
    std::vector<float> position;
    position.push_back(x);position.push_back(y);position.push_back(posZ);
    return position;
}

bool EUTelPatRecTriplets::getDoublet(EUTelHit & hit1, EUTelHit & hit2, doublets& doublet){
    const double curvX = EUTelNav::_curv.at(0); 
    const double curvY =  EUTelNav::_curv.at(1); 

    float initDis = geo::gGeometry().getInitialDisplacementToFirstPlane();
    //Remove the curvature as a factor between hits. Therefore only the slope will displace the hit position from plane to plane.
    float x1 = hit1.getPositionGlobal()[0] - 0.5*curvX*pow(hit1.getPositionGlobal()[2] - initDis, 2);
    float y1 = hit1.getPositionGlobal()[1] - 0.5*curvY*pow(hit1.getPositionGlobal()[2] - initDis, 2);
    float x2 = hit2.getPositionGlobal()[0] - 0.5*curvX*pow(hit2.getPositionGlobal()[2] - initDis, 2);
    float y2 = hit2.getPositionGlobal()[1] - 0.5*curvY*pow(hit2.getPositionGlobal()[2] - initDis, 2);

    doublets doublet;
    doublet.pos.push_back((x2 + x1)/2.0);
    doublet.pos.push_back((y2 + y1)/2.0);
    doublet.pos.push_back((hit2.getPositionGlobal()[2] + hit1.getPositionGlobal()[2])/2.0);

    doublet.diff.push_back(x2 - x1);
    doublet.diff.push_back(y2 - y1);
    doublet.diff.push_back( hit2.getPositionGlobal()[2] - hit1.getPositionGlobal()[2]);

    streamlog_out(DEBUG1) << "Doublet delta X: "<< std::abs(doublet.diff.at(0)) << " Cut X: "<<_doubletDistCut.at(0) <<" Delta Y: " << doublet.diff.at(1)<< " Cut Y: " << _doubletDistCut.at(1)  << std::endl;

    if(fabs(doublet.diff.at(0)) >  _doubletDistCut.at(0) or fabs(doublet.diff.at(1)) >  _doubletDistCut.at(1) ){
        return false;
    }

    doublet.slope.push_back( doublet.diff.at(0)/doublet.diff.at(2)); 
    doublet.slope.push_back( doublet.diff.at(1)/doublet.diff.at(2)); 
    doublet.push_back(hit1);
    doublet.push_back(hit2);
    return true;
}

/*This creates map between plane ID and hits on that plane.  We also order the map correcly with geometry.
 */
void EUTelPatRecTriplets::setHitsVecPerPlane()
{
	_mapHitsVecPerPlane.clear();
	int numberOfPlanes = _senZOrderToIDWithoutExcPla.size();
	
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
        std::vector<EUTelHit> tempHitsVecPlaneX; 
		for(size_t j=0 ; j<_allHitsVec.size();++j)
		{
			if(Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*>(_allHitsVec[j]) ) ==  _senZOrderToIDWithoutExcPla.at(i))
			{
				tempHitsVecPlaneX.push_back(EUTelHit(_allHitsVec.at(j)));
			}
		}
		_mapHitsVecPerPlane[ _senZOrderToIDWithoutExcPla.at(i)] = 	tempHitsVecPlaneX;
	}	
}

void EUTelPatRecTriplets::testHitsVecPerPlane(){
	if(_mapHitsVecPerPlane.size() !=  _senZOrderToIDWithoutExcPla.size()){
		streamlog_out(ERROR0) <<"The size of the planes with hits " << _mapHitsVecPerPlane.size() <<"  Sensors from Geometry with no excluded planes: "<<  _senZOrderToIDWithoutExcPla.size()<<std::endl;
		throw(lcio::Exception("The number of planes that could contain hits and the number of planes is different. Problem could be the gear file has to planes at the same z position.")); 	

	}
	for(size_t i=0 ;i<_mapHitsVecPerPlane.size();++i){
		int planeID = _senZOrderToIDWithoutExcPla.at(i);
		if(_mapHitsVecPerPlane.at(planeID).size() <= 0){
			streamlog_out(DEBUG1) << "One plane has no hits at all. Is this correct?" << std::endl;
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
