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
    _tripletsVec.clear();
    float omega = -1.0/_beamE;
	const gear::BField& Bfield = geo::gGeometry().getMagneticField();
	gear::Vector3D vectorGlobal(0.1,0.1,0.1);
	const double Bx = (Bfield.at( vectorGlobal ).x());  
	const double By = (Bfield.at( vectorGlobal ).y());
	const double Bz = (Bfield.at( vectorGlobal ).z());
    double curvX = 0.0003*Bx*omega; 
    double curvY = 0.0003*By*omega; 
    std::vector<unsigned int > cenPlane(1,4);

    for(std::vector<unsigned int>::iterator cenPlaneID = cenPlane.begin(); cenPlaneID != cenPlane.end(); ++cenPlaneID) {
        unsigned int cenID = *cenPlaneID; 
        EVENT::TrackerHitVec& hitCentre = _mapHitsVecPerPlane[cenID];
        EVENT::TrackerHitVec& hitCentreLeft = _mapHitsVecPerPlane[cenID - 1];
        EVENT::TrackerHitVec& hitCentreRight = _mapHitsVecPerPlane[cenID + 1];

		EVENT::TrackerHitVec::iterator itHit;
		EVENT::TrackerHitVec::iterator itHitLeft;
		EVENT::TrackerHitVec::iterator itHitRight;
		for ( itHitLeft = hitCentreLeft.begin(); itHitLeft != hitCentreLeft.end(); ++itHitLeft ) {
            const int hitLeftLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHitLeft) );
            double hitLeftPos[] = { (*itHitLeft)->getPosition()[0], (*itHitLeft)->getPosition()[1], (*itHitLeft)->getPosition()[2] };
            double hitLeftPosGlobal[3];
            geo::gGeometry().local2Master(hitLeftLoc ,hitLeftPos,hitLeftPosGlobal);
            for ( itHitRight = hitCentreRight.begin(); itHitRight != hitCentreRight.end(); ++itHitRight ) {
                const int hitRightLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHitRight) );
                double hitRightPos[] = { (*itHitRight)->getPosition()[0], (*itHitRight)->getPosition()[1], (*itHitRight)->getPosition()[2] };
                double hitRightPosGlobal[3];
                geo::gGeometry().local2Master(hitRightLoc ,hitRightPos,hitRightPosGlobal);
                doublets doublet;
                doublet = getDoublet(hitLeftPosGlobal,hitRightPosGlobal,curvX,curvY);
                if(abs(doublet.diff.at(0)) >  _doubletDistCut.at(0) or abs(doublet.diff.at(1)) >  _doubletDistCut.at(1) ){
                    continue;
                }
                //Now loop through all hits on plane between two hits which create doublets. 
                for ( itHit = hitCentre.begin(); itHit != hitCentre.end(); ++itHit ) {
                    const int hitLoc = Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*> (*itHit) );
                    double hitPos[] = { (*itHit)->getPosition()[0], (*itHit)->getPosition()[1], (*itHit)->getPosition()[2] };
                    double hitPosGlobal[3];
                    geo::gGeometry().local2Master(hitLoc ,hitPos,hitPosGlobal);
                    float initDis = geo::gGeometry().getInitialDisplacementToFirstPlane();
                    float x1 = hitPos[0] - 0.5*curvX*pow(hitPos[2] - initDis, 2);
                    float y1 = hitPos[1] - 0.5*curvY*pow(hitPos[2] - initDis, 2);
                    double delX = doublet.pos.at(0) - x1;
                    double delY = doublet.pos.at(1) - y1;
                    if(abs(delX) >  _doubletCenDistCut.at(0) or abs(delY) >  _doubletCenDistCut.at(1) ){
                        continue;
                    }
                    triplets triplet;
                    triplet.hits.push_back(EUTelHit(*itHitLeft));
                    triplet.hits.push_back(EUTelHit(*itHit));
                    triplet.hits.push_back(EUTelHit(*itHitRight));
                    triplet.pos.push_back(doublet.pos.at(0));
                    triplet.pos.push_back(doublet.pos.at(1));
                    triplet.pos.push_back(doublet.pos.at(2));
                    triplet.cenPlane = cenID;
                    triplet.slope.push_back(doublet.slope.at(0));
                    triplet.slope.push_back(doublet.slope.at(1));
                    _tripletsVec.push_back(triplet); 
                }
            }
        }
    }
}
std::vector<EUTelTrack> EUTelPatRecTriplets::getTracks( ){
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
    for(itLeftTriplet = leftTriplets.begin();itLeftTriplet != leftTriplets.end();  itLeftTriplet++){
        for(itRightTriplet = rightTriplets.begin();itRightTriplet != rightTriplets.end();  itRightTriplet++){
            if(abs(itRightTriplet->slope.at(0) - itLeftTriplet->slope.at(0)) > _tripletSlopeCuts.at(0)  or abs(itRightTriplet->slope.at(1) - itLeftTriplet->slope.at(1)) >_tripletSlopeCuts.at(1)  ){
                continue;
            }
            //Do we have a DUT. We will only add one DUT at the moment this must be updated to add multiple DUTs 
            if(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() - 6 != 0 ){
            }else{
                //No DUT use average 
                float aveZPosTrip = (itLeftTriplet->pos.at(2)+ itRightTriplet->pos.at(2))/2.0;
                std::vector<float> posLeftAtZ = getTripPosAtZ(*itLeftTriplet,aveZPosTrip); 
                std::vector<float> posRightAtZ = getTripPosAtZ(*itRightTriplet,aveZPosTrip); 
                if(abs(posLeftAtZ.at(0)- posRightAtZ.at(0)) > _tripletConnectDistCut.at(0) or abs(posLeftAtZ.at(0)- posRightAtZ.at(0)) > _tripletConnectDistCut.at(0)){
                    continue;
                }

            }
        }
    }

}
std::vector<float>  EUTelPatRecTriplets::getTripPosAtZ(triplets trip, float posZ ){
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
    doublet.diff.push_back( doublet.diff.at(1)/doublet.diff.at(2)); 

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
