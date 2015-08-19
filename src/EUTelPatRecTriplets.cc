#include "EUTelPatRecTriplets.h"
#include <IMPL/TrackerPulseImpl.h>
#include "EUTelNav.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
namespace eutelescope {

  EUTelPatRecTriplets::EUTelPatRecTriplets()
  {}

EUTelPatRecTriplets::EUTelPatRecTriplets(AIDA::IHistogram1D * DoubletXseperationHistoRight, AIDA::IHistogram1D * DoubletYseperationHistoRight, AIDA::IHistogram1D * DoubletXseperationHistoLeft,
					   AIDA::IHistogram1D * DoubletYseperationHistoLeft, AIDA::IHistogram1D * TripletXseperationHistoRight, AIDA::IHistogram1D * TripletYseperationHistoRight,
					   AIDA::IHistogram1D * TripletXseperationHistoLeft, AIDA::IHistogram1D * TripletYseperationHistoLeft, AIDA::IHistogram1D * TripletDistCutXHisto,
					   AIDA::IHistogram1D * TripletDistCutYHisto, AIDA::IHistogram1D * TripletSlopeHistoX, AIDA::IHistogram1D * TripletSlopeHistoY, AIDA::IHistogram1D * DUTWindowHisto ):  
_eventNumber(0),
_totalNumberOfHits(0),
_numberTripletsLeft(0),
_numberTripletsRight(0),
_tripletSlopeCuts(0,0),
_beamE(-1.),
_hitNum(6),
_mode(1),
_numberOfTracksTotal(0),
_numberOfTracksDUTTotal(0),
_DoubletXseperationHistoRight(DoubletXseperationHistoRight),
_DoubletYseperationHistoRight(DoubletYseperationHistoRight),
_DoubletXseperationHistoLeft(DoubletXseperationHistoLeft),
_DoubletYseperationHistoLeft(DoubletYseperationHistoLeft),
_TripletXseperationHistoRight(TripletXseperationHistoRight),
_TripletYseperationHistoRight(TripletYseperationHistoRight),
_TripletXseperationHistoLeft(TripletXseperationHistoLeft),
_TripletYseperationHistoLeft(TripletYseperationHistoLeft),
_TripletDistCutXHisto(TripletDistCutXHisto),
_TripletDistCutYHisto(TripletDistCutYHisto),
_TripletSlopeHistoX(TripletSlopeHistoX),
_TripletSlopeHistoY(TripletSlopeHistoY),
_DUTWindowHisto(DUTWindowHisto)
{
    EUTelExcludedPlanes();

}
EUTelPatRecTriplets::~EUTelPatRecTriplets()  
{}

void EUTelPatRecTriplets::setPlaneDimensionsVec(EVENT::IntVec& planeDimensions){
	if(planeDimensions.size() != EUTelExcludedPlanes::_senNoDeadMaterial.size() and std::find(EUTelExcludedPlanes::_senNoDeadMaterial.begin(), EUTelExcludedPlanes::_senNoDeadMaterial.end(), 271) == EUTelExcludedPlanes::_senNoDeadMaterial.end() ){
		streamlog_out(ERROR) << "The size of planesDimensions input is: "<< planeDimensions.size()<<" The size of sensorIDsVec is: " << EUTelExcludedPlanes::_senNoDeadMaterial.size()<< std::endl;
		throw(lcio::Exception( "The input dimension vector not the same as the number of planes!"));
	}
	_planeDimensions.clear();
    unsigned int scatCounter=0;
	for(size_t i=0; i<EUTelExcludedPlanes::_senNoDeadMaterial.size(); ++i){
		const int planeID = EUTelExcludedPlanes::_senNoDeadMaterial.at(i);
		if (_planeDimensions.find(planeID) == _planeDimensions.end() and planeID != 271){//This is to check that we don't try to map the same sensor to two different plane dimensions.
			_planeDimensions[planeID] = planeDimensions.at(i-scatCounter);
		}else{
//			streamlog_out(ERROR5) <<"The z position is : "<< i <<" The plane ID is: " << planeID <<std::endl;
//			throw(lcio::Exception( "You are trying to map the same sensor ID to two different plane dimensions. There is something wrong with you gear file input. Make sure there is some distance between your planes in the gear file!"));
		}
        if(planeID == 271){
            _planeDimensions[planeID] = 2;
            scatCounter++;
        }
	}
}	    
///Public get functions
/*This is the function which should be used in processEvent
*/
std::vector<EUTelTrack> EUTelPatRecTriplets::getTracks(){
    EUTelNav::init(getBeamMomentum());
    setHitsVecPerPlane();
    testHitsVecPerPlane();
    std::vector<EUTelTrack>  tracks;
    if(_mode = 1){
        tracks = getMinFakeTracks();
    }else{
        ///Basic pattern recognition used to find DUT hit.
        /// Can also use basic pattern recognition to find more tracks than triplets.
    }
    for(size_t i = 0 ; i < tracks.size();++i){
        tracks.at(i).print();
    }
    streamlog_out( DEBUG1 ) << "tracks.size() = "<<tracks.size()<<std::endl;
    printTrackQuality(tracks );
    return tracks;
}

//get functions private
EUTelTrack EUTelPatRecTriplets::getTrack(triplets tripLeft,triplets tripRight){
    std::vector<EUTelHit> hits;
	streamlog_out(DEBUG1) << "Fill using left arm..."  << "    Number of states in left arm: " << tripLeft.hits.size() << std::endl;
    for(size_t i = 0 ; i < tripLeft.hits.size() ; ++i){
        hits.push_back(tripLeft.hits.at(i));
    }
	streamlog_out(DEBUG1) << "Fill using right arm..." << "    Number of states in right arm: " << tripRight.hits.size()<< std::endl;
    for(size_t i = 0 ; i < tripRight.hits.size() ; ++i){
        hits.push_back(tripRight.hits.at(i));
    }
    streamlog_out(DEBUG1) << "Got hits from triplet." << std::endl;

    EUTelTrack track = EUTelTrackCreate::getTrackFourHits(hits);
    return track;



	if(_numberOfTracksTotal % 500 == 0){
	  //streamlog_out(MESSAGE5) << "Percentage tracks without DUT hit: " << static_cast<float>(_tracksWithoutHit)/static_cast<float>(_numberOfTracksTotal)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of tracks per event: " << static_cast<float>(_numberOfTracksTotal)/static_cast<float>(getEventNumber() +1)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of left arm triplets per event: " << static_cast<float>(_numberTripletsLeft)/static_cast<float>(getEventNumber() +1)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of right arm triplets per event: " << static_cast<float>(_numberTripletsRight)/static_cast<float>(getEventNumber() +1)<< std::endl;
    }

}
bool EUTelPatRecTriplets::getTriplet(doublets& doublet, EUTelHit& hitCen,triplets& triplet){

    const float initDis = geo::gGeometry().getInitialDisplacementToFirstPlane();
    const float x1 = hitCen.getPositionGlobal()[0] - 0.5*EUTelNav::_curv[0]*pow(hitCen.getPositionGlobal()[2] - initDis, 2);
    const float y1 = hitCen.getPositionGlobal()[1] - 0.5*EUTelNav::_curv[1]*pow(hitCen.getPositionGlobal()[2] - initDis, 2);
    std::vector<float> pos = getDoubPosAtZ(doublet, hitCen.getPositionGlobal()[2]);
    const double delX = pos.at(0) - x1;
    const double delY = pos.at(1) - y1;
    streamlog_out(DEBUG1) << "Doublet delta X to centre hit: "<< delX << " Cut X: "<<_doubletCenDistCut.at(0) <<" Distance Y: " << delY<< " Cut Y: " << _doubletCenDistCut.at(1)  << std::endl;


		if(hitCen.getLocation() == 1){
		  _TripletXseperationHistoLeft->fill(delX);
		  _TripletYseperationHistoLeft->fill(delY);
		}

		else if(hitCen.getLocation() == 4){
		  _TripletXseperationHistoRight->fill(delX);
		  _TripletYseperationHistoRight->fill(delY);
		}

    if(fabs(delX) >  _doubletCenDistCut.at(0) or fabs(delY) >  _doubletCenDistCut.at(1) ){
        return false;
    }


    streamlog_out(DEBUG1) << "PASS 2!! " << std::endl;

    triplet.matches = 0;
    triplet.fitID = -999;
    triplet.pos.push_back(doublet.pos.at(0));
    triplet.pos.push_back(doublet.pos.at(1));
    triplet.pos.push_back(doublet.pos.at(2));
    ///Used to find final tracks.
    triplet.cenPlane = hitCen.getLocation();
    triplet.slope.push_back(doublet.slope.at(0));
    triplet.slope.push_back(doublet.slope.at(1));
    ///Initialised vector problems here!!
    std::vector<EUTelHit>  hits;
    hits.push_back(doublet.hits.at(0));
    hits.push_back(hitCen);
    hits.push_back(doublet.hits.at(1));
    triplet.hits = hits;
    streamlog_out(DEBUG1) << "triplet created." << std::endl;
    return true;

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
	    getDoublet(*itHitLeft,*itHitRight,doublet);
	    
	    if (i ==0){
	      _DoubletXseperationHistoLeft->fill(doublet.diff.at(0));
	      _DoubletYseperationHistoLeft->fill(doublet.diff.at(1));
	    }
	    
	    if(i==1){
	      _DoubletXseperationHistoRight->fill(doublet.diff.at(0));
	      _DoubletYseperationHistoRight->fill(doublet.diff.at(1));
	    }
	    
	    if(fabs(doublet.diff.at(0)) >  _doubletDistCut.at(0) or fabs(doublet.diff.at(1)) >  _doubletDistCut.at(1) ){
	      continue;
	    }
	    streamlog_out(DEBUG1) << "PASS 1!! " << std::endl;
	    
	    //Now loop through all hits on plane between two hits which create doublets. 
	    for ( itHit = hitCentre.begin(); itHit != hitCentre.end(); ++itHit ) {
	      triplets triplet;
	      bool pass2 = getTriplet(doublet,*itHit, triplet);
	      if(!pass2) continue;
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
        if(cenID.size()-1 == i and tripletsVec.size() == 0){//If not triplets found on first plane end search.
            streamlog_out(DEBUG1) << "FOUND NO TRIPLETS FOR EVENT: " << getEventNumber()  << std::endl;
            break;
        }
    }
    return tripletsVec;

}
std::vector<std::vector<EUTelHit> > EUTelPatRecTriplets::getTrackHitsFromTriplets(std::vector<EUTelPatRecTriplets::triplets>& tripletsVec){
    streamlog_out(DEBUG1) << "Set track states and hits... " << std::endl;
    std::vector<std::vector<EUTelHit> >  tracksHits;
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
           _TripletSlopeHistoX ->fill(itRightTriplet->slope.at(0) - itLeftTriplet->slope.at(0));
           _TripletSlopeHistoY ->fill(itRightTriplet->slope.at(1) - itLeftTriplet->slope.at(1));

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

	    _TripletDistCutXHisto->fill(posLeftAtZ.at(0)- posRightAtZ.at(0));
	    _TripletDistCutYHisto->fill(posLeftAtZ.at(1)- posRightAtZ.at(1));

            if(fabs(posLeftAtZ.at(0)- posRightAtZ.at(0)) > _tripletConnectDistCut.at(0) or fabs(posLeftAtZ.at(1)- posRightAtZ.at(1)) > _tripletConnectDistCut.at(1)){
                continue;
            }


            //Pass without DUT
            streamlog_out(DEBUG1) << "PASS 4!! " << std::endl;
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
                    std::vector<EUTelHit> hits;
                    hits.insert( hits.end(), itLeftTriplet->hits.begin(), itLeftTriplet->hits.end() );
                    hits.insert( hits.end(), itRightTriplet->hits.begin(), itRightTriplet->hits.end() );
                    tracksHits.push_back(hits);
                }
            }
        }
    }
    return tracksHits;
}
std::vector<EUTelHit> EUTelPatRecTriplets::getCorrHitOrder(std::vector<EUTelHit> hits ){
  streamlog_out(DEBUG0) << "enters getCorrHitOrder" <<std::endl; 
    std::vector<EUTelHit> finalHits;
    for(unsigned int  i = 0; i < (EUTelExcludedPlanes::_senInc.size()); ++i){
        unsigned int sensorID = EUTelExcludedPlanes::_senInc.at(i);
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


bool EUTelPatRecTriplets::getDoubHitOnTraj(doublets const& doub, std::vector<unsigned int> const & sen,int const & hitNum, std::vector<EUTelHit>& newHits   ){
    for(std::vector<unsigned int> ::const_iterator itID = sen.begin(); itID != sen.end(); ++itID){
        std::vector<EUTelHit> hits =  _mapHitsVecPerPlane.at(*itID);
        float distBest = 10000000;
        EUTelHit hitBest;
//        std::cout<<"Hits Number " << hits.size() << " " << " ID " << *itID  << std::endl;

        for(std::vector<EUTelHit>::iterator itHit = hits.begin(); itHit != hits.end(); ++itHit){
  //          std::cout<<"Hit location " << itHit->getLocation() << std::endl;

            double hitPosX = itHit->getPositionGlobal()[0]; 
            double hitPosY = itHit->getPositionGlobal()[1]; 
            double hitPosZ = itHit->getPositionGlobal()[2]; 

            std::vector<float>  pos = getDoubPosAtZ(doub, hitPosZ);/// Could calculate this once. Might be a bit off for tilted sensors.

            float dist = getDistLocal(itHit, pos);

    //        std::cout<<"Dist " << dist  << std::endl;

            if(itHit == hits.begin()){
                hitBest = *itHit;
                distBest = dist;
      //          std::cout<<"DistBest begin " << distBest  << std::endl;

            }

            if(dist < distBest){
                hitBest = *itHit;
                distBest = dist;
        //        std::cout<<"DistBest " << distBest  << std::endl;
            }
            if(itHit == hits.end()-1){
               _DUTWindowHisto ->fill(distBest);
            }
        }

        if(distBest >  _dutDistCut){
   //         std::cout<<"DistBest Fail!!!!! " << distBest  << std::endl;

            streamlog_out(DEBUG1) << "Doublet cut!! " << distBest <<">"<< _dutDistCut<<std::endl;
            continue;
        }
    //    std::cout<<"Pass "  << std::endl;

        streamlog_out(DEBUG1) << "PASS Doublet cut!! " << distBest <<"<"<< _dutDistCut << std::endl;
        newHits.push_back(hitBest);
    }
    if(newHits.size() < hitNum){
        return false;
    //    std::cout<<"Fail number of hits "  << std::endl;

    }else{
    //    std::cout<<"Pass hit number "  << std::endl;
        return true;
    }
}
float EUTelPatRecTriplets::getDistLocal(std::vector<EUTelHit>::iterator itHit, std::vector<float>& pos){
    double locPos [3];
    const double referencePoint[]	= {pos.at(0),pos.at(1) ,itHit->getPositionGlobal()[2]};
    geo::gGeometry().master2Local( itHit->getLocation(), referencePoint, locPos);
    double dist = 99999;
    Eigen::Matrix3d rot  =  geo::gGeometry().getRotMatrixEig(itHit->getLocation());
    Eigen::Matrix3d rotInv = rot.transpose();
    Eigen::Vector3d preG;  
    Eigen::Vector3d measG; 
    Eigen::Vector3d diff; 
    preG  << pos.at(0), pos.at(1), 0;
    measG  << itHit->getPositionGlobal()[0] ,  itHit->getPositionGlobal()[1], 0;
    diff = preG - measG;
    Eigen::Vector3d local;
    local = rotInv*diff;

    if(_planeDimensions[itHit->getLocation()] == 2){ 
        streamlog_out(DEBUG0) <<"Pixel:  X delta: " << fabs(local[0]) << " Y delta: " << fabs(local[1]) << std::endl;
        dist = sqrt(pow(local[0],2)+pow(local[1],2));
    }else if(_planeDimensions[itHit->getLocation()] == 1){
        streamlog_out(DEBUG0) << "Strip: " <<"X delta: " << fabs(local[0]) << std::endl;
        ///Keep it positive, the distance that is!!!
        dist = sqrt(pow(local[0],2));
    }else{
        throw(lcio::Exception( "This is not a strip or pixel sensor!"));
    }
    return dist;
}




std::vector<EUTelTrack> EUTelPatRecTriplets::getMinFakeTracks(){
    std::vector<EUTelTrack> tracks;
    std::vector<EUTelPatRecTriplets::triplets> tripletVec = getTriplets();
    streamlog_out(DEBUG1) << "Total number of triplets found:  " << tripletVec.size()  << std::endl;
    std::vector<std::vector<EUTelHit> >  tracksHits= getTrackHitsFromTriplets(tripletVec);
    ///Loop over all hits which make up a track. 
    std::vector< std::pair< std::vector<EUTelHit> , std::vector<EUTelHit> > > tracksAndDUTHits;
    for(std::vector<std::vector<EUTelHit> >::iterator itTrack = tracksHits.begin(); itTrack != tracksHits.end();++itTrack){
        std::vector<EUTelHit> newHits;
        if(EUTelExcludedPlanes::_senInc.size() > 6 ){///Only look for DUTs if we have more planes included.
    //        std::cout<<"The # tracks: " << tracksHits.size() << std::endl;
            doublets doub;
            getDoublet(*(itTrack->begin()),*(itTrack->rbegin()),doub);
            std::vector< unsigned int> dut;
            for(size_t i = 0 ; i < EUTelExcludedPlanes::_senInc.size(); ++i){
                if(EUTelExcludedPlanes::_senInc.at(i) > 5){
                    dut.push_back(EUTelExcludedPlanes::_senInc.at(i));
                }
            }
            streamlog_out(DEBUG1) << "Got hit! "  << std::endl;
            int hitNum=0; //Need a minimum of 1 DUT hit to pass track.
            bool pass =  getDoubHitOnTraj(doub, dut,hitNum, newHits);
      //      if(!pass){
      //          continue;
       //     }
            tracksAndDUTHits.push_back(make_pair(*itTrack,newHits));
        }
    }
    if(EUTelExcludedPlanes::_senInc.size() > 6 ){
        // Each track can only be associated to a single DUT hit on each plane. Must also make sure that a hit is not associated to multiple tracks.
        // So you could have two tracks with the same DUT hit attached.
        std::vector< std::pair< std::vector<EUTelHit> , std::vector<EUTelHit> > > tracksAndDUTHitsUnique;
//        std::cout<<" Tracks before unique: " << tracksAndDUTHits.size() << std::endl;
        tracksAndDUTHitsUnique = getUniqueMatches(tracksAndDUTHits);
//        std::cout<<" Tracks after unique " << tracksAndDUTHitsUnique.size() << std::endl;

        for(size_t i =0 ; i < tracksAndDUTHitsUnique.size() ; ++i){
            std::vector<EUTelHit> track = tracksAndDUTHitsUnique.at(i).first;
            std::vector<EUTelHit> dut = tracksAndDUTHitsUnique.at(i).second;

            std::vector<EUTelHit> combineHits;
            combineHits.reserve( tracks.size() + dut.size() ); // preallocate memory
            combineHits.insert( combineHits.end(), track.begin(), track.end() );
            combineHits.insert( combineHits.end(), dut.begin(), dut.end() );
            streamlog_out(DEBUG1) <<" Combined! "  << std::endl;
            ///Hit order does not matter
            tracks.push_back(EUTelTrackCreate::getTrackFourHits(combineHits));
        }
    }else{
  for(size_t i =0 ; i < tracksHits.size() ; ++i){
      tracks.push_back(EUTelTrackCreate::getTrackFourHits(tracksHits.at(i)));
  }
    }
    return tracks;

}
 std::vector< std::pair< std::vector<EUTelHit> , std::vector<EUTelHit> > > EUTelPatRecTriplets::getUniqueMatches( std::vector< std::pair< std::vector<EUTelHit> , std::vector<EUTelHit> > >& tracksAndDUTHits){
    //This really could be improved! Too many for loops!
    std::vector< std::pair< std::vector<EUTelHit> , std::vector<EUTelHit> > > tracksAndDUTHitsUnique;
    for(size_t i =0 ; i < tracksAndDUTHits.size() ; ++i){///Loop through all tracks and vector of DUT hits
        std::vector<EUTelHit> uniHits;
        if(tracksAndDUTHits.at(i).second.size() != 0 ){///If no DUT hit then pass track and don't check for DUT miss matching.
         //   std::cout<<"Hits on track "<< i << "  :" << tracksAndDUTHits.at(i).second.size() << std::endl;
            for(size_t j = 0; j < tracksAndDUTHits.at(i).second.size(); ++j){///Loop through all hits attached to track
                EUTelHit hit = tracksAndDUTHits.at(i).second.at(j);///This is the hit we are comparing to.
          //      std::cout<<"Check hit "<< j <<" from track " << i  << std::endl;

                bool unique = true;
                for(size_t k =i+1 ; k < tracksAndDUTHits.size() ; ++k){///Loop through all after the one compared to.
                    std::vector<EUTelHit>::iterator itHitMatch = std::find(tracksAndDUTHits.at(k).second.begin(),tracksAndDUTHits.at(k).second.end(),hit);
                    if(itHitMatch != tracksAndDUTHits.at(k).second.end()){///If this is not unique remove from this vector and continue search.
        //                std::cout<<"Not unique " << std::endl;
      //                  std::cout<<"Track " << k << "  compare size before remove"<<tracksAndDUTHits.at(k).second.size() << std::endl;
                        tracksAndDUTHits.at(k).second.erase(itHitMatch);
    //                    std::cout<<"Track " << k << "  compare size after remove"<<tracksAndDUTHits.at(k).second.size() << std::endl;

                        unique = false;
                    }else{
  //                      std::cout<<"Unique wrt track "<< k << std::endl;
                    }

                }
                if(unique){///If unique add to vector. 
                    uniHits.push_back(hit);
                }
                if(j == tracksAndDUTHits.at(i).second.size()-1){///If we are on the last hit of this track create track and add new unique hits. 
//                    std::cout<<"Hits on track after search:" << uniHits.size() << std::endl;

                    tracksAndDUTHitsUnique.push_back(std::make_pair(tracksAndDUTHits.at(i).first,uniHits));
                }
            }
        }else{
            tracksAndDUTHitsUnique.push_back(std::make_pair(tracksAndDUTHits.at(i).first,uniHits));
        }
    }
    return tracksAndDUTHitsUnique;
 }




//EUTelTrack EUTelPatRecTriplets::getHitsOnTrack(std::vector<double>& offset, std::vector<double>& trackSlope,double& qOverP, std::vector<unsigned int> sen, std::vector<EUTelHit>& hits){
//    std::vector<double> curvCorr; curvCorr.push_back(qOverP*EUTelNav::_bFac[0]);curvCorr.push_back(qOverP*EUTelNav::_bFac[1]);
//    double posX = offset.at(0) + dz1*trackSlope.at(0) + 0.5*dz1*dz2*curvCorr.at(0);
//    double posY = offset.at(1) + dz1*trackSlope.at(1) + 0.5*dz1*dz2*curvCorr.at(1);
//    for(std::vector<unsigned int> ::iterator itID = sen.begin(); itID != sen.end(); ++itID){
//        std::vector<EUTelHit> hits =  _mapHitsVecPerPlane.at(*itID);
//        float distBest = 10000000;
//        EUTelHit hitBest;
//        for(std::vector<EUTelHit>::iterator itHit = hits.begin(); itHit != hits.end(); ++itHit){
//            double hitPosX = itHit->getPositionGlobal()[0]; 
//            double hitPosY = itHit->getPositionGlobal()[1]; 
//            double hitPosZ = itHit->getPositionGlobal()[2]; 
//            std::vector<float>  pos = getDoubPosAtZ(doub, hitPosZ);/// Could calculate this once. Might be a bit off for tilted sensors.
//            float dist = getDistLocal(itHit, pos);
//            if(itHit == hits.begin()){
//                hitBest = *itHit;
//                distBest = dist;
//            }
//            if(dist < distBest){
//                hitBest = *itHit;
//                distBest = dist;
//            }
//        }
//        if(distBest >  _doubletCenDistCut.at(0)){
//            continue;
//        }
//        streamlog_out(DEBUG1) << "PASS Doublet cut!! " << std::endl;
//        newHits.push_back(hitBest);
//    }
//    ///Number of hits cut here.
//    return true;
//}


void EUTelPatRecTriplets::getDoublet(EUTelHit const & hit1, EUTelHit const & hit2,  doublets& doublet){
    const double curvX = EUTelNav::_curv.at(0); 
    const double curvY =  EUTelNav::_curv.at(1); 

    float initDis = geo::gGeometry().getInitialDisplacementToFirstPlane();
    //Remove the curvature as a factor between hits. Therefore only the slope will displace the hit position from plane to plane.
    float x1 = hit1.getPositionGlobal()[0] - 0.5*curvX*pow(hit1.getPositionGlobal()[2] - initDis, 2);
    float y1 = hit1.getPositionGlobal()[1] - 0.5*curvY*pow(hit1.getPositionGlobal()[2] - initDis, 2);
    float x2 = hit2.getPositionGlobal()[0] - 0.5*curvX*pow(hit2.getPositionGlobal()[2] - initDis, 2);
    float y2 = hit2.getPositionGlobal()[1] - 0.5*curvY*pow(hit2.getPositionGlobal()[2] - initDis, 2);

    doublet.pos.push_back((x2 + x1)/2.0);
    doublet.pos.push_back((y2 + y1)/2.0);
    doublet.pos.push_back((hit2.getPositionGlobal()[2] + hit1.getPositionGlobal()[2])/2.0);

    doublet.diff.push_back(x2 - x1);
    doublet.diff.push_back(y2 - y1);
    doublet.diff.push_back( hit2.getPositionGlobal()[2] - hit1.getPositionGlobal()[2]);

    doublet.slope.push_back( doublet.diff.at(0)/doublet.diff.at(2)); 
    doublet.slope.push_back( doublet.diff.at(1)/doublet.diff.at(2)); 
    doublet.hits.push_back(hit1);
    doublet.hits.push_back(hit2);
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
std::vector<float>  EUTelPatRecTriplets::getDoubPosAtZ(doublets doub, float posZ){
    float dz = posZ - doub.pos.at(2);
    float x = doub.pos.at(0) + doub.slope.at(0)*dz;
    float y = doub.pos.at(1) + doub.slope.at(1)*dz;
    std::vector<float> position;
    position.push_back(x);position.push_back(y);position.push_back(posZ);
    return position;
}



/*This creates map between plane ID and hits on that plane.  We also order the map correcly with geometry.
 */
void EUTelPatRecTriplets::setHitsVecPerPlane()
{
	_mapHitsVecPerPlane.clear();
	int numberOfPlanes = EUTelExcludedPlanes::_senInc.size();
	
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
			if(Utility::getSensorIDfromHit( static_cast<IMPL::TrackerHitImpl*>(_allHitsVec[j]) ) ==  EUTelExcludedPlanes::_senInc.at(i))
			{
				tempHitsVecPlaneX.push_back(EUTelHit(_allHitsVec.at(j)));
			}
		}
		_mapHitsVecPerPlane[  EUTelExcludedPlanes::_senInc.at(i)] = 	tempHitsVecPlaneX;
	}	
}
///Other 
void EUTelPatRecTriplets::printTrackQuality(std::vector<EUTelTrack>&  tracks )
{
    for(size_t i = 0 ; i < tracks.size(); i++){
        EUTelTrack& track = tracks.at(i);
        bool newTrack =true;
        for(size_t j = 0 ; j < track.getStates().size(); j++){
            EUTelState& state = track.getStates().at(j);
            if(state.getLocation() > 5 and state.getStateHasHit() and newTrack){
                _numberOfTracksDUTTotal = _numberOfTracksDUTTotal + 1; 
                newTrack=false;
            }
        }
    }

	_numberOfTracksTotal = _numberOfTracksTotal + tracks.size();
	if(_numberOfTracksTotal % 1000 == 0){
        streamlog_out(MESSAGE5) << "Number of tracks with DUT hit per event: " << static_cast<float>(_numberOfTracksDUTTotal)/static_cast<float>(getEventNumber() +1)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of tracks per event: " << static_cast<float>(_numberOfTracksTotal)/static_cast<float>(getEventNumber() +1)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of left arm triplets per event: " << static_cast<float>(_numberTripletsLeft)/static_cast<float>(getEventNumber() +1)<< std::endl;
        streamlog_out(MESSAGE5) << "Number of right arm triplets per event: " << static_cast<float>(_numberTripletsRight)/static_cast<float>(getEventNumber() +1)<< std::endl;
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

void EUTelPatRecTriplets::testHitsVecPerPlane(){
	if(_mapHitsVecPerPlane.size() !=  EUTelExcludedPlanes::_senInc.size()){
		streamlog_out(ERROR0) <<"The size of the planes with hits " << _mapHitsVecPerPlane.size() <<"  Sensors from Geometry with no excluded planes: "<<  EUTelExcludedPlanes::_senInc.size()<<std::endl;
		throw(lcio::Exception("The number of planes that could contain hits and the number of planes is different. Problem could be the gear file has to planes at the same z position.")); 	

	}
	for(size_t i=0 ;i<_mapHitsVecPerPlane.size();++i){
		int planeID = EUTelExcludedPlanes::_senInc.at(i);
		if(_mapHitsVecPerPlane.at(planeID).size() <= 0){
			streamlog_out(DEBUG1) << "One plane has no hits at all. Is this correct?" << std::endl;
		}
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

} // namespace eutelescope



