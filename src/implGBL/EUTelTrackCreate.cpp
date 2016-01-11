#include "EUTelTrackCreate.h"	

namespace eutelescope 
{
EUTelTrack EUTelTrackCreate::getTrackFourHits(std::vector<EUTelHit> hits ){
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
    double qOverPCorr;
    //Find correction of curvature through slope change.
    const gear::BField& B = geo::gGeometry().getMagneticField();
    const double Bmag = B.at( TVector3(0.,0.,0.) ).r2();
    if ( Bmag < 1.E-6 ){
        qOverPCorr = 0;
    }else{
        qOverPCorr =  EUTelNav::getCorr(hitArmOne1,hitArmOne2,hitArmTwo1,hitArmTwo2);
    }
    double qOverP = -1.0/EUTelNav::_intBeamE+ qOverPCorr;
    //NOW CREATE TRACK CANDIDATE
    std::vector<double> offset;
    std::vector<double> trackSlope; 
    EUTelNav::getTrackAvePara(hitArmOne1, hitArmTwo2, offset, trackSlope);
    EUTelTrack track = getTrack(hits,offset,trackSlope,qOverP);
    return track;

}
EUTelTrack EUTelTrackCreate::getTrack(std::vector<EUTelHit> hits, std::vector<double> offset, std::vector<double> trackSlope,double qOverP){
    EUTelTrack track;
    track.setQOverP(qOverP);
    std::vector<double> curvCorr; curvCorr.push_back(track.getQOverP()*EUTelNav::_bFac[0]);curvCorr.push_back(track.getQOverP()*EUTelNav::_bFac[1]);
    // loop around planes
    // loop around hits -> if hit matched to plane record hit, if not record my intersection
    for(unsigned int  j = 0; j < (EUTelExcludedPlanes::_senInc.size()); ++j){
        unsigned int sensorID = EUTelExcludedPlanes::_senInc.at(j);
        streamlog_out(DEBUG1) << "The Z position " << j << " sensor ID: " << sensorID  <<std::endl;
        bool hitOnPlane=false;
        EUTelState state; //Create a state for each plane included in the fit.
        for(unsigned int i = 0; i < hits.size(); ++i){//Check the list of hits to see if we have one on this plane.
            if(hits.at(i).getLocation()==sensorID){
                hitOnPlane=true;
                Eigen::Vector3d posPred;
                std::vector<double> slopePred;
                EUTelNav::getTrackPredictionFromParam(offset,trackSlope,qOverP, hits.at(i).getPositionGlobal()[2],posPred,slopePred);
                state.setLocation(hits.at(i).getLocation());
                state.setHit(hits.at(i));
                float intersectionPoint[3];
                intersectionPoint[0] = posPred[0];  intersectionPoint[1] = posPred[1]; intersectionPoint[2] = hits.at(i).getPositionGlobal()[2];
                //intersection might not be inside a volume. 
                state.setPositionGlobal(intersectionPoint);
                state.setDirFromGloSlope(slopePred);
                track.setState(state);
            //    std::cout<<"Added hit! " << " z pos: " << i  <<std::endl;
            }
        }
        if(hitOnPlane){streamlog_out(DEBUG1) <<"this should already be recorded as there is a hit on sensor "<<sensorID<<std::endl;}
        else {
        //    std::cout<<"No hit ID " <<sensorID <<std::endl;
            double z_dut=geo::gGeometry().getOffsetVector(sensorID)[2];
            Eigen::Vector3d posPred;
            std::vector<double> slopePred;
            EUTelNav::getTrackPredictionFromParam(offset,trackSlope,qOverP, z_dut,posPred,slopePred);
            state.setDirFromGloSlope(slopePred);
            //   TVector3 hitPosGlo = hit.getPositionGlobal();
            float intersectionPoint[3];
            intersectionPoint[0] = posPred[0];  intersectionPoint[1] = posPred[1]; intersectionPoint[2] = z_dut;
            //   //intersection might not be inside a volume. 
            streamlog_out(DEBUG1)<<"intersection point on sensorID "<<sensorID<<" = "<<	intersectionPoint[0]<<", "<<intersectionPoint[1]<<", "<<intersectionPoint[2]<<std::endl;
            //add explicit check that it intersects with sensor?   but i want edge effects?   
            //add arc length thingy
            state.setLocation(sensorID);
            state.setPositionGlobal(intersectionPoint);
            track.setState(state);

        }//else
    }//loop about planes, j iterator
    track.print();
    return track;
}


} //namespace eutelescope
