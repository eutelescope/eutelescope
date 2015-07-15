#include "EUTelRadCal.h"
using namespace eutelescope;


float EUTelRadCal::setHomoBlocks(TVector3& start , TVector3& end ,std::vector<int>& sen,std::map<int ,Block>& blocks){
    streamlog_out(DEBUG5) << "/////////////////////////////////////////////////////////////////////////////////////////////////// " << std::endl;
    streamlog_out(DEBUG5) << "/////////////////////////////////////////////////////////////////////////////////////////////////// " << std::endl;
    streamlog_out(DEBUG5) << "              CALCULATING THE TOTAL RADIATION LENGTH BETWEEN TWO POINTS.                            " << std::endl;
    streamlog_out(DEBUG5) << "                                 POINTS GO FROM:                            " << std::endl;
    streamlog_out(DEBUG5) << "              ("<< start[0] << ","  << start[1]<<","<<start[2]<<")"<< "-------------------------------------> ("<< end[0] << ","  << end[1]<<","<<end[2]<<")" << std::endl;
    TVector3 diff = end - start;
    diff.Unit();
    static geo::EUTelGeometryTelescopeGeoDescription& geo = geo::gGeometry();
	TGeoManager* gGeoManager = geo._geoManager;
    gGeoManager->InitTrack( start[0]/*mm*/, start[1]/*mm*/, start[2]/*mm*/, diff[0], diff[1], diff[2] ); //Start point and direction
    TGeoNode *node = gGeoManager->GetCurrentNode( );
    while ( node ) {
        ///This will call the current node and get the ID. Silly to access this way but leave for now.
        double rad,dist;
        move(node,rad,dist);
        fillBlocks(geo,rad,dist,blocks);
    }
}
float EUTelRadCal::setInHomoBlocks(std::map<int ,Block>& blocks, double beamEnergy){
    for(std::map<int, Block>::iterator itBl = blocks.begin(); itBl != blocks.end(); ++itBl){
        if(itBl->first < 0){///If block is medium then integrate.
            getDistWeig(beamEnergy,itBl);
            getVarWeig(beamEnergy,itBl);
        }
    }
}

void EUTelRadCal::move(TGeoNode *node, double & rad, double& dist){
    const double mm2cm = 0.1;
    TGeoMedium *med = NULL;
    if ( node ) med = node->GetVolume()->GetMedium();
    else throw(std::string("Passed non node in move EUTelRadCal! "));
    double radlen = med->GetMaterial()->GetRadLen() /*cm*/;
    double lastrad = 1. / radlen * mm2cm; //calculate 1/radiationlength per cm. This will transform radlen to mm
    node = gGeoManager->FindNextBoundaryAndStep( 1000 /*mm*/ ); //Some large number to check for new boundary. 
    dist  = gGeoManager->GetStep() /*mm*/; //This will output the distance traveled by FindNextBoundaryAndStep
    rad = 0; //This is the calculated (rad per distance x distance)
    double delta = 0.01;//This is the minimum block size 
    streamlog_out(DEBUG5)<<std::endl <<std::endl  << "DECISION: Step size over min?  "  <<" Block width: " << dist << " delta: " << delta  << std::endl;
    if(dist < delta){
       streamlog_out(DEBUG5) << "INCREASE TO MINIMUM DISTANCE!" << std::endl;
        dist = delta;
        double pt[3];
        memcpy( pt, gGeoManager->GetCurrentPoint(), 3 * sizeof (double) ); //Get global position
        const double *dir = gGeoManager->GetCurrentDirection();//Direction vector
        for ( Int_t i = 0; i < 3; i++ ) pt[i] += delta * dir[i]; //Move the current point slightly in the direction of motion. 
        node = gGeoManager->FindNode( pt[0], pt[1], pt[2] );//Move to new node where we will begin to look for more radiation length   
        rad=lastrad*dist; //Calculate radiation length for the increased block.
    }else{
        streamlog_out(DEBUG5) << "OVERMAX!" << std::endl;
        rad = lastrad*dist; //This is the calculated (rad per distance x distance)
    }
    streamlog_out(DEBUG5) << "NEW BLOCK:SensorID: " <<" Block width: " << dist << " Radiation total/per unit length: " << rad<<"/"<< std::endl;
}
/// Will take the EUTelescope object and get the sensor ID from the node.  
/// It will then update the radiation length and add the start and end points to the trajectory.
/**
 */

void EUTelRadCal::fillBlocks(geo::EUTelGeometryTelescopeGeoDescription& geo,double & rad, double& dist, std::map<int ,Block>& blocks){
    double pt[3];
    memcpy( pt, gGeoManager->GetCurrentPoint(), 3 * sizeof (double) ); //Get global position
    int senID = geo.getSensorIDFromManager();
    if(senID !=  *(_planes.rbegin()) and blocks.size() == 0   ){ ///Must not have reached first sensor.
        if(std::find(_planes.begin(),  _planes.end(), senID) !=  _planes.end()){//Should include this to sensor.
            try{///Is it a new block? Yes if exception. Else update sensor radiation length
                blocks.at(senID).totalRad = blocks.at(senID).totalRad + rad;
            }catch(std::out_of_range&  e){///Reached new sensor boundary. Add position and initial rad.
                blocks[senID].startPos.push_back(pt[0]);blocks[senID].startPos.push_back(pt[1]);blocks[senID].startPos.push_back(pt[2]);
                blocks[senID].totalRad = rad;
                if(blocks.size() != 0 ){//If this is not the first block it must be the end for thick scatterer block.
                    int id = -1.0*senID +1;///Add end point to the last medium.  
                    blocks[id].endPos.push_back(pt[0]);blocks[id].endPos.push_back(pt[1]);blocks[id].endPos.push_back(pt[2]);
                }
            }
        }else{///Must be medium inbetween sensors.
            int lastSenID =  blocks.rbegin()->first;//Last sensor we intersected.
            int idMat = -1.0*lastSenID;//Material infront of senID is mapped to -senID
            try{///Update If this has been created before.
                blocks.at(idMat).totalRad = blocks.at(idMat).totalRad + rad;
            }catch(std::out_of_range& e){///New medium block boundary. Add this to the end of the sensor block
                blocks[lastSenID].endPos.push_back(pt[0]);blocks[lastSenID].endPos.push_back(pt[1]);blocks[lastSenID].endPos.push_back(pt[2]);
                if(lastSenID != *(_planes.rbegin()) ){///If we have not got the last sensor then this must be saved. If so then do nothing.
                    blocks[idMat].totalRad = rad;
                    blocks[idMat].startPos.push_back(pt[0]);blocks[idMat].startPos.push_back(pt[1]);blocks[idMat].startPos.push_back(pt[2]);
                }else{
                    ///We have the last plane and now in next part of sensor
                    ///Do nothing
                 }
            }
        }
    }
}
void EUTelRadCal::setRad(EUTelTrack& track){
    _planes = track.getPlaIDs();
    streamlog_out(DEBUG3)<<"compute radiation..."<<std::endl;
    TVector3 gPosS =  track.getStates().at(0).getPositionGlobal();
    TVector3 gPosE =  track.getStates().back().getPositionGlobal();
    /// Take the start position just outside the first volume. So we include this in the calculation of radiaiton length.
    TVector3 start(gPosS[0],gPosS[1],gPosS[2]-0.05);
    TVector3 end(gPosE[0],gPosE[1],gPosE[2]+0.05);
    ///NOW WE CALCULATE THE RADIATION LENGTH FOR THE FULL FLIGHT AND THEN SPLIT THESE INTO LINEAR COMPONENTS FOR SCATTERING ESTIMATION. 
    streamlog_out(DEBUG3)<<"Use TGeo now!"<<std::endl;
    std::map<int ,Block> blocks;
    std::vector<int> pla = track.getPlaIDs();
    setHomoBlocks(start ,end ,pla,blocks );
    setInHomoBlocks(blocks, track.getBeamEnergy());
    _blocks =blocks;

}
void EUTelRadCal::getDistWeig(double& beamEnergy, std::map<int, Block>::iterator itBl){
    TVector3 end(itBl->second.endPos.at(0), itBl->second.endPos.at(1), itBl->second.endPos.at(2)); 
    TVector3 start(itBl->second.startPos.at(0), itBl->second.startPos.at(1), itBl->second.startPos.at(2)); 
    TVector3 diff = end  - start;
    diff.Unit();
    static geo::EUTelGeometryTelescopeGeoDescription& geo = geo::gGeometry();
    TGeoManager* gGeoManager = geo._geoManager;
    gGeoManager->InitTrack( start[0]/*mm*/, start[1]/*mm*/, start[2]/*mm*/, diff[0], diff[1], diff[2] ); //Start point and direction
    TGeoNode *node = gGeoManager->GetCurrentNode( );
    double radTotal =0;
    double distTotal=0;
    double weigMean=0;
    while ( node ) {
        double rad,dist;
        move(node,rad,dist);
        radTotal = radTotal +rad; //Used as check at the end.
        distTotal = distTotal + dist;//Keep total distance travelled since start
        double var =  pow(Utility::getThetaRMSHighland( beamEnergy, rad ),2);
        double weigMean = weigMean +  distTotal*var;//Update where mean should be placed.
        double pt[3];
        memcpy( pt, gGeoManager->GetCurrentPoint(), 3 * sizeof (double) ); //Get global position
        TVector3 newPos(pt[0], pt[1],pt[2] ); 
        if((end - newPos).Mag() ==  0){
            break;
        }
    }
    if(itBl->second.totalRad == radTotal){
        double varTot =  pow(Utility::getThetaRMSHighland( beamEnergy, radTotal ),2);
        itBl->second.weigMean = weigMean/varTot;
    }else{
        //Something must be wrong since the total radiation length should be the same as calculated in Homogenous block
        throw(std::string("Radition lengths from iteration do not match"));
    }

}

void EUTelRadCal::getVarWeig(double& beamEnergy, std::map<int, Block>::iterator itBl){
    TVector3 end(itBl->second.endPos.at(0), itBl->second.endPos.at(1), itBl->second.endPos.at(2)); 
    TVector3 start(itBl->second.startPos.at(0), itBl->second.startPos.at(1), itBl->second.startPos.at(2)); 
    TVector3 diff = end  - start;
    diff.Unit();
    static geo::EUTelGeometryTelescopeGeoDescription& geo = geo::gGeometry();
    TGeoManager* gGeoManager = geo._geoManager;
    gGeoManager->InitTrack( start[0]/*mm*/, start[1]/*mm*/, start[2]/*mm*/, diff[0], diff[1], diff[2] ); //Start point and direction
    TGeoNode *node = gGeoManager->GetCurrentNode( );
    double radTotal =0;
    double distTotal=0;
    double weigVar=0;
    while ( node ) {
        double rad,dist;
        move(node,rad,dist);
        radTotal = radTotal +rad; //Used as check at the end.
        distTotal = distTotal + dist;//Keep total distance travelled since start
        double var =  pow(Utility::getThetaRMSHighland( beamEnergy, rad ),2);
        double weigVar = weigVar +  pow(distTotal-itBl->second.weigMean ,2)*var;//Update where mean should be placed.
        double pt[3];
        memcpy( pt, gGeoManager->GetCurrentPoint(), 3 * sizeof (double) ); //Get global position
        TVector3 newPos(pt[0], pt[1],pt[2] ); 
        if((end - newPos).Mag() ==  0){
            break;
        }
    }
    if(itBl->second.totalRad == radTotal){
        double varTot =  pow(Utility::getThetaRMSHighland( beamEnergy, radTotal ),2);

        itBl->second.weigVar = weigVar/varTot;
    }else{
        //Something must be wrong since the total radiation length should be the same as calculated in Homogenous block
        throw(std::string("Radition lengths from iteration do not match"));
    }

}

