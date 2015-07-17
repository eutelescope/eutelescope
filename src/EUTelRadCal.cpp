#include "EUTelRadCal.h"
using namespace eutelescope;


std::map<int ,Block> EUTelRadCal::getRad(EUTelTrack& track){
    std::map<int ,Block> blocks;
    setIncSenBlocks(track,blocks);
//    getThicknessAndRad(track,blocks); 
//    getScatParam(blocks); 
    return blocks;

}
void EUTelRadCal::setIncSenBlocks(EUTelTrack const & track, std::map<int, Block>& blocks){ 
    for(std::vector<EUTelState>::iterator  itSt = track.getStatesCopy().begin(); itSt != track.getStatesCopy().end(); ++itSt){
            double senRad =  geo::gGeometry().siPlaneRadLength(itSt->getLocation());
            double senSize =  geo::gGeometry().siPlaneZSize(itSt->getLocation());
            ///Must check values for this method of access
            blocks[itSt->getLocation()].senRadPer = senSize/senRad;
            blocks.at(itSt->getLocation()).weigVar = 0;
            blocks.at(itSt->getLocation()).weigMean = 0;
            blocks.at(itSt->getLocation()).medRadPer = 0;


    }
}

void EUTelRadCal::getThicknessAndRad(EUTelTrack const & track, std::map<int, Block> & blocks){ 
    std::map<int,std::vector<std::pair<double,double> > > IDToThickAndRad;
    /// to air = 37.15 g cm/2   
    /// Density:0.001225 g/cm3
    /// 30326 radiaiton length cm-1
    ///Should get this from TGeo. How??
    double radAir = 30326;
    for(std::map<int,Block>::iterator itBl = blocks.begin(); itBl != blocks.end(); ++itBl){//Loop over included sensors.
        int zOrdStartInc =  geo::gGeometry().sensorIDtoZOrder( itBl->first );
        int zOrdEndInc =  geo::gGeometry().sensorIDtoZOrder( (itBl++)->first );
        itBl--; //This is silly but map is not random access iterator. How should I do this?
        int diff= zOrdEndInc - zOrdStartInc - 1; 
        std::vector<std::pair<double,double> > thicknessAndRad; 
        int ordIncrease=0;
        for(int i =0 ; i < diff; ++i){ ///For excluded sensors: collect air then sensor, repeat ......
            int startAirID  = geo::gGeometry().sensorIDsVec().at(zOrdStartInc + i);
            int endAirID = geo::gGeometry().sensorIDsVec().at(zOrdStartInc +i +1);
            double startZPosAir =  geo::gGeometry().siPlaneZPosition(startAirID);
            double endZPosAir =  geo::gGeometry().siPlaneZPosition(endAirID);
            double sizeAir = endZPosAir - startZPosAir;
            thicknessAndRad.push_back(std::make_pair(sizeAir,radAir));
            double excSize =  geo::gGeometry().siPlaneZSize(endAirID);
            double excRad =  geo::gGeometry().siPlaneRadLength(endAirID);
            thicknessAndRad.push_back(std::make_pair(excSize,excRad));
            ordIncrease++;
        }
        ///Add the last block of air before the included sensor.
        double sizeAir =  geo::gGeometry().siPlaneZPosition((itBl++)->first) - geo::gGeometry().siPlaneZPosition( geo::gGeometry().sensorIDsVec().at(zOrdStartInc + ordIncrease));
        itBl--;
        thicknessAndRad.push_back(std::make_pair(sizeAir,radAir));
        itBl->second.thicknessAndRad = thicknessAndRad;
    }
}

void EUTelRadCal::getScatParam(std::map<int, Block> & blocks){ 
    for(std::map<int,Block>::iterator itBl = blocks.begin(); itBl != blocks.end(); ++itBl){
        setMeanWeight(itBl->second);
        setVarWeight(itBl->second);
    }
}
void EUTelRadCal::setMeanWeight(Block & block){ 
    double mean = 0;
    double norm = 0;
    ///Add each homogeneous scatter as a new term in the for loop.
    for(std::vector<std::pair<double, double > >::iterator itThiRad = block.thicknessAndRad.begin (); itThiRad != block.thicknessAndRad.end() ; ++itThiRad){
        mean = mean + 0.5*pow(itThiRad->first,2)/itThiRad->second;        
        norm = norm + itThiRad->first/itThiRad->second; 
    }
    block.weigMean = mean/norm;
}
void EUTelRadCal::setVarWeight(Block & block){ 
    double var = 0;
    double norm = 0;
    for(std::vector<std::pair<double, double > >::iterator itThiRad = block.thicknessAndRad.begin (); itThiRad != block.thicknessAndRad.end() ; ++itThiRad){
        var = var + ((1.0/3.0)*pow(itThiRad->first,3) - pow(itThiRad->first,2)*block.weigMean + pow(block.weigMean,2)*itThiRad->first)/itThiRad->second;
        norm = norm + itThiRad->first/itThiRad->second; 
    }
    block.weigMean = var/norm;
}





