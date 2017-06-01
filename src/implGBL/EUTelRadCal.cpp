#include "EUTelRadCal.h"
using namespace eutelescope;


void EUTelRadCal::setRad(EUTelTrack& track, int& mode){
//    std::cout<<"setRad!!!" <<std::endl; 
    setIncSenBlocks(track);
    if(mode == 0){
        getVarForSensorScatterersOnly(track);
    }else{
        getThicknessAndRad(track); 
        setMeanWeight(track);
        setVarWeight(track);
        setPosVar(track);
    }

}
void EUTelRadCal::setIncSenBlocks(EUTelTrack & track){ 

    for(std::vector<EUTelState>::iterator  itSt = track.getStates().begin(); itSt != track.getStates().end(); ++itSt){

            double senRad =  geo::gGeometry().siPlaneRadLength(itSt->getLocation());
            double senSize =  geo::gGeometry().siPlaneZSize(itSt->getLocation());
            double radPer = geo::gGeometry().planeRadLengthGlobalIncidence(itSt->getLocation(), itSt->getDirGlobalEig()); 
            Block block;
            block.senRadPer = radPer;
            block.weigVar = 0;
            block.weigMean = 0;
            block.medRadPer = 0;
            itSt->block = block;
            track.setRadPerTotal(track.getRadPerTotal() + radPer);//VARIANCE MUST BE CALCULATED FROM THE TOTAL RADIATION LENGTH. Highland formula is non linear under addition.

    }

}

void EUTelRadCal::getThicknessAndRad(EUTelTrack & track){ 
    for(std::vector<EUTelState>::iterator itSt = track.getStates().begin(); itSt != (track.getStates().end()-1); ++itSt){//Loop over included sensors.
        //simplest way without z order to ID?
        int zOrdStartInc =find(geo::gGeometry().sensorIDsVec().begin(), geo::gGeometry().sensorIDsVec().end(), itSt->getLocation() ) - geo::gGeometry().sensorIDsVec().begin();
        int zOrdEndInc =  find(geo::gGeometry().sensorIDsVec().begin(), geo::gGeometry().sensorIDsVec().end(), (itSt+1)->getLocation() ) - geo::gGeometry().sensorIDsVec().begin();
        int diff= zOrdEndInc - zOrdStartInc - 1; 
        std::vector<std::pair<double,double> > thicknessAndRad; 
        int ordIncrease=0;
        for(int i =0 ; i < diff; ++i){ ///For excluded sensors: collect air then sensor, repeat ......
            int startAirID  = geo::gGeometry().sensorIDsVec().at(zOrdStartInc + i);
            int endAirID = geo::gGeometry().sensorIDsVec().at(zOrdStartInc +i +1);
            double thickness;
            double radAir = geo::gGeometry().airBetweenPlanesRadLengthGlobalIncidence(startAirID, endAirID , itSt->getDirGlobalEig() ,thickness);
            track.setRadPerTotal(track.getRadPerTotal() + radAir);
            itSt->block.medRadPer =  itSt->block.medRadPer + radAir;
            thicknessAndRad.push_back(std::make_pair(thickness,radAir));
            double excRad = geo::gGeometry().planeRadLengthGlobalIncidence(endAirID, itSt->getDirGlobalEig()); 
            double excSize =  geo::gGeometry().siPlaneZSize(endAirID);
            track.setRadPerTotal(track.getRadPerTotal() + excRad);
            itSt->block.medRadPer =  itSt->block.medRadPer + excRad;
            thicknessAndRad.push_back(std::make_pair(excSize,excRad));
            ordIncrease++;///This is needed so we can get the last air block after this loop.
        }
        const int startIDLastAirBlock = geo::gGeometry().sensorIDsVec().at(zOrdStartInc + ordIncrease); //This is the ID of the plane which starts the last air block
        double thickness;
        double radAir = geo::gGeometry().airBetweenPlanesRadLengthGlobalIncidence(startIDLastAirBlock, (itSt+1)->getLocation() , itSt->getDirGlobalEig() ,thickness);
        track.setRadPerTotal(track.getRadPerTotal() + radAir);
        itSt->block.medRadPer =  itSt->block.medRadPer + radAir;
        thicknessAndRad.push_back(std::make_pair(thickness,radAir));
        itSt->block.thicknessAndRad = thicknessAndRad;
    }
}

void EUTelRadCal::setMeanWeight(EUTelTrack & track){ 
    ///Add each homogeneous scatter as a new term in the for loop.
    for(std::vector<EUTelState>::iterator itSt = track.getStates().begin(); itSt != track.getStates().end()-1; ++itSt){
        double mean = 0;
        double norm = 0;
        std::vector<std::pair<double,double> > thRa = itSt->block.thicknessAndRad;
        for(std::vector<std::pair<double,double> >::iterator itTR = thRa.begin(); itTR != thRa.end(); ++itTR){
            mean = mean + 0.5*pow(itTR->first ,2)/itTR->second;        
            norm = norm + itTR->first/itTR->second; 
       //     std::cout<< "mean/norm "<< mean << " " << norm << " ID: " << itSt->getLocation() <<std::endl;
        }
        itSt->block.weigMean = mean/norm;
    }
}
void EUTelRadCal::setVarWeight(EUTelTrack & track){ 
    for(std::vector<EUTelState>::iterator itSt = track.getStates().begin(); itSt != track.getStates().end()-1; ++itSt){
        double var = 0;
        double norm = 0;
        std::vector<std::pair<double,double> > thRa = itSt->block.thicknessAndRad;
        for(std::vector<std::pair<double,double> >::iterator itTR = thRa.begin(); itTR != thRa.end(); ++itTR){

            var = var + ((1.0/3.0)*pow(itTR->first,3) - pow(itTR->first,2)*itSt->block.weigMean + pow(itSt->block.weigMean,2)*itTR->first)/itTR->second;
            norm = norm + itTR->first/itTR->second; 
        }
        itSt->block.weigVar = var/norm;
    }
}
void EUTelRadCal::setPosVar(EUTelTrack & track){ 
    for(std::vector<EUTelState>::iterator itSt = track.getStates().begin(); itSt != track.getStates().end(); ++itSt){
        getRelativePosOfScatterers(*itSt);
    }
    getVarForAllScatters(track);
}

void EUTelRadCal::getRelativePosOfScatterers(EUTelState & state){ 
    double firstScatRelPos =  0.51*geo::gGeometry().siPlaneZSize(state.getLocation());
  //  std::cout<<"mean " <<  state.block.weigMean << " var " << state.block.weigVar <<std::endl;
    double secondScatRelPos =  state.block.weigMean + state.block.weigVar/( state.block.weigMean - firstScatRelPos);
    state.block.scatPos.push_back(firstScatRelPos);state.block.scatPos.push_back(secondScatRelPos);
}

void EUTelRadCal::getVarForAllScatters(EUTelTrack & track ){ 
    ///Get the total variance of the system.
    double varTotal = pow(Utility::getThetaRMSHighland( track.getBeamEnergy(),track.getRadPerTotal()  ),2); 
    for(std::vector<EUTelState>::iterator itSt = track.getStates().begin(); itSt != track.getStates().end(); ++itSt){
        ///Included sensor fraction of variance.
        itSt->block.senVar = varTotal*(itSt->block.senRadPer/track.getRadPerTotal()); 
        ///Medium fraction of variance.
        double incMedVar = varTotal*(itSt->block.medRadPer/track.getRadPerTotal()); 
        double denom = itSt->block.weigVar + pow(itSt->block.weigMean - itSt->block.scatPos.at(0),2);
        itSt->block.scatVar.push_back(incMedVar*(itSt->block.weigVar/denom));
        itSt->block.scatVar.push_back(incMedVar*(pow(itSt->block.weigMean - itSt->block.scatPos.at(0),2)/denom));
    }
}
void EUTelRadCal::getVarForSensorScatterersOnly(EUTelTrack & track ){ 
    ///Get the total variance of the system.
    double varTotal = pow(Utility::getThetaRMSHighland( track.getBeamEnergy(),track.getRadPerTotal()  ),2); 
//    std::cout<<"varTotal "<< varTotal << " Beam energy :" <<  track.getBeamEnergy() << " RadPer " << track.getRadPerTotal() <<std::endl; 
    for(std::vector<EUTelState>::iterator itSt = track.getStates().begin(); itSt != track.getStates().end(); ++itSt){
        ///Included sensor fraction of variance.
        itSt->block.senVar = varTotal*(itSt->block.senRadPer/track.getRadPerTotal()); 
    }
}







