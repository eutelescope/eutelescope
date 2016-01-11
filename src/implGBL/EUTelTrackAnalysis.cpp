#include "EUTelTrackAnalysis.h"
using namespace eutelescope;

EUTelTrackAnalysis::EUTelTrackAnalysis(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY, std::map< int,  AIDA::IHistogram2D*> mapFromSensorIDHitMap, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToEfficiencyX, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToEfficiencyY, std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToGloIncXZ,std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToGloIncYZ, 		 std::map< int,   AIDA::IHistogram1D *> mapKinksX, std::map< int,   AIDA::IHistogram1D *> mapKinksY, std::map< int,  AIDA::IProfile2D* > mapFromSensorKinksMap, AIDA::IHistogram1D * beamEnergy){
setSensorIDTo2DHitMap(mapFromSensorIDHitMap);
    _mapKinksX =mapKinksX;
    _mapKinksY =mapKinksY;

    _mapFromSensorKinksMap = mapFromSensorKinksMap;

    setSensorIDTo2DResidualHistogramX(mapFromSensorIDToHistogramX);
    setSensorIDTo2DResidualHistogramY(mapFromSensorIDToHistogramY);
    setSensorIDTo2DResidualEfficiencyX(mapFromSensorIDToEfficiencyX);
    setSensorIDTo2DResidualEfficiencyY(mapFromSensorIDToEfficiencyY);
    setSensorIDToIncidenceAngleXZ(mapFromSensorIDToGloIncXZ);
    setSensorIDToIncidenceAngleYZ(mapFromSensorIDToGloIncYZ);
    setBeamEnergy(beamEnergy);


    for(unsigned int  j = 0; j < (geo::gGeometry().sensorIDsVec().size()); ++j){
        unsigned int sensorID = geo::gGeometry().sensorIDsVec().at(j);
        _senResTotX[sensorID] = 0.0;
        _senResTotY[sensorID] = 0.0;
        _senResTotZ[sensorID] = 0.0;
        _hitNum[sensorID] = 0;

    }



} 
void EUTelTrackAnalysis::plotKinksVsPosition(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotResidualVsPosition------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		TVector3 statePositionGlobal = state.getPositionGlobal();
		typedef std::map<int ,AIDA::IProfile2D*  >::iterator it_type;
		for(it_type iterator = _mapFromSensorKinksMap.begin(); iterator != _mapFromSensorKinksMap.end(); iterator++) {
		  streamlog_out(DEBUG2) << "	state.getLocation() = "<<	state.getLocation()<< ", iterator->first = "<<iterator->first<<std::endl;	
		  if(iterator->first == state.getLocation()){
            //  std::cout<<"ID: " << state.getLocation() << " pos  " << statePositionGlobal[0] << " " <<  statePositionGlobal[1] << " kink: " <<  state.getKinks()[0] <<std::endl;
           //   if(state.getKinks()[0] < 0 ){
                _mapFromSensorKinksMap[ state.getLocation() ]  -> fill( statePositionGlobal[0], statePositionGlobal[1], state.getKinks()[0] );
           //   }
			break;
            }
		}
	} 
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotResidualVsPosition------------------------------END"<< std::endl;
}


void EUTelTrackAnalysis::plotResidualVsPosition(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotResidualVsPosition------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		if(!state.getStateHasHit()){
			continue;
		}
		EUTelHit hit = state.getHit();	
		const double* statePosition = state.getPosition();
		TVector3 statePositionGlobal = state.getPositionGlobal();
		const double* hitPosition = hit.getPosition();
		float residual[2];
		if(statePositionGlobal[0]>8)streamlog_out(DEBUG2) << "*** HELEN ***** statePositionGlobal[0]>8 ********"<< std::endl;
		if(statePosition[0]<-5)streamlog_out(DEBUG2) << "*** HELEN ***** statePosition[0]<-5 ********"<< std::endl;
		streamlog_out(DEBUG2) << "state.getLocation() = "<<	state.getLocation()<< std::endl;
		streamlog_out(DEBUG2) << "State position: " << statePosition[0]<<","<<statePosition[1]<<","<<statePosition[2]<< std::endl;
streamlog_out(DEBUG2) << "State position glo: " << statePositionGlobal[0]<<","<<statePositionGlobal[1]<<","<<statePositionGlobal[2]<< std::endl;
		streamlog_out(DEBUG2) << "Hit position: " << hitPosition[0]<<","<<hitPosition[1]<<","<<hitPosition[2]<< std::endl;

		typedef std::map<int ,AIDA::IProfile2D*  >::iterator it_type;
		for(it_type iterator = _mapFromSensorIDToHistogramX.begin(); iterator != _mapFromSensorIDToHistogramX.end(); iterator++) {
		  streamlog_out(DEBUG2) << "	state.getLocation() = "<<	state.getLocation()<< ", iterator->first = "<<iterator->first<<std::endl;	
		  if(iterator->first == state.getLocation()){
			 
			residual[0]=std::abs(statePosition[0]-hitPosition[0]);
			streamlog_out(DEBUG2) << "Add residual X : " << residual[0]<< std::endl;

			_mapFromSensorIDToHistogramX[ state.getLocation() ]  -> fill( statePositionGlobal[0], statePositionGlobal[1], residual[0], 1 );
			break;
			}
		}
		for(it_type iterator = _mapFromSensorIDToHistogramY.begin(); iterator != _mapFromSensorIDToHistogramY.end(); iterator++) {
			if(iterator->first == state.getLocation()){

			residual[1]=std::abs(statePosition[1]-hitPosition[1]);
			streamlog_out(DEBUG2) << "Add residual Y : " << residual[1]<< std::endl;
			_mapFromSensorIDToHistogramY[ state.getLocation() ]  -> fill( statePositionGlobal[0], statePositionGlobal[1], residual[1], 1 );
			break;
			}
		}
	} 
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotResidualVsPosition------------------------------END"<< std::endl;
}

void EUTelTrackAnalysis::plotHitMap(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotHitMap------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		const double* statePosition = state.getPosition();
		TVector3 statePositionGlobal = state.getPositionGlobal();
	streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotHitMap------------------------------still here"<< std::endl;

		typedef std::map<int ,AIDA::IHistogram2D*  >::iterator it_type;
		for(it_type iterator = _mapFromSensorIDHitMap.begin(); iterator != _mapFromSensorIDHitMap.end(); iterator++) {
		  streamlog_out(DEBUG2) << "	state.getLocation() = "<<	state.getLocation()<< ", iterator->first = "<<iterator->first<<std::endl;	
		  if(iterator->first == state.getLocation()){

			_mapFromSensorIDHitMap[ state.getLocation() ]  -> fill(statePositionGlobal[0], statePositionGlobal[1],1.);
			break;
			}
		}
		
	} 
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotHitMap------------------------------END"<< std::endl;
}
void EUTelTrackAnalysis::plotEfficiencyVsPosition(EUTelTrack track, IntVec sensorIDs){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotEfficiencyVsPosition------------------------------BEGIN"<< std::endl;
  //need function to retrieve maximum and minimum Xand Y positions of hit co-ordinates
  //need array for 8 locations
  std::map< int, float > minhitx;
  std::map< int, float > maxhitx;
  std::map< int, float > minhity;
  std::map< int, float > maxhity;
  for (size_t i = 0; i < sensorIDs.size() ; ++i){
    minhitx.insert(std::make_pair(sensorIDs.at(i),0.));
    maxhitx.insert(std::make_pair(sensorIDs.at(i),0.));
    minhity.insert(std::make_pair(sensorIDs.at(i),0.));
    maxhity.insert(std::make_pair(sensorIDs.at(i),0.));
  }
  //this is in the wrong place. need it once per job - not per event
  //separate function? calculate here to start. 
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
			const double* statePosition = state.getPosition();
		const double* hitPosition = NULL;
		streamlog_out(DEBUG0) << " In state loop, "<<i<<" out of "<< states.size()<<" state.getStateHasHit() = "<<state.getStateHasHit()<<std::endl;
		EUTelHit hit;
		//statePosition = state.getPosition();
		TVector3 statePositionGlobal = state.getPositionGlobal();
		if(state.getStateHasHit()){
		  hit = state.getHit();	
		  //streamlog_out(DEBUG2)<<"hit on location "<<state.getLocation()<<" has time "<<hit.getTime()<<std::endl;
		  hitPosition = hit.getPosition();
		  streamlog_out(DEBUG2) << " HELEN , hitPosition[0] = "<<hitPosition[0]<<"minhitx[state.getLocation()] = "<< minhitx[state.getLocation()]<<std::endl;
		  if(hitPosition[0]<minhitx[state.getLocation()]){minhitx[state.getLocation()]=hitPosition[0];}
		  if(hitPosition[0]>maxhitx[state.getLocation()]){maxhitx[state.getLocation()]=hitPosition[0];}
		  if(hitPosition[1]<minhity[state.getLocation()]){minhity[state.getLocation()]=hitPosition[1];}
		  if(hitPosition[1]>maxhity[state.getLocation()]){maxhity[state.getLocation()]=hitPosition[1];}
 // 		  if(hitPosition[0]>maxhitx->at(state.getLocation())){maxhitx->at(state.getLocation())=hitPosition[0];}
 // 		  if(hitPosition[1]<minhity->at(state.getLocation())){minhity->at(state.getLocation())=hitPosition[1];}
 // 		  if(hitPosition[1]>maxhity->at(state.getLocation())){maxhity->at(state.getLocation())=hitPosition[1];}
		  streamlog_out(DEBUG2) << " HELEN , hitPosition[0] = "<<hitPosition[0]<<"minhitx[state.getLocation()] = "<< minhitx[state.getLocation()]<<std::endl;
		  streamlog_out(DEBUG2) << " HELEN , hitPosition[0] = "<<hitPosition[0]<<"maxhitx[state.getLocation()] = "<< maxhitx[state.getLocation()]<<std::endl;
		  streamlog_out(DEBUG2) << " HELEN , hitPosition[1] = "<<hitPosition[1]<<"minhity[state.getLocation()] = "<< minhity[state.getLocation()]<<std::endl;
		  streamlog_out(DEBUG2) << " HELEN , hitPosition[1] = "<<hitPosition[1]<<"maxhity[state.getLocation()] = "<< maxhity[state.getLocation()]<<std::endl;
		//	continue;
		}
		else streamlog_out(DEBUG2) << " state.getStateHasHit() = "<<state.getStateHasHit()<<std::endl;
		float residual[2];
		if(statePositionGlobal[0]>8)streamlog_out(DEBUG2) << "*** HELEN ***** statePositionGlobal[0]>8 ********"<< std::endl;
		if(statePosition[0]<-5)streamlog_out(DEBUG2) << "*** HELEN ***** statePosition[0]<-5 ********"<< std::endl;
		streamlog_out(DEBUG2) << "Looking at location = "<<state.getLocation()<< std::endl;
		streamlog_out(DEBUG2) << "State position loc: " << statePosition[0]<<","<<statePosition[1]<<","<<statePosition[2]<< std::endl;
		streamlog_out(DEBUG2) << "State position glo: " << statePositionGlobal[0]<<","<<statePositionGlobal[1]<<","<<statePositionGlobal[2]<< std::endl;
		if(state.getStateHasHit()){
		streamlog_out(DEBUG2) << "Hit position: " << hitPosition[0]<<","<<hitPosition[1]<<","<<hitPosition[2]<< std::endl;
		}
		else 	streamlog_out(DEBUG2) << "No hit"<<std::endl;

		typedef std::map<int ,AIDA::IProfile2D*  >::iterator it_type;
		for(it_type iterator = _mapFromSensorIDToEfficiencyX.begin(); iterator != _mapFromSensorIDToEfficiencyX.end(); iterator++) {


		  bool hasMatchedXHit = false;
			if(iterator->first == state.getLocation()){
			  //check if hit near this location
			  if(state.getStateHasHit()){

			    if(std::abs(statePosition[0]-hitPosition[0])<0.2/*&&std::abs(statePosition[1]-hitPosition[1])<0.5*/){
			      hasMatchedXHit=true;
			      streamlog_out(DEBUG0) << "then we have hit!"<< std::endl;
			       }
			  }else streamlog_out(DEBUG0) << "!state.getStateHasHit()"<<std::endl;
			
			//residual[0]=statePosition[0]-hitPosition[0];
			streamlog_out(DEBUG0) << "Add efficiency X : " << residual[0]<< std::endl;

			_mapFromSensorIDToEfficiencyX[ state.getLocation() ]  -> fill( statePositionGlobal[0], statePositionGlobal[1], hasMatchedXHit, 1 );
			//	if(//if safely inside DUT
			//   _EffX[state.getLocation() ]-> fill(hasMatchedXHit);
			break;
			}
		}
		for(it_type iterator = _mapFromSensorIDToEfficiencyY.begin(); iterator != _mapFromSensorIDToEfficiencyY.end(); iterator++) {
		  bool hasMatchedYHit = false;
			if(iterator->first == state.getLocation()){
			  if(state.getStateHasHit()){
			    if(fabs(statePosition[0]-hitPosition[0])<0.5/*&&fabs(statePosition[1]-hitPosition[1])<0.5*/){
			      hasMatchedYHit=true;
			      //then we have hit!
			       }
			  }
			  //	residual[1]=statePosition[1]-hitPosition[1];
			  streamlog_out(DEBUG0) << "Add efficeincy Y : " << residual[1]<< std::endl;
			  _mapFromSensorIDToEfficiencyY[ state.getLocation() ]  -> fill( statePosition[0], statePosition[1], hasMatchedYHit, 1 );
			  break;
			}
		}
	} 
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotEfficiencyVsPosition------------------------------END"<< std::endl;
}
void EUTelTrackAnalysis::plotBeamEnergy(EUTelTrack track){
    streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotBeamEnergy------------------------------BEGIN"<< std::endl;
//    std::cout << "Track energy: " << track.getBeamEnergy() << std::endl;
    _beamEnergy-> fill(track.getBeamEnergy() );
    streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotBeamEnergy------------------------------END"<< std::endl;
}
void EUTelTrackAnalysis::plotPValueVsBeamEnergy(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueVsBeamEnergy------------------------------BEGIN"<< std::endl;
	float pValue = calculatePValueForChi2(track);
	
	_pValueVsBeamEnergy->fill(track.getBeamEnergy(), pValue);
	streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueVsBeamEnergy------------------------------END pvalue = "<< pValue<< std::endl;
}


void EUTelTrackAnalysis::plotIncidenceAngles(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotIncidenceAngles------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		float incidenceXZ = state.getSlopeXGlobal();
		typedef std::map<int , AIDA::IHistogram1D * >::iterator it_type;
		for(it_type iterator =_mapFromSensorIDToIncidenceXZ.begin(); iterator != _mapFromSensorIDToIncidenceXZ.end(); iterator++) {
			if(iterator->first == state.getLocation()){
			streamlog_out(DEBUG2) << "Add incidence XZ : " << incidenceXZ  << std::endl;
			_mapFromSensorIDToIncidenceXZ[ state.getLocation() ] ->fill(incidenceXZ);
			break;
			}
		}
	} 
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		float incidenceYZ = state.getSlopeYGlobal();
		typedef std::map<int , AIDA::IHistogram1D * >::iterator it_type;
		for(it_type iterator =_mapFromSensorIDToIncidenceYZ.begin(); iterator != _mapFromSensorIDToIncidenceYZ.end(); iterator++) {
			if(iterator->first == state.getLocation()){
			streamlog_out(DEBUG2) << "Add incidence YZ : " << incidenceYZ  << std::endl;
			_mapFromSensorIDToIncidenceYZ[ state.getLocation() ] ->fill(incidenceYZ);
			break;
			}
		}
	}
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotIncidenceAngles------------------------------END"<< std::endl;
}
void EUTelTrackAnalysis::plotKinks(EUTelTrack& track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotIncidenceAngles------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		float kinksX = state.getKinks()[0];
		float kinksY = state.getKinks()[1];
		typedef std::map<int , AIDA::IHistogram1D * >::iterator it_type;
		for(it_type iterator =_mapKinksX.begin(); iterator != _mapKinksX.end(); iterator++) {
			if(iterator->first == state.getLocation()){
                _mapKinksX[ state.getLocation() ] ->fill(kinksX);
                _mapKinksY[ state.getLocation() ] ->fill(kinksY);
                break;
			}
		}
	} 
}

void EUTelTrackAnalysis::plotPValueWithIncidenceAngles(EUTelTrack track){
	streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithIncidenceAngles------------------------------BEGIN"<< std::endl;
	float pValue = calculatePValueForChi2(track);
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		float incidenceXZ = state.getSlopeXGlobal();
		typedef std::map<int , AIDA::IProfile1D * >::iterator it_type;
		for(it_type iterator =_mapFromSensorIDToPValuesVsIncidenceXZ.begin(); iterator != _mapFromSensorIDToPValuesVsIncidenceXZ.end(); iterator++) {
			if(iterator->first == state.getLocation()){
			  streamlog_out(DEBUG2) << "Add incidence XZ : " << incidenceXZ  << "pValue = "<<pValue<<std::endl;
			_mapFromSensorIDToPValuesVsIncidenceXZ[ state.getLocation() ] ->fill(incidenceXZ, pValue);
			break;
			}
		}
	} 
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		float incidenceYZ = state.getSlopeYGlobal();
		typedef std::map<int , AIDA::IProfile1D * >::iterator it_type;
		for(it_type iterator =_mapFromSensorIDToPValuesVsIncidenceYZ.begin(); iterator != _mapFromSensorIDToPValuesVsIncidenceYZ.end(); iterator++) {
			if(iterator->first == state.getLocation()){
			streamlog_out(DEBUG2) << "Add incidence YZ : " << incidenceYZ  << std::endl;
			_mapFromSensorIDToPValuesVsIncidenceYZ[ state.getLocation() ] ->fill(incidenceYZ, pValue);
			break;
			}
		}
	} 
	streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithIncidenceAngles------------------------------END"<< std::endl;
}


void EUTelTrackAnalysis::plotPValueWithPosition(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithPosition------------------------------BEGIN"<< std::endl;
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithPosition------------------------------BEGIN"<< std::endl;
	float pValue = calculatePValueForChi2(track);
	std::vector<EUTelState> states = track.getStates(); streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithPosition------------------------------have states"<< std::endl;
	for(size_t i=0; i<states.size();++i){
	  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithPosition------------------------------inside loop i="<<i<< std::endl;
		EUTelState state  = states.at(i);
		state.print();

		TVector3 statePosition = state.getPositionGlobal();
		streamlog_out(DEBUG2) << "State position: " << statePosition[0]<<","<<statePosition[1]<<","<<statePosition[2]<< std::endl;

		typedef std::map<int ,AIDA::IProfile2D*  >::iterator it_type;
		for(it_type iterator = _mapFromSensorIDTo2DPValuesWithPosition.begin(); iterator != _mapFromSensorIDTo2DPValuesWithPosition.end(); iterator++) {
			if(iterator->first == state.getLocation()){
				_mapFromSensorIDTo2DPValuesWithPosition[ state.getLocation() ]  -> fill( statePosition[0], statePosition[1],pValue , 1 );
				break;
			}
		} 
	}
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithPosition------------------------------END"<< std::endl;
}
float EUTelTrackAnalysis::calculatePValueForChi2(EUTelTrack track){//std::cout<<"pigeon"<<std::endl;
//  std::cout<<"pigeon = "<<track.getNdf()<<std::endl;
  if(track.getNdf()==0) return 1;
  boost::math::chi_squared mydist(track.getNdf());//std::cout<<"pigeon2"<<std::endl;
  float pValue = 1 - boost::math::cdf(mydist,track.getChi2());//std::cout<<"pigeon3"<<std::endl;
    return pValue;
}
void EUTelTrackAnalysis::print(){
    streamlog_out(MESSAGE9) << "Analysis Results: " << std::endl;
    for(size_t  j = 0; j < (geo::gGeometry().sensorIDsVec().size()); ++j){
        unsigned int sensorID = geo::gGeometry().sensorIDsVec().at(j);
        ///Output all planes even if excluded. Add one to hit count so we do not divide by zero.
        streamlog_out(MESSAGE9) <<"Sensor " << sensorID << " Absolute residual average X " <<  _senResTotX[sensorID]/static_cast<float>(_hitNum[sensorID]+1) << " residual average Y " <<  _senResTotY[sensorID]/static_cast<float>(_hitNum[sensorID]+1)<< " hit number " << _hitNum[sensorID]<<std::endl;;

    }
}



void EUTelTrackAnalysis::setTotNum(EUTelTrack& track){
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		if(state.getStateHasHit()){
            state.print();
            int ID = state.getLocation();
            float resX = state.getPositionGlobal()[0] - state.getHit().getPositionGlobal()[0];
            float resY = state.getPositionGlobal()[1] - state.getHit().getPositionGlobal()[1];
            float resZ = state.getPositionGlobal()[2] - state.getHit().getPositionGlobal()[2];
            _senResTotX[ID] = _senResTotX[ID] + pow(resX,2); 
            _senResTotY[ID] = _senResTotY[ID] + pow(resY,2); 
            _senResTotZ[ID] = _senResTotZ[ID] + pow(resZ,2); 
            _hitNum[ID] = _hitNum[ID] + 1; 
        }
	} 
}


//FLOAT EUTelTrackAnalysis::calculatePValueForChi2(EUTelTrack track){
//  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::calculatePValueForChi2------------------------------BEGIN"<< std::endl;
//	float chi2Float=track.getChi2();
//	int   ndfInt = track.getNdf();
//	std::string chi2 = numberToString(chi2Float);
//	std::string ndf = numberToString(ndfInt);
//	float pValue=0;
////	std::cout << "Chi2: " << chi2 <<" and ndf " << ndf <<std::endl;
//	const std::string command = "calculatePValue.pyc " + chi2 + " " + ndf;
//	redi::ipstream pValueStream( command.c_str( ));
//	if ( !pValueStream.is_open( )){
//		throw(lcio::Exception("Could not open the pValue file. "));
//	}else{
//		std::string str;
//		while (pValueStream >> str) {
////		std::cout << str << std::endl;
//		pValue = std::atof(str.c_str());
//		}
//	}
//	pValueStream.close( );
//  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::calculatePValueForChi2------------------------------END"<< std::endl;
////	std::cout<<"Here is the p-value: " <<pValue <<std::endl;
//	return pValue;
//}
//template<typename T>
//std::string EUTelTrackAnalysis::numberToString(T number){
//	std::string Result;        
//	std::ostringstream convert;
//	convert << number;   
//	Result = convert.str();
////	return Result;
//	return 1.0;
//}



//void EUTelTrackAnalysis::plotPValueWithPosition(EUTelTrack track){
//    streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithPosition------------------------------BEGIN"<< std::endl;
//    float pValue = calculatePValueForChi2(track);
//    std::vector<EUTelState> states = track.getStates();
//    for(size_t i=0; i<states.size();++i){
//        EUTelState state  = states.at(i);
//        state.print();
//
//        float* statePosition = state.getPosition();
//        streamlog_out(DEBUG2) << "State position: " << statePosition[0]<<","<<statePosition[1]<<","<<statePosition[2]<< std::endl;
//
//        typedef std::map<int ,AIDA::IProfile2D*  >::iterator it_type;
//        for(it_type iterator = _mapFromSensorIDTo2DPValuesWithPosition.begin(); iterator != _mapFromSensorIDTo2DPValuesWithPosition.end(); iterator++) {
//            if(iterator->first == state.getLocation()){
//                _mapFromSensorIDTo2DPValuesWithPosition[ state.getLocation() ]  -> fill( statePosition[0], statePosition[1],pValue , 1 );
//                break;
//            }
//        } 
//    }
//    streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithPosition------------------------------END"<< std::endl;
//}
//
