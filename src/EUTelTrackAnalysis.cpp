#include "EUTelTrackAnalysis.h"
using namespace eutelescope;
EUTelTrackAnalysis::EUTelTrackAnalysis(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY, std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkXZ,std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkYZ,  AIDA::IHistogram1D * beamEnergy){
setSensorIDTo2DResidualHistogramX(mapFromSensorIDToHistogramX);
setSensorIDTo2DResidualHistogramY(mapFromSensorIDToHistogramY);
setSensorIDToIncidenceAngleXZ(mapFromSensorIDToKinkXZ);
setSensorIDToIncidenceAngleYZ(mapFromSensorIDToKinkYZ);
setBeamEnergy(beamEnergy);

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
		const float* statePosition = state.getPosition();
		const double* hitPosition = hit.getPosition();
		float residual[2];
		streamlog_out(DEBUG2) << "State position: " << statePosition[0]<<","<<statePosition[1]<<","<<statePosition[2]<< std::endl;
		streamlog_out(DEBUG2) << "Hit position: " << hitPosition[0]<<","<<hitPosition[1]<<","<<hitPosition[2]<< std::endl;

		typedef std::map<int ,AIDA::IProfile2D*  >::iterator it_type;
		for(it_type iterator = _mapFromSensorIDToHistogramX.begin(); iterator != _mapFromSensorIDToHistogramX.end(); iterator++) {
			if(iterator->first == state.getLocation()){
			residual[0]=statePosition[0]-hitPosition[0];
			streamlog_out(DEBUG2) << "Add residual X : " << residual[0]<< std::endl;

			_mapFromSensorIDToHistogramX[ state.getLocation() ]  -> fill( hitPosition[0], hitPosition[1], residual[0], 1 );
			break;
			}
		}
		for(it_type iterator = _mapFromSensorIDToHistogramY.begin(); iterator != _mapFromSensorIDToHistogramY.end(); iterator++) {
			if(iterator->first == state.getLocation()){

			residual[1]=statePosition[1]-hitPosition[1];
			streamlog_out(DEBUG2) << "Add residual Y : " << residual[1]<< std::endl;
			_mapFromSensorIDToHistogramY[ state.getLocation() ]  -> fill( hitPosition[0], hitPosition[1], residual[1], 1 );
			break;
			}
		}
	} 
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotResidualVsPosition------------------------------END"<< std::endl;
}

void EUTelTrackAnalysis::plotBeamEnergy(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotBeamEnergy------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	EUTelState state  = states.at(0);
	state.print();
	float omega = -1.0/state.getMomLocal().Mag();	
	_beamEnergy-> fill(-1.0/omega );
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotBeamEnergy------------------------------END"<< std::endl;
}
void EUTelTrackAnalysis::plotPValueVsBeamEnergy(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueVsBeamEnergy------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	EUTelState state  = states.at(0);
	state.print();
	float omega = -1.0/state.getMomLocal().Mag();	
	float pValue = calculatePValueForChi2(track);
	_pValueVsBeamEnergy->fill(-1.0/omega, pValue);

  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueVsBeamEnergy------------------------------END"<< std::endl;
}


void EUTelTrackAnalysis::plotIncidenceAngles(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotIncidenceAngles------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		TVectorD stateVec = state.getStateVec();
		float incidenceXZ = stateVec[1];
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
		TVectorD stateVec = state.getStateVec();
		float incidenceYZ = stateVec[2];
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
void EUTelTrackAnalysis::plotPValueWithIncidenceAngles(EUTelTrack track){
	streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotPValueWithIncidenceAngles------------------------------BEGIN"<< std::endl;
	float pValue = calculatePValueForChi2(track);
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		TVectorD stateVec = state.getStateVec();
		float incidenceXZ = stateVec[1];
		typedef std::map<int , AIDA::IProfile1D * >::iterator it_type;
		for(it_type iterator =_mapFromSensorIDToPValuesVsIncidenceXZ.begin(); iterator != _mapFromSensorIDToPValuesVsIncidenceXZ.end(); iterator++) {
			if(iterator->first == state.getLocation()){
			streamlog_out(DEBUG2) << "Add incidence XZ : " << incidenceXZ  << std::endl;
			_mapFromSensorIDToPValuesVsIncidenceXZ[ state.getLocation() ] ->fill(incidenceXZ, pValue);
			break;
			}
		}
	} 
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		TVectorD stateVec = state.getStateVec();
		float incidenceYZ = stateVec[2];
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
	float pValue = calculatePValueForChi2(track);
	std::vector<EUTelState> states = track.getStates();
	for(size_t i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();

		const float* statePosition = state.getPosition();
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
float EUTelTrackAnalysis::calculatePValueForChi2(EUTelTrack track){
//    boost::math::chi_squared mydist(track.getNdf());
//    float pValue = 1 - boost::math::cdf(mydist,track.getChi2());
    return 1.0;
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
