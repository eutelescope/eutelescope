#include "EUTelTrackAnalysis.h"
using namespace eutelescope;
EUTelTrackAnalysis::EUTelTrackAnalysis(map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX, map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY){
setSensorIDTo2DResidualHistogramX(mapFromSensorIDToHistogramX);
setSensorIDTo2DResidualHistogramY(mapFromSensorIDToHistogramY);
} 

void EUTelTrackAnalysis::plotResidualVsPosition(EUTelTrack track){
  streamlog_out(DEBUG2) << " EUTelTrackAnalysis::plotResidualVsPosition------------------------------BEGIN"<< std::endl;
	std::vector<EUTelState> states = track.getStates();
	for(int i=0; i<states.size();++i){
		EUTelState state  = states.at(i);
		state.print();
		EVENT::TrackerHit* hit = state.getHit();	
		float* statePosition = state.getPosition();
		const double* hitPosition = hit->getPosition();
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


