#include "EUTelUtility.h"
#include <fstream>      // std::ifstream, std::ofstream
#include <iostream>
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelHistogramManager.h"
#include "EUTelMillepede.h"
// AIDA
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IProfile2D.h>
#endif // MARLIN_USE_AIDA
#include "EUTelHit.h"
//#include <boost/math/distributions/chi_squared.hpp> 

#include <map>

namespace eutelescope {

	class  EUTelTrackAnalysis {
        
	public:
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramX;
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramY;
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDTo2DPValuesWithPosition;
		std::map< int,   AIDA::IHistogram1D *> _mapFromSensorIDToIncidenceXZ;
		std::map< int,  AIDA::IHistogram1D * > _mapFromSensorIDToIncidenceYZ;
		std::map< int,   AIDA::IProfile1D *> _mapFromSensorIDToPValuesVsIncidenceXZ;
		std::map< int,  AIDA::IProfile1D * > _mapFromSensorIDToPValuesVsIncidenceYZ;

    EUTelTrackAnalysis(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY, std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkXZ,std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkYZ, AIDA::IHistogram1D * beamEnergy );

		template<typename T>
		std::string numberToString(T number);
		void plotResidualVsPosition(EUTelTrack);
		void plotBeamEnergy(EUTelTrack);
		void plotIncidenceAngles(EUTelTrack track);
		void plotPValueWithPosition(EUTelTrack track);
		void plotPValueWithIncidenceAngles(EUTelTrack track);
		void plotPValueVsBeamEnergy(EUTelTrack track);
		void setBeamEnergy(AIDA::IHistogram1D *  beamEnergy){ _beamEnergy = beamEnergy; }
		void setPValueBeamEnergy(AIDA::IProfile1D *  pValueVsBeamEnergy){ _pValueVsBeamEnergy = pValueVsBeamEnergy; }
		void setSensorIDTo2DResidualHistogramX(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX){_mapFromSensorIDToHistogramX=mapFromSensorIDToHistogramX;}
		void setSensorIDTo2DResidualHistogramY(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY){_mapFromSensorIDToHistogramY=mapFromSensorIDToHistogramY;}
		void setSensorIDTo2DPValuesWithPosition(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDTo2DPValuesWithPosition){_mapFromSensorIDTo2DPValuesWithPosition=mapFromSensorIDTo2DPValuesWithPosition;}
		void setSensorIDToIncidenceAngleXZ( std::map< int,  AIDA::IHistogram1D * > mapFromSensorIDToKinkXZ){_mapFromSensorIDToIncidenceXZ=mapFromSensorIDToKinkXZ;}
		void setSensorIDToIncidenceAngleYZ( std::map< int,  AIDA::IHistogram1D * > mapFromSensorIDToKinkYZ){_mapFromSensorIDToIncidenceYZ=mapFromSensorIDToKinkYZ;}
		void setSensorIDToPValuesVsIncidenceAngleXZ( std::map< int,  AIDA::IProfile1D * > mapFromSensorIDToPValuesVsIncidenceXZ){_mapFromSensorIDToPValuesVsIncidenceXZ=mapFromSensorIDToPValuesVsIncidenceXZ;}
		void setSensorIDToPValuesVsIncidenceAngleYZ( std::map< int,  AIDA::IProfile1D * > mapFromSensorIDToPValuesVsIncidenceYZ){_mapFromSensorIDToPValuesVsIncidenceYZ=mapFromSensorIDToPValuesVsIncidenceYZ;}
//		void setHistName(std::string histName){  _histoInfoFileName = histName; }

		AIDA::IHistogram1D * _beamEnergy;
		AIDA::IProfile1D   * _pValueBeamEnergy;
		AIDA::IProfile1D * _pValueVsBeamEnergy;
		float calculatePValueForChi2(EUTelTrack track);
 //       std::string _histoInfoFileName;


	};

}
