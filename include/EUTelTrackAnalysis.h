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
#include <boost/math/distributions/chi_squared.hpp> 

#include <map>

namespace eutelescope {

	class  EUTelTrackAnalysis {
        
	public:
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramX;
		std::map< int,  AIDA::IProfile2D* > _mapFromSensorKinksMap;
		std::map< int, AIDA::IHistogram2D* > _mapFromSensorIDHitMap;
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramY;
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToEfficiencyX;
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToEfficiencyY;
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDTo2DPValuesWithPosition;
		std::map< int,   AIDA::IHistogram1D *> _mapFromSensorIDToIncidenceXZ;
		std::map< int,  AIDA::IHistogram1D * > _mapFromSensorIDToIncidenceYZ;
		std::map< int,   AIDA::IProfile1D *> _mapFromSensorIDToPValuesVsIncidenceXZ;
		std::map< int,  AIDA::IProfile1D * > _mapFromSensorIDToPValuesVsIncidenceYZ;
		std::map< int,   AIDA::IHistogram1D *> _mapKinksX;
		std::map< int,  AIDA::IHistogram1D * > _mapKinksY;

    EUTelTrackAnalysis( std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY, std::map< int,  AIDA::IHistogram2D*> mapFromSensorIDHitMap, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToEfficiencyX, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToEfficiencyY, std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToGloIncXZ,std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToGloIncYZ,  std::map< int,   AIDA::IHistogram1D *> mapKinksX, std::map< int,   AIDA::IHistogram1D *> mapKinksY,		std::map< int,  AIDA::IProfile2D* > mapFromSensorKinksMap, AIDA::IHistogram1D * beamEnergy);

		template<typename T>
		std::string numberToString(T number);
		void plotResidualVsPosition(EUTelTrack);
        void plotKinksVsPosition(EUTelTrack track);
		void plotHitMap(EUTelTrack);
		void plotEfficiencyVsPosition(EUTelTrack,IntVec );
		void plotBeamEnergy(EUTelTrack);
		void plotIncidenceAngles(EUTelTrack track);
		void plotPValueWithPosition(EUTelTrack track);
		void plotPValueWithIncidenceAngles(EUTelTrack track);
		void plotPValueVsBeamEnergy(EUTelTrack track);
        void plotKinks(EUTelTrack & track);

		void setBeamEnergy(AIDA::IHistogram1D *  beamEnergy){ _beamEnergy = beamEnergy; }
		void setPValueBeamEnergy(AIDA::IProfile1D *  pValueVsBeamEnergy){ _pValueVsBeamEnergy = pValueVsBeamEnergy; }
		void setSensorIDTo2DHitMap(std::map< int,  AIDA::IHistogram2D*> mapFromSensorIDHitMap){_mapFromSensorIDHitMap=mapFromSensorIDHitMap;}
		void setSensorIDTo2DResidualHistogramX(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX){_mapFromSensorIDToHistogramX=mapFromSensorIDToHistogramX;}
		void setSensorIDTo2DResidualHistogramY(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY){_mapFromSensorIDToHistogramY=mapFromSensorIDToHistogramY;}
		void setSensorIDTo2DResidualEfficiencyX(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToEfficiencyX){_mapFromSensorIDToEfficiencyX=mapFromSensorIDToEfficiencyX;}
		void setSensorIDTo2DResidualEfficiencyY(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToEfficiencyY){_mapFromSensorIDToEfficiencyY=mapFromSensorIDToEfficiencyY;}
		void setSensorIDTo2DPValuesWithPosition(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDTo2DPValuesWithPosition){_mapFromSensorIDTo2DPValuesWithPosition=mapFromSensorIDTo2DPValuesWithPosition;}
		void setSensorIDToIncidenceAngleXZ( std::map< int,  AIDA::IHistogram1D * > mapFromSensorIDToGloIncXZ){_mapFromSensorIDToIncidenceXZ=mapFromSensorIDToGloIncXZ;}
		void setSensorIDToIncidenceAngleYZ( std::map< int,  AIDA::IHistogram1D * > mapFromSensorIDToGloIncYZ){_mapFromSensorIDToIncidenceYZ=mapFromSensorIDToGloIncYZ;}
		void setSensorIDToPValuesVsIncidenceAngleXZ( std::map< int,  AIDA::IProfile1D * > mapFromSensorIDToPValuesVsIncidenceXZ){_mapFromSensorIDToPValuesVsIncidenceXZ=mapFromSensorIDToPValuesVsIncidenceXZ;}
		void setSensorIDToPValuesVsIncidenceAngleYZ( std::map< int,  AIDA::IProfile1D * > mapFromSensorIDToPValuesVsIncidenceYZ){_mapFromSensorIDToPValuesVsIncidenceYZ=mapFromSensorIDToPValuesVsIncidenceYZ;}
//		void setHistName(std::string histName){  _histoInfoFileName = histName; }
        void setTotNum(EUTelTrack & track);
        void print();

        protected:
        std::map<int, float > _senResTotX;
        std::map<int, float > _senResTotY;
        std::map<int, float > _senResTotZ;
        std::map<int, int >  _hitNum;
		AIDA::IHistogram1D * _beamEnergy;
		AIDA::IProfile1D   * _pValueBeamEnergy;
		AIDA::IProfile1D * _pValueVsBeamEnergy;
		float calculatePValueForChi2(EUTelTrack track);
 //       std::string _histoInfoFileName;


	};

}
