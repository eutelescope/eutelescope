#include "EUTelUtility.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelHistogramManager.h"
#include "boost/math/distributions/chi_squared.hpp" //TO DO: Does cmake install this automatically?
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

#include <map>

namespace eutelescope {

	class  EUTelTrackAnalysis {
        
    public:

    EUTelTrackAnalysis(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX, std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY, std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkXZ,std::map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkYZ);

		void plotResidualVsPosition(EUTelTrack);
		void plotIncidenceAngles(EUTelTrack track);
		void setSensorIDTo2DResidualHistogramX(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX){_mapFromSensorIDToHistogramX=mapFromSensorIDToHistogramX;}
		void setSensorIDTo2DResidualHistogramY(std::map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY){_mapFromSensorIDToHistogramY=mapFromSensorIDToHistogramY;}
		void setSensorIDToIncidenceAngleXZ( std::map< int,  AIDA::IHistogram1D * > mapFromSensorIDToKinkXZ){_mapFromSensorIDToIncidenceXZ=mapFromSensorIDToKinkXZ;}
		void setSensorIDToIncidenceAngleYZ( std::map< int,  AIDA::IHistogram1D * > mapFromSensorIDToKinkYZ){_mapFromSensorIDToIncidenceYZ=mapFromSensorIDToKinkYZ;}

		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramX;
		std::map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramY;
		std::map< int,   AIDA::IHistogram1D *> _mapFromSensorIDToIncidenceXZ;
		std::map< int,  AIDA::IHistogram1D * > _mapFromSensorIDToIncidenceYZ;


	};

}
