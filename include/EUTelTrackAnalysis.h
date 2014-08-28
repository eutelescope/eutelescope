#include "EUTelUtility.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelHistogramManager.h"

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

using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

	class  EUTelTrackAnalysis {

  	private:
        
    public:

    EUTelTrackAnalysis(map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX, map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY, map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkXZ,map< int,   AIDA::IHistogram1D *> mapFromSensorIDToKinkYZ);

		void plotResidualVsPosition(EUTelTrack);
		void plotIncidenceAngles(EUTelTrack track);
		void setSensorIDTo2DResidualHistogramX(map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramX){_mapFromSensorIDToHistogramX=mapFromSensorIDToHistogramX;}
		void setSensorIDTo2DResidualHistogramY(map< int,  AIDA::IProfile2D*> mapFromSensorIDToHistogramY){_mapFromSensorIDToHistogramY=mapFromSensorIDToHistogramY;}
		void setSensorIDToIncidenceAngleXZ( map< int,  AIDA::IHistogram1D * > mapFromSensorIDToKinkXZ){_mapFromSensorIDToIncidenceXZ=mapFromSensorIDToKinkXZ;}
		void setSensorIDToIncidenceAngleYZ( map< int,  AIDA::IHistogram1D * > mapFromSensorIDToKinkYZ){_mapFromSensorIDToIncidenceYZ=mapFromSensorIDToKinkYZ;}

		map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramX;
		map< int, AIDA::IProfile2D* > _mapFromSensorIDToHistogramY;
		map< int,   AIDA::IHistogram1D *> _mapFromSensorIDToIncidenceXZ;
		map< int,  AIDA::IHistogram1D * > _mapFromSensorIDToIncidenceYZ;


	};

}
