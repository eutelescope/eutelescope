/*
 * EUTelTrackSelection.h 
 * 
 * Created on: April 19th 2015 
 *     author:Alexander Morton 
 * 
 *          This class creates a track based selection object which will take a vector of tracks and output a subset of these based on an a combined series of cuts. 
 *          The cuts are track based, by which all the cuts are on track wide parameters. The do not vary from plane to plane. 
 *          Further development: Any new track wide selection should be placed here and implemented in any new EUTelescope processor.
 *          Plane to plane variations of course include the state fit parameters which must be investigated to gain optimum alignment. 
 *          One such parameter is the non-gaussian tail in the measurement of kink angle and the final fit. 
 * 
 * 
 *
 */


#include "EUTelUtility.h"
#include <fstream>      
#include <iostream>
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelHistogramManager.h"
#include "EUTelMillepede.h"

//system includes
#include <string>
#include <vector>
#include <memory>
//#include <iostream>
#include <cmath>

namespace eutelescope {

	class  EUTelTrackSelection {
        
	public:
    EUTelTrackSelection();
    bool removeTracksWithHitsOnPlane(EUTelTrack track,std::vector<int> sensors); 
	};

}
