#ifndef EUTELRADCAL_H
#define	EUTELRADCAL_H

#include "EUTelUtility.h"
#include <iostream>
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGBLFitter.h"
#include "TMath.h"
#include "EUTelBlock.h"

#include <map>

namespace eutelescope {
  //! This class performs all radiation length calculations. 
  /*! This class will use EUTelTrack objects directly. 
   *  Radiation length is stored in terms of blocks. All sensors included in fit are assumed a thin scatterer.
   *  The radiation length in front of each sensor is assumed thick. Integration over the propagated trajectory for thick scatterers is needed.
   *  After this integration the weighted position and variance of the arclength is determined. This is used to model the thick scatterer by two scattering points.
   *  The "integration" is done by getting the needed radiation length and distance for each homogeneous volume.  
   *  These volumes are combined if dead material to describe the full scattering by only two scattering points.
   *
   */ 


	class  EUTelRadCal {
        public:
            void setMeanWeight(EUTelTrack & );
            void setVarWeight(EUTelTrack & );

            //! This will add the scattering information to each state on a track. This comes in the form of blocks. 
            /*!  
             *  @param [in] track
             */

            void setRad(EUTelTrack& track, int& mode);
            void setIncSenBlocks(EUTelTrack & track); 
            void getThicknessAndRad(EUTelTrack & track);
            void getRelativePosOfScatterers(EUTelState & state);
            void setPosVar(EUTelTrack & track); 
            void getVarForAllScatters(EUTelTrack & track );
            void getVarForSensorScatterersOnly(EUTelTrack & track );



	};
}
#endif	/* EUTELRADCAL */
