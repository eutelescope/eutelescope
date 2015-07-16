#ifndef EUTELRADCAL_H
#define	EUTELRADCAL_H

#include "EUTelUtility.h"
#include <iostream>
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGBLFitter.h"
#include "TMath.h"

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
            //! Struct object
            /*! This stores the basic information for each block of radiation length.
             *  Each block contains the information to model the scattering for the sensor and the material infront of it.  
             *
             */

            struct Block{
                bool isSen;
                double senRadPer;
                double medRadPer;
                double weigMean;
                double weigVar;
                std::vector<std::pair<double,double> > thicknessAndRad; 

            };
            void setMeanWeight(Block & block);
            void setVarWeight(Block & block);

            //! The function will take a track and return the radiation length blocks 
            /*! The tracks are needed to determine which planes should be excluded. 
             *  @param [in] track
             */

            void getRad(EUTelTrack& track);
            void setIncSenBlocks(EUTelTrack const & track, std::map<int, Block>& blocks); 
            void getThicknessAndRad(EUTelTrack const & track, std::map<int, Block> & blocks);
            void getScatParam(std::map<int, Block> & blocks);



            void setBlockWithoutTrackInfo(EUTelTrack const  & track, std::map<int, Block>& blocks);


        protected:


	};

}
#endif	/* EUTELRADCAL */
