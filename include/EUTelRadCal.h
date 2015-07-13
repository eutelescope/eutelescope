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

	class  EUTelRadCal {
        public:
            struct Block{
                bool isSen;
                double totalRad;
                double weigMean;
                double weigVar;
                std::vector<double> startPos; 
                std::vector<double> endPos;
            };
    /// Find Radiation length with TGeo Volumes 
    /// This will deterimine the radiation length between the start and end point. 
    /// The radiation length of the track must be calculated in full and the variance determined by the Highland formula
    /// The variance is then split between each object (plane/air) using each objects radition length relative to the total
    /// The next step is translating the variance for each object to scattering points (See GBL scattering method). 
    ///
    /// This function will calculate the radiation length for the planes and all the material in front of it
    /// Calculation is done according to the eq. (27.23)
    /// pdg.lbl.gov/2006/reviews/passagerpp.pdf

    /**
     * \param [in] globalPosStart Begin trajectory here. 
     * \param [out] globalPosFinish This is the end of the trajectory 
     * \param [out] sensors links between ID->Rad length 
     * \para, [out] air links material after certain sensor ID->Rad length after. (Do not access last sensor of air! Does not exist)
     */



            /// Find Radiation length with TGeo Volumes 
            /// This will deterimine the radiation length between the start and end point. 
            /// The radiation length of the track must be calculated in full and the variance determined by the Highland formula
            /// The variance is then split between each object (plane/air) using each objects radition length relative to the total
            /// The next step is translating the variance for each object to scattering points (See GBL scattering method). 
            ///
            /// Planes side by side will produce a map which can only be applied to certain tracks. 
            /// This information is best associated to a track. Since each track has its own scattering variance.  
            ///
            /// To construct the radiation length maps the calcuation only begins when the first detector plane is impacted on the trajectory.
            /// So the global hit position passed to this function should be taken just before the first plane to include all of it's radiation length. 
            /// The final hit position will include all the radiation length bu default. Since TGeo will propagate through the final volume found.
            /// Note the the start and end position only defines the direction. The full detector system is propagated and returned inside the map.
            /// How this is distributed throughout the system is done in the GBL fitter.
            ///
            /// This function will calculate the radiation length for the planes and all the material in front of it
            /// Calculation is done according to the eq. (27.23)
            /// pdg.lbl.gov/2006/reviews/passagerpp.pdf

            /**
             * \param [in] globalPosStart Begin trajectory here. 
             * \param [out] globalPosFinish This is the end of the trajectory 
             * \param [out] sensors links between ID->Rad length 
             * \param [out] air links material after certain sensor ID->Rad length after. (Do not access last sensor of air! Does not exist)
             * \return perRad Percentage radiation length is returned for the full detector system.
             */

            float setHomoBlocks(TVector3& start , TVector3& end ,std::vector<int>& sen,std::map<int ,Block>& blocks);
            void setRad(EUTelTrack& track);

            void move(TGeoNode *node, double& rad, double& dist);

        private:
            void fillBlocks(geo::EUTelGeometryTelescopeGeoDescription& geo,double & rad, double& dist,std::map<int ,Block>& blocks);

        protected:
            std::vector<int> _planes;


	};

}
#endif	/* EUTELRADCAL */
