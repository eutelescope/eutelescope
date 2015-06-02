#ifndef EUTELHIT_H
#define	EUTELHIT_H

#include "EUTelUtility.h"
// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif
#include "EUTelGeometryTelescopeGeoDescription.h"

namespace eutelescope {

	class  EUTelHit{
		public: 
			EUTelHit();
			EUTelHit(EUTelHit* hit);
			EUTelHit(EVENT::TrackerHit* hit);
            void setPosition(const double * position);
            void setID(int id);
            void setTrackFromLCIOVec(std::vector<double> input);
			void setLocation(int location);

            //get
            int getID() const;
            const double* getPosition() const; 
            TVector3 getPositionGlobal() const; 
			int	getLocation() const;


            //print
            void print();


		    double _position[3];	
            int _location; 
            int _locationKnown;
            int _id; //This is used to keep a track of all the hits for track removal.
            std::vector<double> getLCIOOutput(); 
			bool operator==(const EUTelHit compareHit ) const;

  	private:
	};

}
#endif
