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

namespace eutelescope {

	class  EUTelHit{
		public: 
			EUTelHit();
			EUTelHit(EUTelHit* hit);
            void setPosition(const double * position);
            void setID(int id);
            void setTrackFromLCIOVec(std::vector<double> input);

            //get
            int getID() const;
            const double* getPosition() const; 

		    double _position[3];	
            int _id; //This is used to keep a track of all the hits for track removal.
            std::vector<double> getLCIOOutput(); 
  	private:
	};

}
#endif
