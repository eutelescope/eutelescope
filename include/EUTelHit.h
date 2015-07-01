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
            /// Used to set the covariance from processor input.
            void setCov(const std::vector<double> &);
            void setCov(const std::vector<float> &);
            /// Set the covariance matrix from one state to another
            void setCov(TMatrixD );
            void setPulse( EVENT::LCObjectVec&);

            //get
            void getCov(double (&)[4]) const;
            TMatrixD getCov() const {return _cov;};
            TVector3 getPositionGlobal() const; 
			int	getLocation() const;
            //HIT PARAMETERS
            const double* getPosition() const; 
            int getID() const;
            ///This is the tracker pulse collection which you use to link the clustering information.
            EVENT::LCObjectVec getPulse() const; 
            //END HIT PARAMETERS
            //print
            void print();

            std::vector<double> getLCIOOutput(); 
			bool operator==(const EUTelHit compareHit ) const;

  	protected:
		    double _position[3];	
            int _location; 
            int _locationKnown;
            int _id; //This is used to keep a track of all the hits for track removal.
			TMatrixD _cov;
            EVENT::LCObjectVec _pulse;


	};

}
#endif
