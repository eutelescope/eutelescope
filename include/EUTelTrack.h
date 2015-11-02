#ifndef EUTELTRACK_H
#define	EUTELTRACK_H

#include "EUTelUtility.h"
// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif
//lcio
#include "IMPL/TrackImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelState.h"

namespace eutelescope {

	class  EUTelTrack{
		public: 
			EUTelTrack();
			EUTelTrack( const EUTelTrack& track);
			EUTelTrack( const EUTelTrack& track,bool);
			//getters//TRACK PARAMETERS
            float getChi2() const ;
            float getNdf() const;
			double getRadPerTotal() const { return _radPerTotal;}
            float getQOverP() const { return _qOverP; }
            float getBeamEnergy() const { return -1.0/_qOverP; }
            std::vector<int>  getPlaIDs() const;
            std::vector<int>  getPlaIDDUTs() const; 

            //END OF TRACK PARAMETERS
            //Must return reference to change the contents.
			unsigned int getNumberOfHitsOnTrack() const;
			std::vector<EUTelState>& getStates();
            std::vector<EUTelState> getStatesCopy() const;
            std::vector<double> getLCIOOutput();
			//setters
            void setState(EUTelState state);
            void setStates(std::vector<EUTelState> states);
            std::vector<EUTelHit> getHitsCopy() const; 

			void setRadPerTotal(double radPer){this->_radPerTotal = radPer;}
            void setChi2(float chi2);
            void setNdf(float nDF);
            void setTrackFromLCIOVec(std::vector<double> input);
            void setQOverP(double qOverP ){ _qOverP = qOverP; }
			void setTrackUsingCorrection(TVectorD corrections);
			//print
			void print();
  	private:
            std::vector<EUTelState> _states;
            double _qOverP;
            double _radPerTotal;
            float _chi2;
            float _nDF;

	};

}
#endif
