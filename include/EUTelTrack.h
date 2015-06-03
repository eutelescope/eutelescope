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
			//getters
            float getChi2() const ;
            float getNdf() const;
			float getTotalVariance() const { return _var;}
			unsigned int getNumberOfHitsOnTrack() const;
            //Must return reference to change the contents.
			std::vector<EUTelState>& getStates();
            std::vector<EUTelState> getStatesCopy() const;
            std::vector<double> getLCIOOutput();
			//setters
            void setState(EUTelState state);
            void setStates(std::vector<EUTelState> states);
			void setTotalVariance(double rad);
            void setChi2(float chi2);
            void setNdf(float nDF);
            void setTrackFromLCIOVec(std::vector<double> input);

			//print
			void print();
            //
            std::vector<EUTelState> _states;
            double _var;
            float _chi2;
            float _nDF;

  	private:
	};

}
#endif
