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


using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

	class  EUTelTrack : public IMPL::TrackImpl{
		public: 
			EUTelTrack();
			EUTelTrack( const EUTelTrack& track);
			EUTelTrack( const EUTelTrack& track,bool);
			//getters
			int getNumberOfHitsOnTrack() const;
			std::vector<EUTelState> getStates();
			std::vector<EUTelState*> getStatesPointers();
			//print
			void print();

  	private:
	};

}
#endif
