#ifndef EUTELREADERGENERICLCIO_H
#define	EUTELREADERGENERICLCIO_H

#include "EUTelUtility.h"
// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif
#include "EUTelTrack.h"
// LCIO
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackImpl.h>
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/LCRelationImpl.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <IMPL/LCCollectionVec.h>
// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"











namespace eutelescope {

	class  EUTelReaderGenericLCIO{
		public: 
			EUTelReaderGenericLCIO();
            void getColVec( std::vector<EUTelTrack> tracks,LCEvent* evt,std::string colName );
            std::vector<EUTelTrack> getTracks( LCEvent* evt, std::string colName);

  	private:
	};

}
#endif
