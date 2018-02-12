/*
 * Rewritten by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef EUTelTrackingHelixTrackSearch_h
#define EUTelTrackingHelixTrackSearch_h 1

// C++
#include <string>

// LCIO
#include "lcio.h"
#include "IMPL/TrackerHitImpl.h"

// MARLIN
#include "marlin/Processor.h"

// AIDA
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#endif

// EUTELESCOPE
#include "EUTelTrackFitter.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelUtility.h"
#include "EUTelTrackImpl.h"

using namespace lcio;
using namespace marlin;
using namespace std;

class EUTelTrackImpl;

namespace eutelescope
{

    class EUTelProcessorTrackingHelixTrackSearch : public Processor
    {
	private:

	    DISALLOW_COPY_AND_ASSIGN ( EUTelProcessorTrackingHelixTrackSearch )

	    double _residualsRMax;

	    double _eBeam;

	    double _qBeam;

	    double _eBeamUncertatinty;

	    EVENT::FloatVec _beamSpread;

	    int _maxMissingHitsPerTrackCand;

	    int _maxNTracks;

	    string _tgeoFileName;

	public:

	    virtual Processor* newProcessor ( )
	    {
		return new EUTelProcessorTrackingHelixTrackSearch;
	    }

	    EUTelProcessorTrackingHelixTrackSearch ( );

	    #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

	    map < string, AIDA::IHistogram1D* > _aidaHistoMap1D;

	    #endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

	    virtual void check ( LCEvent * evt );

	    virtual void end ( );

	    virtual void init ( );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void processRunHeader ( LCRunHeader* run );

	    void addTrackCandidateToCollection ( LCEvent* evt, std::vector < IMPL::TrackImpl* > & );

	    void bookHistograms ( );

	    void FillHits ( LCCollection*, EVENT::TrackerHitVec& ) const;

	protected:

	    EUTelTrackFitter* _trackFitter;

	    int _nProcessedRuns;

	    int _nProcessedEvents;

	    string _hitInputCollectionName;

	    string _trackCandidateHitsOutputCollectionName;

	    string _fitpointcollectionname;

    };

    EUTelProcessorTrackingHelixTrackSearch gEUTelProcessorTrackingHelixTrackSearch;

} // eutelescope

#endif // EUTelTrackingExhaustiveTrackSearch_h
