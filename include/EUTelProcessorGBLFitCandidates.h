#ifndef EUTELESCOPEPROCESSORGBLFITCANDIDATES_H
#define	EUTELESCOPEPROCESSORGBLFITCANDIDATES_H

#ifdef USE_GBL

// C++
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>

// LCIO
#include <EVENT/LCCollection.h>
#include "lcio.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

// AIDA
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IProfile2D.h>
#endif // MARLIN_USE_AIDA

//EUTelescope
#include "EUTelUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackFitter.h"
#include "EUTelGBLFitter.h"
#include "EUTelExceptions.h"
#include "EUTelEventImpl.h"
#include "EUTelTrackStateImpl.h"
#include "EUTelTrackImpl.h"

using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

 class EUTelProcessorGBLFitCandidates : public Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorGBLFitCandidates);     // prevent users from making (default) copies of processors
        
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorGBLFitCandidates;
        }

        EUTelProcessorGBLFitCandidates();
        
    public:
        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init();

        /** Called for every run.
         */
        virtual void processRunHeader(LCRunHeader* run);

        /** Called for every event - the working horse.
         */
        virtual void processEvent(LCEvent * evt);

        virtual void check(LCEvent * evt);

        /** Called after data processing for clean up. **/
	
			  virtual void end();

    protected:

        // Statistics counters

        /** Number of events processed */
        int _nProcessedRuns;
        /** Number of runs processed */
        int _nProcessedEvents;

        /** Input TrackerHit collection name */
        string _trackCandidatesInputCollectionName;

        /** Output Tracks collection name */
        string _tracksOutputCollectionName;

				string _milleBinaryFilename;

        /** Track fitter */
        EUTelGBLFitter *_trackFitter;

				//This specifies what rotations etc will be taken into account. This directly effects the jacobian that transforms alignment parameters to changes in measured hits
				int _alignmentMode;

        /** Beam charge in [e] */
        double _beamQ;

				//Beam energy. 
				double _eBeam;

				//This is the maximum chi2 of a track that will be used in the millepede alignment fit
				double _maxChi2Cut;

        /** Outlier downweighting option */
        std::string _mEstimatorType;

				void CreateEUTrackandStates(TrackImpl* trackimpl, EUTelTrackImpl*);

};

    /** A global instance of the processor */
    EUTelProcessorGBLFitCandidates gEUTelProcessorGBLFitCandidates;

}


#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H */
