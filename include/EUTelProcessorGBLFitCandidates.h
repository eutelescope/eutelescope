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
//#include "LCIOSTLTypes.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

//GBL
#include "include/GblTrajectory.h"

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
#include "EUTelHistogramManager.h"

using namespace lcio;
using namespace marlin;
using namespace std;
namespace eutelescope {

 class EUTelProcessorGBLFitCandidates : public Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorGBLFitCandidates)      // prevent users from making (default) copies of processors
        
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

        /** Track fitter */
        EUTelGBLFitter *_trackFitter;


        /** Beam charge in [e] */
        double _beamQ;

				//Beam energy. 
				double _eBeam;

				//This is the maximum chi2 of a track that will be used in the millepede alignment fit
				double _maxChi2Cut;
				std::vector<float> _chi2NdfVec;

        /** Outlier downweighting option */
        std::string _mEstimatorType;

        /** Histogram info file name */
			std::string _histoInfoFileName;

        /** x Resolution of planes in PlaneIds */
        FloatVec _SteeringxResolutions;
 
        /** y Resolution of planes in PlaneIds */
        FloatVec _SteeringyResolutions;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        /** AIDA histogram map
         *  Instead of putting several pointers to AIDA histograms as
         *  class members, histograms are booked in the init() method and
         *  their pointers are inserted into this map keyed by their
         *  names.
         *  The histogram filling can proceed recalling an object through
         *  its name
         */
	std::map< std::string, AIDA::IHistogram1D* > _aidaHistoMap1D;

        /** Names of histograms */
        struct _histName {
						static string _chi2CandidateHistName;
						static string  _fitsuccessHistName;
						static string _residGblFitHistNameX0;
						static string _residGblFitHistNameX1;
						static string _residGblFitHistNameX2;
						static string _residGblFitHistNameX3;
						static string _residGblFitHistNameX4;
						static string _residGblFitHistNameX5;
						static string _residGblFitHistNameY0;
						static string _residGblFitHistNameY1;
						static string _residGblFitHistNameY2;
						static string _residGblFitHistNameY3;
						static string _residGblFitHistNameY4;
						static string _residGblFitHistNameY5;

        };

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

				void CreateEUTrackandStates(TrackImpl* trackimpl, EUTelTrackImpl*);

				void outputLCIO(LCEvent* evt, std::vector< EUTelTrack >& tracks);

				void bookHistograms();

				void plotResidual( map< int, map< float, float > >  & SensorResidualError, bool & first_time);
				bool _first_time=true;
				


};

    /** A global instance of the processor */
    EUTelProcessorGBLFitCandidates gEUTelProcessorGBLFitCandidates;

}


#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORGBLFITCANDIDATES_H */
