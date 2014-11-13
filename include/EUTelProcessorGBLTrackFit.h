#ifndef EUTELESCOPEPROCESSORGBLTRACKFIT_H
#define	EUTELESCOPEPROCESSORGBLTRACKFIT_H

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

	class EUTelProcessorGBLTrackFit : public Processor {

		private:
			DISALLOW_COPY_AND_ASSIGN(EUTelProcessorGBLTrackFit)      // prevent users from making (default) copies of processors
						
		public:

			virtual Processor* newProcessor() {
					return new EUTelProcessorGBLTrackFit;
			}

			EUTelProcessorGBLTrackFit();

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

			bool _first_time;
			/** Number of events processed */
			int _nProcessedRuns;
			/** Number of runs processed */
			int _nProcessedEvents;

			/** Beam charge in [e] */
			double _beamQ;

			//Beam energy. 
			double _eBeam;

			//This is the maximum chi2 of a track that will be used in the millepede alignment fit
			double _maxChi2Cut;

			std::vector<float> _chi2NdfVec;

			/** Input TrackerHit collection name */
			string _trackCandidatesInputCollectionName;

			/** Output Tracks collection name */
			string _tracksOutputCollectionName;

			/** Outlier downweighting option */
			std::string _mEstimatorType;

        /** Histogram info file name */
			std::string _histoInfoFileName;

			/** x Resolution of planes in PlaneIds */
			FloatVec _SteeringxResolutions;
 
			/** y Resolution of planes in PlaneIds */
			FloatVec _SteeringyResolutions;

			/** Track fitter */
			EUTelGBLFitter *_trackFitter;
			//Function defined now for the processor////////////////////////////
			void outputLCIO(LCEvent* evt, std::vector< EUTelTrack >& tracks);

			void bookHistograms();

			void plotResidual(map< int, map<float, float > >  & sensorResidual, map< int, map<float, float > >  & sensorResidualError, bool &first_time);
				
//TO DO: Fix all this histogramming stuff.
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
	std::map< int, AIDA::IHistogram1D* > _mapSensorIDToHistogramCorrection0;
	std::map< int, AIDA::IHistogram1D* > _mapSensorIDToHistogramCorrection1;
	std::map< int, AIDA::IHistogram1D* > _mapSensorIDToHistogramCorrection2;
	std::map< int, AIDA::IHistogram1D* > _mapSensorIDToHistogramCorrection3;
	std::map< int, AIDA::IHistogram1D* > _mapSensorIDToHistogramCorrection4;
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

						static string _residGblFitHistNameX0p;
						static string _residGblFitHistNameX1p;
						static string _residGblFitHistNameX2p;
						static string _residGblFitHistNameX3p;
						static string _residGblFitHistNameX4p;
						static string _residGblFitHistNameX5p;
						static string _residGblFitHistNameY0p;
						static string _residGblFitHistNameY1p;
						static string _residGblFitHistNameY2p;
						static string _residGblFitHistNameY3p;
						static string _residGblFitHistNameY4p;
						static string _residGblFitHistNameY5p;

        };

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

	};

/** A global instance of the processor */
EUTelProcessorGBLTrackFit gEUTelProcessorGBLTrackFit;

}

#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORGBLFITCANDIDATES_H */
