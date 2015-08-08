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
#include "EUTelHistogramManager.h"
#include "EUTelReaderGenericLCIO.h"

namespace eutelescope {

	class EUTelProcessorGBLTrackFit : public marlin::Processor {

		private:
			DISALLOW_COPY_AND_ASSIGN(EUTelProcessorGBLTrackFit)
						
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

			/** Called after data processing for clean up. **/

			virtual void end();

    protected:

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
			//Pointer to access millepede object..
			EUTelMillepede* _Mille;


			/** Input TrackerHit collection name */
			std::string _trackCandidatesInputCollectionName;

			/** Output Tracks collection name */
			std::string _tracksOutputCollectionName;

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

			void plotResidual(std::map< int, std::map<float, float > >  & sensorResidual, std::map< int, std::map<float, float > >  & sensorResidualError);
				
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
						static std::string _chi2CandidateHistName;
						static std::string  _fitsuccessHistName;
						static std::string _residGblFitHistNameX0;
						static std::string _residGblFitHistNameX1;
						static std::string _residGblFitHistNameX2;
						static std::string _residGblFitHistNameX3;
						static std::string _residGblFitHistNameX4;
						static std::string _residGblFitHistNameX5;
	                                        static std::string _residGblFitHistNameXDut1;
	  	                                static std::string _residGblFitHistNameXDut2;
	  					static std::string _residGblFitHistNameY0;
						static std::string _residGblFitHistNameY1;
						static std::string _residGblFitHistNameY2;
						static std::string _residGblFitHistNameY3;
						static std::string _residGblFitHistNameY4;
						static std::string _residGblFitHistNameY5;
	                                        static std::string _residGblFitHistNameYDut1;
	                                        static std::string _residGblFitHistNameYDut2;

						static std::string _residGblFitHistNameX0p;
						static std::string _residGblFitHistNameX1p;
						static std::string _residGblFitHistNameX2p;
						static std::string _residGblFitHistNameX3p;
						static std::string _residGblFitHistNameX4p;
						static std::string _residGblFitHistNameX5p;
	                                        static std::string _residGblFitHistNameXDut1p;
	                                        static std::string _residGblFitHistNameXDut2p;
						static std::string _residGblFitHistNameY0p;
						static std::string _residGblFitHistNameY1p;
						static std::string _residGblFitHistNameY2p;
						static std::string _residGblFitHistNameY3p;
						static std::string _residGblFitHistNameY4p;
						static std::string _residGblFitHistNameY5p;
	                                        static std::string _residGblFitHistNameYDut1p;
	                                        static std::string _residGblFitHistNameYDut2p;

        };

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

	};

/** A global instance of the processor */
EUTelProcessorGBLTrackFit gEUTelProcessorGBLTrackFit;

}

#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORGBLFITCANDIDATES_H */
