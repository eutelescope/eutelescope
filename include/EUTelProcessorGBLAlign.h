#ifndef EUTELPROCESSORGBLALIGN_H
#define	EUTELPROCESSORGBLALIGN_H


// C++
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <iostream>

// LCIO
#include <EVENT/LCCollection.h>
#include "lcio.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

//EUTelescope
#include "EUTelUtility.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelExceptions.h"
#include "EUTelGBLFitter.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelEventImpl.h"
#include "EUTelMillepede.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelReaderGenericLCIO.h"

namespace eutelescope {

 class  EUTelProcessorGBLAlign : public marlin::Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorGBLAlign)     // prevent users from making (default) copies of processors
        
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorGBLAlign;
        }

        EUTelProcessorGBLAlign();
        
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

        /** Called after data processing for clean up. **/
	
			  virtual void end();
				void printPointsInformation(std::vector<gbl::GblPoint>& pointList);
				double printSize(const std::string& address);

		protected: 

				std::string _milleBinaryFilename;
				std::string _milleSteeringFilename;
				std::string _milleResultFileName;
				std::string _gear_aligned_file;

        /** Number of events processed */
        int _nProcessedRuns;
        /** Number of runs processed */
        int _nProcessedEvents;
				int _totalTrackCount;

				int _alignmentMode;

        /** Beam charge in [e] */
        double _beamQ;
        //Beam energy. 
        double _eBeam;
        bool _createBinary;
        int _incMed;
        /** Outlier downweighting option */
        std::string _mEstimatorType;

        /** Track fitter */
        EUTelGBLFitter *_trackFitter;

        /** Input TrackerHit collection name */
        std::string _trackCandidatesInputCollectionName;

        /** Output Tracks collection name */
        std::string _tracksOutputCollectionName;

        /** Allows user-added commands in the pede steering file */
				lcio::StringVec _pedeSteerAddCmds;

				EUTelMillepede* _Mille;

        /** Alignment constants file name */
				std::string _alignmentConstantLCIOFile;

        /** Alignment X shift plane ids to be fixed */
	lcio::IntVec _fixedAlignmentXShfitPlaneIds;
        
        /** Alignment Y shift plane ids to be fixed */
	lcio::IntVec _fixedAlignmentYShfitPlaneIds;
        
        /** Alignment Z shift plane ids to be fixed */
	lcio::IntVec _fixedAlignmentZShfitPlaneIds;
        
        /** Alignment X rotation plane ids to be fixed */
	lcio::IntVec _fixedAlignmentXRotationPlaneIds;
        
        /** Alignment Y rotation plane ids to be fixed */
	lcio::IntVec _fixedAlignmentYRotationPlaneIds;
        
        /** Alignment Z rotation plane ids to be fixed */
	lcio::IntVec _fixedAlignmentZRotationPlaneIds;
        
        /** Planes ids to be excluded from refit */
	lcio::IntVec _excludePlanesFromFit;

        /** Alignment plane ids*/
	lcio::IntVec _alignmentPlaneIds;

        /** x Resolution of planes in PlaneIds */
        FloatVec _SteeringxResolutions;
 
        /** y Resolution of planes in PlaneIds */
        FloatVec _SteeringyResolutions;

	IntVec _excludePlanes;
};

    /** A global instance of the processor */
    EUTelProcessorGBLAlign gEUTelProcessorGBLAlign;
}

#endif	/* EUTELPROCESSORMILLEALIGN_H */
