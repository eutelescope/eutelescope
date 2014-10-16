
#ifndef EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H
#define	EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H


// C++
#include <string>

// LCIO
#include "lcio.h"

#include "marlin/Processor.h"

#include "IMPL/TrackerHitImpl.h"

#ifdef USE_GBL

// GBL
#include "include/MilleBinary.h"

// EUTELESCOPE
#include "EUTelTrackFitter.h"
#include "EUTelGBLFitter.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelUtility.h"

namespace eutelescope {

    /** @class EUTelProcessorAlignmentGBLMille 
     * 
     *  @see EUTelTrackFitter
     */
  class EUTelProcessorAlignmentGBLMille : public marlin::Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorAlignmentGBLMille)     // prevent users from making (default) copies of processors
        
    public:

    virtual marlin::Processor* newProcessor() {
            return new EUTelProcessorAlignmentGBLMille;
        }

        EUTelProcessorAlignmentGBLMille();
        
    public:
        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init();

        /** Called for every run.
         */
        virtual void processRunHeader(lcio::LCRunHeader* run);

        /** Called for every event - the working horse.
         */
        virtual void processEvent(lcio::LCEvent * evt);

        virtual void check(lcio::LCEvent * evt);

        /** Called after data processing for clean up.
         */
        virtual void end();

    private:
        /** Prepare MILLEPEDE steering file */
        void writeMilleSteeringFile();

        /** Rename MILLEPEDE result file */
        void moveMilleResultFile( const std::string&, const std::string& );
        
        /** Parse MILLEPEDE result file
         * and create output collection
         */
        bool parseMilleOutput( const std::string& );
        
        /** Run pede executable to solve alignment problem */
        void runPede();

        /** Generate MILLEPEDE parameter labels */
        void fillMilleParametersLabels();
       
    public:
        
        // Processor parameters

    public:

        // Necessary parameters

        /** Beam charge in [e] */
        double _qBeam;

        /** Beam energy in [GeV] */
        double _eBeam;


        // Optional parameters

        /** Alignment mode */
        int _alignmentMode;
      
        /** GBL M-estimator option */
	std::string _mEstimatorType;
        
        /** Parameter ids */
	std::map<int,int> _xShiftsMap;
        
        /** Parameter ids */
	std::map<int,int> _yShiftsMap;
        
        /** Parameter ids */
	std::map<int,int> _zShiftsMap;
        
        /** Parameter ids */
	std::map<int,int> _xRotationsMap;
        
        /** Parameter ids */
	std::map<int,int> _yRotationsMap;
        
        /** Parameter ids */
	std::map<int,int> _zRotationsMap;
        
        /** Mille binary filename */
	std::string _milleBinaryFilename;

        /** Mille steering filename */
	std::string _milleSteeringFilename;

        /** Mille result filename */
	std::string _milleResultFileName;
        
        /** GEAR new filename */
	std::string _gear_aligned_file;


        /** Allows user-added commands in the pede steering file */
	lcio::StringVec _pedeSteerAddCmds;
        
        /** Alignment plane ids*/
	lcio::IntVec _alignmentPlaneIds;
        
        /** plane ids*/
	lcio::IntVec _planeIds;
 
        /** x Resolution of planes in AlignmentPlaneIds */
	lcio::FloatVec _SteeringxResolutions;
 
        /** y Resolution of planes in AlignmentPlaneIds */
	lcio::FloatVec _SteeringyResolutions;

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
	lcio::FloatVec _excludePlanes;
        
        /** Automatic pede run flag*/
        bool _runPede;
        
        /** Alignment constants file name */
	std::string _alignmentConstantLCIOFile;
        
        /** Maximum value of track chi2 for millipede */
        double _maxMilleChi2Cut;

        /** TGeo geometry file name */
	std::string _tgeoFileName;
        
        /** Histogram info file name */
	std::string _histoInfoFileName;

    protected:

        // Input/Output collections of the processor

        /** Input TrackerHit collection name */
	std::string _trackCandidatesInputCollectionName;


    protected:

        /** Track fitter */
        EUTelTrackFitter *_trackFitter;

        /** Mille */
        gbl::MilleBinary * _milleGBL;

    private:

        struct SeedAlignmentConstants {
            SeedAlignmentConstants() : _xResiduals(), _nxResiduals(), _yResiduals(), _nyResiduals() {};
            std::map< int, double > _xResiduals;   //! sum all x residuals for given plane id
            std::map< int, int >    _nxResiduals;  //! number of residuals used to calculate mean for given plane id
            std::map< int, double > _yResiduals;   //! sum all y residuals for given plane id
            std::map< int, int >    _nyResiduals;  //! number of residuals used to calculate mean for given plane id
        };

        /** Initial alignment constants */
        SeedAlignmentConstants _seedAlignmentConstants;


    private:


    protected:

        // Statistics counters

        /** Number of events processed */
        int _nProcessedRuns;
        /** Number of runs processed */
        int _nProcessedEvents;

    private:
        bool _flag_nohistos;

    public:


    };

    /** A global instance of the processor */
    EUTelProcessorAlignmentGBLMille gEUTelProcessorAlignmentGBLMille;

} // eutelescope

#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H */

