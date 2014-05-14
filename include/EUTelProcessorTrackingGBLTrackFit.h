
#ifndef EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H
#define	EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H


// C++
#include <string>

// LCIO
#include "lcio.h"

#include "marlin/Processor.h"

#include "IMPL/TrackerHitImpl.h"

// AIDA
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#endif

#ifdef USE_GBL

// GBL
#include "include/MilleBinary.h"

// EUTELESCOPE
#include "EUTelTrackFitter.h"
#include "EUTelGBLFitter.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelUtility.h"

namespace eutelescope {

    /** @class EUTelProcessorTrackingGBLTrackFit 
     * 
     *  @see EUTelTrackFitter
     */
  class EUTelProcessorTrackingGBLTrackFit : public marlin::Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackingGBLTrackFit)     // prevent users from making (default) copies of processors
        
    public:

    virtual marlin::Processor* newProcessor() {
            return new EUTelProcessorTrackingGBLTrackFit;
        }

        EUTelProcessorTrackingGBLTrackFit();
        
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
        /** Histogram booking */
        void bookHistograms();

        
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
	lcio::IntVec _excludePlanesFromFit;
        
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
	std::string _trackCandidateHitsInputCollectionName;

        /** Output Tracks collection name */
	std::string _tracksOutputCollectionName;

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

    public:
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
        std::map< std::string, AIDA::IHistogram2D* > _aidaHistoMap2D;
        std::map< std::string, AIDA::IProfile2D* >   _aidaProfileMap2D;

        /** Names of histograms */
        struct _histName {
            static std::string _orchi2GblFitHistName;
            static std::string _orchi2ndfGblFitHistName;
            static std::string _orprobGblFitHistName;
            static std::string _chi2GblFitHistName;
            static std::string _chi2ndfGblFitHistName;
            static std::string _probGblFitHistName;
            static std::string _momentumGblFitHistName;
            static std::string _residGblFitHistName;
            static std::string _normResidGblFitHistName;
            static std::string _residGblFitHistNameX;
            static std::string _residGblFitHistNameY;
            static std::string _resid2DGblFitHistNameX;
            static std::string _resid2DGblFitHistNameY;
            static std::string _normResidGblFitHistNameX;
            static std::string _normResidGblFitHistNameY;
            static std::string _resid2DGblFitHistNameXvsX;
            static std::string _resid2DGblFitHistNameXvsY;
            static std::string _resid2DGblFitHistNameYvsX;
            static std::string _resid2DGblFitHistNameYvsY;
            static std::string _normResid2DGblFitHistNameXvsX;
            static std::string _normResid2DGblFitHistNameXvsY;
            static std::string _normResid2DGblFitHistNameYvsX;
            static std::string _normResid2DGblFitHistNameYvsY;
            static std::string _kinkGblFitHistNameX;
            static std::string _kinkGblFitHistNameY;
        };

#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    };

    /** A global instance of the processor */
    EUTelProcessorTrackingGBLTrackFit gEUTelProcessorTrackingGBLTrackFit;

} // eutelescope

#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H */

