// Version: $Id: EUTelMilleGBL.h 2285 2013-01-18 13:46:44Z hperrey $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_GBL

#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelTrackFinder.h"
#include "EUTelTrackFitter.h"
#include "EUTelUtility.h"

// lcio includes <.h>
#include "IMPL/TrackerHitImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "mille/Mille.h"

// GBL
#include "include/MilleBinary.h"

// gear includes <.h>
#include <gear/SiPlanesLayerLayout.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRunHeader.h>
#include "LCIOSTLTypes.h"


// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#endif

// system includes <>
#include <map>
#include <string>
#include <vector>

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TMath.h>
#include <TMinuit.h>
#include <TSystem.h>
#endif


namespace eutelescope {

    class EUTelMilleGBL : public marlin::Processor {
    public:

        EUTelTrackFinder* _theFinder;
        EUTelTrackFitter* _theFitter;

        //! Returns a new instance of EUTelMille

        /*! This method returns a new instance of this processor.  It is
         *  called by Marlin execution framework and it shouldn't be
         *  called/used by the final user.
         *
         *  @return a new EUTelMille.
         */
        virtual Processor * newProcessor() {
            return new EUTelMilleGBL;
        }

        //! Default constructor
        EUTelMilleGBL();

        //! Called at the job beginning.
        /*! This is executed only once in the whole execution. It prints
         *  out the processor parameters and check that the GEAR
         *  environment is properly set up and accessible from Marlin.
         */
        virtual void init();

        //! Called for every run.
        /*! It is called for every run, and consequently the run counter
         *  is incremented. The geometry ID of the file is compared with
         *  the one provided by the GEAR geometry description. In case the
         *  two are different, the user is asked to decide to quit or to
         *  continue with a description that might be wrong.
         *
         *  @param run the LCRunHeader of the this current run
         */
        virtual void processRunHeader(LCRunHeader * run);

        //! Called every event
        /*! This is called for each event in the file. Each element of the
         *  pulse collection is scanned and the center of the cluster is
         *  translated into the external frame of reference thanks to the
         *  GEAR geometry description.
         *
         *  The cluster center might be calculate using a standard linear
         *  charge center of gravity algorithm or applying a more
         *  sophisticated non linear eta function. This behaviour is
         *  regulated by the user from the steering file.
         *
         *  @throw UnknownDataTypeException if the cluster type is unknown
         *
         *  @param evt the current LCEvent event as passed by the
         *  ProcessMgr
         */
        virtual void processEvent(LCEvent * evt);


        //! Called after data processing.
        /*! This method is called when the loop on events is
         *  finished.
         */
        virtual void end();


        //! Histogram booking
        /*! Some control histograms are filled during this procedure in
         *  order to be able to perform easy check on the quality of the
         *  output hits and also to understand if the frame of reference
         *  conversion has been properly done. Of course this method is
         *  effectively doing something only in the case MARLIN_USE_AIDA.
         */
        void bookHistos();

        virtual int guessSensorID(double* hit);

        virtual inline int getAllowedMissingHits() {
            return _allowedMissingHits;
        }

        virtual inline int getMaxTrackCandidates() {
            return _maxTrackCandidates;
        }

        virtual inline int getMimosa26ClusterChargeMin() {
            return _mimosa26ClusterChargeMin;
        }


    protected:


        //! Ordered sensor ID
        /*! Within the processor all the loops are done up to _nPlanes and
         *  according to their position along the Z axis (beam axis).
         *
         *  This vector is containing the sensorID sorted according to the
         *  same rule.
         */
        IntVec _orderedSensorID;
        IntVec _orderedSensorID_wo_excluded;


        //! reference HitCollection name 
        /*!
         */
        std::string _referenceHitCollectionName;
        bool _applyToReferenceHitCollection;
        LCCollectionVec* _referenceHitVec;

        //! TrackerHit collection name
        /*! Input collection with hits.
         */
        StringVec _hitCollectionName;

        //! TRACK collection name
        /*! Output collection with fitted tracks.
         */
        std::string _trackCollectionName;

        //! Hot pixel collection name.
        /*! 
         * this collection is saved in a db file to be used at the clustering level
         */
        std::string _hotPixelCollectionName;

        //! Vector of map arrays, keeps record of hit pixels 
        /*! The vector elements are sorted by Detector ID
         *  For each Detector unique ID element a map of pixels is created. 
         *  first level key   sensor unique 
         *              value sensor map
         *  sensor map key    unique row number
         *             value  vector of column numbers.
         */

        std::map<std::string, bool > _hotPixelMap;

        //! Sensor ID vector
        IntVec _sensorIDVec;

        //! Sensor ID map (inverse sensorIDVec) 
        std::map< int, int > _sensorIDVecMap;
        //! Sensor ID vector, 
        /*! it's position along Z axis
         */
        IntVec _sensorIDVecZOrder;
        //! sensor ID to position along Z id
        /*!
         */
        std::map<int, int> _sensorIDtoZOrderMap;


        // parameters

        float _distanceMax;
        FloatVec _distanceMaxVec;
        std::vector<unsigned int > _excludePlanes; //only for internal usage
        IntVec _excludePlanes_sensorIDs; //this is going to be
        //set by the user.
        IntVec _FixedPlanes; //only for internal usage
        IntVec _FixedPlanes_sensorIDs; //this is going to be
        //set by the user.



        int _maxTrackCandidates;
        int _maxTrackCandidatesTotal;

        std::string _binaryFilename;     

        float _telescopeResolution;
        int _onlySingleHitEvents;
        int _onlySingleTrackEvents;
        int _alignMode;
        int _useResidualCuts;

        FloatVec _residualsXMin;
        FloatVec _residualsYMin;
        FloatVec _residualsXMax;
        FloatVec _residualsYMax;

        FloatVec _resolutionX;
        FloatVec _resolutionY;
        FloatVec _resolutionZ;

        IntVec _FixParameter;


        int _generatePedeSteerfile;
        std::string _pedeSteerfileName;
        std::string _pedeConstraintsFilename;
        int _runPede;
        int _usePedeUserStartValues;
        FloatVec _pedeUserStartValuesX;
        FloatVec _pedeUserStartValuesY;
        FloatVec _pedeUserStartValuesZ;

        FloatVec _pedeUserStartValuesAlpha;
        FloatVec _pedeUserStartValuesBeta;
        FloatVec _pedeUserStartValuesGamma;

        int _inputMode;
        int _allowedMissingHits;
        int _mimosa26ClusterChargeMin;

        float _testModeSensorResolution;
        float _testModeXTrackSlope;
        float _testModeYTrackSlope;

        FloatVec _testModeSensorZPositions;

        FloatVec _testModeSensorXShifts;
        FloatVec _testModeSensorYShifts;
        FloatVec _testModeSensorGamma;
        FloatVec _testModeSensorAlpha;
        FloatVec _testModeSensorBeta;

        std::string _alignmentConstantLCIOFile;
        std::string _alignmentConstantCollectionName;

        std::vector<int> _useSensorRectangular;

    private:

        //! Run number
        int _iRun;

        //! Event number
        int _iEvt;

        // Excluded planes
        int _nExcludePlanes;

        // Statistics
        int _nMilleDataPoints;
        int _nMilleTracks;

        // Mille
        gbl::MilleBinary * _milleGBL;

        //! Conversion ID map.
        /*! In the data file, each cluster is tagged with a detector ID
         *  identify the sensor it belongs to. In the geometry
         *  description, there are along with the sensors also "passive"
         *  layers and other stuff. Those are identify by a layerindex. So
         *  we need a conversion table to go from the detectorID to the
         *  layerindex.
         */
        std::map< int, int > _conversionIdMap;

        //! Silicon planes parameters as described in GEAR
        /*! This structure actually contains the following:
         *  @li A reference to the telescope geoemtry and layout
         *  @li An integer number saying if the telescope is w/ or w/o DUT
         *  @li An integer number saying the number of planes in the
         *  telescope.
         *
         *  This object is provided by GEAR during the init() phase and
         *  stored here for local use.
         */
        gear::SiPlanesParameters * _siPlanesParameters;

        //! Silicon plane layer layout
        /*! This is the real geoemetry description. For each layer
         *  composing the telescope the relevant information are
         *  available.
         *
         *  This object is taken from the _siPlanesParameters during the
         *  init() phase and stored for local use
         */
        gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        //! AIDA histogram map
        /*! Instead of putting several pointers to AIDA histograms as
         *  class members, histograms are booked in the init() method and
         *  their pointers are inserted into this map keyed by their
         *  names.
         *  The histogram filling can proceed recalling an object through
         *  its name
         */
        std::map<std::string, AIDA::IHistogram1D * > _aidaHistoMap1D;
        std::map<std::string, AIDA::IHistogram2D * > _aidaHistoMap2D;

        //! Names of histograms
        static std::string _numberTracksCandidatesHistName;
        static std::string _chi2GblFitHistName;
        static std::string _probGblFitHistName;
        static std::string _residGblFitHistName;
        static std::string _residGblFitHistNameX;
        static std::string _residGblFitHistNameY;
        static std::string _resid2DGblFitHistNameXvsX;
        static std::string _resid2DGblFitHistNameXvsY;
        static std::string _resid2DGblFitHistNameYvsX;
        static std::string _resid2DGblFitHistNameYvsY;
        static std::string _kinkGblFitHistNameX;
        static std::string _kinkGblFitHistNameY;

#endif

        size_t _nPlanes;

        std::vector<std::vector<double> > _xPos;
        std::vector<std::vector<double> > _yPos;
        std::vector<std::vector<double> > _zPos;

        double * _xPosHere;
        double * _yPosHere;
        double * _zPosHere;
        double * _waferResidX;
        double * _waferResidY;
        double * _waferResidZ;
        double * _telescopeResolX;
        double * _telescopeResolY;
        double * _telescopeResolZ;
        double * _xFitPos;
        double * _yFitPos;

        DoubleVec _siPlaneZPosition;

        //! Fill histogram switch
        /*! Only for debug reason
         */
        bool _histogramSwitch;

        //! Limits the pixels on each sensor-plane to a sub-rectangular  
        Utility::RectangularArray _rect;

    };

    //! A global instance of the processor
    EUTelMilleGBL gEUTelMilleGBL;

}
#endif
#endif

#endif
