// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELAPIXKALMAN_H
#define EUTELAPIXKALMAN_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "mille/Mille.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

//APIX includes
#include "APIXFitter.h"

namespace eutelescope {
  class EUTelAPIXKalman : public marlin::Processor {
  public:
    //! Returns a new instance of EUTelAPIXKalman
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelAPIXKalman.
     */
    virtual Processor * newProcessor() {
      return new EUTelAPIXKalman;
    }

    //! Default constructor
    EUTelAPIXKalman ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. Each element of the
     *  pulse collection is scanned and the center of the cluster is
     *  translated into the external frame of reference thanks to the
     *  GEAR geometry description.
     *
     *  The cluster center might be calculate using a standard linear
     *  charge center of gravity algortihm or applying a more
     *  sophisticated non linear eta function. This behaviour is
     *  regulated by the user from the steering file.
     *
     *  @throw UnknownDataTypeException if the cluster type is unknown
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


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


  protected:
    //Counters
    int cutAngleX, cutAngleY, cutChi2, cutResids, cutInTime;

    //! TrackerHit collection name
    /*! Input collection with hits.*/
    std::vector<std::string > _hitCollectionName;


    //! TRACK collection name
    /*! Output collection with fitted tracks.
     */
    std::string _trackCollectionName;

    // Define output track collection
    LCCollectionVec     * _fittrackvec;

    // parameters
    std::vector<int > _excludePlanes;
    std::vector<int > _fixedPlanes;
    std::vector<int > _fixedTranslations;
    std::vector<int > _fixedZRotations;
    std::vector<int > _fixedScales;
    std::vector<int> _inTimeCheck;
    bool _fixedX, _fixedY;
    std::string _binaryFilename;

    std::vector<float> _telescopeResolution;
    bool _useResidualCuts;

    std::vector<float > _residualsXMin;
    std::vector<float > _residualsYMin;
    std::vector<float > _residualsXMax;
    std::vector<float > _residualsYMax;

    std::vector<float> _shiftsX;
    std::vector<float> _shiftsY;
    std::vector<float> _scalesX;
    std::vector<float> _scalesY;

    int _generatePedeSteerfile;
    std::string _pedeSteerfileName;
    bool _useHitResol;
    bool _runPede;
    bool _addToLCIO;
    bool _usePedeUserStartValues;
    bool _doScatter;
    std::vector<float > _pedeUserStartValuesX;
    std::vector<float > _pedeUserStartValuesY;
    std::vector<float > _pedeUserStartValuesGamma;

    std::string _alignmentConstantLCIOFile;
    std::string _alignmentConstantCollectionName;

    float _normalizedResidualsMax;
    double _eBeam;
    int _nSkipMax;

    double _minDxDz;
    double _maxDxDz;
    double _minDyDz;
    double _maxDyDz;
    double _maxChi2;
    
  private:
    int getPlaneIndex(double zPos);
    void fitPermutations(int layerIndex, APIXFitter::FitPlane* pl, APIXFitter::TrackEstimate* estimate, int nSkipped);
    void finalizeTrack();
    void addToMille();
    void plotResiduals(double chi2, int ndof);
    void tryFill(string name, double value);
    void tryBook(string name, string title, int bins, double min, double max);
    void readHitCollection(LCEvent* event);
    std::vector<double> initAlignParams();
    bool generatePedeSteeringFile(std::vector<double> &params, bool shifts, bool rotate, bool scale);
    void runPede(std::vector<double> &params);
    void addToLCIO(double chi2, int ndof);
    void addPlaneHit(APIXFitter::FitPlane* pl, TrackerHitImpl*);
    double getScatterCov(int index);
    bool inTimeGood(APIXFitter::FitPlane* pl);
    bool goodResiduals(APIXFitter::FitPlane* pl);
    APIXFitter::APIXKalman* _fitter;
    std::map<double, int> _zSort;
    
    std::vector< std::vector<TrackerHitImpl*> > _planeHits;

    bool _checkInTime;

    int _nTracks;
    int _expectedTracks;
    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

    // Statistics
    int _nMilleDataPoints;
    int _nMilleTracks;

    // Mille
    Mille * _mille;

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
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    static std::string _numberTracksLocalname;
    static std::string _skippedLocalname;

    static std::string _chi2Localname;
    static std::string _ndofLocalname;

    static std::string _dxdzLocalname;
    static std::string _dydzLocalname;

    static std::string _residualXLocalname;
    static std::string _residualYLocalname;


#endif
    std::vector<double> _siPlaneZPosition;

    //! Fill histogram switch
    /*! Only for debug reason
     */
    bool _histogramSwitch;
  };
  //! A global instance of the processor
  EUTelAPIXKalman gEUTelAPIXKalman;
}
#endif
#endif
