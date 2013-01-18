// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif


// system includes <>
#include <string>
#include <vector>
#include <map>

// define a specific verbosity level for Multilinefitter
namespace streamlog {

  DEFINE_STREAMLOG_LEVEL( MULTIFITTERMESSAGE,  "MULTIFITTERMESSAGE", message_base_level - 50 + 2 , STREAMLOG_MESSAGE_ACTIVE ) 

}


namespace eutelescope {

  //! Straight line fit processor
  /*!
   *
   */

  class EUTelMultiLineFit : public marlin::Processor {

  public:

    //! Variables for hit parameters
    class HitsInPlane {
    public:
      double measuredX;
      double measuredY;
      double measuredZ;
      double seedCharge;
      double clusterCharge;
    };

    virtual void FitTrack(int nPlanesFitter, double xPosFitter[], double yPosFitter[], double zPosFitter[], double xResFit[], double yResFit[], double chi2Fit[2], double residXFit[], double residYFit[], double angleFit[2]);

    //! Returns a new instance of EUTelMultiLineFit
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelMultiLineFit.
     */
    virtual Processor * newProcessor() {
      return new EUTelMultiLineFit;
    }

    //! Default constructor
    EUTelMultiLineFit ();

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

    //! TrackerHit collection name
    /*! Input collection with hits.
     */
    std::string _hitCollectionName;

    //! TRACK collection name
    /*! Output collection with fitted tracks.
     */
    std::string _outputTrackColName;

    //! TRACKERHIT collection name
    /*! Output collection with hits from fitted tracks.
     */
    std::string _correctedHitColName;
    std::string _outputHitColName;

    int _alignmentMode;

    std::vector<float > _alignmentConstantsSecondLayer;
    std::vector<float > _alignmentConstantsThirdLayer;
    std::vector<float > _alignmentConstantsFourthLayer;
    std::vector<float > _alignmentConstantsFifthLayer;
    std::vector<float > _alignmentConstantsSixthLayer;

    std::vector<float > _millepedeConstantsFirstLayer;
    std::vector<float > _millepedeConstantsSecondLayer;
    std::vector<float > _millepedeConstantsThirdLayer;
    std::vector<float > _millepedeConstantsFourthLayer;
    std::vector<float > _millepedeConstantsFifthLayer;
    std::vector<float > _millepedeConstantsSixthLayer;

    float _distanceMax;
    float _distanceDUTMax;
    float _chi2XMax;
    float _chi2YMax;
    int _excludePlane;
    int _maxTrackCandidates;
    int _maxHitsPlane;
    float _hitDistanceXMax;
    float _hitDistanceYMax;

  private:

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

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
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    static std::string _numberTracksLocalname;

    static std::string _chi2XLocalname;
    static std::string _chi2YLocalname;

    static std::string _angleXLocalname;
    static std::string _angleYLocalname;

    static std::string _residualXLocalname;
    static std::string _residualYLocalname;

    static std::string _positionXYLocalname;

    static std::string _seedChargeLocalname;
    static std::string _clusterChargeLocalname;

    static std::string _hitDistanceLocalname;
#endif

    int _nPlanes;

    double ** _xPos;
    double ** _yPos;
    double ** _zPos;
    double ** _seedCharge;
    double ** _clusterCharge;

    double * _xPosHere;
    double * _yPosHere;
    double * _zPosHere;
    double * _waferResidX;
    double * _waferResidY;
    double * _intrResolX;
    double * _intrResolY;
    double * _xFitPos;
    double * _yFitPos;

    //! Fill histogram switch
    /*! Only for debug reason
     */
    bool _histogramSwitch;

    // arrays for efficiency calculation
    int _iHitDUT;
    int _iHitDUTDistanceCut[10];

  };

  //! A global instance of the processor
  EUTelMultiLineFit gEUTelMultiLineFit;

}
#endif
#endif
