// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELALIGN_H
#define EUTELALIGN_H

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

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMinuit.h"
class TMinuit;
#endif


namespace eutelescope {

  class EUTelAlign : public marlin::Processor {

  public:

    static void Chi2Function(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

    //! Variables for hit parameters
    class HitsForFit {
    public:
      double firstLayerMeasuredX;
      double firstLayerMeasuredY;
      double firstLayerMeasuredZ;
      double secondLayerPredictedX;
      double secondLayerPredictedY;
      double secondLayerPredictedZ;
      double secondLayerMeasuredX;
      double secondLayerMeasuredY;
      double secondLayerMeasuredZ;
      double firstLayerResolution;
      double secondLayerResolution;
    };

    class HitsInFirstBox {
    public:
      double measuredX;
      double measuredY;
      double measuredZ;
    };

    virtual void FitTrack(int nPlanesFit, double xPosFit[], double yPosFit[], double zPosFit[], double xResFit[], double yResFit[], double chi2Fit[2], double&, double&, double);

    //! Returns a new instance of EUTelAlign
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelAlign.
     */
    virtual Processor * newProcessor() {
      return new EUTelAlign;
    }

    //! Default constructor
    EUTelAlign ();

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

    static std::vector<HitsForFit> _hitsForFit;

    //! TrackerHit collection name
    /*! Input collection with measured hits.
     */
    std::string _measHitCollectionName;

    std::vector<float > _alignmentConstantsSecondLayer;
    std::vector<float > _alignmentConstantsThirdLayer;

    // Parameters

    int _alignedPlane;
    int _alignedBox;
    int _nPlanesFirstBox;
    int _referencePlane;
    double _resolution;
    double _chi2Cut;
    double _distanceMin;
    double _distanceMax;
    int _nHitsMax;

    std::vector<float > _startValuesForAlignment;

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

    static std::string _distanceLocalname;

    static std::string _residualXSimpleLocalname;
    static std::string _residualYSimpleLocalname;

    static std::string _residualXLocalname;
    static std::string _residualYLocalname;
#endif

    int _nPlanes;
    double * _xMeasPos;
    double * _yMeasPos;
    double * _zMeasPos;
    double * _intrResol;
    double _xMeas, _yMeas;
    double _xPred, _yPred;

    double ** _xPos;
    double ** _yPos;
    double ** _zPos;

    double * _xPosHere;
    double * _yPosHere;
    double * _zPosHere;
    double * _waferResidX;
    double * _waferResidY;
    double * _intrResolX;
    double * _intrResolY;
    double * _xFitPos;
    double * _yFitPos;

  };

  //! A global instance of the processor
  EUTelAlign gEUTelAlign;
}
#endif
#endif
