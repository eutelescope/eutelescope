// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
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

// system includes <>
#include <string>
#include <vector>
#include <map>

#include "TMinuit.h"
class TMinuit;



namespace eutelescope {

  //! Alignment processor
  /*! 
   *
   */

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

    // Parameters

    int _alignedPlane;
    double _chi2Cut;

    std::vector<float > _startValuesForAlignment;

    //! Output collection name
    /*! Output collection with fitted tracks.
     */ 
    //    std::string _outputTrackColName;

    //! TRACKERHIT collection name
    /*! Output collection with hits from fitted tracks.
     */ 
    //    std::string _outputHitColName;

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
    
    int _nPlanes;
    double * _xMeasPos;
    double * _yMeasPos;
    double * _zMeasPos;
    double * _intrResol;
    double _xMeas, _yMeas;
    double _xPred, _yPred;

  };

  //! A global instance of the processor
  EUTelAlign gEUTelAlign;      



}
#endif
#endif
