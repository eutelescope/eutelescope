/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_GEAR

#ifndef EUTELXCORRELATOR_H
#define EUTELXCORRELATOR_H

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram2D.h>
#endif


// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {

  //! Hit cross correlator
  /*! This processor can be used to study the correlation among the
   *  telescope hits and the one detected by the DUT.
   *
   *  The user can select if to make the correlation among all
   *  possible planes in the telescope and the DUT or just among a
   *  specified plane and the DUT.
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Telescope hit collection</b>: A collection containing the
   *  telescope hits
   *
   *  <b>DUT hit collection</b>: A collection containing the DUT hits
   *
   *  <h3>Warning</h3>
   *  When doing the cross correlation after the fitting and then
   *  using as a hit collection one of the two produced by the fitter
   *  itself, <b>NEVER</b> use the @a fithit collection because this
   *  contains also the extra- / inter- polated position at the DUT
   *  plane. Instead use the @a corrfit or any other hit collection
   *  available before the fitter.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   *
   */
  class EUTelExampleProcessorCorrelator : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelExampleProcessorCorrelator
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelExampleProcessorCorrelator.
     */
    virtual Processor * newProcessor() {
      return new EUTelExampleProcessorCorrelator;
    }

    //! Default constructor
    EUTelExampleProcessorCorrelator ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! For each event, we loop over both the telescope and DUT input
     *  hit collections and for each clusters we will be filling in a
     *  2D correlation histogram.
     *
     *  Histograms are booked in the first event.
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
    /*! This method is used to book all the correlation histograms. It
     *  is called by processEvent when processing the first event.
     */
    void bookHistos();


  protected:

    //! Telescope input hit collection name
    std::string _inputTelescopeCollectionName;

    //! DUT input hit collection name
    std::string _inputDUTCollectionName;

    //! A function to guess the sensorID of a hit
    /*! It is checking against the distance of each plane assuming
     *  that this hit is belonging to the plane at the closest distant.
     */
    int guessSensorID( TrackerHitImpl * hit ) ;

  private:

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

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
    /*! This is the real geometry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    //! An array with the Z position of planes
    double * _siPlaneZPosition;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the processEvent
     *  during the first event and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    //! Correlation histogram matrix
    /*! This is used to store the pointers of each histogram
     */
    std::map< unsigned int , AIDA::IHistogram2D* >  _hitXCorrelationMatrix;
    std::map< unsigned int , AIDA::IHistogram2D* >  _hitYCorrelationMatrix;

    //! Base name of the correlation histogram
    static std::string _hitXCorrelationHistoName;
    static std::string _hitYCorrelationHistoName;
#endif

  };

  //! A global instance of the processor
  EUTelExampleProcessorCorrelator gEUTelExampleProcessorCorrelator;


}

#endif
#endif
