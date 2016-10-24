
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelFitHistograms_h
#define EUTelFitHistograms_h 1

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {


  //! Fitted track analysis processor for EUDET Telescope
  /*! This processor was designed for checking track fitting results.
   *  Based on measured and fitted hit positions in telescope planes
   *  alignment of the telescope is also verified.
   *
   * \par Geometry description
   * Geometry information is taken from GEAR.
   *
   * \par Input
   * \c Track collection is taken as an input.
   * For each track referenced \c TrackerHit entries are analysed.
   * Measured and fitted hits can be distinguished by looking into
   * hit type (type <=31 for measured hits, type >=32 for fitted).
   *
   * \param InputCollectionName  Name of the input  Tracker collection
   *
   * \param BeamReferenceID ID of the layer used to check alignment
   *        relative to the beam direction
   *
   * \param TelescopeReferenceIDs IDs of two layers used to check
   *        internal telescope alignment
   *
   * \param HistoInfoFileName Name of the histogram information file.
   *        Using this file histogram parameters can be changed without
   *        recompiling the code.
   * \param DebugEventCount      Print out debug and information
   * messages only for one out of given number of events. If zero, no
   * debug information is printed.
   *
   * \author A.F.Zarnecki, University of Warsaw
   * @version $Id$
   * \date 2007.09.10
   *
   */


  class EUTelFitHistograms : public marlin::Processor {

  public:



    //! Returns a new instance of EUTelFitHistograms
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelFitHistograms
     */
    virtual Processor*  newProcessor() { return new EUTelFitHistograms ; }

    //! Default constructor
    EUTelFitHistograms() ;

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution.
     *
     */
    virtual void init() ;

    //! Called for every run.
    /*!
     * @param run the LCRunHeader of the current run
     */
    virtual void processRunHeader( LCRunHeader* run ) ;

    //! Called every event
    /*! This is called for each event in the file.
     *
     *  @param evt the current LCEvent event
     */
    virtual void processEvent( LCEvent * evt ) ;

    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. It can be used to fill check
     *  plots. For the time being there is nothing to check and do in
     *  this slot.
     *
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check( LCEvent * evt ) ;


    //! Book histograms
    /*! This method is used to books all required
     *  histograms. Histogram pointers are stored into
     *  _aidaHistoMap so that they can be recalled and filled
     * from anywhere in the code.
     */
    void bookHistos();


    //! Called after data processing for clean up.
    /*! Used to release memory allocated in init() step
     */
    virtual void end() ;

  protected:

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

    //! The histogram information file
    /*! This string contain the name of the histogram information
     *  file. This is selected by the user in the steering file.
     *
     *  @see eutelescope::EUTelHistogramManager
     *  @see eutelescope::EUTelHistogramInfo
     */
    std::string _histoInfoFileName;

    //! Input \c Track collection name
    std::string _inputColName ;

    //! Flag for alignment control histograms
    bool _alignCheckHistograms;

    //! ID of the layer used for beam based alignment check
    int _BeamReferenceID;

    //! Local index of this layer
    int _beamID;

    //! IDs of two layers used to check internal telescope alignment
    std::vector<int > _TelescopeReferenceIDs;

    //! Local indexes of the two reference layers
    int _referenceID0;
    int _referenceID1;

    //!  Debug print out for one out of given number of events.
    int _debugCount ;

    // Setup description

    int _nTelPlanes;
    int _iDUT;

    int * _planeSort;
    int * _planeID;
    double * _planePosition;
    bool   * _isActive;

    bool   * _isMeasured;
    double * _measuredX;
    double * _measuredY;
    double * _measuredQ;

    bool   * _isFitted;
    double * _fittedX;
    double * _fittedY;

    // Internal processor variables
    // ----------------------------

    int _nRun ;
    int _nEvt ;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Used to refer to histograms by their names, i.e. to recall
     *  a histogram pointer using histogram name.
     */

    std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;
    
    static std::string _ShiftXvsYHistoName;
    static std::string _ShiftYvsXHistoName;

    static std::string _MeasuredXHistoName;
    static std::string _MeasuredYHistoName;
    static std::string _MeasuredXYHistoName;

    static std::string _FittedXHistoName;
    static std::string _FittedYHistoName;
    static std::string _FittedXYHistoName;

    static std::string _ResidualXHistoName;
    static std::string _ResidualYHistoName;
    static std::string _ResidualXYHistoName;

    static std::string _ScatXHistoName;
    static std::string _ScatYHistoName;
    static std::string _ScatXYHistoName;

    static std::string _AngleXHistoName;
    static std::string _AngleYHistoName;
    static std::string _AngleXYHistoName;

    static std::string _beamShiftXHistoName;
    static std::string _beamShiftYHistoName;
    static std::string _beamShiftXYHistoName;

    static std::string _clusterSignalHistoName;
    static std::string _meanSignalXHistoName;
    static std::string _meanSignalYHistoName;
    static std::string _meanSignalXYHistoName;

    static std::string _beamRotXHistoName;
    static std::string _beamRotYHistoName;
    static std::string _beamRotX2DHistoName;
    static std::string _beamRotY2DHistoName;
    static std::string _beamRot2XHistoName;
    static std::string _beamRot2YHistoName;

    static std::string _relShiftXHistoName;
    static std::string _relShiftYHistoName;
    static std::string _relRotXHistoName;
    static std::string _relRotYHistoName;
    static std::string _relRotX2DHistoName;
    static std::string _relRotY2DHistoName;

#endif

  } ;


  //! A global instance of the processor.
  EUTelFitHistograms aEUTelFitHistograms ;


}

#endif



