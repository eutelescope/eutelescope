
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef CMSPixelClusterTuple_h
#define CMSPixelClusterTuple_h 1

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/ITuple.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {


  //! Processor for building n-tuple with clustering results
  /*! This processor prepared for easier analysis of clustering results.
   *
   * \par Geometry description
   * Geometry information is taken from GEAR. However, it is possible
   * to select DUT layer ID manually.
   *
   * \par Input
   * \c Track collection is taken as an input.
   * For each track referenced \c TrackerHit entries are analysed.
   * Measured and fitted hits can be distinguished by looking into
   * hit type (type <=31 for measured hits, type >=32 for fitted).
   * DUT measuterents are taken directly from \c TrackerHit
   * collection. Hits belonging to DUT are identified by Z position.
   *
   * \param InputCollectionName  Name of the input  Tracker collection
   *
   * \param DebugEventCount      Print out debug and information
   *        messages only for one out of given number of events. If
   *        zero, no debug information is printed.
   *
   * \author U. Grundler, National Tawain University
   * Based on EUTelFitTuple
   *
   */


  class CMSPixelClusterTuple : public marlin::Processor {

  public:



    //! Returns a new instance of CMSPixelClusterTuple
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new CMSPixelClusterTuple
     */
    virtual Processor*  newProcessor() { return new CMSPixelClusterTuple ; }

    //! Default constructor
    CMSPixelClusterTuple() ;

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


    //! Input \c cluster collection name
    std::string _clusterColName ;

    //!  Debug print out for one out of given number of events.
    int _debugCount ;

    //!  Value to be used for missing measurements
    double _missingValue;

    // Setup description

    int _nTelPlanes;

    int * _planeSort;
    int * _planeID;
    double * _planePosition;
    bool   * _isActive;

    bool   * _isMeasured;
    double * _measuredX;
    double * _measuredY;
    double * _measuredZ;
    double * _measuredQ;

    bool   * _isFitted;
    double * _fittedX;
    double * _fittedY;

    // Internal processor variables
    // ----------------------------

    int _nRun ;
    int _nEvt ;
    int _runNr;
    int _evtNr;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)


    static std::string _ClusterTupleName;

    AIDA::ITuple * _ClusterTuple;

#endif

  } ;


  //! A global instance of the processor.
  CMSPixelClusterTuple aCMSPixelClusterTuple ;


}

#endif



