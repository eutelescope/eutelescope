// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelDUTHistograms_h
#define EUTelDUTHistograms_h 1

// eutelescope includes ".h"
//#include "TrackerHitImpl2.h"
#include "IMPL/TrackerHitImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

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


  //! DUT analysis processor for EUDET Telescope
  /*! This processor was designed for analysis of DUT performancs
   *  based on the analytic track fitting results.
   *
   * \par Geometry description
   * Geometry information is taken from GEAR.
   *
   * \par Input
   * \c Track collection with fit results and \c TrackerHit collection
   * with DUT hits are taken as an input.
   *
   * \param InputTrackCollectionName  Name of the input Track collection
   *
   * \param InputHitCollectionName  Name of the input TrackerHit collection,
   *  from which DUT hits are taken
   *
   * \param UseManualDUT Flag for manual DUT selection
   *                      i.e. ignoring GEAR definition
   *
   * \param ManualDUTid  Id of telescope layer which should be used as DUT
   *
   * \param DUTalignment Alignment corrections for DUT: shift in X, Y
   *                     and rotation around Z
   *
   * \param DistMax Maximum allowed distance between fit and matched
   *                 DUT hit.
   *
   * \param DUTpitchX Sensor pitch size in X
   *
   * \param DUTpitchY Sensor pitch size in Y
   *
   * \param HistoInfoFileName Name of the histogram information file.
   *        Using this file histogram parameters can be changed without
   *        recompiling the code.
   *
   * \param DebugEventCount      Print out debug and information
   * messages only for one out of given number of events. If zero, no
   * debug information is printed.
   *
   * \author A.F.Zarnecki, University of Warsaw
   * @version $Id: EUTelDUTHistograms.h,v 1.7 2009-07-15 17:21:28 bulgheroni Exp $
   *
   */


  class EUTelDUTHistograms : public marlin::Processor {

  public:



    //! Returns a new instance of EUTelDUTHistograms
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelDUTHistograms
     */
    virtual Processor*  newProcessor() { return new EUTelDUTHistograms ; }

    //! Default constructor
    EUTelDUTHistograms() ;

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

    virtual int guessSensorID( const double* hit);
    virtual int getClusterSize(int sensorID, TrackerHitImpl * hit, int& sizeX, int& sizeY, int& subMatrix );
    virtual int getSubMatrix(int sensorID, float xlocal);

    //! Called after data processing for clean up.
    /*! Used to release memory allocated in init() step
     */
    virtual void end() ;

  protected:

    //! reference HitCollection name 
    /*!
     */
    std::string      _referenceHitCollectionName;
    bool             _applyToReferenceHitCollection;
    LCCollectionVec* _referenceHitVec;    

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
    std::string _inputTrackColName ;

    //! Input \c TrackerHit collection name
    std::string _inputHitColName ;

    //! Flag for manual DUT selection

    bool _useManualDUT;

    //! Id of telescope layer which should be used as DUT

    int _manualDUTid;

    //!  Debug print out for one out of given number of events.
    int _debugCount ;



    // Internal processor variables
    // ----------------------------


    int _nRun ;
    int _nEvt ;

    int _iDUT;
    double _zDUT;
    double _distMax;

    double _pitchX;
    double _pitchY;

    std::vector<double> _localX;
    std::vector<double> _localY;

    std::vector<int> _clusterSizeX;
    std::vector<int> _clusterSizeY;
    std::vector<int> _subMatrix;

    int  _maptrackid; 
    std::map< int, std::vector<double> >  _trackhitposX;   
    std::map< int, std::vector<double> >  _trackhitposY;
    std::map< int, std::vector<int> >     _trackhitsizeX;
    std::map< int, std::vector<int> >     _trackhitsizeY;
    std::map< int, std::vector<int> >     _trackhitsubM;
    std::map< int, std::vector<int> >     _trackhitsensorID;
 
 
    int _cluSizeXCut ;
    int _cluSizeYCut ;

    int _trackNCluXCut;
    int _trackNCluYCut;



    std::vector<double> _measuredX;
    std::vector<double> _measuredY;

    std::vector<double> _bgmeasuredX;
    std::vector<double> _bgmeasuredY;

    std::map< int, std::vector<double> >  _fittedX;   
    std::map< int, std::vector<double> >  _fittedY;   
 
    std::map< int, std::vector<double> >  _bgfittedX;   
    std::map< int, std::vector<double> >  _bgfittedY;   
 
//obs.igor.280812    std::vector<double> _fittedX;
//obs.igor.280812    std::vector<double> _fittedY;

//obs.igor.280812    std::vector<double> _bgfittedX;
//obs.igor.280812    std::vector<double> _bgfittedY;

    std::vector<float > _DUTalign;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Used to refer to histograms by their names, i.e. to recall
     *  a histogram pointer using histogram name.
     */

    std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;

    static std::string _ClusterSizeXHistoName  ;
    static std::string _ClusterSizeYHistoName  ;
    static std::string _ClusterSizeXYHistoName ;

    static std::string _ClusterSizeXAHistoName  ;
    static std::string _ClusterSizeYAHistoName  ;
    static std::string _ClusterSizeXYAHistoName ;

    static std::string _ClusterSizeXBHistoName  ;
    static std::string _ClusterSizeYBHistoName  ;
    static std::string _ClusterSizeXYBHistoName ;

    static std::string _ClusterSizeXCHistoName  ;
    static std::string _ClusterSizeYCHistoName  ;
    static std::string _ClusterSizeXYCHistoName ;

    static std::string _ClusterSizeXDHistoName  ;
    static std::string _ClusterSizeYDHistoName  ;
    static std::string _ClusterSizeXYDHistoName ;

    static std::string _MeasuredXHistoName;
    static std::string _MeasuredYHistoName;
    static std::string _MeasuredXYHistoName;

    static std::string _MatchedXHistoName;
    static std::string _MatchedYHistoName;
    static std::string _MatchedXYHistoName;

    static std::string _UnMatchedXHistoName;
    static std::string _UnMatchedYHistoName;
    static std::string _UnMatchedXYHistoName;

    static std::string _FittedXHistoName;
    static std::string _FittedYHistoName;
    static std::string _FittedXYHistoName;

    static std::string _EfficiencyXHistoName;
    static std::string _EfficiencyYHistoName;
    static std::string _EfficiencyXYHistoName;

    static std::string _BgEfficiencyXHistoName;
    static std::string _BgEfficiencyYHistoName;
    static std::string _BgEfficiencyXYHistoName;

    static std::string _NoiseXHistoName;
    static std::string _NoiseYHistoName;
    static std::string _NoiseXYHistoName;

    static std::string _ShiftXHistoName;
    static std::string _ShiftYHistoName;
    static std::string _ShiftXYHistoName;

// submatrix A
    static std::string _ShiftXAHistoName;
    static std::string _ShiftYAHistoName;
    static std::string _ShiftXYAHistoName;
    // cluster size 1
    static std::string _ShiftXA1HistoName;
    static std::string _ShiftYA1HistoName;
    static std::string _ShiftXYA1HistoName;
    // cluster size 2
    static std::string _ShiftXA2HistoName;
    static std::string _ShiftYA2HistoName;
    static std::string _ShiftXYA2HistoName;
    // cluster size 3
    static std::string _ShiftXA3HistoName;
    static std::string _ShiftYA3HistoName;
    static std::string _ShiftXYA3HistoName;
    // cluster size 4
    static std::string _ShiftXA4HistoName;
    static std::string _ShiftYA4HistoName;
    static std::string _ShiftXYA4HistoName;
    // cluster size 5
    static std::string _ShiftXA5HistoName;
    static std::string _ShiftYA5HistoName;
    static std::string _ShiftXYA5HistoName;
    // cluster size 6
    static std::string _ShiftXA6HistoName;
    static std::string _ShiftYA6HistoName;
    static std::string _ShiftXYA6HistoName;
    // cluster size 7
    static std::string _ShiftXA7HistoName;
    static std::string _ShiftYA7HistoName;
    static std::string _ShiftXYA7HistoName;


// submatrix B
    static std::string _ShiftXBHistoName;
    static std::string _ShiftYBHistoName;
    static std::string _ShiftXYBHistoName;
    // cluster size 1
    static std::string _ShiftXB1HistoName;
    static std::string _ShiftYB1HistoName;
    static std::string _ShiftXYB1HistoName;
    // cluster size 2
    static std::string _ShiftXB2HistoName;
    static std::string _ShiftYB2HistoName;
    static std::string _ShiftXYB2HistoName;
    // cluster size 3
    static std::string _ShiftXB3HistoName;
    static std::string _ShiftYB3HistoName;
    static std::string _ShiftXYB3HistoName;
    // cluster size 4
    static std::string _ShiftXB4HistoName;
    static std::string _ShiftYB4HistoName;
    static std::string _ShiftXYB4HistoName;
    // cluster size 5
    static std::string _ShiftXB5HistoName;
    static std::string _ShiftYB5HistoName;
    static std::string _ShiftXYB5HistoName;
    // cluster size 6
    static std::string _ShiftXB6HistoName;
    static std::string _ShiftYB6HistoName;
    static std::string _ShiftXYB6HistoName;
    // cluster size 7
    static std::string _ShiftXB7HistoName;
    static std::string _ShiftYB7HistoName;
    static std::string _ShiftXYB7HistoName;


// submatrix C
    static std::string _ShiftXCHistoName;
    static std::string _ShiftYCHistoName;
    static std::string _ShiftXYCHistoName;
    // cluster size 1
    static std::string _ShiftXC1HistoName;
    static std::string _ShiftYC1HistoName;
    static std::string _ShiftXYC1HistoName;
    // cluster size 2
    static std::string _ShiftXC2HistoName;
    static std::string _ShiftYC2HistoName;
    static std::string _ShiftXYC2HistoName;
    // cluster size 3
    static std::string _ShiftXC3HistoName;
    static std::string _ShiftYC3HistoName;
    static std::string _ShiftXYC3HistoName;
    // cluster size 4
    static std::string _ShiftXC4HistoName;
    static std::string _ShiftYC4HistoName;
    static std::string _ShiftXYC4HistoName;
    // cluster size 5
    static std::string _ShiftXC5HistoName;
    static std::string _ShiftYC5HistoName;
    static std::string _ShiftXYC5HistoName;
    // cluster size 6
    static std::string _ShiftXC6HistoName;
    static std::string _ShiftYC6HistoName;
    static std::string _ShiftXYC6HistoName;
    // cluster size 7
    static std::string _ShiftXC7HistoName;
    static std::string _ShiftYC7HistoName;
    static std::string _ShiftXYC7HistoName;


// submatrix D
    static std::string _ShiftXDHistoName;
    static std::string _ShiftYDHistoName;
    static std::string _ShiftXYDHistoName;
    // cluster size 1
    static std::string _ShiftXD1HistoName;
    static std::string _ShiftYD1HistoName;
    static std::string _ShiftXYD1HistoName;
    // cluster size 2
    static std::string _ShiftXD2HistoName;
    static std::string _ShiftYD2HistoName;
    static std::string _ShiftXYD2HistoName;
    // cluster size 3
    static std::string _ShiftXD3HistoName;
    static std::string _ShiftYD3HistoName;
    static std::string _ShiftXYD3HistoName;
    // cluster size 4
    static std::string _ShiftXD4HistoName;
    static std::string _ShiftYD4HistoName;
    static std::string _ShiftXYD4HistoName;
    // cluster size 5
    static std::string _ShiftXD5HistoName;
    static std::string _ShiftYD5HistoName;
    static std::string _ShiftXYD5HistoName;
    // cluster size 6
    static std::string _ShiftXD6HistoName;
    static std::string _ShiftYD6HistoName;
    static std::string _ShiftXYD6HistoName;
    // cluster size 7
    static std::string _ShiftXD7HistoName;
    static std::string _ShiftYD7HistoName;
    static std::string _ShiftXYD7HistoName;




    static std::string _BgShiftXHistoName;
    static std::string _BgShiftYHistoName;
    static std::string _BgShiftXYHistoName;

    static std::string _ShiftXvsYHistoName;
    static std::string _ShiftYvsXHistoName;
    static std::string _ShiftXvsY2DHistoName;
    static std::string _ShiftYvsX2DHistoName;

    static std::string _ShiftXvsXHistoName;
    static std::string _ShiftYvsYHistoName;
    static std::string _ShiftXvsX2DHistoName;
    static std::string _ShiftYvsY2DHistoName;

    static std::string _EtaXHistoName;
    static std::string _EtaYHistoName;
    static std::string _EtaX2DHistoName;
    static std::string _EtaY2DHistoName;
    static std::string _EtaX3DHistoName;
    static std::string _EtaY3DHistoName;

    static std::string _PixelEfficiencyHistoName    ;
    static std::string _PixelResolutionXHistoName    ;
    static std::string _PixelResolutionYHistoName    ;
    static std::string _PixelChargeSharingHistoName ;

#endif

  } ;


  //! A global instance of the processor.
  EUTelDUTHistograms aEUTelDUTHistograms ;


}

#endif



