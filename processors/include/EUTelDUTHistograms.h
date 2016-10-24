
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
#include "EUTELESCOPE.h"

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
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>

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
   * \author A.F.Zarnecki, University of Warsaw
   * @version $Id$
   *
   */


  class EUTelDUTHistograms : public marlin::Processor {

  public:

    virtual int  read_track(LCEvent *event); 
    virtual int  read_track_from_collections(LCEvent *event); 


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

    virtual int getClusterSize(int sensorID, TrackerHit * hit, int& sizeX, int& sizeY, int& subMatrix );
    virtual int getSubMatrix(int sensorID, float xlocal);

    //! Called after data processing for clean up.
    /*! Used to release memory allocated in init() step
     */
    virtual void end() ;

private:
  DISALLOW_COPY_AND_ASSIGN(EUTelDUTHistograms)
  
protected:
    //! reference HitCollection name 
    /*!
     */
    std::string      _referenceHitCollectionName;
    bool             _useReferenceHitCollection;
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
    std::string _inputRecHitColName ;
    std::string _inputFitHitColName ;

    //! Id of telescope layer which should be used as DUT
    int _iDUT;

    // Internal processor variables
    // ----------------------------


    int _nRun ;

    double _zDUT;
    double _distMax;

    double _pitchX;
    double _pitchY;


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

    std::map< int, std::vector<double> >  _localX;   
    std::map< int, std::vector<double> >  _localY;   

    std::map< int, std::vector<double> >  _fittedX;   
    std::map< int, std::vector<double> >  _fittedY;   
 
    std::map< int, std::vector<double> >  _bgfittedX;   
    std::map< int, std::vector<double> >  _bgfittedY;   
 
    std::vector<float > _DUTalign;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram maps
    /*! Used to manage to histograms pointers in a flexible and
      ! efficient way.  
     */

    // ! Histograms for full detector and for sub matrices are handled by enums to make the distiction clear in the code.
    enum detMatrix {SubMatrixA, SubMatrixB, SubMatrixC, SubMatrixD, FullDetector};

    // ! Some histograms are filled for a particular cluster size; this determines the max. cluster size considered
    static const int HistoMaxClusterSize = 7;

    enum projAxis {projX, projY, projXY};

    std::map< projAxis, std::map< detMatrix, AIDA::IBaseHistogram*> > _ClusterSizeHistos;
    std::map< projAxis, std::map< detMatrix, std::map< int, AIDA::IBaseHistogram*> > > _ShiftHistos; // w/ cluster size studies

    std::map< projAxis, AIDA::IBaseHistogram*> _MeasuredHistos;
    std::map< projAxis, AIDA::IBaseHistogram*> _MatchedHistos;
    std::map< projAxis, AIDA::IBaseHistogram*> _UnMatchedHistos;
    std::map< projAxis, AIDA::IBaseHistogram*> _FittedHistos;
    std::map< projAxis, AIDA::IBaseHistogram*> _EfficiencyHistos;
    std::map< projAxis, AIDA::IBaseHistogram*> _BgEfficiencyHistos;
    std::map< projAxis, AIDA::IBaseHistogram*> _NoiseHistos;
    std::map< projAxis, AIDA::IBaseHistogram*> _BgShiftHistos;

    AIDA::IProfile1D* _ShiftXvsYHisto;
    AIDA::IProfile1D* _ShiftYvsXHisto;
    AIDA::IProfile1D* _ShiftXvsXHisto;
    AIDA::IProfile1D* _ShiftYvsYHisto;

    AIDA::IHistogram2D* _ShiftXvsY2DHisto;
    AIDA::IHistogram2D* _ShiftYvsX2DHisto;
    AIDA::IHistogram2D* _ShiftXvsX2DHisto;
    AIDA::IHistogram2D* _ShiftYvsY2DHisto;

    AIDA::IProfile1D* _EtaXHisto;
    AIDA::IProfile1D* _EtaYHisto;
    AIDA::IHistogram2D* _EtaX2DHisto;
    AIDA::IHistogram2D* _EtaY2DHisto;
    AIDA::IProfile2D* _EtaX3DHisto;
    AIDA::IProfile2D* _EtaY3DHisto;

    AIDA::IProfile2D* _PixelEfficiencyHisto    ;
    AIDA::IProfile2D* _PixelResolutionXHisto   ;
    AIDA::IProfile2D* _PixelResolutionYHisto   ;
    AIDA::IProfile2D* _PixelChargeSharingHisto ;

#endif

  } ;


  //! A global instance of the processor.
  EUTelDUTHistograms aEUTelDUTHistograms ;


}

#endif



