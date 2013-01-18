
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelAPIXHistograms_h
#define EUTelAPIXHistograms_h 1

#include "EUTelAlignmentConstant.h"

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

#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <IMPL/TrackImpl.h>

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
   * @version $Id$
   *
   */


  class EUTelAPIXHistograms : public marlin::Processor {

  public:



    //! Returns a new instance of EUTelAPIXHistograms
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelAPIXHistograms
     */
    virtual Processor*  newProcessor() { return new EUTelAPIXHistograms ; }

    //! Default constructor
    EUTelAPIXHistograms() ;

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

bool hasMatchedHit(Track* fittrack);
//
void revertAlignment(double & x, double & y, double & z,  std::string	collectionName, LCEvent * ev);


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

    std::vector<double> _measuredX;
    std::vector<double> _measuredY;
    std::vector<double> _measuredZ;


    std::vector<double> _bgmeasuredX;
    std::vector<double> _bgmeasuredY;

    std::vector<double> _fittedX;
    std::vector<double> _fittedY;

    std::vector<double> _bgfittedX;
    std::vector<double> _bgfittedY;

    std::vector<float > _DUTalign;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Used to refer to histograms by their names, i.e. to recall
     *  a histogram pointer using histogram name.
     */

    std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;

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



	// some APIX stuff (histograms & arrays)
	// libov@mail.desy.de 05 August 2010
	static std::string	_totAPIXmatchedHistoName;
	static std::string	_totAPIXunmatchedHistoName;
	static std::string	_lv1APIXmatchedHistoName;
	static std::string	_lv1APIXunmatchedHistoName;

	std::vector<double> 	_totAPIX;
	std::vector<float>  	_lv1APIX;

	// some other relevant histograms
	static std::string	_hitsInMatchedClusterHistoName;
	static std::string	_hitsInUnmatchedClusterHistoName;

	std::vector<int>	  	_hitsInCluster;

	static std::string	_maxDifflv1MatchedHistoName;
	static std::string	_maxDifflv1UnmatchedHistoName;
	std::vector<float>	_maxDifflv1;

	// ---  V. Libov 13 July (17 August)----
	double _distMaxReference;
	int _referencePlaneID;
	int _referencePlaneIndex;
	bool _onlyIntimeTracks;
	// Implementing DUT plots in its local FoR, 18 August
    static std::string _EfficiencyXLOCALHistoName;
    static std::string _EfficiencyYLOCALHistoName;
    static std::string _EfficiencyXYLOCALHistoName;

	// measured hits in the local FoR of the DUT
    std::vector<double> _measuredXLOCAL;
    std::vector<double> _measuredYLOCAL;

	int		_indexDUT;
	double	_xPitch, _yPitch, _rot00, _rot01, _rot10, _rot11;

	bool	_noHitYet;
	int		_eventsWithNoHit;
	
	double	_transShiftX, _transShiftY;

	void	getTransformationShifts();
	double a,b,c,d;

	// -- more DUT plots, 24 August 2010, libov@mail.desy.de
    static std::string _MeasuredXLOCALHistoName;
    static std::string _MeasuredYLOCALHistoName;
    static std::string _MeasuredXYLOCALHistoName;

    static std::string _FittedXLOCALHistoName;
    static std::string _FittedYLOCALHistoName;
    static std::string _FittedXYLOCALHistoName;

	static std::string	_MatchedClusterSizeXHistoName;
	static std::string	_MatchedClusterSizeYHistoName;

	static std::string	_UnmatchedClusterSizeXHistoName;
	static std::string	_UnmatchedClusterSizeYHistoName;

	static std::string	_ChargeSharingProbXHistoName;
	static std::string	_ChargeSharingProbYHistoName;
	static std::string	_ChargeSharingProbXYHistoName;

	std::vector<int>	  	_clusterSizeX;
	std::vector<int>	  	_clusterSizeY;

	// and even more (06 September 2010 libov@mail.desy.de)
	static std::string	_totSinglePixelClustersXHistoName;
	static std::string	_totSinglePixelClustersYHistoName;
	static std::string	_totSinglePixelClustersXYHistoName;

	static std::string	_totAllClustersXYHistoName;

	std::vector<float>	_allDUTsmeasuredX;
	std::vector<float>	_allDUTsmeasuredY;
	std::vector<float>	_allDUTslv1;
	std::vector<int>		_allDUTssensorID;

	int		_lv1;

	void fillAPIXhits(LCCollection* hitcol, LCCollection* original_zsdata);

	static std::string	_MatchedHitsHistoName;
	static std::string	_MatchedHitsVSEventHistoName;
	static std::string	_NumberOfFittedTracksHistoName;

	static std::string	_NumberOfMatchedReferencesHistoName;
	static std::string	_Lv1OfFirstMatchedReferenceHistoName;

	// 18 January 2011
	StringVec		_alignmentCollectionNames;
	StringVec		_alignmentCollectionSuffixes;
	double			_beta;
	void TransformToLocalFrame(double & x, double & y, double & z, LCEvent * ev);

	// 21 january
	// corrected for non-normal sensors
	std::vector<double> _fittedXcorr;
	std::vector<double> _fittedYcorr;
	std::vector<double> _fittedZcorr;
	void getTrackImpactPoint(double & x, double & y, double & z, Track * tr, LCEvent * ev);
	int	_indexDUTneighbour;
	double	_zDUTneighbour;


#endif

  } ;


  //! A global instance of the processor.
  EUTelAPIXHistograms aEUTelAPIXHistograms ;


}

#endif



