// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelTestFitter_h
#define EUTelTestFitter_h 1

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h> 
#include "lcio.h"

// AIDA includes <.h>
#ifdef MARLIN_USE_AIDA
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {


  //! Analytical track fitting processor for EUDET Telescope
  /*! This processor was designed for fitting tracks to hits reconstructed in
   * the telescope sensor planes. Analytical approach is used, taking
   * into account multiple scattering in the telescope planes.
   * 
   *  \par Method
   *  Track fitting is performed separately in XZ and YZ planes (Z
   *  is defined along the beam axis direction). Track position in each
   *  telescope plane is found by solving matrix equation resulting from
   *  Chi^2 minimum condition. The following approximation is used: 
   *    \li all telescope planes are parallel to each other
   *    \li the incoming beam is perpendicular to the telescope planes
   *    \li the incoming beam has a small angular spread 
   *    \li particle scattering angles in subsequent telescope layers
   *        are also small
   *    \li thicknesses of all material layers are very small compared
   *        to the distances between planes 
   *    \li particle energy losses in telescope layers can be neglected
   * 
   *
   * \par Algorithm  
   * \li Read measured track points from input \c TrackerHit collection
   *     and copy to local tables
   * \li Prepare lists of hits for each active sensor plane, apply
   *     cuts on hit positions, if required
   * \li Count hit numbers, return if not enough planes fired
   * \li Calculate number of fit hypothesis (including missing hit possibility)
   * \li Search the list of fit hypotheses to find the one with best
   *     Chi^2 (including ``penalties'' for missing hits or skipped planes)
   * \li Accept the fit if Chi^2 is below threshold
   * \li Write fitted track to output \c Track collection; measured
   * particle positions corrected for alignment and fitted positions
   * are also written out as \c TrackerHit collections 
   * \li Remove accepted track hits from hit list and repeat procedure 
   *
   * \par Geometry description
   * This version of the processor does use GEAR input!
   * However, corrections to the geometry description (alignment,
   * removing layers from the fit) can be applied with dedicated parameters
   * (see below).
   *
   * 
   * \par Output
   * Fitted particle positions in all telescope planes are stored as 
   * \c TrackerHit collection. If required, measured
   * particle positions corrected for alignment can also be stored as
   * a separate  \c TrackerHit collection. In addition fit results are
   * written in a \c Track collection. Following \c Track variables
   * are filled:  
   *  \li Chi2 of the fit 
   *  \li number of measured hits used in the track fit (as Ndf)
   *  \li reconstructed position at DUT  (as a track reference point)
   *  \li vector of hits (fitted particle positions in all planes)  
   * 
   * \par Main algorithm parameters
   * \param InputCollectionName  Name of the input TrackerHit collection
   * \param CorrectedHitCollectionName Name of the collection for storing
   *        corrected particle positions in telescope planes (hits),
   *        i.e. positions after alignment, as used in the fit
   * \param OutputHitCollectionName Name of the output collection of
   *        fitted particle positions in telescope planes (hits)
   * \param OutputTrackCollectionName Name of the output Track collection
   * \param InputHitsInTrack Flag for storing input (measured) hits in track.
   * \param OutputHitsInTrack Flag for storing output (fitted) hits in track.  
   *        Input and output hits can be distinguished by looking into
   *        hit type (type <=31 for measured hits, type >=32 for fitted).
   *
   * \param AllowMissingHits Allowed number of hits missing in the track
   *        (sensor planes without hits or with hits removed from
   *        given track) 
   * \param MissingHitPenalty  "Penalty" added to track Chi^2 for each
   *        missing hit (no hits in given layer).
   * \param AllowSkipHits Allowed number of hits removed from the track
   *        (because of large Chi^2 contribution)
   * \param SkipHitPenalty  "Penalty" added to track Chi^2 for each hit
   *        removed from the track because of large Chi^2 contribution.
   * \param AllowAmbiguousHits Allow same hit to be used in more than one.
   *        Significantly improves algorithm performance.
   *
   * \param MaxPlaneHits Maximum number of hits considered per
   *        plane. The algorithm slows down if this number is
   *        too large. However, the real limitation comes from
   *        numerical precision. Maximum number is 34 for 6 planes
   *        used in the fit, 72 for 5 planes, 214 for 4 planes.
   *
   * \param UseNominalResolutio Flag for using nominal sensor resolution
   *        (as given in geometry description) instead of hit position
   *        errors. 
   *
   * \param UseDUT Flag for including DUT measurement in the track fit.
   *
   * \param Ebeam Beam energy in [GeV], needed to estimate multiple
   *        scattering. 
   *
   * \param UseBeamConstraint Flag for using beam direction constraint
   *        in the fit. Can improve the fit, if beam angular spread is
   *        small. 
   * \param BeamSpread Assumed angular spread of the beam [rad]
   *
   * \param SearchMultipleTracks Flag for searching multiple tracks in
   *        events with multiple hits 
   *
   * \param Chi2Max Maximum Chi2 for accepted track fit.
   *
   * \par Performance control parameters
   * \param DebugEventCount      Print out debug and information
   *        messages only for one out of given number of events. If
   *        zero, no debug information is printed. 
   *
   * \param SkipLayerIDs Ids of layers which are described in GEAR but
   *        should not be included in the fit. Can be used to remove
   *        layers in front of and behind the telescope, which do not
   *        influence the fit, but can slow down the algorithm
   *        (increase fit matrix size). 
   *
   * \param AlignLayerIDs Ids of layers for which alignment corrections
   *        should be applied
   * \param AlignLayerShiftX Shifts in X, which should be applied to
   *        correct alignment of these layers.
   * \param AlignLayerShiftY Shifts in Y, which should be applied to
   *        correct alignment of these layers.
   * \param AlignLayerRotZ Rotation around Z (beam) axis, which should 
   *        be applied to correct alignment of these layers.
   *
   * \param WindowLayerIDs Ids of layers for which position cuts are
   *        defined. Only hits inside the defined "window" are accepted
   * \param WindowMinX   Lower window edge in X
   * \param WindowMaxX   Upper window edge in X
   * \param WindowMinY   Lower window edge in Y
   * \param WindowMaxY   Upper window edge in Y
   *
   * \param MaskLayerIDs Ids of layers for which position cuts are
   *        defined. Only hits outside the defined "mask" are accepted
   * \param MaskMinX   Lower window edge in X
   * \param MaskMaxX   Upper window edge in X
   * \param MaskMinY   Lower window edge in Y
   * \param MaskMaxY   Upper window edge in Y
   *
   * \param HistoInfoFileName Name of the histogram information file.
   *        Using this file histogram parameters can be changed without 
   *        recompiling the code.
   *
   * \todo
   *  \li Interface to LCCD (alignment)
   *
   * \author A.F.Zarnecki, University of Warsaw
   * @version $Id: EUTelTestFitter.h,v 1.12 2008-01-27 22:55:31 zarnecki Exp $
   * \date 2007.10.30
   *
   */ 


  class EUTelTestFitter : public marlin::Processor {
  
  public:

  
     
    //! Returns a new instance of EUTelTestFitter
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelTestFitter
     */
    virtual Processor*  newProcessor() { return new EUTelTestFitter ; }
  
    //! Default constructor 
    EUTelTestFitter() ;
  
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
    // Fitting functions

    //! Find track in XZ and YZ 
    /*! Fit track in two planes (XZ and YZ) by solving two matrix
     * equations and calculate Chi^2
     */
    double MatrixFit();

    //! Find track in XZ and YZ assuming nominal errors 
    /*! Fit track in two planes: XZ and YZ. When nominal position errors
     * are used, only one matrix equation has to be solved and the
     * inverse matrix can be applied to the second equation.
     */
    double SingleFit();


    //! Find track in all planes assuming nominal errors 
    /*! Fit track in two planes: XZ and YZ. When nominal position errors
     * are assumed and hits are found in all sensor planes, same inverse
     * matrix can be used for all events.
     */
    double NominalFit();

    //! Fit particle track in one plane (XZ or YZ) 
    int DoAnalFit(double * pos, double *err);

    //! Calculate Chi^2 of the fit
    /*! Calculate Chi^2 of the fit taking into account measured particle
     *  positions in X and Y and fitted scattering angles in XZ and YZ
     *  planes 
     */
    double GetFitChi2();

    //! Solve matrix equation
    int GaussjSolve(double * alfa, double * beta, int n);


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


    // Global processor parameters
    // Parameter documentation is already included above

    int _debugCount ;

    std::string _inputColName ;

    std::string _outputTrackColName ;

    std::string _correctedHitColName ;

    std::string _outputHitColName ;

    bool _InputHitsInTrack;

    bool _OutputHitsInTrack;

    std::vector<int >   _SkipLayerIDs;

    std::vector<int >   _AlignLayerIDs;
    std::vector<float > _AlignLayerShiftX;
    std::vector<float > _AlignLayerShiftY;
    std::vector<float > _AlignLayerRotZ;

    std::vector<int >   _WindowLayerIDs;
    std::vector<float > _WindowMinX;
    std::vector<float > _WindowMaxX;
    std::vector<float > _WindowMinY;
    std::vector<float > _WindowMaxY;

    std::vector<int >   _MaskLayerIDs;
    std::vector<float > _MaskMinX;
    std::vector<float > _MaskMaxX;
    std::vector<float > _MaskMinY;
    std::vector<float > _MaskMaxY;

    // Parameters of hit selection algorithm

    int _allowMissingHits;
    int _allowSkipHits;
    int _maxPlaneHits;

    bool _searchMultipleTracks;

    bool _allowAmbiguousHits;

    // Parameters of fitting algorithm

    double _missingHitPenalty;
    double _skipHitPenalty;
    double _chi2Max ;

    bool   _useNominalResolution ;

    bool   _useDUT ;

    bool   _useBeamConstraint ;
    double _beamSpread;

    double _eBeam ;

    // Setup description

    int _nTelPlanes;
    int _nActivePlanes;
    int _iDUT;
  
    int * _planeSort;
    int * _planeID;
    double * _planeShiftX;
    double * _planeShiftY;
    double * _planeRotZ;
    double * _planePosition;
    double * _planeThickness;
    double * _planeX0;
    double * _planeResolution;
    bool   * _isActive;

    std::vector<int> * _planeWindowIDs;
    std::vector<int> * _planeMaskIDs;

    // Internal processor variables
    // ----------------------------

    int _nRun ;
    int _nEvt ;

    // Arrays for selecting different hit combinations

    int * _planeHits;
    int * _planeChoice;
    int * _planeMod;

    // Fitting algorithm arrays

    double * _planeX  ;
    double * _planeEx ;
    double * _planeY  ;
    double * _planeEy ;

    double * _planeDist ;
    double * _planeScat ;

    double * _fitX  ;
    double * _fitEx ;
    double * _fitY  ;
    double * _fitEy ;

    double * _fitArray ;
    double * _nominalFitArray ;
    double * _nominalError ;

 #ifdef MARLIN_USE_AIDA
    //! AIDA histogram map
    /*! Used to refer to histograms by their names, i.e. to recall 
     *  a histogram pointer using histogram name.
     */

    std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;

    // Chi2 histogram names
    static std::string _linChi2HistoName;
    static std::string _logChi2HistoName;
    static std::string _firstChi2HistoName;
    static std::string _bestChi2HistoName;
    static std::string _fullChi2HistoName;

    // Number of reconstructed tracks histogram name
    static std::string _nTrackHistoName;

    // Number of hits histogram names
    static std::string _nAllHitHistoName;
    static std::string _nAccHitHistoName;

    static std::string _nHitHistoName;
    static std::string _nBestHistoName;

    static std::string _hitAmbiguityHistoName;

#endif 

 } ;

  
  //! A global instance of the processor.
  EUTelTestFitter aEUTelTestFitter ;


}

#endif



