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
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

// define a specific verbosity level for Millepede
namespace streamlog {

  DEFINE_STREAMLOG_LEVEL( TESTFITTERMESSAGE,  "TESTFITTERMESSAGE", message_base_level - 50 + 3 , STREAMLOG_MESSAGE_ACTIVE ) 

}

namespace eutelescope {


  //! Analytical track fitting processor for EUDET Telescope
  /*! This processor was designed for fitting tracks to hits
   * reconstructed in the telescope sensor planes. Analytical approach
   * is used, taking into account multiple scattering in the telescope
   * planes. However, as there are usually multiple hits in each
   * plane, main task of the processor is to look for the best track
   * candidate considering all fit possibilities.
   *
   *  \par Fit method
   *  Track fitting is performed separately in XZ and YZ planes (Z
   *  is defined along the beam axis direction). Track position in each
   *  telescope plane is found by solving matrix equation resulting from
   *  \f$ \chi^{2} \f$ minimum condition. The following approximation is used:
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
   * \par Track finding algorithm
   * \li Read measured track points from input \c TrackerHit collection
   *     and copy to local tables
   * \li Prepare lists of hits for each active sensor plane, apply
   *     cuts on hit positions, if required
   * \li Count hit numbers, return if not enough planes fired
   * \li Calculate number of fit hypothesis (including missing hit possibility)
   * \li Search the list of fit hypotheses to find the one with best
   *     \f$ \chi^{2} \f$ (including ``penalties'' for missing hits or
   *     skipped planes)
   * \li Accept the fit if \f$ \chi^{2} \f$ is below threshold
   * \li Write fitted track to output \c Track collection; measured
   *     particle positions corrected for alignment and fitted positions
   *     are also written out as \c TrackerHit collections
   * \li Remove accepted track hits from hit list and repeat procedure
   *
   * \par Geometry description
   * The processor does use GEAR input for telescope layers and DUT
   * description. Corrections or modification to the geometry
   * description (alignment, removing layers from the fit, etc.) can
   * be applied with dedicated parameters (see below).
   *
   *
   * \par Output
   * Fitted particle positions in all telescope planes are stored as
   * \c TrackerHit collection. If required, measured
   * particle positions corrected for alignment can also be stored as
   * a separate  \c TrackerHit collection. In addition fit results are
   * written in a \c Track collection. Following \c Track variables
   * are filled:
   *  \li \f$ \chi^{2} \f$ of the fit
   *  \li number of measured hits used in the track fit (as Ndf)
   *  \li reconstructed position at DUT  (as a track reference point)
   *  \li vector of output hits (fitted particle positions in all planes)
   *  \li vector of input hits (corrected particle positions in fired planes)
   *
   *
   * \section parameters Control parameters
   * Below, parameters defining performance of the algorithm are
   * described. Some suggestion for the optimal choice of parameters
   * are given in the next section.
   *
   *
   * \par Main algorithm parameters
   *  Following parameters define input and output from the processor.
   *
   * \param InputCollectionName  Name of the input TrackerHit collection
   * \param OutputTrackCollectionName Name of the output Track
   *        collection (main output of the processor).
   * \param InputHitsInTrack Flag for storing input (measured,
   *        corrected for alignment) hits together with the track.
   * \param CorrectedHitCollectionName Name of the collection for storing
   *        corrected particle positions in telescope planes (hits),
   *        i.e. positions after alignment, as used in the fit
   * \param OutputHitsInTrack Flag for storing output (fitted) hits
   *        together with the track.
   *        Input and output hits can be distinguished by looking into
   *        hit type (type <=31 for measured hits, type >=32 for fitted).
   * \param OutputHitCollectionName Name of the output collection of
   *        fitted particle positions in telescope planes (hits)
   * \param DebugEventCount      Print out debug and information
   *        messages only for one out of given number of events. If
   *        zero, no debug information is printed.
   * \param HistoInfoFileName Name of the histogram information file.
   *        Using this file histogram parameters can be changed without
   *        recompiling the code.
   * \param Ebeam Beam energy in [GeV], needed to estimate multiple
   *        scattering.
   *
   *
   *
   * \par Geometry description parameters
   * Following parameters can be used to adjust geometry description
   *
   * \param SkipLayerIDs Ids of layers which are described in GEAR but
   *        should not be included in the fit. Can be used to remove
   *        layers in front of and behind the telescope, which do not
   *        influence the fit, but can slow down the algorithm
   *        (increase fit matrix size).
   *
   * \param PassiveLayerIDs Ids of layers which are described as
   *        active layers in GEAR but should be treated as passive in
   *        the fit (their data should be ignored).
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
   *
   * \param WindowLayerIDs Ids of layers for which position cuts are
   *        defined. Only hits inside the defined "window" are accepted
   * \param WindowMinX   Lower window edge in X
   * \param WindowMaxX   Upper window edge in X
   * \param WindowMinY   Lower window edge in Y
   * \param WindowMaxY   Upper window edge in Y
   *
   *
   * \param MaskLayerIDs Ids of layers for which position cuts are
   *        defined. Only hits outside the defined "mask" are accepted
   * \param MaskMinX   Lower window edge in X
   * \param MaskMaxX   Upper window edge in X
   * \param MaskMinY   Lower window edge in Y
   * \param MaskMaxY   Upper window edge in Y
   *
   * \par Fit performance parameters
   * \param MaxPlaneHits Maximum number of hits considered per
   *        plane. The algorithm slows down if this number is
   *        too large. However, the real limitation comes from
   *        numerical precision. To find the best track fit hypothesis
   *        have to be numbered and long integer number is used for this
   *        purpose. However this allows for maximum of \f$ 2^{31} \f$
   *        hypothesis only. To avoid this limit the number of hits in
   *        single plane has to be constrained. Maximum number is 34
   *        for 6 planes used in the fit, 72 for 5 planes, 214 for 4 planes.
   *        Never use higher values!
   *
   * \param Chi2Max Maximum \f$ \chi^{2} \f$ for accepted track fit.
   *
   * \param SearchMultipleTracks Flag for searching multiple tracks in
   *        events with multiple hits. If false, only best (lowest
   *        \f$ \chi^{2} \f$) track is taken.
   *
   * \param AllowAmbiguousHits Allow same hit to be used in more than
   *        one track. Significantly improves algorithm
   *        performance. However, this option can only be used when
   *        missing hits are not allowed (\e  AllowMissingHits set to 0 )
   *
   * \param UseNominalResolution Flag for using nominal sensor resolution
   *        (as given in geometry description) instead of hit position
   *        errors. Improves tracking performance.
   *
   * \param UseDUT Flag for including DUT measurement in the track fit.
   *
   * \param UseBeamConstraint Flag for using beam direction constraint
   *        in the fit. Can improve the fit, if beam angular spread is
   *        small. Improves track searching for multiple hits.
   *
   * \param BeamSpread Assumed angular spread of the beam [rad]
   *
   * \param AllowMissingHits Allowed number of hits missing in the track
   *        (sensor planes without hits or with hits removed from
   *        given track)
   * \param MissingHitPenalty  "Penalty" added to track
   *        \f$ \chi^{2} \f$ for each
   *        missing hit (when no hit is left in active layer).
   * \param AllowSkipHits Allowed number of hits removed from the track
   *        (because of large \f$ \chi^{2} \f$ contribution)
   * \param SkipHitPenalty  "Penalty" added to track
   *        \f$ \chi^{2} \f$ for each hit removed from the track
   *        because of large \f$ \chi^{2} \f$ contribution.
   *
   * \section performance Performance issues
   * As described above, if multiple hits are found in telescope
   * layers or hit rejection is allowed, the algorithm checks all hits
   * selection possibilities (all track hypothesis). This task is
   * optimized to a large extent, but still track finding can be slow
   * for large multiplicities. Here are some suggestions on how to
   * improve performance.
   *
   *  \li Remove addition layers from geometry description. At high
   *      energies, when multiple scattering can be neglected, only
   *      active telescope planes and DUTs are relevant. At low
   *      energies you can still remove planes which are in front of
   *      the first active layer and behind the last one. You can also
   *      consider removing thin passive layers inside the
   *      telescope. Number of layers determines the order of matrix
   *      equation which has to be solved for each track.
   *
   *  \li Use nominal plane resolutions instead of cluster position
   *      errors (set \e UseNominalResolution to \e true ). For full
   *      tracks (without missing hits) matrix inversion is done only
   *      once and not for each track hypothesis.
   *
   *  \li Use beam constraint (set \e UseBeamConstraint to \e true ),
   *      even if beam spread is large. With beam
   *      constraint first two hits are sufficent to recognize bad track
   *      hypothesis. Without beam constraint at least 3 hits are needed.
   *      However, to use beam constraint telescope layers have to be
   *      aligned w.r.t. the beam direction (beam has to be perpendicular to
   *     telescope layers).
   *
   * \li Do not allow for missing hits (set \e AllowMissingHits to 0).
   *     This reduces number of fit hypothesis and improves fit performance.
   *
   * \li Allow same hit to be used in more than one track (set
   *    \e AllowAmbiguousHits to \e true ). If each hit can be used in
   *    one track only, track finding has to be repeated many times,
   *    each time finding the best track hypothesis and then removing
   *    corresponding hits before looking for the next track. If
   *    ambiguity is allowed only one loop over all hypothesis is
   *    needed. Therefor this option can significantly improve algorithm
   *    performance at large multiplicities. The influence on the fit
   *    results is negligible, as probability of more than one track
   *    matching given hit is very, very small.
   *
   * \li Limit number of hits per plane. This should \b not be done by
   *    using \e MaxPlaneHits parameter, as it would bias plane
   *    efficiency calculation. Best way is to define position window
   *    in each active plane (see geometry description parameters
   *    above).
   *
   * \todo
   *  \li Interface to LCCD (alignment)
   *
   * \author A.F.Zarnecki, University of Warsaw, zarnecki@fuw.edu.pl
   * @version $Id: EUTelTestFitter.h,v 1.16 2008-10-04 12:28:53 bulgheroni Exp $
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
     * equations and calculate \f$ \chi^{2} \f$
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

    //! Calculate \f$ \chi^{2} \f$ of the fit
    /*! Calculate \f$ \chi^{2} \f$ of the fit taking into account measured particle
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
    std::vector<int >   _PassiveLayerIDs;

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

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
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



