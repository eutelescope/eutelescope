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
#include "lcio.h"
#include <string>

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
   * \li Read measured track points from input TrackerHit collection
   *     and copy to local tables
   * \li Prepare lists of hits for each active sensor plane
   * \li Count hit numbers, return if not enough planes fired
   * \li Calculate number of fit hypothesis (including missing hit possibility)
   * \li Fit each hypotheses and calculate Chi^2 (including
   * ``penalties'' for missing hits or skipped planes)
   * \li Select the Chi^2
   * \li Write fitted track position at DUT to output TrackerHit collection
   * \li Remove best track hits from hit list and repeat procedure 
   *
   * \par Geometry description
   * This version of the processor does not use GEAR input yet!
   * Needed geometry information are stored in a dedicated ASCII file.
   * The format of the file is as follows. \n
   * First line consists of two int numbers: 
   * \li number of telescope layers N (sensors and passive layers taken
   * into account in the fit) and 
   * \li position of DUT  i_DUT ( i_DUT = 1 ... N ). i_DUT=0 corresponds
   * to telescope without DUT.
   * 
   * Following N lines include description of telescope planes.
   * Planes have to be ordered in position along the beam line ! 
   * For each plane following details have to be given:
   * \li plane alignment correction in horizontal direction in mm (float)
   * \li plane alignment correction in vertical direction in mm (float)
   * \li position of the plane along beam axis in mm (float)
   * \li thickness of the plane in mm (float)
   * \li radiation length in the plane material in mm (float)
   * \li flag indicating sensitive planes (int; 1 for sensitive, 0 for passive plane)
   * \li nominal position resolution of sensitive planes in mm (float)
   *
   * 
   * \par Output
   * Fitted particle positions in all telescope planes are stored as 
   * \c TrackerHit collection. In addition fit results are written in a
   * \c Track collection. Following \c Track variables are filled: 
   *  \li Chi2 of the fit 
   *  \li number of measured hits used in the track fit (as Ndf)
   *  \li reconstructed position at DUT  (as a track reference point)
   *  \li vector of hits (fitted particle positions in all planes)  
   *
   *
   * \param InputCollectionName  Name of the input  TrackerHit collection
   * \param OutputHitCollectionName Name of the output collection of
   *        fitted particle positions in telescope planes (hits)
   * \param OutputTrackCollectionName Name of the output Track collection
   * \param DebugEventCount      Print out debug and information
   * messages only for one out of given number of events. If zero, no
   * debug information is printed. 
   * \param GeometryFileName     Name of the geometry description
   * file. This version of the processor does not use GEAR input yet!
   * Needed geometry information are stored in a dedicated ASCII file.
   * \param AllowMissingHits Allowed number of hits missing in the track
   * (sensor planes without hits or with hits removed from given track)
   * \param AllowSkipHits Allowed number of hits removed from the track
   * (because of large Chi^2 contribution)
   * \param MaxPlaneHits Maximum number of hits considered per
   * plane. The algorithm becomes very slow if this number is too large.
   * \param MissingHitPenalty  "Penalty" added to track Chi^2 for each
   * missing hit (no hits in given layer).
   * \param SkipHitPenalty  "Penalty" added to track Chi^2 for each hit
   * removed from the track because of large Chi^2 contribution.
   * \param Chi2Max Maximum Chi2 for accepted track fit.
   * \param UseNominalResolutio Flag for using nominal sensor resolution
   * (as given in geometry description) instead of hit position errors.
   * \param UseDUT Flag for including DUT measurement in the track fit.
   * \param Ebeam Beam energy in [GeV], needed to estimate multiple scattering.
   * \param UseBeamConstraint Flag for using beam direction constraint
   * in the fit. Can improve the fit, if beam angular spread is small.
   * \param BeamSpread Assumed angular spread of the beam [rad]
   * \param SearchMultipleTracks Flag for searching multiple tracks in
   * events with multiple hits 
   *
   * \warning
   * Algorithm is very fast when there are only single (or double) hits
   * in each sensor layer. With increasing number of hits per layer (HpL)
   * fitting time increases like HpL^N, where N is number of sensors.
   * To protect from "freezing"  \c MaxPlaneHits parameter is introduced
   * (default value is 5 and should rather not be increased).
   *
   * \todo
   *  \li Interface to GEAR geometry description
   *  \li More detailed track parameter output (currently only reconstructed
   *      hit positions are stored)
   *  \li Implement new track selection algorithm, which would avoid
   * testing all track hypothesis (to improve efficiency for events with
   * many hits)
   *
   * \author A.F.Zarnecki, University of Warsaw
   * @version $Id: EUTelTestFitter.h,v 1.2 2007-06-26 16:19:19 zarnecki Exp $
   * \date 2007.06.04
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

    // Global processor parameters
    // Parameter documentation is already included above

    int _debugCount ;
   
    std::string _inputColName ;

    std::string _outputTrackColName ;

    std::string _outputHitColName ;

    // Parameters of hit selection algorithm

    int _allowMissingHits;
    int _allowSkipHits;
    int _maxPlaneHits;

    bool _searchMultipleTracks;


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

    std::string _geometryFileName;

    int _nTelPlanes;
    int _nActivePlanes;
    int _iDUT;
  
    double * _planeShiftX;
    double * _planeShiftY;
    double * _planePosition;
    double * _planeThickness;
    double * _planeX0;
    double * _planeResolution;
    bool   * _isActive;

    // Internal processor variables
    // ----------------------------

    int _nRun ;
    int _nEvt ;

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
  } ;

  
  //! A global instance of the processor.
  EUTelTestFitter aEUTelTestFitter ;


}

#endif



