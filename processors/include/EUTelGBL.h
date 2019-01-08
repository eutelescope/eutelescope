/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELGBL_H
#define EUTELGBL_H

// built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelUtility.h"
#include "EUTelTripletGBLUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

//for gbl::MilleBinary
#include "include/MilleBinary.h"

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <limits>

namespace eutelescope {

  class EUTelGBL : public marlin::Processor {

    public:
      //! Returns a new instance of EUTelGBL
      /*! This method returns a new instance of this processor.  It is
       *  called by Marlin execution framework and it shouldn't be
       *  called/used by the final user.
       *
       *  @return a new EUTelGBL.
       */
      virtual Processor * newProcessor() { return new EUTelGBL; }

      //! Default constructor
      EUTelGBL ();

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
      void bookHistos(std::vector<int> const & );

    protected:
      static int const NO_PRINT_EVENT_COUNTER = 3;
    
      //! Ordered sensor ID
      /*! Within the processor all the loops are done up to _nPlanes and
       *  according to their position along the Z axis (beam axis).
       */
      std::vector<int> _sensorIDVec;

      //! TrackerHit collection name
      /*! Input collection with hits.
      */
      std::vector<std::string> _hitCollectionName;

      //general parameters
      EUTelTripletGBLUtility gblutil;
      double _zMid;
      double _eBeam;
      double _kappa;

      //plane parameters
      std::vector<double> _planePosition;
      std::vector<double> _planeRadLength;
      std::vector<Eigen::Vector2d> _planeWscatSi;
      std::vector<Eigen::Vector2d> _planeWscatAir;
      std::vector<Eigen::Vector2d> _planeMeasPrec;
      std::vector<float> _xResolutionVec;
      std::vector<float> _yResolutionVec;
      int _SUT_ID;
      std::vector<int>_DUT_IDs;
      std::vector<float> _dutCuts;
      std::vector<int> _excludedPlanes;
      std::vector<int> _fixedPlanes;
      int _requiredPlane;
        
      //alignment: alignModes and fixed settings
      Utility::alignMode _alignMode;
      std::string _alignModeString;
      std::vector<int> _FixedXShift;
      std::vector<int> _FixedYShift;
      std::vector<int> _FixedZShift;
      std::vector<int> _FixedXRot;
      std::vector<int> _FixedYRot;
      std::vector<int> _FixedZRot;
      
      //triplet parameters
      double _upstreamTriplet_ResCut;
      double _downstreamTriplet_ResCut;
      double _upDownTripletMatchCut;
      double _upstreamTriplet_SlopeCut;
      double _downstreamTriplet_SlopeCut;
      std::vector<int> _upstreamTriplet_IDs;
      std::vector<int> _downstreamTriplet_IDs;
      int _lastUpstreamSensorID;
      std::map<int, bool> _isSensorUpstream;

      //track parameters
      double _chi2Cut;   
      int _maxTrackCandidates;
      int _maxTrackCandidatesTotal;

      //MILLEPEDE
      std::string _binaryFilename;
      std::string _pedeSteerfileName;
      std::unique_ptr<gbl::MilleBinary> milleAlignGBL;

      //statistics
      int _iRun;
      int _iEvt;
      int _printEventCounter;
      size_t _nMilleDataPoints;
      size_t _nMilleTracks;
      size_t _nPlanes;

      //flags
      int _performAlignment;
      int _suggestAlignmentCuts;
      int _dumpTracks;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //histograms: hits, triplets and tracks
    AIDA::IHistogram1D * hist1D_nTelescopeHits;
    AIDA::IHistogram1D * hist1D_nUpstreamTriplets;
    AIDA::IHistogram1D * hist1D_nDownstreamTriplets;
    AIDA::IHistogram1D * hist1D_nTracksPerEvent;

    //histograms: GBLFit
    AIDA::IHistogram1D * hist1D_gblNdfAlign;
    AIDA::IHistogram1D * hist1D_gblChi2Align;
    AIDA::IHistogram1D * hist1D_gblProbAlign;
    std::vector<AIDA::IHistogram1D*> hist1D_gblAngleX;
    std::vector<AIDA::IHistogram1D*> hist1D_gblAngleY;
    std::vector<AIDA::IHistogram1D*> hist1D_gblResidX;
    std::vector<AIDA::IHistogram1D*> hist1D_gblResidY;
    std::vector<AIDA::IHistogram1D*> hist1D_gblPullX;
    std::vector<AIDA::IHistogram1D*> hist1D_gblPullY;
    std::vector<AIDA::IHistogram1D*> hist1D_gblKinkX;
    std::vector<AIDA::IHistogram1D*> hist1D_gblKinkY;

	//histograms: SUT specific
    AIDA::IProfile2D* profile2D_gblSUTKinkXvsXY;
    AIDA::IProfile2D* profile2D_gblSUTKinkYvsXY;
#endif
  };	

  //! A global instance of the processor
  EUTelGBL gEUTelGBL;
}
#endif
#endif
