// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/TrackerHitImpl.h>

//for gbl::MilleBinary
#include "include/MilleBinary.h"

#include "EUTelUtility.h"
#include "EUTelTripletGBLUtility.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {

  class EUTelAlignGBL : public marlin::Processor {

    public:
      //! Returns a new instance of EUTelAlignGBL
      /*! This method returns a new instance of this processor.  It is
       *  called by Marlin execution framework and it shouldn't be
       *  called/used by the final user.
       *
       *  @return a new EUTelAlignGBL.
       */
      virtual Processor * newProcessor() {
	return new EUTelAlignGBL;
      }

      //! Default constructor
      EUTelAlignGBL ();

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

      //! Called for first event per run
      /*! Reads hotpixel information from hotPixelCollection into hotPixelMap
       * to be used in the sensor exclusion area logic 
       */
      //DP virtual void  FillHotPixelMap(LCEvent *event);

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

      // parameters
      std::vector<int> _excludedPlanes;
      std::vector<int> _FixedPlanes;
      int _requiredPlane;

	  // Fixed directions for alignMode XYZShiftsXYZRot
      std::vector<int> _FixedXShift;
      std::vector<int> _FixedYShift;
      std::vector<int> _FixedZShift;
      std::vector<int> _FixedXRot;
      std::vector<int> _FixedYRot;
      std::vector<int> _FixedZRot;
      
      double _zMid;
      double _eBeam;

      double _upTriResCut;
      double _downTriResCut;
      double _upDownTripletMatchCut;
      double _upSlopeCut;
      double _downSlopeCut;
      std::vector<float> _DUTCuts;
      double _chi2cut;
      
      double _kappa;

      int _maxTrackCandidates;
      int _maxTrackCandidatesTotal;

      std::string _binaryFilename;

      Utility::alignMode _alignMode;
      std::string _alignModeString;

      std::vector<int> _FixParameter;

      int _performAlignment;
      std::string _pedeSteerfileName;

      EUTelTripletGBLUtility gblutil;

      //! Run number
      int _iRun;

      //! Event number
      int _iEvt;

      //! counter for printed events (for debugging)
      int _printEventCounter;

      // Statistics
      size_t _nMilleDataPoints;
      size_t _nMilleTracks;
      size_t _nPlanes;


      //UpstreamTriplet
      std::vector<int> _upstream_triplet_ids;
      //DownstreamTriplet
      std::vector<int> _downstream_triplet_ids;
      //LastUpstreamSensor
      int _last_upstream_sensor;
      //ResolutionX
      std::vector<float> _x_resolution_vec;
      //ResolutionY
      std::vector<float> _y_resolution_vec;

      std::map<int, bool> _is_sensor_upstream;

      std::vector<int>_dut_ids;

      std::vector<double> _planePosition;
      std::vector<double> _planeRadLength;
      std::vector<Eigen::Vector2d> _planeWscatSi;
      std::vector<Eigen::Vector2d> _planeWscatAir;
      std::vector<Eigen::Vector2d> _planeMeasPrec;
      std::vector<int> indexconverter;
      std::unique_ptr<gbl::MilleBinary>  milleAlignGBL; // for producing MillePede-II binary file

    // definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    AIDA::IHistogram1D * nTelHits;
    AIDA::IHistogram1D * nUpstreamTriplets;
    AIDA::IHistogram1D * nDownstreamTriplets;

    //All Tracks
    AIDA::IHistogram1D * selxHistGBLAlign;
    AIDA::IHistogram1D * selyHistGBLAlign;
    AIDA::IHistogram1D * selaxHistGBLAlign;
    AIDA::IHistogram1D * selayHistGBLAlign;
    AIDA::IHistogram1D * seldxHistGBLAlign;
    AIDA::IHistogram1D * seldyHistGBLAlign;
    AIDA::IHistogram1D * selkxHistGBLAlign;
    AIDA::IHistogram1D * selkyHistGBLAlign;
    std::vector<AIDA::IHistogram1D*> seldxSensorHistGBLAlign;
    std::vector<AIDA::IHistogram1D*> seldySensorHistGBLAlign;

    //Track Fit
    AIDA::IHistogram1D * gblndfHistGBLAlign;
    AIDA::IHistogram1D * gblchi2HistGBLAlign;
    AIDA::IHistogram1D * gblprbHistGBLAlign;

    std::vector<AIDA::IHistogram1D*> gblAxHist;
    std::vector<AIDA::IHistogram1D*> gblAyHist;

    std::vector<AIDA::IHistogram1D*> gblRxHist;
    std::vector<AIDA::IHistogram1D*> gblRyHist;

    std::vector<AIDA::IHistogram1D*> gblPxHist;
    std::vector<AIDA::IHistogram1D*> gblPyHist;

    std::vector<AIDA::IHistogram1D*> gblDxHist;
    std::vector<AIDA::IHistogram1D*> gblDyHist;

    std::vector<AIDA::IHistogram1D*> gblKinkXHist;
    std::vector<AIDA::IHistogram1D*> gblKinkYHist;

    AIDA::IHistogram1D * ntracksperevent;
#endif

  };

  //! A global instance of the processor
  EUTelAlignGBL gEUTelAlignGBL;

}
#endif
#endif
