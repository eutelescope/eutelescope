// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelGBLFitter_h
#define EUTelGBLFitter_h 1

#include "EUTelTripletGBLUtility.h"

#include <memory>
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
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

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// ROOT includes
#include <TMatrixD.h>
#include "TH1D.h"

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <deque>

//Eigen
#include <Eigen/Core>

namespace eutelescope {

  class EUTelGBLFitter : public marlin::Processor {

  public:

    //! Returns a new instance of EUTelGBLFitter
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelGBLFitter
     */
    virtual Processor*  newProcessor() { return new EUTelGBLFitter; }
    EUTelGBLFitter(const EUTelGBLFitter&); 
    void operator=(EUTelGBLFitter const&); 


    //! Default constructor
    EUTelGBLFitter();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution.
     *
     */
    virtual void init();

    //! Called for every run.
    /*!
     * @param run the LCRunHeader of the current run
     */
    virtual void processRunHeader( LCRunHeader* run );

    //! Called every event
    /*! This is called for each event in the file.
     *
     *  @param evt the current LCEvent event
     */
    virtual void processEvent( LCEvent * evt );

    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. It can be used to fill check
     *  plots. For the time being there is nothing to check and do in
     *  this slot.
     *
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check( LCEvent * evt );

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
    virtual void end();

  private:

    // Calculate Point-To-Point Jacobian Transport Matrix for distance "ds"
    //TMatrixD JacobianPointToPoint( double ds );

    //! Fill the telescope plane correlation plots:
    void TelescopeCorrelationPlots(std::vector<EUTelTripletGBLUtility::hit> &telescopehits);

    //! Find hit triplets from three telescope planes
    /*! This runs over all hits in the planes of the telescope and
     * tries to match triplets by comparing with the middle planes.
     * Two cut criteria can be set:
     * @param triplet_residual_cut Cut on the hit residual in the middle plane with respect to the triplet defined by first and last plane
     * @param triplet_angle_cut Cut on the triplet track angle
     *
     * @return a vector of found triplets among the given set of hits.
     */
    //void FindTriplets(std::vector<EUTelTripletGBLUtility::hit> &hits, unsigned int plane0, unsigned int plane1, unsigned int plane2, std::vector<EUTelGBLFitterUtility::triplet> &trip);

    //! Match the upstream and downstream triplets to tracks
    //void MatchTriplets(std::vector<EUTelGBLFitterUtility::triplet> &up, std::vector<EUTelGBLFitterUtility::triplet> &down, double z_match, std::vector<EUTelGBLFitterUtility::track> &track);

    //! Check isolation of triplet within vector of triplets
    //bool IsTripletIsolated(std::vector<EUTelGBLFitterUtility::triplet>::iterator it, std::vector<EUTelGBLFitterUtility::triplet> &trip, double z_match, double isolation = 0.3);

  
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

    std::string _inputCollectionTelescope;

    //Analysis parameters
    bool _isFirstEvent;
    int _printEventCounter;
    double _eBeam;
    int _nRun;
    int _nEvt;
    size_t _nPlanes;
    int _ngbl;
    int _dut_plane;
    double _eff_radius;
    double _kappa;
    double _aluthickum;
    bool _CSswitch;
    
    // Cuts for matching:
    double _track_match_cut;
    double _slope_cut;
    double _triplet_res_cut;
    double _probchi2_cut;

    // Partly outdated GEAR readings:
    std::vector<int> _sensorIDVec;
    std::vector<double>  _planePosition;
    std::vector<double> _planeRadLength;
    std::map<size_t, std::vector<double>> _planeResolutionX;
    std::map<size_t, std::vector<double>> _planeResolutionY;

	  std::vector<Eigen::Vector2d> _planeWscatSi;
	  std::vector<Eigen::Vector2d> _planeWscatAir;

    FloatVec _resolution;
    FloatVec _thickness;

    EUTelTripletGBLUtility gblutil;


      // definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  AIDA::IHistogram1D * nAllHitHisto;

  // Correlation plots for telescope planes
  AIDA::IHistogram1D * dx01Histo, * dy01Histo, * du01Histo, * dx02Histo, * dx03Histo, * dx04Histo, * dx05Histo, * dx12Histo, * dy12Histo, * du12Histo, * dx23Histo, * dy23Histo, * du23Histo, * dx34Histo, * dy34Histo, * du34Histo, * dx45Histo, * dy45Histo, * du45Histo;


  // triplets 0-1-2:
  AIDA::IHistogram1D * da02Histo;
  AIDA::IHistogram1D * db02Histo;

  AIDA::IProfile2D * dzcvsxy;
  AIDA::IProfile2D * z3vsxy;

  AIDA::IHistogram1D * tridxHisto;
  AIDA::IHistogram1D * tridyHisto;

  AIDA::IProfile1D * tridxvsx;
  AIDA::IProfile1D * tridxvsy;
  AIDA::IProfile1D * tridxvstx;
  AIDA::IProfile1D * tridxvsty;
  AIDA::IProfile1D * tridyvsx;
  AIDA::IProfile1D * tridyvsy;
  AIDA::IProfile1D * tridyvstx;
  AIDA::IProfile1D * tridyvsty;

  AIDA::IHistogram1D * tridx1Histo, * tridy1Histo, * tridx3Histo, * tridy3Histo, * tridx3bHisto, * tridy3bHisto, * tridx4Histo, * tridy4Histo, * tridx4bHisto, * tridy4bHisto, * tridx5Histo, * tridy5Histo, * tridx5bHisto, * tridy5bHisto, * trixHisto, * triyHisto, * tritxHisto, * trityHisto;
  AIDA::IHistogram2D * trixyHisto;

  AIDA::IHistogram1D * trixdutHisto; // at DUT
  AIDA::IHistogram1D * triydutHisto;
  AIDA::IHistogram2D * trixydutHisto;

  // driplets 3-4-5:

  AIDA::IHistogram1D * dx35Histo;
  AIDA::IHistogram1D * dy35Histo;

  AIDA::IHistogram1D * dridxHisto;
  AIDA::IHistogram1D * dridyHisto;
  AIDA::IHistogram1D * drixHisto;
  AIDA::IHistogram1D * driyHisto;
  AIDA::IHistogram2D * drixyHisto;
  AIDA::IHistogram1D * dritxHisto;
  AIDA::IHistogram1D * drityHisto;

  AIDA::IProfile1D * dridxvsx;
  AIDA::IProfile1D * dridxvsy;
  AIDA::IProfile1D * dridxvstx;
  AIDA::IProfile1D * dridxvsty;
  AIDA::IProfile1D * dridyvsx;
  AIDA::IProfile1D * dridyvsy;
  AIDA::IProfile1D * dridyvstx;
  AIDA::IProfile1D * dridyvsty;

  AIDA::IHistogram1D * bacsxaHisto;
  AIDA::IHistogram1D * bacdyaHisto;
  AIDA::IHistogram1D * bacsxcHisto;
  AIDA::IHistogram1D * bacdycHisto;
  AIDA::IHistogram1D * bacsxcqHisto;
  AIDA::IHistogram1D * bacdycqHisto;

  AIDA::IProfile1D * effix0;
  AIDA::IProfile1D * effiy0;
  AIDA::IProfile1D * effix1;
  AIDA::IProfile1D * effiy1;
  AIDA::IProfile1D * effix2;
  AIDA::IProfile1D * effiy2;
  AIDA::IProfile1D * effix3;
  AIDA::IProfile1D * effiy3;
  AIDA::IProfile1D * effix4;
  AIDA::IProfile1D * effiy4;
  AIDA::IProfile1D * effix5;
  AIDA::IProfile1D * effiy5;

  AIDA::IHistogram1D * ntriHisto;
  AIDA::IHistogram1D * ndriHisto;
  AIDA::IHistogram1D * nsixHisto;

  AIDA::IProfile2D * kinkpixvsxy;

  AIDA::IHistogram1D * clustersize0;
  AIDA::IHistogram1D * sixx0Histo;
  AIDA::IHistogram1D * sixy0Histo;
  AIDA::IHistogram1D * clustersize1;
  AIDA::IHistogram1D * sixx1Histo;
  AIDA::IHistogram1D * sixy1Histo;
  AIDA::IHistogram1D * clustersize2;
  AIDA::IHistogram1D * sixx2Histo;
  AIDA::IHistogram1D * sixy2Histo;
  AIDA::IHistogram1D * clustersize3;
  AIDA::IHistogram1D * sixx3Histo;
  AIDA::IHistogram1D * sixy3Histo;
  AIDA::IHistogram1D * clustersize4;
  AIDA::IHistogram1D * sixx4Histo;
  AIDA::IHistogram1D * sixy4Histo;
  AIDA::IHistogram1D * clustersize5;
  AIDA::IHistogram1D * sixx5Histo;
  AIDA::IHistogram1D * sixy5Histo;

  AIDA::IHistogram2D * sixxylkHisto;

  AIDA::IHistogram1D * derxtiltHisto;
  AIDA::IHistogram1D * derytiltHisto;
  AIDA::IHistogram1D * derxturnHisto;
  AIDA::IHistogram1D * deryturnHisto;

  AIDA::IHistogram1D * selxHisto;
  AIDA::IHistogram1D * selyHisto;
  AIDA::IHistogram1D * selaxHisto;
  AIDA::IHistogram1D * selayHisto;
  AIDA::IHistogram1D * seldxHisto;
  AIDA::IHistogram1D * seldyHisto;
  AIDA::IHistogram1D * selkxHisto;
  AIDA::IHistogram1D * selkyHisto;

  AIDA::IHistogram1D * seldx1Histo;
  AIDA::IHistogram1D * seldy1Histo;
  AIDA::IHistogram1D * seldx3Histo;
  AIDA::IHistogram1D * seldy3Histo;
  AIDA::IHistogram1D * seldx4Histo;
  AIDA::IHistogram1D * seldy4Histo;
  AIDA::IHistogram1D * seldx5Histo;
  AIDA::IHistogram1D * seldy5Histo;
  //AIDA::IHistogram1D * seldx6Histo;
  //AIDA::IHistogram1D * seldy6Histo;

  AIDA::IHistogram1D * gblndfHisto;
  AIDA::IHistogram1D * gblchi2aHisto;
  AIDA::IHistogram1D * gblchi2bHisto;
  AIDA::IHistogram1D * gblprbHisto;
  AIDA::IHistogram2D * gblprbxHisto;
  AIDA::IHistogram2D * gblprbyHisto;

  AIDA::IHistogram1D * badxHisto;
  AIDA::IHistogram1D * badyHisto;
  AIDA::IHistogram1D * badaxHisto;
  AIDA::IHistogram1D * badayHisto;
  AIDA::IHistogram1D * baddxHisto;
  AIDA::IHistogram1D * baddyHisto;
  AIDA::IHistogram1D * badkxHisto;
  AIDA::IHistogram1D * badkyHisto;

  AIDA::IHistogram1D * baddx1Histo;
  AIDA::IHistogram1D * baddy1Histo;
  AIDA::IHistogram1D * baddx3Histo;
  AIDA::IHistogram1D * baddy3Histo;
  AIDA::IHistogram1D * baddx4Histo;
  AIDA::IHistogram1D * baddy4Histo;
  AIDA::IHistogram1D * baddx5Histo;
  AIDA::IHistogram1D * baddy5Histo;
  //AIDA::IHistogram1D * baddx6Histo;
  //AIDA::IHistogram1D * baddy6Histo;

  AIDA::IHistogram1D * goodxHisto;
  AIDA::IHistogram1D * goodyHisto;
  AIDA::IHistogram1D * goodx1Histo;
  AIDA::IHistogram1D * goody1Histo;
  //AIDA::IHistogram1D * goodx6Histo;
  //AIDA::IHistogram1D * goody6Histo;

  AIDA::IHistogram1D * gblax0Histo;
  AIDA::IHistogram1D * gbldx0Histo;
  AIDA::IHistogram1D * gbldx01Histo;
  AIDA::IHistogram1D * gblrx0Histo;
  AIDA::IHistogram1D * gblry0Histo;
  AIDA::IHistogram1D * gblpx0Histo;
  AIDA::IHistogram1D * gblpx0_unbHisto;
  AIDA::IHistogram1D * gblpy0Histo;
  AIDA::IHistogram1D * gblpy0_unbHisto;
  AIDA::IHistogram1D * gblqx0Histo;
  AIDA::IProfile1D * gblrxvsx0;
  AIDA::IProfile1D * gblryvsy0;
  AIDA::IProfile1D * gblrxvsx01;
  AIDA::IProfile1D * gblryvsy01;
  AIDA::IProfile1D * gblrxvsxpix0;
  AIDA::IProfile1D * gblryvsypix0;
  AIDA::IProfile1D * gblrxvsxpix01;
  AIDA::IProfile1D * gblryvsypix01;
  AIDA::IProfile1D * gblrxvsxpix01_CS1;
  AIDA::IProfile1D * gblryvsypix01_CS1;
  AIDA::IProfile1D * gblrxvsxpix01_CS2;
  AIDA::IProfile1D * gblryvsypix01_CS2;
  AIDA::IProfile1D * gblrxvsxpix01_CS3;
  AIDA::IProfile1D * gblryvsypix01_CS3;
  AIDA::IProfile1D * gblrxvsxpix01_CS4;
  AIDA::IProfile1D * gblryvsypix01_CS4;
  AIDA::IProfile1D * gblrxvsxpix51_CS1;
  AIDA::IProfile1D * gblryvsypix51_CS1;
  AIDA::IProfile1D * gblrxvsxpix51_CS2;
  AIDA::IProfile1D * gblryvsypix51_CS2;
  AIDA::IProfile1D * gblrxvsxpix51_CS3;
  AIDA::IProfile1D * gblryvsypix51_CS3;
  AIDA::IProfile1D * gblrxvsxpix51_CS4;
  AIDA::IProfile1D * gblryvsypix51_CS4;
  AIDA::IProfile1D * gblrxvsxpix31;
  AIDA::IProfile1D * gblryvsypix31;
  AIDA::IProfile1D * gblrxvsxpix3cs1;
  AIDA::IProfile1D * gblryvsypix3cs1;
  AIDA::IProfile1D * gblrxvsxpix3cs2;
  AIDA::IProfile1D * gblryvsypix3cs2;
  AIDA::IProfile1D * gblrxvsxpix3cs3;
  AIDA::IProfile1D * gblryvsypix3cs3;
  AIDA::IProfile1D * gblrxvsxpix3cs4;
  AIDA::IProfile1D * gblryvsypix3cs4;
  AIDA::IProfile1D * gblrxvsxpix3cs5;
  AIDA::IProfile1D * gblryvsypix3cs5;
  AIDA::IProfile1D * gblrxvsxpix3cs6;
  AIDA::IProfile1D * gblryvsypix3cs6;


  AIDA::IHistogram1D * gblax1Histo;
  AIDA::IHistogram1D * gbldx1Histo;
  AIDA::IHistogram1D * gbldx11Histo;
  AIDA::IHistogram1D * gblrx1Histo;
  AIDA::IHistogram1D * gblry1Histo;
  AIDA::IHistogram1D * gblpx1Histo;
  AIDA::IHistogram1D * gblpx1_unbHisto;
  AIDA::IHistogram1D * gblpy1Histo;
  AIDA::IHistogram1D * gblpy1_unbHisto;
  AIDA::IHistogram1D * gblqx1Histo;
  AIDA::IHistogram1D * gblsx1Histo;
  AIDA::IHistogram1D * gbltx1Histo;

  AIDA::IHistogram1D * gblax2Histo;
  AIDA::IHistogram1D * gbldx2Histo;
  AIDA::IHistogram1D * gbldx21Histo;
  AIDA::IHistogram1D * gblrx2Histo;
  AIDA::IHistogram1D * gblry2Histo;
  AIDA::IHistogram1D * gblpx2Histo;
  AIDA::IHistogram1D * gblpx2_unbHisto;
  AIDA::IHistogram1D * gblpy2Histo;
  AIDA::IHistogram1D * gblpy2_unbHisto;
  AIDA::IHistogram1D * gblqx2Histo;
  AIDA::IHistogram1D * gblsx2Histo;
  AIDA::IHistogram1D * gbltx2Histo;

  AIDA::IHistogram1D * gblax3Histo;
  AIDA::IHistogram1D * gbldx3Histo;
  AIDA::IHistogram1D * gbldx31Histo;
  AIDA::IHistogram1D * gblrx3Histo;
  AIDA::IHistogram1D * gblry3Histo;
  AIDA::IHistogram1D * gblrx3_cs1Histo;
  AIDA::IHistogram1D * gblry3_cs1Histo;
  AIDA::IHistogram1D * gblrx3_cs2Histo;
  AIDA::IHistogram1D * gblry3_cs2Histo;
  AIDA::IHistogram1D * gblrx3_cs3Histo;
  AIDA::IHistogram1D * gblry3_cs3Histo;
  AIDA::IHistogram1D * gblrx3_cs4Histo;
  AIDA::IHistogram1D * gblry3_cs4Histo;
  AIDA::IHistogram1D * gblrx3_cs5Histo;
  AIDA::IHistogram1D * gblry3_cs5Histo;
  AIDA::IHistogram1D * gblrx3_cs6Histo;
  AIDA::IHistogram1D * gblry3_cs6Histo;
  AIDA::IHistogram1D * gblpx3Histo;
  AIDA::IHistogram1D * gblpy3Histo;
  AIDA::IHistogram1D * gblpx3_cs1Histo;
  AIDA::IHistogram1D * gblpy3_cs1Histo;
  AIDA::IHistogram1D * gblpx3_cs2Histo;
  AIDA::IHistogram1D * gblpy3_cs2Histo;
  AIDA::IHistogram1D * gblpx3_cs3Histo;
  AIDA::IHistogram1D * gblpy3_cs3Histo;
  AIDA::IHistogram1D * gblpx3_cs4Histo;
  AIDA::IHistogram1D * gblpy3_cs4Histo;
  AIDA::IHistogram1D * gblpx3_cs5Histo;
  AIDA::IHistogram1D * gblpy3_cs5Histo;
  AIDA::IHistogram1D * gblpx3_cs6Histo;
  AIDA::IHistogram1D * gblpy3_cs6Histo;
  AIDA::IHistogram1D * gblpx3_cs7Histo;
  AIDA::IHistogram1D * gblpy3_cs7Histo;
  AIDA::IHistogram1D * gblpx3_unbHisto;
  AIDA::IHistogram1D * gblpy3_unbHisto;
  AIDA::IHistogram1D * gblqx3Histo;
  AIDA::IHistogram1D * gblsx3Histo;
  AIDA::IHistogram1D * gbltx3Histo;

  AIDA::IHistogram1D * gblax4Histo;
  AIDA::IHistogram1D * gbldx4Histo;
  AIDA::IHistogram1D * gbldx41Histo;
  AIDA::IHistogram1D * gblrx4Histo;
  AIDA::IHistogram1D * gblry4Histo;
  AIDA::IHistogram1D * gblpx4Histo;
  AIDA::IHistogram1D * gblpx4_unbHisto;
  AIDA::IHistogram1D * gblpy4Histo;
  AIDA::IHistogram1D * gblpy4_unbHisto;
  AIDA::IHistogram1D * gblqx4Histo;
  AIDA::IHistogram1D * gblsx4Histo;
  AIDA::IHistogram1D * gbltx4Histo;

  AIDA::IHistogram1D * gblax5Histo;
  AIDA::IHistogram1D * gbldx5Histo;
  AIDA::IHistogram1D * gbldx51Histo;
  AIDA::IHistogram1D * gblrx5Histo;
  AIDA::IHistogram1D * gblry5Histo;
  AIDA::IHistogram1D * gblpx5Histo;
  AIDA::IHistogram1D * gblpx5_unbHisto;
  AIDA::IHistogram1D * gblpy5Histo;
  AIDA::IHistogram1D * gblpy5_unbHisto;
  AIDA::IHistogram1D * gblqx5Histo;

  /*AIDA::IHistogram1D * gblax6Histo;
  AIDA::IHistogram1D * gbldx6Histo;
  AIDA::IHistogram1D * gbldy6Histo;
  AIDA::IHistogram1D * gblrx6Histo;
  AIDA::IHistogram1D * gblry6Histo;
  AIDA::IHistogram1D * gblpx6Histo;
  AIDA::IHistogram1D * gblpy6Histo;
  AIDA::IHistogram1D * gblqx6Histo;
  AIDA::IHistogram1D * gblsx6Histo;
  AIDA::IHistogram1D * gbltx6Histo;*/

  AIDA::IHistogram1D * gblkxCentreHisto;
  AIDA::IHistogram1D * gblkxCentre1Histo;
  AIDA::IHistogram1D * gblkx1Histo;
  AIDA::IHistogram1D * gblkx2Histo;
  AIDA::IHistogram1D * gblkx3Histo;
  AIDA::IHistogram1D * gblkx4Histo;
  AIDA::IHistogram1D * gblkx5Histo;
  //AIDA::IHistogram1D * gblkx6Histo;

  AIDA::IHistogram1D * sixzx3Histo;
  AIDA::IHistogram1D * sixzy3Histo;
  AIDA::IHistogram1D * sixzx2Histo;
  AIDA::IHistogram1D * sixzy2Histo;
  AIDA::IHistogram1D * sixzx1Histo;
  AIDA::IHistogram1D * sixzy1Histo;

  AIDA::IHistogram1D * sixkxzyHisto;
  AIDA::IHistogram1D * sixkyzxHisto;
  AIDA::IHistogram1D * sixkxzxHisto;
  AIDA::IHistogram1D * sixkyzyHisto;

  AIDA::IHistogram1D * hIso;

  AIDA::IProfile2D * gblnxy_plane0;
  AIDA::IHistogram2D * gblnxy1_plane0;
  AIDA::IProfile2D * gblnxy_plane3;
  AIDA::IHistogram2D * gblnxy1_plane3;

  AIDA::IHistogram2D * gblnCS1xy_plane0;
  AIDA::IHistogram2D * gblnCS2xy_plane0;
  AIDA::IHistogram2D * gblnCS3xy_plane0;
  AIDA::IHistogram2D * gblnCS4xy_plane0;
  AIDA::IHistogram2D * gblnCS5xy_plane0;
  AIDA::IHistogram2D * gblnCS6xy_plane0;
  AIDA::IHistogram2D * gblnCS7xy_plane0;

  AIDA::IHistogram2D * gblnCS1xy_plane3;
  AIDA::IHistogram2D * gblnCS2xy_plane3;
  AIDA::IHistogram2D * gblnCS3xy_plane3;
  AIDA::IHistogram2D * gblnCS4xy_plane3;
  AIDA::IHistogram2D * gblnCS5xy_plane3;
  AIDA::IHistogram2D * gblnCS6xy_plane3;
  AIDA::IHistogram2D * gblnCS7xy_plane3;

  AIDA::IHistogram2D * gblnCS1xy1_plane3;
  AIDA::IHistogram2D * gblnCS2xy1_plane3;
  AIDA::IHistogram2D * gblnCS3xy1_plane3;
  AIDA::IHistogram2D * gblnCS4xy1_plane3;
  AIDA::IHistogram2D * gblnCS5xy1_plane3;
  AIDA::IHistogram2D * gblnCS6xy1_plane3;
  AIDA::IHistogram2D * gblnCS7xy1_plane3;

#endif

  };

  //! A global instance of the processor:
  //EUTelGBLFitter aEUTelGBLFitter;


}

#endif
