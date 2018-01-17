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
    void TelescopeCorrelationPlots(std::vector<EUTelTripletGBLUtility::hit> const & telescopehits);

    void fillTrackhitHisto(EUTelTripletGBLUtility::hit const & hit, int ipl);
  protected:
    std::string _inputCollectionTelescope;
  
    std::map<size_t, bool> _excludedSensorMap;

    //Analysis parameters
    bool _isFirstEvent;
    int _printEventCounter;
    double _eBeam;
    int _nRun;
    int _nEvt;
    size_t _nPlanes;
    int _ngbl;
    EVENT::IntVec _excluded_planes;
    double _eff_radius;
    double _kappa;
    
    // Cuts for matching:
    double _track_match_cut;
    double _slope_cut;
    double _triplet_res_cut;
    double _probchi2_cut;

    // Partly outdated GEAR readings:
    std::vector<int> _sensorIDVec;
    std::vector<double>  _planePosition;
    std::vector<double> _planeRadLength;
    std::map<size_t, std::vector<float>> _planeResolutionX;
    std::map<size_t, std::vector<float>> _planeResolutionY;

	  std::vector<Eigen::Vector2d> _planeWscatSi;
	  std::vector<Eigen::Vector2d> _planeWscatAir;

    FloatVec _telResolution;
    FloatVec _dutResolutionX;
    FloatVec _dutResolutionY;

    EUTelTripletGBLUtility gblutil;

      // definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  AIDA::IHistogram1D * nAllTelHitHisto;
  AIDA::IHistogram1D * nAllDUTHitHisto;

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

  std::vector<AIDA::IHistogram1D*> clustersizeTotal;
  std::vector<AIDA::IHistogram1D*> clustersizeX;
  std::vector<AIDA::IHistogram1D*> clustersizeY;

  std::vector<AIDA::IHistogram1D*> sixXHistos;
  std::vector<AIDA::IHistogram1D*> sixYHistos;

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

  AIDA::IHistogram1D * goodxHisto;
  AIDA::IHistogram1D * goodyHisto;
  AIDA::IHistogram1D * goodx1Histo;
  AIDA::IHistogram1D * goody1Histo;
  
  std::vector<AIDA::IHistogram1D*> gblaxHistos;
  std::vector<AIDA::IHistogram1D*> gbldxHistos;
  std::vector<AIDA::IHistogram1D*> gblrxHistos;
  std::vector<AIDA::IHistogram1D*> gblryHistos;
  std::vector<AIDA::IHistogram1D*> gblpxHistos;
  std::vector<AIDA::IHistogram1D*> gblpyHistos;
  std::vector<AIDA::IHistogram1D*> gblqxHistos;

  std::vector<AIDA::IProfile1D*> gblrxvsx;
  std::vector<AIDA::IProfile1D*> gblryvsy;
  std::vector<AIDA::IProfile1D*> gblrxvsx1;
  std::vector<AIDA::IProfile1D*> gblryvsy1;

  std::vector<AIDA::IProfile1D*> gblrxvsxpix;
  std::vector<AIDA::IProfile1D*> gblryvsypix;
  std::vector<AIDA::IProfile1D*> gblrxvsxpix1;
  std::vector<AIDA::IProfile1D*> gblryvsypix1;

  std::vector<std::vector<AIDA::IProfile1D*>> gblrxvsxpix1CS;
  std::vector<std::vector<AIDA::IProfile1D*>> gblryvsypix1CS;

  std::vector<AIDA::IHistogram1D*> gblkxHistos;

  AIDA::IHistogram1D * gblkxCentreHisto;
  AIDA::IHistogram1D * gblkxCentre1Histo;

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

  std::vector<AIDA::IProfile2D*> gblnxy;
  std::vector<AIDA::IHistogram2D*> gblnxy1;
  std::vector<AIDA::IHistogram2D*> gblcluxvscluy;

  AIDA::IProfile2D * gblnxy_plane0;
  AIDA::IHistogram2D * gblnxy1_plane0;
  AIDA::IProfile2D * gblnxy_plane3;
  AIDA::IHistogram2D * gblnxy1_plane3;

  std::vector<std::vector<AIDA::IHistogram2D*>> gblnCSxy_tot;
  std::vector<std::vector<AIDA::IHistogram2D*>> gblnCSxy_x;
  std::vector<std::vector<AIDA::IHistogram2D*>> gblnCSxy_y;

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
  EUTelGBLFitter aEUTelGBLFitter;
}
#endif
