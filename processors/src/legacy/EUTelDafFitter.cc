/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

/*!
  * This is a track fitting processor for the Eutelescope package.
 *
 * It performs track finding and fitting on a supplied hit collection.
 *
 * The track finder works by propagating all hits to plane 0, currently assuming
 * straight line fits, then running a cluster finder. Hit clusters above some set value
 * are considered track candidates.
 *
 * This track candidate is then fitted using a implementation of a Deterministic
 * Annealing Filter (DAF), that in short is a Kalman Filter running iteratively over 
 * a set of weighted measurements, reweighting the measurements after each fit based 
 * on the residuals and a supplied chi2 cut off.
 *
 * This package uses the Eigen library for linear algebra. This package is very
 * quick when compiled properly, but very slow when compiled for debugging. Make sure 
 * to compile properly before running productions.
 *
 * Running 'cmake -i' inside the build folder, and then when it asks
 * Variable Name: CMAKE_CXX_FLAGS_RELEASE
 * Description: Flags used by the compiler during release builds (/MD /Ob1 /Oi
 * /Ot /Oy /Gs will produce slightly less optimized but smaller files).
 *
 * enter:
 * New Value (Enter to keep current value): -O3 -msse2 -ftree-vectorize -DNDEBUG
 *
 * When it asks
 * Variable Name: CMAKE_BUILD_TYPE
 * enter:
 * New Value (Enter to keep current value): Release
 *
 * If your cpu supports it, you could try -msse4 or -msse3 as well.
 */

// built only if GEAR is used
#if defined(USE_GEAR)

// eutelescope includes ".h"
#include "EUTelDafFitter.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelVirtualCluster.h"

// marlin includes ".h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#include <marlin/AIDAProcessor.h>
#endif

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/TrackerPulse.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IO/LCWriter.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>

// system includes <>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelDafFitter::EUTelDafFitter() : EUTelDafBase("EUTelDafFitter") {
  // Child specific params and description
  dafParams();
}

void EUTelDafFitter::dafParams() {

  _description = "EUTelDafFitter performs track reconstruction. The tracks are "
                 "as final track fit for analysis.";

  registerOutputCollection(LCIO::TRACK, 
			   "TrackCollectionName",
                           "Collection name for fitted tracks",
                           _trackCollectionName, 
			   std::string("fittracks"));

  registerOptionalParameter("AddToLCIO",
			    "Should plots be made and filled?",
                            _addToLCIO,
			    true);

  registerOptionalParameter("FitDuts",
			    "Set this to true if you want DUTs to be included in the track fit.",
			    _fitDuts,
			    false);  
}

void EUTelDafFitter::dafInit() {
  if (_fitDuts) {
    for (size_t ii = 0; ii < _system.planes.size(); ii++) {
      if (find(_dutPlanes.begin(), _dutPlanes.end(),
               _system.planes.at(ii).getSensorID()) != _dutPlanes.end()) {
        _system.planes.at(ii).include();
      }
    }
  }
}

void EUTelDafFitter::dafEvent(LCEvent *event) {

  // Prepare track collection
  if (_addToLCIO) {
    // Define output track and hit collections
    _fittrackVec = new LCCollectionVec(LCIO::TRACK);
    _fitpointVec = new LCCollectionVec(LCIO::TRACKERHIT);
    // Set flag for storing track hits in track collection
    LCFlagImpl flag(_fittrackVec->getFlag());
    flag.setBit(LCIO::TRBIT_HITS);
    _fittrackVec->setFlag(flag.getFlag());
  }

  // Check found tracks
  for (size_t ii = 0; ii < _system.getNtracks(); ii++) {
    // run track fitte
    _nCandidates++;
    // prepare track for DAF fit
    _system.fitPlanesInfoDaf(_system.tracks.at(ii));
    // check resids, intime, angles
    if (not checkTrack(_system.tracks.at(ii))) {
      continue;
    };
    int inTimeHits = checkInTime(_system.tracks.at(ii));
    if (inTimeHits < _nDutHits) {
      continue;
    }

    // Fill plots
    if (_histogramSwitch) {
      fillPlots(_system.tracks.at(ii));
      fillDetailPlots(_system.tracks.at(ii));
    }
    // Dump to LCIO
    if (_addToLCIO) {
      addToLCIO(_system.tracks.at(ii), _fitpointVec);
    }
    _nTracks++;
  }

  // Add track collection
  if (_addToLCIO) {
    event->addCollection(_fittrackVec, _trackCollectionName);
    std::string sfitpoints = "";

    for (int i = 0; i < 2000; i++) {
      sfitpoints = "fitpoints" + i;
      try {
        dynamic_cast<LCCollectionVec *>(event->getCollection(sfitpoints));
      } catch (...) {
        break;
      }
    }
    event->addCollection(_fitpointVec, sfitpoints);
  }
}

void EUTelDafFitter::addToLCIO(daffitter::TrackCandidate<float, 4> &track,
                               LCCollectionVec *lcvec) {
  TrackImpl *fittrack = new TrackImpl();
  // Impact parameters are useless and set to 0
  fittrack->setD0(0.);        // impact parameter of the track in (r-phi)
  fittrack->setZ0(0.);        // impact parameter of the track in (r-z)
  fittrack->setTanLambda(0.); // dip angle of the track at reference point

  daffitter::TrackEstimate<float, 4> &est = track.estimates.at(0);

  // No good way of storing the track angles, so
  fittrack->setOmega(est.getXdz()); // Storing dxdz as Omega
  fittrack->setPhi(est.getYdz());   // Storing dx/dy as phi

  fittrack->setChi2(track.chi2);
  fittrack->setNdf(int(round(track.ndof)));
  // prepare an encoder for the hit collection to store properties
  CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING, lcvec);

  float refpoint[3];

  for (size_t plane = 0; plane < _system.planes.size(); plane++) {
    daffitter::FitPlane<float> &pl = _system.planes.at(plane);
    daffitter::TrackEstimate<float, 4> &estim = track.estimates.at(plane);
    TrackerHitImpl *fitpoint = new TrackerHitImpl();
    // encode and store sensorID
    int sensorID = _system.planes.at(plane).getSensorID();
    idHitEncoder["sensorID"] = sensorID;
    // set local/global bit flag property AND FittedHit property for the hit
    idHitEncoder["properties"] = kHitInGlobalCoord | kFittedHit;
    double pos[3];
    pos[0] = estim.getX() / 1000.0;
    pos[1] = estim.getY() / 1000.0;
    pos[2] = pl.getMeasZ() / 1000.0;

    fitpoint->setPosition(pos);
    // Covariance matrix of the fitted position
    // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).
    float cov[TRKHITNCOVMATRIX];
    cov[0] = estim.cov(0, 0);
    cov[1] = estim.cov(0, 1);
    cov[2] = estim.cov(1, 1);
    // Error of z position of fit is unclear to me, this would be a systematic
    // alignment error. Set to 0 along with all covariances.
    cov[3] = cov[4] = cov[5] = 0.;
    fitpoint->setCovMatrix(cov);
    // store values
    idHitEncoder.setCellID(fitpoint);
    _fitpointVec->push_back(fitpoint);
    fittrack->addHit(fitpoint);

    streamlog_out(DEBUG3) << " hit : sensorID " << idHitEncoder["sensorID"]
                          << " properties: " << idHitEncoder["properties"]
                          << std::endl;

    if (plane == 0) {
      refpoint[0] = pos[0];
      refpoint[1] = pos[1];
      refpoint[2] = pos[2];
    }
    // Also add measurement point
    for (size_t mm = 0; mm < pl.meas.size(); mm++) {
      if (track.weights.at(plane)(mm) < 0.5f) {
        continue;
      }
      if (pl.isExcluded()) {
        continue;
      }
      TrackerHitImpl *meashit = static_cast<TrackerHitImpl *>(
          _hitCollection->getElementAt(pl.meas.at(mm).getIden()));
      fittrack->addHit(meashit);
    }
  }
  fittrack->setReferencePoint(refpoint);
  _fittrackVec->addElement(fittrack);
}

void EUTelDafFitter::dafEnd() {
}
#endif // USE_GEAR
