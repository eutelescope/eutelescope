/* 
 * File:   EUTelGBLFitter.cc
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Created on January 25, 2013, 2:53 PM
 */

// eutelescope includes ".h"
#include "EUTelGBLFitter.h"
#include "EUTelTrackFitter.h"
#include "EUTELESCOPE.h"

// marlin util includes
#include "mille/Mille.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// GBL
#include "include/GblTrajectory.h"
#include "include/GblPoint.h"
#include "include/GblData.h"
#include "include/BorderedBandMatrix.h"
#include "include/MilleBinary.h"
#include "include/VMatrix.h"

// LCIO includes <.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include "lcio.h"

// system includes <>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>

namespace eutelescope {

    EUTelGBLFitter::EUTelGBLFitter() : EUTelTrackFitter("GBLTrackFitter"),
    _alignmentMode(kXYShift),
    _parPropJac(5, 5) {
        _trackCandidates.clear();
        _gblTrackCandidates.clear();

        //
        _fittrackvec = new IMPL::LCCollectionVec(EVENT::LCIO::TRACK);
        //        LCCollectionVec *fitpointvec = new LCCollectionVec(LCIO::TRACKERHIT);

    }

    EUTelGBLFitter::EUTelGBLFitter(std::string name) : EUTelTrackFitter(name),
    _alignmentMode(kXYShift),
    _parPropJac(5, 5) {
        _trackCandidates.clear();
        _gblTrackCandidates.clear();

        _fittrackvec = new IMPL::LCCollectionVec(EVENT::LCIO::TRACK);
    }

    EUTelGBLFitter::~EUTelGBLFitter() {
        delete _fittrackvec;
    }

    TMatrixD EUTelGBLFitter::PropagatePar(double ds) {
        /* for GBL:
           Jacobian for straight line track
           track = q/p, x', y', x, y
                    0,   1,  2, 3, 4
         */
        _parPropJac.UnitMatrix();
        _parPropJac[3][1] = ds; // x = x0 + xp * ds
        _parPropJac[4][2] = ds; // y = y0 + yp * ds
        return _parPropJac;
    }

    void EUTelGBLFitter::SetTrackCandidates(std::map< int, EVENT::TrackerHitVec >& trackCandidates) {
        this->_trackCandidates = trackCandidates;
        return;
    }

    //! Calculates rough estimate of track's intersection point

    //    double* EUTelGBLFitter::GetTrackOffset(const Utility::HitsPVec& trackCand) const {
    //        double a[3] = {(*trackCand.begin())->X(), (*trackCand.begin())->Y(), (*trackCand.begin())->Z()};
    //
    //        return a;
    //    }
    //
    //    double* EUTelGBLFitter::GetTrackSlope(const Utility::HitsPVec& trackCand) const {
    //        TVector3 r1Vec((*trackCand.begin())->X(), (*trackCand.begin())->Y(), (*trackCand.begin())->Z());
    //        TVector3 r2Vec((*trackCand.end())->X(), (*trackCand.end())->Y(), (*trackCand.end())->Z());
    //        double r12Dist = TVector3(r1Vec - r2Vec).Mag();
    //
    //        TVector3 kVec = r2Vec - r1Vec;
    //        kVec *= 1 / r12Dist;
    //
    //        double k[3] = {0., 0., 0.};
    //        k[0] = kVec.X();
    //        k[1] = kVec.Y();
    //        k[2] = kVec.Z();
    //
    //        return k;
    //    }
    //
    //    //! Predict track hit in X direction using simplified model
    //

    double EUTelGBLFitter::InterpolateTrackX(const EVENT::TrackerHitVec& trackCand, const double z) const {
        double x0 = (*trackCand.begin())->getPosition()[0];
        double z0 = (*trackCand.begin())->getPosition()[0];

        double xLast = trackCand[trackCand.size() - 1]->getPosition()[0];
        double zLast = trackCand[trackCand.size() - 1]->getPosition()[1];

        double x = x0 - (x0 - xLast) * ((z0 - z) / (z0 - zLast));

        return x;
    }

    //! Predict track hit in Y direction using simplified model

    double EUTelGBLFitter::InterpolateTrackY(const EVENT::TrackerHitVec& trackCand, const double z) const {
        double y0 = (*trackCand.begin())->getPosition()[0];
        double z0 = (*trackCand.begin())->getPosition()[2];

        double yLast = trackCand[trackCand.size() - 1]->getPosition()[0];
        double zLast = trackCand[trackCand.size() - 1]->getPosition()[2];

        double y = y0 - (y0 - yLast) * ((z0 - z) / (z0 - zLast));

        return y;
    }

    void EUTelGBLFitter::Reset() {
        if (!_gblTrackCandidates.empty()) {
            std::map< int, gbl::GblTrajectory* >::iterator it;
            for (it = _gblTrackCandidates.begin(); it != _gblTrackCandidates.begin(); ++it) delete it->second;
        }
        _gblTrackCandidates.clear();
        _fittrackvec->clear();
        _gblTrackPoints.clear();
    }

    void EUTelGBLFitter::FitTracks() {
        Reset(); // 

        streamlog_out(DEBUG1) << " EUTelGBLFitter::FitTracks() " << std::endl;
        streamlog_out(DEBUG1) << " N track candidates:" << (int) _trackCandidates.size() << std::endl;

        streamlog_out(DEBUG1) << "Following track candidates supplied:" << std::endl;
        std::map< int, EVENT::TrackerHitVec >::const_iterator itTrkCandHelp = _trackCandidates.begin();
        EVENT::TrackerHitVec::const_iterator itHitHelp;
        for (; itTrkCandHelp != _trackCandidates.end(); ++itTrkCandHelp) {
            streamlog_out(DEBUG1) << std::setfill('-') << std::setw(10) << itTrkCandHelp->first << std::setw(10) << ":" << std::setfill(' ') << std::endl;
            itHitHelp = itTrkCandHelp->second.begin();
            for (; itHitHelp != itTrkCandHelp->second.end(); ++itHitHelp) {
                streamlog_out(DEBUG0) << *itHitHelp << std::endl;
                streamlog_out(DEBUG0) << std::setw(10) << (*itHitHelp)->getPosition()[0]
                        << std::setw(10) << (*itHitHelp)->getPosition()[1]
                        << std::setw(10) << (*itHitHelp)->getPosition()[2] << std::endl;
            }
        }

        TMatrixD jacPointToPoint(5, 5);
        TMatrixD proL2m(2, 2);
        TVectorD meas(2);

        const double resx = 21; // [um] telescope initial resolution
        const double resy = 20; // [um] telescope initial resolution

        TVectorD measPrec(2); // precision = 1/resolution^2
        measPrec[0] = 1.0 / resx / resx;
        measPrec[1] = 1.0 / resy / resy;

        TVectorD scat(2);
        scat[0] = 0.;
        scat[1] = 0.;

        double _eBeam = 3.; // Run 5063
        double p = _eBeam; // beam momentum
        double X0Si = 65e-3 / 94; // Si 
        double X0Kap = 60e-3 / 286; // Kapton
        double X0Air = 50 / 304e3; // Air
        double tetSi = (0.0136 * sqrt(X0Si) / p * (1 + 0.038 * std::log(X0Si)));
        double tetKap = (0.0136 * sqrt(X0Kap) / p * (1 + 0.038 * std::log(X0Kap)));
        double tetAir = 0.; //(0.0136 * sqrt(X0Air) / p * ( 1 + 0.038*std::log(X0Air) ));

        TVectorD scatPrec(2);
        scatPrec[0] = 1.0 / (tetSi * tetSi + tetKap * tetKap);
        scatPrec[1] = 1.0 / (tetSi * tetSi + tetKap * tetKap);

        TMatrixD alDer( 2, 3 ); // alignment derivatives
	alDer[0][0] = 1.0; // dx/dx
        alDer[0][1] = 0.0; // dx/dy
	alDer[1][0] = 0.0; // dy/dx
	alDer[1][1] = 1.0; // dy/dy
        
        std::map< int, EVENT::TrackerHitVec >::const_iterator itTrkCand = _trackCandidates.begin();
        EVENT::TrackerHitVec::const_iterator itHit;
        for (; itTrkCand != _trackCandidates.end(); ++itTrkCand) {

            IMPL::TrackImpl * fittrack = new IMPL::TrackImpl();

            //GBL trajectory construction
            std::vector< gbl::GblPoint > pointList;
            jacPointToPoint.UnitMatrix();
            proL2m.UnitMatrix();

            int iLabel = 0;
            itHit = itTrkCand->second.begin();
            double step = 0.;
            double zprev = (*itHit)->getPosition()[2];
            int iPlane = 0;
            for (; itHit != itTrkCand->second.end(); ++itHit) {
                
                std::vector<int> globalLabels(3);
                globalLabels[0] = 10 + iPlane; // dx
                globalLabels[1] = 20 + iPlane; // dy
                globalLabels[2] = 40 + iPlane; // rot
                
                ++iPlane;
                
                double xPred = InterpolateTrackX(itTrkCand->second, (*itHit)->getPosition()[2]);
                double yPred = InterpolateTrackY(itTrkCand->second, (*itHit)->getPosition()[2]);
                
                alDer[0][2] = -yPred; // dx/rot
                alDer[1][2] = xPred; // dy/rot
                
                fittrack -> addHit(*itHit);

                const double* hitpos = (*itHit)->getPosition();
                step = hitpos[2] - zprev;
                gbl::GblPoint point(jacPointToPoint);
                jacPointToPoint = PropagatePar(step);
                zprev = hitpos[2];
                meas[0] = (*itHit)->getPosition()[0] - xPred; // [um]
                meas[1] = (*itHit)->getPosition()[1] - yPred;
                streamlog_out(DEBUG0) << "Residuals:" << std::endl;
                streamlog_out(DEBUG0) << "X:" << std::setw(10) << meas[0] << std::setw(10) << resx << std::endl;
                streamlog_out(DEBUG0) << "Y:" << std::setw(10) << meas[1] << std::setw(10) << resy << std::endl;

                point.setLabel(++iLabel);
                point.addMeasurement(proL2m, meas, measPrec);
                point.addScatterer(scat, scatPrec);
                point.addGlobals(globalLabels,alDer);
                pointList.push_back(point);
            }

            gbl::GblTrajectory* traj = new gbl::GblTrajectory(pointList, false);
            _gblTrackPoints.push_back(pointList);
            double chi2 = 0.;
            double loss = 0.;
            int ndf = 0;
            int ierr = 0;
            const std::string mEstOpt = "T";
            ierr = traj->fit(chi2, ndf, loss, mEstOpt);

            if( !ierr ) traj->milleOut( *_mille );
            
            _gblTrackCandidates.insert(std::make_pair(_fittrackvec->getNumberOfElements(), traj));

            //            // Write fit result out
            //
            //
            //            // Following parameters are not used for Telescope
            //            // and are set to zero (just in case)
            fittrack->setOmega(0.); // curvature of the track
            fittrack->setD0(0.); // impact paramter of the track in (r-phi)
            fittrack->setZ0(0.); // impact paramter of the track in (r-z)
            fittrack->setPhi(0.); // phi of the track at reference point
            fittrack->setTanLambda(0.); // dip angle of the track at reference point
            //
            //            // Used class members
            //
            fittrack->setChi2(chi2); // Chi2 of the fit (including penalties)
            fittrack->setNdf(ndf); // Number of planes fired (!)
            //
            //            // Store track reference point.
            //            
            float refpoint[3] = {0., 0., 0.};
            fittrack->setReferencePoint(refpoint);
            //
            //            
            //            // Store track
            //            
            _fittrackvec->addElement(fittrack);
            //
            //	    delete fittrack;

            streamlog_out(DEBUG1) << "Fit results:" << std::endl;
            streamlog_out(DEBUG1) << "Chi2:" << std::setw(10) << chi2 << std::endl;
            streamlog_out(DEBUG1) << "NDF:" << std::setw(10) << ndf << std::endl;
            streamlog_out(DEBUG1) << "M-estimator loss:" << std::setw(10) << loss << std::endl;
        }


        return;
    }
}
