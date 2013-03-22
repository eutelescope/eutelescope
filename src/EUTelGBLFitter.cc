/* 
 * File:   EUTelGBLFitter.cc
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Created on January 25, 2013, 2:53 PM
 */

#ifdef USE_GBL
 
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
        double x0 = trackCand.front()->getPosition()[0];
        double z0 = trackCand.front()->getPosition()[2];

        double xLast = trackCand.back()->getPosition()[0];
        double zLast = trackCand.back()->getPosition()[2];

        double x = x0 - (x0 - xLast) * ((z0 - z) / (z0 - zLast));

        return x;
    }

    //! Predict track hit in Y direction using simplified model

    double EUTelGBLFitter::InterpolateTrackY(const EVENT::TrackerHitVec& trackCand, const double z) const {
        double y0 = trackCand.front()->getPosition()[1];
        double z0 = trackCand.front()->getPosition()[2];

        double yLast = trackCand.back()->getPosition()[1];
        double zLast = trackCand.back()->getPosition()[2];

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

        TVectorD meas(2);

        const double resx = 0.017/sqrt(12.); // [um] telescope initial resolution
        const double resy = 0.017/sqrt(12.); // [um] telescope initial resolution
        //const double resx = 18/sqrt(12.); // [um] telescope initial resolution
        //const double resy = 18/sqrt(12.); // [um] telescope initial resolution

        TVectorD measPrec(2); // precision = 1/resolution^2
        measPrec[0] = 1.0 / resx / resx;
        measPrec[1] = 1.0 / resy / resy;

	TVectorD scat(2);
	scat.Zero();
	
        double _eBeam = 3.; // Run 5063
//        double _eBeam = 5.; // Run 5080
        double p = _eBeam; // beam momentum
        double X0Si = 50e-3 / 94; // Si 
        double X0Kap = 60e-3 / 286; // Kapton
        double X0Air = 150 / 304e3; // Air
        double tetSi = (0.0136 * sqrt(X0Si) / p * (1 + 0.038 * std::log(X0Si)));
        double tetKap = (0.0136 * sqrt(X0Kap) / p * (1 + 0.038 * std::log(X0Kap)));
        double tetAir = (0.0136 * sqrt(X0Air) / p * ( 1 + 0.038*std::log(X0Air)));

        TVectorD scatPrecSensor(2);
        scatPrecSensor[0] = 1.0 / (tetSi * tetSi + tetKap * tetKap);
        scatPrecSensor[1] = 1.0 / (tetSi * tetSi + tetKap * tetKap);

        TVectorD scatPrecAir(2);
        scatPrecAir[0] = 1.0 / (tetAir * tetAir);
        scatPrecAir[1] = 1.0 / (tetAir * tetAir);

        TMatrixD alDer( 2, 3 ); // alignment derivatives
        alDer.Zero();
	alDer[0][0] = 1.0; // dx/dx
        alDer[0][1] = 0.0; // dx/dy
	alDer[1][0] = 0.0; // dy/dx
	alDer[1][1] = 1.0; // dy/dy
        
        std::map< int, EVENT::TrackerHitVec >::const_iterator itTrkCand = _trackCandidates.begin();
        EVENT::TrackerHitVec::const_iterator itHit;
        for (; itTrkCand != _trackCandidates.end(); ++itTrkCand) {

            if( itTrkCand->second.size() > 6 ) continue;
            
            IMPL::TrackImpl * fittrack = new IMPL::TrackImpl();

            //GBL trajectory construction
            std::vector< gbl::GblPoint > pointList;

            int iLabel = 0;
            int iPlane = 0;
            TMatrixD jacPointToPoint(5, 5);
            jacPointToPoint.UnitMatrix();
            double step = 0.;
            double zprev = itTrkCand->second.front()->getPosition()[2];
            itHit = itTrkCand->second.begin();
            for (; itHit != itTrkCand->second.end(); ++itHit) {
                const double* hitpos = (*itHit)->getPosition();
//                step = hitpos[2] - zprev;			// commented out beacause propagation is done in air's sub-block
//                streamlog_out(DEBUG0) << "Step size:" << step << std::endl;	// commented out beacause propagation is done in air's sub-block
//                jacPointToPoint = PropagatePar(step);		// commented out beacause propagation is done in air's sub-block
                //jacPointToPoint.Print();
//                
                std::vector<int> globalLabels(3);
//                std::vector<int> globalLabels(2);	// fit w/o rotations
                globalLabels[0] = 1 + iPlane*10; // dx
                globalLabels[1] = 2 + iPlane*10; // dy
                globalLabels[2] = 3 + iPlane*10; // rot
//                
//                
                double xPred = InterpolateTrackX(itTrkCand->second, hitpos[2]);
                double yPred = InterpolateTrackY(itTrkCand->second, hitpos[2]);
//                
                alDer[0][2] = -yPred; // dx/rot
                alDer[1][2] = xPred; // dy/rot
                
                fittrack -> addHit(*itHit);

                gbl::GblPoint point(jacPointToPoint);
                meas[0] = hitpos[0] - xPred; // [um]
                meas[1] = hitpos[1] - yPred;
                streamlog_out(DEBUG0) << "Residuals:" << std::endl;
                streamlog_out(DEBUG0) << "X:" << std::setw(20) << meas[0] << std::setw(20) << resx << std::endl;
                streamlog_out(DEBUG0) << "Y:" << std::setw(20) << meas[1] << std::setw(20) << resy << std::endl;
                TMatrixD proL2m(2, 2);
                proL2m.UnitMatrix();
                ++iLabel;
                point.setLabel(iLabel);
                point.addMeasurement(proL2m, meas, measPrec);
                point.addScatterer(scat, scatPrecSensor);
                point.addGlobals(globalLabels,alDer);
                pointList.push_back(point);

		// construct effective scatterers for air
		// the scatters must be at (Z(plane i) + Z(plane i+1))/2. +/- (Z(plane i) - Z(plane i+1))/sqrt(12)
		if ( iPlane < itTrkCand->second.size()-1 ) {
			const double planeSpacing = hitpos[2] - (*(itHit+1))->getPosition()[2];		// works only for 6 hits tracks
			// propagate parameters into air gap
			step = planeSpacing/2. - planeSpacing/sqrt(12.);
			//step = planeSpacing/2.;
			streamlog_out(DEBUG0) << "Step size in the air gap:" << step << std::endl;
			jacPointToPoint = PropagatePar(step);
			// point with scatterer
			gbl::GblPoint pointInAir1(jacPointToPoint);
			//pointInAir1.setLabel(1000+iLabel);
			pointInAir1.setLabel(++iLabel);
			pointInAir1.addScatterer(scat, scatPrecAir);
			pointList.push_back(pointInAir1);

			step = 2.*planeSpacing/sqrt(12.);
			streamlog_out(DEBUG0) << "Step size in the air gap:" << step << std::endl;
			jacPointToPoint = PropagatePar(step);
			// point with scatterer
			gbl::GblPoint pointInAir2(jacPointToPoint);
			//pointInAir2.setLabel(2000+iLabel);
			pointInAir2.setLabel(++iLabel);
			pointInAir2.addScatterer(scat, scatPrecAir);
			pointList.push_back(pointInAir2);

			step = planeSpacing/2. - planeSpacing/sqrt(12.); 
			//step = planeSpacing/2.;
			streamlog_out(DEBUG0) << "Step size in the air gap:" << step << std::endl;
			jacPointToPoint = PropagatePar(step);
		}

                ++iPlane;
                zprev = hitpos[2];
            }

            gbl::GblTrajectory* traj = new gbl::GblTrajectory(pointList, false);
            _gblTrackPoints.push_back(pointList);
            double chi2 = 0.;
            double loss = 0.;
            int ndf = 0;
            int ierr = 0;
            const std::string mEstOpt = "T";
            ierr = traj->fit(chi2, ndf, loss, mEstOpt);
//            ierr = traj->fit(chi2, ndf, loss);

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
//            delete fittrack;

            streamlog_out(DEBUG1) << "Fit results:" << std::endl;
            streamlog_out(DEBUG1) << "Chi2:" << std::setw(10) << chi2 << std::endl;
            streamlog_out(DEBUG1) << "NDF:" << std::setw(10) << ndf << std::endl;
            streamlog_out(DEBUG1) << "M-estimator loss:" << std::setw(10) << loss << std::endl;
        }


        return;
    }
}

#endif
