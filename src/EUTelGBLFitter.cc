/* 
 * File:   EUTelGBLFitter.cc
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Created on January 25, 2013, 2:53 PM
 */

#ifdef USE_GBL

// eutelescope includes ".h"
#include "EUTelGeometryTelescopeGeoDescription.h"
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
#include "LCIOTypes.h"
#include "lcio.h"

// system includes <>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <iterator>
#include <algorithm>

namespace eutelescope {

    EUTelGBLFitter::EUTelGBLFitter() : EUTelTrackFitter("GBLTrackFitter"),
    _trackCandidates(),
    _gblTrackCandidates(),
    _fittrackvec(0),
    _parPropJac(5, 5),
    _eBeam(-1.),
    _hitId2GblPointLabel(),
    _alignmentMode(Utility::XYShiftXYRot),
    _mEstimatorType(),
    _mille(0),
    _paramterIdXShiftsMap(),
    _paramterIdYShiftsMap(),
    _paramterIdZShiftsMap(),
    _paramterIdXRotationsMap(),
    _paramterIdYRotationsMap(),
    _paramterIdZRotationsMap(),
    _chi2cut(1000.) {}

    EUTelGBLFitter::EUTelGBLFitter(std::string name) : EUTelTrackFitter(name),
    _trackCandidates(),
    _gblTrackCandidates(),
    _fittrackvec(0),
    _parPropJac(5, 5),
    _eBeam(-1.),
    _hitId2GblPointLabel(),
    _alignmentMode(Utility::XYShiftXYRot),
    _mEstimatorType(),
    _mille(0),
    _paramterIdXShiftsMap(),
    _paramterIdYShiftsMap(),
    _paramterIdZShiftsMap(),
    _paramterIdXRotationsMap(),
    _paramterIdYRotationsMap(),
    _paramterIdZRotationsMap(),
    _chi2cut(1000.) {}

    EUTelGBLFitter::~EUTelGBLFitter() {
    }
    
    
    void EUTelGBLFitter::setParamterIdXRotationsMap( const std::map<int, int>& map ) {
        _paramterIdXRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdYRotationsMap( const std::map<int, int>& map ) {
        _paramterIdYRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdZRotationsMap( const std::map<int, int>& map ) {
        _paramterIdZRotationsMap = map;
    }

    void EUTelGBLFitter::setParamterIdZShiftsMap( const std::map<int, int>& map ) {
        _paramterIdZShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdYShiftsMap( const std::map<int, int>& map ) {
        _paramterIdYShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdXShiftsMap( const std::map<int, int>& map ) {
        _paramterIdXShiftsMap = map;
    }
    
    std::map<int, int> EUTelGBLFitter::getParamterIdXRotationsMap() const {
        return _paramterIdXRotationsMap;
    }
    
    std::map<int, int> EUTelGBLFitter::getParamterIdYRotationsMap() const {
        return _paramterIdYRotationsMap;
    }
    
    std::map<int, int> EUTelGBLFitter::getParamterIdZRotationsMap() const {
        return _paramterIdZRotationsMap;
    }

    std::map<int, int> EUTelGBLFitter::getParamterIdZShiftsMap() const {
        return _paramterIdZShiftsMap;
    }

    std::map<int, int> EUTelGBLFitter::getParamterIdYShiftsMap() const {
        return _paramterIdYShiftsMap;
    }

    std::map<int, int> EUTelGBLFitter::getParamterIdXShiftsMap() const {
        return _paramterIdXShiftsMap;
    }

    void EUTelGBLFitter::setMEstimatorType( const std::string& mEstimatorType ) {
        std::string mEstimatorTypeLowerCase = mEstimatorType;
        std::transform( mEstimatorType.begin(), mEstimatorType.end(), mEstimatorTypeLowerCase.begin(), ::tolower);
        
        if ( mEstimatorType.size() != 1 ) {
            streamlog_out( WARNING1 ) << "More than one character supplied as M-estimator option" << std::endl;
            streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
            return;
        }
        
        if ( mEstimatorType.compare("t") == 0 ||
             mEstimatorType.compare("h") == 0 ||
             mEstimatorType.compare("c") == 0   ) this->_mEstimatorType = _mEstimatorType;
        else {
            streamlog_out( WARNING1 ) << "M-estimator option " << mEstimatorType << " was not recognized" << std::endl;
            streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
        }
    }

    std::string EUTelGBLFitter::getMEstimatorType( ) const {
        return _mEstimatorType;
    }

    std::map<int, int> EUTelGBLFitter::getHitId2GblPointLabel( ) const {
        return _hitId2GblPointLabel;
    }

    TMatrixD EUTelGBLFitter::propagatePar(double ds) {
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
    
    /** Propagation jacobian for solenoidal magnetic field
     * 
     *  This is C++ port of C. Kleinwort's python code b2thlx
     * 
     * @param ds        (3D) arc length to endpoint
     * @param phi       azimuthal direction at starting point
     * @param bfac      Magnetic field strength
     * @return 
     */
    /*
    TMatrixD EUTelGBLFitter::PropagatePar( double ds, double phi, double bfac ) {
        // for GBL:
        //   Jacobian for helical track
        //   track = q/p, 
        //            0,   1,  2, 3, 4
        //
        
        // TODO: undefined curvature. see klnwrt
        const double curvature          = 1.;
        // TODO: undefined dzds. see klnwrt
        const double dzds               = 1.;
        
        // at starting point
        const double cosLamStart        = 1./sqrt(1. + dzds*dzds);
        const double sinLamStart        = cosLamStart*dzds;
        
        // at end point
        const double cosLamEnd          = cosLamStart;
        const double sinLamEnd          = sinLamStart;
        const double invCosLamEnd       = 1./cosLamEnd;
        
        // direction of magnetic field
        TVector3 BFieldDir(0.,0.,1.);
        
        // track direction vectors
        const double phiEnd = phi * ds * cosLamStart * curvature;
        TVector3 t1(cosLamStart*cos(phi), cosLamStart*sin(phi), sinLamStart);
        TVector3 t2(cosLamEnd*cos(phiEnd), cosLamEnd*sin(phiEnd), sinLamEnd);
        
        
        const double q                  = curvature*cosLamStart;
        const double theta              = q * ds;
        const double sinTheta           = sin(theta);
        const double cosTheta           = cos(theta);
        const double gamma              = BFieldDir.Dot(t2);    // (B,T)
        
        TVector3 an1                    = BFieldDir.Cross(t1);  // B x T0
        TVector3 an2                    = BFieldDir.Cross(t2);  // B x T
        
        // U0, V0
        const double au1                = 1./t1.Perp2();
        TVector3 u1( -au1*t1[1], au1*t1[0], 0. );
        TVector3 v1( t1[2]*u1[2], t1[2]*u1[0], t1[0]*u1[1] - t1[1]*u1[0] );
        
        // U, V
        const double au2                = 1./t2.Perp2();
        TVector3 u2( -au2*t2[1], au1*t2[0], 0. );
        TVector3 v2( t2[2]*u2[2], t2[2]*u2[0], t2[0]*u2[1] - t2[1]*u2[0] );
        
        const double qp                 = -bfac;
        const double pav                = qp / cosLamStart / curvature;
        const double anv                = -BFieldDir.Dot(u2);
        const double anu                = -BFieldDir.Dot(v2);
        const double omCosTheta         = 1. - cosTheta;
        const double tmSinTheta         = theta - sinTheta;
        
        
        _parPropJac.Zero();
             // q/p
        _parPropJac[0][0] = 1.;
        _parPropJac[0][0] = 1.;     

        return _parPropJac;
    }
        */
    void EUTelGBLFitter::SetTrackCandidates(const std::vector< EVENT::TrackerHitVec >& trackCandidates) {
        this->_trackCandidates = trackCandidates;
        return;
    }

    double EUTelGBLFitter::interpolateTrackX(const EVENT::TrackerHitVec& trackCand, const double z) const {
        double x0 = trackCand.front()->getPosition()[0];
        double z0 = trackCand.front()->getPosition()[2];

        double x = x0 - getTrackSlopeX(trackCand) * (z0 - z);

        return x;
    }

    //! Predict track hit in Y direction using simplified model

    double EUTelGBLFitter::interpolateTrackY(const EVENT::TrackerHitVec& trackCand, const double z) const {
        double y0 = trackCand.front()->getPosition()[1];
        double z0 = trackCand.front()->getPosition()[2];

        double y = y0 - getTrackSlopeY(trackCand) * (z0 - z);

        return y;
    }

    double EUTelGBLFitter::getTrackSlopeX(const EVENT::TrackerHitVec& trackCand) const {
        double x0 = trackCand.front()->getPosition()[0];
        double z0 = trackCand.front()->getPosition()[2];

        double xLast = trackCand.back()->getPosition()[0];
        double zLast = trackCand.back()->getPosition()[2];

        double kx = (x0 - xLast) / (z0 - zLast);

        return kx;
    }

    double EUTelGBLFitter::getTrackSlopeY(const EVENT::TrackerHitVec& trackCand) const {
        double y0 = trackCand.front()->getPosition()[1];
        double z0 = trackCand.front()->getPosition()[2];

        double yLast = trackCand.back()->getPosition()[1];
        double zLast = trackCand.back()->getPosition()[2];

        double ky = (y0 - yLast) / (z0 - zLast);

        return ky;
    }

    void EUTelGBLFitter::Reset() {
        if (!_gblTrackCandidates.empty()) {
            std::map< int, gbl::GblTrajectory* >::iterator it;
            for (it = _gblTrackCandidates.begin(); it != _gblTrackCandidates.begin(); ++it) delete it->second;
        }
        _gblTrackCandidates.clear();
//        _fittrackvec->clear();
        _hitId2GblPointLabel.clear();
    }

    /** Add a measurement to GBL point
     * 
     * @param point
     * @param meas measuremet vector (residuals) to be calculated in this routine
     * @param measPrec residuals weights (1/unc^2) to be calculated in this routine
     * @param hitpos hit position
     * @param xPred predicted by hit x-position (first approximation)
     * @param yPred predicted by hit y-position (first approximation)
     * @param hitcov hit covariance matrix
     * @param proL2m projection matrix from track coordinate system onto measurement system
     */
    void EUTelGBLFitter::addMeasurementsGBL(gbl::GblPoint& point, TVectorD& meas, TVectorD& measPrec, const double* hitpos,
            double xPred, double yPred, const EVENT::FloatVec& hitcov, TMatrixD& proL2m) {
        meas[0] = hitpos[0] - xPred;
        meas[1] = hitpos[1] - yPred;
        measPrec[0] = 1. / hitcov[0];
        measPrec[1] = 1. / hitcov[2];

        streamlog_out(DEBUG0) << "Residuals:" << std::endl;
        streamlog_out(DEBUG0) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] << std::endl;
        streamlog_out(DEBUG0) << "Y:" << std::setw(20) << meas[1] << std::setw(20) << measPrec[1] << std::endl;

        point.addMeasurement(proL2m, meas, measPrec);
    }

    /** Add a acatterer to GBL point
     * Add Si + Kapton thin scatterer to a point
     * 
     * @param point
     * @param scat average scattering angle
     * @param scatPrecSensor 1/RMS^2 of the scattering angle (ususally determined from Highland's formula)
     * @param iPlane plane id
     * @param p momentum of the particle
     */
    void EUTelGBLFitter::addScattererGBL(gbl::GblPoint& point, TVectorD& scat, TVectorD& scatPrecSensor, int planeID, double p) {
        const int iPlane = geo::gGeometry().sensorIDtoZOrder(planeID);
        const double radlenSi           = geo::gGeometry()._siPlanesLayerLayout->getSensitiveRadLength(iPlane);
        const double radlenKap          = geo::gGeometry()._siPlanesLayerLayout->getLayerRadLength(iPlane);
        const double thicknessSi        = geo::gGeometry()._siPlanesLayerLayout->getSensitiveThickness(iPlane);
        const double thicknessKap       = geo::gGeometry()._siPlanesLayerLayout->getLayerThickness(iPlane);

        const double X0Si = thicknessSi / radlenSi; // Si 
        const double X0Kap = thicknessKap / radlenKap; // Kapton                

        const double tetSi = Utility::getThetaRMSHighland(p, X0Si);
        const double tetKap = Utility::getThetaRMSHighland(p, X0Kap);

        scatPrecSensor[0] = 1.0 / (tetSi * tetSi + tetKap * tetKap);
        scatPrecSensor[1] = 1.0 / (tetSi * tetSi + tetKap * tetKap);

        point.addScatterer(scat, scatPrecSensor);
    }

    // @TODO iplane, xPred, yPred must not be here. consider refactoring

    /** Add alignment derivative necessary for MILLIPEDE
     * 
     * @param point
     * @param alDer matrix of global parameter (alignment constants) derivatives 
     * @param globalLabels vector of alignment parameters ids
     * @param iPlane plane id
     * @param xPred predicted by hit x-position (first approximation)
     * @param yPred predicted by hit y-position (first approximation)
     * @param xSlope predicted x-slope
     * @param ySlope predicted y-slope
     */
    void EUTelGBLFitter::addGlobalParametersGBL(gbl::GblPoint& point, TMatrixD& alDer, std::vector<int>& globalLabels, int iPlane,
            double xPred, double yPred, double xSlope, double ySlope) {
        alDer[0][0] = 1.0; // dx/dx
        alDer[0][1] = 0.0; // dx/dy
        alDer[1][0] = 0.0; // dy/dx
        alDer[1][1] = 1.0; // dy/dy
        globalLabels[0] = _paramterIdXShiftsMap[iPlane]; // dx
        globalLabels[1] = _paramterIdYShiftsMap[iPlane]; // dy
        if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][2] = -yPred; // dx/rot
            alDer[1][2] = xPred; // dy/rot
            globalLabels[2] = _paramterIdZRotationsMap[iPlane]; // rot z
        }
        if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][3] = xSlope; // dx/dz
            alDer[1][3] = ySlope; // dy/dz
            globalLabels[3] = _paramterIdZShiftsMap[iPlane]; // dz
        }
        if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][4] = xPred*xSlope; // dx/rot y
            alDer[1][4] = xPred*ySlope; // dy/rot y
            globalLabels[4] = _paramterIdYRotationsMap[iPlane]; // drot y
            alDer[0][5] = yPred*xSlope; // dx/rot x
            alDer[1][5] = yPred*ySlope; // dy/rot x
            globalLabels[5] = _paramterIdXRotationsMap[iPlane]; // drot x
        }
        if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
            alDer[0][3] = xPred*xSlope; // dx/rot y
            alDer[1][3] = xPred*ySlope; // dy/rot y
            globalLabels[3] = _paramterIdYRotationsMap[iPlane]; // drot y
        }
        if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
            alDer[0][3] = yPred*xSlope; // dx/rot x
            alDer[1][3] = yPred*ySlope; // dy/rot x
            globalLabels[3] = _paramterIdXRotationsMap[iPlane]; // drot x
        }
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = xPred*xSlope; // dx/rot y
            alDer[1][3] = xPred*ySlope; // dy/rot y
            globalLabels[3] = _paramterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = yPred*xSlope; // dx/rot x
            alDer[1][4] = yPred*ySlope; // dy/rot x
            globalLabels[4] = _paramterIdXRotationsMap[iPlane]; // drot x
        }
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = xPred*xSlope; // dx/rot y
            alDer[1][3] = xPred*ySlope; // dy/rot y
            globalLabels[3] = _paramterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = yPred*xSlope; // dx/rot x
            alDer[1][4] = yPred*ySlope; // dy/rot x
            globalLabels[4] = _paramterIdXRotationsMap[iPlane]; // drot x
        }

        point.addGlobals(globalLabels, alDer);
    }

    void EUTelGBLFitter::pushBachPoint( std::vector< gbl::GblPoint >& pointList, const gbl::GblPoint& point, int hitid ) {
        pointList.push_back(point);
        
        // store point's GBL label for future reference
        _hitId2GblPointLabel[ hitid ] = static_cast<int>(pointList.size());
    }
    /**
     * Set track omega, D0, Z0, Phi, tan(Lambda), Chi2, NDF parameters
     * and add it to the LCIO collection of fitted tracks
     * 
     * @param fittrack pointer to track object to be stored
     * @param chi2     Chi2 of the track fit
     * @param ndf      NDF of the track fit
     */
    void EUTelGBLFitter::prepareLCIOTrack( IMPL::TrackImpl* fittrack, gbl::GblTrajectory* gblTraj, const EVENT::TrackerHitVec& trackCandidate, 
                                          double chi2, int ndf, double omega, double d0, double z0, double phi, double tanlam ) {
        // prepare output collection
        try {
            _fittrackvec = new IMPL::LCCollectionVec( EVENT::LCIO::TRACK );
            IMPL::LCFlagImpl flag( _fittrackvec->getFlag( ) );
            flag.setBit( lcio::LCIO::TRBIT_HITS );
            _fittrackvec->setFlag( flag.getFlag( ) );
        } catch ( ... ) {
            streamlog_out( WARNING2 ) << "Can't allocate output collection" << std::endl;
        }
        
        // prepare track
        float refpoint[3] = { 0., 0., 0. };
        fittrack->setReferencePoint( refpoint );
        
        fittrack->setChi2 ( chi2 );      // Chi2 of the fit (including penalties)
        fittrack->setNdf  ( ndf );        // Number of planes fired (!)
        fittrack->setOmega( omega );       // curvature of the track
        fittrack->setD0   ( d0 );          // impact paramter of the track in (r-phi)
        fittrack->setZ0   ( z0 );          // impact paramter of the track in (r-z)
        fittrack->setPhi  ( phi );         // phi of the track at reference point
        fittrack->setTanLambda( tanlam );   // dip angle of the track at reference point
        
        unsigned int numData;
        TVectorD residual(200);
        TVectorD measErr(200);
        TVectorD residualErr(200);
        TVectorD downWeight(200);

        EVENT::TrackerHitVec::const_iterator itrHit;
        for ( itrHit = trackCandidate.begin(); itrHit != trackCandidate.end(); ++itrHit ) {
            
            const int hitGblLabel = _hitId2GblPointLabel[ (*itrHit)->id() ];
            IMPL::TrackerHitImpl* hit = new IMPL::TrackerHitImpl();
                  
            gblTraj->getMeasResults( hitGblLabel, numData, residual, measErr, residualErr, downWeight);
            
            // retrieve original hit coordinates
            double pos[3] = { (*itrHit)->getPosition()[0], (*itrHit)->getPosition()[1], (*itrHit)->getPosition()[2] };
            
            // correct original values to the fitted ones
            pos[0] -= residual[0];
            pos[1] -= residual[1];
            
            hit -> setPosition(pos);
            
            // prepare covariance matrix
            float cov[TRKHITNCOVMATRIX] = {0.,0.,0.,0.,0.,0.};
            cov[0]=measErr[0]*measErr[0];
            cov[2]=measErr[1]*measErr[1];
            
            hit -> setCovMatrix(cov);
            
            fittrack->addHit( hit );
        } // loop over track's hits

        // add track to LCIO collection vector
        _fittrackvec->push_back( fittrack );
    }
    
    void EUTelGBLFitter::FitTracks() {
        Reset(); // 

        streamlog_out(DEBUG2) << " EUTelGBLFitter::FitTracks() " << std::endl;
        streamlog_out(DEBUG1) << " N track candidates:" << static_cast<int>(_trackCandidates.size()) << std::endl;

        TVectorD meas(2);
        TVectorD measPrec(2); // precision = 1/resolution^2

        TVectorD scat(2);
        scat.Zero();

        TVectorD scatPrecSensor(2);
        TVectorD scatPrecAir(2);

        TMatrixD alDer; // alignment derivatives
        std::vector<int> globalLabels;
        if (_alignmentMode == Utility::XYShift) {
            globalLabels.resize(2);
            alDer.ResizeTo(2, 2);
        } else if (_alignmentMode == Utility::XYShiftXYRot) {
            globalLabels.resize(3);
            alDer.ResizeTo(2, 3);
        } else if (_alignmentMode == Utility::XYZShiftXYRot) {
            globalLabels.resize(4);
            alDer.ResizeTo(2, 4);
        } else if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
            globalLabels.resize(4);
            alDer.ResizeTo(2, 4);
        } else if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
            globalLabels.resize(4);
            alDer.ResizeTo(2, 4);
        } else if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            globalLabels.resize(5);
            alDer.ResizeTo(2, 5);
        } else if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            globalLabels.resize(6);
            alDer.ResizeTo(2, 6);
        }
        alDer.Zero();

        double p = _eBeam; // beam momentum
       
        std::vector< EVENT::TrackerHitVec >::const_iterator itTrkCand = _trackCandidates.begin();
        EVENT::TrackerHitVec::const_iterator itHit;
        for (; itTrkCand != _trackCandidates.end(); ++itTrkCand) {

            // sanity check. Mustn't happen in principle.
            if (itTrkCand->size() > geo::gGeometry().nPlanes()) continue;

            IMPL::TrackImpl * fittrack = new IMPL::TrackImpl();

            //GBL trajectory construction
            std::vector< gbl::GblPoint > pointList;

            TMatrixD jacPointToPoint(5, 5);
            jacPointToPoint.UnitMatrix();
            double step = 0.;
            itHit = itTrkCand->begin();
            for (; itHit != itTrkCand->end(); ++itHit) {
                const double* hitpos = (*itHit)->getPosition();
                const EVENT::FloatVec hitcov = (*itHit)->getCovMatrix();
                const int planeID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itHit) );

                double xPred = interpolateTrackX(*itTrkCand, hitpos[2]);
                double yPred = interpolateTrackY(*itTrkCand, hitpos[2]);

                double xSlope = getTrackSlopeX(*itTrkCand);
                double ySlope = getTrackSlopeY(*itTrkCand);

                gbl::GblPoint point(jacPointToPoint);
                TMatrixD proL2m(2, 2);
                proL2m.UnitMatrix();
                addMeasurementsGBL(point, meas, measPrec, hitpos, xPred, yPred, hitcov, proL2m);
                addScattererGBL(point, scat, scatPrecSensor, planeID, p);
                if (_alignmentMode != Utility::noAlignment) {
                    addGlobalParametersGBL( point, alDer, globalLabels, planeID, xPred, yPred, xSlope, ySlope );
                }
                pushBachPoint( pointList, point, (*itHit)->id() );

                // construct effective scatterers for air
                // the scatters must be at (Z(plane i) + Z(plane i+1))/2. +/- (Z(plane i) - Z(plane i+1))/sqrt(12)
                if ( itHit != (itTrkCand->end() - 1) ) {
                    const TVector3 startPoint( hitpos[0], hitpos[1], hitpos[2] );
                    const double* nexthitpos = (*(itHit + 1))->getPosition();
                    const TVector3 endPoint( nexthitpos[0], nexthitpos[1], nexthitpos[2] );
                    const TVector3 distBeteenHits = endPoint - startPoint;
                    const double hitSpacing = distBeteenHits.Mag(); 
//                    const double hitSpacing = fabs( hitpos[2] - (*(itHit + 1))->getPosition()[2] ); 

                    double X0Air =  hitSpacing / 304e3; // Air
                    double tetAir = Utility::getThetaRMSHighland(p, X0Air);
                    scatPrecAir[0] = 1.0 / (tetAir * tetAir);
                    scatPrecAir[1] = 1.0 / (tetAir * tetAir);

                    // propagate parameters into the air gap between consecutive planes 

                    {   // downstream air scatterer
                        step = hitSpacing / 2. - hitSpacing / sqrt(12.);
                        jacPointToPoint = propagatePar(step);
                        gbl::GblPoint pointInAir1(jacPointToPoint);
                        pointInAir1.addScatterer(scat, scatPrecAir);
                        pointList.push_back(pointInAir1);
                    }
                    
                    {   // upstream air scatterer
                        step = 2. * hitSpacing / sqrt( 12. );
                        jacPointToPoint = propagatePar( step );
                        gbl::GblPoint pointInAir2( jacPointToPoint );
                        pointInAir2.addScatterer( scat, scatPrecAir );
                        pointList.push_back( pointInAir2 );
                    }
                    
                    // propagate to the next hit
                    {
                        step = hitSpacing / 2. - hitSpacing / sqrt( 12. );
                        jacPointToPoint = propagatePar( step );
                    }
                } // if not the last hit
            } // loop over hits

            double loss = 0.;
            double chi2 = 0.;
            int ndf = 0;
            // perform GBL fit
            {
                gbl::GblTrajectory* traj = new gbl::GblTrajectory( pointList, false );

                int ierr = 0;
                if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( chi2, ndf, loss, _mEstimatorType );
                else ierr = traj->fit( chi2, ndf, loss );

                if ( chi2 < _chi2cut ) {
                    if ( !ierr ) traj->milleOut( *_mille );
                }
                
                _gblTrackCandidates.insert( std::make_pair(0, traj ) );
                
                // Write fit result
                {
                    prepareLCIOTrack( fittrack, traj, (*itTrkCand), chi2, ndf, 0., 0., 0., 0., 0. );
                }
            }
        } // loop over supplied track candidates

        return;
    } // EUTelGBLFitter::FitTracks()
} // namespace eutelescope

#endif
