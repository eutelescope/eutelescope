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
#include "EUTelUtilityRungeKutta.h"
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
    _fithitsvec(0),
    _parPropJac(5, 5),
    _beamQ(-1),
    _eBeam(-1.),
    _hitId2GblPointLabel(),
    _hitId2GblPointLabelMille(),
    _alignmentMode(Utility::XYShiftXYRot),
    _mEstimatorType(),
    _mille(0),
    _paramterIdXShiftsMap(),
    _paramterIdYShiftsMap(),
    _paramterIdZShiftsMap(),
    _paramterIdXRotationsMap(),
    _paramterIdYRotationsMap(),
    _paramterIdZRotationsMap(),
    _excludeFromFit(),
    _chi2cut(1000.),
    _eomIntegrator( new EUTelUtilityRungeKutta() ),
    _eomODE( 0 )
    {
                // Initialise ODE integrators for eom and jacobian
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
    }

    EUTelGBLFitter::EUTelGBLFitter(std::string name) : EUTelTrackFitter(name),
    _trackCandidates(),
    _gblTrackCandidates(),
    _fittrackvec(0),
    _fithitsvec(0),
    _parPropJac(5, 5),
    _beamQ(-1),
    _eBeam(-1.),
    _hitId2GblPointLabel(),
    _hitId2GblPointLabelMille(),
    _alignmentMode(Utility::XYShiftXYRot),
    _mEstimatorType(),
    _mille(0),
    _paramterIdXShiftsMap(),
    _paramterIdYShiftsMap(),
    _paramterIdZShiftsMap(),
    _paramterIdXRotationsMap(),
    _paramterIdYRotationsMap(),
    _paramterIdZRotationsMap(),
    _excludeFromFit(),
    _chi2cut(1000.),
    _eomIntegrator( new EUTelUtilityRungeKutta() ),
    _eomODE( 0 )
    {
                // Initialise ODE integrators for eom and jacobian
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
    }

    EUTelGBLFitter::~EUTelGBLFitter() {
    }
    
    void EUTelGBLFitter::setParamterIdPlaneVec( const std::vector<int>& vector)
    {
      _paramterIdPlaneVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdXResolutionVec( const std::vector<float>& vector)
    {
      _paramterIdXResolutionVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdYResolutionVec( const std::vector<float>& vector)
    {
      _paramterIdYResolutionVec = vector;
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
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdXRotationsMap() const {
        return _paramterIdXRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdYRotationsMap() const {
        return _paramterIdYRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdZRotationsMap() const {
        return _paramterIdZRotationsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdZShiftsMap() const {
        return _paramterIdZShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdYShiftsMap() const {
        return _paramterIdYShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdXShiftsMap() const {
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

    const std::map<long, int>& EUTelGBLFitter::getHitId2GblPointLabel( ) const {
        return _hitId2GblPointLabel;
    }

    void EUTelGBLFitter::setExcludeFromFitPlanes( const std::vector<int>& excludedPlanes ) {
        this->_excludeFromFit = excludedPlanes;
    }

    std::vector<int> EUTelGBLFitter::getExcludeFromFitPlanes( ) const {
        return _excludeFromFit;
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
        
    /** Propagation jacobian in inhomogeneous magnetic field
     * 
     * 
     * @param ds        Z - Z0 propagation distance
     * 
     * @return          transport matrix
     */
    
    TMatrixD EUTelGBLFitter::PropagatePar( double ds, double invP, double tx0, double ty0, double x0, double y0, double z0 ) {
        // for GBL:
        //   Jacobian for helical track
        //   track = q/p, x'  y'  x  y
        //            0,   1,  2, 3, 4
        //
        
        streamlog_out( DEBUG2 ) << "EUTelGBLFitter::PropagatePar()" << std::endl;
	// The formulas below are derived from equations of motion of the particle in
        // magnetic field under assumption |dz| small. Must be valid for |dz| < 20-30 cm

        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
        const double Bx         = B.at( vectorGlobal ).x();
        const double By         = B.at( vectorGlobal ).y();
        const double Bz         = B.at( vectorGlobal ).z();
	const double mm = 1000.;
	const double k = 0.299792458/mm;

        const double sqrtFactor = sqrt( 1. + tx0*tx0 + ty0*ty0 );

	const double Ax = sqrtFactor * (  ty0 * ( tx0 * Bx + Bz ) - ( 1. + tx0*tx0 ) * By );
	const double Ay = sqrtFactor * ( -tx0 * ( ty0 * By + Bz ) + ( 1. + ty0*ty0 ) * Bx );

	// Partial derivatives
	const double dAxdty0 = ty0 * Ax / (sqrtFactor*sqrtFactor) + sqrtFactor*( tx0*Bx + Bz );
	const double dAydtx0 = tx0 * Ay / (sqrtFactor*sqrtFactor) + sqrtFactor*( -ty0*By - Bz );

	const double dxdtx0 = ds;
	const double dxdty0 = 0.5 * invP * k * ds*ds * dAxdty0;

	const double dydtx0 = 0.5 * invP * k * ds*ds * dAydtx0;
	const double dydty0 = ds;

	const double dtxdty0 = invP * k * ds * dAxdty0;
	const double dtydtx0 = invP * k * ds * dAydtx0;

	const double dxdinvP0 = 0.5 * k * ds*ds * Ax;
	const double dydinvP0 = 0.5 * k * ds*ds * Ay;

	const double dtxdinvP0 = k * ds * Ax;
	const double dtydinvP0 = k * ds * Ay;

	// Fill-in matrix elements
	_parPropJac.UnitMatrix();
	_parPropJac[1][0] = dtxdinvP0;  _parPropJac[1][2] = dtxdty0;	
	_parPropJac[2][0] = dtydinvP0;  _parPropJac[2][1] = dtydtx0;	
	_parPropJac[3][0] = dxdinvP0;   _parPropJac[3][1] = dxdtx0;	_parPropJac[3][2] = dxdty0;	
	_parPropJac[4][0] = dydinvP0;   _parPropJac[4][1] = dydtx0;	_parPropJac[4][2] = dydty0;	
        
        if ( streamlog_level(DEBUG0) ){
             streamlog_out( DEBUG0 ) << "Propagation jacobian: " << std::endl;
            _parPropJac.Print();
        }
	
        if( streamlog_level(DEBUG4) ) 
        {
             streamlog_out(DEBUG4) << "step:" << ds ;
             streamlog_out(DEBUG4) << " to point x: " << std::setw(8) << x0 ;
             streamlog_out(DEBUG4) << " y: " << std::setw(8) << y0 ;
             streamlog_out(DEBUG4) << " z: " << std::setw(8) << z0 << std::endl;
             streamlog_out(DEBUG4) << " B=0:                  " ;
             streamlog_out(DEBUG4) << "                              and genB>0 : " << std::endl;
             for(int i =0;i<5;i++)
             {
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][0] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][1] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][2] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][3] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][4] ;
             streamlog_out(DEBUG4) << "                 " ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][0] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][1] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][2] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][3] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][4] <<std::endl;
             }
        }

        streamlog_out( DEBUG2 ) << "-----------------------------EUTelGBLFitter::PropagatePar()-------------------------------" << std::endl;   

        return _parPropJac;
    }

    /**
     * Get extrapolated position of the track in global coordinate system
     * Formulas are also valid for vanishing magnetic field.
     * Calculation is based on numerical integration of equations of motion
     * 
     * @param ts track state
     * @param dz propagation distance along z
     * @return vector of track parameters in the global coordinate system
     */
    TVectorD EUTelGBLFitter::getXYZfromDzNum( double invP, double tx, double ty, double x0, double y0, double z0, double dz ) const {
        streamlog_out(DEBUG2) << "EUTelGBLFitter::getXYZfromDzNum()" << std::endl;
        
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
        const double bx         = B.at( vectorGlobal ).x();
        const double by         = B.at( vectorGlobal ).y();
        const double bz         = B.at( vectorGlobal ).z();

        TVector3 hVec(bx,by,bz);
        
        // Setup the equation
        TVectorD trackStateVec(5);
	trackStateVec[0] = x0;
	trackStateVec[1] = y0;
	trackStateVec[2] = tx;
	trackStateVec[3] = ty;
	trackStateVec[4] = invP;
        static_cast< eom::EOMODE* >(_eomODE)->setInitValue( trackStateVec );
        static_cast< eom::EOMODE* >(_eomODE)->setBField( hVec );
        
        // Integrate
        TVectorD result = _eomIntegrator->integrate( dz );
        
        streamlog_out(DEBUG0) << "Result of the integration:" << std::endl;
        streamlog_message( DEBUG0, result.Print();, std::endl; );
        
        streamlog_out(DEBUG2) << "---------------------------------EUTelGBLFitter::getXYZfromDzNum()------------------------------------" << std::endl;
        
        streamlog_out(DEBUG4) << "getXYZfromDzNum   dz:" << dz ;
        streamlog_out(DEBUG4) << " to point x0: " << std::setw(8) << x0 ;
        streamlog_out(DEBUG4) << " y0: " << std::setw(8) << y0 ;
        streamlog_out(DEBUG4) << " z0: " << std::setw(8) << z0 ;
        streamlog_out(DEBUG4) << " tx: " << std::setw(8) << tx ;
        streamlog_out(DEBUG4) << " ty: " << std::setw(8) << ty << std::endl;
       
        return result;
    }
        
    void EUTelGBLFitter::SetTrackCandidates(const EVENT::TrackVec& trackCandidates) {

        this->_trackCandidates = trackCandidates;
        return;
    }

    double EUTelGBLFitter::interpolateTrackX(const EVENT::TrackerHitVec& trackCand, const double z) const {
        const int planeIDStart = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(trackCand.front()) );
        const double* hitPointLocalStart = trackCand.front()->getPosition();
        double hitposStart[] = {hitPointLocalStart[0],hitPointLocalStart[1],hitPointLocalStart[2]};
        double hitPointGlobalStart[] = {0.,0.,0.};
        geo::gGeometry().local2Master(planeIDStart,hitposStart,hitPointGlobalStart);
        
        double x0 = hitPointGlobalStart[0];
        double z0 = hitPointGlobalStart[2];

        double x = x0 - getTrackSlopeX(trackCand) * (z0 - z);

        return x;
    }
    

    double EUTelGBLFitter::interpolateTrackY(const EVENT::TrackerHitVec& trackCand, const double z) const {
        const int planeIDStart = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(trackCand.front()) );
        const double* hitPointLocalStart = trackCand.front()->getPosition();
        double hitposStart[] = {hitPointLocalStart[0],hitPointLocalStart[1],hitPointLocalStart[2]};
        double hitPointGlobalStart[] = {0.,0.,0.};
        geo::gGeometry().local2Master(planeIDStart,hitposStart,hitPointGlobalStart);
        
        double y0 = hitPointGlobalStart[1];
        double z0 = hitPointGlobalStart[2];

        double y = y0 - getTrackSlopeY(trackCand) * (z0 - z);

        return y;
    }
    
    double EUTelGBLFitter::getTrackSlopeX(const EVENT::TrackerHitVec& trackCand) const {
        const int planeIDStart = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(trackCand.front()) );
        const double* hitPointLocalStart = trackCand.front()->getPosition();
        double hitposStart[] = {hitPointLocalStart[0],hitPointLocalStart[1],hitPointLocalStart[2]};
        double hitPointGlobalStart[] = {0.,0.,0.};
        geo::gGeometry().local2Master(planeIDStart,hitposStart,hitPointGlobalStart);
        
        double x0 = hitPointGlobalStart[0];
        double z0 = hitPointGlobalStart[2];

        const int planeIDFinish = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(trackCand.back()) );
        const double* hitPointLocalFinish = trackCand.back()->getPosition();
        double hitposFinish[] = {hitPointLocalFinish[0],hitPointLocalFinish[1],hitPointLocalFinish[2]};
        double hitPointGlobalFinish[] = {0.,0.,0.};
        geo::gGeometry().local2Master(planeIDFinish,hitposFinish,hitPointGlobalFinish);
        
        double xLast = hitPointGlobalFinish[0];
        double zLast = hitPointGlobalFinish[2];

        double kx = (x0 - xLast) / (z0 - zLast);

        return kx;
    }
    
    double EUTelGBLFitter::getTrackSlopeY(const EVENT::TrackerHitVec& trackCand) const {
        const int planeIDStart = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(trackCand.front()) );
        const double* hitPointLocalStart = trackCand.front()->getPosition();
        double hitposStart[] = {hitPointLocalStart[0],hitPointLocalStart[1],hitPointLocalStart[2]};
        double hitPointGlobalStart[] = {0.,0.,0.};
        geo::gGeometry().local2Master(planeIDStart,hitposStart,hitPointGlobalStart);
        
        double y0 = hitPointGlobalStart[1];
        double z0 = hitPointGlobalStart[2];

        const int planeIDFinish = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(trackCand.back()) );
        const double* hitPointLocalFinish = trackCand.back()->getPosition();
        double hitposFinish[] = {hitPointLocalFinish[0],hitPointLocalFinish[1],hitPointLocalFinish[2]};
        double hitPointGlobalFinish[] = {0.,0.,0.};
        geo::gGeometry().local2Master(planeIDFinish,hitposFinish,hitPointGlobalFinish);
        
        double yLast = hitPointGlobalFinish[1];
        double zLast = hitPointGlobalFinish[2];

        double ky = (y0 - yLast) / (z0 - zLast);

        return ky;
    }

    void EUTelGBLFitter::Reset() {
        std::map< int, gbl::GblTrajectory* >::iterator it;
        for (it = _gblTrackCandidates.begin(); it != _gblTrackCandidates.end(); ++it) delete it->second;
        _gblTrackCandidates.clear();
        
        _hitId2GblPointLabel.clear();
        _hitId2GblPointLabelMille.clear();       
    }

    /** Add a measurement to GBL point
     * 
     * @param point
     * @param meas measuremet vector (residuals) to be calculated in this routine
     * @param measPrec residuals weights (1/unc^2) to be calculated in this routine
     * @param hitpos hit position
     * @param predpos predicted by hit x-position (approximation)
     * @param hitcov hit covariance matrix
     * @param proL2m projection matrix from track coordinate system onto measurement system
     */
    void EUTelGBLFitter::addMeasurementsGBL(gbl::GblPoint& point, TVectorD& meas, TVectorD& measPrec, const double* hitpos,
            const double* predpos, const EVENT::FloatVec& hitcov, TMatrixD& proL2m) {
     
        streamlog_out(DEBUG3) << " addMeasurementsGBL " << std::endl;
 
        meas[0] = hitpos[0] - predpos[0];
        meas[1] = hitpos[1] - predpos[1];
        measPrec[0] = 1. / hitcov[0];	// cov(x,x)
        measPrec[1] = 1. / hitcov[2];	// cov(y,y)

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
    void EUTelGBLFitter::addSiPlaneScattererGBL(gbl::GblPoint& point, TVectorD& scat, TVectorD& scatPrecSensor, int planeID, double p) {
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
            const double* predpos, double xSlope, double ySlope) {

        streamlog_out(MESSAGE1) << " addGlobalPArametersGBL " << std::endl;

        alDer[0][0] = -1.0; // dx/dx
        alDer[0][1] =  0.0; // dx/dy
        alDer[1][0] =  0.0; // dy/dx
        alDer[1][1] = -1.0; // dy/dy
        globalLabels[0] = _paramterIdXShiftsMap[iPlane]; // dx
        globalLabels[1] = _paramterIdYShiftsMap[iPlane]; // dy


        if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][2] = -predpos[1]; // dx/rot
            alDer[1][2] =  predpos[0]; // dy/rot
            globalLabels[2] = _paramterIdZRotationsMap[iPlane]; // rot z
        }


        if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][3] =   -xSlope; // dx/dz
            alDer[1][3] =   -ySlope; // dy/dz
            globalLabels[3] = _paramterIdZShiftsMap[iPlane]; // dz
        }

        if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][4] =   predpos[0]*xSlope; // dx/rot y
            alDer[1][4] =   predpos[0]*ySlope; // dy/rot y
            globalLabels[4] = _paramterIdYRotationsMap[iPlane]; // drot y  - actually X?
            alDer[0][5] =  -predpos[1]*xSlope; // dx/rot x          
            alDer[1][5] =  -predpos[1]*ySlope; // dy/rot x         
            globalLabels[5] = _paramterIdXRotationsMap[iPlane]; // drot x  - actually Y?
        }


// partial alignment 
        if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _paramterIdYRotationsMap[iPlane]; // drot y
        }
        if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
            alDer[0][3] = predpos[1]*xSlope; // dx/rot x
            alDer[1][3] = predpos[1]*ySlope; // dy/rot x
            globalLabels[3] = _paramterIdXRotationsMap[iPlane]; // drot x
        }
 
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _paramterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _paramterIdXRotationsMap[iPlane]; // drot x
        }
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _paramterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _paramterIdXRotationsMap[iPlane]; // drot x
        }

        point.addGlobals(globalLabels, alDer);
 
        streamlog_out(MESSAGE1) << " addGlobalPArametersGBL over " << std::endl;
   }

    void EUTelGBLFitter::pushBackPoint( std::vector< gbl::GblPoint >& pointList, const gbl::GblPoint& point, int hitid ) {
        pointList.push_back(point);
        
        // store point's GBL label for future reference
        _hitId2GblPointLabel[ hitid ] = static_cast<int>(pointList.size());
    }

    void EUTelGBLFitter::pushBackPointMille( std::vector< gbl::GblPoint >& pointList, const gbl::GblPoint& point, int hitid ) {
        pointList.push_back(point);
        
        // store point's GBL label for future reference
        _hitId2GblPointLabelMille[ hitid ] = static_cast<int>(pointList.size());
    }
    /**
     * Set track omega, D0, Z0, Phi, tan(Lambda), Chi2, NDF parameters
     * and add it to the LCIO collection of fitted tracks
     * 
     * @param fittrack pointer to track object to be stored
     * @param chi2     Chi2 of the track fit
     * @param ndf      NDF of the track fit
     */
    void EUTelGBLFitter::prepareLCIOTrack( gbl::GblTrajectory* gblTraj, const EVENT::TrackerHitVec& trackCandidate, 
                                          double chi2, int ndf, double omega, double d0, double z0, double phi, double tanlam ) {
        
        IMPL::TrackImpl * fittrack = new IMPL::TrackImpl();        
         
        unsigned int numData;
        TVectorD corrections(5);
        TMatrixDSym correctionsCov(5);
        
        TVectorD residual(2);
        TVectorD measErr(2);
        TVectorD residualErr(2);
        TVectorD downWeight(2);
       
        float cov[TRKHITNCOVMATRIX] = {0.,0.,0.,0.,0.,0.};
 
        float trackRefPoint[3] = { 0., 0., 0. };
        double invP = 1./_eBeam;

        TVectorD prevState = getXYZfromDzNum( invP, 0., 0., 0., 0., 0., 0. );
//        double prevZ = refPoint[2];

        EVENT::TrackerHitVec::const_iterator itrHit;
        for ( itrHit = trackCandidate.begin(); itrHit != trackCandidate.end(); ++itrHit ) {
            const int planeID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itrHit) );
             
            bool excludeFromFit = false;
            if ( std::find( _excludeFromFit.begin(), _excludeFromFit.end(), planeID ) != _excludeFromFit.end() ) excludeFromFit = true;
            if(excludeFromFit)
            {
               streamlog_out(DEBUG3) << "this plane is marked as the one to be skipped -> ID: " << planeID << std::endl;
//               continue;
            }  

            const int hitGblLabel = _hitId2GblPointLabel[ (*itrHit)->id() ];
                  
            gblTraj->getResults( hitGblLabel, corrections, correctionsCov );
            
            streamlog_out(DEBUG2) << "Corrections: " << std::endl;
            streamlog_message( DEBUG2, corrections.Print();, std::endl; );
/// get measured hit position:
// in local:
            const double* hitPointLocal = (*itrHit)->getPosition();
            double hitPointGlobal[] = {0.,0.,0.};
// in global:
            geo::gGeometry().local2Master(planeID,hitPointLocal,hitPointGlobal);

///
            
            // retrieve original hit coordinates
//            double pos[3] = { (*itrHit)->getPosition()[0], (*itrHit)->getPosition()[1], (*itrHit)->getPosition()[2] };
 
            streamlog_out(DEBUG2) << "hit " << hitGblLabel  << " at " << hitPointGlobal[0] << " " << hitPointGlobal[1] << " " << hitPointGlobal[2] << std::endl;
 
            double trackPointLocal[] = { hitPointLocal[0], hitPointLocal[1], hitPointLocal[2] };
            double trackPointGlobal[] = { 0., 0., 0. };
          
            if( !excludeFromFit)
            {
              gblTraj->getMeasResults( hitGblLabel, numData, residual, measErr, residualErr, downWeight);
              // correct original values to the fitted ones
              trackPointLocal[0] -= residual[0];
              trackPointLocal[1] -= residual[1];
              // prepare covariance matrix
              cov[0]=residualErr[0]*residualErr[0];
              cov[2]=residualErr[1]*residualErr[1];
            }

            geo::gGeometry().local2Master( planeID, trackPointLocal, trackPointGlobal );



            streamlog_out(DEBUG2) << "fit " << hitGblLabel  << " at " << trackPointGlobal[0] << " " << trackPointGlobal[1] << " " << trackPointGlobal[2] << std::endl;

// check if it s not the last one, to get distance to the next one and propagate from the track measurement.
            if ( itrHit != (trackCandidate.end() - 1) ) {
            
// here we assume that trackCandidate does not contain hits from same plane ! should be checked explicitely somewhere ...
                    const int nextPlaneID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*(itrHit + 1)) );
                    const double* nextHitPoint = (*(itrHit + 1))->getPosition();
                    double nextHitPointGlobal[] = { 0., 0., 0. };
                    geo::gGeometry().local2Master(nextPlaneID,nextHitPoint,nextHitPointGlobal);
                    
// correct only for linear tracks ?
                    const double hitSpacing = sqrt( (nextHitPointGlobal[0]- trackPointGlobal[0])*(nextHitPointGlobal[0]-trackPointGlobal[0])+
                                                           (nextHitPointGlobal[1]-trackPointGlobal[1])*(nextHitPointGlobal[1]-trackPointGlobal[1])+
                                                           (nextHitPointGlobal[2]-trackPointGlobal[2])*(nextHitPointGlobal[2]-trackPointGlobal[2]) ); 
                    const double hitSpacingDz = -nextHitPointGlobal[2] +trackPointGlobal[2];
 
	            double trackDirLocal[] = { corrections[1], corrections[2], 0. };
		    double trackDirGlobal[] = { 0., 0., 0. };
	            geo::gGeometry().local2MasterVec( planeID, trackDirLocal, trackDirGlobal);

                    prevState = getXYZfromDzNum( invP, trackDirGlobal[0], trackDirGlobal[1], trackPointGlobal[0], trackPointGlobal[1], trackPointGlobal[2], hitSpacingDz );
                    streamlog_out(DEBUG2) << "forward "  << prevState[0] << " "  << prevState[1] << " " << trackPointGlobal[2] << " " << hitSpacingDz << std::endl;

//          jacPointToPoint = PropagatePar( step, invP, corrections[3], corrections[4], corrections[1], corrections[2], hitPointGlobal[2] );
            }
           
            if ( itrHit == trackCandidate.begin() ) {
//              double hitPosGlobal[] = {0.,0.,0.};
//                geo::gGeometry().local2Master( Utility::GuessSensorID(*itrHit), pos, hitPosGlobal);
//                geo::gGeometry().local2Master( planeID, pos, hitPosGlobal);
                trackRefPoint[0] = hitPointGlobal[0];     trackRefPoint[1] = hitPointGlobal[1];     trackRefPoint[2] = hitPointGlobal[2];
                omega +=  corrections[0];
            }

            IMPL::TrackerHitImpl* hit = new IMPL::TrackerHitImpl();            
            hit -> setPosition( trackPointLocal );
            EVENT::LCObjectVec originalHit;
            originalHit.push_back( *itrHit );
            hit -> rawHits() = originalHit;
           
            hit -> setCovMatrix(cov);
            
            fittrack->addHit( hit );
	    _fithitsvec->addElement( hit ); // fitted hit needs to be stored in a collection separately
        } // loop over track hits

        // prepare track
        fittrack->setReferencePoint( trackRefPoint );
        fittrack->setChi2 ( chi2 );      // Chi2 of the fit (including penalties)
        fittrack->setNdf  ( ndf );        // Number of planes fired (!)
        fittrack->setOmega( omega );       // curvature of the track
        fittrack->setD0   ( d0 );          // impact parameter of the track in (r-phi)
        fittrack->setZ0   ( z0 );          // impact parameter of the track in (r-z)
        fittrack->setPhi  ( phi );         // phi of the track at reference point
        fittrack->setTanLambda( tanlam );   // dip angle of the track at reference point
        
        // add track to LCIO collection vector
        _fittrackvec->addElement( fittrack );
    }
    
    void EUTelGBLFitter::FitTracks() {
        Reset(); // 

        // prepare output collection
        try {
            _fittrackvec = new IMPL::LCCollectionVec( EVENT::LCIO::TRACK );
            _fithitsvec = new IMPL::LCCollectionVec( EVENT::LCIO::TRACKERHIT );
            IMPL::LCFlagImpl flag( _fittrackvec->getFlag( ) );
            flag.setBit( lcio::LCIO::TRBIT_HITS );
            _fittrackvec->setFlag( flag.getFlag( ) );
        } catch ( ... ) {
            streamlog_out( ERROR2 ) << "Can't allocate output collection" << std::endl;
        }
        
        streamlog_out(DEBUG2) << " EUTelGBLFitter::FitTracks() " << std::endl;
        streamlog_out(DEBUG1) << " N track candidates:" << static_cast<int>(_trackCandidates.size()) << std::endl;

        TVectorD meas(2);
        TVectorD measPrec(2); // precision = 1/resolution^2

        TVectorD scat(2);
        scat.Zero();

        TVectorD scatPrecSensor(2);
        //TMatrixDSym scatPrec(2,2);
        TVectorD scatPrec(2);

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
       

        EVENT::TrackVec::const_iterator itTrkCand;

        for ( itTrkCand = _trackCandidates.begin(); itTrkCand != _trackCandidates.end(); ++itTrkCand) {
            // sanity check. Mustn't happen in principle.
            if ((*itTrkCand)->getTrackerHits().size() > geo::gGeometry().nPlanes()) continue;

            //GBL trajectory construction
            std::vector< gbl::GblPoint > pointList;

            TMatrixD jacPointToPoint(5, 5);
            jacPointToPoint.UnitMatrix();

            double step = 0.;

	    // Z axis points along beam direction.
	    double pt = ( 1./(*itTrkCand)->getOmega() ) * _beamQ;
	    double px = (*itTrkCand)->getTanLambda() * pt;
	    double py = pt * sin( (*itTrkCand)->getPhi() );
	    double pz = pt * cos( (*itTrkCand)->getPhi() );

	    double tx   = px / pz;
	    double ty   = py / pz;
	    double invP = _beamQ / sqrt( px*px + pt*pt );
 
            streamlog_out(DEBUG4) << "FitTracks   ";
            streamlog_out(DEBUG4) << " omega: " << std::setw(8) << (*itTrkCand)->getOmega()  ;
            streamlog_out(DEBUG4) << " Q: " << std::setw(8) << _beamQ ;
            streamlog_out(DEBUG4) << " px: " << std::setw(8) << px ;
            streamlog_out(DEBUG4) << " py: " << std::setw(8) << py ;
            streamlog_out(DEBUG4) << " pz: " << std::setw(8) << pt << std::endl;
 
	    // Reference point is the last hit on a track
	    const float *refPoint = (*itTrkCand)->getReferencePoint();

	    TVectorD prevState = getXYZfromDzNum( invP, tx, ty, refPoint[0], refPoint[1], refPoint[2], 0. );
	    double prevZ = refPoint[2];
            
            // Loop over hits on a track candidate
            const EVENT::TrackerHitVec& hits = (*itTrkCand)->getTrackerHits();
            EVENT::TrackerHitVec::const_reverse_iterator itHit;
            for ( itHit = hits.rbegin(); itHit != hits.rend(); ++itHit) {
                const int planeID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itHit) );
		if ( planeID < 0 ) streamlog_out( WARNING2 ) << "Can't guess sensor ID. Check supplied hits." << std::endl;

                // Go to global coordinates
                const double* hitPointLocal = (*itHit)->getPosition();
                double hitPointGlobal[] = {0.,0.,0.};
                geo::gGeometry().local2Master(planeID,hitPointLocal,hitPointGlobal);


//                const EVENT::FloatVec hitcov = (*itHit)->getCovMatrix();
                EVENT::FloatVec hitcov(4);
                hitcov[0]=0.01;
                hitcov[1]=0.00;
                hitcov[2]=0.01;
                hitcov[3]=0.00;

                // check Processor Parameters for plane resolution // Denys
                if( _paramterIdXResolutionVec.size() > 0 && _paramterIdYResolutionVec.size() > 0 )
                {
                  for(int izPlane=0;izPlane<_paramterIdPlaneVec.size();izPlane++)
                  {
                    if( _paramterIdPlaneVec[izPlane] == planeID )
                    {  
                      hitcov[0] =  _paramterIdXResolutionVec[izPlane];
                      hitcov[2] =  _paramterIdYResolutionVec[izPlane];
 
                      hitcov[0] *= hitcov[0]; // squared !
                      hitcov[2] *= hitcov[2]; // squared !
                      break;
                    } 
                  }
                }
                else
                {  
                  hitcov = (*itHit)->getCovMatrix();            
                  
                }

                streamlog_out(DEBUG0) << "Hit covariance matrix: [0: "  << hitcov[0] << "] [1: "  << hitcov[1] << "] [2: "  << hitcov[2] << "] [3: " << hitcov[3]  << std::endl;

//                double xPred = interpolateTrackX(hits, hitPointGlobal[2]);
//                double yPred = interpolateTrackY(hits, hitPointGlobal[2]);

		const double dz = hitPointGlobal[2] - prevZ;
//                prevZ = hitPointGlobal[2];


//		double hitPointLocal2[] = {0.,0.,0.};
//                geo::gGeometry().master2Local( hitPointGlobal, hitPointLocal2 );
//		streamlog_out(MESSAGE4) << "hitl2= " << hitPointLocal2[0] << " " << hitPointLocal2[1] << " " << hitPointLocal2[2] << std::endl;

		streamlog_out(DEBUG4) << "hitl= " << hitPointLocal[0] << " " << hitPointLocal[1] << " " << hitPointLocal[2] << std::endl;
		streamlog_out(DEBUG4) << "hitm= " << hitPointGlobal[0] << " " << hitPointGlobal[1] << " " << hitPointGlobal[2] << std::endl;
		streamlog_out(DEBUG4) << "dz= " << dz << std::endl;
//		prevState.Print();
		
		// Propagate track parameters to the next hit
// rubinskiy		TVectorD trackParamPrediction = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], refPoint[2], dz );
                TVectorD trackParamPrediction = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ, dz );
//		trackParamPrediction.Print();
                prevZ = hitPointGlobal[2];

                prevState[0] =   trackParamPrediction[0];
                prevState[1] =   trackParamPrediction[1];
                prevState[2] =   trackParamPrediction[2];

		// Perform transformation to the sensor local coordinates
		double trackPointGlobal[] = { trackParamPrediction[0], trackParamPrediction[1], hitPointGlobal[2] };
		double trackPointLocal[] = { 0., 0., 0. };
                geo::gGeometry().master2Local( trackPointGlobal, trackPointLocal );

		double trackDirGlobal[] = { trackParamPrediction[2], trackParamPrediction[3], 1. };
		double trackDirLocal[] = { 0., 0., 0. };
		geo::gGeometry().master2LocalVec( planeID, trackDirGlobal, trackDirLocal );

                // print track parameters:
		streamlog_out(DEBUG4) << "dir trkg=  " << trackDirGlobal[0] << " " << trackDirGlobal[1] << " " << trackDirGlobal[2] << std::endl;
		streamlog_out(DEBUG4) << "dir trkl=  " << trackDirLocal [0] << " " << trackDirLocal [1] << " " << trackDirLocal [2] << std::endl;

		streamlog_out(DEBUG4) << "pos trkg=  " << trackPointGlobal[0] << " " << trackPointGlobal[1] << " " << trackPointGlobal[2] << std::endl;
		streamlog_out(DEBUG4) << "pos trkl=  " << trackPointLocal [0] << " " << trackPointLocal [1] << " " << trackPointLocal [2] << std::endl;


                gbl::GblPoint point(jacPointToPoint);
//                jacPointToPoint.Print();

		// Calculate projection matrix
                TMatrixD proL2m(2, 2);

		const TGeoHMatrix* globalH = geo::gGeometry().getHMatrix( hitPointGlobal );
		const TGeoHMatrix& globalHInv = globalH->Inverse();
		const double* rotation = globalHInv.GetRotationMatrix();

		proL2m[0][0] = rotation[0]; // x projection, xx
		proL2m[0][1] = rotation[1]; // y projection, xy
		proL2m[1][0] = rotation[3]; // x projection, yx
		proL2m[1][1] = rotation[4]; // y projection, yy

//		proL2m.UnitMatrix();
//		proL2m.Print();
                
                bool excludeFromFit = false;
                if ( std::find( _excludeFromFit.begin(), _excludeFromFit.end(), planeID ) != _excludeFromFit.end() ) excludeFromFit = true;
 
                if ( !excludeFromFit )
                {
                     addMeasurementsGBL(point, meas, measPrec, hitPointLocal, trackPointLocal, hitcov, proL2m);
//                   addMeasurementsGBL(point, meas, measPrec, hitPointGlobal, trackPointGlobal, hitcov, proL2m);
                  if (_alignmentMode != Utility::noAlignment) {
//                    if ( !excludeFromFit ) addGlobalParametersGBL( point, alDer, globalLabels, planeID, trackPointLocal, trackDirLocal[0], trackDirLocal[1] );
//                      addGlobalParametersGBL( point, alDer, globalLabels, planeID, trackPointLocal, trackDirLocal[0], trackDirLocal[1] );
                  }
                }

                addSiPlaneScattererGBL(point, scat, scatPrecSensor, planeID, p);
               
                pushBackPoint( pointList, point, (*itHit)->id() );

                // construct effective scatters for air
                // the scatters must be at (Z(plane i) + Z(plane i+1))/2. +/- (Z(plane i) - Z(plane i+1))/sqrt(12)
                if ( itHit != (hits.rend() - 1) ) {

                    // Go to global coordinates
                    const int nextPlaneID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*(itHit + 1)) );
                    const double* nextHitPoint = (*(itHit + 1))->getPosition();
                    double nextHitPointGlobal[] = { 0., 0., 0. };
                    geo::gGeometry().local2Master(nextPlaneID,nextHitPoint,nextHitPointGlobal);
                    
                    const double hitSpacing = sqrt( (nextHitPointGlobal[0]-hitPointGlobal[0])*(nextHitPointGlobal[0]-hitPointGlobal[0])+
                                                           (nextHitPointGlobal[1]-hitPointGlobal[1])*(nextHitPointGlobal[1]-hitPointGlobal[1])+
                                                           (nextHitPointGlobal[2]-hitPointGlobal[2])*(nextHitPointGlobal[2]-hitPointGlobal[2]) ); 
                    const double hitSpacingDz = nextHitPointGlobal[2]-hitPointGlobal[2];
                    
		    // Check if start and finish points are inside silicon sensor planes. Sometimes it happens because of numerical rounding.
		    double rad = geo::gGeometry().findRadLengthIntegral( hitPointGlobal, nextHitPointGlobal, true );
                    streamlog_out(DEBUG0) << "Rad length( " << planeID << "-" << nextPlaneID << "): " << rad << std::endl;
		    if ( rad < 1e-10 ) {
			streamlog_out(WARNING3) << "Rad length is too small for thick GBL scatterer( " << planeID << "-" << nextPlaneID << "): " << rad << std::endl;
			streamlog_out(WARNING3) << "Check the calculations" << std::endl;
			streamlog_out(WARNING3) << "Setting radiation length to " << 1e-10 << std::endl;
			rad = 1e-10;
		    }
                    
                    double sigmaTheta = Utility::getThetaRMSHighland(p, rad);

		    // Calculate scalar products of local offset directions and track direction
		    double uDir[] = { 1., 0., 0. };
		    double uDirGlobal[] = { 0., 0., 0. };
		    geo::gGeometry().local2Master(planeID,uDir,uDirGlobal);
		    double vDir[] = { 0., 1., 0. };
		    double vDirGlobal[] = { 0., 0., 0. };
		    geo::gGeometry().local2Master(planeID,vDir,vDirGlobal);

		    // Calculate px, py, pz at current point
//		    const double px = p*trackParamPrediction[2] / sqrt( 1. + trackParamPrediction[2]*trackParamPrediction[2] +
//									     trackParamPrediction[3]*trackParamPrediction[3] );
//                    const double py = p*trackParamPrediction[3] / sqrt( 1. + trackParamPrediction[2]*trackParamPrediction[2] + 
//									     trackParamPrediction[3]*trackParamPrediction[3] );
//                    const double pz = p    / sqrt( 1. + trackParamPrediction[2]*trackParamPrediction[2] + trackParamPrediction[3]*trackParamPrediction[3] );			  const double p = sqrt ( px*px + py*py + pz*pz );
//
//		    double c1, c2;
//		    c1 = uDirGlobal[0]*px + uDirGlobal[1]*py + uDirGlobal[2]*pz;
//		    c1 /= p;
//		    c2 = vDirGlobal[0]*px + vDirGlobal[1]*py + vDirGlobal[2]*pz;
//		    c2 /= p;
//
//                    scatPrec[0][0] = 1. - c1*c1;
//                    scatPrec[0][1] = -c1*c2;
//                    scatPrec[1][0] = -c1*c2;
//                    scatPrec[1][1] = 1. - c2*c2;
//		    scatPrec *= 1.0 / (sigmaTheta * sigmaTheta) * (1. - c1*c1 - c2*c2);

		    scatPrec[0] = 1.0 / (sigmaTheta * sigmaTheta);
		    scatPrec[1] = 1.0 / (sigmaTheta * sigmaTheta);

                    streamlog_out(DEBUG0) << "Scattering between hits:" << std::endl;
//                    scatPrec.Print();
                    
                    // propagate parameters into the air gap between consecutive planes 

                    {   // downstream air scatterer
                        step = hitSpacingDz / 2. - hitSpacingDz / sqrt(12.);
// rubinsky                        prevZ += step;
                        jacPointToPoint = PropagatePar( step, invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ );
//rubinsky                       jacPointToPoint = PropagatePar( step, invP, trackParamPrediction[2], trackParamPrediction[3], trackParamPrediction[0], trackParamPrediction[1], hitPointGlobal[2] );
//                        jacPointToPoint.Print();
// rubinskiy           prevState = getXYZfromDzNum( invP, trackParamPrediction[2], trackParamPrediction[3], trackParamPrediction[0], trackParamPrediction[1], hitPointGlobal[2], step );
                       prevState = getXYZfromDzNum( invP, trackParamPrediction[2], trackParamPrediction[3], trackParamPrediction[0], trackParamPrediction[1], prevZ, step );
                       prevZ += step;

                        gbl::GblPoint pointInAir1(jacPointToPoint);
                        pointInAir1.addScatterer(scat, scatPrec);
                        pointList.push_back(pointInAir1);
                    }
                    
                    {   // upstream air scatterer
                        step = 2*hitSpacingDz / sqrt( 12. ); // rubinskiy
// rubinsky                        prevZ += step;
// rubinsky                         jacPointToPoint = PropagatePar( step, invP, prevState[2], prevState[3], prevState[0], prevState[1], hitPointGlobal[2] );
                        jacPointToPoint = PropagatePar( step, invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ );
//                        jacPointToPoint.Print();
//rubinskiy             prevState = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], hitPointGlobal[2], step );
                        prevState = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ, step );
                        prevZ += step;
                        gbl::GblPoint pointInAir2( jacPointToPoint );
                        pointInAir2.addScatterer( scat, scatPrec );
                        pointList.push_back( pointInAir2 );
                    }
                    
                    // propagate to the next hit
                    {
                        step = hitSpacingDz / 2. - hitSpacingDz / sqrt( 12. );
// not needed after this step anymore                         prevZ += step;
// rubinskiy            jacPointToPoint = PropagatePar( step, invP, prevState[2], prevState[3], prevState[0], prevState[1], hitPointGlobal[2] );
                        jacPointToPoint = PropagatePar( step, invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ );
//                        jacPointToPoint.Print();
                    }
                } // if not the last hit
            } // loop over hits

            double loss = 0.;
            double chi2 = 0.;
            int ndf = 0;
            // perform GBL fit
            {
                // check magnetic field at (0,0,0)
                const gear::BField&   B = geo::gGeometry().getMagneticFiled();
                const double Bmag       = B.at( TVector3(0.,0.,0.) ).r2();
                
                gbl::GblTrajectory* traj = 0;
                if ( Bmag < 1.E-6 ) {
                   traj = new gbl::GblTrajectory( pointList, false );
                } else {
                   traj = new gbl::GblTrajectory( pointList, true );
                }
                
//                if ( !traj.isValid() ) {
//                        streamlog_out(WARNING1) << "Not valid GBL trajectory. Check the input." << std::endl;
//                        streamlog_message( MESSAGE0, traj->printData(); traj->printPoints(1); traj->printTrajectory(1);, std::endl;); 
//                    continue;
//                }
                
                int ierr = 0;

                if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( chi2, ndf, loss, _mEstimatorType );
                else ierr = traj->fit( chi2, ndf, loss );

              if ( chi2 < _chi2cut ) 
                {
                    if ( ierr )
                    {
			traj->printTrajectory(1);
			traj->printData();
			traj->printPoints(1);
		    }
                } 
                
                 EVENT::TrackVec::const_iterator begin = _trackCandidates.begin();
                _gblTrackCandidates.insert( std::make_pair( std::distance( begin, itTrkCand ), traj ) );
                
                // Write fit result
                {
                    prepareLCIOTrack( traj, (*itTrkCand)->getTrackerHits(), chi2, ndf, invP, 0., 0., 0., 0. );
                }

                // Prepare and write Mille Out
                if (_alignmentMode != Utility::noAlignment) 
                {
//                     prepareMilleOut( traj, (*itTrkCand)->getTrackerHits(), chi2, ndf, invP, 0., 0., 0., 0. );
                   prepareMilleOut( traj, itTrkCand, chi2, ndf, invP, 0., 0., 0., 0. );
                }
            }
        } // loop over supplied track candidates

        return;
    } // EUTelGBLFitter::FitTracks()


    void EUTelGBLFitter::prepareMilleOut( gbl::GblTrajectory* gblTraj, const EVENT::TrackVec::const_iterator& itTrkCand, 
                                          double chi2, int ndf, double omega, double d0, double z0, double phi, double tanlam ) {

        const EVENT::TrackerHitVec trackCandidate = (*itTrkCand)->getTrackerHits();

// now get starting point:
          TMatrixD jacPointToPoint(5, 5);
          jacPointToPoint.UnitMatrix();

          double step = 0.;

          // Z axis points along beam direction.
          double pt = ( 1./(*itTrkCand)->getOmega() ) * _beamQ;
          double px = (*itTrkCand)->getTanLambda() * pt;
	  double py = pt * sin( (*itTrkCand)->getPhi() );
	  double pz = pt * cos( (*itTrkCand)->getPhi() );

	  double tx   = px / pz;
	  double ty   = py / pz;
	  double invP = _beamQ / sqrt( px*px + pt*pt );
 
          const float *refPoint = (*itTrkCand)->getReferencePoint();

          TVectorD prevState = getXYZfromDzNum( invP, tx, ty, refPoint[0], refPoint[1], refPoint[2], 0. );
          double prevZ = refPoint[2];

          streamlog_out(DEBUG4) << "FitTracks   ";
          streamlog_out(DEBUG4) << " omega: " << std::setw(8) << (*itTrkCand)->getOmega()  ;
          streamlog_out(DEBUG4) << " Q: " << std::setw(8) << _beamQ ;
          streamlog_out(DEBUG4) << " px: " << std::setw(8) << px ;
          streamlog_out(DEBUG4) << " py: " << std::setw(8) << py ;
          streamlog_out(DEBUG4) << " pz: " << std::setw(8) << pt << std::endl;
 
          streamlog_out(DEBUG3) << "refPoint "  << refPoint[0] << " "  << refPoint[1] << " " << refPoint[2]  << std::endl;

          std::vector< gbl::GblPoint > pointList;

//        TMatrixD jacPointToPoint(5, 5);
//        jacPointToPoint.UnitMatrix();

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


         
        unsigned int numData;
        TVectorD corrections(5);
        TMatrixDSym correctionsCov(5);
        
        TVectorD residual(2);
        TVectorD measErr(2);
        TVectorD residualErr(2);
        TVectorD downWeight(2);

        TVectorD scat(2);
        scat.Zero();
        TVectorD scatPrecSensor(2);

        double p = _eBeam; // beam momentum 
 
        if( abs(_eBeam)<1e-12) streamlog_out(WARNING) << " provided beam energy is too low, _eBeam = " << _eBeam << " check inputs!" << std::endl;
//        double invP = 1./_eBeam;

        EVENT::TrackerHitVec::const_reverse_iterator itrHit;
        for ( itrHit = trackCandidate.rbegin(); itrHit != trackCandidate.rend(); ++itrHit ) 
        {
            const int planeID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itrHit) );

            const int hitGblLabel = _hitId2GblPointLabel[ (*itrHit)->id() ];
                  
            gblTraj->getResults( hitGblLabel, corrections, correctionsCov );
            
            streamlog_out(MESSAGE0) << "Corrections: " << std::endl;
            streamlog_message( MESSAGE0, corrections.Print();, std::endl; );
            
            gblTraj->getMeasResults( hitGblLabel, numData, residual, measErr, residualErr, downWeight);
            
            // retrieve original hit coordinates
            double pos[3] = { (*itrHit)->getPosition()[0], (*itrHit)->getPosition()[1], (*itrHit)->getPosition()[2] };
 
            streamlog_out(DEBUG1) << "meas point : " << hitGblLabel << std::endl;
	    streamlog_out(DEBUG1) << "p0 : " << pos[0] << " p1: " << pos[1] << " p2: " << pos[2] << std::endl;

	    double hitPointLocal[] = {pos[0], pos[1], pos[2]};
            double hitPointGlobal[] = {0.,0.,0.};
            geo::gGeometry().local2Master( planeID, hitPointLocal , hitPointGlobal);
            streamlog_out(MESSAGE1) << "hitg2= " << hitPointGlobal[0] << " " << hitPointGlobal[1] << " " << hitPointGlobal[2] << std::endl;
            streamlog_out(MESSAGE1) << "hitl2= " << hitPointLocal[0] << " " << hitPointLocal[1] << " " << hitPointLocal[2] << std::endl;

           
            // correct original values to the fitted ones
            pos[0] -= residual[0] ; // residual = meas - fitted -> to get fitted from measured;
            pos[1] -= residual[1] ;


	    double trackPointLocal[] = { pos[0], pos[1], pos[2] };
	    double trackPointGlobal[] = { 0., 0., 0. };
            geo::gGeometry().local2Master( planeID, trackPointLocal, trackPointGlobal );

//--------------- get the track slope...
            const double dz = trackPointGlobal[2] - prevZ;
//                step = dz; 
//                prevZ = hitPointGlobal[2];

//                TVectorD trackParamPrediction = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ, dz );
//                TVectorD trackParamPrediction = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], refPoint[2], dz );
                TVectorD trackParamPrediction = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ, dz );
                streamlog_message( DEBUG2,                trackParamPrediction.Print();, std::endl; );
                prevZ = trackPointGlobal[2];

                prevState[0] =   trackParamPrediction[0];
                prevState[1] =   trackParamPrediction[1];
                prevState[2] =   trackParamPrediction[2];
//---------------

	    double trackDirGlobal[] = { trackParamPrediction[2], trackParamPrediction[3], 1.};
            double trackDirLocal[] = { 0., 0., 0. };
//           double trackDirLocal[] = { corrections[1], corrections[2], 1. };

//	    geo::gGeometry().local2MasterVec( planeID, trackDirLocal, trackDirGlobal );
	    geo::gGeometry().master2LocalVec( planeID, trackDirGlobal, trackDirLocal );

            streamlog_out(MESSAGE1) << "fitg2= " << trackPointGlobal[0] << " " << trackPointGlobal[1] << " " << trackPointGlobal[2] << std::endl;
            streamlog_out(MESSAGE1) << "fitl2= " << trackPointLocal[0] << " " << trackPointLocal[1] << " " << trackPointLocal[2] << std::endl;
 
            streamlog_out(MESSAGE1) << "fitDIRg2= " << trackDirGlobal[0] << " " << trackDirGlobal[1] << " " << trackDirGlobal[2] << std::endl;
            streamlog_out(MESSAGE1) << "fitDIRl2= " << trackDirLocal[0] << " " << trackDirLocal[1] << " " << trackDirLocal[2] << std::endl;
   


                EVENT::FloatVec hitcov(4);
                hitcov[0]=0.01;
                hitcov[1]=0.00;
                hitcov[2]=0.01;
                hitcov[3]=0.00;

                // check Processor Parameters for plane resolution // Denys
                if( _paramterIdXResolutionVec.size() > 0 && _paramterIdYResolutionVec.size() > 0 )
                {
                  for(int izPlane=0;izPlane<_paramterIdPlaneVec.size();izPlane++)
                  {
                    if( _paramterIdPlaneVec[izPlane] == planeID )
                    {  
                      hitcov[0] =  _paramterIdXResolutionVec[izPlane];
                      hitcov[2] =  _paramterIdYResolutionVec[izPlane];
 
                      hitcov[0] *= hitcov[0]; // squared !
                      hitcov[2] *= hitcov[2]; // squared !
                      break;
                    } 
                  }
                }
                else
                {  
                  hitcov = (*itrHit)->getCovMatrix();            
                  
                }




            gbl::GblPoint point(jacPointToPoint);
            // Calculate projection matrix

 		// Calculate projection matrix
                TMatrixD proL2m(2, 2);

		const TGeoHMatrix* globalH = geo::gGeometry().getHMatrix( hitPointGlobal );
		const TGeoHMatrix& globalHInv = globalH->Inverse();
		const double* rotation = globalHInv.GetRotationMatrix();

		proL2m[0][0] = rotation[0]; // x projection, xx
		proL2m[0][1] = rotation[1]; // y projection, xy
		proL2m[1][0] = rotation[3]; // x projection, yx
		proL2m[1][1] = rotation[4]; // y projection, yy

// add measurment (residuals) in the measurement system (module 2D coordinates)
                addMeasurementsGBL( point, residual, measErr, hitPointLocal, trackPointLocal, hitcov, proL2m);

// add global derivatives derived from the track parameters after the track fit (coordinate system?)
                addGlobalParametersGBL( point, alDer, globalLabels, planeID, trackPointLocal, trackDirLocal[0], trackDirLocal[1] );
//              addGlobalParametersGBL( point, alDer, globalLabels, planeID, trackPointLocal, trackDirGlobal[0], trackDirGlobal[1] );

// add scatterrers
                addSiPlaneScattererGBL(point, scat, scatPrecSensor, planeID, p);
 
                if ( itrHit != ( trackCandidate.rend() -1) )
                {
                    // Go to global coordinates
                    const int nextPlaneID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*(itrHit + 1)) );
                    const double* nextHitPoint = (*(itrHit + 1))->getPosition();
                    double nextHitPointGlobal[] = { 0., 0., 0. };
                    geo::gGeometry().local2Master(nextPlaneID,nextHitPoint,nextHitPointGlobal);
                    double step = hitPointGlobal[2] - nextHitPointGlobal[2];
                    streamlog_out(MESSAGE1) << "nexg2= " << nextHitPointGlobal[0] << " " << nextHitPointGlobal[1] << " " << nextHitPointGlobal[2] << std::endl;

                    jacPointToPoint = PropagatePar( step, invP, corrections[3], corrections[4], corrections[1], corrections[2], hitPointGlobal[2] );
                }

                pushBackPointMille( pointList, point, (*itrHit)->id() );
         
        } // loop over track hits

        gbl::GblTrajectory* traj;
        traj = new gbl::GblTrajectory( pointList, false );
 
        if ( chi2 < _chi2cut ) 
        {
             traj->milleOut( *_mille );
        }
        
        delete traj; 
 
   }// Method end prepareMilleOut

} // namespace eutelescope

#endif
