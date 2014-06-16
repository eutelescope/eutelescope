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
    _parameterIdXShiftsMap(),
    _parameterIdYShiftsMap(),
    _parameterIdZShiftsMap(),
    _parameterIdXRotationsMap(),
    _parameterIdYRotationsMap(),
    _parameterIdZRotationsMap(),
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
    _parameterIdXShiftsMap(),
    _parameterIdYShiftsMap(),
    _parameterIdZShiftsMap(),
    _parameterIdXRotationsMap(),
    _parameterIdYRotationsMap(),
    _parameterIdZRotationsMap(),
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
      _parameterIdPlaneVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdXResolutionVec( const std::vector<float>& vector)
    {
      _parameterIdXResolutionVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdYResolutionVec( const std::vector<float>& vector)
    {
      _parameterIdYResolutionVec = vector;
    }
       

    void EUTelGBLFitter::setParamterIdXRotationsMap( const std::map<int, int>& map ) {
        _parameterIdXRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdYRotationsMap( const std::map<int, int>& map ) {
        _parameterIdYRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdZRotationsMap( const std::map<int, int>& map ) {
        _parameterIdZRotationsMap = map;
    }

    void EUTelGBLFitter::setParamterIdZShiftsMap( const std::map<int, int>& map ) {
        _parameterIdZShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdYShiftsMap( const std::map<int, int>& map ) {
        _parameterIdYShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdXShiftsMap( const std::map<int, int>& map ) {
        _parameterIdXShiftsMap = map;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdXRotationsMap() const {
        return _parameterIdXRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdYRotationsMap() const {
        return _parameterIdYRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdZRotationsMap() const {
        return _parameterIdZRotationsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdZShiftsMap() const {
        return _parameterIdZShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdYShiftsMap() const {
        return _parameterIdYShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdXShiftsMap() const {
        return _parameterIdXShiftsMap;
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
 
    void EUTelGBLFitter::SetTrackCandidates( const vector<IMPL::TrackImpl*>& trackCandidatesVec) {

        this->_trackCandidatesVec = trackCandidatesVec	;
        return;
    }
       
    void EUTelGBLFitter::SetTrackCandidates( const EVENT::TrackVec& trackCandidates) {

        this->_trackCandidates = trackCandidates;
        return;
    }


    void EUTelGBLFitter::Clear() {
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
     
        streamlog_out(DEBUG4) << " addMeasurementsGBL " << std::endl;
 
        meas[0] = hitpos[0] - predpos[0];
        meas[1] = hitpos[1] - predpos[1];
        measPrec[0] = 1. / hitcov[0];	// cov(x,x)
        measPrec[1] = 1. / hitcov[2];	// cov(y,y)

        streamlog_out(DEBUG4) << "Residuals:" << std::endl;
        streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] << std::endl;
        streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20) << measPrec[1] << std::endl;

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

        const double X0Si  = thicknessSi / radlenSi; // Si 
        const double X0Kap = thicknessKap / radlenKap; // Kapton                

        const double tetSi  = Utility::getThetaRMSHighland(p, X0Si);
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

        streamlog_out(MESSAGE0) << " addGlobalParametersGBL " << std::endl;

        alDer[0][0] = -1.0; // dx/dx
        alDer[0][1] =  0.0; // dx/dy
        alDer[1][0] =  0.0; // dy/dx
        alDer[1][1] = -1.0; // dy/dy
        globalLabels[0] = _parameterIdXShiftsMap[iPlane]; // dx
        globalLabels[1] = _parameterIdYShiftsMap[iPlane]; // dy


        if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][2] = -predpos[1]; // dx/rot
            alDer[1][2] =  predpos[0]; // dy/rot
            globalLabels[2] = _parameterIdZRotationsMap[iPlane]; // rot z
        }


        if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][3] =   -xSlope; // dx/dz
            alDer[1][3] =   -ySlope; // dy/dz
            globalLabels[3] = _parameterIdZShiftsMap[iPlane]; // dz
        }

        if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][4] =   predpos[0]*xSlope; // dx/rot y
            alDer[1][4] =   predpos[0]*ySlope; // dy/rot y
            globalLabels[4] = _parameterIdYRotationsMap[iPlane]; // drot y  - actually X?
            alDer[0][5] =  -predpos[1]*xSlope; // dx/rot x          
            alDer[1][5] =  -predpos[1]*ySlope; // dy/rot x         
            globalLabels[5] = _parameterIdXRotationsMap[iPlane]; // drot x  - actually Y?
        }


// partial alignment 
        if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
        }
        if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
            alDer[0][3] = predpos[1]*xSlope; // dx/rot x
            alDer[1][3] = predpos[1]*ySlope; // dy/rot x
            globalLabels[3] = _parameterIdXRotationsMap[iPlane]; // drot x
        }
 
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _parameterIdXRotationsMap[iPlane]; // drot x
        }
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _parameterIdXRotationsMap[iPlane]; // drot x
        }

        point.addGlobals(globalLabels, alDer);
 
        streamlog_out(MESSAGE0) << " addGlobalPArametersGBL over " << std::endl;
   }

    void EUTelGBLFitter::pushBackPoint( std::vector< gbl::GblPoint >& pointListTrack, const gbl::GblPoint& pointTrack, int hitid ) {
        pointListTrack.push_back(pointTrack);
       
        streamlog_out(MESSAGE0) << endl << "pushBackPoint: " << hitid << std::endl;
        // store point's GBL label for future reference
        _hitId2GblPointLabel[ hitid ] = static_cast<int>(pointListTrack.size());
    }

    void EUTelGBLFitter::pushBackPointMille( std::vector< gbl::GblPoint >& pointListMille, const gbl::GblPoint& pointMille, int hitid ) {
        pointListMille.push_back(pointMille);
        
        streamlog_out(MESSAGE0) << endl << "pushBackPointMille: " << hitid << std::endl;
        // store point's GBL label for future reference
        _hitId2GblPointLabelMille[ hitid ] = static_cast<int>(pointListMille.size());
    }

    /**
     * merging too sources of information: a copy of a Track Candidate (*itTrkCand) gets updated with gblTraj
     * 
     * @param fittrack pointer to track object to be stored
     * @param chi2     Chi2 of the track fit
     * @param ndf      NDF of the track fit
     */
     void EUTelGBLFitter::prepareLCIOTrack( gbl::GblTrajectory* gblTraj, const vector<IMPL::TrackImpl*>::const_iterator& itTrkCand, 
                                          double chi2, int ndf, double omega, double d0, double z0, double phi, double tanlam ) {
 
// output  track

        IMPL::TrackImpl * fittrack = new IMPL::TrackImpl( **itTrkCand); 
 
        unsigned int numData;
        TVectorD corrections(5);
        TMatrixDSym correctionsCov(5);
 
        TVectorD residual(2);
        TVectorD measErr(2);
        TVectorD residualErr(2);
        TVectorD downWeight(2);
 
          streamlog_out(MESSAGE1) << endl; 
 
          int nstates = fittrack->getTrackStates().size();
          streamlog_out(MESSAGE1) << "states " << nstates << "    " << endl;

          for(int i=0;i < nstates; i++) 
            {
                const IMPL::TrackStateImpl* const_trkState = static_cast <const IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i) ) ;

                double fitPointLocal[] = {0.,0.,0.};
                fitPointLocal[0] = const_trkState->getReferencePoint()[0] ;
                fitPointLocal[1] = const_trkState->getReferencePoint()[1] ;
                fitPointLocal[2] = const_trkState->getReferencePoint()[2] ;

                double bd0        = const_trkState->getD0() ;
 	        double bphi       = const_trkState->getPhi() ;
                double bomega     = const_trkState->getOmega() ;
	        double btanlambda = const_trkState->getTanLambda() ;
	        double bz0        = const_trkState->getZ0() ;

                const int hitGblLabel = _hitId2GblPointLabel[ const_trkState->id() ];

                streamlog_out(MESSAGE1) << hitGblLabel << " corr: " << const_trkState->id() 
                                        << " [d0]" << std::setw(6) << bd0 << ":" 
                                        << " [phi]" << std::setw(6) << bphi << ":" 
                                        << " [ome]" << std::setw(6) << bomega << ":" 
                                        << " [tanl]" << std::setw(6) << btanlambda << ":" 
                                        << " [z0]" << std::setw(6) << bz0 << ":" << " ["<< setw(3) << const_trkState->getLocation() <<"]" 
                                        << " points:"  << setw(8) << fitPointLocal[0] << setw(8) << fitPointLocal[1] << setw(8) << fitPointLocal[2] << " "  ;   

               streamlog_out(MESSAGE1) << endl; 
            }
 
          streamlog_out(MESSAGE1) << endl; 
  
        const EVENT::TrackerHitVec& ihits = (*itTrkCand)->getTrackerHits();
        int itrk = (*itTrkCand)->id();
        int nhits =  ihits.size( ) ;
        int expec = _parameterIdPlaneVec.size();
        streamlog_out(MESSAGE1) <<  " track itrk:" <<  itrk  << " with " << nhits << " at least " << expec << "       " ;//std::endl;
        for (int i = 0; i< ihits.size(); i++ ) 
        { 
            EVENT::TrackerHit* ihit = ihits[i];
            int ic = ihit->id();
            streamlog_out(MESSAGE0) <<  ic << " ";
        }
        streamlog_out(MESSAGE1) << std::endl;
 
          for(int i=0;i < nstates; i++) 
            {
                const IMPL::TrackStateImpl* const_trkState = static_cast <const IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i) ) ;
                      IMPL::TrackStateImpl* trkState = static_cast <      IMPL::TrackStateImpl*> ( fittrack->getTrackStates().at(i) ) ;


                float fitPointLocal[] = {0.,0.,0.};
                fitPointLocal[0] = const_trkState->getReferencePoint()[0] ;
                fitPointLocal[1] = const_trkState->getReferencePoint()[1] ;
                fitPointLocal[2] = const_trkState->getReferencePoint()[2] ;

                double bd0        = const_trkState->getD0() ;
 	        double bphi       = const_trkState->getPhi() ;
                double bomega     = const_trkState->getOmega() ;
	        double btanlambda = const_trkState->getTanLambda() ;
	        double bz0        = const_trkState->getZ0() ;

                trkState->setD0(        const_trkState->getD0() + corrections[0] ) ;
 	        trkState->setPhi(       const_trkState->getPhi() + corrections[1] ) ;
                trkState->setOmega(     const_trkState->getOmega() + corrections[2] ) ;
	        trkState->setTanLambda( const_trkState->getTanLambda() + corrections[3] ) ;
	        trkState->setZ0(        const_trkState->getZ0() + corrections[4] ) ;

                double ed0        = trkState->getD0() ;
 	        double ephi       = trkState->getPhi()  ;
                double eomega     = trkState->getOmega()  ;
	        double etanlambda = trkState->getTanLambda()  ;
	        double ez0        = trkState->getZ0()  ;


                const int hitGblLabel = _hitId2GblPointLabel[ const_trkState->id() ];

                _hitId2GblPointLabel.insert( std::make_pair(   trkState->id() , hitGblLabel  ) );

                gblTraj->getResults( hitGblLabel, corrections, correctionsCov );

                streamlog_out(MESSAGE1) << hitGblLabel << " corr: " << trkState->id() 
                                        << " [d0]" << std::setw(6)  << ed0 <<  ":"
                                        << " [phi]" << std::setw(6)  << ephi<<  ":"
                                        << " [ome]" << std::setw(6)  << eomega<<  ":"
                                        << " [tanl]" << std::setw(6)  << etanlambda<<  ":"
                                        << " [z0]" << std::setw(6)  << ez0 <<  ":" << " ["<< setw(3) << trkState->getLocation() <<"]" ;      

                if( trkState->getLocation() >= 0 ) {
                  gblTraj->getMeasResults( hitGblLabel, numData, residual, measErr, residualErr, downWeight);
                  // correct original values to the fitted ones

                  fitPointLocal[0] += residual[0];
                  fitPointLocal[1] += residual[1];
                  trkState->setReferencePoint(fitPointLocal);

                  double PointLocal[] = {0.,0.,0.};
                  PointLocal[0] = trkState->getReferencePoint()[0] ;
                  PointLocal[1] = trkState->getReferencePoint()[1] ;
                  PointLocal[2] = trkState->getReferencePoint()[2] ;

                  streamlog_out(MESSAGE1) << " points:"  << setw(8) << PointLocal[0] << setw(8)  << PointLocal[1] << setw(8) << PointLocal[2] << " "  ;   
                  streamlog_out(MESSAGE1) << "  re[0]:"  << std::setw(7) << residual[0] << " re[1]:"  << std::setw(7) << residual[1] ;
                }

                streamlog_out(MESSAGE1) << endl;      

            }
            streamlog_out(MESSAGE1) << std::endl;


        const EVENT::TrackerHitVec& chits = fittrack->getTrackerHits();
         itrk = fittrack->id();
         nhits =  chits.size( ) ;
         expec =  _parameterIdPlaneVec.size();
		streamlog_out(MESSAGE1) <<  " track itrk:" <<  itrk  << " with " << nhits << " at least " << expec  << "       ";//std::endl;
        for (int i = 0; i< chits.size(); i++ ) 
        { 
            EVENT::TrackerHit* ihit = chits[i];
            int ic = ihit->id();
            streamlog_out(MESSAGE0) <<  ic << " ";
        }
        streamlog_out(MESSAGE1) << std::endl;
 
         // prepare track
        fittrack->setChi2 ( chi2 );      // Chi2 of the fit (including penalties)
        fittrack->setNdf  ( ndf );        // Number of planes fired (!)
        
        // add track to LCIO collection vector
        _fittrackvec->addElement( fittrack );
 
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
            
            streamlog_out(DEBUG2) << "LCIO - Corrections: " << std::endl;
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
                    const int nextPlaneID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itrHit+1) );
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
   
    void EUTelGBLFitter::FitSingleTrackCandidate(EVENT::TrackVec::const_iterator& itTrkCand)
    {
            // relies on sane itTrkCand -> sanity checks ?

	    // Z axis points along beam direction.
	    double pt = ( 1./(*itTrkCand)->getOmega() ) * _beamQ;
	    double px = (*itTrkCand)->getTanLambda() * pt;
	    double py = pt * sin( (*itTrkCand)->getPhi() );
	    double pz = pt * cos( (*itTrkCand)->getPhi() );

	    double tx   = px / pz;
	    double ty   = py / pz;
	    double invP = _beamQ / sqrt( px*px + pt*pt );
 
            streamlog_out(DEBUG4) << "::FitSingleTrackCandidate FitTracks   ";
            streamlog_out(DEBUG4) << " omega: " << std::setw(8) << (*itTrkCand)->getOmega()  ;
            streamlog_out(DEBUG4) << " Q: " << std::setw(8) << _beamQ ;
            streamlog_out(DEBUG4) << " px: " << std::setw(8) << px ;
            streamlog_out(DEBUG4) << " py: " << std::setw(8) << py ;
            streamlog_out(DEBUG4) << " pz: " << std::setw(8) << pt << std::endl;
 
	    // Reference point is the last hit on a track
	    const float *refPoint = (*itTrkCand)->getReferencePoint();

            // remember last hit-point from the track candidate below
            double start[3] = { 0.,0.,-500.}; // 
 
            // Loop over hits on a track candidate
            const EVENT::TrackerHitVec& hits = (*itTrkCand)->getTrackerHits();
            EVENT::TrackerHitVec::const_reverse_iterator itHit;

      const map< int, int > sensorMap = geo::gGeometry().sensorZOrdertoIDs();
      int planeID     = sensorMap.at(0); // the first first plane in the array of the planes according to z direction. // assume not tilted plane. 
      const int    iPlane             = geo::gGeometry().sensorIDtoZOrder(planeID);
      const double thicknessSen       = geo::gGeometry()._siPlanesLayerLayout->getSensitiveThickness(iPlane );
      const double thicknessLay       = geo::gGeometry()._siPlanesLayerLayout->getLayerThickness(iPlane );
      start[2]    = geo::gGeometry().siPlaneZPosition(planeID) - thicknessSen - thicknessLay; //initial z position to the most-first plane
      double dir[3]   = {0.,0.,1.};  // as good as any other value along z axis.
      float dpoint[3] = {0.,0.,0.};  // initialise output point-vector (World frame)
      float lpoint[3] = {0.,0.,0.};  // initialise output point-vector (sensor frame)

      streamlog_out ( DEBUG0 ) << "::FitSingleTrackCandidate starting point: "  <<  start[0] << " "  << start[1] << " " << start[2] << " dir:"<< dir[0]<<" "<<dir[1]<<" "<<dir[2] << endl;
 
      map< int, int >::const_iterator iter = sensorMap.begin();
      while ( iter != sensorMap.end() ) {
        bool found = false;
        if( iter->first > -999 )
          {
            int getErrorID = geo::gGeometry( ).findNextPlaneEntrance(  start,  dir, iter->second, dpoint ) ;

            // get to local frame for sensor newSensorID:
//            geo::gGeometry().master2Local(dpoint, lpoint);


            if( getErrorID < 0 ) 
            {
              streamlog_out ( DEBUG4 ) << "no Entrance: " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] << " err:"<< getErrorID << endl;
            } else {     
              streamlog_out ( DEBUG0 ) << "identified NextPlane Entrance: " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] <<  " err:"<< getErrorID << endl;
              start[0] = dpoint[0];
              start[1] = dpoint[1];
              start[2] = dpoint[2];
              for ( itHit = hits.rbegin(); itHit != hits.rend(); ++itHit) {

                 const int planeID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itHit) );

                 if(planeID == getErrorID)
                 {
                   streamlog_out ( DEBUG4 ) << "match sensor found: extrapl@ " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] <<  " @"<< getErrorID << endl;
                  //
                  const double* hitPointLocal = (*itHit)->getPosition();
                  double hitPointGlobal[] = {0.,0.,0.};
                  geo::gGeometry().local2Master(planeID,hitPointLocal,hitPointGlobal);
   		  streamlog_out(DEBUG4) << "hitm= " << hitPointGlobal[0] << " " << hitPointGlobal[1] << " " << hitPointGlobal[2] << std::endl;
                  found = true;
                 }
              }
              if(!found) 
              {
                streamlog_out ( DEBUG4 ) << "identified NextPlane UNMATCHED " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] <<  " err:"<< getErrorID << endl;
              } 
            }
          }
        ++iter;
      } 
    }
 
    void EUTelGBLFitter::CalculateProjMatrix(TMatrixD& proL2m, double* hitPointGlobal )
    {  
		// Calculate projection matrix

		const TGeoHMatrix* globalH = geo::gGeometry().getHMatrix( hitPointGlobal );
		const TGeoHMatrix& globalHInv = globalH->Inverse();
		const double* rotation = globalHInv.GetRotationMatrix();

		proL2m[0][0] = rotation[0]; // x projection, xx
		proL2m[0][1] = rotation[1]; // y projection, xy
		proL2m[1][0] = rotation[3]; // x projection, yx
		proL2m[1][1] = rotation[4]; // y projection, yy

    }

// GBL Trajectory treatment ::  Fit and dump into LCIO
    void EUTelGBLFitter::PerformFitGBLTrajectory( gbl::GblTrajectory* traj, vector<IMPL::TrackImpl*>::const_iterator& itTrkCand, double invP  ) {
                
            streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::PerformFitGBLTrajectory -- starting " << endl;

            double loss = 0.;
            double chi2 = 0.;
            int ndf = 0;
            // perform GBL fit

                int ierr = 0;

		        if ( streamlog_level(MESSAGE0) ){
	        	  std::cout << "pre fit FitTrack - trajectory: " << std::endl;
//		 	  traj->printTrajectory(1);
//			  traj->printData();
//	 		  traj->printPoints(1);
	        	  std::cout << "pre fit FitTrack - trajectory: end;" << std::endl;

		        }
	
                if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( chi2, ndf, loss, _mEstimatorType );
                else ierr = traj->fit( chi2, ndf, loss );

                streamlog_out(MESSAGE0) << "ierr : "<< ierr << " and chi2: " << chi2 << std::endl;

                if ( ierr  )
                    {
		        if ( streamlog_level(MESSAGE0) ){
	        	  std::cout << "after fit FitTrack - trajectory: " << std::endl;
//		 	  traj->printTrajectory(1);
//			  traj->printData();
//	 		  traj->printPoints(1);
	        	  std::cout << "after fit FitTrack - trajectory: " << std::endl;
  		        }
		    }
  
                    // for some reason (??) need to keep the same numbering for trajectories as for the Track Candidates
                    vector<IMPL::TrackImpl*>::const_iterator begin = _trackCandidatesVec.begin();
                    _gblTrackCandidates.insert( std::make_pair( std::distance( begin, itTrkCand ), traj ) );
                
                    // Write fit result
                    prepareLCIOTrack( traj, itTrkCand, chi2, ndf, invP, 0., 0., 0., 0. );

                    // Prepare and write Mille Out only when enabled (and chi2 is below _chi2cut)
                    if (_alignmentMode != Utility::noAlignment) 
                    {
                       prepareMilleOut( traj, itTrkCand) ;//, chi2, ndf, invP, 0., 0., 0., 0. );
                    }
 
           streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::PerformFitGBLTrajectory -- finished " << endl;
           
    }

// convert input TrackCandidates and TrackStates into a GBL Trajectory
void EUTelGBLFitter::FillInformationToGBLPointObject(IMPL::TrackImpl* trackimpl){
	// sanity check. Mustn't happen in principle.
  if ( trackimpl->getTrackerHits().size() > geo::gGeometry().nPlanes() ){
  	streamlog_out(ERROR) << "Sanity check. This should not happen in principle. Number of hits is greater then number of planes" << std::endl;
   	return;
  }



}


    void EUTelGBLFitter::TrackCandidatesToGBLTrajectory( vector<IMPL::TrackImpl*>::const_iterator& itTrkCand) {

            // sanity check. Mustn't happen in principle.
            if ( (*itTrkCand)->getTrackerHits().size() > geo::gGeometry().nPlanes() )
            {
              streamlog_out(ERROR) << "Sanity check. This should not happen in principle. Number of hits is greater then number of planes" << std::endl;
              return;
            }



            // Z axis points along beam direction.
            double pt = ( 1./(*itTrkCand)->getOmega() ) * _beamQ;
            double px = (*itTrkCand)->getTanLambda() * pt;
   	    double py = pt * sin( (*itTrkCand)->getPhi() );
	    double pz = pt * cos( (*itTrkCand)->getPhi() );

	    double tx   = px / pz;
	    double ty   = py / pz;
//	    double invP = _beamQ / sqrt( px*px + pt*pt );
 
            double p = _eBeam; // beam momentum
            double invP = 1./p;

            //GBL trajectory construction
            // lsit of GBL points that goes eventually into the trajectory
            std::vector< gbl::GblPoint > pointList;

            // introduce new jacobian from point to point // not sure why we need it here ? 
            //
            TMatrixD jacPointToPoint(5, 5);
            jacPointToPoint.UnitMatrix();

            // get measurement points:
            const EVENT::TrackerHitVec& hits = (*itTrkCand)->getTrackerHits();
            EVENT::TrackerHitVec::const_reverse_iterator itHit;
            int imatch=0;
            streamlog_out( MESSAGE0 ) << "list requested planes: " ;
            for(int izPlane=0; izPlane<_parameterIdPlaneVec.size(); izPlane++) {
              int kPlaneID =  _parameterIdPlaneVec[izPlane];
              streamlog_out( MESSAGE0 ) << " [" << imatch << ":" << kPlaneID ;
              bool ifound = false;
              for ( itHit = hits.rbegin(); itHit != hits.rend(); ++itHit) {
              const int planeID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itHit) );
               if( kPlaneID == planeID )
                {  
                  imatch++;
                  streamlog_out( MESSAGE0 ) << " yes "   ;
                  ifound = true;
                  break;
                }
              }     
              if( !ifound)  streamlog_out( MESSAGE0 ) << " not "   ;      
              streamlog_out( MESSAGE0 ) << "] "   ;
            }
            streamlog_out( MESSAGE0 ) << std::endl;

            streamlog_out( MESSAGE1 ) << "from the list of requested planes and available hit on a track candidate: " << imatch << "/" <<_parameterIdPlaneVec.size() << " found" << std::endl;

            if( imatch != _parameterIdPlaneVec.size() ) 
            { 
              streamlog_out(MESSAGE1) << " Number of hits does not correspond to the number of planes selected for tracking" << std::endl;
              return;
            }



            TVectorD meas(2);
            TVectorD measPrec(2); // precision = 1/resolution^2

            TVectorD scat(2);
            scat.Zero();

            TVectorD scatPrecSensor(2);
            TVectorD scatPrec(2);
 

            // loop through scattering planes (including measurement planes) ::
            int nstates = (*itTrkCand)->getTrackStates().size();
            for(int i=0;i < nstates; i++) 
            {

                streamlog_out(MESSAGE1) << "state: at " << i << " of " << nstates ;
                IMPL::TrackStateImpl* trk = static_cast < IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i) ) ;
                int trkVolumeID =  trk->getLocation();

                if ( trkVolumeID < 0 )
                {
                   streamlog_out(DEBUG0) << "Sanity check. SensorID can not be negative. Negative TrackStates are kept for the beginning and end of the TrackStates on Pattern Recognition. skip this one.";
                   streamlog_out(MESSAGE1) << std::endl; // to be consistent with the MESSAGE1 level printouts in the for-loop
                   continue;
                } 

                streamlog_out(MESSAGE1) << "  [ " << trk->id() << " ] " ;
 
                double fitPointLocal[] = {0.,0.,0.};
                fitPointLocal [0] = trk->getReferencePoint()[0] ;
                fitPointLocal [1] = trk->getReferencePoint()[1] ;
                fitPointLocal [2] = trk->getReferencePoint()[2] ;

                double fitPointGlobal[] = {0.,0.,0.};
                geo::gGeometry().local2Master( trkVolumeID, fitPointLocal, fitPointGlobal);

                // fill next GBL point:
                gbl::GblPoint point(jacPointToPoint);


                if(  trkVolumeID == -1 )
                {
                   // begining of the track - does not exist normally
                }
                else if(  trkVolumeID == -2 )
                {
                   // end of the track - can ignore
                   streamlog_out(MESSAGE1) << std::endl;
                }
                if(  trkVolumeID >= 0 )
                {


                 // look for the measurement point corresponding to the TrackState Location (fitted hit on the volume surface)
                 EVENT::TrackerHitVec::const_iterator itHit;
                 int imatch=0;
                 for ( itHit = hits.begin(); itHit != hits.end(); itHit++) {
                  const int planeID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(*itHit) );
                  if( trkVolumeID == planeID )
                  { 
                    EVENT::FloatVec hitcov(4);
                    hitcov[0]=0.01;
                    hitcov[1]=0.00;
                    hitcov[2]=0.01;
                    hitcov[3]=0.00;

                    // check Processor Parameters for plane resolution // Denys
                    if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 )
                    {
                      for(int izPlane=0;izPlane<_parameterIdPlaneVec.size();izPlane++)
                      {
                        if( _parameterIdPlaneVec[izPlane] == planeID )
                        {  
                          hitcov[0] =  _parameterIdXResolutionVec[izPlane];
                          hitcov[2] =  _parameterIdYResolutionVec[izPlane];

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

 
                    const double* hitPointLocal = (*itHit)->getPosition();
                    double hitPointGlobal[] = {0.,0.,0.};
                    geo::gGeometry().local2Master(planeID,hitPointLocal,hitPointGlobal);
/*
                     streamlog_out(MESSAGE0) << " | " <<  setw(3)  << right 
                                               <<   planeID  << " " << fixed << setw(11) << setprecision (3) << noshowpos << setw(7) << " " <<
                                                                      hitPointGlobal[0] << " : " <<  fitPointGlobal[0] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointGlobal[1] << " : " <<  fitPointGlobal[1] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointGlobal[2] << " : " <<  fitPointGlobal[2] ;
*/
                     streamlog_out(MESSAGE1) << " | " <<  setw(3)  << right 
                                               <<   planeID  << " " << fixed << setw(11) << setprecision (3) << noshowpos << setw(7) << " " <<
                                                                      hitPointLocal[0] << " : " <<  fitPointLocal[0] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointLocal[1] << " : " <<  fitPointLocal[1] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointLocal[2] << " : " <<  fitPointLocal[2] << " glo " << setw(8) << fitPointGlobal[2] ;

                 
                     // Calculate projection matrix
                     // GBL language "Local" -> our language "Telescope"="global" system
                     // GLB language "measurement" -> our language "measurement"="detector" system
                     TMatrixD proL2m(2, 2);
                     CalculateProjMatrix(proL2m, hitPointGlobal);
                
                     bool excludeFromFit = false;
                     if ( std::find( _excludeFromFit.begin(), _excludeFromFit.end(), planeID ) != _excludeFromFit.end() ) excludeFromFit = true;
 
                     if ( !excludeFromFit )
                     {
                       addMeasurementsGBL(point, meas, measPrec, hitPointLocal, fitPointLocal, hitcov, proL2m);
                     }

                  }                          
                 }

                 addSiPlaneScattererGBL(point, scat, scatPrecSensor, trkVolumeID, p);
                 pushBackPoint( pointList, point, trk->id() );
 
// get jacobian to arrive to next point (if not at last point already):: 
                 if( i < nstates-1 ){
                   // get next point
                   IMPL::TrackStateImpl* trk = static_cast < IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i+1) ) ;
                   const float *fitPoint  = trk->getReferencePoint();
                   const int nextVolumeID = trk->getLocation();
                   double hitPointLocal[] = {static_cast<double>(fitPoint[0]), static_cast<double>(fitPoint[1]), static_cast<double>(fitPoint[2]) };
                   double hitPointGlobal[] = {0.,0.,0.};
                   geo::gGeometry().local2Master( nextVolumeID, hitPointLocal, hitPointGlobal);

                   double ds = hitPointGlobal[2] - fitPointGlobal[2]; // actually it should be a method which returns Jacobian matrix (from state1 to state2)
                   jacPointToPoint = propagatePar(ds);
                   streamlog_out(MESSAGE1) << " step[" << i << "] = " << ds << std::endl;
                 } else {
                   streamlog_out(MESSAGE1) << std::endl;
                 }

                }

           }

            // check magnetic field at (0,0,0)
            const gear::BField&   B = geo::gGeometry().getMagneticFiled();
            const double Bmag       = B.at( TVector3(0.,0.,0.) ).r2();
                
            gbl::GblTrajectory* traj = 0;
            if ( Bmag < 1.E-6 ) {
                traj = new gbl::GblTrajectory( pointList, false );
            } else {
                traj = new gbl::GblTrajectory( pointList, true );
            }
 
            PerformFitGBLTrajectory( traj, itTrkCand, invP  );

    }

// convert input TrackCandidates and TrackStates into a GBL Trajectory
    void EUTelGBLFitter::TrackCandidatesToGBLTrajectories( ) {
  	//      Clear(); //  This should be done explictly outside the class. So it is not hidden away.

        ////////////////////////////////////////////////////////////////////////// prepare output collection!!!!!!!!!!NOT HERE HOW STUPID! Unless I am missing something
     //   try {
       //     _fittrackvec = new IMPL::LCCollectionVec( EVENT::LCIO::TRACK );
         //   _fithitsvec = new IMPL::LCCollectionVec( EVENT::LCIO::TRACKERHIT );
      //      IMPL::LCFlagImpl flag( _fittrackvec->getFlag( ) );
       //     flag.setBit( lcio::LCIO::TRBIT_HITS );
       //     _fittrackvec->setFlag( flag.getFlag( ) );
      //  } catch ( ... ) {
       //     streamlog_out( ERROR2 ) << "Can't allocate output collection" << std::endl;
      //  }
				////////////////////////////////////////////////////////////////////////////////////////////////////////
 
//
        vector<IMPL::TrackImpl*>::const_iterator itTrkCand;
        int trackcounter = 0;
        for ( itTrkCand = _trackCandidatesVec.begin(); itTrkCand != _trackCandidatesVec.end(); ++itTrkCand) {
            streamlog_out(MESSAGE) << " track candidate # " << trackcounter << " of " << _trackCandidatesVec.size() << std::endl; 

           // take care of a track candidate
           TrackCandidatesToGBLTrajectory( itTrkCand  );         
           trackcounter++;
        }
 
    }

    void EUTelGBLFitter::prepareMilleOut( gbl::GblTrajectory* gblTraj, const vector<IMPL::TrackImpl*>::const_iterator& itTrkCand) {



        IMPL::TrackImpl * fittrack =  *itTrkCand;

        unsigned int numData;
        TVectorD corrections(5);
        TMatrixDSym correctionsCov(5);
 
        TVectorD residual(2);
        TVectorD measErr(2);
        TVectorD residualErr(2);
        TVectorD downWeight(2);
        
          streamlog_out(MESSAGE1) << endl; 
 
          int nstates = fittrack->getTrackStates().size();
          streamlog_out(MESSAGE1) << "states " << nstates << "    " << endl;

// now get starting point:  // needed by GBL point
          TMatrixD jacPointToPoint(5, 5);
          jacPointToPoint.UnitMatrix();

          // Z axis points along beam direction.
          double pt = ( 1./(*itTrkCand)->getOmega() ) * _beamQ;
          double px = (*itTrkCand)->getTanLambda() * pt;
	  double py = pt * sin( (*itTrkCand)->getPhi() );
	  double pz = pt * cos( (*itTrkCand)->getPhi() );

	  double tx   = px / pz;
	  double ty   = py / pz;
	  double invP = _beamQ / sqrt( px*px + pt*pt );
 
          std::vector< gbl::GblPoint > pointList;


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
 
 
          for(int i=0;i < nstates; i++) 
            {
                const IMPL::TrackStateImpl* const_trkState = static_cast <const IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i) ) ;

                double fitPointLocal[] = {0.,0.,0.};
                fitPointLocal[0] = const_trkState->getReferencePoint()[0] ;
                fitPointLocal[1] = const_trkState->getReferencePoint()[1] ;
                fitPointLocal[2] = const_trkState->getReferencePoint()[2] ;

                double bd0        = const_trkState->getD0() ;
 	        double bphi       = const_trkState->getPhi() ;
                double bomega     = const_trkState->getOmega() ;
	        double btanlambda = const_trkState->getTanLambda() ;
	        double bz0        = const_trkState->getZ0() ;

                const int hitGblLabel = _hitId2GblPointLabel[ const_trkState->id() ];

                streamlog_out(MESSAGE1) << hitGblLabel << " corr: " << const_trkState->id() 
                                        << " [d0]" << std::setw(6) << bd0 << ":" 
                                        << " [phi]" << std::setw(6) << bphi << ":" 
                                        << " [ome]" << std::setw(6) << bomega << ":" 
                                        << " [tanl]" << std::setw(6) << btanlambda << ":" 
                                        << " [z0]" << std::setw(6) << bz0 << ":" << " ["<< setw(3) << const_trkState->getLocation() <<"]" 
                                        << " points:"  << setw(8) << fitPointLocal[0] << setw(8) << fitPointLocal[1] << setw(8) << fitPointLocal[2] << " "  ;   

               streamlog_out(MESSAGE1) << endl; 
            }
 
        streamlog_out(MESSAGE1) << endl; 

        TVectorD scat(2);
        scat.Zero();
        TVectorD scatPrecSensor(2);

        double p = _eBeam; // beam momentum 
 
        if( abs(_eBeam)<1e-12) streamlog_out(WARNING) << " provided beam energy is too low, _eBeam = " << _eBeam << " check inputs!" << std::endl;


        gbl::GblTrajectory* traj;
        traj = new gbl::GblTrajectory( pointList, false );
 
        if ( streamlog_level(MESSAGE0) ){
          std::cout << "MilleOut - trajectory: " << std::endl;
 	  traj->printTrajectory(1);
	  traj->printData();
 	  traj->printPoints(1);
        }

        traj->milleOut( *_mille );
        
        delete traj; 

    }

    void EUTelGBLFitter::prepareMilleOut( gbl::GblTrajectory* gblTraj, const EVENT::TrackVec::const_iterator& itTrkCand ) {
//                                          double chi2, int ndf, double omega, double d0, double z0, double phi, double tanlam ) {

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
 
          const float *refPoint = (*itTrkCand)->getReferencePoint(); // will fail miserably now with state in LOCAL frame

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
            const int planeID = Utility::GuessSensorID(static_cast< IMPL::TrackerHitImpl* >(*itrHit));

            const int hitGblLabel = _hitId2GblPointLabel[ (*itrHit)->id() ];
                  
            gblTraj->getResults( hitGblLabel, corrections, correctionsCov );
            
            streamlog_out(MESSAGE0) << "MilleOut - Corrections: " << std::endl;
            streamlog_message( MESSAGE0, corrections.Print();, std::endl; );
            
            gblTraj->getMeasResults( hitGblLabel, numData, residual, measErr, residualErr, downWeight);
            
            // retrieve original hit coordinates
            double pos[3] = { (*itrHit)->getPosition()[0], (*itrHit)->getPosition()[1], (*itrHit)->getPosition()[2] };
 
            streamlog_out(DEBUG1) << "meas point : " << hitGblLabel << std::endl;
	    streamlog_out(DEBUG1) << "p0 : " << pos[0] << " p1: " << pos[1] << " p2: " << pos[2] << std::endl;

	    double hitPointLocal[]  = {pos[0], pos[1], pos[2]};
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

            TVectorD trackParamPrediction = getXYZfromDzNum( invP, prevState[2], prevState[3], prevState[0], prevState[1], prevZ, dz );
            streamlog_message( DEBUG2,                trackParamPrediction.Print();, std::endl; );
            prevZ = trackPointGlobal[2];

            prevState[0] =   trackParamPrediction[0];
            prevState[1] =   trackParamPrediction[1];
            prevState[2] =   trackParamPrediction[2];
//---------------

	    double trackDirGlobal[] = { trackParamPrediction[2], trackParamPrediction[3], 1.};
            double trackDirLocal[] = { 0., 0., 0. };

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
                if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 )
                {
                  for(int izPlane=0;izPlane<_parameterIdPlaneVec.size();izPlane++)
                  {
                    if( _parameterIdPlaneVec[izPlane] == planeID )
                    {  
                      hitcov[0] =  _parameterIdXResolutionVec[izPlane];
                      hitcov[2] =  _parameterIdYResolutionVec[izPlane];
 
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
                TMatrixD proL2m(2, 2);
                CalculateProjMatrix(proL2m, hitPointGlobal);
 
                bool excludeFromFit = false;
                if ( std::find( _excludeFromFit.begin(), _excludeFromFit.end(), planeID ) != _excludeFromFit.end() ) excludeFromFit = true;
 
                if ( !excludeFromFit )
                {
// add measurment (residuals) in the measurement system (module 2D coordinates)
                    addMeasurementsGBL( point, residual, measErr, hitPointLocal, trackPointLocal, hitcov, proL2m);
                }
// add scatterrers
                addSiPlaneScattererGBL(point, scat, scatPrecSensor, planeID, p);

// add global derivatives derived from the track parameters after the track fit (coordinate system?) 
// this one needed only for alignment with millepede
                addGlobalParametersGBL( point, alDer, globalLabels, planeID, trackPointLocal, trackDirLocal[0], trackDirLocal[1] );
 
                if ( itrHit != ( trackCandidate.rend() -1) )
                {
                    // Go to global coordinates
   		    const int nextPlaneID = Utility::GuessSensorID(static_cast< IMPL::TrackerHitImpl* >(*(itrHit+1)) );
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
 
        if ( streamlog_level(MESSAGE0) ){
          std::cout << "MilleOut - trajectory: " << std::endl;
 	  traj->printTrajectory(1);
	  traj->printData();
 	  traj->printPoints(1);
        }

        traj->milleOut( *_mille );
        
        delete traj; 
 
   }// Method end prepareMilleOut

} // namespace eutelescope

#endif
