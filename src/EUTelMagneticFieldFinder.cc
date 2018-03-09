/* 
 * File:   EUTelMagneticFieldFinder.cc
 *
 * Created on July 2, 2013, 12:59 PM
 */

#include "EUTelMagneticFieldFinder.h"
#include <UTIL/CellIDEncoder.h>
#include "streamlog/streamlog.h"

#include "EUTelGeometryTelescopeGeoDescription.h"

#include "EUTelUtilityRungeKutta.h"

#include "gear/gearimpl/Vector3D.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif

//LCIO
#include "IMPL/TrackImpl.h"

#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <math.h>

using namespace eutelescope;

    namespace {

        struct hasSensorID
        {
            int _requiredSensorID;

            bool operator( )( MeasurementLayer* layer ) {
                return ( layer->sensorID( ) == _requiredSensorID );
            }
        } ;
    }
    
    EUTelKalmanFilter::EUTelKalmanFilter() : EUTelTrackFitter( "KalmanFilter" ), 
            _tracksCartesian(), 
            _trackStates(), 
            _allHits(),
            _allMeasurements(),
            _isHitsOK(false),
            _isReady(false),
            _allowedMissingHits(0),
            _maxTrackCandidates(0),
            _beamE(-1.),
            _beamQ(-1.),
            _beamEnergyUncertainty(0.),
            _beamAngularSpread(2,-1.),
	    _jacobianF(5,5),
            _trkParamCovCkkm1(5,5),
            _processNoiseQ(5,5),
            _gainK(5,2),
            _residualCovR(2,2),
            _eomIntegrator( new EUTelUtilityRungeKutta() ),
            _jacobianIntegrator( new EUTelUtilityRungeKutta() ),
            _eomODE( 0 ),
            _jacobianODE( 0 ){
                // Initialise ODE integrators for eom and jacobian       
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }

                {
                    _jacobianIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
            }

    EUTelKalmanFilter::EUTelKalmanFilter( std::string name ) : EUTelTrackFitter( name ),
            _tracksCartesian(),
            _trackStates(), 
            _allHits(),
            _allMeasurements(),
            _isHitsOK(false),
            _isReady(false),
            _allowedMissingHits(0),
            _maxTrackCandidates(0),
            _beamE(-1.),
            _beamQ(-1.),
            _beamEnergyUncertainty(0.),
            _beamAngularSpread(2,-1.),
	    _jacobianF(5,5),
            _trkParamCovCkkm1(5,5),
            _processNoiseQ(5,5),
            _gainK(5,2),
            _residualCovR(2,2),
            _eomIntegrator( new EUTelUtilityRungeKutta() ),
            _jacobianIntegrator( new EUTelUtilityRungeKutta() ),
            _eomODE( 0 ),
            _jacobianODE( 0 ){
                // Initialise ODE integrators for eom and jacobian       
                {
                    _eomODE = new eom::EOMODE(5);
                    // _eomIntegrator integrator becomes the owner of _eomODE and ButcherTableauDormandPrince
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }

                {
                    _jacobianIntegrator->setRhs( _jacobianODE );
                    _jacobianIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
            }

    EUTelKalmanFilter::~EUTelKalmanFilter() { 
        delete _eomIntegrator;
        delete _jacobianIntegrator;
    }

    

    
    

    /** Propagate track state by dz 
     * 
     * @param ts track state
     * @param dz propagation distance
     */
    void EUTelKalmanFilter::propagateTrackRefPoint( EUTelTrackStateImpl* ts, double dz ) {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::propagateTrackRefPoint()" << std::endl;
          // The formulas below are derived from equations of motion of the particle in 
          // magnetic field under assumption |dz| small. Must be valid for |dz| < 10 cm

	  // Get track parameters
	  const double invP = ts->getInvP();
	  const double x0 = ts->getX();
	  const double y0 = ts->getY();
	  //const double x0 = ts->getReferencePoint()[0];
	  //const double y0 = ts->getReferencePoint()[1];
	  const double z0 = ts->getReferencePoint()[2];
	  const double tx0 = ts->getTx();
	  const double ty0 = ts->getTy();

          // Get magnetic field vector
          gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field
	  const gear::BField&   B = geo::gGeometry().getMagneticField();
	  const double Bx         = B.at( vectorGlobal ).x();
	  const double By         = B.at( vectorGlobal ).y();
	  const double Bz         = B.at( vectorGlobal ).z();
	  const double mm = 1000.;
	  const double k = 0.299792458/mm;
          
	  const double sqrtFactor = sqrt( 1. + tx0*tx0 + ty0*ty0 );

	  const double Ax = sqrtFactor * (  ty0 * ( tx0 * Bx + Bz ) - ( 1. + tx0*tx0 ) * By );
	  const double Ay = sqrtFactor * ( -tx0 * ( ty0 * By + Bz ) + ( 1. + ty0*ty0 ) * Bx );

	  double x = x0 + tx0 * dz + 0.5 * k * invP * Ax * dz*dz;
	  double y = y0 + ty0 * dz + 0.5 * k * invP * Ay * dz*dz;
          
//	  const double tx = tx0 + invP * k * Ax * dz;
//	  const double ty = ty0 + invP * k * Ay * dz;

	  streamlog_out( DEBUG0 ) << "Old track state (x,y,tx,ty,invP):" << std::endl;
	  streamlog_out( DEBUG0 ) << std::setw(15) << x0
				  << std::setw(15) << y0
				  << std::setw(15) << tx0
				  << std::setw(15) << ty0
				  << std::setw(15) << invP << std::endl;
	  streamlog_out( DEBUG0 ) << "Old track ref point (x,y,z):" << std::endl;
	  streamlog_out( DEBUG0 ) << std::setw(15) << ts->getReferencePoint()[0]
				  << std::setw(15) << ts->getReferencePoint()[1]
				  << std::setw(15) << ts->getReferencePoint()[2] << std::endl;

	  const float newPos[] = {static_cast<float>(x), static_cast<float>(y), static_cast<float>(z0+dz)};
	  ts->setLocation( EUTelTrackStateImpl::AtOther );
	  ts->setReferencePoint( newPos );
          
	  streamlog_out( DEBUG0 ) << "New track ref point (x,y,z):" << std::endl;
	  streamlog_out( DEBUG0 ) << std::setw(15) << ts->getReferencePoint()[0]
				  << std::setw(15) << ts->getReferencePoint()[1]
				  << std::setw(15) << ts->getReferencePoint()[2] << std::endl;
                                  
          streamlog_out(DEBUG2) << "-------------------EUTelKalmanFilter::propagateTrackRefPoint()-------------------" << std::endl;
    }
    
    /** Find the hit closest to the intersection of a track with given sensor
     * 
     * @param ts track state
     * @return hit closest to the intersection of the track with the sensor plane
     * 
     */
    const EVENT::TrackerHit* EUTelKalmanFilter::findClosestHit( const EUTelTrackStateImpl* ts, int sensorID ) {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::findClosestHit()" << std::endl;
        hasSensorID needSensorID;
        needSensorID._requiredSensorID = sensorID;
        std::vector< MeasurementLayer* >::iterator itLayer = std::find_if (_allMeasurements.begin(), _allMeasurements.end(), needSensorID );
        
        if ( itLayer == _allMeasurements.end() ) {
            streamlog_out(DEBUG0) << "No hits in intersected plane " << sensorID << std::endl;
            return NULL;
        }
        
        EVENT::TrackerHitVec& hitInPlane = (*itLayer)->getHits();
        
        double maxDistance = std::numeric_limits<double>::max();
        EVENT::TrackerHitVec::const_iterator itClosestHit;
        
        EVENT::TrackerHitVec::const_iterator itHit;
        streamlog_out(DEBUG0) << "Hits in plane vector " << &hitInPlane << std::endl;
        streamlog_out(DEBUG0) << "N hits in plane " << sensorID << ": " << hitInPlane.size() << std::endl;
        for ( itHit = hitInPlane.begin(); itHit != hitInPlane.end(); ++itHit ) {
	    const double distance = getResidual( ts, *itHit ).Norm2Sqr();
            streamlog_out(DEBUG0) << "Distance^2 between hit and track intersection: " << distance << std::endl;
            if ( distance < maxDistance ) {
		itClosestHit = itHit;
		maxDistance = distance;
	    }
        }
        streamlog_out(DEBUG0) << "Minimal distance^2 between hit and track intersection: " << maxDistance << std::endl;
        streamlog_out(DEBUG2) << "----------------------EUTelKalmanFilter::findClosestHit()------------------------" << std::endl;

	return *itClosestHit;
    }
    
  double EUTelKalmanFilter::getXYPredictionPrecision( const EUTelTrackStateImpl* ts ) const {
      streamlog_out(DEBUG2) << "EUTelKalmanFilter::getXYPredictionPrecision()" << std::endl;
      
      TMatrixDSym Ckkm1 = getTrackStateCov(ts);
      //double xyPrec = _residualsRMax;//sqrt( Ckkm1[0][0]*Ckkm1[0][0] + Ckkm1[1][1]*Ckkm1[1][1] );
      double xyPrec = 1.0 / sqrt( Ckkm1[0][0]*Ckkm1[0][0] + Ckkm1[1][1]*Ckkm1[1][1] );
      
      streamlog_out(DEBUG0) << "Minimal combined UV resolution : " << xyPrec << std::endl;
      streamlog_out(DEBUG2) << "----------------------EUTelKalmanFilter::getXYPredictionPrecision()------------------------" << std::endl;

      return xyPrec;
    }

    /** Calculate track parameters propagation jacobian for given track state
     *  and propagation distance. The expressions were derived in parabolic approximation
     *  valid for small values of propagation distance |dz| < 10cm. Can be iterated if necessary.
     * 
     * @param ts track state
     * @param dz propagation distance
     * @return 
     */
    const TMatrixD& EUTelKalmanFilter::getPropagationJacobianF( const EUTelTrackStateImpl* ts, double dz ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getPropagationJacobianF()" << std::endl;
	// The formulas below are derived from equations of motion of the particle in
        // magnetic field under assumption |dz| small. Must be valid for |dz| < 10 cm

	const double mm = 1000.;
	const double k = 0.299792458/mm;

	// Get track parameters
	const double invP = ts->getInvP();
	const double x0 = ts->getX();
        const double y0 = ts->getY();
        const double z0 = ts->getReferencePoint()[2];
        const double tx0 = ts->getTx();
        const double ty0 = ts->getTy();

        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field
	const gear::BField&   B = geo::gGeometry().getMagneticField();
        const double Bx         = B.at( vectorGlobal ).x();
        const double By         = B.at( vectorGlobal ).y();
        const double Bz         = B.at( vectorGlobal ).z();
        
        const double sqrtFactor = sqrt( 1. + tx0*tx0 + ty0*ty0 );

	const double Ax = sqrtFactor * (  ty0 * ( tx0 * Bx + Bz ) - ( 1. + tx0*tx0 ) * By );
	const double Ay = sqrtFactor * ( -tx0 * ( ty0 * By + Bz ) + ( 1. + ty0*ty0 ) * Bx );

	// Partial derivatives
	//const double dAxdtx0 = tx0 * Ax / (sqrtFactor*sqrtFactor) + sqrtFactor*( ty0*Bx - 2. * tx0 * By );
	const double dAxdty0 = ty0 * Ax / (sqrtFactor*sqrtFactor) + sqrtFactor*( tx0*Bx + Bz );
	const double dAydtx0 = tx0 * Ay / (sqrtFactor*sqrtFactor) + sqrtFactor*( -ty0*By - Bz );
	//const double dAydty0 = ty0 * Ay / (sqrtFactor*sqrtFactor) + sqrtFactor*( -tx0*By + 2. * ty0 * Bx );

	const double dxdtx0 = dz;
	const double dxdty0 = 0.5 * invP * k * dz*dz * dAxdty0;

	const double dydtx0 = 0.5 * invP * k * dz*dz * dAydtx0;
	const double dydty0 = dz;

	const double dtxdty0 = invP * k * dz * dAxdty0;
	const double dtydtx0 = invP * k * dz * dAydtx0;

	const double dxdinvP0 = 0.5 * k * dz*dz * Ax;
	const double dydinvP0 = 0.5 * k * dz*dz * Ay;

	const double dtxdinvP0 = k * dz * Ax;
	const double dtydinvP0 = k * dz * Ay;

	// Fill-in matrix elements
	_jacobianF.UnitMatrix();
	_jacobianF[0][2] = dxdtx0;	_jacobianF[0][3] = dxdty0;	_jacobianF[0][4] = dxdinvP0;
	_jacobianF[1][2] = dydtx0;	_jacobianF[1][3] = dydty0;	_jacobianF[1][4] = dydinvP0;
	_jacobianF[2][3] = dtxdty0;	_jacobianF[2][4] = dtxdinvP0;
	_jacobianF[3][2] = dtydtx0;	_jacobianF[3][4] = dtydinvP0;
        
        if ( streamlog_level(DEBUG0) ){
             streamlog_out( DEBUG0 ) << "Propagation jacobian: " << std::endl;
            _jacobianF.Print();
        }
	
        streamlog_out( DEBUG2 ) << "-----------------------------EUTelKalmanFilter::getPropagationJacobianF()-------------------------------" << std::endl;
        
	return _jacobianF;
    }
    
    /**
     * Get extrapolated position of the track in global coordinate system
     * Extrapolation performed along the helix. Formulas are also valid for vanishing magnetic field.
     * 
     * @param ts track state
     * @param s arc length
     * @return 3d vector of coordinates in the global coordinate system
     */
    TVector3 EUTelKalmanFilter::getXYZfromArcLength1( const EUTelTrackStateImpl* ts, double s ) const {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::getXYZfromArcLength()" << std::endl;
        
        // Get starting track position
        const float* x = ts->getReferencePoint();
        const double x0 = x[0];
        const double y0 = x[1];
        const double z0 = x[2];
        
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticField();
        const double bx         = B.at( vectorGlobal ).x();
        const double by         = B.at( vectorGlobal ).y();
        const double bz         = B.at( vectorGlobal ).z();
        TVector3 hVec(bx,by,bz);
               
        TVector3 pVec = getPfromCartesianParameters( ts );

        const double p = pVec.Mag();
	const double mm = 1000.;
        const double k = -0.299792458/mm*_beamQ*hVec.Mag();
        const double rho = k/p;
               
        // Calculate end track position
	TVector3 pos;
	if ( fabs( k ) > 1.E-6  ) {
		// Non-zero magnetic field case
        	pos.SetX( x0 + 1./p * pVec.X() * s );
		pos.SetY( y0 + 1./k * pVec.Y() * sin( rho*s ) + 1./k * pVec.Z() * ( 1. - cos( rho*s ) ) );
                pos.SetZ( z0 + 1./k * pVec.Z() * sin( rho*s ) - 1./k * pVec.Y() * ( 1. - cos( rho*s ) ) );
        } else {
		// Vanishing magnetic field case
		const double cosA = cosAlpha( ts );
		const double cosB = cosBeta( ts );
		pos.SetX( x0 + cosA * s );
		pos.SetY( y0 + cosB * s );
		pos.SetZ( z0 + 1./p * pVec.Z() * s );
	}
        
        streamlog_out(DEBUG2) << "---------------------------------EUTelKalmanFilter::getXYZfromArcLength()------------------------------------" << std::endl;
        
        return pos;
    }

    /**
     * Get extrapolated position of the track in global coordinate system
     * Extrapolation performed along the helix. Formulas are also valid for vanishing magnetic field.
     * 
     * @param ts track state
     * @param s arc length
     * @return 3d vector of coordinates in the global coordinate system
     */
    TVector3 EUTelKalmanFilter::getXYZfromArcLength( const EUTelTrackStateImpl* ts, double s ) const {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::getXYZfromArcLength()" << std::endl;
        
        // Get starting track position
        const float* x = ts->getReferencePoint();
        const double x0 = x[0];
        const double y0 = x[1];
        const double z0 = x[2];
        
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticField();
        const double bx         = B.at( vectorGlobal ).x();
        const double by         = B.at( vectorGlobal ).y();
        const double bz         = B.at( vectorGlobal ).z();
        TVector3 hVec(bx,by,bz);
               
        TVector3 pVec = getPfromCartesianParameters( ts );

	const double H = hVec.Mag();
        const double p = pVec.Mag();
	const double mm = 1000.;
        const double k = -0.299792458/mm*_beamQ*H;
        const double rho = k/p;
        
        // Calculate end track position
	TVector3 pos( x0, y0, z0 );
	if ( fabs( k ) > 1.E-6  ) {
		// Non-zero magnetic field case
		TVector3 pCrossH = pVec.Cross(hVec.Unit());
		TVector3 pCrossHCrossH = pCrossH.Cross(hVec.Unit());
		const double pDotH = pVec.Dot(hVec.Unit());
		TVector3 temp1 = pCrossHCrossH;	temp1 *= ( -1./k * sin( rho * s ) );
		TVector3 temp2 = pCrossH;       temp2 *= ( -1./k * ( 1. - cos( rho * s ) ) );
		TVector3 temp3 = hVec;          temp3 *= ( pDotH / p * s );
		pos += temp1;
		pos += temp2;
		pos += temp3;
        } else {
		// Vanishing magnetic field case
		const double cosA = cosAlpha( ts );
		const double cosB = cosBeta( ts );
		pos.SetX( x0 + cosA * s );
		pos.SetY( y0 + cosB * s );
		pos.SetZ( z0 + 1./p * pVec.Z() * s );
	}
        
        streamlog_out(DEBUG2) << "---------------------------------EUTelKalmanFilter::getXYZfromArcLength()------------------------------------" << std::endl;
        
        return pos;
    }
    
    /**
     * Get extrapolated position of the track in global coordinate system
     * Formulas are also valid for vanishing magnetic field.
     * Calculation is based on numerical integration of equations of motion
     * 
     * @param ts track state
     * @param dz propagation distance along z
     * @return 3d vector of coordinates in the global coordinate system
     */
    TVector3 EUTelKalmanFilter::getXYZfromDzNum( const EUTelTrackStateImpl* ts, double dz ) const {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::getXYZfromDzNum()" << std::endl;
        
        // Get starting track position
        
        TVectorD trackStateVec = getTrackStateVec( ts );
	const float* x = ts->getReferencePoint();
        const double x0 = x[0];
        const double y0 = x[1];
        const double z0 = x[2];
        
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticField();
        const double bx         = B.at( vectorGlobal ).x();
        const double by         = B.at( vectorGlobal ).y();
        const double bz         = B.at( vectorGlobal ).z();
        TVector3 hVec(bx,by,bz);
        
        // Setup the equation
	trackStateVec[0] = x0;
	trackStateVec[1] = y0;
        static_cast< eom::EOMODE* >(_eomODE)->setInitValue( trackStateVec );
        static_cast< eom::EOMODE* >(_eomODE)->setBField( hVec );
        
        // Integrate
        TVectorD result = _eomIntegrator->integrate( dz );
        
        TVector3 pos(result[0],result[0],z0+dz);
        
        streamlog_out(DEBUG0) << "Result of the integration:" << std::endl;
        streamlog_message( DEBUG0, result.Print();, std::endl; );
        
        streamlog_out(DEBUG2) << "---------------------------------EUTelKalmanFilter::getXYZfromDzNum()------------------------------------" << std::endl;
        
        return pos;
    }

    /** Calculate cos of the angle between Z(beam) and X(solenoid field axis)
     *  from track parameters
     * 
     * @param ts track state
     * @return cos(alpha)
     */
    double EUTelKalmanFilter::cosAlpha( const EUTelTrackStateImpl* ts ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::cosAlpha()" << std::endl;
        const double tx = ts->getTx( );
        const double ty = ts->getTy( );
        const double cos = tx / sqrt( 1. + tx * tx + ty * ty );
        
        streamlog_out( DEBUG0 ) << "cosAlpha= " << cos << std::endl;
        
        return cos;
    }

    /** Calculate cos of the angle between Z(beam) and Y
     *  from track parameters
     * 
     * @param ts track state
     * @return cos(beta)
     */
    double EUTelKalmanFilter::cosBeta( const EUTelTrackStateImpl* ts ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::cosBeta()" << std::endl;
        const double tx = ts->getTx( );
        const double ty = ts->getTy( );
        const double cos = ty / sqrt( 1. + tx * tx + ty * ty );
        
        streamlog_out( DEBUG0 ) << "cosBeta= " << cos << std::endl;
        
        return ty / sqrt( 1. + tx * tx + ty * ty );
    }

    /** Convert track state to the vector object. Useful for matrix operations
     * 
     * @param ts track stare
     * @return vector of parameters
     */
    TVectorD EUTelKalmanFilter::getTrackStateVec( const EUTelTrackStateImpl* ts ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getTrackStateVec()" << std::endl;
        TVectorD x(5);
        x[0] = ts->getX();
        x[1] = ts->getY();
        x[2] = ts->getTx();
        x[3] = ts->getTy();
        x[4] = ts->getInvP();
        
        if ( streamlog_level(DEBUG0) ){
            streamlog_out( DEBUG0 ) << "Track state:" << std::endl;
            x.Print();
        }
        
        return x;
    }
    
    /** Convert track state parameter covariances to the matrix object. Useful for matrix operations
     * 
     * @param ts track state
     * @return covariance matrix
     */
    TMatrixDSym EUTelKalmanFilter::getTrackStateCov( const EUTelTrackStateImpl* ts ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getTrackStateVec()" << std::endl;
        TMatrixDSym C(5);
        
        const EVENT::FloatVec& trkCov = ts->getCovMatrix();
        
        C.Zero();
            
        C[0][0] = trkCov[0]; 
        C[1][0] = trkCov[1];  C[1][1] = trkCov[2]; 
        C[2][0] = trkCov[3];  C[2][1] = trkCov[4];  C[2][2] = trkCov[5]; 
        C[3][0] = trkCov[6];  C[3][1] = trkCov[7];  C[3][2] = trkCov[8];  C[3][3] = trkCov[9]; 
        C[4][0] = trkCov[10]; C[4][1] = trkCov[11]; C[4][2] = trkCov[12]; C[4][3] = trkCov[13]; C[4][4] = trkCov[14]; 
        
        if ( streamlog_level(DEBUG0) ){
            streamlog_out( DEBUG0 ) << "Track state covariance matrix:" << std::endl;
            C.Print();
        }
        
        return C;
    }
        
    /**
     * Propagate track state k-1 -> k
     * @param ts track state to update
     */
    void EUTelKalmanFilter::propagateTrackState( EUTelTrackStateImpl* ts ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::propagateTrackState()" << std::endl;
        TVectorD xkm1 = getTrackStateVec( ts );
        TVectorD xkkm1 = _jacobianF * xkm1;
        
        streamlog_message( DEBUG0, xkkm1.Print();, std::endl; );
        
        ts->setX( xkkm1[0] );
        ts->setY( xkkm1[1] );
        ts->setTx( xkkm1[2] );
        ts->setTy( xkkm1[3] );
        ts->setInvP( xkkm1[4] );
        
        streamlog_out( DEBUG2 ) << "-----------------------------------EUTelKalmanFilter::propagateTrackState()----------------------------------" << std::endl;
    }
    
    /** Retrieve hit covariance matrix from hit object. Useful for matrix operations
     * 
     * @param hit
     * @return hit covariance matrix
     */
    TMatrixDSym EUTelKalmanFilter::getHitCov( const EVENT::TrackerHit* hit ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getHitCov()" << std::endl;
        const EVENT::FloatVec uvcov = hit->getCovMatrix();
        TMatrixDSym V(2);
        V[0][0] = uvcov[0];                             //cov(x,x)
        V[1][0] = uvcov[1];   V[1][1] = uvcov[2];       //cov(y,x), cov(y,y)
        
        if ( streamlog_level(DEBUG0) ){
            streamlog_out( DEBUG0 ) << "Hit covariance matrix:" << std::endl;
            V.Print();
        }
        
        streamlog_out( DEBUG2 ) << "--------------------------------------EUTelKalmanFilter::getHitCov()-----------------------------------------" << std::endl;
        
        return V;
    }
    
    /** Calculate residual vector between given track and hit
     * 
     * @param ts track state
     * @param hit hit
     * @return 
     */
    TVectorD EUTelKalmanFilter::getResidual( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getResidual()" << std::endl;

        const double* uvpos = hit->getPosition();
        TVectorD mk(2); 
        mk[0] = uvpos[0];          mk[1] = uvpos[1];
        
	TVector3 trkPointVec = getXYZfromArcLength( ts, 0. );
        double trkPoint[] = { trkPointVec.X(), trkPointVec.Y(), trkPointVec.Z() };
	double trkPointLocal[] = {0.,0.,0.};
		    const int sensorID = Utility::GuessSensorID( hit );
	geo::gGeometry().master2Local(sensorID,trkPoint,trkPointLocal);
        
	TVectorD prediction(2);
	prediction[0] = trkPointLocal[0];	prediction[1] = trkPointLocal[1];    
        
        streamlog_out( DEBUG0 ) << "Hit (id=" << hit->id() << ") local(u,v) coordinates: (" << uvpos[0] << "," << uvpos[1] << ")" << std::endl;
        streamlog_out( DEBUG0 ) << "Prediction for hit (id=" << hit->id() << ") local(u,v) coordinates: (" 
				<< trkPointLocal[0] << "," << trkPointLocal[1] << ")" << std::endl;

/*
        TMatrixD Hk = getH(ts);
        if ( streamlog_level(DEBUG0) ){
            streamlog_out( DEBUG0 ) << "Projection Matrix Hk:" << std::endl;
            Hk.Print();
        }

        TVectorD xkkm1 = getTrackStateVec( ts );       // track state must be updated to x(k,k-1) not initial xk        
        xkkm1[0] -= xc;
        xkkm1[1] -= yc;

        TVectorD rk(2);
        rk = mk;        rk -= Hk * xkkm1;
*/        

        TVectorD rk(2);
        rk = mk - prediction;
        
        if ( streamlog_level(DEBUG0) ){
            streamlog_out( DEBUG0 ) << "Residual vector rk:" << std::endl;
            rk.Print();
        }
        
        streamlog_out( DEBUG2 ) << "----------------------------------EUTelKalmanFilter::getResidual()------------------------------------" << std::endl;
        
        return rk;
    }
    
    /** Retrieve residuals covariance
     * 
     * @param ts track state
     * @param hit
     * @return 
     */
    TMatrixDSym EUTelKalmanFilter::getResidualCov( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getResidualCov()" << std::endl;
        TMatrixD Hk = getH(ts);
        
        _processNoiseQ.Zero();
        TMatrixDSym Ckm1 = getTrackStateCov( ts );
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF );        //Ckkm1 += _processNoiseQ;       
        TMatrixDSym Rkkm1 = Ckkm1.Similarity(Hk);
        Rkkm1 += getHitCov(hit);
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Residual covariance matrix:" << std::endl;
            Rkkm1.Print();
        }
        
        streamlog_out( DEBUG2 ) << "-----------------------------------EUTelKalmanFilter::getResidualCov()------------------------------------" << std::endl;
        
        return Rkkm1;
    }
    
    /** Retrieve Kalman gain matrix
     * 
     * @param ts track state
     * @param hit
     * 
     * @return Gain matrix K
     */
    const TMatrixD& EUTelKalmanFilter::updateGainK( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::updateGainK()" << std::endl;
        _processNoiseQ.Zero();
        TMatrixDSym Ckm1 = getTrackStateCov( ts );
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF );        //Ckkm1 += _processNoiseQ;
        TMatrixD Ht(5,2);     Ht = Ht.Transpose( getH(ts) );
        
        _gainK = Ckkm1 * Ht * getResidualCov( ts, hit ).Invert();
//        _gainK *= Ht;
//        _gainK *= Ckkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Gain matrix:" << std::endl;
            _gainK.Print();
        }
        
        streamlog_out( DEBUG2 ) << "----------------------------------------EUTelKalmanFilter::updateGainK()------------------------------------" << std::endl;
        
        return _gainK;
    }
    
    /** Update track state given new hit
     * 
     * @param ts track state
     * @param hit
     */
    double EUTelKalmanFilter::updateTrackState( EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::updateTrackState()" << std::endl;
        TVectorD rkkm1(2);      rkkm1 = getResidual( ts, hit );
        TVectorD xk(5);         xk = getTrackStateVec( ts );
        
        TMatrixD Kk = updateGainK( ts, hit );
        xk += Kk * rkkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Updated track parameters:" << std::endl;
            xk.Print();
        }
        
        TMatrixD Hk = getH( ts );
        TMatrixD I(5,5);     I.UnitMatrix();
        I -= Kk*Hk;

        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Track parameters projection matrix Hk:" << std::endl;
            Hk.Print();
            streamlog_out( DEBUG0 ) << "Gain matrix Kk:" << std::endl;
            Kk.Print();
        }
        
        _processNoiseQ.Zero();
        TMatrixDSym Ckm1 = getTrackStateCov( ts );
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF );        //Ckkm1 += _processNoiseQ;
        
        TMatrixD Ck(5,5);
        Ck = I*Ckkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Updated track covariance matrix:" << std::endl;
            Ck.Print();
        }
        
        ts -> setX( xk[0] );
        ts -> setY( xk[1] );
        ts -> setTx( xk[2] );
        ts -> setTy( xk[3] );
        ts -> setInvP( xk[4] );
        
        float trkCov[15] = { static_cast<float>(Ck[0][0]), static_cast<float>(Ck[1][0]), static_cast<float>(Ck[1][1]), 
                             static_cast<float>(Ck[2][0]), static_cast<float>(Ck[2][1]), static_cast<float>(Ck[2][2]),
                             static_cast<float>(Ck[3][0]), static_cast<float>(Ck[3][1]), static_cast<float>(Ck[3][2]),
                             static_cast<float>(Ck[3][3]), static_cast<float>(Ck[4][0]), static_cast<float>(Ck[4][1]),
                             static_cast<float>(Ck[4][2]), static_cast<float>(Ck[4][3]), static_cast<float>(Ck[4][4]) };

        ts->setCovMatrix( trkCov );
        ts->setLocation( EUTelTrackStateImpl::AtOther );

        const float newPos[] = {static_cast<float>(xk[0]), static_cast<float>(xk[1]), static_cast<float>(ts->getReferencePoint()[2])};
        ts->setReferencePoint( newPos );
        
        double chi2 = 0.;
        
        TVectorD rk(2);
        TMatrixD Ir(2,2);        Ir.UnitMatrix();
        Ir -= Hk * Kk;
        rk = Ir * rkkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Residual vector rk:" << std::endl;
            rk.Print();
        }
        
        TMatrixD Rk(2,2);       Rk = getHitCov( hit );
        TMatrixD HkT(5,2);      HkT.Transpose( Hk );
        TMatrixD HkCkHkT(2,2);  HkCkHkT = Hk * Ck * HkT;
        
        Rk -= HkCkHkT;

        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Residual covariance matrix Rk:" << std::endl;
            Rk.Print();
        }
        
        chi2 = Rk.Invert().Similarity(rk);
        
        streamlog_out( DEBUG0 ) << "Chi2 contribution of the hit: " << chi2 << std::endl;
        
        streamlog_out( DEBUG2 ) << "------------------------------------------EUTelKalmanFilter::updateTrackState()----------------------------------" << std::endl;
        
        return chi2;
    }
    
    /** Retrieve track state projection onto measurement space matrix
     * 
     * @param ts track state
     * @return 
     */
    TMatrixD EUTelKalmanFilter::getH( const EUTelTrackStateImpl* ts ) const {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::getH()" << std::endl;
        TMatrixD H(2,5);
        H.Zero();
        TVector3 trkPointVec = getXYZfromArcLength( ts, 0. );
        double trkPoint[] = { trkPointVec.X(), trkPointVec.Y(), trkPointVec.Z() };
        const TGeoHMatrix* globalH = geo::gGeometry().getHMatrix( trkPoint );
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Local to global transformation matrix:" << std::endl;
            globalH->Print();
        }
        
        const TGeoHMatrix& globalHInv = globalH->Inverse();
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Global to local transformation matrix:" << std::endl;
            globalHInv.Print();
        }
        
//        const double* shift = globalHInv.GetTranslation();
        const double* rotation = globalHInv.GetRotationMatrix();

        // Fill necessary components
        H[0][0] = rotation[0]; // x projection, xx
        H[0][1] = rotation[1]; // y projection, xy
        H[1][0] = rotation[3]; // x projection, yx
        H[1][1] = rotation[4]; // y projection, yy

        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Matrix H:" << std::endl;
            H.Print();
        }
        
        return H;
    }

    /**
     * Convert EUTelTrackImpl representation to LCIO's TrackImpl class
     * 
     * @param track track to be converted
     * 
     * @return converted track object
     */
    IMPL::TrackImpl* EUTelKalmanFilter::cartesian2LCIOTrack( EUTelTrackImpl* track ) const {
        IMPL::TrackImpl* LCIOtrack = new IMPL::TrackImpl;

        const EUTelTrackStateImpl* ts = track->getTrackState( EUTelTrackStateImpl::AtLastHit );
        const float* r = ts->getReferencePoint();
        const double rx = r[0];
        const double ry = r[1];
        const double rz = r[2];
        
        LCIOtrack->setReferencePoint( r );
        // Assign hits to LCIO TRACK
        const EVENT::TrackerHitVec& trkcandhits = track->getTrackerHits();
        EVENT::TrackerHitVec::const_iterator itrHit;
        for ( itrHit = trkcandhits.begin(); itrHit != trkcandhits.end(); ++itrHit ) {
             LCIOtrack->addHit( *itrHit );
        }
                
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( rx, ry, rz );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticField();
        const double bx         = B.at( vectorGlobal ).x();
        const double by         = B.at( vectorGlobal ).y();
        const double bz         = B.at( vectorGlobal ).z();
        TVector3 h(bx,by,bz);
               
        TVector3 p = getPfromCartesianParameters( ts );

	//const double H = h.Mag();
	//const double mm = 1000.;
        //const double k = -0.299792458/mm*_beamQ*H;
        
        const double pt  = sqrt(p.Y()*p.Y() + p.Z()*p.Z());
        
        const double om  =    _beamQ/( pt );
        const double d0  =    0.;               // d0 in the plane transverse to magnetic field vector
        const double z0  =    0.;               // z0 along the magnetic field
        const double phi =    atan2( p.Y(), p.Z() );
        const double tanlam = p.X()/pt;
        
        const double chi2 = track->getChi2();
        LCIOtrack->setChi2( chi2 );
        LCIOtrack->setNdf( trkcandhits.size() ); // NDF contains number of hits used for the fit
        
        LCIOtrack->setOmega( om );
        LCIOtrack->setD0( d0 );
        LCIOtrack->setZ0( z0 );
        LCIOtrack->setPhi( phi );
        LCIOtrack->setTanLambda( tanlam );

        return LCIOtrack;

    }
 
    
    MeasurementLayer::MeasurementLayer() : _id(-1), _allHits() {}
    
    MeasurementLayer::MeasurementLayer( int id ) : _id(id), _allHits() {}
    
    MeasurementLayer::~MeasurementLayer() { }
    
    void MeasurementLayer::addHit( EVENT::TrackerHit* hit ) {
        _allHits.push_back( hit );
    }    



void EUTelKalmanFilter::setHits ( EVENT::TrackerHitVec& hits )
{
    streamlog_out ( DEBUG5 ) << "KalmanFilter: Setting hits" << std::endl;
    if ( hits.empty ( ) )
    {
	streamlog_out ( DEBUG5 ) << "Empty hit vector!" << std::endl;
	_isHitsOK = false;
	return;
    }
    _allHits = hits;
    _isHitsOK = sortHitsByMeasurementLayers ( _allHits );

}


bool EUTelKalmanFilter::sortHitsByMeasurementLayers ( const EVENT::TrackerHitVec& hits )
{

    streamlog_out ( DEBUG5 ) << "KalmanFilter: Sorting hits" << std::endl;
    bool isReady = true;

    if ( hits.empty ( ) )
    {
	streamlog_out ( WARNING2 ) << "No hits supplied in this event!" << std::endl;
	isReady = false;
	return isReady;
    }

    // clear information leftovers from previous events
    std::vector < MeasurementLayer* > ::iterator itLayer;
    for ( itLayer = _allMeasurements.begin ( ); itLayer != _allMeasurements.end ( ); ++itLayer )
    {
	delete ( *itLayer );
    }
    _allMeasurements.clear ( );

    // distribute hits to the measurement layers
    std::map < int, MeasurementLayer* > measLayers;
    std::map < int, MeasurementLayer* > ::iterator itMeasLayer;

    int sensorID = -1;
    EVENT::TrackerHitVec::const_iterator itr;
    for ( itr = hits.begin ( ); itr != hits.end ( ); ++itr )
    {
	sensorID = Utility::GuessSensorID ( *itr );
	itMeasLayer = measLayers.find ( sensorID );
	if ( itMeasLayer != measLayers.end ( ) )
	{
	    itMeasLayer -> second -> addHit ( *itr );
	}
	else
	{
	    measLayers[sensorID] = new MeasurementLayer ( sensorID );
	    measLayers[sensorID] -> addHit ( *itr );
	}
    }

    // sort measurement layers so that layers encountered by the track are first
    _allMeasurements = std::vector < MeasurementLayer* > ( geo::gGeometry ( ) .nPlanes ( ), NULL );
    int numberAlongZ = -1;
    std::vector < int > ::const_iterator itSensorID;
    for ( itSensorID = geo::gGeometry ( ) .sensorIDsVec ( ) .begin ( ); itSensorID != geo::gGeometry ( ) .sensorIDsVec ( ) .end ( ); ++itSensorID )
    {
	sensorID = ( *itSensorID );

	itMeasLayer = measLayers.find ( sensorID );
	if ( itMeasLayer != measLayers.end ( ) )
	{
	    numberAlongZ = geo::gGeometry ( ) .sensorIDtoZOrder ( sensorID );
	    _allMeasurements.at ( numberAlongZ ) = itMeasLayer -> second;
	}
    }

    // remove elements without MeasurementLayer assigned
    for ( itLayer = _allMeasurements.begin ( ); itLayer != _allMeasurements.end ( ); )
    {
	if ( (*itLayer) == NULL )
	{
	    streamlog_out ( DEBUG5 ) << "Removing a layer, no measurement attached!" << std::endl;
	    itLayer = _allMeasurements.erase ( itLayer );
	}
	else
	{
	    ++itLayer;
	}
    }

    return isReady;
}


bool EUTelKalmanFilter::initialise ( )
{
    streamlog_out ( DEBUG5 ) << "KalmanFilter: Initialisation" << std::endl;
    _isReady = true;

    // check the validity of the supplied beam energy 
    if ( _beamE < 1.E-6 )
    {
	streamlog_out ( ERROR5 ) << "Beam energy was set incorrectly!" << std::endl;
	_isReady = false;
	return _isReady;
    }

    if ( _beamEnergyUncertainty < 0.0 )
    {
	streamlog_out ( ERROR5 ) << "Beam uncertainty is negative. Check supplied values!" << std::endl;
	_isReady = false;
	return _isReady;
    }

    // Check validity of supplied hits
    if ( !_isHitsOK )
    {
	streamlog_out ( WARNING1 ) << "Hit collection is empty! Check supplied hits!" << std::endl;
	_isReady = false;
	return _isReady;
    }

    return _isReady;
}


void EUTelKalmanFilter::FitTracks ( )
{
    streamlog_out ( DEBUG5 ) << "KalmanFilter: Fitting tracks" << std::endl;

    // check if the fitter was configured correctly
    if ( !_isReady )
    {
	streamlog_out ( ERROR5 ) << "Can't fit. The fitter is not initialised!" << std::endl;
	return;
    }

    // flush data from previous event
    reset ( );

    // initialise Kalman filter states
    initialiseSeeds ( );

    // start Kalman filter
    std::vector < EUTelTrackImpl* > ::iterator itTrk;
    for ( itTrk = _tracksCartesian.begin ( ); itTrk != _tracksCartesian.end ( ); )
    {
	bool isGoodTrack = true;
	EUTelTrackStateImpl* state = const_cast < EUTelTrackStateImpl* > ( ( *itTrk ) -> getTrackState ( EUTelTrackStateImpl::AtFirstHit ) );

	LCCollectionVec *fitpointvec = new LCCollectionVec ( LCIO::TRACKERHIT );
	const float* xnew0 = ( state ) -> getReferencePoint ( );

	CellIDEncoder < TrackerHitImpl > fitHitEncoder ( EUTELESCOPE::HITENCODING, fitpointvec );
	TrackerHitImpl * fitpoint0 = new TrackerHitImpl;
	const int firstsensorid = geo::gGeometry ( ) .getSensorID ( xnew0 );
	fitHitEncoder["sensorID"] =  firstsensorid;
	fitHitEncoder["properties"] = kFittedHit;
	fitHitEncoder.setCellID ( fitpoint0 );
	double fitpos0[3] = { xnew0[0], xnew0[1], xnew0[2] };
	fitpoint0 -> setPosition ( fitpos0 );
	fitpointvec -> push_back ( fitpoint0 );
	setFitHits ( fitpointvec );
	( *itTrk ) -> addHit ( fitpoint0 );

	streamlog_out ( DEBUG1 ) << "Fit point: " << xnew0[0] << " " << xnew0[1] << " " << xnew0[2] << std::endl;

	double dz = findIntersection ( state );
	if ( dz < 0 )
	{
	    isGoodTrack = false;
	    delete ( *itTrk );
	    itTrk = _tracksCartesian.erase ( itTrk );
	}

	size_t hitcount = 1;

	// iterate until the particle flies out of the detector volume
	while ( dz > 0 )
	{

	    propagateTrackRefPoint ( state, dz );

	    // calculate track's new position
	    const float* xnew = ( state ) -> getReferencePoint ( );

	    streamlog_out ( DEBUG1 ) << "Fit point: " << xnew[0] << " " << xnew[1] << " " << xnew[2] << std::endl;

	    // determine id of the sensor in which track reference point is located
	    bool findhit = true;
	    const int newSensorID = geo::gGeometry ( ) .getSensorID ( xnew );
	    if ( newSensorID != -999 )
	    {
		findhit = true;
	    }
	    else
	    {
		streamlog_out ( DEBUG5 ) << "New point is outside of any sensor volume, removing this track candidate from further consideration." << std::endl;
		delete ( *itTrk );
		itTrk = _tracksCartesian.erase ( itTrk );
		findhit = false;
		isGoodTrack = false;
		break;
	    }

	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    // @TODO Allow missing hits
	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    EVENT::TrackerHit* closestHit = const_cast < EVENT::TrackerHit* > ( findClosestHit ( state, newSensorID ) );
	    if ( closestHit )
	    {
		const double distance = getResidual ( state, closestHit ) .Norm2Sqr ( );
		const double distanceCut = dz / 100.0 * getXYPredictionPrecision ( state );
		if ( distance > distanceCut )
		{
		    streamlog_out ( DEBUG5 ) << "Closest hit is outside of search window, skipping current plane." << std::endl;
		    //itTrk = _tracks.erase( itTrk );
		    findhit = false;
		    //isGoodTrack = false;
		    //break;
		}
		TVector3 test = getXYZfromArcLength(state, dz);
	    }
	    else
	    {
		streamlog_out ( DEBUG5 ) << "There are no hits in the closest plane." << std::endl;
		findhit = false;
	    }

	    if ( findhit )
	    {
		getPropagationJacobianF ( state, dz );
		propagateTrackState ( state );
		double chi2 = ( *itTrk ) -> getChi2 ( ) + updateTrackState ( state, closestHit );
		( *itTrk ) -> setChi2 ( chi2 );
		( *itTrk ) -> addHit ( closestHit );
		hitcount++;

		TrackerHitImpl * fitpoint1 = new TrackerHitImpl;
		fitHitEncoder["sensorID"] =  newSensorID;
		fitHitEncoder["properties"] = kFittedHit;
		fitHitEncoder.setCellID ( fitpoint1 );
		double fitpos[3] = { xnew[0], xnew[1], xnew[2] };
		fitpoint1 -> setPosition ( fitpos );
		fitpointvec -> push_back ( fitpoint1 );
		setFitHits ( fitpointvec );
		( *itTrk ) -> addHit ( fitpoint1 );
	    }
	    else
	    {
		getPropagationJacobianF ( state, dz );
		propagateTrackState( state );
	    }

	    // find next intersection
	    dz = findIntersection ( state );

	}

	if ( isGoodTrack && hitcount < geo::gGeometry ( ) .nPlanes ( ) - _allowedMissingHits)
	{
	    streamlog_out ( DEBUG5 ) << "Track candidate has " << hitcount << " hits, system has " << geo::gGeometry ( ) .nPlanes ( ) << " planes, " << _allowedMissingHits << " missing hits allowed." << std::endl;
	    streamlog_out ( DEBUG5) << "Too many missing hits, removing this track candidate from further consideration." << std::endl;
	    delete ( *itTrk );
	    itTrk = _tracksCartesian.erase ( itTrk );
	    isGoodTrack = false;
	}

	if ( isGoodTrack )
	{
	    streamlog_out ( DEBUG5 ) << "Found good track candidate!" << std::endl;
	    state -> setLocation ( EUTelTrackStateImpl::AtLastHit );
	    _tracks.push_back ( cartesian2LCIOTrack ( *itTrk ) );
	    delete ( *itTrk );
	    ++itTrk;
	}

    }

}


void EUTelKalmanFilter::reset ( )
{
    streamlog_out ( DEBUG5 ) << "KalmanFilter: resetting" << std::endl;
    _tracks.clear ( );
    _tracksCartesian.clear ( );
    _trackStates.clear ( );
}


void EUTelKalmanFilter::initialiseSeeds ( )
{
    streamlog_out ( DEBUG5 ) << "KalmanFilter: initialising seeds" << std::endl;

    if ( _allMeasurements.empty ( ) )
    {
	streamlog_out ( WARNING2 ) << "Can't initialise track seeds for the finder. No hits in this event." << std::endl;
	return;
    }

    const EVENT::TrackerHitVec& hitFirstLayer = _allMeasurements[0] -> getHits ( );
    streamlog_out ( DEBUG5 ) << "Number of hits in first non-empty layer: " << hitFirstLayer.size ( ) << std::endl;

    EVENT::TrackerHitVec::const_iterator itHit;
    for ( itHit = hitFirstLayer.begin ( ); itHit != hitFirstLayer.end ( ); ++itHit )
    {

	EUTelTrackStateImpl* state = new EUTelTrackStateImpl;

	// use beam direction as a seed track state
	const double* uvpos = ( *itHit ) -> getPosition ( );
	const int sensorID = Utility::GuessSensorID ( *itHit );
	double temp[] = { 0.0, 0.0, 0.0 };

	geo::gGeometry ( ) .local2Master ( sensorID, uvpos, temp );
	float posGlobal[] = { static_cast < float > ( temp[0] ), static_cast < float > ( temp[1] ), static_cast < float > ( temp[2] ) };

	gear::Vector3D vectorGlobal ( temp[0], temp[1], temp[2] );
	const double q = _beamQ;
	const double invp = q / _beamE;

	// fill track parameters covariance matrix
	float trkCov[15] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	const EVENT::FloatVec uvcov = ( *itHit ) -> getCovMatrix ( );

	// initial X,Y covariances are defined by seeding hit uncertainty multiplied by a huge factor
	trkCov[0] = uvcov[0] * 10000.0; // cov ( x, x )
	trkCov[1] = uvcov[1] * 10000.0; // cov ( y, x )
	trkCov[2] = uvcov[2] * 10000.0; // cov ( y, y )

	// TX,TY covariances are defined by the beam spread at the first plane
	if ( _beamAngularSpread[0] > 0.0 )
	{
	    trkCov[5] = _beamAngularSpread[0] * _beamAngularSpread[0]; // cov ( tx, x ) = 0, cov ( tx, y ) = 0, cov ( tx, tx ) = 0
	}
	else
	{
	    trkCov[5] = 100000.0;
	}

	if ( _beamAngularSpread[1] > 0.0 )
	{
	    trkCov[9] = _beamAngularSpread[1] * _beamAngularSpread[1]; // cov ( ty, x ) = 0, cov ( ty, y ) = 0, cov ( ty, tx ) = 0, cov ( ty, ty ) = 0
	}
	else
	{
	    trkCov[9] = 100000.0;
	}

	trkCov[14] = ( _beamEnergyUncertainty * _beamE ) * ( _beamEnergyUncertainty * _beamE ); // cov ( q/p, x ) = 0, cov ( q/p, y ) = 0, cov ( q/p, tx ) = 0, cov ( q/p, ty ) = 0, cov ( q/p, q/p ) = 0

	state -> setCovMatrix ( trkCov );
	state -> setReferencePoint ( posGlobal );
	state -> setLocation ( EUTelTrackStateImpl::AtFirstHit );
	state -> setTx ( 0.0 ); // seed with 0 at first hit given by beam direction
	state -> setTy ( 0.0 ); // seed with 0 at first hit given by beam direction
	state -> setX ( posGlobal[0] ); // 0.0 at first hit
	state -> setY ( posGlobal[1] ); // 0.0 at first hit
	state -> setInvP ( invp ); // independent of reference point

	EUTelTrackImpl* track = new EUTelTrackImpl;
	track -> addHit ( ( *itHit ) );
	track -> addTrackState ( state );
	_tracksCartesian.push_back ( track );
    }

    streamlog_out ( DEBUG5 ) << "KalmanFilter: done initialising seeds" << std::endl;
}


double EUTelKalmanFilter::findIntersection ( EUTelTrackStateImpl* ts )
{
    streamlog_out ( DEBUG5 ) << "KalmanFilter: finding intersection" << std::endl;

    // get track position
    const float* x = ts -> getReferencePoint ( );
    const double x0 = x[0];
    const double y0 = x[1];
    const double z0 = x[2];
    TVector3 trkVec ( x0, y0, z0 );

    // get magnetic field vector
    gear::Vector3D vectorGlobal ( x0, y0, z0 );
    const gear::BField& B = geo::gGeometry ( ) .getMagneticField ( );
    const double bx = B.at ( vectorGlobal ) .x ( );
    const double by = B.at ( vectorGlobal ) .y ( );
    const double bz = B.at ( vectorGlobal ) .z ( );
    TVector3 hVec ( bx, by, bz );
    const double H = hVec.Mag ( );

    // calculate track momentum from track parameters
    TVector3 pVec = getPfromCartesianParameters ( ts );
    const double p = pVec.Mag ( );
    const double mm = 1000.0;
    const double k = -0.299792458 / mm * _beamQ * H;
    const double rho = k / p;

    // determine id of the sensor in which track reference point is located
    int sensorID = geo::gGeometry ( ) .getSensorID ( x );
    int sensorZorder = geo::gGeometry ( ) .sensorIDtoZOrder ( sensorID );

    streamlog_out ( DEBUG5 ) << "SensorID: " << sensorID << std::endl;
    streamlog_out ( DEBUG5 ) << "Sensor Z order: " << sensorZorder << std::endl;

    if ( sensorID < 0 )
    {
	streamlog_out ( DEBUG5 ) << "Track interseciton was not found, sensorID < 0" << std::endl;
	return -999.0;
    }

    // get plane normals in global reference system
    std::map < int, TVector3 > planesNorm;
    EVENT::IntVec sensID = geo::gGeometry ( ) .sensorIDsVec ( );
    EVENT::IntVec::const_iterator itPlaneId;
    for ( itPlaneId = sensID.begin ( ); itPlaneId != sensID.end ( ); ++itPlaneId )
    {
	TVector3 norm = geo::gGeometry ( ) .siPlaneNormal ( *itPlaneId );
	if ( norm.Mag2 ( ) > 1e-6 )
	{
	    planesNorm[*itPlaneId] = norm;
	}
	else
	{
	    streamlog_out( ERROR5 ) << "Wrong sensor normal vector: ID: " << *itPlaneId << " x: " << norm.X ( ) << " y: " << norm.Y ( ) << " z: " << norm.Z ( ) << std::endl;
	}
    }

    // next plane downstream
    const int nextPlaneId = geo::gGeometry ( ) .sensorZOrderToID ( sensorZorder + 1 );
    if ( nextPlaneId > 0 )
    {
	itPlaneId = std::find ( sensID.begin ( ), sensID.end ( ), nextPlaneId ); 
    }
    else
    {
	streamlog_out ( DEBUG5 ) << "No following plane found!" << std::endl;
	return -999.0;
    }

    TVector3 sensorCenter ( geo::gGeometry ( ) .siPlaneXPosition ( *itPlaneId ), geo::gGeometry ( ) .siPlaneYPosition ( *itPlaneId ), geo::gGeometry ( ) .siPlaneZPosition ( *itPlaneId ) );
    TVector3 delta = trkVec - sensorCenter;
    TVector3 pVecCrosH = pVec.Cross ( hVec.Unit ( ) );

    streamlog_out ( DEBUG5 ) << "Current point (X,Y,Z): " << std::setw ( 15 ) << x0 << std::setw ( 15 ) << y0 << std::setw ( 15 ) << z0 << std::endl;
    streamlog_out ( DEBUG5 ) << "PlaneID: " << *itPlaneId << std::endl;
    streamlog_out ( DEBUG5 ) << "Normal vector: " << std::endl;
    //planesNorm[*itPlaneId].Print ( );
    streamlog_out ( DEBUG5 ) << "P x H vector: " << std::endl;
    //pVecCrosH.Print ( );
    streamlog_out ( DEBUG5 ) << "Rho: " << rho << std::endl;
    streamlog_out ( DEBUG5 ) << "P: " << p << std::endl;

    const double a = -0.5 * rho * ( planesNorm[*itPlaneId].Dot ( pVecCrosH ) ) / p;
    const double b = planesNorm[*itPlaneId].Dot ( pVec ) / p;
    const double c = planesNorm[*itPlaneId].Dot ( delta );

    // solutions are sorted in ascending order
    std::vector < double > sol = Utility::solveQuadratic (a, b, c);

    TVector3 newPos[2];
    newPos[0] = getXYZfromArcLength ( ts, sol[0] );
    newPos[1] = getXYZfromArcLength ( ts, sol[1] );

    streamlog_out ( DEBUG5 ) << "Next intersected volume can be: " << *itPlaneId << std::endl;
    streamlog_out ( DEBUG5 ) << "Plane: " << *itPlaneId << std::endl;
    streamlog_out ( DEBUG5 ) << "Solutions for arc length: " << std::setw ( 15 ) << sol[0] << std::setw ( 15 ) << sol[1] << std::endl;
    streamlog_out ( DEBUG5 ) << "First solution (X,Y,Z):  " << std::setw ( 15 ) << newPos[0].X ( ) << std::setw ( 15 ) << newPos[0].Y ( ) << std::setw ( 15 ) << newPos[0].Z ( ) << std::endl;
    streamlog_out ( DEBUG5 ) << "Second solution (X,Y,Z): " << std::setw ( 15 ) << newPos[1].X ( ) << std::setw ( 15 ) << newPos[1].Y ( ) << std::setw ( 15 ) << newPos[1].Z ( ) << std::endl;

    // choose solution with minimal positive arc length. It will correspond to the closest point along the helix
    double solution = ( sol[0] > 0. ) ? sol[0] : ( ( sol[0] < 0. && sol[1] > 0. ) ? sol[1] : -1.0 );
    int solutionNum = ( sol[0] > 0. ) ? 0 : ( ( sol[0] < 0. && sol[1] > 0. ) ? 1 : -1 );
    double dz = ( solution > 0. ) ? newPos[ solutionNum ].Z ( ) - trkVec.Z ( ) : -1.0 ;

    if ( dz < 1.E-6 )
    {
	streamlog_out ( DEBUG5 ) << "Track intersection was not found" << std::endl;
	return -999.0;
    }

    return dz;
}


TVector3 EUTelKalmanFilter::getPfromCartesianParameters ( const EUTelTrackStateImpl* ts ) const
{
    streamlog_out ( DEBUG5 ) << "KalmanFilter getting p from cartesian parameters" << std::endl;
    const double p  = 1.0 / ts -> getInvP ( ) * _beamQ;
    const double tx = ts -> getTx ( );
    const double ty = ts -> getTy ( );
    const double px = p * tx / sqrt ( 1.0 + tx * tx + ty * ty );
    const double py = p * ty / sqrt ( 1.0 + tx * tx + ty * ty );
    const double pz = p / sqrt ( 1.0 + tx * tx + ty * ty );

    return TVector3 ( px, py, pz );
}
