/* 
 * File:   EUTelMagneticFieldFinder.cc
 *
 * Created on July 2, 2013, 12:59 PM
 */

#include "EUTelMagneticFieldFinder.h"

#include "streamlog/streamlog.h"

#include "EUTelGeometryTelescopeGeoDescription.h"

#include "gear/gearimpl/Vector3D.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif

#include <iostream>
#include <functional>
#include <algorithm>

namespace eutelescope {

    EUTelKalmanFilter::EUTelKalmanFilter() : EUTelTrackFitter( "KalmanFilter" ), 
            _tracks(), 
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
            _gainK(5,5),
            _residualCovR(2,2){}

    EUTelKalmanFilter::EUTelKalmanFilter( std::string name ) : EUTelTrackFitter( name ),
            _tracks(),
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
            _gainK(5,5),
            _residualCovR(2,2){}

    EUTelKalmanFilter::~EUTelKalmanFilter() {}
    
    /** Perform Kalman filter track search and track fit */
    void EUTelKalmanFilter::FitTracks() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::FitTracks()" << std::endl;
        // Check if the fitter was configured correctly
        if ( !_isReady ) {
            streamlog_out(ERROR1) << _name << ": Can't fit. The fitter is not initialised!" << std::endl;
            return;
        }
        
        // Initialise Kalman filter states
        initialiseSeeds();

	// 
	findIntersection( _trackStates[0] );
	//
	propagateTrackState( _trackStates[0], 750. );
	//
	getPropagationJacobianF( _trackStates[0], 750. );
	_jacobianF.Print();
	//
	
    }
    
    /** Check validity of the user input */
    bool EUTelKalmanFilter::initialise() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::initialise()" << std::endl;
        
        _isReady = true;
        
        // Check the validity of supplied beam energy 
        if ( _beamE < 1.E-6 ) {
            streamlog_out(ERROR1) << "Beam direction was set incorrectly" << std::endl;
            _isReady = false;
            return _isReady;
        }
        
        if ( _beamEnergyUncertainty < 0 ) {
            streamlog_out(ERROR1) << "Beam uncertainty is negative. Check supplied values" << std::endl;
            _isReady = false;
            return _isReady;
        }
        
        // Check validity of supplied hits
        if ( !_isHitsOK ) {
            streamlog_out(ERROR1) << "Collection of hits is empty. Check supplied hits" << std::endl;
            _isReady = false;
            return _isReady;
        }
        
        streamlog_out(DEBUG2) << "Initialisation successfully completed" << std::endl;
        return _isReady;
    }
    
    void EUTelKalmanFilter::reset() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::reset()" << std::endl;
    }
    
    void EUTelKalmanFilter::setHits( EVENT::TrackerHitVec& hits ) {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::setHits()" << std::endl;
        if ( hits.empty() ) {
            _isHitsOK = false;
            return;
        }
        _allHits = hits;
        _isHitsOK = sortHitsByMeasurementLayers(_allHits);

    }
    
    /** Generate seed track states necessary to
     * start Kalman filter
     * 
     * Generate as many starting states as number of hits in first telescope plane
     * plus one additional state for missing hit in the fist plane if allowed
     * 
     */
    void EUTelKalmanFilter::initialiseSeeds() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::initialiseSeeds()" << std::endl;

        if ( _allMeasurements.empty() ) {
            streamlog_out(WARNING1) << "Can't initialise track seeds for the finder. No hits in this event." << std::endl;
            return;
        }
        
        const EVENT::TrackerHitVec& hitFirstLayer = _allMeasurements[0]->getHits();
        streamlog_out(DEBUG1) << "N hits in first non-empty layer: " << hitFirstLayer.size() << std::endl;
        EVENT::TrackerHitVec::const_iterator itHit;
        for ( itHit = hitFirstLayer.begin(); itHit != hitFirstLayer.end(); ++itHit ) {
            // The object will be owned by LCIO collection. Must not to free.
            EUTelTrackStateImpl* state = new EUTelTrackStateImpl;
//
//            // Use beam direction as a seed track state
            const double* uvpos = (*itHit)->getPosition();
            const int sensorID = Utility::GuessSensorID( *itHit );
            double temp[] = {0.,0.,0.};
            geo::gGeometry().local2Master( sensorID, uvpos, temp);
            float posGlobal[] = { temp[0], temp[1], temp[2] };
            
//            std::cout << std::fixed;
//            std::cout << "Local coordinates ( sensorID =  " << sensorID << " ) : " << std::endl;
//            std::cout << std::setw(10) << std::setprecision(5) << uvpos[0] << std::setw(10) << std::setprecision(5) << uvpos[1] << std::setw(10) << std::setprecision(5) << uvpos[2] << std::endl;
//            std::cout << "Global coordinates ( sensorID =  " << sensorID << " ) : " << std::endl;
//            std::cout << std::setw(10) << std::setprecision(5) << xxx[0] << std::setw(10) << std::setprecision(5) << xxx[1] << std::setw(10) << std::setprecision(5) << xxx[2] << std::endl;
            
            gear::Vector3D vectorGlobal( temp[0], temp[1], temp[2] );
            const double q          = _beamQ;      // assume electron beam
            const double invp       = fabs(q)/_beamE;
            
            // Fill track parameters covariance matrix
            float trkCov[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            
            const EVENT::FloatVec uvcov = (*itHit)->getCovMatrix();
            // X,Y covariances are defined by seeding hit uncertainty
            trkCov[0] = uvcov[0];                               //cov(x,x)
            trkCov[1] = uvcov[1];   trkCov[2] = uvcov[2];       //cov(y,x), cov(y,y)
            // TX,TY covariances are defined by beam spread at the first plane
            
            if ( _beamAngularSpread[0] > 0. ) trkCov[5] = _beamAngularSpread[0] * _beamAngularSpread[0];          //cov(tx,x)=0, cov(tx,y)=0, cov(tx,tx)
            else trkCov[5] = 1.E5;
            if ( _beamAngularSpread[1] > 0. ) trkCov[9] = _beamAngularSpread[1] * _beamAngularSpread[1];          //cov(ty,x)=0, cov(ty,y)=0, cov(ty,tx)=0, cov(ty,ty)
            else trkCov[9] = 1.E5;
            trkCov[14] = ( _beamEnergyUncertainty * _beamE ) * ( _beamEnergyUncertainty * _beamE );               //cov(q/p,x)=0, cov(q/p,y)=0, cov(q/p,tx)=0, cov(q/p,ty)=0, cov(q/p,q/p)
            
            state->setCovMatrix(trkCov);          
            state->setReferencePoint( posGlobal );
            state->setLocation( EUTelTrackStateImpl::AtFirstHit );
            state->setTx(0.);         // seed with 0 at first hit. Given by beam direction.
            state->setTy(0.);         // seed with 0 at first hit. Given by beam direction.
            state->setX(posGlobal[0]);          // 0. at first hit
            state->setY(posGlobal[1]);          // 0. at first hit
            state->setInvP(invp);            // independent of reference point
            _trackStates.push_back( state );
        }
    }
    
    /** Calculate track momentum from track parameters 
     * @param ts track state with specified tx,ty,x,y,invP
     * @return 3-vector of momentum
     */
    TVector3 EUTelKalmanFilter::getPfromCartesianParameters( const EUTelTrackStateImpl* ts ) const {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::getPfromCartesianParameters()" << std::endl;
        const double p  = 1. / ts->getInvP() * fabs( _beamQ );
        const double tx = ts->getTx();
        const double ty = ts->getTy();
        const double px = p*tx / sqrt( 1. + tx*tx + ty*ty );
        const double py = p*ty / sqrt( 1. + tx*tx + ty*ty );
        const double pz = p    / sqrt( 1. + tx*tx + ty*ty );
        
        return TVector3(px,py,pz);
    }
    
    std::vector< double > EUTelKalmanFilter::findIntersection( EUTelTrackStateImpl* ts ) const {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::findIntersection()" << std::endl;
        
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( 0.,0.,0. );        // assuming uniform magnetic field running along X direction
        const double bx         = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).x();
        const double by         = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).y();
        const double bz         = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).z();
        TVector3 hVec(bx,by,bz);
        
        // Calculate track momentum from track parameters
        TVector3 pVec = getPfromCartesianParameters( ts );
        const double p = pVec.Mag();
        const double k = -0.299792458*_beamQ*hVec.Mag();
        const double rho = k/p;
        
        // Calculate track position
        const float* x = ts->getReferencePoint();
        const double x0 = x[0];
        const double y0 = x[1];
        const double z0 = x[2];
        TVector3 trkVec(x0,y0,z0);
        
        // Get planes normals in global reference system
        std::map< int, TVector3 > planesNorm;
        EVENT::IntVec sensID = geo::gGeometry().sensorIDsVec();
        EVENT::IntVec::const_iterator itPlaneId;
        for ( itPlaneId = sensID.begin(); itPlaneId != sensID.end(); ++itPlaneId ) {
            TVector3 norm = geo::gGeometry().siPlaneNormal( *itPlaneId );
            if ( norm.Mag2() > 1e-6 ) planesNorm[ *itPlaneId ] = norm;
            else {
                streamlog_out( ERROR0 ) << "Wrong sensor normal vector:\n";
                streamlog_out( ERROR0 ) << "ID: " << *itPlaneId << "\n";
                streamlog_out( ERROR0 ) << "x: " << norm.X() << " y: " << norm.Y() << " z: " << norm.Z()<< std::endl;
            }
        }
        
        // Construct quadratic equation
        for ( itPlaneId = sensID.begin(); itPlaneId != sensID.end(); ++itPlaneId ) {
            TVector3 sensorCenter( geo::gGeometry().siPlaneXPosition( *itPlaneId ),
                                   geo::gGeometry().siPlaneYPosition( *itPlaneId ),
                                   geo::gGeometry().siPlaneZPosition( *itPlaneId ) );
            TVector3 delta = trkVec - sensorCenter;
            TVector3 pVecCrosH = pVec.Cross( hVec );

	     streamlog_out (DEBUG0) << "-------------------------------------------" << std::endl;
	     streamlog_out (DEBUG0) << "Current point (X,Y,Z): " << std::setw(15) << x0 
								 << std::setw(15) << y0 
								 << std::setw(15) << z0 << std::endl;
	     streamlog_out (DEBUG0) << "PlaneID: " << *itPlaneId << std::endl;
	     streamlog_out (DEBUG0) << "Normal vector" << std::endl;
	     planesNorm[*itPlaneId].Print();
	     streamlog_out (DEBUG0) << "P x H vector" << std::endl;
	     pVecCrosH.Print();
	     streamlog_out (DEBUG0) << "Rho: " << rho << std::endl;
	     streamlog_out (DEBUG0) << "P: " << p << std::endl;

            const double a = -0.5 * rho * ( planesNorm[*itPlaneId].Dot( pVecCrosH ) ) / p;
            const double b = planesNorm[*itPlaneId].Dot( pVec ) / p;
            const double c = planesNorm[*itPlaneId].Dot( delta );
            
            std::vector< double > sol = Utility::solveQuadratic(a,b,c);
            
            TVector3 newPos1 = getXYZfromArcLenght( ts, sol[0] );
            TVector3 newPos2 = getXYZfromArcLenght( ts, sol[1] );
	    
	    streamlog_out (DEBUG0) << "Plane: " << *itPlaneId << std::endl;
	    streamlog_out (DEBUG0) << "Solutions for arc length: " << std::setw(15) << sol[0] << std::setw(15) << sol[1] << std::endl;
	    streamlog_out (DEBUG0) << "First solution (X,Y,Z): " << std::setw(15) << newPos1.X() 
								 << std::setw(15) << newPos1.Y() 
								 << std::setw(15) << newPos1.Z() << std::endl;
	    streamlog_out (DEBUG0) << "Second solution (X,Y,Z): " << std::setw(15) << newPos2.X() 
								 << std::setw(15) << newPos2.Y() 
								 << std::setw(15) << newPos2.Z() << std::endl;

//            TVector3 newPos = ( newPos1.Z() > 0 ) ? newPos1 : newPos2;
        }
	   std::vector< double > newPos;
	   newPos.push_back(1.); newPos.push_back(1.); newPos.push_back(1.);
           return newPos; 
    }

    /** Propagate track state by dz */
    void EUTelKalmanFilter::propagateTrackState( EUTelTrackStateImpl* ts, double dz ) {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::propagateTrackState()" << std::endl;
          // The formulas below are derived from equations of motion of the particle in 
          // magnetic field under assumption |dz| small. Must be valid for |dz| < 10 cm

          // Get magnetic field vector
          gear::Vector3D vectorGlobal( 0.,0.,0. );        // assuming uniform magnetic field running along X direction
          const double Bx = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).x();
          const double By = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).y();
          const double Bz = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).z();
	  const double k = 0.299792458;

	  // Get track parameters
	  const double invP = ts->getInvP();
	  const double x0 = ts->getX();
	  const double y0 = ts->getY();
	  const double z0 = ts->getReferencePoint()[2];
	  const double tx0 = ts->getTx();
	  const double ty0 = ts->getTy();
	  const double sqrtFactor = sqrt( 1. + tx0*tx0 + ty0*ty0 );

	  const double Ax = sqrtFactor * (  ty0 * ( tx0 * Bx + Bz ) - ( 1. + tx0*tx0 ) * By );
	  const double Ay = sqrtFactor * ( -tx0 * ( ty0 * By + Bz ) + ( 1. + ty0*ty0 ) * Bx );

	  const double x = x0 + tx0 * dz + 0.5 * k * invP * Ax * dz*dz;
	  const double y = y0 + ty0 * dz + 0.5 * k * invP * Ay * dz*dz;

	  const double tx = tx0 + invP * k * Ax * dz;
	  const double ty = ty0 + invP * k * Ay * dz;

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

	  const float newPos[] = {x, y, z0+dz};
	  ts->setLocation( EUTelTrackStateImpl::AtOther );
	  ts->setReferencePoint( newPos );
	  ts->setX(x);
	  ts->setY(y);
	  ts->setTx(tx);
	  ts->setTy(ty);

	  streamlog_out( DEBUG0 ) << "New track state (x,y,tx,ty,invP):" << std::endl;
	  streamlog_out( DEBUG0 ) << std::setw(15) << x
				  << std::setw(15) << y
				  << std::setw(15) << tx
				  << std::setw(15) << ty
				  << std::setw(15) << invP << std::endl;
	  streamlog_out( DEBUG0 ) << "New track ref point (x,y,z):" << std::endl;
	  streamlog_out( DEBUG0 ) << std::setw(15) << ts->getReferencePoint()[0]
				  << std::setw(15) << ts->getReferencePoint()[1]
				  << std::setw(15) << ts->getReferencePoint()[2] << std::endl;
    }

    const TMatrixD& EUTelKalmanFilter::getPropagationJacobianF( const EUTelTrackStateImpl* ts, double dz ) {
	  // The formulas below are derived from equations of motion of the particle in
          // magnetic field under assumption |dz| small. Must be valid for |dz| < 10 cm

          // Get magnetic field vector
          gear::Vector3D vectorGlobal( 0.,0.,0. );        // assuming uniform magnetic field running along X direction
          const double Bx = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).x();
          const double By = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).y();
          const double Bz = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).z();
	  const double k = 0.299792458;

	  // Get track parameters
	  const double invP = ts->getInvP();
	  //const double x0 = ts->getX();
          //const double y0 = ts->getY();
          //const double z0 = ts->getReferencePoint()[2];
          const double tx0 = ts->getTx();
          const double ty0 = ts->getTy();
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
	
	return _jacobianF;
    }
    
    TVector3 EUTelKalmanFilter::getXYZfromArcLenght( const EUTelTrackStateImpl* ts, double s ) const {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::getXYZfromArcLenght()" << std::endl;
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( 0.,0.,0. );        // assuming uniform magnetic field running along X direction
        const double bx         = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).x();
        const double by         = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).y();
        const double bz         = geo::gGeometry().getMagneticFiled().at( vectorGlobal ).z();
        TVector3 hVec(bx,by,bz);
               
        TVector3 pVec = getPfromCartesianParameters( ts );

        const double p = pVec.Mag();
        const double k = -0.299792458*_beamQ*hVec.Mag();
        const double rho = k/p;
        
        // Calculate starting track position
        const float* x = ts->getReferencePoint();
        const double x0 = x[0];
        const double y0 = x[1];
        const double z0 = x[2];
        
        // Calculate end track position
	TVector3 pos;
	if ( fabs( k ) > 1.E-6  ) {
		// Non-zero magnetic field case
        	pos.SetX( x0 + 1./p * pVec.X() * s);
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
        return pos;
    }

	double EUTelKalmanFilter::cosAlpha( const EUTelTrackStateImpl* ts) const {
		const double tx = ts->getTx();
		const double ty = ts->getTy();
		return tx / sqrt( 1. + tx*tx + ty*ty );
	}

	double EUTelKalmanFilter::cosBeta( const EUTelTrackStateImpl* ts) const {
		const double tx = ts->getTx();
		const double ty = ts->getTy();
		return ty / sqrt( 1. + tx*tx + ty*ty );
	}

    TVectorD EUTelKalmanFilter::getTrackStateVec( const EUTelTrackStateImpl* ts ) const {
        TVectorD x(5);
        x[0] = ts->getX();
        x[1] = ts->getY();
        x[2] = ts->getTx();
        x[3] = ts->getTy();
        x[4] = ts->getInvP();
        
        return x;
    }
    
    TMatrixDSym EUTelKalmanFilter::getTrackStateCov( const EUTelTrackStateImpl* ts ) const {
        TMatrixDSym C(5);
        
        const EVENT::FloatVec& trkCov = ts->getCovMatrix();
        
        C.Zero();
            
        C[0][0] = trkCov[0]; 
        C[1][0] = trkCov[1];  C[1][1] = trkCov[2]; 
        C[2][0] = trkCov[3];  C[2][1] = trkCov[4];  C[2][2] = trkCov[5]; 
        C[3][0] = trkCov[6];  C[3][1] = trkCov[7];  C[3][2] = trkCov[8];  C[3][3] = trkCov[9]; 
        C[4][0] = trkCov[10]; C[4][1] = trkCov[11]; C[4][2] = trkCov[12]; C[4][3] = trkCov[13]; C[4][4] = trkCov[14]; 
        
        return C;
    }
        
    /**
     * Propagate track state k-1 -> k
     * @param ts track state to update
     */
    void EUTelKalmanFilter::propagateTrackState( EUTelTrackStateImpl* ts ) {
        TVectorD xkm1 = getTrackStateVec( ts );
        TVectorD xkkm1 = _jacobianF * xkm1;
        
        _processNoiseQ.Zero();
        TMatrixDSym Ck = getTrackStateCov( ts );
        TMatrixDSym Ckkm1 = Ck.Similarity( _jacobianF );        Ckkm1 += _processNoiseQ;
        
    }
    
    TMatrixD EUTelKalmanFilter::getHitCov( const EVENT::TrackerHit* hit ) const {
        const EVENT::FloatVec uvcov = hit->getCovMatrix();
        TMatrixDSym V(2);
        V[0][0] = uvcov[0];                             //cov(x,x)
        V[1][0] = uvcov[1];   V[1][1] = uvcov[2];       //cov(y,x), cov(y,y)
        
        return V;
    }
    
    TVectorD EUTelKalmanFilter::getResidual( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) const {
        const double* uvpos = hit->getPosition();
        const int sensorID = Utility::GuessSensorID( hit );
        double temp[] = {0.,0.,0.};
        geo::gGeometry().local2Master( sensorID, uvpos, temp );
        float posGlobal[] = { temp[0], temp[1], temp[2] };
        
        const TVectorD xkkm1 = getTrackStateVec( ts );
        TVectorD mk(2); 
        mk[0] = hit->getPosition()[0];          mk[1] = hit->getPosition()[1];
    
        TMatrixD Hk = getH(ts);
        
        TVectorD rk(2);
        rk = mk;        rk -= Hk * xkkm1;
        
        return rk;
    }
    
    TMatrixD EUTelKalmanFilter::getResidualCov( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) const {
        TMatrixD Hk = getH(ts);
        TMatrixDSym Ck = getTrackStateCov( ts );
        TMatrixDSym Ckkm1 = Ck.Similarity( _jacobianF );        Ckkm1 += _processNoiseQ;
        
        Ckkm1.Similarity(Hk);
        
        TMatrixD Rkkm1 = getHitCov(hit);
        Rkkm1 += Ckkm1;
        return Rkkm1;
    }
    
    
    TMatrixD EUTelKalmanFilter::getH( const EUTelTrackStateImpl* ts ) const {
        TMatrixD H(2,5);
        H.Zero();
        TVector3 trkPointVec = getXYZfromArcLenght( ts, 0. );
        double trkPoint[] = { trkPointVec.X(), trkPointVec.Y(), trkPointVec.Z() };
        TGeoHMatrix* globalH = geo::gGeometry().getHMatrix( trkPoint );
        const TGeoHMatrix& globalHInv = globalH->Inverse();
        
        const double* shift = globalHInv.GetTranslation();
        const double* rotation = globalHInv.GetRotationMatrix();

        // Fill necessary components
        H[0][0] = rotation[0]; // x projection
        H[0][1] = rotation[1]; // y projection
        H[1][0] = rotation[3]; // x projection
        H[1][1] = rotation[4]; // y projection
        
        return H;
    }
    
    /**
     * Apply filtering to given state and measurement
     * 
     * @param ts track state
     * @param r  measurement residual
     */
    void EUTelKalmanFilter::filterTrackState( EUTelTrackStateImpl* ts, TVectorD& rkkm1 ) {
        return;       
    }
        
    /**
     * Distributes hits among measurement layers in the order
     * seen by track traversing the telescope planes 
     * 
     * @param hits Vector of hits to be assigned to the measurement layers of Klaman filter
     */
    bool EUTelKalmanFilter::sortHitsByMeasurementLayers( const EVENT::TrackerHitVec& hits ) {
        bool isReady = false;
        
        if ( hits.empty() ) {
            streamlog_out(WARNING1) << "No hits supplied this event." << std::endl;
            return isReady;
        }
        
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::sortHitsByMeasurementLayers()" << std::endl;
       
        // Clear information leftover from previous events
        std::vector< MeasurementLayer* >::iterator itLayer;
        for ( itLayer = _allMeasurements.begin(); itLayer != _allMeasurements.end(); ++itLayer) {
            delete (*itLayer);
        }
        _allMeasurements.clear();
        
        // Distribute all hits among measurement layers
        std::map< int, MeasurementLayer* > measLayers;
        std::map< int, MeasurementLayer* >::iterator itMeasLayer;
        
        int sensorID = -1;
        EVENT::TrackerHitVec::const_iterator itr;
        for( itr = hits.begin(); itr != hits.end(); ++itr ) {
            sensorID = Utility::GuessSensorID( *itr );
            itMeasLayer = measLayers.find( sensorID );
            if ( itMeasLayer != measLayers.end() ) itMeasLayer->second->addHit( *itr );
            else {
                measLayers[ sensorID ] = new MeasurementLayer( sensorID );
                measLayers[ sensorID ]->addHit( *itr );
            }
        }
        
        // Sort measurement layers such that layers encountered by track first
        // are in the front of array
        _allMeasurements = std::vector< MeasurementLayer* >( geo::gGeometry().nPlanes(), 0 );   // flush vector
        int numberAlongZ = -1;
        std::vector< int >::const_iterator itSensorID;
        for ( itSensorID = geo::gGeometry().sensorIDsVec().begin(); itSensorID != geo::gGeometry().sensorIDsVec().end(); ++itSensorID ) {
            sensorID = *itSensorID;
            itMeasLayer = measLayers.find( sensorID );
            if ( itMeasLayer != measLayers.end() ) {
                numberAlongZ = geo::gGeometry().sensorIDtoZOrder( sensorID );
                _allMeasurements.at( numberAlongZ ) = itMeasLayer->second;
            }
        }
        
        // remove elements without MeasurementLayer assigned
        for ( itLayer = _allMeasurements.begin(); itLayer != _allMeasurements.end(); ) {
            if ( !(*itLayer) ) itLayer = _allMeasurements.erase(itLayer);
            else ++itLayer;
        }
        
        isReady = true;
        return isReady;
    }
    
    MeasurementLayer::MeasurementLayer() : _id(-1), _allHits() {}
    
    MeasurementLayer::MeasurementLayer( int id ) : _id(id), _allHits() {}
    
    MeasurementLayer::~MeasurementLayer() { }
    
    void MeasurementLayer::addHit( EVENT::TrackerHit* hit ) {
        _allHits.push_back( hit );
    }
    
} // namespace eutelescope
