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
            _beamQ(-1.){}

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
            _beamQ(-1.){}

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
    }
    
    /** Check validity of the user input */
    bool EUTelKalmanFilter::initialise() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::initialise()" << std::endl;
        
        _isReady = true;
        
        // Check the validity of supplied beam energy 
        if ( _beamE <= 0. ) {
            streamlog_out(ERROR1) << "Beam direction was set incorrectly" << std::endl;
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
        
        const EVENT::TrackerHitVec& hitFirstLayer = _allMeasurements[1]->getHits();
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
        const double p  = ts->getInvP() / fabs( _beamQ );
        const double tx = ts->getTx();
        const double ty = ts->getTy();
        const double px = p*ty / sqrt( tx*tx + ty*ty +tx*tx*ty*ty );
        const double py = p*tx / sqrt( tx*tx + ty*ty +tx*tx*ty*ty );
        const double pz = p*tx*ty / sqrt( tx*tx + ty*ty +tx*tx*ty*ty );
        
        return TVector3(px,py,pz);
    }
    
    std::vector< double > EUTelKalmanFilter::findIntersection( EUTelTrackStateImpl* ts ) const {
        
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
            
            const double a = -0.5 * rho * ( planesNorm[*itPlaneId].Dot( pVecCrosH ) ) / p;
            const double b = planesNorm[*itPlaneId].Dot( pVec ) / p;
            const double c = planesNorm[*itPlaneId].Dot( delta );
            
            std::vector< double > sol = Utility::solveQuadratic(a,b,c);
            
//            TVector3 newPos1 = getXYZfromArcLenght( ts, sol[0] );
//            TVector3 newPos2 = getXYZfromArcLenght( ts, sol[1] );
//            TVector3 newPos = ( newPos1.Z() > 0 ) ? newPos1 : newPos2;
            
        }
    }
    
    TVector3 EUTelKalmanFilter::getXYZfromArcLenght( const EUTelTrackStateImpl* ts, double s ) {
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
        TVector3 pos( x0 + 1./p * pVec.X() * s,
                      y0 + 1./k * pVec.Y() * sin( rho*s ) + 1./k * pVec.Z() * ( 1. - cos( rho*s ) ),
                      z0 + 1./k * pVec.Z() * sin( rho*s ) - 1./k * pVec.Y() * ( 1. - cos( rho*s ) ) );
        
        return pos;
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
