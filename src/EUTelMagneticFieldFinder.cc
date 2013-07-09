/* 
 * File:   EUTelMagneticFieldFinder.cc
 *
 * Created on July 2, 2013, 12:59 PM
 */

#include "EUTelMagneticFieldFinder.h"

#include "streamlog/streamlog.h"

#include <iostream>
#include <functional>
#include <algorithm>

namespace eutelescope {

    EUTelKalmanFilter::EUTelKalmanFilter() : EUTelTrackFitter( "KalmanFilter" ), 
            _tracks(), 
            _trackStates(), 
            _allHits(),
            _isReady(false),
            _allowedMissingHits(0),
            _maxTrackCandidates(0),
            _beamDir(-1.,0.,0.,0.){}

    EUTelKalmanFilter::EUTelKalmanFilter( std::string name ) : EUTelTrackFitter( name ),
            _tracks(),
            _trackStates(), 
            _allHits(),
            _isReady(false),
            _allowedMissingHits(0),
            _maxTrackCandidates(0),
            _beamDir(-1.,0.,0.,0.){}

    EUTelKalmanFilter::~EUTelKalmanFilter() {}
    
    /** Perform Kalman filter track search and track fit */
    void EUTelKalmanFilter::FitTracks() {
        
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
        _isReady = true;
        
        // Check the validity of supplied beam energy 
        const double beamEnergy = _beamDir.Energy();
        if ( beamEnergy <= 0. ) {
            streamlog_out(ERROR1) << "Beam direction was set incorrectly" << std::endl;
            _isReady = false;
            return _isReady;
        }
        
        return _isReady;
    }
    
    /** Generate seed track states necessary to
     * start Kalman filter
     * 
     * Generate as many starting states as number of hits in first telescope plane
     * plus one additional state for missing hit in the fist plane if allowed
     * 
     */
    void EUTelKalmanFilter::initialiseSeeds() {

        // The object will be owned by LCIO collection. Must not to free.
        IMPL::TrackStateImpl* state = new IMPL::TrackStateImpl;
        
        // Use beam direction as a seed track state
        state->setLocation( IMPL::TrackStateImpl::AtFirstHit );
//        state->
        
        _trackStates.push_back( state );
    }
      
    /**
     * Distributes hits among measurement layers in the order
     * seen by track traversing the telescope planes 
     * 
     * @param hits Vector of hits to be assigned to the measurement layers of Klaman filter
     */
    void EUTelKalmanFilter::sortHitsByMeasurementLayers( const EVENT::TrackerHitVec& hits ) {
       
        // Distribute all hits among measurement layers
        std::map< int, MeasurementLayer* > measLayers;
        std::map< int, MeasurementLayer* >::iterator itMeasLayer;
        
        int sensorID = -1;
        EVENT::TrackerHitVec::const_iterator itr;
        for( itr = hits.begin(); itr != hits.end(); ++itr ) {
            sensorID = Utility::GuessSensorID( *itr );
            itMeasLayer = measLayers.find( sensorID );
            if ( itMeasLayer != measLayers.end() ) itMeasLayer->second->addHit( *itr );
            else measLayers[ sensorID ] = new MeasurementLayer( sensorID );
        }
        
        // Sort measurement layers such that layers encountered by a track first are
        // in the front of array
        
    }
    
    MeasurementLayer::MeasurementLayer() : _id(-1), _allHits() {}
    
    MeasurementLayer::MeasurementLayer( int id ) : _id(id), _allHits() {}
    
    MeasurementLayer::~MeasurementLayer() {}
    
    void MeasurementLayer::addHit( EVENT::TrackerHit* hit ) {
        _allHits.push_back( hit );
    }
    
}
