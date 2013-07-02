/* 
 * File:   EUTelMagneticFieldFinder.cc
 *
 * Created on July 2, 2013, 12:59 PM
 */

#include "EUTelMagneticFieldFinder.h"

namespace eutelescope {

    EUTelKalmanFilter::EUTelKalmanFilter() : EUTelTrackFitter( "KalmanFilter" ) {}

    EUTelKalmanFilter::EUTelKalmanFilter( std::string name ) : EUTelTrackFitter( name ) {}

    EUTelKalmanFilter::~EUTelKalmanFilter() {}
    
    void EUTelKalmanFilter::FitTracks() {}
}
