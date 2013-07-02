/* 
 * File:   EUTelMagneticFieldFinder.h
 *
 * Created on July 2, 2013, 12:53 PM
 */

#ifndef EUTELMAGNETICFIELDFINDER_H
#define	EUTELMAGNETICFIELDFINDER_H

// system includes <>
#include <string>
#include <map>

// EUTELESCOPE
#include "EUTelUtility.h"
#include "EUTelTrackFitter.h"

namespace eutelescope {
    
    class EUTelKalmanFilter : public EUTelTrackFitter {
        
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelKalmanFilter)        // prevent users from making (default) copies of processors
                
    public:
        EUTelKalmanFilter( );
        
        explicit EUTelKalmanFilter( std::string name );

        virtual ~EUTelKalmanFilter();
    
        /** Fit supplied hits */
        virtual void FitTracks();
    };
}

#endif	/* EUTELMAGNETICFIELDFINDER_H */

