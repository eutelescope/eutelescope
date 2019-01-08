/*
 * Created by Eda Yildirim
 *  (2015 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 * Updated by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef EUTELMISSINGCOORDINATEESTIMATOR_H
#define EUTELMISSINGCOORDINATEESTIMATOR_H

// eutelescope includes ".h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <set>

namespace eutelescope
{

    //! Missing Coordinate Estimator
    /*! As the name suggests this processor finds the estimated
     *  position of the missing coordinate on your strip sensor.
     *  How it works is: it gets the hits from two specified reference sensors,
     *  finds the closest hit pairs, make a straight line and finds
     *  the estimated position in one axis on your sensor.
     *  No promises that this will work with tilted sensors and/or with magnetic field
     *  One needs to used this with merged hits and after pre-alignment
     */

    class EUTelMissingCoordinateEstimator : public marlin::Processor
    {

	private:
	    DISALLOW_COPY_AND_ASSIGN ( EUTelMissingCoordinateEstimator )

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new EUTelMissingCoordinateEstimator;
	    }

	    EUTelMissingCoordinateEstimator ();

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void end ( );

	    void bookHistos ( );

	protected:

	    std::string _inputHitCollectionName;

	    std::string _outputHitCollectionName;

	    EVENT::IntVec _referencePlanes;

	    EVENT::IntVec _dutPlanes;

	    std::string _missingCoordinate;

	    float _maxResidual;

	    TrackerHitImpl* cloneHit ( TrackerHitImpl *inputHit );

	private:

	    bool _multihitmode;

	    unsigned int _missingHitPos;

	    unsigned int _knownHitPos;

	    unsigned int _nDutHits;

	    unsigned int _nDutHitsCreated;

	    unsigned int _maxExpectedCreatedHitPerDUTHit;

	    unsigned int _nNoReferenceCount;

	    unsigned int _nResidualFailCount;

	    unsigned int _numberOfCreatedHitsPerDUTHit[10];
    };

    EUTelMissingCoordinateEstimator gEUTelMissingCoordinateEstimator;

}

#endif
