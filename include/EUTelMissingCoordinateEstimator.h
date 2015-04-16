/*
 * Created by Eda Yildirim
 *  (2015 DESY)
 *
 *  email:eda.yildirim@cern.ch
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

namespace eutelescope {
    
    //! Missing Coordinate Estimator
    /*! As the name suggest this processor is finds the estimated
     *  position of the missing coordinate on your sensor
     *  How it works is simple, it gets the hits from specified two sensors
     *  finds the closest hit pairs, make a straight line out of it and find
     *  the estimated position in one axis on your sensor you want.
     *  No promises that this will work with tilted sensors and/or with magnetic field
     *  One needs to used this with merged hits and after pre-alignment
     */
    
    class EUTelMissingCoordinateEstimator : public marlin::Processor {
        
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelMissingCoordinateEstimator)
        
    public:
        
        
        //! Returns a new instance of EUTelMissingCoordinateEstimator
        /*! This method returns a new instance of this processor.  It is
         *  called by Marlin execution framework and it shouldn't be
         *  called/used by the final user.
         *
         *  @return a new EUTelMissingCoordinateEstimator.
         */
        virtual Processor * newProcessor() {
            return new EUTelMissingCoordinateEstimator;
        }
        
        //! Default constructor
        EUTelMissingCoordinateEstimator ();
        
        //! Called at the job beginning.
        /*! This is executed only once in the whole execution. It prints
         *  out the processor parameters and check that the GEAR
         *  environment is properly set up and accessible from Marlin.
         */
        virtual void init ();
        
        //! Called for every run.
        /*! It is called for every run, and consequently the run counter
         *  is incremented. The geometry ID of the file is compared with
         *  the one provided by the GEAR geometry description. In case the
         *  two are different, the user is asked to decide to quit or to
         *  continue with a description that might be wrong.
         *
         *  @param run the LCRunHeader of the this current run
         */
        virtual void processRunHeader (LCRunHeader * run);
        
        //! Called every event
        /*! This is called for each event in the file. Each element of the
         *  pulse collection is scanned and the center of the cluster is
         *  translated into the external frame of reference thanks to the
         *  GEAR geometry description.
         *
         *  The cluster center might be calculate using a standard linear
         *  charge center of gravity algortihm or applying a more
         *  sophisticated non linear eta function. This behaviour is
         *  regulated by the user from the steering file.
         *
         *  @throw UnknownDataTypeException if the cluster type is unknown
         *
         *  @param evt the current LCEvent event as passed by the
         *  ProcessMgr
         */
        virtual void processEvent (LCEvent * evt);
        
        
        //! Called after data processing.
        /*! This method is called when the loop on events is
         *  finished.
         */
        virtual void end();
        
        
        //! Histogram booking
        /*! Does nothing
         */
        void bookHistos();
        
        
    protected:
        
        //! Input TrackerHit collection name
        /*! This is the name the user wants to give to the input hit
         *  collection.
         */
        std::string _inputHitCollectionName;
        
        //! Output TrackerHit collection name
        /*! This is the name the user wants to give to the output hit
         *  collection.
         */
        std::string _outputHitCollectionName;
        
        //! Reference planes
        /*! This is the list of sensorIDs that their hits will be used 
         *  to estimate the missing coordinate on your DUT. You have to 
         *  give exactly 2 sensorIDs. For better results use the ones 
         *  that are closest to your DUT
         */
        EVENT::IntVec _referencePlanes;

        //! DUT Planes
        /*! This is the list of sensorIDs that missing coordinate of their 
         *  hits needs to be found. Notice that if the specified coordinate 
         *  already exists it will be overwritten
         */
        EVENT::IntVec _dutPlanes;
        
        
        //! Missing Coordinate
        /*! The coordinate axis that needs to be estimated. 
         *  You have to set this to either X or Y.
         */
        std::string _missingCoordinate;
        
        //! Max Residual
        /*! This processor will look for a closest hits (in known coordinate)
         *  to determine if the hits are correlated
         *  The hits will be considered as correlated if the residual is smaller
         *  than MaxResidual
         */
        float _maxResidual;
        
        //! Clone Hit
        /*! This method is used to clone TrackerHitImpl object
         */
        TrackerHitImpl* cloneHit(TrackerHitImpl *inputHit);
        
        
    private:
        
        //! Run number
        int _iRun;
        
        //! Event number
        int _iEvt;
        
        //! Missing hit position
        unsigned int _missingHitPos;

        //! Known hit position
        unsigned int _knownHitPos;
        
        //! Number of DUT hits
        unsigned int _nDutHits;
        
        //! Number of Dut hits created
        unsigned int _nDutHitsCreated;
	       
 	//! Number of Expected created hit per DUT hit
 	unsigned int _maxExpectedCreatedHitPerDUTHit;
	
	//! Count number of created hit per DUT Hit
	std::vector<unsigned int> _numberOfCreatedHitPerDUTHit; 
    };
    
    //! A global instance of the processor
    EUTelMissingCoordinateEstimator gEUTelMissingCoordinateEstimator;
    
}
#endif
