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

//LCIO
#include "IMPL/TrackImpl.h"

#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <math.h>

namespace eutelescope {

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
						_AllowedSharedHitsOnTrackCandidate(0),
            _maxTrackCandidates(0),
            _beamE(-1.),
            _beamQ(-1.),
            _beamEnergyUncertainty(0.),
            _beamAngularSpread(2,-1.),
	    _jacobianF(5,5),
            _trkParamCovCkkm1(5,5),
            _processNoiseQ(5,5),
            _residualCovR(2,2),
            _eomIntegrator( new EUTelUtilityRungeKutta() ),
            _jacobianIntegrator( new EUTelUtilityRungeKutta() ),
            _eomODE( 0 ),
            _jacobianODE( 0 ),
            _planesForPR(1) 
            {
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
            _residualCovR(2,2),
            _eomIntegrator( new EUTelUtilityRungeKutta() ),
            _jacobianIntegrator( new EUTelUtilityRungeKutta() ),
            _eomODE( 0 ),
            _jacobianODE( 0 ),
            _planesForPR( 1 ) 
            {
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
 
    /** */
    int EUTelKalmanFilter::findNextPlaneEntrance(  EUTelTrackStateImpl* ts , int sensorID  ){

// how to make sure that ts is in localframe system??

      const double x0 = ts->getReferencePoint()[0];
      const double y0 = ts->getReferencePoint()[1];
      const double z0 = ts->getReferencePoint()[2];
      int tsPlaneID   = ts->getLocation();
	

      // get the very first sensor : why?
      // remember last hit-point from the track candidate below
      double lpoint[3]   = { x0, y0, z0}; // 
      double gpoint[3]   = {0.,0.,0.};  // initialise output point-vector global frame (World)

      geo::gGeometry().local2Master( tsPlaneID, lpoint, gpoint );
      streamlog_out ( DEBUG2 ) << "findNextPlaneEntrance: " << tsPlaneID << " " <<  lpoint[0] << " " <<  lpoint[1] << " " << lpoint[2] <<  " :"
                                << gpoint[0] << " "    << gpoint[1] << " "   << gpoint[2] << " " << endl; 

      gpoint[2] += -100.; 
//      float start[3]     = {static_cast<float> (gpoint[0]),static_cast<float> (gpoint[1]),static_cast<float> (gpoint[2]) };  

      double dir[3]      = {0.,0.,1.};  // as good as any other value along z axis.

      float nextPoint[3] = {0.,0.,0.};  // initialise output point-vector local frame (sensor)

            int newSensorID   = geo::gGeometry( ).findNextPlaneEntrance( gpoint,  dir, sensorID, nextPoint ) ;

// get to local frame for sensor newSensorID:
            double dGpoint[3] = {static_cast<double> (nextPoint[0]),static_cast<double> (nextPoint[1]),static_cast<double> (nextPoint[2]) };
            double dLpoint[3] = {0.,0.,0.};  // initialise output point-vector local frame (sensor)
            geo::gGeometry().master2Local( dGpoint, dLpoint);

            if( newSensorID < 0 ) 
            {
              streamlog_out ( DEBUG4 ) << "no Entrance: " <<  lpoint[0] << " " <<  lpoint[1] << " " << lpoint[2] << " err:"<< newSensorID << endl;
            } else {     
              streamlog_out ( DEBUG2 ) << "identified NextPlane Entrance: " <<  lpoint[0] << " " <<  lpoint[1] << " " << lpoint[2] <<  " at : "<< newSensorID << endl;
              const float opoint[3] = { lpoint[0], lpoint[1], lpoint[2] };
              ts->setReferencePoint( opoint );  
              ts->setLocation( newSensorID );
           }

      return newSensorID;
    }

    /** */



    void EUTelKalmanFilter::propagateFromRefPoint( 	std::vector< EUTelTrackImpl* >::iterator &itTrk    ){
			
	//Set up the geometry//////////////////////////////////////////////////////////////////////// 
      	const map< int, int > sensorMap = geo::gGeometry().sensorZOrdertoIDs();
      	int planeID     = sensorMap.at(0); // the first first plane in the array of the planes according to z direction. // assume not tilted plane. 
	//      const int    iPlane             = geo::gGeometry().sensorIDtoZOrder(planeID);
      	const double thicknessSen       = geo::gGeometry().siPlaneZSize( planeID );
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

      	EUTelTrackStateImpl* state = const_cast<EUTelTrackStateImpl*>((*itTrk)->getFirstTrackState( ));
	state->setbeamQ(_beamQ); //Set the beam charge here. This is not perfect I think. Since we could set it as a static variable. However how this should be used in other processor I am unsure???????/
	streamlog_out ( DEBUG0 ) << "Memory location of initial state: " << state << endl; 

      	if( state == 0 )
      	{
        	streamlog_out ( WARNING0 ) << "track _tracksCartesian return a NULL state, skip this one " << endl; 
        	return ;
      	}
	
      	// loop through all known sensors (local, defined above):
      	map< int, int >::const_iterator iter = sensorMap.begin();
      	while ( iter->first < sensorMap.size()-1 ) { // do not iterate through the very last plane
        if( iter->first < 0 ) continue;
	
        float dpoint[3];
        int next = iter->first + 1;
        map< int, int >::const_iterator nextSensor = sensorMap.find( next );
        int newSensorID = state->findIntersectionWithCertainID( nextSensor->second, dpoint );
        int sensorIntersection = geo::gGeometry( ).getSensorID(dpoint);
        streamlog_out ( DEBUG5 ) << " propagateFromRefPoint: pos:" <<  iter->first << " from ID= " << iter->second << " to ID=" <<  nextSensor->second << " and got: " << newSensorID << " sensorIntersection= " << sensorIntersection << std::endl;

	if( newSensorID < 0 or sensorIntersection < 0 ){ //If there was no intersection of infinite plane or if the intersection is not in the sensor.
		streamlog_out ( DEBUG4 ) << "Point (" <<  dpoint[0] << ", " <<  dpoint[1] << ", " << dpoint[2]  << ")" << "found on no sensor. New sensor ID (Was there intersection?):"  << newSensorID <<"SensorID (Was the intersection on the plane?)"<< sensorIntersection  << std::endl;
	}
	else{ 
		//So we found a intersection so we need to determine that new state.
		streamlog_out ( DEBUG3 ) << "Intersection on a plane " << newSensorID << std::endl;
		EUTelTrackStateImpl* state_new =  new EUTelTrackStateImpl(); //This creates a new state object
		state_new->setHit(NULL); //By default it should contain no hit
		state_new->setbeamQ(_beamQ); //Set the beam charge here. This is not perfect I think. Since we could set it as a static variable. However how this should be used in other processor I am unsure???????/
		
		//Here we fill the state with its new approximate new hit position. Nothing else is filled yet since this will depend on if hit information is there.
		TMatrixD jacobian(5,5);
		jacobian.Zero();
		jacobian = state->getPropagationJacobianF((dpoint[2] - state->getZParameter())); //Find all the relations between state variables at a particular z parameter dpoint[2] 
        	
  		//jacobian.Print();
		nextStateUsingJacobianFinder(state, state_new, jacobian); //Here we determine the new state position and CovMatrix using the jacobian. This might not need to be done now but would involve changing closestHit()????
		state_new->setZParameter( dpoint[2] ); //Set this here since it is not a state variable but a parameter
		state_new->setLocation( newSensorID );

		double global[] = { state_new->getX(),state_new->getY(),state_new->getZParameter() };
		double local[3];
		geo::gGeometry().master2Localtwo( state_new->getLocation(), global, local );
		const float localRef[3] = {local[0], local[1], local[2]}; 
		state_new->setReferencePoint(localRef);
			
		streamlog_out ( DEBUG5 ) << "Both FindIntersection and Jacobian should be the same" << std::endl;
		streamlog_out ( DEBUG5 ) << "Point (" <<  dpoint[0] << ", " <<  dpoint[1] << ") from findIntersection" << std::endl;
		streamlog_out ( DEBUG5 ) << "Point (" <<  state_new->getX() << ", " <<  state_new->getY() << ") from jacobian to sensor : " << state_new->getLocation()  << std::endl<<std::endl;
						
		////////////Find next closest hit and determine if it is within window. If both fill new state with information about hit. Otherwise fill without it/////////////////////////////// 
         	EVENT::TrackerHit* closestHit = const_cast< EVENT::TrackerHit* > ( findClosestHit( state_new, newSensorID ) ); //This will look for the closest hit but not if it is within the excepted range		
 		if ( closestHit ){ //Just check if the closestHit exist 
         	
			const double* uvpos = closestHit->getPosition(); //Get that hits position
        		const double distance = getResidual( state_new, closestHit ).Norm2Sqr( ); //Determine the residual to it. //Distance is in mm.
            		const double DCA = 10000;  //DCA since is mm^2 since distance is.                             
			//getXYPredictionPrecision( state_new ); // basically RMax cut at the moment. MUST FIX!!!! HOW SHOULD THIS BE TREAT?;;
	           	streamlog_out ( DEBUG1 ) << "NextPlane " << newSensorID << " " << uvpos[0] << " " << uvpos[1] << " " << uvpos[2] << " resid:" << distance << " ResidualCut: " << DCA << endl;
           		if ( distance > DCA ) {
             	           	streamlog_out ( DEBUG1 ) << "Closest hit is outside of search window." << std::endl;
              	           	streamlog_out ( DEBUG1 ) << "Skipping current plane. Covariant Matrix and position already updated to this point " << std::endl;
	  	  	   	streamlog_out ( DEBUG1 ) << " Distance between them: "<< distance << endl;
            		}
			else{
            			streamlog_out (DEBUG1 ) << "NextPlane MATCHED. Position of Hit (Local) " <<  uvpos[0] << " " <<  uvpos[1] << " " << uvpos[2] 
                                                        <<" Position of state (Global) " << state_new->getX() 
                                                        << "," << state_new->getY()
                                                        <<","<< state_new->getZParameter() 
                                                        <<" Distance between them: " << distance 
                                                        << " Sensor ID:" << state_new->getLocation() 
                                                        << " Seed we began at: " << (*itTrk) << endl;
	              		streamlog_out ( DEBUG1 ) << "Will now alter Cov matrix and state variables using hit information " << std::endl;
				TMatrixD HMatrix = state_new->getH(); //We need to be able to move from the measurement to the state space
				TMatrixD GainMatrix = updateGainK( state_new, closestHit ); //This is a matrix that tells you how much the state should be changed with the information from the hit
				UpdateStateUsingHitInformation( state_new ,closestHit, jacobian, GainMatrix, HMatrix); //Update the state on the track
				UpdateTrackUsingHitInformation( state_new , closestHit, (*itTrk), jacobian, GainMatrix,HMatrix); //Update the track itself.									
            		}

		}
                if( streamlog_level(DEBUG1) ) state_new->Print();
  		(*itTrk)->addTrackState( new EUTelTrackStateImpl(*state_new) ); //New memory allocation here should be deleted by LCIO memory management. I think...
		streamlog_out ( DEBUG1 ) << "Memory location of initial state after all allocation (Should be the same): " << state << endl; 
		streamlog_out ( DEBUG1 ) << "Memory location of state_new: " << state_new << endl; 
						
	 	state = state_new; //MUST DOUBLE CHECK THIS MEMORY MANAGEMENT DOES NOT LOOK RIGHT TO ME
		streamlog_out ( DEBUG1 ) << "Memory location of initial state after it was made equal to state_new: " << state << endl; 		

          }
	  ++iter;	
       
 
    }

}

    //Print the list of tracks given in _collection
    void EUTelKalmanFilter::Print( std::string Name, std::vector< EUTelTrackImpl*> & _collection) 
    {
       int itrk = 0; 
       int size_itTrk  = _collection.size();

       std::vector< EUTelTrackImpl* >::iterator itTrk;
 
       for ( itTrk = _collection.begin(); itTrk != _collection.end(); itTrk++, itrk++ ) 
       {
          if( (*itTrk) == 0 ) 
          {         
             streamlog_out(WARNING1) << " Track vector: " << Name.c_str() << " contains no information at track = " << itrk << "." << std::endl;
             continue; 
          }
          const EVENT::TrackerHitVec& ihits = (*itTrk)->getTrackerHits();
          int nhits =  ihits.size( ) ;
          int expec =  geo::gGeometry( ).nPlanes( ) - _allowedMissingHits;
          streamlog_out(DEBUG5) << " Track vector " << Name.c_str() << " at track number " <<  itrk << " of size " << size_itTrk << " with " << nhits << " and expecting at least " << expec << " hits " << std::endl;
          (*itTrk)->Print();
       }
    }
 
    void EUTelKalmanFilter::Prune( std::vector< EUTelTrackImpl*> &_collection, std::vector< EUTelTrackImpl*> &_collection_to_delete ) 
    {
      int itrk = 0 ;
      std::vector< EUTelTrackImpl* >::iterator itTrk;

      for (  itTrk = _collection.begin(); itTrk != _collection.end();) 
      {
        // obsolte code from very old times ....
        // check that it's not a fitted hit (type =32)
	//     if( (*itTrk)->getType() > 31 ) continue;            Is this the best way to deal with type? Should it no be in LCIO bit field?
 
        bool iend = std::find( _collection_to_delete.begin(), _collection_to_delete.end(), (*itTrk) ) == _collection_to_delete.end(); 
        if( iend )
        {
           streamlog_out(DEBUG3) << " track " << itrk << " at " << (*itTrk) << " NOT found  in _collection_to_delete " << std::endl;
itTrk++;
	   itrk++ ;
        } else {
           streamlog_out(DEBUG3) << " track " << itrk << " at " << (*itTrk) << "     found  in _collection_to_delete, deleting ... " << std::endl;
           delete (*itTrk);
           itTrk = _collection.erase(itTrk);
        } 
      }
    }


    /** Perform Kalman filter track search and track fit */
    void EUTelKalmanFilter::SearchTrackCandidates() {

      streamlog_out(MESSAGE1) << "EUTelKalmanFilter::SearchTrackCandidates()" << std::endl;

      // Check if the fitter was configured correctly
        if ( !_isReady ) {
            streamlog_out(ERROR1) << _name << ": Can't fit. The fitter is not initialised!" << std::endl;
            return;
        }
        
        // Flush data from previous event
        reset();
        
        // Initialise Kalman filter states
        initialiseSeeds();
        
	// Start Kalman filter
	std::vector< EUTelTrackImpl* >::iterator itTrk;
        int local_itTrk = 0;
        int size_itTrk = _tracksCartesian.size();
				//LOOP over list of all tracks. At this point they are only the hits on the first and second plane as track objects.
        for ( itTrk = _tracksCartesian.begin(); itTrk != _tracksCartesian.end(); itTrk++ ) {
            streamlog_out(MESSAGE0) << "beginning now at : " << local_itTrk << " of " << size_itTrk << " itTrk: " << (*itTrk) << endl;
            bool isGoodTrack = true; 
         
             streamlog_out ( DEBUG5 ) << "The track seed before propagation to find other states and hits on other surfaces. " << endl; 
            (*itTrk)->Print();

						//Get the first state of the track. Remember at this point we only have single hits on a track
            EUTelTrackStateImpl* state = const_cast<EUTelTrackStateImpl*>((*itTrk)->getFirstTrackState( ));
            
						//Check is the state of the hit on the track exists
            if( state == 0 )
            {
               streamlog_out ( WARNING0 ) << "track _tracksCartesian return a NULL state, skip this one " << endl; 
               isGoodTrack = false;
               delete (*itTrk); //Is this really needed?
               itTrk = _tracksCartesian.erase( itTrk ); itTrk--; //Must itTrk-- since the erase will take us to the next hit automatically 
               continue;
            }
						//This will propagate the initial seed to the other planes. On the new plane it will add hits and change the state depending on relation to hit.
            propagateFromRefPoint(itTrk);
 
						//Check the number of hit on the track after propagation and collecting hits is over the minimum
            if ( isGoodTrack && ( *itTrk )->getHitsOnTrack().size( ) < geo::gGeometry( ).nPlanes( ) - _allowedMissingHits ) {
		streamlog_out ( MESSAGE2 ) << "Number of hits on the track: " << ( *itTrk )->getTrackerHits( ).size( ) <<" Number needed: " << geo::gGeometry( ).nPlanes( ) - _allowedMissingHits << std::endl;
               	streamlog_out ( DEBUG5 ) << "Track candidate has to many missing hits." << std::endl;
               	streamlog_out ( DEBUG5 ) << "Removing this track candidate from further consideration." << std::endl;
           	(*itTrk)->Print();
               	delete (*itTrk);  //Is this really needed?
               	isGoodTrack = false;
               	itTrk = _tracksCartesian.erase( itTrk ); itTrk--; 
            }

            if ( isGoodTrack ) {
//               state->setLocation( EUTelTrackStateImpl::AtLastHit );
               	int nstates = (*itTrk)->getTrackStates().size();
	
               	streamlog_out(DEBUG5) << "'Tracks' looped through. I.e initial hit seed as a track. (Ignore states with no hits): "  << local_itTrk << ". Number of seeds: " << size_itTrk << ". At seed: " << *itTrk << " after propagation with: " << ( *itTrk )->getTrackerHits( ).size( ) << " hits collected on track.  Expecting at least " << geo::gGeometry( ).nPlanes( ) - _allowedMissingHits << " The states. I.e Including planes with no hits: "<< nstates << std::endl << std::endl << std::endl << std::endl;

             	streamlog_out ( DEBUG5 ) << "Successful track state after propagation: " << endl; 
            	(*itTrk)->Print();
            }

            local_itTrk++;
        }

        streamlog_out(MESSAGE0) << "Finished looping through all seeds! Looped at total of : " << local_itTrk << " after seeds with no state are deduced. Total seeds were: " << size_itTrk << endl;
  
       streamlog_out(MESSAGE1) << "------------------------------EUTelKalmanFilter::SearchTrackCandidates()---------------------------------" << std::endl;

    }
 

    // Perform track pruning this removes tracks that have the same hits used to create the track on some planes
    void EUTelKalmanFilter::PruneTrackCandidates() {

      	streamlog_out(MESSAGE1) << "EUTelKalmanFilter::PruneTrackCandidates()" << std::endl;
 
	//   std::vector< IMPL::TrackImpl* >::iterator itTrk;
	std::vector< EUTelTrackImpl* >::iterator itTrk;
      	std::vector< EUTelTrackImpl* >::iterator jtTrk;
      	std::vector< EUTelTrackImpl* >  _tracks_to_delete;

      	int itrk = 0;
      	int jtrk = 0;

	//Loop over all tracks 
	for ( itTrk = _tracksCartesian.begin(); itTrk != _tracksCartesian.end(); itTrk++, itrk++ ) 
	{
        	streamlog_out(MESSAGE1) <<  "Loop at track number:" <<  itrk << " at  " << ( *itTrk) << std::endl;
   
	        const EVENT::TrackerHitVec ihits = (*itTrk)->getTrackerHits(); //get the hits contained within this track object
 
		//Now loop through all tracks one ahead of the original track itTrk. This is done since we want to compare all the track to each other to if they have similar hits     
	        for ( jtrk = itrk+1; jtrk < _tracksCartesian.size(); jtrk++ ) 
	        {
          		int hitscount=0;
 
          		EUTelTrackImpl* jtTrack = _tracksCartesian[jtrk];
          		const EVENT::TrackerHitVec jhits = jtTrack->getTrackerHits();
			// cross check every track to all following ones in the _track collection
			int outer_loop = 0;
			int inner_loop = 0;
		        for(int i=0;i<ihits.size();i++)
          		{ // loop through hits in itTrk candidate
				streamlog_out(MESSAGE1) <<  "The number of hits to loop through (outer) :" << ihits.size()<< " Beginning loop (outer) : " <<outer_loop<< std::endl;
				outer_loop++; 
            			EVENT::TrackerHit* ihit = ihits[i];
            			int ic = ihit->id();
				inner_loop = 0;
            
				for(int j=0;j<jhits.size();j++)
            			{ // loop through hits in jtTrk candidate 
					streamlog_out(MESSAGE1) <<  "The number of hits to loop through (compare, inner)  :" << jhits.size()<< " Beginning loop (inner) : " <<inner_loop<< std::endl; 
					inner_loop++;
               				EVENT::TrackerHit* jhit = jhits[j];
               				int jc = jhit->id();
               				if(ic == jc ){
                 				hitscount++; 
                 				streamlog_out(MESSAGE1) <<  "Hit number on track you are comparing all other to :" << i << ". Hit ID: " << ic << ". Hit number of comparison: " << j << ". Hit ID of this comparison : " << jc << ". Number of common hits: " << hitscount << std::endl; 
               				}
            			}
          		} 
    			//Here we fill up a vector with tracks we want to delete.
          		if(hitscount > _AllowedSharedHitsOnTrackCandidate) {   
            			if( std::find(_tracks_to_delete.begin(), _tracks_to_delete.end(), *itTrk ) == _tracks_to_delete.end() || _tracks_to_delete.size() == 0 ) {
               				streamlog_out(DEBUG5) <<  "Track " << itrk << " and " << jtrk <<"Have similar hits. Remove track " <<*itTrk << std::endl;
              				_tracks_to_delete.push_back( *itTrk );
            			}
            			continue;
          		}
        	}

      	}

      	Print( "_tracksCartesian ", _tracksCartesian);   

      	Print( "Tracks that are to be deleted: ", _tracks_to_delete);   
 
      	Prune( _tracksCartesian, _tracks_to_delete);   

      	Print( "After deletion the tracks left: ", _tracksCartesian);   
 
      	streamlog_out(MESSAGE1) << "------------------------------EUTelKalmanFilter::PruneTrackCandidates()---------------------------------" << std::endl;
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
            streamlog_out(WARNING1) << "Collection of hits is empty. Check supplied hits" << std::endl;
            _isReady = false;
            return _isReady;
        }
       
        streamlog_out(DEBUG4) << "Initialise fitted hit vector, size = " << hitFittedVec.size() << ", resetting ..." << std::endl;
        hitFittedVec.clear();
        streamlog_out(DEBUG4) << "Initialise fitted hit vector, size = " << hitFittedVec.size()  << std::endl;
                
     
        streamlog_out(DEBUG2) << "Initialisation successfully completed" << std::endl;
        streamlog_out(DEBUG2) << "------------------------------EUTelKalmanFilter::initialise()------------------------------" << std::endl;
        return _isReady;
    }
    
    void EUTelKalmanFilter::reset() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::reset()" << std::endl;
        _tracksCartesian.clear();
        _trackStates.clear();
        streamlog_out(DEBUG2) << "-------------------------------EUTelKalmanFilter::reset()-----------------------------------" << std::endl;
    }
    
    void EUTelKalmanFilter::setHits( EVENT::TrackerHitVec& hits ) {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::setHits()" << std::endl;
        if ( hits.empty() ) {
            _isHitsOK = false;
            return;
        }
        _allHits = hits;
        _isHitsOK = sortHitsByMeasurementLayers(_allHits);
        streamlog_out(DEBUG2) << "-------------------------------EUTelKalmanFilter::setHits()------------------------------------" << std::endl;
    }
 

    /** Prune seed track states necessary to
     * start Kalman filter
     * 
     * 
     */
    void EUTelKalmanFilter::pruneSeeds() {
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::pruneSeeds()" << std::endl;

        if ( _allMeasurements.empty() ) {
            streamlog_out(WARNING1) << "Can't initialise track seeds for the finder. No hits in this event." << std::endl;
            return;
        }


	// Start Kalman filter
	std::vector< EUTelTrackImpl* >::iterator itTrk;
	for ( itTrk = _tracksCartesian.begin(); itTrk != _tracksCartesian.end(); itTrk++ ) {
            bool isDuplTrack = false; 
             
            std::vector< EUTelTrackImpl* >::iterator itTrk_in;
	    for ( itTrk_in = _tracksCartesian.begin(); itTrk_in != _tracksCartesian.end(); itTrk_in++ ) {
              bool isDuplTrack = false; 

            }
            
        }
 
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

//
// loop over the first 2 planes assuming the numbering to go 0,1,2,3,.. etc.
// building track seeds from every hit found in the planes. could do from all planes, but perhaps an overkill - parameter control?
//          
        int maxPlanes =  std::min( static_cast<int> (_allMeasurements.size()) , static_cast<int> (_planesForPR) );
        for( int iplane = 0; iplane < maxPlanes ; iplane++) {

          streamlog_out(MESSAGE1) << "_allMeasurments plane : " << iplane << " of " << maxPlanes  << std::endl;
      
        const EVENT::TrackerHitVec& hitFirstLayer = _allMeasurements[iplane]->getHits();
        streamlog_out(DEBUG1) << "N hits in first non-empty layer: " << hitFirstLayer.size() << std::endl;
        EVENT::TrackerHitVec::const_iterator itHit;
        for ( itHit = hitFirstLayer.begin(); itHit != hitFirstLayer.end(); ++itHit ) {
            // The object will be owned by LCIO collection. Must not to free.
            EUTelTrackStateImpl* state = new EUTelTrackStateImpl;
//
//            // Use beam direction as a seed track state
            const double* uvpos = (*itHit)->getPosition();
            float posLocal[] =  { static_cast<float>(uvpos[0]), static_cast<float>(uvpos[1]), static_cast<float>(uvpos[2]) };

            const int sensorID = Utility::getSensorIDfromHit( *itHit );

            double temp[] = {0.,0.,0.};
            geo::gGeometry().local2Master( sensorID, uvpos, temp);
            float posGlobal[] = { static_cast<float>(temp[0]), static_cast<float>(temp[1]), static_cast<float>(temp[2]) };
 
            streamlog_out ( DEBUG2 ) << "pick next hit" << sensorID << " " <<  posLocal[0] << " " << posLocal[1] << " " << posLocal[2] <<  " :"
                                << posGlobal[0] << " "    << posGlobal[1] << " "   << posGlobal[2] << " " << endl; 

                       
//            gear::Vector3D vectorGlobal( temp[0], temp[1], temp[2] );
            const double q          = _beamQ;      // assume electron beam
            //const double invp       = fabs(q)/_beamE;
            const double invp       = q/_beamE;
            
            // Fill track parameters covariance matrix
            float trkCov[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            
            const EVENT::FloatVec uvcov = (*itHit)->getCovMatrix();
            // initial X,Y covariances are defined by seeding hit uncertainty multiplied by a huge factor
            trkCov[0] = uvcov[0]*1.E4;                               //cov(x,x)
            trkCov[1] = uvcov[1]*1.E4;   trkCov[2] = uvcov[2]*1.E4;       //cov(y,x), cov(y,y)
            // TX,TY covariances are defined by beam spread at the first plane
            
            if ( _beamAngularSpread[0] > 0. ) trkCov[5] = _beamAngularSpread[0] * _beamAngularSpread[0];          //cov(tx,x)=0, cov(tx,y)=0, cov(tx,tx)
            else trkCov[5] = 1.E5;
            if ( _beamAngularSpread[1] > 0. ) trkCov[9] = _beamAngularSpread[1] * _beamAngularSpread[1];          //cov(ty,x)=0, cov(ty,y)=0, cov(ty,tx)=0, cov(ty,ty)
            else trkCov[9] = 1.E5;
            trkCov[14] = ( _beamEnergyUncertainty * _beamE ) * ( _beamEnergyUncertainty * _beamE );               //cov(q/p,x)=0, cov(q/p,y)=0, cov(q/p,tx)=0, cov(q/p,ty)=0, cov(q/p,q/p)
            
						//NOTE THESE PARAMETERS ARE WHAT DEFINE THE MOVEMENT OF A PARTICLE THROUGH A SURFACE. 
						//NOTE THIS IS GLOBAL TELESCOPE FRAME. Also note that there is no SetZ. Since this state will be moved along the z axis by the jacobian calculated later. So it is only a parameter
            state->setCovMatrix( trkCov );          
            state->setReferencePoint( posLocal ); //Is this needed if we store the global position of the hit?
            state->setLocation( sensorID ) ; //  EUTelTrackStateImpl::AtFirstHit );
            state->setTx(0.);         // seed with 0 at first hit. Given by beam direction.
            state->setTy(0.);         // seed with 0 at first hit. Given by beam direction.
            state->setX(posGlobal[0]);          // 0. at first hit
            state->setY(posGlobal[1]);          // 0. at first hit
            state->setInvP(invp);            // independent of reference point
						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    state->setZParameter(posGlobal[2]);

            if( streamlog_level(DEBUG1) ){            
              state->Print();
            }
            EUTelTrackImpl* track = new EUTelTrackImpl;
						state->setHit(*itHit);
            track->addTrackState( state );

            _tracksCartesian.push_back( track );

          }// loop over for the iplane

        }// loop over for the planes to extract seed hits

        streamlog_out(DEBUG2) << "--------------------------------EUTelKalmanFilter::initialiseSeeds()---------------------------" << std::endl;
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
      
      TMatrixDSym Ckkm1 = ts->getTrackStateCov();
      double xyPrec = getWindowSize();   //sqrt( Ckkm1[0][0]*Ckkm1[0][0] + Ckkm1[1][1]*Ckkm1[1][1] );
      
      streamlog_out(DEBUG0) << "Minimal combined UV resolution : " << xyPrec << std::endl;
      streamlog_out(DEBUG2) << "----------------------EUTelKalmanFilter::getXYPredictionPrecision()------------------------" << std::endl;

      return xyPrec;
    }

        
    /**
     * Propagate track state k-1 -> k
     * @param ts track state to update
     */
    void EUTelKalmanFilter::propagateTrackState( EUTelTrackStateImpl* ts ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::propagateTrackState()" << std::endl;
        TVectorD xkm1 = ts->getTrackStateVec();
        TVectorD xkkm1 = _jacobianF * xkm1;
        
        streamlog_message( DEBUG0, xkkm1.Print();, std::endl; );
        
        ts->setX( xkkm1[0] );
        ts->setY( xkkm1[1] );
        ts->setTx( xkkm1[2] );
        ts->setTy( xkkm1[3] );
        ts->setInvP( xkkm1[4] );
        
        streamlog_out( DEBUG2 ) << "-----------------------------------EUTelKalmanFilter::propagateTrackState()----------------------------------" << std::endl;
    }

	void EUTelKalmanFilter::nextStateUsingJacobianFinder(EUTelTrackStateImpl* input, EUTelTrackStateImpl* output, TMatrixD& jacobian){
        	streamlog_out( DEBUG5 ) << "EUTelKalmanFilter::nextStateUsingJacobianFinder()" << std::endl;
		/////////////////////////////////////////////////////////////////////////////////////////////////////Here we update the new global position of the track
		TVectorD xkm1 = input->getTrackStateVec();
	        streamlog_out( DEBUG4 ) << "Before transformation" << std::endl;
	        streamlog_message( DEBUG4, xkm1.Print();, std::endl; );
	        TVectorD xkkm1 = jacobian * xkm1; 
	        streamlog_out( DEBUG4 ) << "After transformation" << std::endl;
        	streamlog_message( DEBUG4, xkkm1.Print();, std::endl; );
        
	        output->setX( xkkm1[0] );
	        output->setY( xkkm1[1] );
	        output->setTx( xkkm1[2] );
	        output->setTy( xkkm1[3] );
	        output->setInvP( xkkm1[4] );

		///////////////////////////////////////////////////////////Here we update the Covariant matrix
       		TMatrixDSym TrackCovInput = input->getTrackStateCov();
	        TMatrixDSym TrackCovOutput  = TrackCovInput.Similarity( jacobian ); 

        	float trkCov[15] = { static_cast<float>(TrackCovOutput[0][0]), static_cast<float>(TrackCovOutput[1][0]), static_cast<float>(TrackCovOutput[1][1]), 
                	             static_cast<float>(TrackCovOutput[2][0]), static_cast<float>(TrackCovOutput[2][1]), static_cast<float>(TrackCovOutput[2][2]),
                        	     static_cast<float>(TrackCovOutput[3][0]), static_cast<float>(TrackCovOutput[3][1]), static_cast<float>(TrackCovOutput[3][2]),
	                             static_cast<float>(TrackCovOutput[3][3]), static_cast<float>(TrackCovOutput[4][0]), static_cast<float>(TrackCovOutput[4][1]),
        	                     static_cast<float>(TrackCovOutput[4][2]), static_cast<float>(TrackCovOutput[4][3]), static_cast<float>(TrackCovOutput[4][4]) };

	        input->setCovMatrix( trkCov );    
				
        
        	streamlog_out( DEBUG5 ) << "-----------------------------------EUTelKalmanFilter::nextStateUsingJacobianFinder()----------------------------------" << std::endl;
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

	//Determined hit on plane ins local coordinates 
        const double* uvpos = hit->getPosition();
        TVectorD mk(3); 
        mk[0] = uvpos[0];          
	mk[1] = uvpos[1];
	
        streamlog_out( DEBUG3 ) << "Hit (id=" << hit->id() << ") local(u,v) coordinates of hit: (" << mk[0] << "," << mk[1] <<","<<mk[2] << ")" << std::endl;
	/////////////////////////////////////////////////////////////////////////

	//Determine local coordinates of state prediction////////////////////////
	TVectorD prediction(3);
	double localState[3];
	const double input[3] = {ts->getX(),ts->getY(),ts->getZParameter()};

       	streamlog_out( DEBUG3 ) << "	Global(u,v,z) coordinates before transform to local  (" << input[0] << "," << input[1] << "," << input[2] << ")" << " planeID: " <<  ts->getLocation() << std::endl;

	geo::gGeometry().master2Localtwo( ts->getLocation(), input, localState );
	prediction[0] = localState[0];	prediction[1] = localState[1]; prediction[2] = localState[2];
	
        streamlog_out( DEBUG3 ) << "	Prediction for hit (id=" << hit->id() << ") local(u,v) coordinates of state: ("  << prediction[0] << "," << prediction[1] <<","<<prediction[2] << ")" << std::endl;
	//////////////////////////////////////////////////////////////////////////

        
      
        TVectorD rk(2);
        rk[0] = mk[0] - prediction[0];
        rk[1] = mk[1] - prediction[1];
				
				
        
        if ( streamlog_level(DEBUG2) ){
            streamlog_out( DEBUG2 ) << "	Residual vector rk: (" <<rk[0] <<","<<rk[1]<<")"<< std::endl;
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
        TMatrixD Hk = ts->getH();
        
        _processNoiseQ.Zero();
        TMatrixDSym Ckm1 = ts->getTrackStateCov( ); // 5x5
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF ); // 5x5        //Ckkm1 += _processNoiseQ;       
        TMatrixDSym Rkkm1 = Ckkm1.Similarity(Hk); // 5x5
        Rkkm1 += getHitCov(hit); // hitCov is 2x2 ??
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Residual covariance matrix:" << std::endl;
            Rkkm1.Print();
        }
        
        streamlog_out( DEBUG2 ) << "-----------------------------------EUTelKalmanFilter::getResidualCov()------------------------------------" << std::endl;
        
        return Rkkm1;
    }
    
    /** Retrieve Kalman gain matrix. This is a 5x2 matrix so takes a state 5x1 as an argument.
     * 
     * @param ts track state
     * @param hit
     * 
     * @return Gain matrix K
     */
    TMatrixD EUTelKalmanFilter::updateGainK( const EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::updateGainK()" << std::endl;

	TMatrixD gainK(5,2);
        _processNoiseQ.Zero(); //Not sure the need for this?
        TMatrixDSym Ckm1 = ts->getTrackStateCov( );
        TMatrixDSym Ckkm1 = Ckm1.Similarity( _jacobianF ); //Transform covariant matrix from one position to the next      
        TMatrixD Ht(5,2);     Ht = Ht.Transpose( ts->getH() );//This matrix transforms from local to global.
        
        gainK = Ckkm1 * Ht * getResidualCov( ts, hit ).Invert();//take residual of hit and track in local frame. Transform to global. Then multiply by global state covariant matrix
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Gain matrix:" << std::endl;
            gainK.Print();
        }
        
        streamlog_out( DEBUG2 ) << "----------------------------------------EUTelKalmanFilter::updateGainK()------------------------------------" << std::endl;
        
        return gainK;
    }
    
    /** Update track state given new hit
     * 
     * @param ts track state
     * @param hit
     */


		void EUTelKalmanFilter::UpdateStateUsingHitInformation(EUTelTrackStateImpl* input,EVENT::TrackerHit* hit, const TMatrixD& jacobian, TMatrixD & KGain, TMatrixD & HMatrix){
       streamlog_out( DEBUG2 ) << "-----------------------EUTelKalmanFilter::UpdateStateUsingHitInformation()-------------------------------START" << std::endl;
				//Get the residual of the hit and the track and the state vector/////////
        TVectorD residual = getResidual( input, hit ); //This is just the components of distance in x and y
        TVectorD state = input->getTrackStateVec( );       
				///////////////////////////////////////////////////////////////////////

				///////////////////////////////////////////////////////First the state. Note that if the hit is accurate you want to change the state to that position
        state += KGain * residual;   //If the hit is very certain KGain will be unity so the new state position will be the hit position. If the hit is very uncertain then KGain is very small and has no effect. Note matrix (5x1)=+(5x2)*(2x1). 
        input -> setX( state[0] );
        input -> setY( state[1] );
        input -> setTx( state[2] );
        input -> setTy( state[3] );
        input -> setInvP( state[4] );
				///////////////////////////////////////////////////////////////////////////
				
		/////////////////////////////////////////////////////////////////Now the Cov Matrix. Note that if the hit is very accurate it will reduce the Cov matrix
		TMatrixDSym CovMatrix = input->getTrackStateCov();
		TMatrixD I(5,5);
		I.UnitMatrix();
		I -= KGain*HMatrix;  //If we have low uncertainty in the hit then I will be 0. So the covariant will be 0 at this state // Matrix    (5x5)=(5x2)*(2x5). Simply transform gain from
        	TMatrixD newCovMatrix = I*CovMatrix;

	        float trkCov[15] = { static_cast<float>(newCovMatrix[0][0]), static_cast<float>(newCovMatrix[1][0]), static_cast<float>(newCovMatrix[1][1]), 
                             static_cast<float>(newCovMatrix[2][0]), static_cast<float>(newCovMatrix[2][1]), static_cast<float>(newCovMatrix[2][2]),
                             static_cast<float>(newCovMatrix[3][0]), static_cast<float>(newCovMatrix[3][1]), static_cast<float>(newCovMatrix[3][2]),
                             static_cast<float>(newCovMatrix[3][3]), static_cast<float>(newCovMatrix[4][0]), static_cast<float>(newCovMatrix[4][1]),
                             static_cast<float>(newCovMatrix[4][2]), static_cast<float>(newCovMatrix[4][3]), static_cast<float>(newCovMatrix[4][4]) };

        input->setCovMatrix( trkCov );
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				input->setHit(hit);
        streamlog_out( DEBUG2 ) << "-----------------EUTelKalmanFilter::UpdateStateUsingHitInformation()-------------------------------END" << std::endl;

	}

	
	void EUTelKalmanFilter::UpdateTrackUsingHitInformation( EUTelTrackStateImpl* input,const EVENT::TrackerHit* hit, EUTelTrackImpl* track, const TMatrixD& jacobian, TMatrixD & KGain, TMatrixD & HMatrix){
				streamlog_out( DEBUG2 ) << "-----------------EUTelKalmanFilter::UpdateTrackUsingHitInformation()-------------------------------BEGIN" << std::endl;			
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Now determine the chi2 of the track.
				TMatrixDSym newCovMatrix = input->getTrackStateCov( );
        TVectorD residual(2);   residual = getResidual( input, hit ); //This is just the components of distance in x and y

        TMatrixD HMatrixTranspose(5,2);      HMatrixTranspose.Transpose( HMatrix );
        TMatrixD newCovMatrixMeas(2,2);  newCovMatrixMeas = HMatrix * newCovMatrix * HMatrixTranspose; //This is the new state covariant matrix in measurements space.
				TMatrixD hitCov(2,2);       hitCov = getHitCov( hit );
				hitCov -= newCovMatrixMeas;
				TMatrixD I(2,2);
				I.UnitMatrix();
				I -= HMatrix*KGain;
				TVectorD change_residual(2); change_residual =  I*residual;

				double chi2 = hitCov.Invert().Similarity(change_residual);
        track->setChi2( chi2 );
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				streamlog_out( DEBUG5 ) << "Add hit: " << hit << " to track " << track << std::endl;
				streamlog_out( DEBUG2 ) << "-----------------EUTelKalmanFilter::UpdateTrackUsingHitInformation()-------------------------------END" << std::endl;
	
	}



    double EUTelKalmanFilter::updateTrackState( EUTelTrackStateImpl* ts, const EVENT::TrackerHit* hit ) {
        streamlog_out( DEBUG2 ) << "EUTelKalmanFilter::updateTrackState()" << std::endl;
        TVectorD rkkm1(2);      rkkm1 = getResidual( ts, hit );
        TVectorD xk(5);         xk = ts->getTrackStateVec( );
        
        TMatrixD Kk = updateGainK( ts, hit );
        xk += Kk * rkkm1;
        
        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Updated track parameters:" << std::endl;
            xk.Print();
        }
        
        TMatrixD Hk = ts->getH();
        TMatrixD I(5,5);     I.UnitMatrix();
        I -= Kk*Hk;

        if ( streamlog_level(DEBUG0) ) {
            streamlog_out( DEBUG0 ) << "Track parameters projection matrix Hk:" << std::endl;
            Hk.Print();
            streamlog_out( DEBUG0 ) << "Gain matrix Kk:" << std::endl;
            Kk.Print();
        }
        
        _processNoiseQ.Zero();
        TMatrixDSym Ckm1 = ts->getTrackStateCov( );
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
    
    IMPL::TrackImpl* EUTelKalmanFilter::cartesian2LCIOTrack( EUTelTrackImpl* track ) const {

        IMPL::TrackImpl* LCIOtrack = new IMPL::TrackImpl;


        int nstates =  track->getTrackStates().size();
        for(int i=0;i < nstates; i++) 
        {
//           EUTelTrackStateImpl* nexttrackstate = (track->getTrackStates().at(i)) ;  
          EUTelTrackStateImpl* nexttrackstate = new EUTelTrackStateImpl( *(track->getTrackStates().at(i)) ); // memory leak source ?
          IMPL::TrackStateImpl* implstate     = static_cast <IMPL::TrackStateImpl*> (nexttrackstate );
          LCIOtrack->addTrackState( implstate );
        }

        // Assign hits to LCIO TRACK
        const EVENT::TrackerHitVec& trkcandhits = track->getTrackerHits();
        EVENT::TrackerHitVec::const_iterator itrHit;
        for ( itrHit = trkcandhits.begin(); itrHit != trkcandhits.end(); ++itrHit ) {
             LCIOtrack->addHit( *itrHit );
        }

        return LCIOtrack;

    }
    
    /**
     * Distributes hits among measurement layers in the order
     * seen by track traversing the telescope planes 
     * 
     * @param hits Vector of hits to be assigned to the measurement layers of Klaman filter
     */
    bool EUTelKalmanFilter::sortHitsByMeasurementLayers( const EVENT::TrackerHitVec& hits ) {
        bool isReady = true;
        
        if ( hits.empty() ) {
            streamlog_out(WARNING1) << "No hits supplied this event." << std::endl;
            isReady = false;
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

            const int sensorID = Utility::getSensorIDfromHit( *itr );

            itMeasLayer = measLayers.find( sensorID );
            if ( itMeasLayer != measLayers.end() ) itMeasLayer->second->addHit( *itr );
            else {
                measLayers[ sensorID ] = new MeasurementLayer( sensorID );
                measLayers[ sensorID ]->addHit( *itr );
            }
        }
        
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::sortHitsByMeasurementLayers()  start loop through the .sensorIDsVec()"  << std::endl;
        int nplanes        =  geo::gGeometry().nPlanes();
        int nplanesVecSize =  geo::gGeometry().sensorIDsVec().size();
        streamlog_out(DEBUG2) << "EUTelKalmanFilter::sortHitsByMeasurementLayers()  nplanes : " << nplanes << " vectorSize: " << nplanesVecSize << std::endl;

        // Sort measurement layers such that layers encountered by track first
        // are in the front of array
        _allMeasurements = std::vector< MeasurementLayer* >( nplanesVecSize , NULL );   // flush vector
        int numberAlongZ = -1;
        std::vector< int >::const_iterator itSensorID;
        for ( itSensorID = geo::gGeometry().sensorIDsVec().begin(); itSensorID != geo::gGeometry().sensorIDsVec().end(); ++itSensorID ) {
            sensorID = (*itSensorID);
            streamlog_out(DEBUG1) << " sensorID  : " << sensorID << std::endl;
            itMeasLayer = measLayers.find( sensorID );
            if ( itMeasLayer != measLayers.end() ) {
                numberAlongZ = geo::gGeometry().sensorIDtoZOrder( sensorID );
                _allMeasurements.at( numberAlongZ ) = itMeasLayer->second;
            }
        }
        
        // remove elements without MeasurementLayer assigned
        for ( itLayer = _allMeasurements.begin(); itLayer != _allMeasurements.end(); ) {
            if ( (*itLayer) == NULL ) itLayer = _allMeasurements.erase(itLayer);
            else ++itLayer;
        }
        
        return isReady;
    }
    
    MeasurementLayer::MeasurementLayer() : _id(-1), _allHits() {}
    
    MeasurementLayer::MeasurementLayer( int id ) : _id(id), _allHits() {}
    
    MeasurementLayer::~MeasurementLayer() { }
    
    void MeasurementLayer::addHit( EVENT::TrackerHit* hit ) {
        _allHits.push_back( hit );
    }    
    
} // namespace eutelescope
