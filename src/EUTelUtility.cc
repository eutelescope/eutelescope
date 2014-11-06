/* 
 * File:   EUTelUtility.cpp
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Created on January 22, 2013, 11:54 AM
 */

// eutelescope includes ".h"
#include "EUTelUtility.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelFFClusterImpl.h"

// lcio includes <.h>
#include <EVENT/LCEvent.h>

#include <cstdio>

using namespace std;

namespace eutelescope {

    namespace Utility {


			std::string outputColourString(std::string inputString, std::string colour){
				std::string outputString;
				if(colour == "RED"){
					outputString = "\033[31m" +  inputString + "\033[39m";
					return	outputString;
				}
				if(colour =="BLUE"){
					outputString = "\033[34m" +  inputString + "\033[39m";
					return	outputString;
				}
				if(colour =="GREEN"){
					outputString = "\033[32m" +  inputString + "\033[39m";
					return	outputString;
				}
				if(colour =="YELLOW"){
					outputString = "\033[33m" +  inputString + "\033[39m";
					return	outputString;
				}
	
			}
		 


        /**
         * Fills indices of not excluded planes
         * 
         * 
         *      *               *
         * |    |       ...     |       |
         * 0    1       ...     k      k+1
         * 
         * Plane marked with (*) to be excluded from consideration
         * 
         *      *                               *       
         * |    |       |       |       ...     |       |
         * 0   -1       1       2       ...    -1      k-2
         * 
         * Array of indices of not excluded planes
        & * 
         * @param indexconverter 
         *              returned indices of not excluded planes
         * 
         * @param excludePlanes
         *              array of plane ids to be excluded
         * 
         * @param nPlanes
         *              total number of planes
         */
				
        void FillNotExcludedPlanesIndices( std::vector<int>& indexconverter, const std::vector<unsigned int >& excludePlanes,  unsigned int nPlanes ) {
            int icounter = 0;
            int nExcludePlanes = static_cast< int >( excludePlanes.size());
            for( unsigned int i = 0; i < nPlanes; i++) {
                int excluded = 0; //0 - not excluded, 1 - excluded
                if ( nExcludePlanes > 0) {
                    for( int j = 0; j < nExcludePlanes; j++ ) {
                        if ( i == excludePlanes[ j ] ) {
                            excluded = 1;
                            break;
                        }
                    }
                }
                if ( excluded == 1 )
                    indexconverter[ i ] = -1;
                else {
                    indexconverter[ i ] = icounter;
                    icounter++;
                }
            }
            streamlog_out( DEBUG ) << "FillNotExcludedPlanesIndices" << std::endl;
        }
        
        bool HitContainsHotPixels( const IMPL::TrackerHitImpl* hit, const std::map<std::string, bool >& hotPixelMap ) {
            bool skipHit = false;

            try {
                try {
                    LCObjectVec clusterVector = hit->getRawHits();

                    EUTelVirtualCluster * cluster = NULL;

                    if (hit->getType() == kEUTelSparseClusterImpl) {

                        TrackerDataImpl * clusterFrame = dynamic_cast<TrackerDataImpl*> (clusterVector[0]);
                        if (clusterFrame == 0) {
                            // found invalid result from cast
                            throw UnknownDataTypeException("Invalid hit found in method hitContainsHotPixels()");
                        }

                        eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel > *cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelGenericSparsePixel > (clusterFrame);
                        int sensorID = cluster->getDetectorID();

                        for (unsigned int iPixel = 0; iPixel < cluster->size(); iPixel++) {
                            EUTelGenericSparsePixel m26Pixel;
                            cluster->getSparsePixelAt(iPixel, &m26Pixel);
                            {
                                char ix[100];
                                sprintf(ix, "%d,%d,%d", sensorID, m26Pixel.getXCoord(), m26Pixel.getYCoord());
                                std::map < std::string, bool >::const_iterator z = hotPixelMap.find(ix);
                                if (z != hotPixelMap.end() && hotPixelMap.at(ix) == true) {
                                    skipHit = true;
                                    streamlog_out(DEBUG3) << "Skipping hit as it was found in the hot pixel map." << std::endl;
                                    break;
//                                    delete cluster;
//                                    return true; // if TRUE  this hit will be skipped
                                } else {
                                    skipHit = false;
                                }
                            }
                        }
                        delete cluster;
                    } else if (hit->getType() == kEUTelBrickedClusterImpl) {

                        // fixed cluster implementation. Remember it
                        //  can come from
                        //  both RAW and ZS data

                        cluster = new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                        delete cluster;
                    } else if (hit->getType() == kEUTelDFFClusterImpl) {

                        // fixed cluster implementation. Remember it can come from
                        // both RAW and ZS data
                        cluster = new EUTelDFFClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                        delete cluster;
                    } else if (hit->getType() == kEUTelFFClusterImpl) {

                        // fixed cluster implementation. Remember it can come from
                        // both RAW and ZS data
                        cluster = new EUTelFFClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                        delete cluster;
                    }
                    
//                if ( cluster != 0 ) delete cluster;
                    
                } catch (lcio::Exception e) {
                    // catch specific exceptions
                    streamlog_out(ERROR) << "Exception occured in hitContainsHotPixels(): " << e.what() << std::endl;
                }
            } catch (...) {
                // if anything went wrong in the above return FALSE, meaning do not skip this hit
                return 0;
            }

            // if none of the above worked return FALSE, meaning do not skip this hit
            return skipHit;

        }
        
        /**
         * Provides access to raw cluster information for given hit
         * Constructed object is owned by caller. Cluster must be destroyed by caller.
         * 
         * @param hit 
         * @return raw data cluster information
         */
	std::auto_ptr<EUTelVirtualCluster> GetClusterFromHit( const IMPL::TrackerHitImpl* hit ) 
	{
            
            LCObjectVec clusterVector = hit->getRawHits();
            
            if (hit->getType() == kEUTelBrickedClusterImpl) 
	    {
            	return std::auto_ptr<EUTelVirtualCluster>( new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl*>(clusterVector[0])) );
            } 
	    else if (hit->getType() == kEUTelDFFClusterImpl) 
	    {
		    return std::auto_ptr<EUTelVirtualCluster>( new EUTelDFFClusterImpl(static_cast<TrackerDataImpl*>(clusterVector[0])) );
            } 
	    else if (hit->getType() == kEUTelFFClusterImpl) 
	    {
		    return std::auto_ptr<EUTelVirtualCluster>( new EUTelFFClusterImpl(static_cast<TrackerDataImpl*>(clusterVector[0])) );
	    } 
	    else if (hit->getType() == kEUTelSparseClusterImpl) 
	    {
		    return std::auto_ptr<EUTelVirtualCluster>( new EUTelSparseClusterImpl<EUTelGenericSparsePixel>(static_cast<TrackerDataImpl*>(clusterVector[0])) );
            }
	    else 
	    {
		    streamlog_out(WARNING2) << "Unknown cluster type: " << hit->getType() << std::endl;
                    return std::auto_ptr<EUTelVirtualCluster>();
	    }
        }

        /**
         * Determine hit's plane id
         * @param hit 
         * @return plane id
         */
        int getSensorIDfromHit( EVENT::TrackerHit* hit ) {
            if ( hit == NULL ) {
                streamlog_out(ERROR) << "getSensorIDfromHit:: An invalid hit pointer supplied! will exit now\n" << std::endl;
                return -1;
            }

            try {

                UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );

                int sensorID = hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit))["sensorID"];
                return sensorID;

            } catch (...) {
                streamlog_out(ERROR) << "getSensorIDfromHit() produced an exception!" << std::endl;
            }

            return -1;
        }     
 
        std::map<std::string, bool > FillHotPixelMap( EVENT::LCEvent *event, const std::string& hotPixelCollectionName ) {
            
            std::map < std::string, bool > hotPixelMap;
            
            LCCollectionVec *hotPixelCollectionVec = 0;
            try {
                hotPixelCollectionVec = static_cast<LCCollectionVec*> (event->getCollection(hotPixelCollectionName));
            } catch (...) {
                streamlog_out( MESSAGE4 ) << "hotPixelCollectionName " << hotPixelCollectionName.c_str() << " not found" << std::endl;
                return hotPixelMap;
            }

            CellIDDecoder<TrackerDataImpl> cellDecoder(hotPixelCollectionVec);

            for (int i = 0; i < hotPixelCollectionVec -> getNumberOfElements(); i++) {
                TrackerDataImpl* hotPixelData = dynamic_cast<TrackerDataImpl*> (hotPixelCollectionVec->getElementAt(i));
                SparsePixelType type = static_cast<SparsePixelType> (static_cast<int> (cellDecoder(hotPixelData)["sparsePixelType"]));

                int sensorID = static_cast<int> (cellDecoder(hotPixelData)["sensorID"]);

                if (type == kEUTelGenericSparsePixel) {
                    std::auto_ptr< EUTelSparseClusterImpl< EUTelGenericSparsePixel > > m26Data(new EUTelSparseClusterImpl< EUTelGenericSparsePixel > (hotPixelData));
                    std::vector< EUTelGenericSparsePixel* > m26PixelVec;
                    EUTelGenericSparsePixel m26Pixel;

                    //Push all single Pixels of one plane in the m26PixelVec
                    for (unsigned int iPixel = 0; iPixel < m26Data->size(); iPixel++) {
                        std::vector<int> m26ColVec();
                        m26Data->getSparsePixelAt(iPixel, &m26Pixel);
                        streamlog_out(DEBUG0) << iPixel << " of " << m26Data->size() << " HotPixelInfo:  " << m26Pixel.getXCoord() << " " << m26Pixel.getYCoord() << " " << m26Pixel.getSignal() << std::endl;
                        try {
                            char ix[100];
                            sprintf(ix, "%d,%d,%d", sensorID, m26Pixel.getXCoord(), m26Pixel.getYCoord());
                            hotPixelMap[ix] = true;
                        } catch (...) {
                            std::cout << "can not add pixel " << std::endl;
                            std::cout << sensorID << " " << m26Pixel.getXCoord() << " " << m26Pixel.getYCoord() << " " << std::endl;
                        }
                    }
                }
            }
            return hotPixelMap;
        }

        /** Highland's formula for multiple scattering 
         * @param p momentum of the particle [GeV/c]
         * @param x thickness of the material in units of radiation lenght
         */
        double getThetaRMSHighland( double p, double x ) {
            double tet = (0.0136 * sqrt(x) / p * (1 + 0.038 * std::log(x)));
            return tet;
        }
        
        /** Calculate median 
         * Sort supplied vector and determine the median
         */
        double getMedian( std::vector<double>& vec ) {
            std::sort( vec.begin( ), vec.end( ) );
            double median = -999.;
            size_t size = vec.size( );
            if ( size % 2 == 0 ) {
                median = ( vec[size / 2 - 1] + vec[size / 2] ) / 2;
            } else {
                median = vec[size / 2];
            }
            return median;
        }
        
        /**
         * Calculate 2D curvature of the track with given pt and charge q
         * in solenoidal magnetic field with strength B
         * @param pt transverse momentum of the track [GeV/c]
         * @param B  magnetic field strength [T]
         * @param q  particle charge [e]
         * @return 1/R 2D curvature of the track [1/m]
         */
        double getCurvature( double pt, double B, double q ) {
            double rho = 0.;
            if ( pt > 0. ) rho = 0.299792458 * q * B / pt;
            
            return rho;
        }
        
        /**
         *  Solves quadratic equation a*x^2 + b*x + c = 0
         * @param a 
         * @param b
         * @param c
         * @return vector of solution sorted in descending order
         */
        vector< double > solveQuadratic( double a, double b, double c) {
		streamlog_out( DEBUG1 ) << "Solving quadratic equation with coefficients:\na: " 
		<< a << "\nb: " << b << "\nc:" << c << std::endl;
                //Solutions
                vector< double > X;              //initialise with two doubles equal 0.

                if( fabs( a ) > 1.E-10 )
                {
                        //The equation has the form
                        // a*x^2 + b*x + c = 0
                        double disc2 =  b*b - 4.*a*c ;
                        if( disc2 < 0. )
                        {
                                cout << " Quadratic equation solution is imaginary! WARNING! disc2 < 0: " << disc2 << endl;
                                return X;
                        }
                        double disc = sqrt( disc2 );
                        double denom = 2.*a;
                        double num1 = -b + disc;
                        double num2 = -b - disc;

			X.push_back( num1 / denom );	// larger root
			X.push_back( num2 / denom );	// smaller root
                        //X[0] = num1 / denom;            // bigger root
                        //X[1] = num2 / denom;            // lower root
                }
                else
                {
                        //Degenerate case, when a = 0.
                        //The linear equation has the form
                        // b*x + c = 0

			X.push_back( -c/b );
			X.push_back( -c/b );
                        //X[0] = -c/b;
                        //X[1] = -c/b;
                }

                return X;
        }

        /** get cluster size in X and Y for TrackerHit derived hits:
         *  includes only certain type of data
         *  to-do: replace with generic hit and cluster structure not requiring explicit declaration of data types
         */

        void getClusterSize(const IMPL::TrackerHitImpl * hit, int& sizeX, int& sizeY ) {
        // rewrite from EUTelAXPITbTrackTuple:
         if(hit==0)
         {
            streamlog_out( ERROR5 ) << "An invalid hit pointer supplied! will exit now\n" << endl;
            return ;
         }

	  std::auto_ptr<EUTelVirtualCluster> cluster = GetClusterFromHit( hit ) ;

	       if(cluster.get() != NULL )  cluster->getClusterSize(sizeX, sizeY);
 
        }
  
 	void copyLCCollectionHitVec(  IMPL::LCCollectionVec* input, LCCollectionVec* output ) {


		// Prepare output collection
//	 	LCFlagImpl flag( input->getFlag() );
//	  	flag.setBit( LCIO::TRBIT_HITS );
//	  	output->setFlag( flag.getFlag( ) );


	   	// deep copy of all elements  - requires clone of original elements
	   	//
	   	int nElements = input->getNumberOfElements() ;

           	streamlog_out( DEBUG4) << "HIT : copy of n= " << nElements << " element for collection : " << input->getTypeName() << std::endl;

	   	for(int i=0; i< nElements ; i++){
	        	IMPL::TrackerHitImpl *hit = static_cast<IMPL::TrackerHitImpl *> ( input->getElementAt(i) );
                	streamlog_out( DEBUG4) << " i= " << i << " type : " << hit->getType() << std::endl;
    			output->push_back(  hit ) ;
	   	}
        
        }


	//Create once per event     !!
	void copyLCCollectionTrackVec(  IMPL::LCCollectionVec* input,  LCCollectionVec* output) {

		// Prepare output collection
  		LCFlagImpl flag( input->getFlag() );
  		flag.setBit( LCIO::TRBIT_HITS );
  		output->setFlag( flag.getFlag( ) );

	   	// deep copy of all elements  - requires clone of original elements
	   	//
	   	int nElements = input->getNumberOfElements() ;

                streamlog_out( DEBUG4) << "TRACK: copy of n= " << nElements << " element for collection : " << input->getTypeName() << std::endl;
   
	   	for(int i=0; i< nElements ; i++){
                        IMPL::TrackImpl *trk = static_cast<IMPL::TrackImpl *> ( input->getElementAt(i) );
                        streamlog_out( DEBUG4) << " i= " << i << 
                                                    " type : " << trk->getType() << 
                                                    " nstates: " << trk->getTrackStates().size() << 
                                                    " nhits: " << trk->getTrackerHits().size() << 
                                                    std::endl;
	     		output->push_back(  trk ) ;
//	     		output->push_back(  trk->clone() ) ;
	   	}
        
        }

        float DoubleToFloat(double a){
                return static_cast<float> (a);  
        }

        float* toFloatN(double* a, int N){ 
           float *vec = new float[N];
           for(int i=0;i<N;i++) {
             vec[i] = DoubleToFloat( a[i]) ;  
           }
           return vec;            
        }

        const float* HitCDoubleShiftCFloat(const double* hitPosition, TVectorD& residual ){

             double hit[3];
             std::copy(hitPosition, hitPosition+3, hit);             

             hit[0] -= residual[0]; // why minus?
             hit[1] -= residual[1];

             streamlog_out(MESSAGE1)  << " " << hit[0]  << " " << hit[1] << " " << hit[2]  << std::endl;

             float *opoint = new float[3];
                         
 //c++0x     std::ransform( hitPointLocal, hitPointLocal+3, const_cast<float*>( temp), [](double hitPointLocal){return (float) hitPointLocal;}  );
             std::transform( hit, hit+3, opoint, DoubleToFloat  );
	  return opoint;
        }

  }
}
