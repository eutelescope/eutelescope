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
#include "EUTelAPIXSparsePixel.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelFFClusterImpl.h"

// lcio includes <.h>
#include "IMPL/TrackerHitImpl.h"

#include <EVENT/LCEvent.h>

#include <cstdio>

namespace eutelescope {

    namespace Utility {

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
         * 
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

                    EUTelVirtualCluster * cluster;

                    if (hit->getType() == kEUTelSparseClusterImpl) {

                        TrackerDataImpl * clusterFrame = dynamic_cast<TrackerDataImpl*> (clusterVector[0]);
                        if (clusterFrame == 0) {
                            // found invalid result from cast
                            throw UnknownDataTypeException("Invalid hit found in method hitContainsHotPixels()");
                        }

                        eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > *cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > (clusterFrame);
                        int sensorID = cluster->getDetectorID();

                        for (unsigned int iPixel = 0; iPixel < cluster->size(); iPixel++) {
                            EUTelSimpleSparsePixel m26Pixel;
                            cluster->getSparsePixelAt(iPixel, &m26Pixel);
                            {
                                char ix[100];
                                sprintf(ix, "%d,%d,%d", sensorID, m26Pixel.getXCoord(), m26Pixel.getYCoord());
                                std::map < std::string, bool >::const_iterator z = hotPixelMap.find(ix);
                                if (z != hotPixelMap.end() && hotPixelMap.at(ix) == true) {
                                    skipHit = true;
                                    streamlog_out(DEBUG3) << "Skipping hit as it was found in the hot pixel map." << std::endl;
                                    return true; // if TRUE  this hit will be skipped
                                } else {
                                    skipHit = false;
                                }
                            }
                        }

                    } else if (hit->getType() == kEUTelBrickedClusterImpl) {

                        // fixed cluster implementation. Remember it
                        //  can come from
                        //  both RAW and ZS data

                        cluster = new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));

                    } else if (hit->getType() == kEUTelDFFClusterImpl) {

                        // fixed cluster implementation. Remember it can come from
                        // both RAW and ZS data
                        cluster = new EUTelDFFClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                    } else if (hit->getType() == kEUTelFFClusterImpl) {

                        // fixed cluster implementation. Remember it can come from
                        // both RAW and ZS data
                        cluster = new EUTelFFClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                    }
                    else if (hit->getType() == kEUTelAPIXClusterImpl) {
                        TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> (clusterVector[0]);

                        cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > (clusterFrame);

                        eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > (clusterFrame);

                        int sensorID = apixCluster->getDetectorID();

                        for (unsigned int iPixel = 0; iPixel < apixCluster->size(); ++iPixel) {
                            EUTelAPIXSparsePixel apixPixel;
                            apixCluster->getSparsePixelAt(iPixel, &apixPixel);
                            {
                                char ix[100];
                                sprintf(ix, "%d,%d,%d", sensorID, apixPixel.getXCoord(), apixPixel.getYCoord());
                                std::map < std::string, bool >::const_iterator z = hotPixelMap.find(ix);
                                if (z != hotPixelMap.end() && hotPixelMap.at(ix) == true) {
                                    skipHit = true;
                                    streamlog_out(DEBUG3) << "Skipping hit as it was found in the hot pixel map." << std::endl;
                                    return true; // if TRUE  this hit will be skipped
                                } else {
                                    skipHit = false;
                                }
                            }
                        }

                        return skipHit; // if TRUE  this hit will be skipped
                    }

//                    delete cluster;
                    
                } catch (lcio::Exception e) {
                    // catch specific exceptions
                    streamlog_out(ERROR) << "Exception occured in hitContainsHotPixels(): " << e.what() << std::endl;
                }
            } catch (...) {
                // if anything went wrong in the above return FALSE, meaning do not skip this hit
                return 0;
            }

            // if none of the above worked return FALSE, meaning do not skip this hit
            return 0;

        }
        
        EUTelVirtualCluster* GetClusterFromHit( const IMPL::TrackerHitImpl* hit ) {
            
            LCObjectVec clusterVector = hit->getRawHits();
            
            if (hit->getType() == kEUTelBrickedClusterImpl) {
                        return new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                    } else if (hit->getType() == kEUTelDFFClusterImpl) {
                        return new EUTelDFFClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                    } else if (hit->getType() == kEUTelFFClusterImpl) {
                        return new EUTelFFClusterImpl(static_cast<TrackerDataImpl *> (clusterVector[0]));
                    } else if (hit->getType() == kEUTelAPIXClusterImpl) {
                        return new EUTelSparseClusterImpl< EUTelAPIXSparsePixel >(static_cast<TrackerDataImpl *> (clusterVector[0]));
                    } else if (hit->getType() == kEUTelSparseClusterImpl) {
			return new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > (static_cast<TrackerDataImpl *> (clusterVector[0]));
                    } else {
                        throw UnknownDataTypeException("Unknown cluster type");
                    }
            
        }

        /**
         * Determine hit's plane id
         * @param hit 
         * @return plane id
         */
        int GuessSensorID( const EVENT::TrackerHit* hit ) {
            if ( hit == NULL ) {
                streamlog_out(ERROR) << "An invalid hit pointer supplied! will exit now\n" << std::endl;
                return -1;
            }

            try {
                EUTelVirtualCluster * cluster = GetClusterFromHit( static_cast< const IMPL::TrackerHitImpl*> (hit) );

                if ( cluster != NULL ) {
                    int sensorID = cluster->getDetectorID();
		    delete cluster;
                    return sensorID;
                }
            } catch (...) {
                streamlog_out(ERROR) << "guessSensorID() produced an exception!" << std::endl;
            }

            return -1;
        }     
        
        std::map<std::string, bool > FillHotPixelMap( EVENT::LCEvent *event, const std::string& hotPixelCollectionName ) {
            
            std::map < std::string, bool > hotPixelMap;
            
            LCCollectionVec *hotPixelCollectionVec = 0;
            try {
                hotPixelCollectionVec = static_cast<LCCollectionVec*> (event->getCollection(hotPixelCollectionName));
            } catch (...) {
                streamlog_out(MESSAGE) << "hotPixelCollectionName " << hotPixelCollectionName.c_str() << " not found" << std::endl;
                return hotPixelMap;
            }

            CellIDDecoder<TrackerDataImpl> cellDecoder(hotPixelCollectionVec);

            for (int i = 0; i < hotPixelCollectionVec -> getNumberOfElements(); i++) {
                TrackerDataImpl* hotPixelData = dynamic_cast<TrackerDataImpl*> (hotPixelCollectionVec->getElementAt(i));
                SparsePixelType type = static_cast<SparsePixelType> (static_cast<int> (cellDecoder(hotPixelData)["sparsePixelType"]));

                int sensorID = static_cast<int> (cellDecoder(hotPixelData)["sensorID"]);

                if (type == kEUTelAPIXSparsePixel) {
                    std::auto_ptr< EUTelSparseDataImpl< EUTelAPIXSparsePixel > > apixData(new EUTelSparseDataImpl< EUTelAPIXSparsePixel > (hotPixelData));
                    std::vector< EUTelAPIXSparsePixel* > apixPixelVec;
                    EUTelAPIXSparsePixel apixPixel;

                    //Push all single Pixels of one plane in the apixPixelVec
                    for (unsigned int iPixel = 0; iPixel < apixData->size(); iPixel++) {
                        std::vector<int> apixColVec();
                        apixData->getSparsePixelAt(iPixel, &apixPixel);
                        streamlog_out(MESSAGE) << iPixel << " of " << apixData->size() << " HotPixelInfo:  " << apixPixel.getXCoord() << " " << apixPixel.getYCoord() << " " << apixPixel.getSignal() << std::endl;
                        try {
                            char ix[100];
                            sprintf(ix, "%d,%d,%d", sensorID, apixPixel.getXCoord(), apixPixel.getYCoord());
                            hotPixelMap[ix] = true;
                        } catch (...) {
                            streamlog_out(ERROR) << "problem adding pixel to hotpixel map! " << std::endl;
                        }
                    }

                } else if (type == kEUTelSimpleSparsePixel) {
                    std::auto_ptr< EUTelSparseClusterImpl< EUTelSimpleSparsePixel > > m26Data(new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > (hotPixelData));
                    std::vector< EUTelSimpleSparsePixel* > m26PixelVec;
                    EUTelSimpleSparsePixel m26Pixel;

                    //Push all single Pixels of one plane in the m26PixelVec
                    for (unsigned int iPixel = 0; iPixel < m26Data->size(); iPixel++) {
                        std::vector<int> m26ColVec();
                        m26Data->getSparsePixelAt(iPixel, &m26Pixel);
                        streamlog_out(MESSAGE) << iPixel << " of " << m26Data->size() << " HotPixelInfo:  " << m26Pixel.getXCoord() << " " << m26Pixel.getYCoord() << " " << m26Pixel.getSignal() << std::endl;
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
                //Solutions
                vector< double > X(2, 0.);              //initialise with two doubles equal 0.

                if( fabs( a ) > 0. )
                {
                        //The equation has the form
                        // a*x^2 + b*x + c = 0
                        double disc2 =  b*b - 4.*a*c ;
                        if( disc2 < 0. )
                        {
                                cout << "WARNING! disc2 < 0: " << disc2 << endl;
                                return X;
                        }
                        double disc = sqrt( disc2 );
                        double denom = 2.*a;
                        double num1 = -b + disc;
                        double num2 = -b - disc;

                        X[0] = num1 / denom;            //bigger root
                        X[1] = num2 / denom;            //lower root
                }
                else
                {
                        //Degenerate case, when a = 0.
                        //The linear equation has the form
                        // b*x + c = 0

                        X[0] = -c/b;
                        X[1] = -c/b;
                }

                return X;
        }
        
    }
}
