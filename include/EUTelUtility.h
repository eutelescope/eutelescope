/* 
 * File:   EUTelUtility.h
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Various routines and classes useful for EUTelescope data analysis
 *
 * Created on January 22, 2013, 11:35 AM
 */

#ifndef EUTELUTILITY_H
#define	EUTELUTILITY_H

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// lcio includes <.h>
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackImpl.h"
#include "EVENT/LCEvent.h"

// ROOT
#include "TVectorD.h" 

// system includes <>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#ifndef DISALLOW_COPY_AND_ASSIGN
//Following #define stops the accidental creation of a copy or assignment operator by causing a link error.
//Copy and Assignment operators not allowed because they are unnecessary and the cause of many bugs
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
TypeName(const TypeName&); \
void operator=(const TypeName&);
#endif

namespace eutelescope {

    namespace Utility {

				std::string outputColourString(std::string inputString, std::string colour);
			
        class Hit;

        /**
         * Vector of Hits in plane
         */
        typedef std::vector< Hit > HitsVec;
        typedef std::vector< Hit* > HitsPVec;

        /**
         * Represents geometric properties of hits
         */
        class Hit {
        public:

            Hit() : _PlaneID(-1), _x(0.0), _y(0.0), _z(0.0) {
            }

            Hit(int id,
                    double x,
                    double y,
                    double z) : _PlaneID(id), _x(x), _y(y), _z(z) {
            }

            bool operator<(const Hit& b) const {
                return ( _z < b.Z());
            }

            inline void SetZ(double z) {
                this->_z = z;
            }

            inline double Z() const {
                return _z;
            }

            inline void SetY(double y) {
                this->_y = y;
            }

            inline double Y() const {
                return _y;
            }

            inline void SetX(double x) {
                this->_x = x;
            }

            inline double X() const {
                return _x;
            }

            void SetPlaneID(int _PlaneID) {
                this->_PlaneID = _PlaneID;
            }

            int GetPlaneID() const {
                return _PlaneID;
            }

            void Print() {
                streamlog_out(DEBUG0) << "Hit:" << std::endl;
                streamlog_out(DEBUG0) << "PlaneID" << std::setw(10) << setiosflags(std::ios::right) << _PlaneID << std::endl;
                streamlog_out(DEBUG0) << "X:" << std::setw(10) << std::setprecision(6) << setiosflags(std::ios::right) << _x << std::endl;
                streamlog_out(DEBUG0) << "Y:" << std::setw(10) << std::setprecision(6) << setiosflags(std::ios::right) << _y << std::endl;
                streamlog_out(DEBUG0) << "Z:" << std::setw(10) << std::setprecision(6) << setiosflags(std::ios::right) << _z << std::endl;
            }

        protected:
            int _PlaneID;

            double _x;
            double _y;
            double _z;
        };

        /**
         * Specify a Rectangular in a sensor
         */
        class SensorRectangular {
        protected:
            int sensor; // SensorID
            int A; // lowest pixel in X direction
            int B; // lowest pixel in Y direction
            int C; // highest pixel in X direction
            int D; // highest pixel in Y direction
        public:

            SensorRectangular(int s, int a, int b, int c, int d) : sensor(s), A(a), B(b), C(c), D(d) {
            };

            SensorRectangular() : sensor(0), A(0), B(0), C(0), D(0) {
            };

            int getSensor() const {
                return sensor;
            }

            //look if x and y are inside the foreseen rectangular

            bool isInside(int x, int y) const {
                return (x >= A && x <= C && y >= B && y <= D);
            }

            void print() {
                streamlog_out(MESSAGE4) << "Sensor: " << sensor << ": (" << A << "|" << B << ") to (" << C << "|" << D << ")" << std::endl;
            }

        };

        class RectangularArray {
        protected:
            std::map<int, SensorRectangular > _rect;

        public:

            void addRectangular(SensorRectangular &s) {
                _rect[s.getSensor() ] = s;
            }

            bool isInside(int s, int x, int y) {
                std::map<int, SensorRectangular >::iterator it = _rect.find(s);
                if (it == _rect.end()) { // not in the map means no limit on this sensor -> always true
                    return true;
                }
                SensorRectangular cSensor = _rect[s];
                return cSensor.isInside(x, y);
            }

        };

        /*!
         * Fills indices of not excluded planes
         */

        void FillNotExcludedPlanesIndices(
                std::vector< int >&,
                const std::vector< unsigned int >&,
                unsigned int = 0);

        /*!
         * Checks if a hit belongs to a hot-pixel-map
         */

        //! Called for first event per run
        /*! Reads hotpixel information from hotPixelCollection into hotPixelMap
         * to be used in the sensor exclusion area logic 
         */

        std::map<std::string, bool > FillHotPixelMap(EVENT::LCEvent *event, const std::string& hotPixelCollectionName);

        bool HitContainsHotPixels(const IMPL::TrackerHitImpl * hit, const std::map<std::string, bool >& hotPixelMap);

	std::auto_ptr<EUTelVirtualCluster> GetClusterFromHit(const IMPL::TrackerHitImpl*);

        int getSensorIDfromHit( EVENT::TrackerHit* hit);

        /** Highland's formula for multiple scattering */
        double getThetaRMSHighland( double, double );
        
        /** Calculate median */
        double getMedian(std::vector<double>& );
        
        /** Calculate track's 2D curvature */
        double getCurvature( double, double, double );
        
        /** Solve quadratic equation a*x^2 + b*x + c = 0 */
        std::vector< double > solveQuadratic( double, double, double );
 
        /** getClusterSize method from TrackerHit object:: assumes known cluster types */
        void getClusterSize(const IMPL::TrackerHitImpl * hit, int& sizeX, int& sizeY ) ;

        /** */
        void copyLCCollectionHitVec( LCCollectionVec*, LCCollectionVec* );
 
        /** */
        void copyLCCollectionTrackVec( LCCollectionVec*, LCCollectionVec* );

        /**  type conversion  */
        float DoubleToFloat(double a);
        float* DoubleToFloatN(double* a, int N);

        /** */ 
        const float* HitCDoubleShiftCFloat(const double* , TVectorD& );

        /** Tokenize string */
                /**
         * String tokenizer
         * @param source
         * @param delimiter
         * @param keepEmpty
         * @return vector of tokens
         */
        std::vector<std::string> inline stringSplit(const std::string& source, const char* delimiter = " ", bool keepEmpty = false ) {
            std::vector<std::string> results;

            size_t prev = 0;
            size_t next = 0;

            while ( ( next = source.find_first_of( delimiter, prev ) ) != std::string::npos ) {
                if ( keepEmpty || ( next - prev != 0 ) ) {
                    results.push_back( source.substr( prev, next - prev ) );
                }
                prev = next + 1;
            }

            if ( prev < source.size( ) ) {
                results.push_back( source.substr( prev ) );
            }

            return results;
        }
        
        /** Possible choices of alignment degrees of freedom */
        enum AlignmentMode {
            noAlignment, XYShift, XYShiftXYRot, XYZShiftXYRot, XYShiftYZRotXYRot, XYShiftXZRotXYRot, XYShiftXZRotYZRotXYRot, XYZShiftXZRotYZRotXYRot
        };

        /** @struct BinomialCombination
         * This structure generates combinations of n objects 
         * taken k at a time. Combinations are composed from numbers 0 .. n-1.
         * 
         * @see D. Knuth "Art of programming. Combinatorial Algorithms" Volume 4A
         * 
         * \b Usage:
         * \code{.cpp}
         *  BinomialCombination com(5,3);
         * 
         *  std::vector < std::vector < int > > combinations;
         *  combinations = com.sampleCombinations();
         * \endcode
         * 
         */
        struct BinomialCombination {

            DISALLOW_COPY_AND_ASSIGN( BinomialCombination )     // prevent users from making (default) copies
            
            // Constructors 
            
            BinomialCombination() : _nObjects(1), _kTaken(1) {
            }

            BinomialCombination(int n, int k) : _nObjects(n), _kTaken(k) {
            }

            /** Sample posible combinations of n objekts taken k at
             * one time. Generated combinations are sorted in lexicografical
             * order.
             * 
             * @return vector of possible combinations
             */
            std::vector < std::vector < int > > sampleCombinations() {

               std::vector < int > statevec;    // statevec stores current combination
                // intialise algorithm
                {
                    for (int j = 1; j <= _kTaken; ++j) statevec.push_back(j - 1);
                    statevec.push_back(_nObjects);
                    statevec.push_back(0);
                }

                std::vector < std::vector < int > > result; // result keeps all generated combinations
                // iterate over combinations
                std::vector < int > currentCombination(_kTaken);
                int j = 1;
                do {
                    std::copy( statevec.begin(), statevec.begin() + _kTaken, 
                               currentCombination.begin() );
                    result.push_back( currentCombination );
                    for ( j = 1; (statevec[ j - 1 ] + 1 == statevec[ j ]); ++j ) {
                        statevec[ j - 1 ] = j - 1;
                    }
                    statevec[ j - 1 ]++;
                } while (j <= _kTaken);
                return result;
            }

        private:
            
            /** Total number of objects */
            int _nObjects;
            
            /** Number of objects taken at one time */
            int _kTaken;
        };

    }
}

#endif	/* EUTELUTILITY_H */
