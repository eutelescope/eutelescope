/* 
 * File:   EUTelGeometryTelescopeGeoDescription.h
 *
 */

#ifndef EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H
#define	EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H

// C++
#include <map>
#include <string>

// LCIO includes
#include "LCIOSTLTypes.h"

// MARLIN
#include "marlin/Global.h"

// GEAR includes
#include "gear/GearMgr.h"
#include "gear/SiPlanesLayerLayout.h"
#include "gear/SiPlanesParameters.h"
#include "gear/BField.h"

// EUTELESCOPE
#include "EUTelUtility.h"
#include "EUTelGenericPixGeoMgr.h"
//#include "EUTelGenericPixGeoDescr.h"

//#ifdef USE_TGEO
// ROOT
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TVector3.h"


//#endif //USE_TGEO

// built only if GEAR is available
#ifdef USE_GEAR

/** @class EUTelGeometryTelescopeGeoDescription
 * This class is supposed to keep globally accesible 
 * telescope geometry description.
 * 
 * It is based on singleton design pattern and furnishes
 * a facade for GEAR description
 */

namespace eutelescope {

    namespace geo {


        class EUTelGeometryTelescopeGeoDescription {
            
        private:
            EUTelGeometryTelescopeGeoDescription();

            DISALLOW_COPY_AND_ASSIGN(EUTelGeometryTelescopeGeoDescription)      // prevent users from making (default) copies of processors


        public:
            /** Retrieves the instanstance of geometry.
             * Performs lazy intialization if necessary.
             * @TODO this routine has to be considered to be constant
             */
            static EUTelGeometryTelescopeGeoDescription& getInstance();
 
           
            /** Number of planes in the setup */
            size_t nPlanes() const;
            
            /** Z coordinates of centers of planes */
            const EVENT::DoubleVec& siPlanesZPositions() const;
            
            /** X coordinate of center of sensor 
             * with given ID in global coordinate frame */
            double siPlaneXPosition( int );
            
            /** Y coordinate of center of sensor 
             * with given ID in global coordinate frame */
            double siPlaneYPosition( int );
            
            /** Z coordinate of center of sensor 
             * with given ID in global coordinate frame */
            double siPlaneZPosition( int );
            
            /** Rotation around X axis of the global coordinate frame */
            double siPlaneXRotation( int );
            
            /** Rotation around Y axis of global coordinate frame */
            double siPlaneYRotation( int );
            
            /** Rotation around Z axis of global coordinate frame */
            double siPlaneZRotation( int );
            
            /** Sensor X side size */
            double siPlaneXSize( int );
            
            /** Sensor Y side size */
            double siPlaneYSize( int );
            
            /** Sensor Z side size */
            double siPlaneZSize( int );
            
            /** Sensor medium radiation length */
            double siPlaneMediumRadLen( int );
            
	    /** Name of pixel geometry library */
	    std::string geoLibName( int );
            
	    /** Plane normal vector (nx,ny,nz) */
            TVector3 siPlaneNormal( int );
            
            
            /** Map from sensor ID to number along Z */
            const std::map<int, int>& sensorZOrdertoIDs() const;
            
            /** Map from sensor ID to number along Z */
            const std::map<int, int>& sensorIDstoZOrder() const;
            
            int sensorIDtoZOrder( int ) const;
            
            int sensorZOrderToID( int ) const;
            
            /** Vector of all sensor IDs */
            const EVENT::IntVec& sensorIDsVec() const;

        public:
            virtual ~EUTelGeometryTelescopeGeoDescription();

        public:
            // TGeo stuff

            /** Initialize TGeo geometry 
             * Establish access to TGeoManager, load geometry description file.
             * 
             * @param tgeofilename name of a .root or .gdml file with valid geometry
             * 
             * @see ROOT TGeoManager::Import
             */
            void initializeTGeoDescription(std::string tgeofilename);
            
            void initializeTGeoDescription( std::string& geomName, bool dumpRoot );

            // Geometry operations
        public:
            float findRadLengthIntegral( const double[], const double[], bool );
            
            int getSensorID( const float globalPos[] ) const;
           
            int getSensorID( const IMPL::TrackerHitImpl* ) const;
 
            void local2Master( int, const double[], double[] );

						void local2masterHit(EVENT::TrackerHit* hit_input, IMPL::TrackerHitImpl* hit_output, LCCollection * hitCollectionOutput);

						void master2localHit(EVENT::TrackerHit* hit_input, IMPL::TrackerHitImpl* hit_output, LCCollection * hitCollectionOutput);
            
            void master2Local( const double[], double[] );

			void local2MasterVec( int, const double[], double[] );
 
			void master2LocalVec( int, const double[], double[] );
            
            const TGeoHMatrix* getHMatrix( const double globalPos[] );
            
            /** Magnetic field */
            const gear::BField& getMagneticFiled() const;

			/** Returns a pointer to the EUTelGenericPixGeoDescr of given plane */
			EUTelGenericPixGeoDescr* getPixGeoDescr( int );

			/** Returns the TGeo path of given plane */
			std::string  getPlanePath( int  );

        public:
            /** Silicon planes parameters as described in GEAR
             * This structure actually contains the following:
             *  @li A reference to the telescope geoemtry and layout
             *  @li An integer number saying if the telescope is w/ or w/o DUT
             *  @li An integer number saying the number of planes in the
             *  telescope.
             *
             *  This object is provided by GEAR during the init() phase and
             *  stored here for local use.
             */
            gear::SiPlanesParameters* _siPlanesParameters;

            /** Silicon plane layer layout
             * This is the real geoemetry description. For each layer
             *  composing the telescope the relevant information are
             *  available.
             *
             *  This object is taken from the _siPlanesParameters during the
             *  init() phase and stored for local use
             */
            gear::SiPlanesLayerLayout* _siPlanesLayerLayout;


        private:
            /** Vector of Sensor IDs */
            EVENT::IntVec _sensorIDVec;

            /** Sensor ID map (inverse sensorIDVec) */
            std::map< int, int > _sensorIDVecMap;

            /** Map from number along the Z axis (beam axis) to sensor ID */
            std::map<int, int> _sensorZOrderToIDMap;

            /** Map from sensor ID to number along Z */
            std::map<int, int> _sensorIDtoZOrderMap;

            /** X coordinate of the sensors centers in global coordinate frame [mm]*/
            EVENT::DoubleVec _siPlaneXPosition;
            
            /** Y coordinate of the sensors centers in global coordinate frame [mm]*/
            EVENT::DoubleVec _siPlaneYPosition;
            
            /** Z coordinate of the sensors centers in global coordinate frame [mm]*/
            EVENT::DoubleVec _siPlaneZPosition;
            
            /** Rotation around X axis of the global coordinate frame [rad]*/
            EVENT::DoubleVec _siPlaneXRotation;
            
            /** Rotation around Y axis of global coordinate frame [rad]*/
            EVENT::DoubleVec _siPlaneYRotation;
            
            /** Rotation around Z axis of global coordinate frame [rad]*/
            EVENT::DoubleVec _siPlaneZRotation;
            
            /** Sensor X side length [mm]*/
            EVENT::DoubleVec _siPlaneSizeX;
            
            /** Sensor Y side length [mm]*/
            EVENT::DoubleVec _siPlaneSizeY;
            
            /** Sensor Z side length [mm]*/
            EVENT::DoubleVec _siPlaneSizeZ;
            
            /** Radiation length of the sensor [mm]*/
            EVENT::DoubleVec _siPlaneRadLength;

	    /** Name of the pixel geometry library for each plane*/
	    EVENT::StringVec _geoLibName;

            /** Number of planes including DUT */
            size_t _nPlanes;

            /** Pointer to the pixel geometry manager */
			EUTelGenericPixGeoMgr* _pixGeoMgr;
            //#ifdef  USE_TGEO

			/** Flag if geoemtry is already initialized */
			bool _isGeoInitialized;

			/** Map containing plane path (string) and corresponding planeID */
			std::map<int, std::string> _planePath;

        public:
            // TGeo stuff
            /** @TODO this must be coupled with GEAR
             * description. No checks of consistency between GEAR and TGeo
             * descriptions are being done, currently.
             */

            /** Geometry manager global object */
            TGeoManager* _geoManager;
            //#endif // USE_TGEO

            int findNextPlaneEntrance(  double* ,  double *, int, float*  );
            int findNextPlane(  double* lpoint,  double* ldir,  float* newpoint );

        };
        
        inline EUTelGeometryTelescopeGeoDescription& gGeometry() {
                return EUTelGeometryTelescopeGeoDescription::getInstance(); 
        }

        
    } // namespace geo
} // namespace eutelescope

#endif  // USE_GEAR

#endif	/* EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H */

