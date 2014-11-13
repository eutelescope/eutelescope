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
#include "gearimpl/SiPlanesLayerLayoutImpl.h"
#include "gearimpl/SiPlanesParametersImpl.h"

#include "gearimpl/TrackerPlanesLayerLayoutImpl.h"
#include "gearimpl/TrackerPlanesParametersImpl.h"

#include "gear/BField.h"

// EUTELESCOPE
#include "EUTelUtility.h"
#include "EUTelGenericPixGeoMgr.h"
//#include "EUTelGenericPixGeoDescr.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMatrixD.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif


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

        // Iterate over registered GEAR objects and construct their TGeo representation
        const Double_t PI     = 3.141592653589793;
        const Double_t DEG    = 180./PI; 
        const Double_t RADIAN = PI/180.; 

        class EUTelGeometryTelescopeGeoDescription {
            
        private:
            /** */
            EUTelGeometryTelescopeGeoDescription();

            /** */ 
            DISALLOW_COPY_AND_ASSIGN(EUTelGeometryTelescopeGeoDescription)      // prevent users from making (default) copies of processors

            /** need only for pede2lcio*/
            gear::GearMgr* _gearManager;


            /** */ 
            bool _siPlanesDefined;

            /** */ 
            bool _telPlanesDefined;

        public:
            /** Retrieves the instanstance of geometry.
             * Performs lazy intialization if necessary.
             * @TODO this routine has to be considered to be constant
             */
            static EUTelGeometryTelescopeGeoDescription& getInstance( gear::GearMgr* _g );
  
            /** */
            void updateGearManager();  
 
            /** */
            unsigned counter() { return _counter++; }
						void setInitialDisplacementToFirstPlane(float initialDisplacement);

            /** needed only for pede2lcio*/ 
            void setGearManager( gear::GearMgr* value ) { _gearManager = value ; }

            /** Number of planes in the setup */
            inline size_t getSiPlanesLayoutID() const { return _siPlanesLayoutID; } ;

             /** Number of planes in the setup */
            void setSiPlanesLayoutID(size_t value) { _siPlanesLayoutID = value; } ;          
                    
            /** Number of planes in the setup */
            size_t nPlanes() const;
            
            /** Z coordinates of centers of planes */
            const EVENT::DoubleVec& siPlanesZPositions() const;
           
            /** set methods */

            /** set X position  */
            void setPlaneXPosition(int sensorID, double value);
 
            /** set Y position  */
            void setPlaneYPosition(int sensorID, double value);
 
            /** set Z position  */
            void setPlaneZPosition(int sensorID, double value);
 
            /** set X rotation  */
            void setPlaneXRotation(int sensorID, double value);
 
            /** set Y rotation  */
            void setPlaneYRotation(int sensorID, double value);
 
            /** set Z rotation  */
            void setPlaneZRotation(int sensorID, double value);
 
            /** set X rotation  */
            void setPlaneXRotationRadians(int sensorID, double value /* in Radians */);
 
            /** set Y rotation  */
            void setPlaneYRotationRadians(int sensorID, double value /* in Radians */);
 
            /** set Z rotation  */
            void setPlaneZRotationRadians(int sensorID, double value /* in Radians */);

            /** */ 
            float siPlaneRotation1(int sensorID);

            /** */ 
            float siPlaneRotation2(int sensorID);

            /** */ 
            float siPlaneRotation3(int sensorID);

            /** */ 
            float siPlaneRotation4(int sensorID);
 
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

             /** Rotation around X axis of the global coordinate frame */
            double siPlaneXRotationRadians( int );
            
            /** Rotation around Y axis of global coordinate frame */
            double siPlaneYRotationRadians( int );
            
            /** Rotation around Z axis of global coordinate frame */
            double siPlaneZRotationRadians( int );

          
            /** Sensor X side size */
            double siPlaneXSize( int );
            
            /** Sensor Y side size */
            double siPlaneYSize( int );
            
            /** Sensor Z side size */
            double siPlaneZSize( int );
 
            /** Sensor X side pixel pitch [mm] */
            double siPlaneXPitch( int );
            
            /** Sensor Y side pixel pitch [mm] */
            double siPlaneYPitch( int );

            /** Sensor X side size in pixels */
            double siPlaneXNpixels( int );
            
            /** Sensor Y side size in pixels */
            double siPlaneYNpixels( int );
 
            /** Sensor X side size in pixels */
            double siPlaneXResolution( int );
            
            /** Sensor Y side size in pixels */
            double siPlaneYResolution( int );
            
            /** Sensor medium radiation length */
            double siPlaneRadLength( int );
            
	    /** Name of pixel geometry library */
	    std::string geoLibName( int );
            
	    /** Plane normal vector (nx,ny,nz) */
            TVector3 siPlaneNormal( int );
            TVector3 siPlaneXAxis( int);
						TVector3 siPlaneYAxis( int );
            void initialisePlanesToExcluded(FloatVec planeIDs );
            /** Map from sensor ID to number along Z */
            const std::map<int, int>& sensorZOrdertoIDs() const;
            
            std::map<int, int>& sensorZOrderToIDWithoutExcludedPlanes(); 
						std::map<int,int>& sensorIDToZOrderWithoutExcludedPlanes();
            /** Map from sensor ID to number along Z */
            const std::map<int, int>& sensorIDstoZOrder() const;
            
            int sensorIDtoZOrder( int ) const;
            
            int sensorZOrderToID( int ) const;

	
            /** Vector of all sensor IDs */
            const EVENT::IntVec& sensorIDsVec() const;

        public:
            virtual ~EUTelGeometryTelescopeGeoDescription();

        private:
            /** reading initial info from gear: part of contructor */
	    void readSiPlanesLayout();

            /** reading initial info from gear: part of contructor */
	    void updateSiPlanesLayout();


            /** reading initial info from gear: part of contructor
              * new GEAR from branch/TelPlanes
              */
	    void readTrackerPlanesLayout(); 
           
            /**  */
	    void updateTrackerPlanesLayout(); 


            /** housing for the above two 
              */    
            void readGear();

            void translateSiPlane2TGeo(TGeoVolume*,int );

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
           
            	void local2Master( int, const double[], double[] );

		void local2masterHit(EVENT::TrackerHit* hit_input, IMPL::TrackerHitImpl* hit_output, LCCollection * hitCollectionOutput);
		
		void master2localHit(EVENT::TrackerHit* hit_input, IMPL::TrackerHitImpl* hit_output, LCCollection * hitCollectionOutput);
            
            	void master2Local( const double[], double[] );

            	void master2Localtwo(int, const double[], double[] );

		void local2MasterVec( int, const double[], double[] );
 
		void master2LocalVec( int, const double[], double[] );

		int findIntersectionWithCertainID( float x0, float y0, float z0, float px, float py, float pz, float beamQ, int nextPlaneID, float outputPosition[],TVector3& outputMomentum, float& arcLength );
		TVector3 getXYZMomentumfromArcLength(TVector3 momentum, TVector3 globalPositionStart, float charge, float  arcLength );
		float getInitialDisplacementToFirstPlane() const;

		TVector3 getXYZfromArcLength( TVector3 pos,TVector3 pVec , float _beamQ, double s) const;
		TMatrixD getPropagationJacobianCurvilinear(float ds, float qbyp, TVector3 t1, TVector3 t2);
		TMatrixD getLocalToCurvilinearTransformMatrix(TVector3 globalMomentum, int  planeID, float charge);


		TMatrix getPropagationJacobianF( float x0, float y0, float z0, float px, float py, float pz, float _beamQ, float dz );

                void CalculateProjMatrix( TMatrixD& proL2m, double* hitPointGlobal )
		{  
		// Calculate projection matrix

		const TGeoHMatrix* globalH = getHMatrix( hitPointGlobal );
		const TGeoHMatrix& globalHInv = globalH->Inverse();
		const double* rotation = globalHInv.GetRotationMatrix();

		proL2m[0][0] = rotation[0]; // x projection, xx
		proL2m[0][1] = rotation[1]; // y projection, xy
		proL2m[1][0] = rotation[3]; // x projection, yx
		proL2m[1][1] = rotation[4]; // y projection, yy

    		}


        	const TGeoHMatrix* getHMatrix( const double globalPos[] );
            
            /** Magnetic field */
            const gear::BField& getMagneticFiled() const;

			/** Returns a pointer to the EUTelGenericPixGeoDescr of given plane */
			EUTelGenericPixGeoDescr* getPixGeoDescr( int );

			/** Returns the TGeo path of given plane */
			std::string  getPlanePath( int  );

        private:
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

            /**
             */
            gear::TrackerPlanesParameters*  _trackerPlanesParameters;
 
            /**
             */
            gear::TrackerPlanesLayerLayout* _trackerPlanesLayerLayout;

// overwrite private to public ::
        private :

            /** */
            size_t _siPlanesLayoutID;
						float _initialDisplacement; 

            /** Vector of Sensor IDs */
            EVENT::IntVec _sensorIDVec;

            /** Sensor ID map (inverse sensorIDVec) */
            std::map< int, int > _sensorIDVecMap;

            /** Map from number along the Z axis (beam axis) to sensor ID */
            std::map<int, int> _sensorZOrderToIDMap;

            /** Map from sensor ID to number along Z */
            std::map<int, int> _sensorIDtoZOrderMap;

						std::map<int,int> _sensorZOrderToIDWithoutExcludedPlanes;
            /** X coordinate of the sensors centers in global coordinate frame [mm]*/

						std::map<int, int> _sensorIDToZOrderWithoutExcludedPlanes;

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
           
            /** deprecated rotaion natrix elements */
            EVENT::DoubleVec _siPlaneRotation1; 

            /** deprecated rotaion natrix elements */
            EVENT::DoubleVec _siPlaneRotation2; 

            /** deprecated rotaion natrix elements */
            EVENT::DoubleVec _siPlaneRotation3; 

            /** deprecated rotaion natrix elements */
            EVENT::DoubleVec _siPlaneRotation4; 

	    /** Sensor X side length [mm]*/
            EVENT::DoubleVec _siPlaneXSize;
            
            /** Sensor Y side length [mm]*/
            EVENT::DoubleVec _siPlaneYSize;
            
            /** Sensor Z side length [mm]*/
            EVENT::DoubleVec _siPlaneZSize;
 
            /** Sensor X side pitch length [mm]*/
            EVENT::DoubleVec _siPlaneXPitch;
            
            /** Sensor Y side pitch length [mm]*/
            EVENT::DoubleVec _siPlaneYPitch;
 
            /** Sensor X side pitch length [pixels]*/
            EVENT::DoubleVec _siPlaneXNpixels;
            
            /** Sensor Y side pitch length [pixels]*/
            EVENT::DoubleVec _siPlaneYNpixels;

            /** Sensor X side pitch length [pixels]*/
            EVENT::DoubleVec _siPlaneXResolution;
            
            /** Sensor Y side pitch length [pixels]*/
            EVENT::DoubleVec _siPlaneYResolution;
            
            /** Radiation length of the sensor [mm]*/
            EVENT::DoubleVec _siPlaneRadLength;

	    /** Name of the pixel geometry library for each plane*/
	    EVENT::StringVec _geoLibName;

            /** Number of planes including DUT */
            size_t _nPlanes;

            /** Pointer to the pixel geometry manager */
            EUTelGenericPixGeoMgr* _pixGeoMgr;
            //#ifdef  USE_TGEO

        private:
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

            /** */
            static unsigned _counter;

        };
        
        inline EUTelGeometryTelescopeGeoDescription& gGeometry( gear::GearMgr* _g = marlin::Global::GEAR ) {
                return EUTelGeometryTelescopeGeoDescription::getInstance( _g ); 
        }

        
    } // namespace geo
} // namespace eutelescope

#endif  // USE_GEAR
 
#endif	/* EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H */

