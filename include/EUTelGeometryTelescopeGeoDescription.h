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

// EUTELESCOPE
#include "EUTelUtility.h"

//#ifdef USE_TGEO
// ROOT
#include "TGeoManager.h"

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

            // Geometry operations
        public:
            void findRad(Double_t x, Double_t y, Double_t z,
                    Double_t theta, Double_t phi, Int_t &nbound, Float_t &length, Float_t &safe, Float_t &rad, Bool_t verbose);

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


        public:
            //    /** Ordered sensor ID
            //     *  This vector contains sensorID sorted according to their 
            //     *  position along the Z axis (beam axis)
            //     */
            //    EVENT::IntVec _orderedSensorID;
            //
            //    /** @TODO: Add description to this variable */
            //    EVENT::IntVec _orderedSensorID_wo_excluded;

            /** Vector of Sensor IDs */
            EVENT::IntVec _sensorIDVec;

            /** Sensor ID map (inverse sensorIDVec) */
            std::map< int, int > _sensorIDVecMap;

            /** Sensor ID vector ordered according to their 
             *  position along the Z axis (beam axis) */
            EVENT::IntVec _sensorIDVecZOrder;

            /** Map from sensor ID to number along Z */
            std::map<int, int> _sensorIDtoZOrderMap;

            // an associative map for getting also the sensorID ordered
            /** @TODO: Add description to this variable */
            std::map< double, int > _sensorIDMap;

            /** @TODO: Add description to this variable */
            EVENT::DoubleVec _siPlaneZPosition;

            /** Number of planes including DUT */
            size_t _nPlanes;


            //#ifdef  USE_TGEO
        public:
            // TGeo stuff
            /** @TODO this must be coupled with GEAR
             * description. No checks of consistency between GEAR and TGeo
             * descriptions are being done, currently.
             */

            /** Geometry manager global object */
            TGeoManager* _geoManager;
            //#endif // USE_TGEO

        };
        
        inline EUTelGeometryTelescopeGeoDescription& gGeometry() {
                return EUTelGeometryTelescopeGeoDescription::getInstance(); 
        }
        
    } // namespace geo
} // namespace eutelescope

#endif  // USE_GEAR

#endif	/* EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H */

