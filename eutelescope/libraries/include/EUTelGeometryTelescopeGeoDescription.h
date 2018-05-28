/*
 * File:   EUTelGeometryTelescopeGeoDescription.h
 *
 */
#ifndef EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H
#define EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H

// C++
#include <array>
#include <map>
#include <memory>
#include <string>

// MARLIN
#include "marlin/Global.h"

// GEAR includes
#include "gear/BField.h"
#include "gear/GearMgr.h"
#include "gearimpl/SiPlanesLayerLayoutImpl.h"
#include "gearimpl/SiPlanesParametersImpl.h"
#include "gearimpl/TrackerPlanesLayerLayoutImpl.h"
#include "gearimpl/TrackerPlanesParametersImpl.h"

// EUTELESCOPE
#include "EUTelGenericPixGeoMgr.h"
#include "EUTelGeoSupportClasses.h"
#include "EUTelUtility.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMatrixD.h"
#else
#error *** You need ROOT to compile this code.  ***
#endif

// Eigen
#include <Eigen/Core>

// ROOT
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TVector3.h"

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

    // Iterate over registered GEAR objects and construct their TGeo
    // representation
    const Double_t PI = 3.141592653589793;
    const Double_t DEG = 180. / PI;
    const Double_t RADIAN = PI / 180.;

    class EUTelGeometryTelescopeGeoDescription {

    private:
      /** */
      EUTelGeometryTelescopeGeoDescription();

      /** */
      DISALLOW_COPY_AND_ASSIGN(
          EUTelGeometryTelescopeGeoDescription) // prevent users from making
                                                // (default) copies of
                                                // processors

      /** need only for pede2lcio*/
      gear::GearMgr *_gearManager;

      /** */
      bool _siPlanesDefined;

      /** */
      bool _telPlanesDefined;

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
      gear::SiPlanesParameters *_siPlanesParameters;

      /** Silicon plane layer layout
       * This is the real geoemetry description. For each layer
       *  composing the telescope the relevant information are
       *  available.
       *
       *  This object is taken from the _siPlanesParameters during the
       *  init() phase and stored for local use
       */
      gear::SiPlanesLayerLayout *_siPlanesLayerLayout;

      /** */
      gear::TrackerPlanesParameters *_trackerPlanesParameters;

      /** */
      gear::TrackerPlanesLayerLayout *_trackerPlanesLayerLayout;

      /** */
      size_t _siPlanesLayoutID;

      float _initialDisplacement;

      /** Vector of Sensor IDs */
      std::vector<int> _sensorIDVec;

      /** Pointer to the pixel geometry manager */
      EUTelGenericPixGeoMgr *_pixGeoMgr;

      /** Flag if geoemtry is already initialized */
      bool _isGeoInitialized;

      /** Map containing plane path (string) and corresponding planeID */
      std::map<int, std::string> _planePath;
	  std::map<int, TGeoMatrix*> _TGeoMatrixMap;
      /** */
      static unsigned _counter;

      std::map<int, TVector3> _planeNormalMap;
      std::map<int, TVector3> _planeXMap;
      std::map<int, TVector3> _planeYMap;
      std::map<int, double> _planeRadMap;

      std::vector<std::unique_ptr<EUTelLayer>> _telescopeLayers;
      std::map<int, EUTelLayer *> _telescopeLayerMap;
      std::map<std::string, EUTelMaterial> _materialMap;
      std::map<int, EUTelActive *> _activeMap;

    public:
      /** Retrieves the instanstance of geometry.
       * Performs lazy intialization if necessary.
       * @TODO this routine has to be considered to be constant
       */
      static EUTelGeometryTelescopeGeoDescription &
      getInstance(gear::GearMgr *_g);

      /** */
      void updateGearManager();

      /** */
      unsigned counter() { return _counter++; }

      auto getActivePlaneMap() const -> decltype(_activeMap) {
        return _activeMap;
      }

      void setInitialDisplacementToFirstPlane(float initialDisplacement) {
        _initialDisplacement = initialDisplacement;
      };

      /** needed only for pede2lcio*/
      void setGearManager(gear::GearMgr *value) { _gearManager = value; }

      /** Number of planes in the setup */
      inline size_t getSiPlanesLayoutID() const { return _siPlanesLayoutID; };

      /** Number of planes in the setup */
      void setSiPlanesLayoutID(size_t value) { _siPlanesLayoutID = value; };

      /** Number of planes in the setup */
      size_t nPlanes() const { return _activeMap.size(); };

      /** set methods */
      /** set X position  */
      /*
              inline void setPlaneOffset(int sensorID, Eigen::Vector3d pos){
         return; this->clearMemoizedValues(); };

              inline void setPlaneAbsolutePosition(int sensorID, Eigen::Vector3d
         pos){ return; this->clearMemoizedValues(); };
              inline void setPlaneAbsolutePosition(int sensorID, double x,
         double y, double z){ return; this->clearMemoizedValues(); };
      */
      inline void alignGlobalPos(int sensorID, Eigen::Vector3d const &pos) {
        std::cout << "Aligning sensor: " << sensorID << " to position: " << pos
                  << std::endl;
        _activeMap.at(sensorID)->alignPos(pos);
        this->clearMemoizedValues();
        return;
      };

      inline void alignGlobalPos(int sensorID, double x, double y, double z) {
        Eigen::Vector3d pos = Eigen::Vector3d(x, y, z);
        std::cout << "Aligning sensor: " << sensorID << " to position: " << pos
                  << std::endl;
        _activeMap.at(sensorID)->alignPos(pos);
        this->clearMemoizedValues();
        return;
      };

      inline void alignGlobalRot(int sensorID, Eigen::Matrix3d const &rot) {
        std::cout << "Aligning sensor: " << sensorID << " to rotation: " << rot
                  << std::endl;
        _activeMap.at(sensorID)->alignRot(rot);
        this->clearMemoizedValues();
        return;
      };

      /*
              inline void setPlaneOffset(int sensorID, double x, double y,
         double z){ return; this->clearMemoizedValues(); };

              inline void setParentLayerPosition(int sensorID, double x, double
         y, double z){
                      _activeMap.at(sensorID)->setParentPosition(x, y, z);
                      this->clearMemoizedValues();
                      return;
              };
      */
      /** set X rotation  */
      //	inline void setPlaneXRotation(int sensorID, double value){
      //_planeSetup[sensorID].alpha = value; this->clearMemoizedValues(); };

      /** set Y rotation  */
      //	inline void setPlaneYRotation(int sensorID, double value){
      //_planeSetup[sensorID].beta = value; this->clearMemoizedValues(); };

      /** set Z rotation  */
      //	inline void setPlaneZRotation(int sensorID, double value){
      //_planeSetup[sensorID].gamma = value; this->clearMemoizedValues(); };

      /** set X rotation in radians */
      //	inline void setPlaneXRotationRadians(int sensorID, double
      //value){ _planeSetup[sensorID].alpha = value*DEG;
      //this->clearMemoizedValues(); };

      /** set Y rotation in radians */
      //	inline void setPlaneYRotationRadians(int sensorID, double
      //value){ _planeSetup[sensorID].beta = value*DEG;
      //this->clearMemoizedValues(); };

      /** set Z rotation in radians */
      //	inline void setPlaneZRotationRadians(int sensorID, double
      //value){ _planeSetup[sensorID].gamma = value*DEG;
      //this->clearMemoizedValues(); };

      inline void setPlanePitch(int sensorID, double xPitch, double yPitch) {
        _activeMap.at(sensorID)->setPitch(xPitch, yPitch);
      }
      inline void setPlaneNoPixels(int sensorID, int xNo, int yNo) {
        _activeMap.at(sensorID)->setNoPixels(xNo, yNo);
      }
      // GETTER
      /** */
      int siPlaneRotation1(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(0, 0);
      };

      /** */
      int siPlaneRotation2(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(1, 0);
      };

      /** */
      int siPlaneRotation3(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(0, 1);
      };

      /** */
      int siPlaneRotation4(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(1, 1);
      };

      auto getPlanePosition(int sensorID) const
          -> decltype(_activeMap.at(sensorID)->getPosition()) {
        return _activeMap.at(sensorID)->getPosition();
      }

      /** X coordinate of center of sensor
       * with given ID in global coordinate frame */
      double siPlaneXPosition(int sensorID) {
        return _activeMap.at(sensorID)->getPosition().coeff(0);
      };

      /** Y coordinate of center of sensor
       * with given ID in global coordinate frame */
      double siPlaneYPosition(int sensorID) {
        return _activeMap.at(sensorID)->getPosition().coeff(1);
      };

      /** Z coordinate of center of sensor
       * with given ID in global coordinate frame */
      double siPlaneZPosition(int sensorID) {
        return _activeMap.at(sensorID)->getPosition().coeff(2);
      };

      /** Rotation around X axis of the global coordinate frame */
      double siPlaneXRotation(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(0)*DEG;
      };

      /** Rotation around Y axis of global coordinate frame */
      double siPlaneYRotation(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(1)*DEG;
      };

      /** Rotation around Z axis of global coordinate frame */
      double siPlaneZRotation(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(2)*DEG;
      };

      /** Rotation around X axis of the global coordinate frame */
      double siPlaneXRotationRadians(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(0);
      };

      /** Rotation around Y axis of global coordinate frame */
      double siPlaneYRotationRadians(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(1);
      };

      /** Rotation around Z axis of global coordinate frame */
      double siPlaneZRotationRadians(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(2);
      };

      /** Sensor X side size */
      double siPlaneXSize(int sensorID) {
        return _activeMap.at(sensorID)->getSize().coeff(0);
      };

      /** Sensor Y side size */
      double siPlaneYSize(int sensorID) {
        return _activeMap.at(sensorID)->getSize().coeff(1);
      };

      /** Sensor Z side size */
      double siPlaneZSize(int sensorID) {
        return _activeMap.at(sensorID)->getSize().coeff(2);
      };

      /** Sensor X side pixel pitch [mm] */
      double siPlaneXPitch(int sensorID) {
        return _activeMap.at(sensorID)->getPitch().first;
      };

      /** Sensor Y side pixel pitch [mm] */
      double siPlaneYPitch(int sensorID) {
        return _activeMap.at(sensorID)->getPitch().second;
      };

      /** Sensor X side size in pixels */
      int siPlaneXNpixels(int sensorID) {
        return _activeMap.at(sensorID)->getNoPixels().first;
      };

      /** Sensor Y side size in pixels */
      int siPlaneYNpixels(int sensorID) {
        return _activeMap.at(sensorID)->getNoPixels().second;
      };

      /** Sensor X side size in pixels */
      double siPlaneXResolution(int sensorID) {
        return _activeMap.at(sensorID)->getResolution().first;
      };

      /** Sensor Y side size in pixels */
      double siPlaneYResolution(int sensorID) {
        return _activeMap.at(sensorID)->getResolution().second;
      };

      /** Sensor medium radiation length */
      double siPlaneRadLength(int sensorID) {
        return _activeMap.at(sensorID)->getRadLength();
      };

      /** Name of pixel geometry library */
      std::string geoLibName(int sensorID) {
        return _activeMap.at(sensorID)->getGeometry();
      };

      /** Plane normal vector (nx,ny,nz) */
      TVector3 siPlaneNormal(int);

      TVector3 siPlaneXAxis(int);

      TVector3 siPlaneYAxis(int);

      /** Vector of all sensor IDs */
      const std::vector<int> &sensorIDsVec() const { return _sensorIDVec; };

      Eigen::Matrix3d rotationMatrixFromAngles(int sensorID);

      Eigen::Vector3d getOffsetVector(int sensorID);

      Eigen::Matrix3i getFlipMatrix(int sensorID);

      Eigen::Vector3d globalXAxis(int sensorID);

      Eigen::Vector3d globalYAxis(int sensorID);

      void writeGEARFile(std::string filename);

      virtual ~EUTelGeometryTelescopeGeoDescription();

      /** Initialize TGeo geometry
       * Establish access to TGeoManager, load geometry description file.
       *
       * @param tgeofilename name of a .root or .gdml file with valid geometry
       *
       * @see ROOT TGeoManager::Import
       */
      void initializeTGeoDescription(std::string tgeofilename);
      void initializeTGeoDescription(std::string const &geomName,
                                     bool dumpRoot);

      double FindRad(Eigen::Vector3d const &startPt,
                     Eigen::Vector3d const &endPt);

      double planeRadLengthGlobalIncidence(int planeID,
                                           Eigen::Vector3d incidenceDir);
      double planeRadLengthLocalIncidence(int planeID,
                                          Eigen::Vector3d incidenceDir);

      void local2Master(int sensorID, std::array<double, 3> const &localPos,
                        std::array<double, 3> &globalPos);
      void master2Local(int sensorID, std::array<double, 3> const &globalPos,
                        std::array<double, 3> &localPos);
      void local2MasterVec(int sensorID, std::array<double, 3> const &localVec,
                           std::array<double, 3> &globalVec);
      void master2LocalVec(int sensorID, std::array<double, 3> const &globalVec,
                           std::array<double, 3> &localVec);

      void local2Master(int, const double[], double[]);
      void master2Local(int, const double[], double[]);
      void local2MasterVec(int, const double[], double[]);
      void master2LocalVec(int, const double[], double[]);

      bool findIntersectionWithCertainID(float x0, float y0, float z0, float px,
                                         float py, float pz, float beamQ,
                                         int nextPlaneID,
                                         float outputPosition[],
                                         TVector3 &outputMomentum,
                                         float &arcLength, int &newNextPlaneID);

      TVector3 getXYZMomentumfromArcLength(TVector3 momentum,
                                           TVector3 globalPositionStart,
                                           float charge, float arcLength);

      // This outputs the total percentage radiation length for the full
      // detector system.
      float calculateTotalRadiationLengthAndWeights(
          const double startD[3], const double endD[3],
          std::map<const int, double> &, std::map<const int, double> &);
      double addKapton(std::map<const int, double> &mapSensor);

      float getInitialDisplacementToFirstPlane() const {
        return _initialDisplacement;
      };

      /** Magnetic field */
      const gear::BField &getMagneticField() const {
        return _gearManager->getBField();
      };

      /** Returns a pointer to the EUTelGenericPixGeoDescr of given plane */
      EUTelGenericPixGeoDescr *getPixGeoDescr(int planeID) {
        return _pixGeoMgr->getPixGeoDescr(planeID);
      };

      /** Returns the TGeo path of given plane */
      std::string getPlanePath(int planeID) {
        return _planePath.find(planeID)->second;
      };

      // TGeo stuff
      /** @TODO this must be coupled with GEAR
       * description. No checks of consistency between GEAR and TGeo
       * descriptions are being done, currently.
       */

      /** Geometry manager global object */
      std::unique_ptr<TGeoManager> _geoManager = nullptr;

    private:
      void updatePlaneInfo(int sensorID);

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

      /** housing for the above two */
      void readGear();

      void translateSiPlane2TGeo(TGeoVolume *, int);

      void clearMemoizedValues() {
        _planeNormalMap.clear();
        _planeXMap.clear();
        _planeYMap.clear();
        _planeRadMap.clear();
      }
    };

    inline EUTelGeometryTelescopeGeoDescription &
    gGeometry(gear::GearMgr *_g = marlin::Global::GEAR) {
      return EUTelGeometryTelescopeGeoDescription::getInstance(_g);
    }
  } // namespace geo
} // namespace eutelescope
#endif // USE_GEAR
#endif /* EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H */
