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

// Eigen
#include <Eigen/Core>

// ROOT
#include "TGeoManager.h"
#include "TGeoMatrix.h"

/** @class EUTelGeometryTelescopeGeoDescription
 * This class is supposed to keep globally accesible
 * telescope geometry description.
 *
 * It is based on singleton design pattern and furnishes
 * a facade for GEAR description
 */
namespace eutelescope {
  namespace geo {
    static const double PI = 3.141592653589793;
    static const double DEG = 180. / PI;
    static const double RADIAN = PI / 180.;

    class EUTelGeometryTelescopeGeoDescription {
    private:
      /** Default constructor */
      EUTelGeometryTelescopeGeoDescription();

      /** need only for pede2lcio*/
      gear::GearMgr *_gearManager;

      /** Flag if the SiPlanesLayerLayout in GEAR is used*/
      bool _siPlanesDefined;

      /** Flag if the rackerPlanesLayerLayout in GEAR is used*/
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

      /** Tracker planes parameters as described in GEAR */
      gear::TrackerPlanesParameters *_trackerPlanesParameters;

      /** Tracker planes layout as described in GEAR */
      gear::TrackerPlanesLayerLayout *_trackerPlanesLayerLayout;

      /** The layout ID from the GEAR files */
      size_t _layoutID;

      /** Vector of Sensor IDs 
       * Ordered according to their global z-position 
       * */
      std::vector<int> _sensorIDVec;

      /** Pointer to the pixel geometry manager */
      std::unique_ptr<EUTelGenericPixGeoMgr> _pixGeoMgr;

      /** Flag if geoemtry is already initialized */
      bool _isGeoInitialized;

      /** Map containing the path to the TGeoNode in ROOT's TGeo framework for each plane (identified by its planeID) */
      std::map<int, std::string> _planePath;

      /** Map holding the transformation matrix for each plane (identified by its planeID) */
	    std::map<int, TGeoMatrix*> _TGeoMatrixMap;

      /** Conter to indicate if instance of this object exists */
      static unsigned _counter;

      /** Map containing the normal vector of each plane */
      std::map<int, Eigen::Vector3d> _planeNormalMap;
      /** Map containing the x-direction vector of each plane */
      std::map<int, Eigen::Vector3d> _planeXMap;
      /** Map containing the y-direction vector of each plane */
      std::map<int, Eigen::Vector3d> _planeYMap;
      /** Map containing the radiation length of each plane */
      std::map<int, double> _planeRadMap;

      /** Map containing all materials defined in GEAR file */
      std::map<std::string, EUTelMaterial> _materialMap;

      //TODO
      std::vector<std::unique_ptr<EUTelLayer>> _telescopeLayers;
      std::map<int, EUTelLayer*> _telescopeLayerMap;
      std::map<int, EUTelActive*> _activeMap;

    public:
      /** Retrieves the instanstance of geometry.
       * Performs lazy intialization if necessary.
       */
      static EUTelGeometryTelescopeGeoDescription &
      getInstance(gear::GearMgr *_g);

      /**TODO */
      void updateGearManager();

      /** Increment the static instance counter */
      unsigned counter() { return _counter++; }

      auto getActivePlaneMap() const -> decltype(_activeMap) {
        return _activeMap;
      }

      /** needed only for pede2lcio*/
      void setGearManager(gear::GearMgr *value) { _gearManager = value; }

      /** Get the GEAR internally used identifier */
      inline size_t getLayoutID() const { return _layoutID; };

      /** Set the GEAR internally used identifier */
      void setLayoutID(size_t value) { _layoutID = value; };

      /** Number of planes in the setup */
      size_t nPlanes() const { return _activeMap.size(); };

      /** Align a given plane (sensorID) to a provided global position (in [mm]) 
       *  As this modifies the geometrical position it is important to clear 
       *  memoized values
       */
      inline void alignGlobalPos(int sensorID, Eigen::Vector3d const &pos) {
        streamlog_out(MESSAGE4) << "Aligning sensor: " << sensorID
                                << " to position: " << pos << std::endl;
        _activeMap.at(sensorID)->alignPos(pos);
        this->clearMemoizedValues();
        return;
      };

      /** Align a given plane (sensorID) to a provided global position (in [mm])
       *  As this modifies the geometrical position it is important to clear 
       *  memoized values
       */
      inline void alignGlobalPos(int sensorID, double const & x, double const & y, double const & z) {
        Eigen::Vector3d pos = Eigen::Vector3d(x, y, z);
        alignGlobalPos(sensorID, pos);
        return;
      };

      /** Align a given plane (sensorID) to provided global rotations 
       *  As this modifies the geometrical position it is important to
       *  clear memoized values
       */
      inline void alignGlobalRot(int sensorID, Eigen::Matrix3d const &rot) {
        streamlog_out(MESSAGE4) << "Aligning sensor: " << sensorID
                                << " to rotation: " << rot << std::endl;
        _activeMap.at(sensorID)->alignRot(rot);
        this->clearMemoizedValues();
        return;
      };

      /** Set the given plane's typical pixel pitch in [mm] 
       *  Note that this in not well defined for irregular sensors, thus this
       *  value should not be used if in doubt - it is intended for dynamic 
       *  histogram binning etc. where a deviation from this value can be
       *  tolerated
       */
      inline void setPlanePitch(int sensorID, double const & xPitch, double const & yPitch) {
        _activeMap.at(sensorID)->setPitch(xPitch, yPitch);
      }

      /** Set the given plane's amoutn of pixels in x- and y-direction
       * I.e. the plane's pixel matrix dimensions
       */
      inline void setPlaneNoPixels(int sensorID, int xNo, int yNo) {
        _activeMap.at(sensorID)->setNoPixels(xNo, yNo);
      }

      /** Get the first flip matrix coefficient for the given plane
       *  Can only be plus or minus one or zero
       */
      int planeFlip1(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(0, 0);
      };

      /** Get the second flip matrix coefficient for the given plane
       *  Can only be plus or minus one or zero
       */
      int planeFlip2(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(1, 0);
      };

      /** Get the third flip matrix coefficient for the given plane
       *  Can only be plus or minus one or zero
       */
      int planeFlip3(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(0, 1);
      };

      /** Get the fourth flip matrix coefficient for the given plane
       *  Can only be plus or minus one or zero
       */
      int planeFlip4(int sensorID) {
        return _activeMap.at(sensorID)->getFlipMatrix().coeff(1, 1);
      };

      /** Returns the given plane's position in global coordinates, in [mm] */ 
      auto getPlanePosition(int sensorID) const
          -> decltype(_activeMap.at(sensorID)->getPosition()) {
        return _activeMap.at(sensorID)->getPosition();
      }

      /** X position of sensor center in the global coordinate frame, in [mm] */
      double getPlaneXPosition(int sensorID) {
        return _activeMap.at(sensorID)->getPosition().coeff(0);
      };

      /** Y position of sensor center in the global coordinate frame, in [mm] */
      double getPlaneYPosition(int sensorID) {
        return _activeMap.at(sensorID)->getPosition().coeff(1);
      };

      /** Z position of sensor center in the global coordinate frame, in [mm] */
      double getPlaneZPosition(int sensorID) {
        return _activeMap.at(sensorID)->getPosition().coeff(2);
      };

      /** Rotation around X axis of the global coordinate frame, in [deg] */
      double getPlaneXRotationDegrees(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(0)*DEG;
      };

      /** Rotation around Y axis of global coordinate frame, in [deg] */
      double getPlaneYRotationDegrees(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(1)*DEG;
      };

      /** Rotation around Z axis of global coordinate frame, in [deg] */
      double getPlaneZRotationDegrees(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(2)*DEG;
      };

      /** Rotation around X axis of the global coordinate frame, in [rad] */
      double getPlaneXRotationRadians(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(0);
      };

      /** Rotation around Y axis of global coordinate frame, in [rad] */
      double getPlaneYRotationRadians(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(1);
      };

      /** Rotation around Z axis of global coordinate frame, in [rad] */
      double getPlaneZRotationRadians(int sensorID) {
        return _activeMap.at(sensorID)->getGlobalRotationAngles().coeff(2);
      };

      /** Sensor X side size, in [mm] */
      double getPlaneXSize(int sensorID) {
        return _activeMap.at(sensorID)->getSize().coeff(0);
      };

      /** Sensor Y side size, in [mm] */
      double getPlaneYSize(int sensorID) {
        return _activeMap.at(sensorID)->getSize().coeff(1);
      };

      /** Sensor Z side size, in [mm] */
      double getPlaneZSize(int sensorID) {
        return _activeMap.at(sensorID)->getSize().coeff(2);
      };

      /** Sensor X side pixel pitch, in [mm] */
      double getPlaneXPitch(int sensorID) {
        return _activeMap.at(sensorID)->getPitch().first;
      };

      /** Sensor Y side pixel pitch, in [mm] */
      double getPlaneYPitch(int sensorID) {
        return _activeMap.at(sensorID)->getPitch().second;
      };

      /** Number of pixels in x-direction */
      int getPlaneNumberOfPixelsX(int sensorID) {
        return _activeMap.at(sensorID)->getNoPixels().first;
      };

      /** Number of pixels in y-direction */
      int getPlaneNumberOfPixelsY(int sensorID) {
        return _activeMap.at(sensorID)->getNoPixels().second;
      };

      /** Resolution of sensor in x-direction, in [mm] */ 
      double getPlaneXResolution(int sensorID) {
        return _activeMap.at(sensorID)->getResolution().first;
      };

      /** Resolution of sensor in y-direction, in [mm] */ 
      double getPlaneYResolution(int sensorID) {
        return _activeMap.at(sensorID)->getResolution().second;
      };

      /** Return the sensor's radiation length in [mm]*/
      double getPlaneRadiationLength(int sensorID) {
        return _activeMap.at(sensorID)->getRadLength();
      };

      /** Name of pixel geometry library */
      std::string geoLibName(int sensorID) {
        return _activeMap.at(sensorID)->getGeometry();
      };

      /** Get the plane's normal vector in global coordinates */
      Eigen::Vector3d getPlaneNormalVector(int);

      /** Get the plane's x-direction vector in global coordinates */
      Eigen::Vector3d getPlaneXVector(int);

      /** Get the plane's y-direction vector in global coordinates */
      Eigen::Vector3d getPlaneYVector(int);

      /** Vector of all sensor IDs */
      const std::vector<int> & sensorIDsVec() const { 
        return _sensorIDVec;
      };

      Eigen::Matrix3d rotationMatrixFromAngles(int sensorID);

      Eigen::Vector3d getOffsetVector(int sensorID);

      Eigen::Matrix3i getFlipMatrix(int sensorID);

      void writeGEARFile(std::string filename);

      virtual ~EUTelGeometryTelescopeGeoDescription();

      /** Initialize TGeo geometry
       * Establish access to TGeoManager, load geometry description file.
       *
       * @param tgeofilename name of a .root or .gdml file with valid geometry
       *
       * @see ROOT TGeoManager::Import
       */
      void initializeTGeoDescription(std::string const & tgeofilename);
      void initializeTGeoDescription(std::string const & geomName, bool dumpRoot);

      /** Retrieve the radiation length 
        * Radiation length between two points (each point in [mm,mm,mm]). 
        * The result is the radiation length in [mm] 
        */
      double getRadiationLengthBetweenPoints(Eigen::Vector3d const &startPt, Eigen::Vector3d const &endPt);

      double planeRadLengthGlobalIncidence(int planeID, Eigen::Vector3d incidenceDir);
      double planeRadLengthLocalIncidence(int planeID, Eigen::Vector3d incidenceDir);

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

      // This outputs the total percentage radiation length for the full
      // detector system.
//      float calculateTotalRadiationLengthAndWeights(
 //         const double startD[3], const double endD[3],
  //        std::map<const int, double> &, std::map<const int, double> &);
   //   double addKapton(std::map<const int, double> &mapSensor);

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
#endif /* EUTELGEOMETRYTELESCOPEGEODESCRIPTION_H */
