#ifndef EUTelGeoSupportClasses_h
#define EUTelGeoSupportClasses_h

#include "EUTelUtility.h"

#include "GEAR.h"
//#include "gearxml/GearXML.h"

#include "gearimpl/TrackerPlanesLayerLayoutImpl.h"
#include "gearimpl/TrackerPlanesParametersImpl.h"
/** Helper classes to represent the geometry reflected in the GEAR file in
 * EUTelescope .The description is based on layers and planes. A plane can 
 * either be an active or passive one - where active derives from passive. 
 * Even for the case of a one plane telescope, you will have one layer and 
 * one active (plane). The layer holds the alignment information for all
 * planes in it. Planes (i.e. active & passive) can further hold
 * corrections to this.
 *
 * The alignment of a layer is stored in two objects: the rotation matrix R
 * and a position vector V. A local coordinate x is transformed into a 
 * global X one by:
 *
 *			X = V + R*x
 *
 * If additional alignment constants are stored for a plane, denoted by V' 
 * and R' the transformation is as follows:
 *
 *			X = (V+V') + (R'*R)x
 *
 * Notice that the individual angles are NOT simply added, but the rotation
 * are carried out subsequently!
 */
namespace eutelescope {
  namespace geo {

    class EUTelActive;
    class EUTelPassive;

    class EUTelMaterial {
    public:
      EUTelMaterial(double A, double Z, double rho)
          : _A(A), _Z(Z), _density(rho) {}
      EUTelMaterial(EUTelMaterial const &other) = default;

      double _A, _Z, _density;
      // Radiation length in [mm]
      double _radLength = 0;
    };

    class EUTelLayer {
    public:
      EUTelLayer(int ID) : _ID(ID) {}

    private:
      int _ID;
      gear::TrackerPlanesLayerImpl *_gearLayerPtr;
      std::string _info = "";
      std::vector<std::unique_ptr<EUTelActive>> _activeVec;
      std::vector<std::unique_ptr<EUTelPassive>> _passiveVec;
      Eigen::Matrix3d _rotMatrix;
      Eigen::Vector3d _posVector;
      Eigen::Vector3d _posUncVector;
      Eigen::Vector3d _angleVector;
      Eigen::Vector3d _angleUncVector;

    public:
      int getID() const { return _ID; }

      void setGearPtr(gear::TrackerPlanesLayerImpl *ptr) {
        _gearLayerPtr = ptr;
      }
      gear::TrackerPlanesLayerImpl *getGearPtr() const { return _gearLayerPtr; }

      void setInfo(std::string const &s) { _info = s; }

      template <typename... Par> void addActivePlane(Par &&... par) {
        _activeVec.push_back(std::forward<Par>(par)...);
      }

      template <typename... Par> void addPassivePlane(Par &&... par) {
        _passiveVec.push_back(std::forward<Par>(par)...);
      }

      std::vector<std::unique_ptr<EUTelActive>> const &getActivePlanes() const {
        return _activeVec;
      }

      std::vector<std::unique_ptr<EUTelPassive>> const &
      getPassivePlanes() const {
        return _passiveVec;
      }

      void setPosition(double x, double y, double z) {
        _posVector = Eigen::Vector3d(x, y, z);
      }

      void setPosition(Eigen::Vector3d const &posVec) { _posVector = posVec; }

      void setRotation(double alpha, double beta, double gamma) {
        _angleVector = Eigen::Vector3d(alpha, beta, gamma);
        _rotMatrix = Utility::rotationMatrixFromAngles(alpha, beta, gamma);
      }

      void setRotationDeg(double alpha, double beta, double gamma) {
		    auto alphaRad = static_cast<double>(alpha*Utility::PI/180.);
		    auto betaRad = static_cast<double>(beta*Utility::PI/180.);
		    auto gammaRad = static_cast<double>(gamma*Utility::PI/180.);
        _angleVector = Eigen::Vector3d(alphaRad, betaRad, gammaRad);
        _rotMatrix = Utility::rotationMatrixFromAngles(alphaRad, betaRad, gammaRad);
      }

      void setPositionUnc(double xUnc, double yUnc, double zUnc) {
        _posUncVector = Eigen::Vector3d(xUnc, yUnc, zUnc);
      }

      void setRotationUnc(double alphaUnc, double betaUnc, double gammaUnc) {
        _angleUncVector = Eigen::Vector3d(alphaUnc, betaUnc, gammaUnc);
      }

      void reseatPosition(Eigen::Vector3d const &posVec);
      void reseatRotation(Eigen::Matrix3d const &rotVec);

      Eigen::Vector3d const &getPosVec() const { return _posVector; }
      Eigen::Vector3d const &getPosUncVec() const { return _posUncVector; }
      Eigen::Matrix3d const &getRotMat() const { return _rotMatrix; }
      Eigen::Vector3d const &getAngleVec() const { return _angleVector; }
    };

    class EUTelPassive {
    public:
      EUTelPassive(int ID, EUTelLayer *parent, EUTelMaterial &mat)
          : _ID(ID), _parentLayer(parent), _material(mat) {}
      virtual ~EUTelPassive() = default;

    protected:
      int _ID;

      std::string _info = "";

      EUTelLayer *_parentLayer;
      EUTelMaterial &_material;

      double _thickness;
      double _sizeX, _sizeY;

      Eigen::Vector3d _offVector;
      Eigen::Vector3d _deltaRotAnglesVector;
      Eigen::Matrix3d _deltaRotMatrix;
      Eigen::Vector3d _totalPosVector;
      Eigen::Matrix3d _totalRotMatrix;
      Eigen::Vector3d _globalRotationAngles;

    public:
      void setOffset(double xOff, double yOff, double zOff) {
        _offVector = Eigen::Vector3d(xOff, yOff, zOff);
        _totalPosVector = _parentLayer->getPosVec() + _offVector;
      }
      void setOffset(Eigen::Vector3d const &offset) {
        _offVector = offset;
        _totalPosVector = _parentLayer->getPosVec() + _offVector;
      }

      void setDeltaRotation(double alpha, double beta, double gamma) {
        _deltaRotAnglesVector = Eigen::Vector3d(alpha, beta, gamma);
        _deltaRotMatrix = Utility::rotationMatrixFromAngles(alpha, beta, gamma);
        _totalRotMatrix = _deltaRotMatrix * _parentLayer->getRotMat();
        _globalRotationAngles =
            Utility::getRotationAnglesFromMatrix(_totalRotMatrix);
      }

      void setDeltaRotation(Eigen::Matrix3d const &deltaRotMat) {
        _deltaRotAnglesVector =
            Utility::getRotationAnglesFromMatrix(deltaRotMat);
        _deltaRotMatrix = deltaRotMat;
        _totalRotMatrix = _deltaRotMatrix * _parentLayer->getRotMat();
        _globalRotationAngles =
            Utility::getRotationAnglesFromMatrix(_totalRotMatrix);
      }

      void setSize(double x, double y, double z) {
        _sizeX = x;
        _sizeY = y;
        _thickness = z;
      }

      void setInfo(std::string const &s) { _info = s; }

      Eigen::Vector3d const &getPosition() const { return _totalPosVector; }

      Eigen::Vector3d const &getOffset() const { return _offVector; }

      Eigen::Vector3d const &getGlobalRotationAngles() const {
        return _globalRotationAngles;
      }

      Eigen::Matrix3d const &getGlobalRotationMatrix() const {
        return _totalRotMatrix;
      }

      Eigen::Vector3d getSize() const {
        return Eigen::Vector3d(_sizeX, _sizeY, _thickness);
      }

      double getRadLength() const { return _material._radLength; }

      int getID() const { return _ID; }

      EUTelLayer const *getParent() const { return _parentLayer; }

      void alignPos(Eigen::Vector3d const &globalPos) {
        _parentLayer->reseatPosition(globalPos);
        _offVector = Eigen::Vector3d(0, 0, 0);
        _totalPosVector = _parentLayer->getPosVec();
      }

      void alignRot(Eigen::Matrix3d const &globalRotMat) {
        _parentLayer->reseatRotation(globalRotMat);
        _deltaRotAnglesVector = Eigen::Vector3d(0, 0, 0);
        _deltaRotMatrix = Eigen::Matrix3d::Identity();
        _totalRotMatrix = globalRotMat;
        _globalRotationAngles =
            Utility::getRotationAnglesFromMatrix(_totalRotMatrix);
      }
    };

    class EUTelActive : public EUTelPassive {
    public:
      EUTelActive(int ID, EUTelLayer *parent, EUTelMaterial &mat)
          : EUTelPassive(ID, parent, mat){};

    protected:
      double _xOffUnc = 0;
      double _yOffUnc = 0;
      double _zOffUnc = 0;

      double _xyDeltaRotUnc = 0;
      double _zxDeltaRotUnc = 0;
      double _zyDeltaRotUnc = 0;

      Eigen::Matrix2i _flipMatrix;

      double _pitchX, _pitchY;
      int _nPixelsX, _nPixelsY;
      double _resX, _resY;

      std::string _geometry = "CAST";

      bool _enabled = true;

      bool _axisFlipped = false;

    public:
      void setOffsetUnc(double xOffUnc, double yOffUnc, double zOffUnc) {
        _xOffUnc = xOffUnc;
        _yOffUnc = yOffUnc;
        _zOffUnc = zOffUnc;
      }

      void setDeltaRotationUnc(double alphaUnc, double betaUnc,
                               double gammaUnc) {
        _xyDeltaRotUnc = gammaUnc;
        _zyDeltaRotUnc = alphaUnc;
        _zxDeltaRotUnc = betaUnc;
      }

      void setFlips(int f1, int f2, int f3, int f4) {
        _flipMatrix << f1, f2, f3, f4;
      }

      void setPitch(double xPitch, double yPitch) {
        _pitchX = xPitch;
        _pitchY = yPitch;
      }

      void setNoPixels(int noX, int noY) {
        _nPixelsX = noX;
        _nPixelsY = noY;
      }

      void setResolution(double xRes, double yRes) {
        _resX = xRes;
        _resY = yRes;
      }

      void setEnabled(bool enabled) { _enabled = enabled; }
      void setGeometry(std::string const &g) { _geometry = g; }

      Eigen::Matrix2i const &getFlipMatrix() const { return _flipMatrix; }

      std::pair<double, double> getResolution() const {
        return std::make_pair(_resX, _resY);
      }

      std::pair<int, int> getNoPixels() const {
        return std::make_pair(_nPixelsX, _nPixelsY);
      }

      std::pair<double, double> getPitch() const {
        return std::make_pair(_pitchX, _pitchY);
      }

      std::string getGeometry() const { return _geometry; }
    };
  }
} // namespaces
#endif
