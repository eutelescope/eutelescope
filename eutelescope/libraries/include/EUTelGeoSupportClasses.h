#ifndef EUTelGeoSupportClasses_h
#define EUTelGeoSupportClasses_h
namespace eutelescope {
namespace geo{

class EUTelActive;
class EUTelPassive;

class EUTelMaterial {


};

class EUTelLayer {
  public:

	EUTelLayer(int ID): _ID(ID){}

	int _ID;

	double _xPos = 0;
	double _yPos = 0;
	double _zPos = 0;

	double _xyRot = 0;
	double _zxRot = 0;
	double _zyRot = 0;

	double _xPosUnc = 0;
	double _yPosUnc = 0;
	double _zPosUnc = 0;

	double _xyRotUnc = 0;
	double _zxRotUnc = 0;
	double _zyRotUnc = 0;

	std::string _info = "";

	std::vector<EUTelActive> _activeVec;
	std::vector<EUTelPassive> _passiveVec;
};

class EUTelPassive {
  public:

	EUTelPassive(int ID, EUTelLayer& parent, EUTelMaterial& mat): _ID(ID), _parentLayer(parent), _material(mat){}
	virtual ~EUTelPassive() = default;

	int _ID;
	double _xOff = 0;
	double _yOff = 0;
	double _zOff = 0;

	double _xyDeltaRot = 0;
	double _zxDeltaRot = 0;
	double _zyDeltaRot = 0;

	std::string _info = "";

	EUTelLayer& _parentLayer;
	EUTelMaterial& _material;

	double _thickness;
	double _sizeX, _sizeY;

};

class EUTelActive: public EUTelPassive {
  public:
	EUTelActive(int ID, EUTelLayer& parent, EUTelMaterial& mat): EUTelPassive(ID, parent, mat){};

	double _xOffUnc = 0;
	double _yOffUnc = 0;
	double _zOffUnc = 0;

	double _xyDeltaRotUnc = 0;
	double _zxDeltaRotUnc = 0;
	double _zyDeltaRotUnc = 0;

	int _flip1 = 1;
	int _flip2 = 0;
	int _flip3 = 0;
	int _flip4 = 1;

	double _pitchX, _pitchY;
	double _nPixelsX, _nPixelsY;
	double _resX, _resY;

	std::string _geometry = "CAST";	

	bool _enabled = true;
};


}} //namespaces
#endif
