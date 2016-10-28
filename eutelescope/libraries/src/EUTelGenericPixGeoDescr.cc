#include "EUTelGenericPixGeoDescr.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

using namespace eutelescope;
using namespace geo;


EUTelGenericPixGeoDescr::EUTelGenericPixGeoDescr(double sizeX, double sizeY, double sizeZ, int minX, int maxX, int minY, int maxY, double radLen): 
				_tGeoManager( gGeometry()._geoManager.get() ),
				_sizeSensitiveAreaX(sizeX),
				_sizeSensitiveAreaY(sizeY),
				_sizeSensitiveAreaZ(sizeZ),
				_minIndexX(minX),
				_minIndexY(minY),
				_maxIndexX(maxX),
				_maxIndexY(maxY),
				_radLength(radLen)
{}

