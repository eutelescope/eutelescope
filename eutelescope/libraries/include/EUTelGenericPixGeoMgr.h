#ifndef EUTELGENERICPIXGEOMGR_H
#define	EUTELGENERICPIXGEOMGR_H

//STL
#include <map>
#include <string>
#include <utility>

//EUTelescope
#include "EUTelGenericPixGeoDescr.h"

//ROOT includes
#include "TGeoManager.h"

namespace eutelescope {
namespace geo {

/** @class EUTelGenericPixGeoMgr
 * This class holds is the interface between the pixel geoemtry
 * description and the actual telescope set up. Meaning that it 
 * provides methods to obtain information about a specific plane
 * and it's geoemtry. The planeID is the crucial index used to 
 * address any plane.
 * This class also provides functionality to call the method of
 * the pixel geometry descriptions which load the geometry into
 * TGeo planes.
 * This class should not be exposed to the user, it should only
 * be used in the EUTelGeometryTelescopeGeoDescription instance.
 */
class EUTelGenericPixGeoMgr {

public:

	/** Default constructor */
	EUTelGenericPixGeoMgr();
	/** Destructor */
	~EUTelGenericPixGeoMgr();

	/** Method to add a plane to the pixel geo manager.
	 *  It is the method to which will call the factory method 
	 *  of a implementation of the EUTelGenericPixGeoDescr.
	 *	It is used by the EUTelGeometryTelescopeGeoDescription
	 *  to add pixel layouts to planes.
	 *
	 *  @param planeID The ID of the plane
	 *
	 *  @param geoName The name of the shared library holding
	 *  the pixel layout in form of a EUTelGenericPixGeoDescr
	 *
	 *  @param planeVolume The path of the plane into which the
	 *  geometry shall be loaded
	 */
	void addPlane(int planeID, std::string geoName, std::string planeVolume);

	void addCastedPlane(int planeID, int xPixel, int yPixel, double xSize, double ySize, double zSize, double radLength, std::string planeVolume);

	/** Method to get the EUTelGenericPixGeoDescr of a plane. 
	 * 
	 *  @param planeID The plane of which the EUTelGenericPixGeoDescr
	 *  is desired.
	 */
	EUTelGenericPixGeoDescr* getPixGeoDescr(int planeID);

protected:	

	/** Map of the geo library name and the actual pointer to the instance of it. */
	std::map<std::string, EUTelGenericPixGeoDescr*> _pixelDescriptions;	

	/** Map of the casted geoemtry and the actual pointer to the instance of it. */
	std::map<std::string, EUTelGenericPixGeoDescr*> _castedDescriptions;	
	
	/** Map of the planeID and corresponding EUTelGenericPixGeoDescr* */
	std::map<int, EUTelGenericPixGeoDescr* > _geoDescriptions;

}; //class EUTelGenericGeoMgr

} //namespace geo
} //namespace eutelescope

#endif	//EUTELGENERICPIXGEOMGR_H
