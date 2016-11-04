#ifndef FEI4Single400uEdge400uEdge_h 
#define FEI4Single400uEdge400uEdge_h	

  /** @class FEI4Single400uEdge
	* This class is the implementation of  @class EUTelGenericPixGeoDescr
	* for a FEI4 layout with edge pixels which are 400 microns long, the
	* other properties are: 80 x 336 pixels, 250 x 50 microns**2 size
	* with exception of the edge pixels (X=0,79) which are longer.
    */

//EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"

//ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

namespace eutelescope {
namespace geo {

class FEI4Single400uEdge : public EUTelGenericPixGeoDescr {
	
	public:
		FEI4Single400uEdge();
		~FEI4Single400uEdge();

		void createRootDescr(char const *);
		std::string getPixName(int, int);
		std::pair<int, int> getPixIndex(char const *);

	protected:
		TGeoMaterial* matSi;
		TGeoMedium* Si;
		TGeoVolume* plane;

};

extern "C"
{
	EUTelGenericPixGeoDescr* maker();
}

} //namespace geo
} //namespace eutelescope

#endif	//FEI4SINGLE_H
