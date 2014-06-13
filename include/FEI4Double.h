#ifndef FEI4DOUBLE_H
#define	FEI4DOUBLE_H

  /** @class FEI4Double
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

class FEI4Double : public EUTelGenericPixGeoDescr {
	
	public:
		FEI4Double();
		~FEI4Double();

		void createRootDescr(char const *);
		std::string getPixName(int, int);
		std::pair<int, int> getPixIndex(char const *);

	protected:
		TGeoMaterial* matSi;
		TGeoMedium* Si;
		TGeoVolume* plane;

};

} //namespace geo
} //namespace eutelescope

#endif	//FEI4DOUBLE_H
