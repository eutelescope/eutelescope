#ifndef FEI4DOUBLE_H
#define	FEI4DOUBLE_H

  /** @class FEI4Double
	* This class is the implementation of  @class EUTelGenericPixGeoDescr
	* for a FEI4 double chip layout with  default edge pixels (250 microns)
	* and two prolonged pixels in between the two single chips of 450 micron
	* size.
	* 160 x 672 pixels, in x-direction: 79 times 250 mu, 2 times 450 mu and 
	* again 79 times 250 mu. In y-direction 336 times 50 mu
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

extern "C"
{
	EUTelGenericPixGeoDescr* maker();
}

} //namespace geo
} //namespace eutelescope

#endif	//FEI4DOUBLE_H
