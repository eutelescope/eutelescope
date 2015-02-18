#ifndef FEI4FOURCHIP_H
#define	FEI4FOURCHIP_H

  /** @class FEI4Single
	* This class is the implementation of  @class EUTelGenericPixGeoDescr
	* for a for chip FEI4 layout.
	* Two double chips make up a four chip module, the two double chips are
	* seperated by 1.56 mm dead area. 
	* One double chip contains 336 pixels in y direction featuring 50 microns
	* pitch. In total 160 pixels are along x the x direction, all except the 
	* two pixels in the centre (79, 80 where numbering starts from from 0) are
	* 250 microns long. Each of the two centre pixels are 450 microns long.
    */

//EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"

//ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

namespace eutelescope {
namespace geo {

class FEI4FourChip : public EUTelGenericPixGeoDescr {
	
	public:
		FEI4FourChip();
		~FEI4FourChip();

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

}//nanespace geo
} //namespace eutelescope

#endif	//FEI4FOURCHIP_H
