#ifndef FEI4Single_h
#define FEI4Single_h

  /** @class FEI4Single
	* This class is the implementation of  @class EUTelGenericPixGeoDescr
	* for a default FEI4 layout 
    */

//EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"

//ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

namespace eutelescope {
namespace geo {

class FEI4Single : public EUTelGenericPixGeoDescr {
	
	public:
		FEI4Single();
		~FEI4Single();

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
