#ifndef EUTELEXCLUDEPLANES_H
#define EUTELEXCLUDEPLANES_H

#include "EUTelGeometryTelescopeGeoDescription.h"

namespace eutelescope 
{

class EUTelExcludedPlanes
{
	public: 
		static std::vector<int> getExcludedPlaneIDVec(std::vector<int> const& planesToExclude);
		static std::vector<int> getExcludedPlaneIDVecZSorted(std::vector<int> const& planesToExclude);
  	
	private:
		static std::vector<int> getRelativeComplementSet(std::vector<int> const& planesToExclude);
		static bool sortSensorIDByZPos(int i, int j);
		EUTelExcludedPlanes();
};

}
#endif
