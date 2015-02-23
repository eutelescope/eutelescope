#include "EUTelExcludedPlanes.h"
#include <algorithm>

namespace eutelescope 
{

std::vector<int> EUTelExcludedPlanes::getExcludedPlaneIDVec(std::vector<int> const& planesToExclude)
{
	std::vector<int> resultVec = getRelativeComplementSet(planesToExclude);
	std::sort(resultVec.begin(), resultVec.end());
	return resultVec;
}



std::vector<int> EUTelExcludedPlanes::getExcludedPlaneIDVecZSorted(std::vector<int> const& planesToExclude)
{
	std::vector<int> resultVec = getRelativeComplementSet(planesToExclude);
	std::sort(resultVec.begin(), resultVec.end(), sortSensorIDByZPos);
	return resultVec;
}

/* Takes a vector of planes to exclude, retrieves the complete sensorID
 * vec from the geometry class and computes the relative complement set.
 * I.e. all sensors which are not excluded */
std::vector<int> EUTelExcludedPlanes::getRelativeComplementSet(std::vector<int> const& planesToExclude)
{
        EVENT::IntVec const& sensorIDsVec = geo::gGeometry().sensorIDsVec();
	std::vector<int> resultVec;

	for(EVENT::IntVec::const_iterator it = sensorIDsVec.begin(); it != sensorIDsVec.end(); it++)
	{
		int currentSensorID = *it;
		std::vector<int>::const_iterator pos = std::find(planesToExclude.begin(), planesToExclude.end(), currentSensorID);
		if(pos == sensorIDsVec.end()) resultVec.push_back(currentSensorID);
	}

	return resultVec;
}

bool EUTelExcludedPlanes::sortSensorIDByZPos(int i, int j)
{
	return geo::gGeometry().siPlaneZPosition(i) < geo::gGeometry().siPlaneZPosition(j);
}

} //namespace eutelescope
