#include "EUTelExcludedPlanes.h"
#include <algorithm>

namespace eutelescope 
{
EUTelExcludedPlanes::EUTelExcludedPlanes(){
    _senNoDeadMaterial = getSensorsNotDeadMaterial();
}


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
    std::vector<int> resultVec;
    const std::vector<int>& sensorIDsVec = geo::gGeometry().sensorIDsVec();

    ///We look for dead material and excluded sensors. Should split this in a better way.
    for(std::vector<int>::const_iterator it = sensorIDsVec.begin(); it != sensorIDsVec.end(); it++)
    {
        int currentSensorID = *it;
        std::vector<int>::const_iterator pos = std::find(planesToExclude.begin(), planesToExclude.end(), currentSensorID);
        if(pos == planesToExclude.end() and currentSensorID <= 99 or currentSensorID ==271){ ///Not all end statements are create equal! :-)
            resultVec.push_back(currentSensorID);
        }
    }

	return resultVec;
}
std::vector<int> EUTelExcludedPlanes::getSensorsNotDeadMaterial()
{
    std::vector<int> resultVec;
    const std::vector<int>& sensorIDsVec = geo::gGeometry().sensorIDsVec();

    for(std::vector<int>::const_iterator it = sensorIDsVec.begin(); it != sensorIDsVec.end(); it++){
        int currentSensorID = *it;
        if(currentSensorID <= 99 or currentSensorID == 271){ 
            resultVec.push_back(currentSensorID);
        }
    }
	return resultVec;
}

void EUTelExcludedPlanes::setRelativeComplementSet(std::vector<int> const& planesToExclude)
{
    _senInc = getRelativeComplementSet(planesToExclude);
}
void EUTelExcludedPlanes::setPlaneInc(std::vector<int> const& plaInc){
    _senInc = plaInc;
}
bool EUTelExcludedPlanes::sortSensorIDByZPos(int i, int j)
{
	return geo::gGeometry().siPlaneZPosition(i) < geo::gGeometry().siPlaneZPosition(j);
}
    std::vector<int>  EUTelExcludedPlanes::_senInc;
    std::vector<int>  EUTelExcludedPlanes::_senNoDeadMaterial;

} //namespace eutelescope
