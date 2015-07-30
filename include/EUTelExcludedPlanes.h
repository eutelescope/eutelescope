#ifndef EUTELEXCLUDEPLANES_H
#define EUTELEXCLUDEPLANES_H

#include "EUTelGeometryTelescopeGeoDescription.h"

namespace eutelescope 
{

class EUTelExcludedPlanes
{
public: 
    EUTelExcludedPlanes();
    static std::vector<int> getExcludedPlaneIDVec(std::vector<int> const& planesToExclude);
    static std::vector<int> getExcludedPlaneIDVecZSorted(std::vector<int> const& planesToExclude);
    static void setRelativeComplementSet(std::vector<int> const& planesToExclude);
    static std::vector<int> getSensorsNotDeadMaterial();
    static void setPlaneInc(std::vector<int> const& plaInc);
    static  std::vector<int>  _senInc;
    static  std::vector<int>  _senNoDeadMaterial;

private:
    static std::vector<int> getRelativeComplementSet(std::vector<int> const& planesToExclude);
    static bool sortSensorIDByZPos(int i, int j);

};

}
#endif
