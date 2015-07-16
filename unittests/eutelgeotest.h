#ifndef EUTELGEOTEST_H
#define EUTELGEOTEST_H

#include "EUTelGeometryTelescopeGeoDescription.h"

#include "gearxml/GearXML.h"
#include "gear/GearMgr.h"
#include "gear/GEAR.h"

class eutelgeotest {

public:
	eutelgeotest();

protected:
	gear::GearMgr* gearManager;
};
#endif
