#include "eutelgeotest.h"

namespace eugeo = eutelescope::geo;

eutelgeotest::eutelgeotest() {
	gear::GearXML gearXML( "unitTestGear1.xml" ) ;
	gearManager = gearXML.createGearMgr();
	std::string name("test.root");
	eugeo::gGeometry( gearManager ).initializeTGeoDescription(name,false);
}
