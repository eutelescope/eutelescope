#include "FEI4Single.h"

namespace eutelescope {
namespace geo {

FEI4Single::FEI4Single(): EUTelGenericPixGeoDescr(	20.00, 16.8, 0.025,	//size X, Y, Z
							0, 79, 0, 335,		//min max X,Y
							9.3660734 )		//rad length					
{
	//Create the material for the sensor
	matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, -_radLength, 45.753206 );
	Si = new TGeoMedium("FEI4Silicon",1, matSi);

	/* Make a box for the sensitive area
	Size is: x=2*400+78*250=20300 microns and y=336*50=16800 microns
	MakeBox takes the half of those values in mm as arguments */
	plane = _tGeoManager->MakeBox( "sns_fei4", Si, 10.0, 8.4, 0.0125 );

	auto row = plane->Divide("row", 2, 336, 0, 1, 0, "N"); 
	row->Divide("col", 1,  80, 0, 1, 0, "N"); 
}

FEI4Single::~FEI4Single()
{
	//delete matSi;
	//delete Si;
}

void  FEI4Single::createRootDescr(char const * planeVolume)
{
	//Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
	TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
	//Finaly add the sensitive area to the plane
	topplane->AddNode(plane, 1);
}

std::string FEI4Single::getPixName(int x , int y)
{
	char buffer [100];
	//since pixel 0|0 is located on the upper left corner we have to correct y by 335-y+1 
	//(one for the offset in TGeo which starts counting at 1)
	snprintf( buffer, 100, "/sns_fei4_1/row_%d/col_%d", 336-y, x+1);
	//Return the full path
	return std::string( buffer ); 
}

//TODO: parse the path to a pixel number!
std::pair<int, int>  FEI4Single::getPixIndex(char const*){return std::make_pair(0,0); }

EUTelGenericPixGeoDescr* maker()
{
	FEI4Single* mPixGeoDescr = new FEI4Single();
	return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
}

} //namespace geo
} //namespace eutelescope

