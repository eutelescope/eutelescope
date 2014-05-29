#include "Mimosa26GeoDescr.h"

namespace eutelescope {
namespace geo {

Mimosa26GeoDescr::Mimosa26GeoDescr(): EUTelGenericPixGeoDescr(	21.2, 10.6, 0.02,		//size X, Y, Z
																0, 1151, 0, 575,		//min max X,Y
																93.660734 )				//rad length
{
	//Create the material for the sensor
	matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, _radLength, 45.753206 );
	Si = new TGeoMedium("MimosaSilicon",1, matSi);

	//Create a plane for the sensitive area
	plane = _tGeoManager->MakeBox( "sensarea_mimosa", Si, 10.6, 5.3, 0.01 );
	//Divide the regions to create pixels
  	TGeoVolume* row = plane->Divide("mimorow", 1 , 1152 , 0 , 1, 0, "N"); 
	row->Divide("mimopixel", 2 , 576, 0 , 1, 0, "N");

}

Mimosa26GeoDescr::~ Mimosa26GeoDescr()
{
	//It appears that ROOT will take ownership and delete that stuff! 
	//delete matSi,
	//delete Si;
}

void  Mimosa26GeoDescr::createRootDescr(char const * planeVolume)
{
	//Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
	TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
	//Add the sensitive area to the plane
//	topplane->AddNode(plane, 1);
     
        float xoffset = 18.4/2.;
        float yoffset = 18.4/2.;
        topplane->AddNode(plane, 1, new TGeoTranslation( 0.00 + xoffset, 0 + yoffset, 0) );
}

std::string Mimosa26GeoDescr::getPixName(int x , int y)
{
	char buffer [100];
	//return path to the pixel, don't forget to shift indices by +1+
	snprintf( buffer, 100, "/sensarea_mimosa_1/mimorow_%d/mimopixel_%d", x+1, y+1);
	return std::string( buffer ); 
}
	/*TODO*/ std::pair<int, int>  Mimosa26GeoDescr::getPixIndex(char const*){return std::make_pair(0,0); }

} //namespace geo
} //namespace eutelescope

