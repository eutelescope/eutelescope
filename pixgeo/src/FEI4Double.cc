#include "FEI4Double.h"

namespace eutelescope {
namespace geo {

FEI4Double::FEI4Double(): EUTelGenericPixGeoDescr(	40.40, 16.8, 0.025,	//size X, Y, Z
							0, 159, 0, 335,		//min max X,Y
							9.3660734 )		//rad length					
{
	//Create the material for the sensor
	matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, -_radLength, 45.753206 );
	Si = new TGeoMedium("FEI4Silicon",1, matSi);

	/* Make a box for the sensitive area
	Size is: x=2*78*250+2*450=40400 microns and y=336*50=16800 microns
	MakeBox takes the half of those values in mm as arguments */
	plane = _tGeoManager->MakeBox( "sensarea_fei4d", Si, 20.20, 8.4, 0.0125 );

	//Create volumes for the centre and edge region(s)
	TGeoVolume* centreregion = _tGeoManager->MakeBox("fei4dcentreregion",Si, 0.45, 8.4, 0.0125);
	TGeoVolume* edgeregion = _tGeoManager->MakeBox("fei4dedgeregion",Si, 9.875, 8.4, 0.0125);

	//Divide the regions to create pixels
	TGeoVolume* edgerow = edgeregion->Divide("fei4dedgerow", 2 , 336, 0 , 1, 0, "N"); 
	edgerow->Divide("fei4dedgepixel", 1, 79, 0, 1, 0, "N");
	TGeoVolume* centrerow = centreregion->Divide("fei4dcentrerow", 2 , 336, 0 , 1, 0, "N"); 
	centrerow->Divide("fei4dcentrepixel", 1 , 2, 0 , 1, 0, "N"); 

	//And place them to make a doublechip
	plane->AddNode(centreregion, 1);
	plane->AddNode(edgeregion, 1, new TGeoTranslation(-10.325,0,0));
	plane->AddNode(edgeregion, 2, new TGeoTranslation(10.325,0,0));
}

FEI4Double::~FEI4Double()
{
	//deletion of medium and material done by root
}

void  FEI4Double::createRootDescr(char const * planeVolume)
{
	//Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
	TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
	//Finaly add the sensitive area to the plane
	topplane->AddNode(plane, 1);
}

std::string FEI4Double::getPixName(int x , int y)
{
	char buffer [200];

	//since pixel 0|0 is located on the upper left corner we have to correct y by 335-y+1 
	//(one for the offset in TGeo which starts counting at 1)
	if (x < 79 )
	{
		snprintf( buffer, 200, "/sensarea_fei4d_1/fei4dedgeregion_1/fei4dedgerow_%d/fei4dedgepixel_%d", 336-y, x+1);
	}

	else if ( x == 79)
	{
		snprintf( buffer, 200, "/sensarea_fei4d_1/fei4dcentreregion_1/fei4dcentrerow_%d/fei4dcentrepixel_1", 336-y);
	}

	else if ( x == 80)
	{
		snprintf( buffer, 200, "/sensarea_fei4d_1/fei4dcentreregion_1/fei4dcentrerow_%d/fei4dcentrepixel_2", 336-y);
	}

	else
	{
	
		snprintf( buffer, 200, "/sensarea_fei4d_1/fei4dedgeregion_2/fei4dedgerow_%d/fei4dedgepixel_%d", 336-y, x-80);
	}
	//Return the full path
	return std::string( buffer ); 
}

//TODO: parse the path to a pixel number!
std::pair<int, int>  FEI4Double::getPixIndex(char const*){return std::make_pair(0,0); }

EUTelGenericPixGeoDescr* maker()
{
	FEI4Double* mPixGeoDescr = new FEI4Double();
	return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
}

} //namespace geo
} //namespace eutelescope

