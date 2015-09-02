#include "FEI4FourChip.h"
#include <iostream>

namespace eutelescope {
namespace geo {

FEI4FourChip::FEI4FourChip(): EUTelGenericPixGeoDescr(	40.4, 35.18, 0.025,	//size X, Y, Z
							0, 159, 0, 671,		//min max X,Y
							9.3660734 )		//rad length					
{
	//Create the material for the sensor
	matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, -_radLength, 45.753206 );
	Si = new TGeoMedium("FEI4Silicon",1, matSi);

	/* Make a box for the sensitive area
	MakeBox takes the half of size in mm as arguments */
	plane = _tGeoManager->MakeBox( "sensarea_fei4four", Si, 20.2, 17.59, 0.0125 );

	//Create volumes for the double chipe
	TGeoVolume* doublechip = _tGeoManager->MakeBox("fei4double",Si, 20.2, 8.4, 0.0125);

	//Create volumes for the different regions
	TGeoVolume* normalRegion = _tGeoManager->MakeBox("fei4normreg",Si, 9.875, 8.4, 0.0125);
	TGeoVolume* centreRegion = _tGeoManager->MakeBox("fei4centreg",Si, 0.45, 8.4, 0.0125);

	//Divide the regions to create pixels
	TGeoVolume* normalCol =  normalRegion->Divide("col", 1 , 79, 0 , 1, 0, "N"); 
	normalCol->Divide("pixel", 2 , 336, 0 , 1, 0, "N");

	TGeoVolume* centreCol = centreRegion->Divide("col", 1 , 2, 0 , 1, 0, "N"); 
	centreCol->Divide("pixel", 2 , 336, 0 , 1, 0, "N"); 

	//And place them to make a doublechip
	doublechip->AddNode(normalRegion, 1,  new TGeoTranslation(-10.325,0,0));
	doublechip->AddNode(centreRegion, 1);
	doublechip->AddNode(normalRegion, 2,  new TGeoTranslation(10.325,0,0));

	//Place two double chips for a four chip module
	plane->AddNode(doublechip, 1, new TGeoTranslation(0,-9.19,0));
	plane->AddNode(doublechip, 2, new TGeoTranslation(0,9.19,0));
}

FEI4FourChip::~FEI4FourChip()
{
	//delete matSi;
	//delete Si;
}

void  FEI4FourChip::createRootDescr(char const * planeVolume)
{
	//Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
	TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
	//Finaly add the sensitive area to the plane
	topplane->AddNode(plane, 1);
}

std::string FEI4FourChip::getPixName(int x, int y)
{
	char buffer [100];
	int doublechip = 1;

	//If y index is larger than 335 then we are on double chip two (the upper one)
	if(y > 335)
	{
		doublechip = 2;
		y -= 336;
	}

	if (x < 79 )
	{
		//in normalRegion_1
		snprintf( buffer, 100, "/sensarea_fei4four_1/fei4double_%d/fei4normreg_1/col_%d/pixel_%d", doublechip, x+1, y+1);
	}

	else if ( x == 79 )
	{
		//in centreRegion_1, col_1
		snprintf( buffer, 100, "/sensarea_fei4four_1/fei4double_%d/fei4centreg_1/col_1/pixel_%d", doublechip, y+1);
	}
	else if( x == 80 )
	{
		//in centreRegion_1, col_2
		snprintf( buffer, 100, "/sensarea_fei4four_1/fei4double_%d/fei4centreg_1/col_2/pixel_%d", doublechip, y+1);
	}
	else
	{
		//in normalRegion_2
		snprintf( buffer, 100, "/sensarea_fei4four_1/fei4double_%d/fei4normreg_2/col_%d/pixel_%d", doublechip, x-80, y+1);
	}
	//Return the full path
	return std::string( buffer ); 
}

//TODO: parse the path to a pixel number!
std::pair<int, int>  FEI4FourChip::getPixIndex(char const*){return std::make_pair(0,0); }

EUTelGenericPixGeoDescr* maker()
{
	FEI4FourChip* mPixGeoDescr = new FEI4FourChip();
	return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
}

} //namespace geo
} //namespace eutelescope

