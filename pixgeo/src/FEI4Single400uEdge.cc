#include "FEI4Single400uEdge.h"

namespace eutelescope {
namespace geo {

FEI4Single400uEdge::FEI4Single400uEdge(): EUTelGenericPixGeoDescr(	20.30, 16.8, 0.025,	//size X, Y, Z
							0, 79, 0, 335,		//min max X,Y
							9.3660734 )		//rad length					
{
	//Create the material for the sensor
	matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, -_radLength, 45.753206 );
	Si = new TGeoMedium("FEI4Silicon",1, matSi);

	/* Make a box for the sensitive area
	Size is: x=2*400+78*250=20300 microns and y=336*50=16800 microns
	MakeBox takes the half of those values in mm as arguments */
	plane = _tGeoManager->MakeBox( "sensarea_fei4", Si, 10.15, 8.4, 0.0125 );

	//Create volumes for the centre and edge regions
	TGeoVolume* centreregion = _tGeoManager->MakeBox("fei4centreregion",Si, 9.75, 8.4, 0.0125);
	TGeoVolume* edgeregion   = _tGeoManager->MakeBox("fei4edgeregion"  ,Si, 0.2 , 8.4, 0.0125);

	//Divide the regions to create pixels
 	edgeregion->Divide("fei4edgepixel",   2, 336, 0, 1, 0, "N"); 
	TGeoVolume* centrerow = centreregion->Divide("fei4centrerow", 2, 336, 0, 1, 0, "N");
	centrerow ->Divide("fei4centrepixel", 1,  78, 0, 1, 0, "N"); 

        //And place them to make a singlechip
	plane->AddNode(centreregion, 1, new TGeoTranslation( 0.00 , 0 , 0) );
	plane->AddNode(edgeregion,   1, new TGeoTranslation(-9.95 , 0 , 0) );
	plane->AddNode(edgeregion,   2, new TGeoTranslation( 9.95 , 0 , 0) );

}

FEI4Single400uEdge::~FEI4Single400uEdge()
{
	//delete matSi;
	//delete Si;
}

void  FEI4Single400uEdge::createRootDescr(char const * planeVolume)
{
	//Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
	TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
	//Finaly add the sensitive area to the plane
	topplane->AddNode(plane, 1);
}

std::string FEI4Single400uEdge::getPixName(int x , int y)
{
	char buffer [100];

	//since pixel 0|0 is located on the upper left corner we have to correct y by 335-y+1 
	//(one for the offset in TGeo which starts counting at 1)
	if (x == 0 )
	{
		snprintf( buffer, 100, "/sensarea_fei4_1/fei4edgeregion_1/fei4edgepixel_%d", 336-y);
	}

	else if ( x == 79 )
	{
		snprintf( buffer, 100, "/sensarea_fei4_1/fei4edgeregion_2/fei4edgepixel_%d", 336-y);
	}
	if(x > 0 && x < 79 )
	{
	
		snprintf( buffer, 100, "/sensarea_fei4_1/fei4centreregion_1/fei4centrerow_%d/fei4centrepixel_%d", 336-y, x);
	}
	//Return the full path
	return std::string( buffer ); 
}

//TODO: parse the path to a pixel number!
std::pair<int, int>  FEI4Single400uEdge::getPixIndex(char const*){return std::make_pair(0,0); }

EUTelGenericPixGeoDescr* maker()
{
	FEI4Single400uEdge* mPixGeoDescr = new FEI4Single400uEdge();
	return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
}

} //namespace geo
} //namespace eutelescope

