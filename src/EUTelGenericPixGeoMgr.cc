//STL
#include <map>
#include <sstream>
#include <stdexcept>

//System
#include <dlfcn.h>

//EUTelescope
#include "EUTelGenericPixGeoMgr.h"

// MARLIN
#include "marlin/Global.h"
#include "marlin/VerbosityLevels.h"

//Geometry implementations
#include "GEARPixGeoDescr.h"

//ROOT includes
#include "TGeoBBox.h"
#include "TObjArray.h"

using namespace eutelescope;
using namespace geo;

//Default constructor
EUTelGenericPixGeoMgr::EUTelGenericPixGeoMgr() {}

//Destructor
EUTelGenericPixGeoMgr::~EUTelGenericPixGeoMgr() 
{
	//Clean up by deleting all EUTelGenericPixGeoDescr 
	std::map<std::string, EUTelGenericPixGeoDescr*>::iterator it;
	for(it = _pixelDescriptions.begin(); it != _pixelDescriptions.end(); ++it )
	{
		streamlog_out( MESSAGE3 ) << "Deleting " << (*it).first << std::endl;
		delete (*it).second;
	}
}

//TODO: comments
void EUTelGenericPixGeoMgr::addCastedPlane(int planeID, int xPixel, int yPixel, double xSize, double ySize, double zSize, double radLength, std::string planeVolume)
{
	EUTelGenericPixGeoDescr* pixgeodescrptr = NULL;
	int xSizeMap  = static_cast<int>(1000*xSize+0.5);  
	int ySizeMap  = static_cast<int>(1000*ySize+0.5); 
	int zSizeMap  = static_cast<int>(1000*zSize+0.5); 

	std::stringstream stream;
	stream << xPixel << yPixel << xSizeMap << ySizeMap << zSizeMap;
	std::string name = stream.str();

        std::map<std::string, EUTelGenericPixGeoDescr*>::iterator it;
	it = _castedDescriptions.find(name);
	
	if( it!=_castedDescriptions.end() )
	{
		//if it is, use it!
		streamlog_out( MESSAGE3 )  << "Found " << name << ", using it" << std::endl;
		pixgeodescrptr = (*it).second;
		streamlog_out( MESSAGE3 ) << "Inserting " << name << " into map" << std::endl;
		_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
	}
	else
	{
		streamlog_out( MESSAGE3 ) << "Didnt find " << name << " yet, thus creating" << std::endl;
		pixgeodescrptr = new GEARPixGeoDescr(xPixel, yPixel, xSize, ySize, zSize, radLength);
		streamlog_out( MESSAGE3 ) << "Inserting " << name << " into map" << std::endl;
		_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
		_castedDescriptions.insert( std::make_pair(name, pixgeodescrptr) );
	}
	
	streamlog_out( MESSAGE3 )  << "Adding plane: " << planeID << " with geoLibName: " << name << " in volume " << planeVolume << std::endl;
	pixgeodescrptr->createRootDescr(planeVolume);
}


void EUTelGenericPixGeoMgr::addPlane(int planeID, std::string geoName, std::string planeVolume)
{
	EUTelGenericPixGeoDescr* pixgeodescrptr = NULL;
	std::map<std::string, EUTelGenericPixGeoDescr*>::iterator it;
	
	//Check if the geoemtry is loaded already (e.g. for some other plane)
	it = _pixelDescriptions.find(geoName);	
	if( it != _pixelDescriptions.end() )
	{
		//if it is, use it!
		streamlog_out( MESSAGE3 )  << "Found " << geoName << ", using it" << std::endl;
		pixgeodescrptr = (*it).second;
		streamlog_out( MESSAGE3 ) << "Inserting " << planeID << " into map" << std::endl;
		_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
	}	
	//Otherwise we load it
	else
	{
		streamlog_out( MESSAGE3 ) << "Didnt find " << geoName << " yet, thus creating" << std::endl;
		
		//CMake will call shared libraried "lib"+LibraryName.so
		std::string libName = std:: string("lib").append(geoName);
		
		//Load shared library, be sure to export the path of the lib to LD_LIBRARY_PATH!
		void *hndl = dlopen(libName.c_str(), RTLD_NOW);
		if(hndl == NULL)
		{
			streamlog_out(ERROR7) << "Loading of " << libName << " failed: " << dlerror() << std::endl;
			throw std::runtime_error("dlopen could not open shared libraray");
		}

		try
		{
			//Use maker() to create the instance of the description
			void *mkr = dlsym(hndl, "maker");
			pixgeodescrptr = reinterpret_cast<EUTelGenericPixGeoDescr*(*)()>(mkr)();

			streamlog_out( MESSAGE3 ) << "Inserting " << planeID << " into map" << std::endl;
			_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
			_pixelDescriptions.insert( std::make_pair(geoName, pixgeodescrptr) );
		}
		
		catch(...)
		{
			streamlog_out( ERROR3 ) << "Could not retrieve pointer to EUTelGenericPixGeoDescr object, did you include a maker() function w/o C++ name mangling?" << std::endl;
			//While runtime errors are not processed in Marlin, they will cause a controlled crash!
			throw std::runtime_error("Retrieval of EUTelGenericPixGeoDescr* failed");
		}
	}

	streamlog_out( MESSAGE3 )  << "Adding plane: " << planeID << " with geoLibName: " << geoName << " in volume " << planeVolume << std::endl;
	
	//Call the factory method to actually load the geoemtry!
	pixgeodescrptr->createRootDescr(planeVolume);
}

EUTelGenericPixGeoDescr* EUTelGenericPixGeoMgr::getPixGeoDescr(int planeID)
{
	std::map<int, EUTelGenericPixGeoDescr*>::iterator returnGeoDescrIt = _geoDescriptions.find(planeID);
	if(returnGeoDescrIt == _geoDescriptions.end())
	{
		std::stringstream ss;
		ss << "Could not retrieve plane: " << planeID << " from PixGeoMgr. Did you include it in your GEAR file?";
		throw std::runtime_error( ss.str() );
	}
	else
	{
		return static_cast<EUTelGenericPixGeoDescr*>( returnGeoDescrIt->second );
	}
}
