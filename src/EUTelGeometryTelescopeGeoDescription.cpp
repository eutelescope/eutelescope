/* 
 * File:   EUTelGeometryTelescopeGeoDescription.cpp
 * 
 */
#include "EUTelGeometryTelescopeGeoDescription.h"

// C++
#include <algorithm>
#include <string>
#include <cstring>
#include <sstream>

// MARLIN
#include "marlin/Global.h"
#include "marlin/VerbosityLevels.h"

//GEAR
#include "GEAR.h" 
#include "gearxml/GearXML.h"

// EUTELESCOPE
#include "EUTelExceptions.h"
#include "EUTelGenericPixGeoMgr.h"
#include "EUTelNav.h"

// ROOT
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoBBox.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TError.h"

using namespace eutelescope;
using namespace geo;

unsigned EUTelGeometryTelescopeGeoDescription::_counter = 0;

/**TODO: Replace me: NOP*/
EUTelGeometryTelescopeGeoDescription& EUTelGeometryTelescopeGeoDescription::getInstance( gear::GearMgr* _g )
{
	static  EUTelGeometryTelescopeGeoDescription instance;
	unsigned i = EUTelGeometryTelescopeGeoDescription::_counter;
	
	//do it only once!
	if( i < 1 )
	{
		instance.setGearManager(_g);
		instance.readGear();
	}
	
	instance.counter();
	return instance;
}

/**TODO: Replace me: NOP*/
int EUTelGeometryTelescopeGeoDescription::sensorIDtoZOrder( int planeID ) const
{
	std::map<int,int>::const_iterator it = _sensorIDtoZOrderMap.find(planeID);
	if( it != _sensorIDtoZOrderMap.end() )
	{
		return it->second;
	}
	else
	{
		std::stringstream ss;
		ss << planeID;
		std::string errMsg = "EUTelGeometryTelescopeGeoDescription::sensorIDtoZOrder: Could not find planeID: " + ss.str(); 
		throw InvalidGeometryException(errMsg);
	}
}

//Note  that to determine these axis we MUST use the geometry class after initialisation. By this I mean directly from the root file create.
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneNormal( int planeID )
{
	std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
	if( it != _sensorIDVec.end() )
	{
		const double zAxisLocal[3]  = {0,0,1};
		double zAxisGlobal[3]; 
		local2MasterVec(planeID,zAxisLocal,zAxisGlobal); 
		TVector3 normVec(zAxisGlobal[0],zAxisGlobal[1],zAxisGlobal[2]);
		return normVec;
	}
	else
	{
		std::stringstream ss;
		ss << planeID;
		std::string errMsg = "EUTelGeometryTelescopeGeoDescription::siPlaneNormal: Could not find planeID: " + ss.str();
		throw InvalidGeometryException(errMsg);
	}
}

/**TODO: Replace me: NOP*/
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneXAxis( int planeID )
{
	std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
	if( it != _sensorIDVec.end() )
	{
		const double xAxisLocal[3]  = {1,0,0};
		double xAxisGlobal[3]; 
		local2MasterVec(planeID,xAxisLocal , xAxisGlobal ); 
		TVector3 xVec(xAxisGlobal[0],xAxisGlobal[1],xAxisGlobal[2]);
		return xVec;
	}
	else
	{
		std::stringstream ss;
		ss << planeID;
		std::string errMsg = "EUTelGeometryTelescopeGeoDescription::siPlaneXAxis: Could not find planeID: " + ss.str();
		throw InvalidGeometryException(errMsg);
	}
}


/**TODO: Replace me: NOP*/
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneYAxis( int planeID )
{
	std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
	if( it != _sensorIDVec.end() )
	{
		const double yAxisLocal[3]  = {0,1,0};
		double yAxisGlobal[3]; 
		local2MasterVec(planeID,yAxisLocal , yAxisGlobal ); 
		TVector3 yVec(yAxisGlobal[0],yAxisGlobal[1],yAxisGlobal[2]);
		return yVec;
	}
	else
	{
		std::stringstream ss;
		ss << planeID;
		std::string errMsg = "EUTelGeometryTelescopeGeoDescription::siPlaneYAxis: Could not find planeID: " + ss.str(); 
		throw InvalidGeometryException(errMsg);
	}
}

/**TODO: Replace me: NOP*/
void EUTelGeometryTelescopeGeoDescription::initialisePlanesToExcluded(IntVec planeIDs)
{
	int counter=0;
	for(size_t i = 0 ; i <_sensorZOrderToIDMap.size(); ++i){
			bool excluded=false;
			for(size_t j =0; j< planeIDs.size(); ++j){
					if(_sensorZOrderToIDMap[i] == planeIDs[j]){
							excluded=true;
							break;
					} 
			}
			if(!excluded){
					_sensorIDToZOrderWithoutExcludedPlanes[_sensorZOrderToIDMap[i]] =  counter;
					_sensorZOrderToIDWithoutExcludedPlanes[counter]=_sensorZOrderToIDMap[i];
					counter++;
			}
	}
	//Check if the number of excluded planes set is the same as (total-number of plane IDs inputed that should be excluded)
	if(geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() != (geo::gGeometry().sensorIDstoZOrder().size()-planeIDs.size())){
			throw(lcio::Exception( "The number of Planes-Excluded is not correct. This could be a problem with geometry."));
	}else{
			streamlog_out(DEBUG0) <<"The correct number of planes have been excluded" << std::endl;
	}
}

/** Sensor ID vector ordered according to their position along the Z axis (beam axis)
 *  Numeration runs from 0 to nPlanes-1 */
int EUTelGeometryTelescopeGeoDescription::sensorZOrderToID( int znumber ) const
{
	std::map<int,int>::const_iterator it = _sensorZOrderToIDMap.find( znumber );
	if( it != _sensorZOrderToIDMap.end() )
	{
		return it->second;
	}
	else
	{
		std::stringstream ss;
		ss << znumber;
		std::string errMsg = "EUTelGeometryTelescopeGeoDescription::sensorZOrderToID: Could not find snumber: " + ss.str(); 
		throw InvalidGeometryException(errMsg);
	}
}

void EUTelGeometryTelescopeGeoDescription::readSiPlanesLayout()
{
	// sensor-planes in geometry navigation:
	_siPlanesParameters = const_cast< gear::SiPlanesParameters*> (&( _gearManager->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast< gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));

	_nPlanes = _siPlanesLayerLayout->getNLayers(); 

	//read the geoemtry names from the "Geometry" StringVec section of the gear file
	lcio::StringVec geometryNameParameters;

	try
	{
		geometryNameParameters  =  _siPlanesParameters->getStringVals("Geometry");
	}
	catch(gear::UnknownParameterException e)
	{
		streamlog_out(MESSAGE6) << "No Geometry field found in GEAR file, assuming CAST for all planes" << std::endl;
		for(size_t i = 0; i < _nPlanes; i++)
		{
			geometryNameParameters.push_back("CAST");
		}
	}

	setSiPlanesLayoutID( _siPlanesParameters->getSiPlanesID() ) ;

	// create an array with the z positions of each layer
	for (size_t iPlane = 0; iPlane < _nPlanes; iPlane++)
	{
		EUTelPlane thisPlane;

		thisPlane.xPos	= _siPlanesLayerLayout->getLayerPositionX(iPlane);
		thisPlane.yPos	= _siPlanesLayerLayout->getLayerPositionY(iPlane);
		thisPlane.zPos	= _siPlanesLayerLayout->getLayerPositionZ(iPlane);

		thisPlane.alpha	= _siPlanesLayerLayout->getLayerRotationZY(iPlane);
		thisPlane.beta	= _siPlanesLayerLayout->getLayerRotationZX(iPlane);
		thisPlane.gamma	= _siPlanesLayerLayout->getLayerRotationXY(iPlane);

		thisPlane.pixGeoName = geometryNameParameters[iPlane];

		thisPlane.r1	= _siPlanesLayerLayout->getSensitiveRotation1(iPlane);
		thisPlane.r2	= _siPlanesLayerLayout->getSensitiveRotation2(iPlane);
		thisPlane.r3	= _siPlanesLayerLayout->getSensitiveRotation3(iPlane);
		thisPlane.r4	= _siPlanesLayerLayout->getSensitiveRotation4(iPlane);

		thisPlane.xSize	= _siPlanesLayerLayout->getSensitiveSizeX(iPlane);
		thisPlane.ySize	= _siPlanesLayerLayout->getSensitiveSizeY(iPlane);
		thisPlane.zSize	= _siPlanesLayerLayout->getSensitiveThickness(iPlane);

		thisPlane.xPixelNo	= _siPlanesLayerLayout->getSensitiveNpixelX(iPlane);
		thisPlane.yPixelNo	= _siPlanesLayerLayout->getSensitiveNpixelY(iPlane);
		thisPlane.xPitch	= _siPlanesLayerLayout->getSensitivePitchX(iPlane);
		thisPlane.yPitch	= _siPlanesLayerLayout->getSensitivePitchY(iPlane); 
		thisPlane.xRes	= _siPlanesLayerLayout->getSensitiveResolution(iPlane); //should be ResolutionX
		thisPlane.yRes	= _siPlanesLayerLayout->getSensitiveResolution(iPlane); //should be ResolutionY

		thisPlane.radLength	= _siPlanesLayerLayout->getSensitiveRadLength(iPlane);

		_planeSetup[_siPlanesLayerLayout->getID(iPlane)] = thisPlane;
	}


	// sort the array with increasing z
	std::sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());

	// clear the sensor ID vector
	_sensorIDVec.clear();

	// clear the sensor ID map
	_sensorIDVecMap.clear();
	_sensorIDtoZOrderMap.clear();

	for(int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++)
	{
		int sensorID = _siPlanesLayerLayout->getID(iPlane);
		_sensorIDVec.push_back(sensorID);
		_sensorIDVecMap.insert(std::make_pair(sensorID, iPlane));

		// count number of the sensors to the left of the current one:
		int sensorsToTheLeft = 0;
		int kposition = _siPlanesLayerLayout->getSensitivePositionZ(iPlane);
		for (int jPlane = 0; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++)
		{
			if(_siPlanesLayerLayout->getSensitivePositionZ(jPlane) + 1e-06 < kposition  )
			{
				sensorsToTheLeft++;
			}
		}
		_sensorZOrderToIDMap.insert(std::make_pair(sensorsToTheLeft, sensorID));        
		_sensorIDtoZOrderMap.insert(std::make_pair(sensorID, sensorsToTheLeft));
	}
	_nPlanes = _siPlanesParameters->getSiPlanesNumber();
}

void EUTelGeometryTelescopeGeoDescription::readTrackerPlanesLayout()
{
	// sensor-planes in geometry navigation:
	_trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&( _gearManager->getTrackerPlanesParameters()));
	_trackerPlanesLayerLayout = const_cast< gear::TrackerPlanesLayerLayout*> (&(_trackerPlanesParameters->getTrackerPlanesLayerLayout()));

	setSiPlanesLayoutID( _trackerPlanesParameters->getLayoutID() ) ;

	// clear the sensor ID vector
	_sensorIDVec.clear();
	// clear the sensor ID map
	_sensorIDVecMap.clear();
	_sensorIDtoZOrderMap.clear();
	
	//should be filled based on the length of the sensor vector after the loop
	_nPlanes = 0; 

	// create an array with the z positions of each layer
	int nLayers = _trackerPlanesLayerLayout->getNLayers();
	for (int iLayer = 0; iLayer < nLayers; iLayer++)
	{
		gear::TrackerPlanesLayerImpl* _trackerPlanesLayerImpl = const_cast<gear::TrackerPlanesLayerImpl*>(_trackerPlanesLayerLayout->getLayer( iLayer));

		int nsensitive = _trackerPlanesLayerImpl->getNSensitiveLayers();
		gear::TrackerPlanesSensitiveLayerImplVec& vector = _trackerPlanesLayerImpl->getSensitiveLayerVec();

		for(int iSensLayer = 0; iSensLayer < nsensitive; iSensLayer++)
		{       

			gear::TrackerPlanesSensitiveLayerImpl& sensitiveLayer = vector.at(iSensLayer);
			int sensorID =   sensitiveLayer.getID();

			EUTelPlane thisPlane;

			thisPlane.xPos	= sensitiveLayer.getPositionX();
			thisPlane.yPos	= sensitiveLayer.getPositionY();
			thisPlane.zPos	= sensitiveLayer.getPositionZ();

			thisPlane.alpha	= sensitiveLayer.getRotationZY();
			thisPlane.beta	= sensitiveLayer.getRotationZX();
			thisPlane.gamma	= sensitiveLayer.getRotationXY();

			thisPlane.pixGeoName = sensitiveLayer.getInfo();

			//TODO
			thisPlane.r1	= 1;//sensitiveLayer.getRotation1();
			thisPlane.r2	= 0;//sensitiveLayer.getRotation2();
			thisPlane.r3	= 0;//sensitiveLayer.getRotation3();
			thisPlane.r4	= 1;//sensitiveLayer.getRotation4();

			thisPlane.xSize	= sensitiveLayer.getSizeX();
			thisPlane.ySize	= sensitiveLayer.getSizeY();
			thisPlane.zSize	= sensitiveLayer.getThickness();

			thisPlane.xPixelNo	= sensitiveLayer.getNpixelX();
			thisPlane.yPixelNo	= sensitiveLayer.getNpixelY();
			thisPlane.xPitch	= sensitiveLayer.getPitchX();
			thisPlane.yPitch	= sensitiveLayer.getPitchY();
			thisPlane.xRes	= sensitiveLayer.getResolutionX();
			thisPlane.yRes	= sensitiveLayer.getResolutionY();

			thisPlane.radLength	= sensitiveLayer.getRadLength();

			_planeSetup[sensorID] = thisPlane;

			_sensorIDVec.push_back(sensorID);
			_sensorIDVecMap.insert(std::make_pair(sensorID, iLayer)); // what if there are more then 1 sensore per layer?
			streamlog_out(DEBUG1) << " iter: " << _sensorIDVec.at( _sensorIDVec.size()-1 ) << " " << sensorID << " " << sensitiveLayer.getInfo() .c_str() << std::endl; 
		}
	}

	_nPlanes =  _sensorIDVec.size(); 

	for(size_t i=0; i< _siPlaneZPosition.size(); i++)
	{
		int sensorsToTheLeft = 0;
		int sensorID = _sensorIDVec.at(i);

		for(size_t j=0; j< _siPlaneZPosition.size(); j++)
		{ 
			if( _siPlaneZPosition.at(j) < _siPlaneZPosition.at(i) - 1e-06 )
			{
				sensorsToTheLeft++;
			}
		}
		_sensorZOrderToIDMap.insert(std::make_pair(sensorsToTheLeft, sensorID));        
		_sensorIDtoZOrderMap.insert(std::make_pair(sensorID, sensorsToTheLeft));
	}
}

EUTelGeometryTelescopeGeoDescription::EUTelGeometryTelescopeGeoDescription() :
_gearManager( marlin::Global::GEAR ),
_siPlanesDefined(false),
_telPlanesDefined(false),
_siPlanesParameters(nullptr),
_siPlanesLayerLayout(nullptr),
_trackerPlanesParameters(nullptr),
_trackerPlanesLayerLayout(nullptr),
_sensorIDVec(),
_sensorIDVecMap(),
_sensorZOrderToIDMap(),
_sensorIDtoZOrderMap(),
_nPlanes(0),
_isGeoInitialized(false),
_geoManager(nullptr)
{
	//Set ROOTs verbosity to only display error messages or higher (so info will not be streamed to stderr)
	gErrorIgnoreLevel =  kError;  

	//Pixel Geometry manager creation
	_pixGeoMgr = new EUTelGenericPixGeoMgr();
}

void EUTelGeometryTelescopeGeoDescription::readGear()
{
    if( _gearManager == nullptr )
	{
        streamlog_out(ERROR2) << "The GearMgr is not available, for an unknown reason." << std::endl;
        throw eutelescope::InvalidGeometryException("GEAR manager is not initialised");
    }

    try
	{
      _siPlanesParameters = const_cast< gear::SiPlanesParameters*> (&(_gearManager->getSiPlanesParameters()));
      streamlog_out(MESSAGE1)  << "gear::SiPlanes : " << _siPlanesParameters << std::endl;
      _siPlanesDefined = true;
    }
	catch(...)
	{
		streamlog_out(WARNING)   << "gear::SiPlanes NOT found " << std::endl;
    }
    try
	{
		_trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&(_gearManager->getTrackerPlanesParameters()));
		streamlog_out(MESSAGE1)  << "gear::TrackerPlanes : " << _trackerPlanesParameters << std::endl;
		_telPlanesDefined = true;
	}
	catch(...)
	{
		streamlog_out(WARNING)   << "gear::TrackerPlanes NOT found "  << std::endl;
    }

    if( _siPlanesDefined )
	{
		readSiPlanesLayout();
    }
    else if( _telPlanesDefined )
	{
		readTrackerPlanesLayout();
    }
	else
	{
		streamlog_out(ERROR5) << "Your GEAR file neither contains SiPlanes nor TrackerPlanes and thus is not valid" << std::endl;
		throw eutelescope::InvalidGeometryException("GEAR file invalid, does not contain SiPlanes nor TrackerPlanes");
	}
}

EUTelGeometryTelescopeGeoDescription::~EUTelGeometryTelescopeGeoDescription()
{
	delete _geoManager;
	delete _pixGeoMgr;
}

/**
 * Initialise ROOT geometry objects from external .root file
 * @param tgeofilename name of .root file
 */
void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( std::string tgeofilename ) {
    
    _geoManager = TGeoManager::Import( tgeofilename.c_str() );
    if( !_geoManager ) {
        streamlog_out( WARNING ) << "Can't read file " << tgeofilename << std::endl;
    }

    _geoManager->CloseGeometry();
}

/**
 *
 */
void EUTelGeometryTelescopeGeoDescription::translateSiPlane2TGeo(TGeoVolume* pvolumeWorld, int SensorId ){
	double xc, yc, zc;   // volume center position 
	double alpha, beta, gamma;
	double rotRef1, rotRef2, rotRef3, rotRef4; // for backward compatibility with previous GEAR. We only need 2 entries from gear file for fast z rotation function in TGeoRotations. 

	std::stringstream strId;
	strId << SensorId;

	// Get sensor center position
	xc = siPlaneXPosition( SensorId );
	yc = siPlaneYPosition( SensorId );
	zc = siPlaneZPosition( SensorId );

	// Get sensor orientation
	alpha = siPlaneXRotation( SensorId ); //  in degrees !
	beta  = siPlaneYRotation( SensorId ); // 
	gamma = siPlaneZRotation( SensorId ); // 

	rotRef1 = siPlaneRotation1( SensorId );
	rotRef2 = siPlaneRotation2( SensorId );
	rotRef3 = siPlaneRotation3( SensorId );
	rotRef4 = siPlaneRotation4( SensorId );

	//We must check that the input is correct. Since this is a combination of initial rotations and reflections the determinate must be 1 or -1
	float	determinant = rotRef1*rotRef4 - rotRef2*rotRef3  ;
	if(determinant==1 or determinant==-1){ 
		streamlog_out(MESSAGE5) << "SensorID: " << SensorId << ". Determinant =  " <<determinant <<"  This is the correct determinate for this transformation." << std::endl;   
	}else{
		streamlog_out(ERROR5) << "SensorID: " << SensorId << ". Determinant =  " <<determinant << std::endl;   
		throw(lcio::Exception("The initial rotation and reflection matrix does not have determinant of 1 or -1. Gear file input must be wrong.")); 	
	}
	//Create spatial TGeoTranslation object.
	std::string stTranslationName = "matrixTranslationSensor";
	stTranslationName.append( strId.str() );
	TGeoTranslation* pMatrixTrans = new TGeoTranslation( stTranslationName.c_str(), xc, yc, zc );
	//ALL clsses deriving from TGeoMatrix are not owned by the ROOT geometry manager, invoking RegisterYourself() transfers
	//ownership and thus ROOT will clean up
	pMatrixTrans->RegisterYourself();      

	//Create TGeoRotation object. 
	//Translations are of course just positional changes in the global frame.
	//Note that each subsequent rotation is using the new coordinate system of the last transformation all the way back to the global frame.
	//The way to think about this is that each rotation is the multiplication of the last rotation matrix by a new one.
	//The order is:
	//Integer Z rotation and reflections.
	//Z rotations specified by in degrees.
	//X rotations 
	//Y rotations
	TGeoRotation * pMatrixRotRefCombined = new TGeoRotation();
	double integerRotationsAndReflections[9]={rotRef1,rotRef2,0,rotRef3,rotRef4,0,0,0,1};
	pMatrixRotRefCombined->SetMatrix(integerRotationsAndReflections);
	pMatrixRotRefCombined->RotateZ(gamma);//Z Rotation (degrees)//This will again rotate a vector around z axis usign the right hand rule.  
	pMatrixRotRefCombined->RotateX(alpha);//X Rotations (degrees)//This will rotate a vector usign the right hand rule round the x-axis
	pMatrixRotRefCombined->RotateY(beta);//Y Rotations (degrees)//Same again for Y axis 
	pMatrixRotRefCombined->RegisterYourself();//We must allow the matrix to be used by the TGeo manager.
	// Combined translation and orientation
	TGeoCombiTrans* combi = new TGeoCombiTrans( *pMatrixTrans, *pMatrixRotRefCombined );
	//This is to print to screen the rotation and translation matrices used to transform from local to global frame.
	streamlog_out(MESSAGE9) << "THESE MATRICES ARE USED TO TAKE A POINT IN THE LOCAL FRAME AND MOVE IT TO THE GLOBAL FRAME."  << std::endl;   
	streamlog_out(MESSAGE9) << "SensorID: " << SensorId << " Rotation/Reflection matrix for this object."  << std::endl;   
	const double* rotationMatrix =  combi->GetRotationMatrix();	
	streamlog_out(MESSAGE9) << std::setw(10) <<rotationMatrix[0]<<"  "<<rotationMatrix[1]<<"   "<<rotationMatrix[2]<< std::endl;
	streamlog_out(MESSAGE9) << std::setw(10) <<rotationMatrix[3]<<"  "<<rotationMatrix[4]<<"   "<<rotationMatrix[5]<< std::endl;
	streamlog_out(MESSAGE9) << std::setw(10) <<rotationMatrix[6]<<"  "<<rotationMatrix[7]<<"   "<<rotationMatrix[8]<< std::endl;

	//streamlog_out(MESSAGE9) << std::setw(10) <<rotationMatrix[0] << std::setw(10) <<rotationMatrix[1]<< std::setw(10) <<rotationMatrix[2]<< std::setw(10)<< std::endl<< std::endl; 
	//streamlog_out(MESSAGE9) << std::setw(10) <<rotationMatrix[3] << std::setw(10) <<rotationMatrix[4]<< std::setw(10) <<rotationMatrix[5]<< std::setw(10)<< std::endl<< std::endl; 
	//streamlog_out(MESSAGE9) << std::setw(10) <<rotationMatrix[6] << std::setw(10) <<rotationMatrix[7]<< std::setw(10) <<rotationMatrix[8]<< std::setw(10)<< std::endl<< std::endl; 
	const double* translationMatrix =  combi->GetTranslation();	
	streamlog_out(MESSAGE9) << "SensorID: " << SensorId << " Translation vector for this object."  << std::endl;   
	streamlog_out(MESSAGE9) << std::setw(10) <<translationMatrix[0] << std::setw(10) <<translationMatrix[1]<< std::setw(10) <<translationMatrix[2]<< std::setw(10)<< std::endl; 

	combi->RegisterYourself();   
	

	// Construction of sensor objects

	// Construct object medium. Required for radiation length determination

	// assume SILICON, though all information except of radiation length is ignored
	double a       = 28.085500;     
	double z       = 14.000000;
	double density = 2.330000;
	double radl    = siPlaneRadLength( SensorId );
	double absl    = 45.753206;
	std::string stMatName = "materialSensor";
	stMatName.append( strId.str() );
	TGeoMaterial* pMat = new TGeoMaterial( stMatName.c_str(), a, z, density, radl, absl );
	pMat->SetIndex( 1 );
	// Medium: medium_Sensor_SILICON
	int numed   = 0;  // medium number
	double par[8];
	par[0]  = 0.000000; // isvol
	par[1]  = 0.000000; // ifield
	par[2]  = 0.000000; // fieldm
	par[3]  = 0.000000; // tmaxfd
	par[4]  = 0.000000; // stemax
	par[5]  = 0.000000; // deemax
	par[6]  = 0.000000; // epsil
	par[7]  = 0.000000; // stmin
	std::string stMedName = "mediumSensor";
	stMedName.append( strId.str() );
	TGeoMedium* pMed = new TGeoMedium( stMedName.c_str(), numed, pMat, par );

	// Construct object shape
	// Shape: Box type: TGeoBBox
	// TGeo requires half-width of box side
	Double_t dx = siPlaneXSize( SensorId ) / 2.;
	Double_t dy = siPlaneYSize( SensorId ) / 2.;
	Double_t dz = siPlaneZSize( SensorId ) / 2.;
	TGeoShape *pBoxSensor = new TGeoBBox( "BoxSensor", dx, dy, dz );
	// Volume: volume_Sensor1

	// Geometry navigation package requires following names for objects that have an ID
	// name:ID
	std::string stVolName = "volume_SensorID:";
	stVolName.append( strId.str() );

	_planePath.insert( std::make_pair(SensorId, "/volume_World_1/"+stVolName+"_1") );

	TGeoVolume* pvolumeSensor = new TGeoVolume( stVolName.c_str(), pBoxSensor, pMed );
	pvolumeSensor->SetVisLeaves( kTRUE );
	pvolumeWorld->AddNode(pvolumeSensor, 1/*(SensorId)*/, combi);

	//this line tells the pixel geometry manager to load the pixel geometry into the plane			
	streamlog_out(DEBUG1) << " sensorID: " << SensorId << " " << stVolName << std::endl;   
	std::string name = geoLibName(SensorId);

	if( name == "CAST" )
	{
	_pixGeoMgr->addCastedPlane( SensorId, siPlaneXNpixels(SensorId), siPlaneYNpixels(SensorId), siPlaneXSize(SensorId), siPlaneYSize(SensorId), siPlaneZSize(SensorId), siPlaneRadLength(SensorId), stVolName);
	}

	else
	{
	_pixGeoMgr->addPlane( SensorId, name, stVolName);
	}
}

/**
 * Initialise ROOT geometry objects from GEAR objects
 * 
 * @param geomName name of ROOT geometry object
 * @param dumpRoot dump automatically generated ROOT geometry file for further inspection
 */
void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( std::string& geomName, bool dumpRoot = false ) {
    // get access to ROOT's geometry manager
    
	if( _isGeoInitialized )
	{
		streamlog_out( WARNING3 ) << "EUTelGeometryTelescopeGeoDescription: Geometry already initialized, using old initialization" << std::endl;
		return;
	}
	else
	{
    		_geoManager = new TGeoManager("Telescope", "v0.1");
	}

	if( !_geoManager )
	{
		streamlog_out( ERROR3 ) << "Can't instantiate ROOT TGeoManager " << std::endl;
		return;
	}
   
    
    // Create top world volume containing telescope geometry
    
    // Create air mixture
    // see http://pdg.lbl.gov/2013/AtomicNuclearProperties/HTML_PAGES/104.html
    double air_density = 1.2e-3;         // g/cm^3
    double air_radlen  = 36.62;          // g/cm^2
    TGeoMixture* pMatAir = new TGeoMixture("AIR",3,air_density);
    pMatAir->DefineElement(0, 14.007, 7.,  0.755267 );     //Nitrogen
    pMatAir->DefineElement(1, 15.999, 8.,  0.231781 );     //Oxygen
    pMatAir->DefineElement(2, 39.948, 18., 0.012827 );     //Argon
    pMatAir->DefineElement(3, 12.011, 6.,  0.000124 );     //Carbon
    pMatAir->SetRadLen( air_radlen );
    // Medium: medium_World_AIR
    TGeoMedium* pMedAir = new TGeoMedium("medium_World_AIR", 3, pMatAir );

    // The World is the 10 x 10m x 10m box filled with air mixture
    Double_t dx,dy,dz;
    dx = 5000.000000; // [mm]
    dy = 5000.000000; // [mm]
    dz = 5000.000000; // [mm]
    TGeoShape *pBoxWorld = new TGeoBBox("Box_World", dx,dy,dz);
    // Volume: volume_World
    TGeoVolume* pvolumeWorld = new TGeoVolume("volume_World",pBoxWorld, pMedAir);
    pvolumeWorld->SetLineColor(4);
    pvolumeWorld->SetLineWidth(3);
    pvolumeWorld->SetVisLeaves(kTRUE);

   // Set top volume of geometry
   gGeoManager->SetTopVolume( pvolumeWorld );
   
   IntVec::const_iterator itrPlaneId;
   for ( itrPlaneId = _sensorIDVec.begin(); itrPlaneId != _sensorIDVec.end(); ++itrPlaneId ) {
       
       translateSiPlane2TGeo(pvolumeWorld, *itrPlaneId );
 
   } // loop over sensorID

    _geoManager->CloseGeometry();
    _isGeoInitialized = true;
    // Dump ROOT TGeo object into file
    if ( dumpRoot ) _geoManager->Export( geomName.c_str() );
    return;
}

Eigen::Matrix3d EUTelGeometryTelescopeGeoDescription::rotationMatrixFromAngles(int sensorID)
{
	return rotationMatrixFromAngles( (long double)siPlaneXRotationRadians(sensorID), (long double)siPlaneYRotationRadians(sensorID), (long double)siPlaneZRotationRadians(sensorID) );
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::globalXAxis(int sensorID)
{
	Eigen::Vector3d xAxis(1,0,0);
	return rotationMatrixFromAngles(sensorID)*xAxis;
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::globalYAxis(int sensorID)
{
	Eigen::Vector3d yAxis(0,1,0);
	return rotationMatrixFromAngles(sensorID)*yAxis;
}

/** Returns the rotation matrix for given angles
 *  It alpha rotation is around initial X axis, then beta around the new Y' axis
 *  and finally the gamam rotation around the new Z'' axis */
Eigen::Matrix3d EUTelGeometryTelescopeGeoDescription::rotationMatrixFromAngles(long double alpha, long double beta, long double gamma)
{
	//Eigen::IOFormat IO(6, 0, ", ", ";\n", "[", "]", "[", "]");
	//std::cout << "alpha, beta, gamma: " << alpha << ", " << beta << ", " << gamma << std::endl;
	long double cosA = cos(alpha);
	long double sinA = sin(alpha);
	long double cosB = cos(beta);
	long double sinB = sin(beta);
	long double cosG = cos(gamma);
	long double sinG = sin(gamma);

	//std::cout << "trig" << cosA << ", " << cosB << ", " << cosG << ", " << sinA << ", " << sinB << ", " << sinG <<  std::endl;

	Eigen::Matrix3d rotMat;
	rotMat <<	(double)(cosB*cosG+sinA*sinB*sinG),	(double)(sinA*sinB*cosG-cosB*sinG),	(double)(cosA*sinB),
	      		(double)(cosA*sinG),			(double)(cosA*cosG),			(double)(-sinA),
			(double)(sinA*cosB*sinG-sinB*cosG),	(double)(sinA*cosB*cosG+sinB*sinG),	(double)(cosA*cosB);
	//std::cout << rotMat.format(IO) << std::endl;
	return rotMat;
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getRotationAnglesFromMatrix(Eigen::Matrix3d rotMat)
{
	long double alpha = asin((long double)(-rotMat(1,2)));
	long double cosA = cos(alpha);

	long double beta = asin((long double)(rotMat(0,2)/cosA));
	long double gamma = asin((long double)(rotMat(1,0)/cosA));

	Eigen::Vector3d vec;

	vec << alpha, beta, gamma;
	return vec;
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getOffsetVector(int sensorID)
{
	Eigen::Vector3d offsetVec;
	offsetVec << siPlaneXPosition(sensorID), siPlaneYPosition(sensorID), siPlaneZPosition(sensorID); 
	return offsetVec;
}

Eigen::Matrix3d EUTelGeometryTelescopeGeoDescription::getFlipMatrix(int sensorID)
{
	Eigen::Matrix3d flipMat;
	flipMat << 	siPlaneRotation1(sensorID),	siPlaneRotation2(sensorID),	0,
			siPlaneRotation3(sensorID), 	siPlaneRotation4(sensorID),	0,
			0,				0,				1;	
	return flipMat;
}

/** Determine id of the sensor in which point is locate
 * 
 * @param globalPos 3D point in global reference frame
 * @return sensorID or -999 if the point in outside of sensor volume
 */
//MUST OUTPUT -999 TO SIGNIFY THAT NO SENSOR HAS BEEN FOUND. SINCE USED IN PATTERN RECOGNITION THIS WAY.
int EUTelGeometryTelescopeGeoDescription::getSensorID( const float globalPos[] ) const {
    streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::getSensorID() " << std::endl;
    
    _geoManager->FindNode( globalPos[0], globalPos[1], globalPos[2] );

    std::vector<std::string> split;
 
    int sensorID = -999;

    const char* volName1 = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
    streamlog_out(DEBUG2) << "init sensorID  : " << sensorID  <<  " " << volName1 << std::endl;

    while( _geoManager->GetLevel() > 0 ) { 
      const char* volName = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
      streamlog_out( DEBUG1 ) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") found in volume: " << volName << " level: " << _geoManager->GetLevel() << std::endl;
      split = Utility::stringSplit( std::string( volName ), "/", false);
      if ( split.size() > 0 && split[0].length() > 16 && (split[0].substr(0,16) == "volume_SensorID:") ) {
         int strLength = split[0].length(); 
         sensorID = strtol( (split[0].substr(16, strLength )).c_str(), NULL, 10 );
         streamlog_out(DEBUG1) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was found at :" << sensorID << std::endl;       
         break;
      }
      _geoManager->CdUp();	////////////////////////////////////////THIS NEEDS TO BE FIXED. If partice falls in the pixel volume and to find sensor ID you need to be on the sensor volume
    }

    const char* volName2 = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
    streamlog_out( DEBUG2 ) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") found in volume: " << volName2 << " no moving around any more" << std::endl;
    
          if( sensorID >= 0 )
	  {
//                sensorID = strtol( (split[0].substr(16, strLength )).c_str(), NULL, 10 );
		streamlog_out(DEBUG5) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was found. sensorID = " << sensorID << std::endl;
          }
	  else
	  {
		streamlog_out(DEBUG5) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was not found inside any sensor! sensorID = " << sensorID << std::endl;
	  }
//        }
//        streamlog_out(DEBUG5) << "sensorID  : " << sensorID  << std::endl;

	return sensorID;
}

/**
 * Coordinate transformation from local reference frame of sensor with a given sensorID
 * to the global coordinate system
 * 
 * @param sensorID Id of the sensor (specifies local coordinate system)
 * @param localPos (x,y,z) in local coordinate system
 * @param globalPos (x,y,z) in global coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::local2Master( int sensorID, const double localPos[], double globalPos[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::local2Master() " << std::endl;
    const double sensorCenterX = siPlaneXPosition( sensorID );
    const double sensorCenterY = siPlaneYPosition( sensorID );
    const double sensorCenterZ = siPlaneZPosition( sensorID );
   
    streamlog_out(DEBUG0) << "Sensor id: " << sensorID << std::endl;
    streamlog_out(DEBUG0) << "Sensor center: " << "(" << sensorCenterX << "," << sensorCenterY << "," << sensorCenterZ << ")" << std::endl;

    _geoManager->cd( _planePath[sensorID].c_str() );
    _geoManager->GetCurrentNode()->LocalToMaster( localPos, globalPos );
 
    const char* volName = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
    streamlog_out( DEBUG0 ) << "sensorCenter (" << sensorCenterX << "," << sensorCenterY << "," << sensorCenterZ << ") found in volume: " << volName << std::endl;

    streamlog_out(DEBUG0) << std::fixed;
    streamlog_out(DEBUG0) << "Local coordinates ( sensorID =  " << sensorID << " ) : " << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << localPos[0] << std::setw(10) << std::setprecision(5) << localPos[1] << std::setw(10) << std::setprecision(5) << localPos[2] << std::endl;
    streamlog_out(DEBUG0) << "Global coordinates ( sensorID =  " << sensorID << " ) : " << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << globalPos[0] << std::setw(10) << std::setprecision(5) << globalPos[1] << std::setw(10) << std::setprecision(5) << globalPos[2] << std::endl;
}

/**
 * Coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param sensorID Id of the sensor (specifies local coordinate system)
 * @param globalPos (x,y,z) in global coordinate system
 * @param localPos (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::master2Local(int sensorID, const double globalPos[], double localPos[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::master2Local()--BEGIN " << std::endl;

    _geoManager->cd( _planePath[sensorID].c_str() );
		_geoManager->GetCurrentNode()->MasterToLocal( globalPos, localPos );
    
    streamlog_out(DEBUG0) << std::fixed;
    streamlog_out(DEBUG0) << "Global coordinates:" << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << globalPos[0] << std::setw(10) << std::setprecision(5) << globalPos[1] << std::setw(10) << std::setprecision(5) << globalPos[2] << std::endl;
    streamlog_out(DEBUG0) << "Local coordinates: " << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << localPos[0] << std::setw(10) << std::setprecision(5) << localPos[1] << std::setw(10) << std::setprecision(5) << localPos[2] << std::endl;
		streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::master2Local()----END " << std::endl;

}

/**
 * Vector coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param globalVec (x,y,z) in global coordinate system
 * @param localVec (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::local2MasterVec( int sensorID, const double localVec[], double globalVec[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::local2masterVec() " << std::endl;
    const double sensorCenterX = siPlaneXPosition( sensorID );
    const double sensorCenterY = siPlaneYPosition( sensorID );
    const double sensorCenterZ = siPlaneZPosition( sensorID );
    
    streamlog_out(DEBUG0) << "Sensor id: " << sensorID << std::endl;
    streamlog_out(DEBUG0) << "Sensor center: " << "(" << sensorCenterX << "," << sensorCenterY << "," << sensorCenterZ << ")" << std::endl;
    
    _geoManager->cd( _planePath[sensorID].c_str() );
    _geoManager->GetCurrentNode()->LocalToMasterVect( localVec, globalVec );

    streamlog_out(DEBUG0) << std::fixed;
    streamlog_out(DEBUG0) << "Global coordinates:" << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << globalVec[0] << std::setw(10) << std::setprecision(5) << globalVec[1] << std::setw(10) << std::setprecision(5) << globalVec[2] << std::endl;
    streamlog_out(DEBUG0) << "Local coordinates: " << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << localVec[0] << std::setw(10) << std::setprecision(5) << localVec[1] << std::setw(10) << std::setprecision(5) << localVec[2] << std::endl;
}


/**
 * Vector coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param globalVec (x,y,z) in global coordinate system
 * @param localVec (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::master2LocalVec( int sensorID, const double globalVec[], double localVec[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::master2LocalVec() " << std::endl;
    const double sensorCenterX = siPlaneXPosition( sensorID );
    const double sensorCenterY = siPlaneYPosition( sensorID );
    const double sensorCenterZ = siPlaneZPosition( sensorID );
    
    streamlog_out(DEBUG0) << "Sensor id: " << sensorID << std::endl;
    streamlog_out(DEBUG0) << "Sensor center: " << "(" << sensorCenterX << "," << sensorCenterY << "," << sensorCenterZ << ")" << std::endl;
    
    _geoManager->cd( _planePath[sensorID].c_str() );
    _geoManager->GetCurrentNode()->MasterToLocalVect( globalVec, localVec );

//    _geoManager->FindNode( sensorCenterX, sensorCenterY, sensorCenterZ );    
//    _geoManager->MasterToLocalVect( globalVec, localVec );
    
    streamlog_out(DEBUG0) << std::fixed;
    streamlog_out(DEBUG0) << "Global coordinates:" << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << globalVec[0] << std::setw(10) << std::setprecision(5) << globalVec[1] << std::setw(10) << std::setprecision(5) << globalVec[2] << std::endl;
    streamlog_out(DEBUG0) << "Local coordinates: " << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << localVec[0] << std::setw(10) << std::setprecision(5) << localVec[1] << std::setw(10) << std::setprecision(5) << localVec[2] << std::endl;
}

/**
 * Local-to-Global coordinate transformation matrix.
 * Corresponding volume is determined automatically.
 * 
 * @param globalPos (x,y,z) in global coordinate system
 * @return 
 */
const TGeoHMatrix* EUTelGeometryTelescopeGeoDescription::getHMatrix( const double globalPos[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getHMatrix()--------BEGIN " << std::endl;
    _geoManager->FindNode( globalPos[0], globalPos[1], globalPos[2] );    
    const TGeoHMatrix* globalH = _geoManager->GetCurrentMatrix();
	//	if(streamlog_out(DEBUG2)){
  //  streamlog_out(DEBUG2) << "Transformation matrix " << std::endl;
	//	globalH->Print();
	//	}
		
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getHMatrix()----END " << std::endl;
    return globalH;
}

/**
 * Calculate effective radiation length traversed by particle traveling between two points
 * along straight line.
 * 
 * Calculation is done according to the eq. (27.23)
 * @see http://pdg.lbl.gov/2006/reviews/passagerpp.pdf
 * 
 * @param globalPosStart starting point in the global coordinate system
 * @param globalPosFinish ending point in the global coordinate system
 * @param skipBoundaryVolumes if true subtract rad length of the volumes containing start and finish points
 * 
 * @return radiation length in units of X0
 */
float EUTelGeometryTelescopeGeoDescription::findRadLengthIntegral( const double globalPosStart[], const double globalPosFinish[], bool skipBoundaryPonitsVolumes ) {

    streamlog_out(DEBUG1) << "EUTelGeometryTelescopeGeoDescription::findRadLengthIntegral()" << std::endl;
    
    float rad = 0.;        // integral of radiation length in units of X0. Set to 0 at the start.
    
    const double mm2cm = 0.1;
    
    /* TGeo uses cm and grams as internal units e.g. in radiation length and density. Telescope/LCIO uses mm. Therefore this routine is full of 
     annoying conversion factors */    
    
		//Calcuate the distance^2 of start and final position
    const double stepLenght2 = ( globalPosFinish[0] - globalPosStart[0] )*( globalPosFinish[0] - globalPosStart[0] ) +
                               ( globalPosFinish[1] - globalPosStart[1] )*( globalPosFinish[1] - globalPosStart[1] ) +
                               ( globalPosFinish[2] - globalPosStart[2] )*( globalPosFinish[2] - globalPosStart[2] );
    
		//This is the distance between the start and final position
    const double stepLenght  = TMath::Sqrt( stepLenght2 );

    // don't need conversion factor to for calculation of directions. This is just the direction of the track.  WE ASSUME A LINEAR MOTION! This should be ok for a approximate solution
    const double xp  = ( globalPosFinish[0] - globalPosStart[0] )/stepLenght;
    const double yp  = ( globalPosFinish[1] - globalPosStart[1] )/stepLenght;
    const double zp  = ( globalPosFinish[2] - globalPosStart[2] )/stepLenght;

    streamlog_out(DEBUG0) << "Start point (x,y,z):" << globalPosStart[0] << "," << globalPosStart[1] << "," << globalPosStart[2] << std::endl;
    streamlog_out(DEBUG0) << "Finish point (x,y,z):" << globalPosFinish[0] << "," << globalPosFinish[1] << "," << globalPosFinish[2] << std::endl;
    streamlog_out(DEBUG0) << "Direction (nx,ny,nz):" << xp << "," << yp << "," << zp << std::endl;
    
    double snext;
    double pt[3], loc[3];
    double epsil = 1.E-2; //1.E-7; //This number must be in the order of 10s of microns to make sure you leave the starting volume everytime.
    double lastrad = 0.;
    int ismall       = 0;
    int nbound       = 0;
    float length     = 0.;
    TGeoMedium *med;
    TGeoShape *shape;
    
    // Get starting node
    gGeoManager->InitTrack( globalPosStart[0]/*mm*/, globalPosStart[1]/*mm*/, globalPosStart[2]/*mm*/, xp, yp, zp ); //This the start point and direction
    TGeoNode *nextnode = gGeoManager->GetCurrentNode( );
    
    double currentStep = stepLenght /*mm*/;  //This is total distance we are going to travel
    // Loop over all, encountered during the propagation, volumes 
    while ( nextnode ) {
        med = NULL;
        
	// Check if current point is inside silicon sensor. Radiation length of silicon sensors is accounted in thin scatterers of GBL.
        bool isBoundaryVolume = false;
        if ( gGeoManager->IsSameLocation( globalPosStart[0], globalPosStart[1], globalPosStart[2] ) ||
             gGeoManager->IsSameLocation( globalPosFinish[0], globalPosFinish[1], globalPosFinish[2] ) ) isBoundaryVolume = true; //This used to ignore the first and last boundary. Why I get the last bu why ignore the first?
        
        if ( nextnode ) med = nextnode->GetVolume()->GetMedium();  //We begin in the current node. Which includes the 'worldvolume' of air. If for some reason you are not in the world at all then return 0
        else return 0.;
        
        shape = nextnode->GetVolume()->GetShape();
        
        // make a step to the next intersection point
        if ( currentStep > 1.e-9 /*mm*/ ) nextnode = gGeoManager->FindNextBoundaryAndStep( currentStep /*mm*/ ); //Now move to the next volume if the step is above 0. If not then return the already calculated rad to that point 
        else return rad;
        
        snext  = gGeoManager->GetStep() /*mm*/; //This will output the distance traveled by FindNextBoundaryAndStep
        
        // Small steps treatment
        if ( snext < 1.e-8 /*mm*/ ) {
            ismall++;
            
            // Terminate calculation if too many small steps done
            if ( ismall > 3 ) {
                streamlog_out( WARNING1 ) << "ERROR: Small steps in: " << gGeoManager->GetPath() << " shape=" << shape->ClassName() << std::endl;
                return rad;
            }

            // increase step size (epsilon) and advance along the particle direction
            memcpy( pt, gGeoManager->GetCurrentPoint(), 3 * sizeof (double) );
            const double *dir = gGeoManager->GetCurrentDirection();
            for ( Int_t i = 0; i < 3; i++ ) pt[i] += epsil * dir[i]; //Move the current point slightly in the direction of motion. 
            snext = epsil; //Now set the change in distance to the epsil. This is done since snext is used to calculate the total raditaion length. NOTE epsil is just a small number since we assume we are within a substance
            length += snext; //Change the length traveled not by epsil.
            
            // Ignore start and finish volumes if required
            if ( skipBoundaryPonitsVolumes && isBoundaryVolume ) {//Now if it is not the last of first boundary
                rad += 0.;
            } else {
                rad += lastrad*snext; //This is the calculated (rad per distance x distance)
            }
            
            gGeoManager->CdTop( ); //This moves the current node to the top one. Which should be the world_volume
            nextnode = gGeoManager->FindNode( pt[0], pt[1], pt[2] );    // Check if particle is crossed the boundary with the new incrementally moved position
            if ( gGeoManager->IsOutside() ) return rad;                 // This checks if the particle is still with the geometry. If not then just return the current calculated radiation length
            TGeoMatrix *mat = gGeoManager->GetCurrentMatrix(); //This matrix is the transform from global to local coordinates              
            mat->MasterToLocal( pt, loc );                     //Now transform the global coordinates of pt to local coordinates loc.
            if ( !gGeoManager->GetCurrentVolume()->Contains( loc ) ) { //If then point is not in the volume then we must not be at the top so move up again and then set that to the current node
                gGeoManager->CdUp();
                nextnode = gGeoManager->GetCurrentNode();               // move to new volume
            }
            continue; //We continue since we dont need to calculate radiation length again since it is the same as before
        } else {
            ismall = 0;
        }//END OF SMALL STEP TREATMENT
        
        // Normal steps case
        nbound++;
        length += snext;
        currentStep -= snext;
        if ( med ) { //If medium is not NULL
            double radlen = med->GetMaterial()->GetRadLen() /*cm*/;
            if ( radlen > 1.e-9 && radlen < 1.e10 ) {
                
                lastrad = 1. / radlen * mm2cm; //calculate 1/radiationlength per cm
                
                // Ignore start and finish volumes if required
                if ( skipBoundaryPonitsVolumes && isBoundaryVolume ) { //Do the same which is done in small volume approximation. Add the radiation length if it is not the first or last volume
                    rad += 0.;
                } else {
                    rad += lastrad*snext; 
                }
                
            } else {
                lastrad = 0.;
            }
            streamlog_out( DEBUG0 ) << "STEP #" << nbound << std::endl;
            streamlog_out( DEBUG0 ) << "   step[mm]=" << snext << "   length[mm]=" << length
                    << " rad[X0]=" << snext * mm2cm / radlen << " " << med->GetName( ) 
                    << " rho[g/cm^3]=" << med->GetMaterial()->GetDensity() <<" radlen[cm]=" << radlen << " Boundary:" << (isBoundaryVolume?"yes":"no")
		    << std::endl;
        }//END OF IF MEDIUM NOT NULL
    }
    
    streamlog_out(DEBUG1) << "--------EUTelGeometryTelescopeGeoDescription::findRadLengthIntegral()--------" << std::endl;
    
    return rad;
}

//
// straight line - shashlyk plane assembler
//
int EUTelGeometryTelescopeGeoDescription::findNextPlane(  double* lpoint,  double* ldir, float* newpoint ){
	if(newpoint==NULL)
	{
		throw(lcio::Exception("You have passed a NULL pointer to findNextPlane.")); 	
	}
	//Here we set the normalised direction and starting point.
	double normdir = TMath::Sqrt(ldir[0]*ldir[0]+ldir[1]*ldir[1]+ldir[2]*ldir[2]); 
	streamlog_out( DEBUG0 ) << "::findNextPlane lpoint: "  << lpoint[0] << " " << lpoint[1] << " "<< lpoint[2] << " " << std::endl;
	ldir[0] = ldir[0]/normdir; 
	ldir[1] = ldir[1]/normdir; 
	ldir[2] = ldir[2]/normdir;
	streamlog_out( DEBUG0 ) << "::findNextPlane ldir  : "  << ldir  [0] << " " << ldir  [1] << " "<< ldir  [2] << " " << std::endl;

	for(int ip=0;ip<3;ip++) 
	{
	 newpoint[ip] = static_cast<float> (lpoint[ip]);
	}  
	int currentSensorID = getSensorID(newpoint); 
	//initialise the track.
	gGeoManager->InitTrack( lpoint, ldir );
	TGeoNode *node = gGeoManager->GetCurrentNode( );

	Int_t inode    = node->GetIndex();
	Int_t i        = 0;

	streamlog_out( DEBUG0 ) << "::findNextPlane look for next node, starting at node: " << node << " id: " << inode  << " currentSensorID: " << currentSensorID << std::endl;

	//   double kStep = 1e-03;
	while(( node = gGeoManager->FindNextBoundaryAndStep() ))
	{
		 inode = node->GetIndex();
		 streamlog_out( DEBUG0 ) << "::findNextPlane found next node: " << node << " id: " << inode << std::endl;
		 const double* point = gGeoManager->GetCurrentPoint();
		 const double* dir   = gGeoManager->GetCurrentDirection();
		 double ipoint[3] ;
		 double idir[3]   ;

		 for(int ip=0;ip<3;ip++) 
		 {
			 ipoint[ip] = point[ip];
			 idir[ip]   = dir[ip];
			 if(ip==2) ipoint[ip]+=0.01 ; // assumption !!! step by one um into the new volume // new volume is thicker than 1 um
			 newpoint[ip] = static_cast<float> (ipoint[ip]);
		 }  
		 int sensorID = getSensorID(newpoint); 
		 i++;     
		
		 gGeoManager->SetCurrentPoint( ipoint);
		 gGeoManager->SetCurrentDirection( idir);

		 streamlog_out( DEBUG0 ) << "::findNextPlane i=" << i  << " " << inode << " " << ipoint[0]  << " " << ipoint[1] << " " << ipoint[2]  << " sensorID:" << sensorID <<  std::endl;
		 if(sensorID >= 0 && sensorID != currentSensorID ) return sensorID;
	}


	return -100;
}
//This will take in a global coordinate and direction and output the new global point on the next sensor. 
bool EUTelGeometryTelescopeGeoDescription::findNextPlaneEntrance(  TVector3 lpoint,  TVector3 ldir, int nextSensorID, float* newpoint ){
	streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::findNextPlaneEntrance()------BEGIN" << std::endl;
	if(newpoint==NULL)
	{
		throw(lcio::Exception("You have passed a NULL pointer to findNextPlane.")); 	
	}
	//initialise direction and location in global telescope coordinates
	double dlPoint[3];
	dlPoint[0] = lpoint[0];	dlPoint[1] = lpoint[1];	dlPoint[2] = lpoint[2];
	double dlDir[3];
	dlDir[0] = ldir[0];	dlDir[1] = ldir[1];	dlDir[2] = ldir[2];

	_geoManager->InitTrack( dlPoint, dlDir );

	TGeoNode *node = _geoManager->GetCurrentNode( ); //Return the volume i.e 'node' that contains that point.
	Int_t inode =  node->GetIndex();
	Int_t stepNumber=0;

	streamlog_out( DEBUG0 ) << "findNextPlaneEntrance node: " << node << " id: " << inode << std::endl;

	//Keep looping until you have left this plane volume and are at another. Note FindNextBoundaryAndStep will only take you to the next volume 'node' it will not enter it.
	while(( node = _geoManager->FindNextBoundaryAndStep() )){
		inode = node->GetIndex();
		const double* point = _geoManager->GetCurrentPoint(); //This will be the new global coordinates after the move
		const double* dir   = _geoManager->GetCurrentDirection(); //This will be the same direction. Since we will only travel in a straight line.  
		double ipoint[3] ;
		double idir[3]   ;
		//Here we set the coordinates and move into the volume in the z direction.
		for(int ip=0;ip<3;ip++){
			ipoint[ip] = point[ip];
			idir[ip]   = dir[ip];
			if(ip==2){ 
				ipoint[ip]+=0.001 ; // assumption !!! step by one um into the new volume // new volume is thicker than 1 um
			}
			newpoint[ip] = static_cast<float> (ipoint[ip]);
		}  
		int sensorID = getSensorID(newpoint); 

		_geoManager->SetCurrentPoint( ipoint);
		_geoManager->SetCurrentDirection( idir);

		streamlog_out( DEBUG0 ) << "Loop number: " << stepNumber  << ". Index of next boundary: " << inode << ". Current global point: " << ipoint[0]  << " " << ipoint[1] << " " << ipoint[2]  << " sensorID: " << sensorID << ". Input of expect next sensor: " << nextSensorID << std::endl;
		streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::findNextPlaneEntrance()------END" << std::endl;

		if( sensorID == nextSensorID ){
			return true;
		}
		//We return false to say we have not found intersection on the plane.
		if(stepNumber == 10){
			return false;
		}
		stepNumber++;     
	}

}

/**
* Find closest surface intersected by the track and return the position
*/
bool EUTelGeometryTelescopeGeoDescription::findIntersectionWithCertainID(	float x0, float y0, float z0, 
										float px, float py, float pz, 
										float beamQ, int nextPlaneID, float outputPosition[],
										TVector3& outputMomentum, float& arcLength, int& newNextPlaneID)
{
	//positions are in mm
	TVector3 trkVec(x0,y0,z0);
	TVector3 pVec(px,py,pz);

	//Assuming uniform magnetic field running along X direction. 
	//Why do we need this assumption? Equations of motion do not seem to dictate this.
	gear::Vector3D vectorGlobal( x0, y0, z0 );      
	const gear::BField& B	= geo::gGeometry().getMagneticField();
	const double bx         = B.at( vectorGlobal ).x();
	const double by         = B.at( vectorGlobal ).y();
	const double bz         = B.at( vectorGlobal ).z();

	//B field is in units of Tesla
	TVector3 hVec(bx,by,bz);
	const double H = hVec.Mag();

	//Calculate track momentum from track parameters and fill some useful variables
	const double p = pVec.Mag();
	//This is a constant used in the derivation of this equation. This is the distance light travels in a nano second    
	const double constant =  -0.299792458; 
	const double mm = 1000;
	const double combineConstantsAndMagneticField = (constant*beamQ*H)/mm;
	const double rho = combineConstantsAndMagneticField/p; 
				      
	//Determine geometry of sensor to be used to determine the point of intersection.//////////////////////////////////////
	TVector3 norm = geo::gGeometry().siPlaneNormal( nextPlaneID  );       
	streamlog_out (DEBUG5) << "The normal of the plane is (x,y,z): "<<norm[0]<<","<<norm[1]<<","<<norm[2]<< std::endl;
	TVector3 sensorCenter( geo::gGeometry().siPlaneXPosition( nextPlaneID  ), geo::gGeometry().siPlaneYPosition( nextPlaneID  ), geo::gGeometry().siPlaneZPosition( nextPlaneID  ) );
	streamlog_out (DEBUG5) << "The sensor centre (x,y,z): "<<sensorCenter[0]<<","<<sensorCenter[1]<<","<<sensorCenter[2] << std::endl;
	TVector3 delta = (trkVec - sensorCenter);//Must be in mm since momentum is.
	TVector3 pVecCrosH = pVec.Cross(hVec.Unit());

/*
	if( streamlog_level(DEBUG5) )
	{
		streamlog_out(DEBUG5) << "-------------------------------------------" << std::endl;
		streamlog_out(DEBUG5) << "Current point (X,Y,Z): " << std::setw(15) << x0  << std::setw(15) << y0 << std::setw(15) << z0 << std::endl;
		streamlog_out(DEBUG5) << "Next PlaneID : " << nextPlaneID << std::endl;
		streamlog_out(DEBUG5) << "Normal vector" << std::endl; norm.Print();
		streamlog_out(DEBUG5) << "Sensor centre" << std::endl; sensorCenter.Print();
      		streamlog_out(DEBUG5) << "P x H vector" << std::endl; pVecCrosH.Print();
		streamlog_out (DEBUG5) << "Rho: " << rho << " P: " << p << " Delta: " << std::endl; delta.Print();
 	}
*/

	//Solution to the plane equation and the curved line intersection will be an 
	//quadratic with the coefficients. The solution is the arc length along the curve
	const double a = -0.5 * rho * ( norm.Dot( pVecCrosH ) ) / p;
	const double b = norm.Dot( pVec ) / p;
	const double c = norm.Dot( delta );

	//Solution must be in femto metres, vector of arc length
	std::vector<double> sol = Utility::solveQuadratic(a,b,c); 
	//solutions are sorted in ascending order, choose minimum arc length
	double solution = ( sol[0] > 0. ) ? sol[0] : ( ( sol[0] < 0. && sol[1] > 0. ) ? sol[1] : -1. );
	if(solution < 0)
	{
		return false;
	}
			
	//Determine the global position from arc length.             
	TVector3 newPos;
	TVector3 newMomentum;
	newPos = EUTelNav::getXYZfromArcLength(trkVec,pVec,beamQ,solution);
	newMomentum = EUTelNav::getXYZMomentumfromArcLength(pVec, trkVec, beamQ, solution);
	outputMomentum[0] = newMomentum[0];
	outputMomentum[1] = newMomentum[1];
	outputMomentum[2] = newMomentum[2];

	//Is the new point within the sensor. If not then we may have to propagate a little bit further to enter.
 	const float pos[3] = {newPos[0],newPos[1],newPos[2]};
	int sensorIDCheck = getSensorID(pos); 
	bool foundIntersection = false;
	if(sensorIDCheck == nextPlaneID){
		streamlog_out( DEBUG3 ) << "INTERSECTION FOUND! " << std::endl;
		foundIntersection = true;
		outputPosition[0] = newPos[0];
		outputPosition[1] = newPos[1];
		outputPosition[2] = newPos[2];
	}
	if(!foundIntersection)
	{
		foundIntersection = findNextPlaneEntrance( newPos,  newMomentum, nextPlaneID, outputPosition);
	}
	if(!foundIntersection)
	{
		foundIntersection = findNextPlaneEntrance( newPos,  -newMomentum, nextPlaneID, outputPosition);
	}

	arcLength = solution;		

	//We have still not found intersection
	if(!foundIntersection)
	{
		streamlog_out( DEBUG3 ) << "FINAL: NO INTERSECTION FOUND. " << std::endl;
		return false;
	}
	else
	{
		newNextPlaneID = nextPlaneID;
		return true;
	}
}

//This function will intake position and direction. Then using the gear file and magnetic field will output position and sensor ID in correct order of intersection. 
//We need to introduce the idea of:
//sensitive volume => data and state to be created
//scatter volume => state only to be created
//volume => This will cause scattering but no state is to be created on this volume
//At the moment everything in gear file is assumed to be sensitive volume

//std::map<int,double> EUTelGeometryTelescopeGeoDescription::UsingStateReturnAllVolumesIntersected(){}

void EUTelGeometryTelescopeGeoDescription::updateSiPlanesLayout() {
 streamlog_out( MESSAGE1 ) << "EUTelGeometryTelescopeGeoDescription::updateSiPlanesLayout() --- START ---- " << std::endl;

    gear::SiPlanesParameters*    siplanesParameters = const_cast< gear::SiPlanesParameters*> (&( _gearManager->getSiPlanesParameters()));
    gear::SiPlanesLayerLayout*  siplanesLayerLayout = const_cast< gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));

    // data member::
    _nPlanes = siplanesLayerLayout->getNLayers(); 
 
    // create an array with the z positions of each layer
    for(size_t iPlane = 0; iPlane < _nPlanes; iPlane++)
    {
        int sensorID =  _sensorIDVec.at(iPlane);
        
	siplanesLayerLayout->setLayerPositionX( iPlane, siPlaneXPosition(sensorID) );
        siplanesLayerLayout->setLayerPositionY(  iPlane, siPlaneYPosition(sensorID) );
        siplanesLayerLayout->setLayerPositionZ(  iPlane, siPlaneZPosition(sensorID) );
        siplanesLayerLayout->setLayerRotationZY( iPlane, siPlaneXRotation(sensorID) );
        siplanesLayerLayout->setLayerRotationZX( iPlane, siPlaneYRotation(sensorID) );
        siplanesLayerLayout->setLayerRotationXY( iPlane, siPlaneZRotation(sensorID) );
    }


    // ------- add to GearMgr ----
    if( _gearManager != 0 )
    {
      _gearManager->setSiPlanesParameters( siplanesParameters ) ;

    }
 streamlog_out( MESSAGE1 ) << "EUTelGeometryTelescopeGeoDescription::updateSiPlanesLayout() --- OVER ---- " << std::endl;
}


void EUTelGeometryTelescopeGeoDescription::updateTrackerPlanesLayout() {

    streamlog_out( MESSAGE1 ) << "EUTelGeometryTelescopeGeoDescription::updateTrackerPlanesLayout() --- START ---- " << std::endl;

    gear::TrackerPlanesParameters*  trackerplanesParameters  = const_cast< gear::TrackerPlanesParameters*>  (&( _gearManager->getTrackerPlanesParameters() ));
    gear::TrackerPlanesLayerLayout* trackerplanesLayerLayout = const_cast< gear::TrackerPlanesLayerLayout*> (&(  trackerplanesParameters->getTrackerPlanesLayerLayout() ));
    
    trackerplanesParameters->setLayoutID( getSiPlanesLayoutID() ) ;
 


    // create an array with the z positions of each layer
    int nLayers = trackerplanesLayerLayout->getNLayers();
    for (int iLayer = 0; iLayer < nLayers; iLayer++) {
        gear::TrackerPlanesLayerImpl*  trackerplanesLayerImpl = const_cast< gear::TrackerPlanesLayerImpl*>  ( trackerplanesLayerLayout->getLayer( iLayer) );
        int nsensitive =  trackerplanesLayerImpl->getNSensitiveLayers() ;

        gear::TrackerPlanesSensitiveLayerImplVec& vector =  trackerplanesLayerImpl->getSensitiveLayerVec();
       
        for (int iSensLayer = 0; iSensLayer < nsensitive; iSensLayer++) {       

            gear::TrackerPlanesSensitiveLayerImpl& sensitiveLayer = vector.at(iSensLayer);
 
            for( size_t iplane = 0; iplane < _sensorIDVec.size(); iplane++ )
	    {
              int sensorID =  _sensorIDVec.at(iplane);
            
	      if( sensitiveLayer.getID() !=  _sensorIDVec.at( iplane) ) continue;  
              
	      sensitiveLayer.setPositionX( siPlaneXPosition(sensorID) );
              sensitiveLayer.setPositionY( siPlaneYPosition(sensorID) );
              sensitiveLayer.setPositionZ( siPlaneZPosition(sensorID) );

              sensitiveLayer.setRotationZY( siPlaneXRotation(sensorID) );
              sensitiveLayer.setRotationZX( siPlaneYRotation(sensorID) );
              sensitiveLayer.setRotationXY( siPlaneZRotation(sensorID) );
            }
        }
    }

    // ------- add to GearMgr ----
    if( _gearManager != 0 )
    {
   	 _gearManager->setTrackerPlanesParameters( trackerplanesParameters ) ;
    }
    streamlog_out( MESSAGE1 ) << "EUTelGeometryTelescopeGeoDescription::updateTrackerPlanesLayout() --- OVER ---- " << std::endl;
}

void EUTelGeometryTelescopeGeoDescription::updateGearManager()
{

	if( _siPlanesDefined )
	{
		updateSiPlanesLayout();
	}
	else if( _telPlanesDefined )
	{
		updateTrackerPlanesLayout();
	}
}

void EUTelGeometryTelescopeGeoDescription::writeGEARFile(std::string filename)
{
	updateGearManager();
	gear::GearXML::createXMLFile(marlin::Global::GEAR, filename);
}
