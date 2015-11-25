/* 
 * File:   EUTelGeometryTelescopeGeoDescription.cpp
 * 
 */
#include "EUTelGeometryTelescopeGeoDescription.h"

// C++
#include <algorithm>
#include <string>
#include <cstring>
#include <cmath>
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

#define  SCATTER_IDENTIFIER 314

using namespace eutelescope;
using namespace geo;

unsigned EUTelGeometryTelescopeGeoDescription::_counter = 0;

/**TODO: Replace me: NOP*/
EUTelGeometryTelescopeGeoDescription& EUTelGeometryTelescopeGeoDescription::getInstance( gear::GearMgr* _g ) {
	static  EUTelGeometryTelescopeGeoDescription instance;
	unsigned i = EUTelGeometryTelescopeGeoDescription::_counter;
	
	//do it only once!
	if( i < 1 ) {
		instance.setGearManager(_g);
		instance.readGear();
	}
	instance.counter();
	return instance;
}

//Note  that to determine these axis we MUST use the geometry class after initialisation. By this I mean directly from the root file create.
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneNormal( int planeID )
{
	std::map<int, TVector3>::iterator mapIt = _planeNormalMap.find(planeID);
	if( mapIt != _planeNormalMap.end() ) {
		return mapIt->second;
	} else {
		std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
		if( it != _sensorIDVec.end() ) {
			std::array<double,3> const zAxisLocal {{0,0,1}};
			std::array<double,3> zAxisGlobal; 
			local2MasterVec(planeID, zAxisLocal, zAxisGlobal); 
			TVector3 normVec(zAxisGlobal.data());
			_planeNormalMap[planeID] = normVec;
			return normVec;
		} else {
			std::stringstream ss;
			ss << planeID;
			std::string errMsg = "EUTelGeometryTelescopeGeoDescription::siPlaneNormal: Could not find planeID: " + ss.str();
			throw InvalidGeometryException(errMsg);
		}
	}
}

/**TODO: Replace me: NOP*/
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneXAxis( int planeID ) {
	std::map<int, TVector3>::iterator mapIt = _planeXMap.find(planeID);
	if( mapIt != _planeXMap.end() ) {
		return mapIt->second;
	} else {
		std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
		if( it != _sensorIDVec.end() ) {
			std::array<double,3> const xAxisLocal {{1,0,0}};
			std::array<double,3> xAxisGlobal; 
			local2MasterVec(planeID, xAxisLocal, xAxisGlobal); 
			TVector3 xVec(xAxisGlobal.data());
			_planeXMap[planeID] = xVec;
			return xVec;
		} else {
			std::stringstream ss;
			ss << planeID;
			std::string errMsg = "EUTelGeometryTelescopeGeoDescription::siPlaneXAxis: Could not find planeID: " + ss.str();
			throw InvalidGeometryException(errMsg);
		}
	}
}

/**TODO: Replace me: NOP*/
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneYAxis( int planeID ) {
	std::map<int, TVector3>::iterator mapIt = _planeYMap.find(planeID);
	if( mapIt != _planeYMap.end() ) {
		return mapIt->second;
	} else {
		std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
		if( it != _sensorIDVec.end() ) {
			std::array<double,3> const yAxisLocal {{0,1,0}};
			std::array<double,3> yAxisGlobal; 
			local2MasterVec(planeID, yAxisLocal, yAxisGlobal); 
			TVector3 yVec(yAxisGlobal.data());
			_planeYMap[planeID] = yVec;
			return yVec;
		} else {
			std::stringstream ss;
			ss << planeID;
			std::string errMsg = "EUTelGeometryTelescopeGeoDescription::siPlaneYAxis: Could not find planeID: " + ss.str(); 
			throw InvalidGeometryException(errMsg);
		}
	}
}

void EUTelGeometryTelescopeGeoDescription::readSiPlanesLayout() {
	// sensor-planes in geometry navigation:
	_siPlanesParameters = const_cast< gear::SiPlanesParameters*> (&( _gearManager->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast< gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));
	_nPlanes = _siPlanesLayerLayout->getNLayers(); 

	//read the geoemtry names from the "Geometry" StringVec section of the gear file
	lcio::StringVec geometryNameParameters;

	try {
		geometryNameParameters  =  _siPlanesParameters->getStringVals("Geometry");
	} catch(gear::UnknownParameterException e) {
		streamlog_out(MESSAGE6) << "No Geometry field found in GEAR file, assuming CAST for all planes" << std::endl;
		for(size_t i = 0; i < _nPlanes; i++) {
			geometryNameParameters.push_back("CAST");
		}
	}

	setSiPlanesLayoutID( _siPlanesParameters->getSiPlanesID() ) ;

	// create an array with the z positions of each layer
	for (size_t iPlane = 0; iPlane < _nPlanes; iPlane++) {
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

		//GEAR uses mm wheras TGeo will use cm
		thisPlane.radLength	= _siPlanesLayerLayout->getSensitiveRadLength(iPlane)/10;

		_planeSetup[_siPlanesLayerLayout->getID(iPlane)] = thisPlane;
	}

	_sensorIDVec.clear();

	for(int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++) {
		int sensorID = _siPlanesLayerLayout->getID(iPlane);
		_sensorIDVec.push_back(sensorID);
	}
	_nPlanes = _siPlanesParameters->getSiPlanesNumber();
	std::sort(_sensorIDVec.begin(), _sensorIDVec.end(), doCompare(*this) );	
}

void EUTelGeometryTelescopeGeoDescription::readTrackerPlanesLayout() {
	// sensor-planes in geometry navigation:
	_trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&( _gearManager->getTrackerPlanesParameters()));
	_trackerPlanesLayerLayout = const_cast< gear::TrackerPlanesLayerLayout*> (&(_trackerPlanesParameters->getTrackerPlanesLayerLayout()));

	setSiPlanesLayoutID( _trackerPlanesParameters->getLayoutID() ) ;

	_sensorIDVec.clear();
	
	//should be filled based on the length of the sensor vector after the loop
	_nPlanes = 0; 

	// create an array with the z positions of each layer
	int nLayers = _trackerPlanesLayerLayout->getNLayers();
	for (int iLayer = 0; iLayer < nLayers; iLayer++) {
		gear::TrackerPlanesLayerImpl* _trackerPlanesLayerImpl = const_cast<gear::TrackerPlanesLayerImpl*>(_trackerPlanesLayerLayout->getLayer( iLayer));

		int nsensitive = _trackerPlanesLayerImpl->getNSensitiveLayers();
		gear::TrackerPlanesSensitiveLayerImplVec& vector = _trackerPlanesLayerImpl->getSensitiveLayerVec();

		for(int iSensLayer = 0; iSensLayer < nsensitive; iSensLayer++) {
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

			//GEAR uses mm wheras TGeo will use cm
			thisPlane.radLength	= sensitiveLayer.getRadLength()/10;

			_planeSetup[sensorID] = thisPlane;

			_sensorIDVec.push_back(sensorID);
			streamlog_out(DEBUG1) << " iter: " << _sensorIDVec.at( _sensorIDVec.size()-1 ) << " " << sensorID << " " << sensitiveLayer.getInfo() .c_str() << std::endl; 
		}
	}
	std::sort(_sensorIDVec.begin(), _sensorIDVec.end(), doCompare(*this) );
	_nPlanes =  _sensorIDVec.size(); 
}

EUTelGeometryTelescopeGeoDescription::EUTelGeometryTelescopeGeoDescription():
_gearManager( marlin::Global::GEAR ),
_siPlanesDefined(false),
_telPlanesDefined(false),
_siPlanesParameters(nullptr),
_siPlanesLayerLayout(nullptr),
_trackerPlanesParameters(nullptr),
_trackerPlanesLayerLayout(nullptr),
_sensorIDVec(),
_nPlanes(0),
_isGeoInitialized(false),
_geoManager(nullptr)
{
	//Set ROOTs verbosity to only display error messages or higher (so info will not be streamed to stderr)
	gErrorIgnoreLevel =  kError;  
	//Pixel Geometry manager creation
	_pixGeoMgr = new EUTelGenericPixGeoMgr();
}

void EUTelGeometryTelescopeGeoDescription::readGear() {
	if( _gearManager == nullptr ) {
		streamlog_out(ERROR2) << "The GearMgr is not available, for an unknown reason." << std::endl;
		throw eutelescope::InvalidGeometryException("GEAR manager is not initialised");
	} 
	try {
		_siPlanesParameters = const_cast< gear::SiPlanesParameters*> (&(_gearManager->getSiPlanesParameters()));
		streamlog_out(MESSAGE1)  << "gear::SiPlanes : " << _siPlanesParameters << std::endl;
		_siPlanesDefined = true;
	} catch(...) {
		streamlog_out(WARNING)   << "gear::SiPlanes NOT found " << std::endl;
	}
	try {
		_trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&(_gearManager->getTrackerPlanesParameters()));
		streamlog_out(MESSAGE1)  << "gear::TrackerPlanes : " << _trackerPlanesParameters << std::endl;
		_telPlanesDefined = true;
	} catch(...) {
		streamlog_out(WARNING)   << "gear::TrackerPlanes NOT found "  << std::endl;
	}

	if( _siPlanesDefined ) {
		readSiPlanesLayout();
	} else if( _telPlanesDefined ) {
		readTrackerPlanesLayout();
	} else {
		streamlog_out(ERROR5) << "Your GEAR file neither contains SiPlanes nor TrackerPlanes and thus is not valid" << std::endl;
		throw eutelescope::InvalidGeometryException("GEAR file invalid, does not contain SiPlanes nor TrackerPlanes");
	}
}

EUTelGeometryTelescopeGeoDescription::~EUTelGeometryTelescopeGeoDescription() {
	delete _geoManager;
	_geoManager = nullptr;
	delete _pixGeoMgr;
	_pixGeoMgr = nullptr;
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
void EUTelGeometryTelescopeGeoDescription::translateSiPlane2TGeo(TGeoVolume* pvolumeWorld, int SensorId ) {
	double xc, yc, zc;   // volume center position 
	double alpha, beta, gamma;
	double rotRef1, rotRef2, rotRef3, rotRef4;

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
	float determinant = rotRef1*rotRef4 - rotRef2*rotRef3  ;
	if(determinant==1 or determinant==-1) { 
		streamlog_out(DEBUG5) << "SensorID: " << SensorId << ". Determinant =  " <<determinant <<"  This is the correct determinate for this transformation." << std::endl;   
	} else {
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
	TGeoRotation* pMatrixRotRefCombined = new TGeoRotation();
	//We have to ensure that we retain a right handed coordinate system, i.e. if we only flip the x or y axis, we have to also flip the z-axis. If we flip both we have to flip twice.	
	double integerRotationsAndReflections[9]={rotRef1,rotRef2,0,rotRef3,rotRef4,0,0,0, determinant};
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
	
	// Construct object medium. Required for radiation length determination
	// assume SILICON, though all information except of radiation length is ignored
	double a       = 28.085500;     
	double z       = 14.000000;
	double density = 2.330000;
	double radl    = siPlaneRadLength( SensorId );
	double absl    = 45.753206;
	std::string stMatName = "materialSensor";
	stMatName.append( strId.str() );
	TGeoMaterial* pMat = new TGeoMaterial( stMatName.c_str(), a, z, density, -radl, absl );
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

	// Geometry navigation package requires following names for objects that have an ID  name:ID
	std::string stVolName = "volume_SensorID:";
	stVolName.append( strId.str() );

	_planePath.insert( std::make_pair(SensorId, "/volume_World_1/"+stVolName+"_1") );

	TGeoVolume* pvolumeSensor = new TGeoVolume( stVolName.c_str(), pBoxSensor, pMed );
	pvolumeSensor->SetVisLeaves( kTRUE );
	pvolumeWorld->AddNode(pvolumeSensor, 1/*(SensorId)*/, combi);

	//this line tells the pixel geometry manager to load the pixel geometry into the plane			
	streamlog_out(DEBUG1) << " sensorID: " << SensorId << " " << stVolName << std::endl;   
	std::string name = geoLibName(SensorId);

	if( name == "CAST" ) {
		_pixGeoMgr->addCastedPlane( SensorId, siPlaneXNpixels(SensorId), siPlaneYNpixels(SensorId), siPlaneXSize(SensorId), siPlaneYSize(SensorId), siPlaneZSize(SensorId), siPlaneRadLength(SensorId), stVolName);
	} else {
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
	if( _isGeoInitialized ) {
		streamlog_out( WARNING3 ) << "EUTelGeometryTelescopeGeoDescription: Geometry already initialized, using old initialization" << std::endl;
		return;
	} else {
    		_geoManager = new TGeoManager("Telescope", "v0.1");
	}

	if( !_geoManager ) {
		streamlog_out( ERROR3 ) << "Can't instantiate ROOT TGeoManager " << std::endl;
		return;
	}
   
    
    // Create top world volume containing telescope geometry
    
    // Create air mixture
    // see http://pdg.lbl.gov/2013/AtomicNuclearProperties/HTML_PAGES/104.html
    double air_density = 1.2e-3;         // g/cm^3
    double air_radlen  = 36.62;          // g/cm^2 //Must be -ve to use this value rather than internal one calculated.

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
   }
    _geoManager->CloseGeometry();
    _isGeoInitialized = true;
    // Dump ROOT TGeo object into file
    if ( dumpRoot ) _geoManager->Export( geomName.c_str() );
    return;
}

Eigen::Matrix3d EUTelGeometryTelescopeGeoDescription::rotationMatrixFromAngles(int sensorID) {
	return rotationMatrixFromAngles( (long double)siPlaneXRotationRadians(sensorID), (long double)siPlaneYRotationRadians(sensorID), (long double)siPlaneZRotationRadians(sensorID) );
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::globalXAxis(int sensorID) {
	Eigen::Vector3d xAxis(1,0,0);
	return rotationMatrixFromAngles(sensorID)*xAxis;
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::globalYAxis(int sensorID) {
	Eigen::Vector3d yAxis(0,1,0);
	return rotationMatrixFromAngles(sensorID)*yAxis;
}

/** Returns the rotation matrix for given angles
 *  It alpha rotation is around initial X axis, then beta around the new Y' axis
 *  and finally the gamam rotation around the new Z'' axis */
Eigen::Matrix3d EUTelGeometryTelescopeGeoDescription::rotationMatrixFromAngles(long double alpha, long double beta, long double gamma) {
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

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getRotationAnglesFromMatrix(Eigen::Matrix3d rotMat) {
	long double alpha = asin((long double)(-rotMat(1,2)));
	long double cosA = cos(alpha);

	long double beta = asin((long double)(rotMat(0,2)/cosA));
	long double gamma = asin((long double)(rotMat(1,0)/cosA));

	Eigen::Vector3d vec;

	vec << alpha, beta, gamma;
	return vec;
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getOffsetVector(int sensorID) {
	Eigen::Vector3d offsetVec;
	offsetVec << siPlaneXPosition(sensorID), siPlaneYPosition(sensorID), siPlaneZPosition(sensorID); 
	return offsetVec;
}

Eigen::Matrix3d EUTelGeometryTelescopeGeoDescription::getFlipMatrix(int sensorID) {
	Eigen::Matrix3d flipMat;
	flipMat << 	siPlaneRotation1(sensorID),	siPlaneRotation2(sensorID),	0,
			siPlaneRotation3(sensorID), 	siPlaneRotation4(sensorID),	0,
			0,				0,				1;	
	return flipMat;
}

int EUTelGeometryTelescopeGeoDescription::getSensorIDFromManager() {
    std::vector<std::string> split;
 
    int sensorID = -999;

  	int levelStart =	geo::gGeometry()._geoManager->GetLevel();
    while( _geoManager->GetLevel() > 0 ) { 
      const char* volName = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );
      split = Utility::stringSplit( std::string( volName ), "/", false);
      if ( split.size() > 0 && split[0].length() > 16 && (split[0].substr(0,16) == "volume_SensorID:") ) {
         int strLength = split[0].length(); 
         sensorID = strtol( (split[0].substr(16, strLength )).c_str(), NULL, 10 );
         break;
      }
      _geoManager->CdUp();	////////////////////////////////////////THIS NEEDS TO BE FIXED. If partice falls in the pixel volume and to find sensor ID you need to be on the sensor volume
    }
  	int levelEnd =	geo::gGeometry()._geoManager->GetLevel();

//		std::cout <<" node level end : " << _geoManager->GetLevel() <<std::endl;

		//Must return the manager pointing to the same node before we looked for the sensorID
		for(int i =0 ; i < (levelStart - levelEnd); i++ ){
			geo::gGeometry()._geoManager->CdDown(0);
		}
//		std::cout <<" node level re : " << _geoManager->GetLevel() <<std::endl;

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
    _geoManager->cd( _planePath[sensorID].c_str() );
    _geoManager->GetCurrentNode()->LocalToMaster( localPos, globalPos );
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
    _geoManager->cd( _planePath[sensorID].c_str() );
    _geoManager->GetCurrentNode()->MasterToLocal( globalPos, localPos );
}

/**
 * Vector coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param globalVec (x,y,z) in global coordinate system
 * @param localVec (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::local2MasterVec( int sensorID, const double localVec[], double globalVec[] ) {
    _geoManager->cd( _planePath[sensorID].c_str() );
    _geoManager->GetCurrentNode()->LocalToMasterVect( localVec, globalVec );
}

/**
 * Vector coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param globalVec (x,y,z) in global coordinate system
 * @param localVec (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::master2LocalVec( int sensorID, const double globalVec[], double localVec[] ) {
    _geoManager->cd( _planePath[sensorID].c_str() );
    _geoManager->GetCurrentNode()->MasterToLocalVect( globalVec, localVec );
}

void EUTelGeometryTelescopeGeoDescription::local2Master( int sensorID, std::array<double,3> const & localPos, std::array<double,3>& globalPos) {
	this->local2Master(sensorID, localPos.data(), globalPos.data());
}
void EUTelGeometryTelescopeGeoDescription::master2Local(int sensorID, std::array<double,3> const & globalPos, std::array<double,3>& localPos) {
	this->master2Local(sensorID, globalPos.data(), localPos.data());
}
void EUTelGeometryTelescopeGeoDescription::local2MasterVec( int sensorID, std::array<double,3> const & localVec, std::array<double,3>& globalVec) {
	this->local2MasterVec(sensorID, localVec.data(), globalVec.data());
}
void EUTelGeometryTelescopeGeoDescription::master2LocalVec( int sensorID, std::array<double,3> const & globalVec, std::array<double,3>& localVec) {
	this->master2LocalVec(sensorID, globalVec.data(), localVec.data());
}

/**
 * Local-to-Global coordinate transformation matrix.
 * Corresponding volume is determined automatically.
 * 
 * @param globalPos (x,y,z) in global coordinate system
 * @return 
 */
const TGeoHMatrix* EUTelGeometryTelescopeGeoDescription::getHMatrix( const double globalPos[] ) {
    _geoManager->FindNode( globalPos[0], globalPos[1], globalPos[2] );    
    const TGeoHMatrix* globalH = _geoManager->GetCurrentMatrix();
	//	if(streamlog_out(DEBUG2)){
  //  streamlog_out(DEBUG2) << "Transformation matrix " << std::endl;
	//	globalH->Print();
	//	}
		
    return globalH;
}

TMatrixD EUTelGeometryTelescopeGeoDescription::getRotMatrix( int sensorID ) {
	streamlog_out(DEBUG0) << "EUTelGeometryTelescopeGeoDescription::getRotMatrix()--------BEGIN " << std::endl;
	std::array<double,3> const local {{0,0,0}};
	std::array<double,3> global;
//	std::cout << "Sensor ID " << sensorID << std::endl;
	TMatrixD TRotMatrix(3,3);
	if(sensorID != SCATTER_IDENTIFIER) {
		local2Master( sensorID, local, global );
		_geoManager->FindNode( global[0], global[1], global[2] );    
		const TGeoHMatrix* globalH = _geoManager->GetCurrentMatrix();
		const double* rotMatrix = globalH->GetRotationMatrix();
		TRotMatrix[0][0] = *rotMatrix; TRotMatrix[0][1] = *(rotMatrix+1);TRotMatrix[0][2] = *(rotMatrix+2);
		TRotMatrix[1][0] = *(rotMatrix+3); TRotMatrix[1][1] = *(rotMatrix+4);TRotMatrix[1][2] = *(rotMatrix+5);
		TRotMatrix[2][0] = *(rotMatrix+6); TRotMatrix[2][1] = *(rotMatrix+7);TRotMatrix[2][2] = *(rotMatrix+8);
	} else {
		TRotMatrix.UnitMatrix();
	}
	//	std::cout<< "Here is the first element of rotation matrix: " << TRotMatrix[0][0]<<std::endl;
	//	std::cout<< "Here is the last element of rotation matrix: " << TRotMatrix[2][2]<<std::endl;

    return TRotMatrix;
    streamlog_out(DEBUG0) << "EUTelGeometryTelescopeGeoDescription::getRotMatrix()----END " << std::endl;

}
int EUTelGeometryTelescopeGeoDescription::getSensorID( float const globalPos[] ) const {
    double pos[3] = {globalPos[0],globalPos[1],globalPos[2]};
    return getSensorID(pos);
}
/** Determine id of the sensor in which point is locate
 *  * 
 *  * @param globalPos 3D point in global reference frame
 *  * @return sensorID or -999 if the point in outside of sensor volume
 *  */
int EUTelGeometryTelescopeGeoDescription::getSensorID( double const globalPos[] ) const {
    streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::getSensorID() " << std::endl;
    const float constPos[3] = {globalPos[0],globalPos[1],globalPos[2]};
    _geoManager->FindNode( constPos[0], constPos[1], constPos[2] );

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
    _geoManager->CdUp();  ////////////////////////////////////////THIS NEEDS TO BE FIXED. If partice falls in the pixel volume and to find sensor ID you need to be on the sensor volume
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
    return sensorID;
}

int EUTelGeometryTelescopeGeoDescription::getSensorID(std::array<double,3> const globalPos) const {
	return this->getSensorID(globalPos.data());
}
int EUTelGeometryTelescopeGeoDescription::getSensorID(std::array<float,3> const globalPos) const {
	return this->getSensorID(globalPos.data());
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
 * 
 * @return radiation length in units of X0
 */

float EUTelGeometryTelescopeGeoDescription::findRad( const std::map<int,int>& sensorIDToZOrderWithoutExcludedPlanes, const double globalPosStart[], const double globalPosFinish[], std::map< const int, double> &sensors, 	std::map< const int, double> &air ){
    streamlog_out(DEBUG5) << "/////////////////////////////////////////////////////////////////////////////////////////////////// " << std::endl;
    streamlog_out(DEBUG5) << "/////////////////////////////////////////////////////////////////////////////////////////////////// " << std::endl;
    streamlog_out(DEBUG5) << "              CALCULATING THE TOTAL RADIATION LENGTH BETWEEN TWO POINTS.                            " << std::endl;
    streamlog_out(DEBUG5) << "                                 POINTS GO FROM:                            " << std::endl;
    streamlog_out(DEBUG5) << "              ("<< globalPosStart[0] << ","  << globalPosStart[1]<<","<<globalPosStart[2]<<")"<< "-------------------------------------> ("<< globalPosFinish[0] << ","  << globalPosFinish[1]<<","<<globalPosFinish[2]<<")" << std::endl;
    const double mm2cm = 0.1;
    bool foundFirstPlane=false;
    double blockEnd=0;
    double total=0;
    int sensorLeftSide=0;

    //get direction from positions.
    const double stepLength =TMath::Sqrt( ( globalPosFinish[0] - globalPosStart[0] )*( globalPosFinish[0] - globalPosStart[0] ) +
                               ( globalPosFinish[1] - globalPosStart[1] )*( globalPosFinish[1] - globalPosStart[1] ) +
                               ( globalPosFinish[2] - globalPosStart[2] )*( globalPosFinish[2] - globalPosStart[2] ));

    const double xp  = ( globalPosFinish[0] - globalPosStart[0] )/stepLength;
    const double yp  = ( globalPosFinish[1] - globalPosStart[1] )/stepLength;
    const double zp  = ( globalPosFinish[2] - globalPosStart[2] )/stepLength;
    //We get the object we are in currently on this track and begin to loop to each plane 
    gGeoManager->InitTrack( globalPosStart[0]/*mm*/, globalPosStart[1]/*mm*/, globalPosStart[2]/*mm*/, xp, yp, zp ); //This the start point and direction
    TGeoNode *nextnode = gGeoManager->GetCurrentNode( );
    while ( nextnode ) {
        int sensorID = getSensorIDFromManager();
        //If not in the first plane then look forward.
        if(sensorID == sensorIDToZOrderWithoutExcludedPlanes.at(0) or foundFirstPlane ){ //We want to make sure to add radiation length from first plane to last plane only;
            foundFirstPlane = true;
        }else{
            //nextnode is the next found and we can get the step using GetStep. stepLength is the max distance to travel before we find another node. 
            nextnode = gGeoManager->FindNextBoundaryAndStep( stepLength /*mm*/ );  
        }
        //Don't do anything until we find the first sensor.
        if(foundFirstPlane){
            TGeoMedium *med = NULL;
            if ( nextnode ) med = nextnode->GetVolume()->GetMedium();
            else return 0.; //We return 0 to get rid of the track but not the event.
            double radlen = med->GetMaterial()->GetRadLen() /*cm*/;
            double lastrad = 1. / radlen * mm2cm; //calculate 1/radiationlength per cm. This will transform radlen to mm
            nextnode = gGeoManager->FindNextBoundaryAndStep( stepLength /*mm*/ );  
            double snext  = gGeoManager->GetStep() /*mm*/; //This will output the distance traveled by FindNextBoundaryAndStep
            double rad = 0; //This is the calculated (rad per distance x distance)
            double delta = 0.01;//This is the minimum block size 
            streamlog_out(DEBUG5)<<std::endl <<std::endl  << "DECISION: Step size over min?  "  <<" Block width: " << snext << " delta: " << delta  << std::endl;
            if(snext < delta){
               streamlog_out(DEBUG5) << "INCREASE TO MINIMUM DISTANCE!" << std::endl;
                snext = delta;
                double pt[3];
                memcpy( pt, gGeoManager->GetCurrentPoint(), 3 * sizeof (double) ); //Get global position
                const double *dir = gGeoManager->GetCurrentDirection();//Direction vector
                for ( Int_t i = 0; i < 3; i++ ) pt[i] += delta * dir[i]; //Move the current point slightly in the direction of motion. 
                nextnode = gGeoManager->FindNode( pt[0], pt[1], pt[2] );//Move to new node where we will begin to look for more radiation length   
                rad=lastrad*snext; //Calculate radiation length for the increased block.
                blockEnd += snext;
           }else{
                streamlog_out(DEBUG5) << "OVERMAX!" << std::endl;
                blockEnd += snext;
                rad = lastrad*snext; //This is the calculated (rad per distance x distance)
           }
            total = total + rad;
            streamlog_out(DEBUG5) << "NEW BLOCK:SensorID: " << sensorID   <<" Block width: " << snext << " Radiation total/per unit length: " << rad<<"/"<< 1.0/lastrad << " Block end position: " << blockEnd  << std::endl;
            streamlog_out(DEBUG5) << "DECISION: Where should block be placed?"  << std::endl;

            //Now we have the block. We place it in the planes or in the air if excluded.
            //Work Flow: Check we are at end. If not then set radiation length to sensor if included.Else attach the radiation length last one included and found. 
            if(sensorID != sensorIDToZOrderWithoutExcludedPlanes.at( sensorIDToZOrderWithoutExcludedPlanes.size()-1 )){
                if(sensorIDToZOrderWithoutExcludedPlanes.find(sensorID) !=  sensorIDToZOrderWithoutExcludedPlanes.end()){
                    sensors[sensorID] = sensors[sensorID] +  rad;
                    sensorLeftSide =sensorID;
                    air[sensorID] = 0 ;
                    streamlog_out( DEBUG5 ) << "Placed in sensor:: " << sensorID <<std::endl;
                }else{
                    air[sensorLeftSide] = air[sensorLeftSide] +  rad;
                    streamlog_out( DEBUG5 ) << "Placed in air after sensor:: " << sensorLeftSide <<std::endl;
                }
            }else{//IF WE HAVE REACHED THE FINAL SENSOR THEN DO NOT LOOK FOR ANYMORE RADIATION LENGTH
                sensors[sensorID] = sensors[sensorID] + rad;
                streamlog_out( DEBUG5 ) << "FINISHED: Sensor: " <<sensorID << " Total Rad: " << total <<std::endl;
                return total;
            }
        }
    }
    return 0; //Return 0 if no nextnode to remove track 
}
//This will output the X/X0 of the the full detector system. This is needed to calculate the for each individual scatter the proper correction. 
//Note we can not determine this correction for each scatterer individually since this correction would introduce a non linear term which would be unphysical. 
float EUTelGeometryTelescopeGeoDescription::calculateTotalRadiationLengthAndWeights(const double start[3]  ,const double end[3],  std::map<const int,double> & mapSensor, std::map<const int ,double> & mapAir){
//	streamlog_out(DEBUG1) << "calculateTotalRadiationLength()------------------------------BEGIN" <<std::endl;
//	std::map<const int, double> sensors; //This will store all the sensor scattering.
//	std::map<const int, double> air; //This will store the air directly infront of one plane
//	int lastPlaneID; // = 	geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at( geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()-1 );
////	std::cout<<"Reached here" <<std::endl;
////	std::cout << " Here is the position: " <<_planeSetup[lastPlaneID].zPos  <<std::endl; 
//	const double endUpdate[] = {end[0],end[1],_planeSetup[lastPlaneID].zPos +0.025};//Must make sure we add all silicon.
//	//THIS WILL RETURN THE TOTAL RADIATION LENGTH FOUND AND THE FRACTION FOR AIR AND PLANES.
//	//This will be for non excluded planes. This is sorted in mapWeightsToSensor(...)
//    const std::map<int,int>& sensorIDToZOrderWithoutExcludedPlanes;
//	float perRad =	findRad( sensorIDToZOrderWithoutExcludedPlanes, start,endUpdate, sensors, air ); 
//	//NOW WE REDUCE EXCLUDED PLANES TO DEAD MATERIAL. THIS IS ABSORBED IN THE AIR OF THE PLANES NOT EXCLUDED.
//	//First two with excluded. The last two without.
//	mapWeightsToSensor(sensors,air, mapSensor,mapAir);
//	perRad = perRad +	addKapton(mapSensor);
//	streamlog_out(DEBUG0) << "X/X0 (TOTAL SYSTEM) : " << perRad <<std::endl;
//	bool pass = testOutput(mapSensor, mapAir);
//	if(!pass){
//		return 0;
//	}
////	std::cout << "here the rad" << perRad << std::endl;
//	return perRad;
//	streamlog_out(DEBUG1) << "calculateTotalRadiationLength()------------------------------END" <<std::endl;
    return 1.0;

}
//This function wil not add kapton to excluded planes.
//TO DO: The found planes will be different from the excluded. Since sometimes we will miss a plane if there are two DUT for example. Must keep a note of these since we will be still adding radiation length incorrectly if not accounted for. These track are removed at the moment since some entries of radiation length will zero.
double EUTelGeometryTelescopeGeoDescription::addKapton(std::map<const int, double> & mapSensor){
//	for(unsigned int i = 0 ; i<geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() ; i++){
//		mapSensor[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)] = mapSensor[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)] + 0.0002;
//	}
	return 1.0;// 2*geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()*0.0001;

}


double EUTelGeometryTelescopeGeoDescription::FindRad(Eigen::Vector3d const & startPt, Eigen::Vector3d const & endPt) {

	Eigen::Vector3d track = endPt-startPt;
	double length = track.norm();
	track.normalize();

	double snext;
	Eigen::Vector3d point;
	Eigen::Vector3d direction;
	double epsil = 0.00001;
	double rad    = 0.;
	double propagatedDistance = 0;
	bool reachedEnd = false;

	TGeoMedium* med;
	gGeoManager->InitTrack(startPt(0), startPt(1), startPt(2), track(0), track(1), track(2));
	TGeoNode* nextnode = gGeoManager->GetCurrentNode();

	while(nextnode && !reachedEnd) {
		med = nullptr;
		if (nextnode) med = nextnode->GetVolume()->GetMedium();

		nextnode = gGeoManager->FindNextBoundaryAndStep(length);
		snext  = gGeoManager->GetStep();

		if( propagatedDistance+snext >= length ) {
			snext = length - propagatedDistance;
			reachedEnd = true;
		}
		//snext gets very small at a transition into a next node, in this case we need to manually propagate a small (epsil)
		//step into the direction of propagation. This introduces a small systematic error, depending on the size of epsil as
	    	if(snext < 1.e-8) {
			const double * currDir = gGeoManager->GetCurrentDirection();
			const double * currPt = gGeoManager->GetCurrentPoint();

			direction(0) = currDir[0]; direction(1) = currDir[1]; direction(2) = currDir[2];
			point(0) = currPt[0]; point(1) = currPt[1]; point(2) = currPt[2];

			point = point + epsil*direction;

			gGeoManager->CdTop();
			nextnode = gGeoManager->FindNode(point(0),point(1),point(2));
			snext = epsil;
		}	
		if(med) {
			//ROOT returns the rad length in cm while we use mm, therefore factor of 10
			double radlen = med->GetMaterial()->GetRadLen();
			if (radlen > 1.e-5 && radlen < 1.e10) {
				rad += snext/(radlen*10);
			} 
		}
		propagatedDistance += snext; 
	}
	return rad;   
}

double EUTelGeometryTelescopeGeoDescription::planeRadLengthGlobalIncidence(int planeID, Eigen::Vector3d incidenceDir) {
	
	incidenceDir.normalize();
	double normRad;
	
	TVector3 planeNormalT = siPlaneNormal(planeID);
	Eigen::Vector3d planeNormal(planeNormalT(0), planeNormalT(1), planeNormalT(2));
	
	std::map<int, double>::iterator mapIt = _planeRadMap.find(planeID);
	if( mapIt != _planeRadMap.end() ) {
		normRad = mapIt->second;
	} else {
		Eigen::Vector3d planePosition(siPlaneXPosition(planeID), siPlaneYPosition(planeID), siPlaneZPosition(planeID));

		//We have to propagate halfway to to front and halfway back + a minor safety margin
		Eigen::Vector3d startPoint = planePosition - 0.51*siPlaneZSize(planeID)*planeNormal;
		Eigen::Vector3d endPoint = planePosition + 0.51*siPlaneZSize(planeID)*planeNormal;

		normRad = FindRad(startPoint, endPoint);
		_planeRadMap[planeID] = normRad;
	}
	double scale = std::abs(incidenceDir.dot(planeNormal));
	return normRad/scale;
}

double EUTelGeometryTelescopeGeoDescription::planeRadLengthLocalIncidence(int planeID, Eigen::Vector3d incidenceDir) {
	
	incidenceDir.normalize();
	double normRad;

	std::map<int, double>::iterator mapIt = _planeRadMap.find(planeID);
	if( mapIt != _planeRadMap.end() ) {
		normRad = mapIt->second;
	} else {
		Eigen::Vector3d planePosition(siPlaneXPosition(planeID), siPlaneYPosition(planeID), siPlaneZPosition(planeID));

		TVector3 planeNormalT = siPlaneNormal(planeID);
		Eigen::Vector3d planeNormal(planeNormalT(0), planeNormalT(1), planeNormalT(2));
	
		//We have to propagate halfway to to front and halfway back + a minor safety margin
		Eigen::Vector3d startPoint = planePosition - 0.51*siPlaneZSize(planeID)*planeNormal;
		Eigen::Vector3d endPoint = planePosition + 0.51*siPlaneZSize(planeID)*planeNormal;

		normRad = FindRad(startPoint, endPoint);
		_planeRadMap[planeID] = normRad;
	}
	double scale = std::abs(incidenceDir(2));
	return normRad/scale;
}

bool EUTelGeometryTelescopeGeoDescription::testOutput(std::map< const int,double> & mapSensor,std::map<const int,double> & mapAir){
//    bool foundRadZero = false;
//
//	//TEST ONE SENSOR ONE SCATTERER AFTER.
//	if((mapSensor.size()-1) != mapAir.size()){ //last plane does not contain any air scattering information.
//		throw(std::string("We did not determine the radiation along the track correctly! Sensor != Air"));
//	}
//	//TEST SAME NUMBER AS EXCLUDED SENSORS
//    if( (geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size()-1) != mapAir.size()){
//		throw(std::string("We do not have a scatterer for each plane included "));
//	}
//    streamlog_out(DEBUG5) << "/////////////////////////////////////////////////////////////////////////////////////////////////// " << std::endl;
//    streamlog_out(DEBUG5) << "                 THIS IS WHAT WE WILL CONSTRUCT THE PLANES AND SCATTERING FROM           " << std::endl;
//    for(unsigned int i = 0 ; i<geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() ; i++){
//        streamlog_out(DEBUG5) << "Sensor ID:   "<< geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i)<< " X/X0: " << mapSensor[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i) ] <<" Mass in front of sensor X/X0: "  << mapAir[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i) ] << std::endl;
//        if(mapSensor[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i) ] == 0 or (mapAir[geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().at(i) ] == 0 and i != geo::gGeometry().sensorZOrderToIDWithoutExcludedPlanes().size() -1 )){
//            foundRadZero = true;
//        }
//    }
//    streamlog_out(DEBUG5) << "/////////////////////////////////////////////////////////////////////////////////////////////////// " << std::endl;
//    streamlog_out(DEBUG5) << "/////////////////////////////////////////////////////////////////////////////////////////////////// " << std::endl;
//    if(foundRadZero){
//        return false;
//    }else{
//        return true;
//    }
       return true;

}

//TO DO: Can not exclude thie first plane from pattern recognition. This could also be a problem in general. So we will need to look into pattern recognition.
void EUTelGeometryTelescopeGeoDescription::mapWeightsToSensor(std::map<const int,double> sensor,std::map<const int,double> air,  std::map< const  int, double > & mapSen, std::map< const  int, double > & mapAir ){
//	std::cout<< "sensor size: " << sensor.size() <<std::endl;
//	std::cout<< "air size: " << air.size() <<std::endl;
//
//	unsigned int j=0;
//	double ExcPlaneScat=0;
//	int beforeExcPos=0;
//	bool addMore=false;
//	for(unsigned int  i=0; i<_sensorZOrderToIDMap.size() ; i++){
//		const int sensorID = _sensorZOrderToIDMap[i];
//
//		//Here we check if we have added all the scatterers due to excluded planes. 
//		if(_sensorZOrderToIDWithoutExcludedPlanes[j] == _sensorZOrderToIDMap[i] and  addMore){
//			const int sensorIDBefore = _sensorZOrderToIDWithoutExcludedPlanes[beforeExcPos];
//			//Do not need to add plane scattering again.
//			mapAir[sensorIDBefore] = mapAir[sensorIDBefore]+ExcPlaneScat ;	//Add the same air before and the new plane and air.
//			addMore=false;
//			ExcPlaneScat=0;
//		}
//
//		if(_sensorZOrderToIDWithoutExcludedPlanes[j] == _sensorZOrderToIDMap[i]){//Add the scatterers as normal.
//			const int sensorID = _sensorZOrderToIDWithoutExcludedPlanes[j];
//			mapSen[sensorID] = sensor[sensorID];	
//			if(_sensorZOrderToIDWithoutExcludedPlanes.size() - 1 != j){ //Do not add the last scatterer since not scattering beyond last plane.
//				mapAir[sensorID] = air[sensorID];	
//			}
//			beforeExcPos=j;
//			j++;
//		}else{
//			ExcPlaneScat=ExcPlaneScat + sensor[sensorID] + air[sensorID]; //Add plane and air for excluded plane.
//			addMore=true;
//		}
//	}
}
//
// straight line - shashlyk plane assembler
//
int EUTelGeometryTelescopeGeoDescription::findNextPlane(  double* lpoint,  double* ldir, float* newpoint ) {
	if( newpoint== nullptr) {
		throw(lcio::Exception("You have passed a NULL pointer to findNextPlane.")); 	
	}
	//Here we set the normalised direction and starting point.
	double normdir = TMath::Sqrt(ldir[0]*ldir[0]+ldir[1]*ldir[1]+ldir[2]*ldir[2]); 
	streamlog_out( DEBUG0 ) << "::findNextPlane lpoint: "  << lpoint[0] << " " << lpoint[1] << " "<< lpoint[2] << " " << std::endl;
	ldir[0] = ldir[0]/normdir; 
	ldir[1] = ldir[1]/normdir; 
	ldir[2] = ldir[2]/normdir;
	streamlog_out( DEBUG0 ) << "::findNextPlane ldir  : "  << ldir  [0] << " " << ldir  [1] << " "<< ldir  [2] << " " << std::endl;

	for(int ip=0;ip<3;ip++) {
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
	while( node = gGeoManager->FindNextBoundaryAndStep() ) {
		 inode = node->GetIndex();
		 streamlog_out( DEBUG0 ) << "::findNextPlane found next node: " << node << " id: " << inode << std::endl;
		 const double* point = gGeoManager->GetCurrentPoint();
		 const double* dir   = gGeoManager->GetCurrentDirection();
		 double ipoint[3] ;
		 double idir[3]   ;

		 for(int ip=0;ip<3;ip++) {
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
	if( newpoint == nullptr ) {
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
	while( node = _geoManager->FindNextBoundaryAndStep() ) {
		inode = node->GetIndex();
		const double* point = _geoManager->GetCurrentPoint(); //This will be the new global coordinates after the move
		const double* dir   = _geoManager->GetCurrentDirection(); //This will be the same direction. Since we will only travel in a straight line.  
		double ipoint[3] ;
		double idir[3]   ;
		//Here we set the coordinates and move into the volume in the z direction.
		for(int ip=0;ip<3;ip++) {
			ipoint[ip] = point[ip];
			idir[ip]   = dir[ip];
			if(ip==2) { 
				ipoint[ip]+=0.001 ; // assumption !!! step by one um into the new volume // new volume is thicker than 1 um
			}
			newpoint[ip] = static_cast<float> (ipoint[ip]);
		}
		int sensorID = getSensorID(newpoint); 

		_geoManager->SetCurrentPoint( ipoint);
		_geoManager->SetCurrentDirection( idir);

		streamlog_out( DEBUG0 ) << "Loop number: " << stepNumber  << ". Index of next boundary: " << inode << ". Current global point: " << ipoint[0]  << " " << ipoint[1] << " " << ipoint[2]  << " sensorID: " << sensorID << ". Input of expect next sensor: " << nextSensorID << std::endl;
		streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::findNextPlaneEntrance()------END" << std::endl;

		if( sensorID == nextSensorID ) {
			return true;
		}
		//We return false to say we have not found intersection on the plane.
		if(stepNumber == 10) {
			return false;
		}
		stepNumber++;     
	}
    return false; //If the correct sensor ID is not found by this point then we have failed to find the intersection
}


//This function will intake position and direction. Then using the gear file and magnetic field will output position and sensor ID in correct order of intersection. 
//We need to introduce the idea of:
//sensitive volume => data and state to be created
//scatter volume => state only to be created
//volume => This will cause scattering but no state is to be created on this volume
//At the moment everything in gear file is assumed to be sensitive volume

//std::map<int,double> EUTelGeometryTelescopeGeoDescription::UsingStateReturnAllVolumesIntersected(){}

void EUTelGeometryTelescopeGeoDescription::updateSiPlanesLayout() {
	gear::SiPlanesParameters* siplanesParameters = const_cast< gear::SiPlanesParameters*> (&( _gearManager->getSiPlanesParameters()));
	gear::SiPlanesLayerLayout* siplanesLayerLayout = const_cast< gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));

	// data member::
	_nPlanes = siplanesLayerLayout->getNLayers(); 

	// create an array with the z positions of each layer
	for(size_t iPlane = 0; iPlane < _nPlanes; iPlane++) {
		int sensorID =  _sensorIDVec.at(iPlane);

		siplanesLayerLayout->setLayerPositionX( iPlane, siPlaneXPosition(sensorID) );
		siplanesLayerLayout->setLayerPositionY(  iPlane, siPlaneYPosition(sensorID) );
		siplanesLayerLayout->setLayerPositionZ(  iPlane, siPlaneZPosition(sensorID) );
		siplanesLayerLayout->setLayerRotationZY( iPlane, siPlaneXRotation(sensorID) );
		siplanesLayerLayout->setLayerRotationZX( iPlane, siPlaneYRotation(sensorID) );
		siplanesLayerLayout->setLayerRotationXY( iPlane, siPlaneZRotation(sensorID) );
	}


	// ------- add to GearMgr ----
	if( _gearManager != nullptr ) {
		_gearManager->setSiPlanesParameters( siplanesParameters ) ;
	}
}

void EUTelGeometryTelescopeGeoDescription::updateTrackerPlanesLayout() {
	gear::TrackerPlanesParameters* trackerplanesParameters  = const_cast< gear::TrackerPlanesParameters*>  (&( _gearManager->getTrackerPlanesParameters() ));
	gear::TrackerPlanesLayerLayout* trackerplanesLayerLayout = const_cast< gear::TrackerPlanesLayerLayout*> (&(  trackerplanesParameters->getTrackerPlanesLayerLayout() ));

	trackerplanesParameters->setLayoutID( getSiPlanesLayoutID() );

	// create an array with the z positions of each layer
	int nLayers = trackerplanesLayerLayout->getNLayers();
	for (int iLayer = 0; iLayer < nLayers; iLayer++) {
		gear::TrackerPlanesLayerImpl*  trackerplanesLayerImpl = const_cast< gear::TrackerPlanesLayerImpl*>  ( trackerplanesLayerLayout->getLayer( iLayer) );
		int nsensitive =  trackerplanesLayerImpl->getNSensitiveLayers() ;
		gear::TrackerPlanesSensitiveLayerImplVec& vector =  trackerplanesLayerImpl->getSensitiveLayerVec();

		for (int iSensLayer = 0; iSensLayer < nsensitive; iSensLayer++) {       
			gear::TrackerPlanesSensitiveLayerImpl& sensitiveLayer = vector.at(iSensLayer);
			for( size_t iplane = 0; iplane < _sensorIDVec.size(); iplane++ ) {
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
	if( _gearManager != nullptr ) {
		_gearManager->setTrackerPlanesParameters( trackerplanesParameters );
	}
}

void EUTelGeometryTelescopeGeoDescription::updateGearManager() {
	if( _siPlanesDefined ) {
		updateSiPlanesLayout();
	} else if( _telPlanesDefined ) {
		updateTrackerPlanesLayout();
	}
}

void EUTelGeometryTelescopeGeoDescription::writeGEARFile(std::string filename) {
	updateGearManager();
	gear::GearXML::createXMLFile(marlin::Global::GEAR, filename);
}
