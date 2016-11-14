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
#include "gearimpl/SimpleMaterialImpl.h"

// EUTELESCOPE
#include "EUTelExceptions.h"
#include "EUTelGenericPixGeoMgr.h"

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

		thisPlane.f1	= _siPlanesLayerLayout->getSensitiveRotation1(iPlane);
		thisPlane.f2	= _siPlanesLayerLayout->getSensitiveRotation2(iPlane);
		thisPlane.f3	= _siPlanesLayerLayout->getSensitiveRotation3(iPlane);
		thisPlane.f4	= _siPlanesLayerLayout->getSensitiveRotation4(iPlane);

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

double getRadLength(int A, int Z) {
	auto nom = 716.4*A;
	auto denom = Z*(Z+1)*log(287/sqrt(Z));
	return nom/denom;
}

void EUTelGeometryTelescopeGeoDescription::readTrackerPlanesLayout() {
	// sensor-planes in geometry navigation:
	_trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&( _gearManager->getTrackerPlanesParameters()));
	_trackerPlanesLayerLayout = const_cast< gear::TrackerPlanesLayerLayout*> (&(_trackerPlanesParameters->getTrackerPlanesLayerLayout()));

	setSiPlanesLayoutID( _trackerPlanesParameters->getLayoutID() ) ;

	_sensorIDVec.clear();
	
	//should be filled based on the length of the sensor vector after the loop
	_nPlanes = 0; 

	auto matNameVec = _gearManager->getMaterialNames();
	std::cout << "Known materials: " << std::endl;
	for(auto& matName: matNameVec) {
		std::cout << "-> " << matName << std::endl;

		//auto& GEARMat = dynamic_cast<gear::SimpleMaterialImpl const&>(_gearManager->getSimpleMaterial(matName));
		auto const& GEARMat = _gearManager->getSimpleMaterial(matName);
		auto A = GEARMat.getA();
		auto Z = GEARMat.getZ();
		auto density = GEARMat.getDensity();
		std::cout << "Radiation length computed to be: " << getRadLength(A,Z) << " g*cm^-2" << '\n';
		std::cout << "Yielding: " << getRadLength(A,Z)/density << " cm" << std::endl;
	}  

	// create an array with the z positions of each layer
	int nLayers = _trackerPlanesLayerLayout->getNLayers();
	for (int iLayer = 0; iLayer < nLayers; iLayer++) {

		auto _trackerPlanesLayerImpl = const_cast<gear::TrackerPlanesLayerImpl*>(_trackerPlanesLayerLayout->getLayer( iLayer));
		auto layerID = _trackerPlanesLayerImpl->getID();

		std::cout << "Got Layer: " << layerID << " .. checking material and sensitive layers ..." << std::endl;
		std::cout << "It has positions: " << _trackerPlanesLayerImpl->getPositionX() << "|" << _trackerPlanesLayerImpl->getPositionY() << "|" << _trackerPlanesLayerImpl->getPositionZ()  << std::endl; 

		int nsensitive = _trackerPlanesLayerImpl->getNSensitiveLayers();
		auto sensitiveLayerVector = _trackerPlanesLayerImpl->getSensitiveLayerVec();
		auto materialLayerVector = _trackerPlanesLayerImpl->getMaterialLayerVec();

		for(auto& sensitiveLayer: sensitiveLayerVector) {
			int sensorID =   sensitiveLayer.getID();

			std::cout << "It has sensitive: " << sensorID << std::endl;

			EUTelPlane thisPlane;

			thisPlane.xPos	= sensitiveLayer.getOffsetX();
			thisPlane.yPos	= sensitiveLayer.getOffsetY();
			thisPlane.zPos	= sensitiveLayer.getOffsetZ();

			thisPlane.xPosUnc	= sensitiveLayer.getOffsetXunc();
			thisPlane.yPosUnc	= sensitiveLayer.getOffsetYunc();
			thisPlane.zPosUnc	= sensitiveLayer.getOffsetZunc();

			thisPlane.alpha	= sensitiveLayer.getDeltaRotationZY();
			thisPlane.beta	= sensitiveLayer.getDeltaRotationZX();
			thisPlane.gamma	= sensitiveLayer.getDeltaRotationXY();

			thisPlane.alphaUnc	= sensitiveLayer.getDeltaRotationZYunc();
			thisPlane.betaUnc	= sensitiveLayer.getDeltaRotationZXunc();
			thisPlane.gammaUnc	= sensitiveLayer.getDeltaRotationXYunc();

			auto pixName = sensitiveLayer.getGeometry();

			thisPlane.pixGeoName = pixName.empty() ? "CAST" : pixName;

			thisPlane.f1	= sensitiveLayer.getFlip1();
			thisPlane.f2	= sensitiveLayer.getFlip2();
			thisPlane.f3	= sensitiveLayer.getFlip3();
			thisPlane.f4	= sensitiveLayer.getFlip4();

			thisPlane.xPixelNo	= sensitiveLayer.getNpixelX();
			thisPlane.yPixelNo	= sensitiveLayer.getNpixelY();
			thisPlane.xPitch	= sensitiveLayer.getPitchX();
			thisPlane.yPitch	= sensitiveLayer.getPitchY();
			thisPlane.xRes	= sensitiveLayer.getResolutionX();
			thisPlane.yRes	= sensitiveLayer.getResolutionY();

			thisPlane.xSize	= thisPlane.xPixelNo*thisPlane.xPitch; 
			thisPlane.ySize	= thisPlane.yPixelNo*thisPlane.yPitch;
			thisPlane.zSize	= sensitiveLayer.getThickness();

			thisPlane.enabled = sensitiveLayer.getEnabled();
			
			//GEAR uses mm wheras TGeo will use cm
			thisPlane.radLength	= 10;//sensitiveLayer.getRadLength()/10;

			_planeSetup[sensorID] = thisPlane;

			_sensorIDVec.push_back(sensorID);

			std::cout << "Pos: " << thisPlane.xPos << "|" << thisPlane.yPos << "|" << thisPlane.zPos << '\n';
			std::cout << "PosUnc: " << thisPlane.xPosUnc << "|" << thisPlane.yPosUnc << "|" << thisPlane.zPosUnc << '\n';

			std::cout << "Rot: " << thisPlane.alpha << "|" << thisPlane.beta << "|" << thisPlane.gamma << '\n';
			std::cout << "RotUnc: " << thisPlane.alphaUnc << "|" << thisPlane.betaUnc << "|" << thisPlane.gammaUnc << '\n';
			
			std::cout << "Flip: " << thisPlane.f1 << "|" << thisPlane.f2 << "|" << thisPlane.f3 << "|" << thisPlane.f4 << '\n';

		}
		for(auto& materialLayer: materialLayerVector) {
			auto sensorID = materialLayer.getID();
			std::cout << "It has material: " << sensorID << std::endl;
		}
	}
	std::sort(_sensorIDVec.begin(), _sensorIDVec.end(), doCompare(*this) );
	_nPlanes =  _sensorIDVec.size(); 

	writeGEARFile("test.xml");
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
	_geoManager.release();
	delete _pixGeoMgr;
	_pixGeoMgr = nullptr;
}

/**
 * Initialise ROOT geometry objects from external .root file
 * @param tgeofilename name of .root file
 */
void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( std::string  tgeofilename ) {
    
    _geoManager = std::unique_ptr<TGeoManager>(TGeoManager::Import(tgeofilename.c_str()));
    _geoManager->SetBit(kCanDelete);
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
	pMatrixRotRefCombined->RotateX(alpha);//X Rotations (degrees)//This will rotate a vector usign the right hand rule round the x-axis
	pMatrixRotRefCombined->RotateY(beta);//Y Rotations (degrees)//Same again for Y axis
	pMatrixRotRefCombined->RotateZ(gamma);//Z Rotation (degrees)//This will again rotate a vector around z axis usign the right hand rule.  
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
void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( std::string const & geomName, bool dumpRoot = false ) {
	if( _isGeoInitialized ) {
		streamlog_out( WARNING3 ) << "EUTelGeometryTelescopeGeoDescription: Geometry already initialized, using old initialization" << std::endl;
		return;
	} else {
    		_geoManager = std::make_unique<TGeoManager>("Telescope", "v0.1");
			_geoManager->SetBit(kCanDelete);
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
	rotMat <<	(double)(cosB*cosG),	(double)(sinA*sinB*cosG - sinG*cosA),	(double)(sinA*sinG + sinB*cosA*cosG),
	      		(double)(sinG*cosB),			(double)(sinA*sinB*sinG + cosA*cosG),			(double)(-sinA*cosG + sinB*sinG*cosA),
			(double)(-sinB),	(double)(sinA*cosB),	(double)(cosA*cosB);//for rotation order X->Y->Z
	/*rotMat <<	(double)(cosB*cosG+sinA*sinB*sinG),	(double)(sinA*sinB*cosG-cosB*sinG),	(double)(cosA*sinB),
	      		(double)(cosA*sinG),			(double)(cosA*cosG),			(double)(-sinA),
			(double)(sinA*cosB*sinG-sinB*cosG),	(double)(sinA*cosB*cosG+sinB*sinG),	(double)(cosA*cosB);*/ //for rotation order Z->X->Y
	//std::cout << rotMat.format(IO) << std::endl;
	return rotMat;
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getRotationAnglesFromMatrix(Eigen::Matrix3d rotMat) {
	//for rotation order X->Y->Z
	long double beta = asin((long double)(-rotMat(2,0)));
	long double cosB = cos(beta);

	long double alpha = asin((long double)(rotMat(2,1)/cosB));
	long double gamma = asin((long double)(rotMat(1,0)/cosB));
	
	//for rotation order Z->X->Y 
	/*long double alpha = asin((long double)(-rotMat(1,2)));
	long double cosA = cos(alpha);

	long double beta = asin((long double)(rotMat(0,2)/cosA));
	long double gamma = asin((long double)(rotMat(1,0)/cosA));*/

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
		gear::TrackerPlanesSensitiveLayerImplVec & vector =  const_cast< gear::TrackerPlanesSensitiveLayerImplVec &>(trackerplanesLayerImpl->getSensitiveLayerVec());

		for (int iSensLayer = 0; iSensLayer < nsensitive; iSensLayer++) {       
			gear::TrackerPlanesSensitiveLayerImpl& sensitiveLayer = vector.at(iSensLayer);
			for( size_t iplane = 0; iplane < _sensorIDVec.size(); iplane++ ) {
				int sensorID =  _sensorIDVec.at(iplane);
				if( sensitiveLayer.getID() !=  _sensorIDVec.at( iplane) ) continue;  
				sensitiveLayer.setOffsetX( siPlaneXPosition(sensorID) );
				sensitiveLayer.setOffsetY( siPlaneYPosition(sensorID) );
				sensitiveLayer.setOffsetZ( siPlaneZPosition(sensorID) );

				sensitiveLayer.setDeltaRotationZY( siPlaneXRotation(sensorID) );
				sensitiveLayer.setDeltaRotationZX( siPlaneYRotation(sensorID) );
				sensitiveLayer.setDeltaRotationXY( siPlaneZRotation(sensorID) );
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
