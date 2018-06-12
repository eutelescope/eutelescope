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
#include "EUTelUtility.h"

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
Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getPlaneNormalVector( int planeID )
{
	auto mapIt = _planeNormalMap.find(planeID);
	if( mapIt != _planeNormalMap.end() ) {
		return mapIt->second;
	} else {
		std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
		if( it != _sensorIDVec.end() ) {
			std::array<double,3> const zAxisLocal {{0,0,1}};
			std::array<double,3> zAxisGlobal; 
			local2MasterVec(planeID, zAxisLocal, zAxisGlobal); 
			Eigen::Vector3d normVec(zAxisGlobal.data());
			_planeNormalMap[planeID] = normVec;
			return normVec;
		} else {
			std::stringstream ss;
			ss << planeID;
			std::string errMsg = "EUTelGeometryTelescopeGeoDescription::getPlaneNormalVector: Could not find planeID: " + ss.str();
			throw InvalidGeometryException(errMsg);
		}
	}
}

/**TODO: Replace me: NOP*/
Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getPlaneXVector( int planeID ) {
	auto mapIt = _planeXMap.find(planeID);
	if( mapIt != _planeXMap.end() ) {
		return mapIt->second;
	} else {
		std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
		if( it != _sensorIDVec.end() ) {
			std::array<double,3> const xAxisLocal {{1,0,0}};
			std::array<double,3> xAxisGlobal; 
			local2MasterVec(planeID, xAxisLocal, xAxisGlobal); 
			Eigen::Vector3d xVec(xAxisGlobal.data());
			_planeXMap[planeID] = xVec;
			return xVec;
		} else {
			std::stringstream ss;
			ss << planeID;
			std::string errMsg = "EUTelGeometryTelescopeGeoDescription::getPlaneXVector: Could not find planeID: " + ss.str();
			throw InvalidGeometryException(errMsg);
		}
	}
}

/**TODO: Replace me: NOP*/
Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getPlaneYVector( int planeID ) {
	auto mapIt = _planeYMap.find(planeID);
	if( mapIt != _planeYMap.end() ) {
		return mapIt->second;
	} else {
		std::vector<int>::iterator it = std::find(_sensorIDVec.begin(), _sensorIDVec.end(), planeID);
		if( it != _sensorIDVec.end() ) {
			std::array<double,3> const yAxisLocal {{0,1,0}};
			std::array<double,3> yAxisGlobal; 
			local2MasterVec(planeID, yAxisLocal, yAxisGlobal); 
			Eigen::Vector3d yVec(yAxisGlobal.data());
			_planeYMap[planeID] = yVec;
			return yVec;
		} else {
			std::stringstream ss;
			ss << planeID;
			std::string errMsg = "EUTelGeometryTelescopeGeoDescription::getPlaneYVector: Could not find planeID: " + ss.str(); 
			throw InvalidGeometryException(errMsg);
		}
	}
}

void EUTelGeometryTelescopeGeoDescription::readSiPlanesLayout() {
	// sensor-planes in geometry navigation:
	_siPlanesParameters = const_cast<gear::SiPlanesParameters*> (&( _gearManager->getSiPlanesParameters()));
	_siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));
	auto nPlanes = static_cast<size_t>(_siPlanesLayerLayout->getNLayers());

	//read the geoemtry names from the "Geometry" StringVec section of the gear file
	lcio::StringVec geometryNameParameters;

	try {
		geometryNameParameters  =  _siPlanesParameters->getStringVals("Geometry");
	} catch(gear::UnknownParameterException& e) {
		streamlog_out(MESSAGE6) << "No Geometry field found in GEAR file, assuming CAST for all planes" << std::endl;
		for(size_t i = 0; i < nPlanes; i++) {
			geometryNameParameters.push_back("CAST");
		}
	}

	setLayoutID( _siPlanesParameters->getSiPlanesID() ) ;

	// create an array with the z positions of each layer
	for (size_t iPlane_sz = 0; iPlane_sz < nPlanes; iPlane_sz++) {

    auto iPlane = static_cast<int>(iPlane_sz);
		auto ID = _siPlanesLayerLayout->getID(iPlane);

		auto thisMat = EUTelMaterial(0, 0, 0);
    //The radiation length is in [mm]
		thisMat._radLength = _siPlanesLayerLayout->getSensitiveRadLength(iPlane);
		_materialMap.insert(std::make_pair(std::to_string(ID),thisMat));

		auto thisLayer = std::make_unique<EUTelLayer>(ID);
			
		auto matIt = _materialMap.find(std::to_string(ID));
		auto activePlane = std::make_unique<EUTelActive>(ID, thisLayer.get(), matIt->second);

		auto xP	= _siPlanesLayerLayout->getLayerPositionX(iPlane);
		auto yP	= _siPlanesLayerLayout->getLayerPositionY(iPlane);
		auto zP	= _siPlanesLayerLayout->getLayerPositionZ(iPlane);
		thisLayer->setPosition(xP, yP, zP);

		auto alpha	= _siPlanesLayerLayout->getLayerRotationZY(iPlane);
		auto beta	= _siPlanesLayerLayout->getLayerRotationZX(iPlane);
		auto gamma	= _siPlanesLayerLayout->getLayerRotationXY(iPlane);
		thisLayer->setRotationDeg(alpha, beta, gamma);	

		auto f1	= _siPlanesLayerLayout->getSensitiveRotation1(iPlane);
		auto f2	= _siPlanesLayerLayout->getSensitiveRotation2(iPlane);
		auto f3	= _siPlanesLayerLayout->getSensitiveRotation3(iPlane);
		auto f4	= _siPlanesLayerLayout->getSensitiveRotation4(iPlane);
    
    auto pmOneOrZero = [](double const & f) -> int {
      if( f == 0.0 ) return 0;
      else if( f == 1.0 ) return 1;
      else if( f == -1.0 ) return -1;
      else throw std::runtime_error("EUTelGeometryTelescopeGeoDescription::readSiPlanesLayout flip matrix exntry is not plus/minus one or zero");
    };

		activePlane->setFlips( pmOneOrZero(f1), pmOneOrZero(f2), pmOneOrZero(f3), pmOneOrZero(f4) ); 
	
		activePlane->setGeometry( geometryNameParameters[iPlane]  );

		auto xSize	= _siPlanesLayerLayout->getSensitiveSizeX(iPlane);
		auto ySize	= _siPlanesLayerLayout->getSensitiveSizeY(iPlane);
		auto zSize	= _siPlanesLayerLayout->getSensitiveThickness(iPlane);
		activePlane->setSize(xSize, ySize, zSize);	

		auto xPixelNo	= _siPlanesLayerLayout->getSensitiveNpixelX(iPlane);
		auto yPixelNo	= _siPlanesLayerLayout->getSensitiveNpixelY(iPlane);
		activePlane->setNoPixels(xPixelNo, yPixelNo);
	
		auto xPitch	= _siPlanesLayerLayout->getSensitivePitchX(iPlane);
		auto yPitch	= _siPlanesLayerLayout->getSensitivePitchY(iPlane); 
		activePlane->setPitch(xPitch, yPitch);
		
		auto xRes	= _siPlanesLayerLayout->getSensitiveResolution(iPlane); //should be ResolutionX
		auto yRes	= _siPlanesLayerLayout->getSensitiveResolution(iPlane); //should be ResolutionY
		activePlane->setResolution(xRes, yRes);

		activePlane->setOffset( 0,0,0);
		activePlane->setDeltaRotation( 0,0,0); 

		_activeMap.insert(std::pair<int, EUTelActive*>(ID, activePlane.get()));
		_sensorIDVec.push_back(ID);
		thisLayer->addActivePlane(std::move(activePlane));
		_telescopeLayerMap[ID] = thisLayer.get(); 
		_telescopeLayers.push_back(std::move(thisLayer));
	}

	std::sort(_sensorIDVec.begin(), _sensorIDVec.end(), [&](int a, int b)-> bool {
		return getPlaneZPosition(a) < getPlaneZPosition(b);
	}); 

	for(auto& layer: _telescopeLayers){
		std::cout << "Si Layer: " << layer->getID() << '\n';
		for(auto& active: layer->getActivePlanes()) {
			std::cout << "\t\t -- active: " << active->getID() << '\n';
		}
		for(auto& passive: layer->getPassivePlanes()) {
			std::cout << "\t\t -- passive: " << passive->getID() << '\n';
		}
	}


}

inline double getRadLength(double A, double Z) {
	auto nom = 716.4*A;
	auto denom = Z*(Z+1)*log(287/sqrt(Z));
	return nom/denom;
}

void EUTelGeometryTelescopeGeoDescription::readTrackerPlanesLayout() {
	// sensor-planes in geometry navigation:
	_trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&( _gearManager->getTrackerPlanesParameters()));
	_trackerPlanesLayerLayout = const_cast< gear::TrackerPlanesLayerLayout*> (&(_trackerPlanesParameters->getTrackerPlanesLayerLayout()));

	setLayoutID( _trackerPlanesParameters->getLayoutID() ) ;

	_sensorIDVec.clear();
	
	auto matNameVec = _gearManager->getMaterialNames();
	std::cout << "Known materials: " << std::endl;
	for(auto& matName: matNameVec) {
		std::cout << "-> " << matName << std::endl;
		auto const & GEARMat = _gearManager->getSimpleMaterial(matName);
		auto mat = EUTelMaterial(GEARMat.getA(), GEARMat.getZ(), GEARMat.getDensity());

		if(GEARMat.getRadLength() == 0) {
			mat._radLength = getRadLength(mat._A,mat._Z)/GEARMat.getDensity()*10;
			std::cout << "Radiation length computed to be: " <<  getRadLength(mat._A,mat._Z) << " g*cm^-2" << '\n';
			std::cout << "Yielding: " <<  getRadLength(mat._A,mat._Z)/GEARMat.getDensity()*10 << " mm" << std::endl;
		} else {
			mat._radLength = GEARMat.getRadLength();
		}
		_materialMap.insert(std::make_pair(matName,mat));
	}  

	// create an array with the z positions of each layer
	auto nLayers = _trackerPlanesLayerLayout->getNLayers();
	for (size_t iLayer = 0; iLayer < nLayers; iLayer++) {

		auto _trackerPlanesLayerImpl = const_cast<gear::TrackerPlanesLayerImpl*>(_trackerPlanesLayerLayout->getLayer( static_cast<int>(iLayer)));
		auto layerID = _trackerPlanesLayerImpl->getID();

		if( _telescopeLayerMap.find(layerID) != _telescopeLayerMap.end()){
			streamlog_out(ERROR5)  << "Layer: " << layerID << " already existing, layer IDs must be unique. Terminating" << std::endl;
			throw std::exception();		
		}

		auto thisLayer = std::make_unique<EUTelLayer>(layerID);

		thisLayer->setPosition(  _trackerPlanesLayerImpl->getPositionX(), _trackerPlanesLayerImpl->getPositionY(), _trackerPlanesLayerImpl->getPositionZ() );
		thisLayer->setRotation( _trackerPlanesLayerImpl->getRotationZY(),  _trackerPlanesLayerImpl->getRotationZX(), _trackerPlanesLayerImpl->getRotationXY() );
		thisLayer->setPositionUnc( _trackerPlanesLayerImpl->getPositionXunc(), _trackerPlanesLayerImpl->getPositionYunc(), _trackerPlanesLayerImpl->getPositionZunc() );
		thisLayer->setRotationUnc( _trackerPlanesLayerImpl->getRotationZYunc(), _trackerPlanesLayerImpl->getRotationZXunc(),_trackerPlanesLayerImpl->getRotationXYunc() );
		thisLayer->setInfo( _trackerPlanesLayerImpl->getInfo() );
		thisLayer->setGearPtr(_trackerPlanesLayerImpl);

		auto sensitiveLayerVector = _trackerPlanesLayerImpl->getSensitiveLayerVec();
		for(auto& sensitiveLayer: sensitiveLayerVector) {		

			int activeID = sensitiveLayer.getID();
			auto materialName = sensitiveLayer.getMaterial();

			auto matIt = _materialMap.find(materialName);
			if(matIt == _materialMap.end()) {
				streamlog_out(ERROR5)	<< "Could not find material: " << materialName << " required for <active> with ID: " << activeID 
										<< " on layer: " << layerID << " in GEAR file. Please make sure you incldued it. Terminating!\n";
				throw std::exception(); 
			}
			auto activePlane = std::make_unique<EUTelActive>(activeID, thisLayer.get(), matIt->second);

			activePlane->setOffset( sensitiveLayer.getOffsetX(), sensitiveLayer.getOffsetY(), sensitiveLayer.getOffsetZ() );
			activePlane->setDeltaRotation( sensitiveLayer.getDeltaRotationZY(), sensitiveLayer.getDeltaRotationZX(), sensitiveLayer.getDeltaRotationXY() );
			activePlane->setOffsetUnc( sensitiveLayer.getOffsetXunc(), sensitiveLayer.getOffsetYunc(), sensitiveLayer.getOffsetZunc() );
			activePlane->setDeltaRotationUnc( sensitiveLayer.getDeltaRotationZYunc(), sensitiveLayer.getDeltaRotationZXunc(), sensitiveLayer.getDeltaRotationXYunc() );
			auto pixName = sensitiveLayer.getGeometry();
			activePlane->setGeometry( pixName.empty() ? "CAST" : pixName );
			activePlane->setFlips( sensitiveLayer.getFlip1(), sensitiveLayer.getFlip2(), sensitiveLayer.getFlip3(), sensitiveLayer.getFlip4() );
			activePlane->setInfo( sensitiveLayer.getInfo() );

			auto sizeX = sensitiveLayer.getSizeX();
			auto sizeY = sensitiveLayer.getSizeY();


			auto pitchX  = sensitiveLayer.getPitchX();
			auto pitchY  = sensitiveLayer.getPitchY();
			auto nPixelsX = sensitiveLayer.getNpixelX();
			auto nPixelsY = sensitiveLayer.getNpixelY();

			if( (sizeX != 0 || sizeY != 0 ) && ( pitchX != 0 || pitchY != 0 ) ) {
				streamlog_out(ERROR5) << "You're messing with GEAR! On plane: " << activeID <<'\n' << "values for size as well as pitch, and noPixels are defined." 
										<< " You must either define the number of pixels and their pitch and not use the external pixel libraries or merely size and nothign else!\n";
			}

			if( pitchX == 0 ) {
				activePlane->setSize(sizeX, sizeY, sensitiveLayer.getThickness());	
			} else {
				activePlane->setSize(pitchX*nPixelsX, pitchY*nPixelsY, sensitiveLayer.getThickness() );
			}
			activePlane->setPitch(pitchX, pitchY );
			activePlane->setNoPixels(nPixelsX, nPixelsY );
			activePlane->setResolution(sensitiveLayer.getResolutionX(), sensitiveLayer.getResolutionY() );

			activePlane->setEnabled( sensitiveLayer.getEnabled() );


			auto previousSameIdIt = _activeMap.find(activeID);
			if(previousSameIdIt == _activeMap.end()) {
				_activeMap.insert(std::pair<int, EUTelActive*>(activeID, activePlane.get()));
			} else {
				streamlog_out(ERROR5)	<< "Trying to add <active> with ID: " << activeID << " from layer: " << layerID << ", but same activeID/sensorID already existing on layer: "
										<< (previousSameIdIt->second)->getParent()->getID() << '\n'; 
										throw std::exception();
			}
			_sensorIDVec.push_back(activeID);
			thisLayer->addActivePlane(std::move(activePlane));
		}	
		
		auto materialLayerVector = _trackerPlanesLayerImpl->getMaterialLayerVec();
		for(auto& materialLayer: materialLayerVector) {
			auto passiveID = materialLayer.getID();
			auto materialName = materialLayer.getMaterial();

			auto matIt = _materialMap.find(materialName);
			if(matIt == _materialMap.end()) {
				streamlog_out(ERROR5)	<< "Could not find material: " << materialName << " required for <passive> with ID: " << passiveID 
										<< " on layer: " << layerID << " in GEAR file. Please make sure you incldued it. Terminating!\n";
				throw std::exception();
			} 
			auto passivePlane = std::make_unique<EUTelPassive>(passiveID, thisLayer.get(), matIt->second);

			passivePlane->setOffset( materialLayer.getOffsetX(), materialLayer.getOffsetY(), materialLayer.getOffsetZ());
			passivePlane->setDeltaRotation( materialLayer.getDeltaRotationZY(), materialLayer.getDeltaRotationZX(), materialLayer.getDeltaRotationXY());

			passivePlane->setSize( materialLayer.getSizeX(), materialLayer.getSizeY(), materialLayer.getThickness() );
			passivePlane->setInfo( materialLayer.getInfo() );

			thisLayer->addPassivePlane(std::move(passivePlane));
		}

	_telescopeLayerMap[layerID] = thisLayer.get(); 
	_telescopeLayers.push_back(std::move(thisLayer));
	}

	for(auto& layer: _telescopeLayers){
		std::cout << "Tracker Layer: " << layer->getID() << '\n';
		for(auto& active: layer->getActivePlanes()) {
			std::cout << "\t\t -- active: " << active->getID() << '\n';
		}
		for(auto& passive: layer->getPassivePlanes()) {
			std::cout << "\t\t -- passive: " << passive->getID() << '\n';
		}
	}


	std::sort(_sensorIDVec.begin(), _sensorIDVec.end(), [&](int a, int b)-> bool {
		return getPlaneZPosition(a) < getPlaneZPosition(b);
	}); 

	std::cout << "Sensor IDs ordered by Z: \n";
	for(auto& ID: _sensorIDVec) {
		std::cout << ID << '\n';
	}
	std::cout << std::endl;
	
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
_isGeoInitialized(false),
_geoManager(nullptr)
{
	//Set ROOTs verbosity to only display error messages or higher (so info will not be streamed to stderr)
	gErrorIgnoreLevel =  kError;  
	//Pixel Geometry manager creation
	_pixGeoMgr = std::make_unique<EUTelGenericPixGeoMgr>();
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
}

/**
 * Initialise ROOT geometry objects from external .root file
 * @param tgeofilename name of .root file
 */
void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( std::string const & tgeofilename ) {
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
	int rotRef1, rotRef2, rotRef3, rotRef4;

	std::stringstream strId;
	strId << SensorId;

	// Get sensor center position
	xc = getPlaneXPosition( SensorId );
	yc = getPlaneYPosition( SensorId );
	zc = getPlaneZPosition( SensorId );

	// Get sensor orientation
	alpha = getPlaneXRotationDegrees( SensorId ); //  in degrees !
	beta  = getPlaneYRotationDegrees( SensorId ); // 
	gamma = getPlaneZRotationDegrees( SensorId ); // 

	rotRef1 = planeFlip1( SensorId );
	rotRef2 = planeFlip2( SensorId );
	rotRef3 = planeFlip3( SensorId );
	rotRef4 = planeFlip4( SensorId );

	//We must check that the input is correct. Since this is a combination of initial rotations and reflections the determinate must be 1 or -1
	int determinant = rotRef1*rotRef4 - rotRef2*rotRef3;
	if(determinant == 1 or determinant == -1) { 
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
	std::array<double, 9> integerRotationsAndReflections = { static_cast<double>(rotRef1), static_cast<double>(rotRef2), 0.,
                                                           static_cast<double>(rotRef3), static_cast<double>(rotRef4), 0.,
                                                           0.                   , 0. , static_cast<double>(determinant)   };
	pMatrixRotRefCombined->SetMatrix(integerRotationsAndReflections.data());
	std::cout << "Rotating plane " << SensorId << " to gamma: " << gamma << std::endl;
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
	double radl    = getPlaneRadiationLength( SensorId );
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
	Double_t dx = getPlaneXSize( SensorId ) / 2.;
	Double_t dy = getPlaneYSize( SensorId ) / 2.;
	Double_t dz = getPlaneZSize( SensorId ) / 2.;
	TGeoShape *pBoxSensor = new TGeoBBox( "BoxSensor", dx, dy, dz );

	std::cout << "Box for sensor: " << SensorId << " is: " << dx << "|" << dy  << "|" << dz << '\n';

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
		_pixGeoMgr->addCastedPlane( SensorId, getPlaneNumberOfPixelsX(SensorId), getPlaneNumberOfPixelsY(SensorId), getPlaneXSize(SensorId), getPlaneYSize(SensorId), getPlaneZSize(SensorId), getPlaneRadiationLength(SensorId), stVolName);
	} else {
		_pixGeoMgr->addPlane( SensorId, name, stVolName);
		updatePlaneInfo(SensorId);
	}
}

void EUTelGeometryTelescopeGeoDescription::updatePlaneInfo(int sensorID) {
	
	auto pixGeoDescr = getPixGeoDescr(sensorID );

	double sizeX, sizeY;
  int minX, minY, maxX, maxY;

	pixGeoDescr->getSensitiveSize(sizeX, sizeY);
	pixGeoDescr->getPixelIndexRange(minX, maxX, minY, maxY);
	
	int noX = maxX - minX + 1;
	int noY = maxY - minY + 1;

	std::cout << "Sensor " << sensorID << " has pitch: " << sizeX/noX << "|" << sizeY/noY << std::endl;

	setPlanePitch(sensorID, sizeX/noX, sizeY/noY);
	setPlaneNoPixels(sensorID, noX, noY);

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

   for(auto& mapEntry: _planePath) {
		auto const & pathName = mapEntry.second;
		auto sensorID = mapEntry.first;
    	_geoManager->cd( pathName.c_str() );
		  _TGeoMatrixMap[sensorID] = _geoManager->GetCurrentNode()->GetMatrix();
	  } 
    return;
}

Eigen::Matrix3d EUTelGeometryTelescopeGeoDescription::rotationMatrixFromAngles(int sensorID) {
	return Utility::rotationMatrixFromAngles( static_cast<long double>(getPlaneXRotationRadians(sensorID)), 
                                            static_cast<long double>(getPlaneYRotationRadians(sensorID)), 
                                            static_cast<long double>(getPlaneZRotationRadians(sensorID)));
}

Eigen::Vector3d EUTelGeometryTelescopeGeoDescription::getOffsetVector(int sensorID) {
	Eigen::Vector3d offsetVec;
	offsetVec << getPlaneXPosition(sensorID), getPlaneYPosition(sensorID), getPlaneZPosition(sensorID); 
	return offsetVec;
}

Eigen::Matrix3i EUTelGeometryTelescopeGeoDescription::getFlipMatrix(int sensorID) {
	Eigen::Matrix3i flipMat;
	flipMat << 	planeFlip1(sensorID),	planeFlip2(sensorID),	0,
			        planeFlip3(sensorID), planeFlip4(sensorID),	0,
			        0                   , 0                   , 1;	
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
	_TGeoMatrixMap[sensorID]->LocalToMaster(localPos, globalPos);
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
	_TGeoMatrixMap[sensorID]->MasterToLocal(globalPos, localPos);
}

/**
 * Vector coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param globalVec (x,y,z) in global coordinate system
 * @param localVec (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::local2MasterVec( int sensorID, const double localVec[], double globalVec[] ) {
	_TGeoMatrixMap[sensorID]->LocalToMasterVect(localVec, globalVec);
}

/**
 * Vector coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param globalVec (x,y,z) in global coordinate system
 * @param localVec (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::master2LocalVec( int sensorID, const double globalVec[], double localVec[] ) {
	_TGeoMatrixMap[sensorID]->MasterToLocalVect(globalVec, localVec);
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

double EUTelGeometryTelescopeGeoDescription::getRadiationLengthBetweenPoints(Eigen::Vector3d const & startPt, Eigen::Vector3d const & endPt) {

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

	TGeoMedium* med = nullptr;
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
	
	Eigen::Vector3d planeNormal = getPlaneNormalVector(planeID);
	
	std::map<int, double>::iterator mapIt = _planeRadMap.find(planeID);
	if( mapIt != _planeRadMap.end() ) {
		normRad = mapIt->second;
	} else {
		Eigen::Vector3d planePosition(getPlaneXPosition(planeID), getPlaneYPosition(planeID), getPlaneZPosition(planeID));

		//We have to propagate halfway to to front and halfway back + a minor safety margin
		Eigen::Vector3d startPoint = planePosition - 0.51*getPlaneZSize(planeID)*planeNormal;
		Eigen::Vector3d endPoint = planePosition + 0.51*getPlaneZSize(planeID)*planeNormal;

		normRad = getRadiationLengthBetweenPoints(startPoint, endPoint);
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
		Eigen::Vector3d planePosition(getPlaneXPosition(planeID), getPlaneYPosition(planeID), getPlaneZPosition(planeID));
		Eigen::Vector3d planeNormal = getPlaneNormalVector(planeID);
	
		//We have to propagate halfway to to front and halfway back + a minor safety margin
		Eigen::Vector3d startPoint = planePosition - 0.51*getPlaneZSize(planeID)*planeNormal;
		Eigen::Vector3d endPoint = planePosition + 0.51*getPlaneZSize(planeID)*planeNormal;

		normRad = getRadiationLengthBetweenPoints(startPoint, endPoint);
		_planeRadMap[planeID] = normRad;
	}
	double scale = std::abs(incidenceDir(2));
	return normRad/scale;
}

void EUTelGeometryTelescopeGeoDescription::updateSiPlanesLayout() {
	auto siplanesParameters = const_cast<gear::SiPlanesParameters*> (&( _gearManager->getSiPlanesParameters()));
	auto siplanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));

	// data member::
	auto nPlanes = static_cast<size_t>(siplanesLayerLayout->getNLayers());

	// create an array with the z positions of each layer
	for(size_t iPlane_sz = 0; iPlane_sz < nPlanes; iPlane_sz++) {

    auto iPlane = static_cast<int>(iPlane_sz);
		int sensorID =  _sensorIDVec.at(iPlane);

		std::cout << "Set layer " << sensorID << " gamma to: " <<  getPlaneZRotationDegrees(sensorID) << std::endl;

		siplanesLayerLayout->setLayerPositionX( iPlane, getPlaneXPosition(sensorID) );
		siplanesLayerLayout->setLayerPositionY(  iPlane, getPlaneYPosition(sensorID) );
		siplanesLayerLayout->setLayerPositionZ(  iPlane, getPlaneZPosition(sensorID) );
		siplanesLayerLayout->setLayerRotationZY( iPlane, getPlaneXRotationDegrees(sensorID) );
		siplanesLayerLayout->setLayerRotationZX( iPlane, getPlaneYRotationDegrees(sensorID) );
		siplanesLayerLayout->setLayerRotationXY( iPlane, getPlaneZRotationDegrees(sensorID) );
	}

	// ------- add to GearMgr ----
	if( _gearManager != nullptr ) {
		_gearManager->setSiPlanesParameters( siplanesParameters ) ;
	}
}

void EUTelGeometryTelescopeGeoDescription::updateTrackerPlanesLayout() {

	for(auto& layer: _telescopeLayers) {
		auto GEARLayerPtr = layer->getGearPtr();
		auto layerPos = layer->getPosVec();
		auto layerPosUnc = layer->getPosUncVec();
		auto layerAngle = layer->getAngleVec();

		std::cout << "Writing back position: " << layerPos << std::endl;

		GEARLayerPtr->setPositionX(layerPos.coeff(0));
		GEARLayerPtr->setPositionY(layerPos.coeff(1));
		GEARLayerPtr->setPositionZ(layerPos.coeff(2));

		GEARLayerPtr->setRotationZY(layerAngle.coeff(0)*DEG);
		GEARLayerPtr->setRotationZX(layerAngle.coeff(1)*DEG);
		GEARLayerPtr->setRotationXY(layerAngle.coeff(2)*DEG);

		GEARLayerPtr->setPositionXunc(layerPosUnc.coeff(0));
		GEARLayerPtr->setPositionYunc(layerPosUnc.coeff(1));
		GEARLayerPtr->setPositionZunc(layerPosUnc.coeff(2));

/*  
    virtual void setRotationXYunc( double value)   = 0;
    virtual void setRotationZXunc( double value)   = 0;
    virtual void setRotationZYunc( double value)   = 0;
*/
	}

	// ------- add to GearMgr ----
//	if( _gearManager != nullptr ) {
//		_gearManager->setTrackerPlanesParameters( trackerplanesParameters );
//	}
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
