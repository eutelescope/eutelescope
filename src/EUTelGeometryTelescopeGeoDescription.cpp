/* 
 * File:   EUTelGeometryTelescopeGeoDescription.cpp
 * 
 */

#include "EUTelGeometryTelescopeGeoDescription.h"

// C++
#include <algorithm>
#include <string>

// MARLIN
#include "marlin/Global.h"
#include "marlin/VerbosityLevels.h"

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

#include <cstring>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;
using namespace geo;
using namespace std;


EUTelGeometryTelescopeGeoDescription& EUTelGeometryTelescopeGeoDescription::getInstance() {
      static  EUTelGeometryTelescopeGeoDescription instance;
      return instance;
}

size_t EUTelGeometryTelescopeGeoDescription::nPlanes( ) const {
    return _nPlanes;
}

const EVENT::DoubleVec& EUTelGeometryTelescopeGeoDescription::siPlanesZPositions( ) const {
    return _siPlaneZPosition;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneXPosition( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneXPosition[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneYPosition( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneYPosition[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneZPosition( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneZPosition[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneXRotation( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneXRotation[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneYRotation( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneYRotation[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneZRotation( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneZRotation[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneXSize( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneSizeX[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneYSize( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneSizeY[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneZSize( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneSizeZ[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneMediumRadLen( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneRadLength[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

std::string EUTelGeometryTelescopeGeoDescription::geoLibName( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _geoLibName[ _sensorIDtoZOrderMap[ planeID ] ];
    return "failed";
}

TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneNormal( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) {
        TVector3 normVec( 0., 0., 1. );
        normVec.RotateX( siPlaneXRotation( planeID) ); // to be in rad
        normVec.RotateY( siPlaneYRotation( planeID) ); // to be in rad
        normVec.RotateZ( siPlaneZRotation( planeID) ); // to be in rad
        return normVec;
    }
    return TVector3(0.,0.,0.);
}

const std::map<int, int>& EUTelGeometryTelescopeGeoDescription::sensorIDstoZOrder( ) const {
    return _sensorIDtoZOrderMap;
}

int EUTelGeometryTelescopeGeoDescription::sensorIDtoZOrder( int planeID ) const {
    std::map<int,int>::const_iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return it->second;
    return -1;
}

/** Sensor ID vector ordered according to their position along the Z axis (beam axis)
 *  Numeration runs from 0 to nPlanes-1 */
int EUTelGeometryTelescopeGeoDescription::sensorZOrderToID( int znumber ) const {
    std::map<int,int>::const_iterator it;
    it = _sensorZOrderToIDMap.find( znumber );
    if ( it != _sensorZOrderToIDMap.end() ) return it->second;
    return -1;
}

 
            /** Map from sensor ID to number along Z */
const std::map<int, int>& EUTelGeometryTelescopeGeoDescription::sensorZOrdertoIDs() const {

return _sensorZOrderToIDMap;
}
            

const EVENT::IntVec& EUTelGeometryTelescopeGeoDescription::sensorIDsVec( ) const {
    return _sensorIDVec;
}

EUTelGeometryTelescopeGeoDescription::EUTelGeometryTelescopeGeoDescription() :
_siPlanesParameters(0),
_siPlanesLayerLayout(0),
_sensorIDVec(),
_sensorIDVecMap(),
_sensorZOrderToIDMap(),
_sensorIDtoZOrderMap(),
_siPlaneXPosition(),
_siPlaneYPosition(),
_siPlaneZPosition(),
_siPlaneXRotation(),
_siPlaneYRotation(),
_siPlaneZRotation(),
_nPlanes(0),
_isGeoInitialized(false),
_geoManager(0)
{
    // Check if the GEAR manager is not corrupted, otherwise stop

    if (!marlin::Global::GEAR) {
        streamlog_out(ERROR2) << "The GearMgr is not available, for an unknown reason." << std::endl;
        throw eutelescope::InvalidGeometryException("GEAR manager is not initialised");
    }

    // sensor-planes in geometry navigation:
    _siPlanesParameters = const_cast<gear::SiPlanesParameters*> (&(marlin::Global::GEAR->getSiPlanesParameters()));
    _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));
    
    //read the geoemtry names from the "Geometry" StringVec section of the gear file
    lcio::StringVec geometryNameParameters =  _siPlanesParameters->getStringVals("Geometry");
    
    // create an array with the z positions of each layer
    for (int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++) {
        _siPlaneXPosition.push_back(_siPlanesLayerLayout->getLayerPositionX(iPlane));
        _siPlaneYPosition.push_back(_siPlanesLayerLayout->getLayerPositionY(iPlane));
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getLayerPositionZ(iPlane));
        _siPlaneXRotation.push_back(_siPlanesLayerLayout->getLayerRotationZY(iPlane));
        _siPlaneYRotation.push_back(_siPlanesLayerLayout->getLayerRotationZX(iPlane));
        _siPlaneZRotation.push_back(_siPlanesLayerLayout->getLayerRotationXY(iPlane));
        
        _siPlaneSizeX.push_back(_siPlanesLayerLayout->getLayerSizeX(iPlane));
        _siPlaneSizeY.push_back(_siPlanesLayerLayout->getLayerSizeY(iPlane));
        _siPlaneSizeZ.push_back(_siPlanesLayerLayout->getLayerThickness(iPlane));
        
        _siPlaneRadLength.push_back(_siPlanesLayerLayout->getLayerRadLength(iPlane));
	_geoLibName.push_back(geometryNameParameters[iPlane]);
    }


    if (_siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT) {
        _siPlaneXPosition.push_back(_siPlanesLayerLayout->getDUTPositionX());
        _siPlaneYPosition.push_back(_siPlanesLayerLayout->getDUTPositionY());
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getDUTPositionZ());
        // WARNING No DUT rotations in GEAR!!!!!!!!!!
        // TODO: Need this in GEAR
        _siPlaneXRotation.push_back(0.);
        _siPlaneYRotation.push_back(0.);
        _siPlaneZRotation.push_back(0.);
    }

    // sort the array with increasing z
    std::sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());

    // clear the sensor ID vector
    _sensorIDVec.clear();

    // clear the sensor ID map
    _sensorIDVecMap.clear();
    _sensorIDtoZOrderMap.clear();

    double* keepZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];

    for (int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++) {
        int sensorID = _siPlanesLayerLayout->getID(iPlane);
        _sensorIDVec.push_back(sensorID);
        _sensorIDVecMap.insert(std::make_pair(sensorID, iPlane));

        // count number of the sensors to the left of the current one:
        int sensorsToTheLeft = 0;
        keepZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
        for (int jPlane = 0; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++)
            if (_siPlanesLayerLayout->getLayerPositionZ(jPlane) + 1e-06 < keepZPosition[ iPlane ]) sensorsToTheLeft++;

        _sensorZOrderToIDMap.insert(std::make_pair(sensorsToTheLeft, sensorID));        
        _sensorIDtoZOrderMap.insert(std::make_pair(sensorID, sensorsToTheLeft));
    }

    delete [] keepZPosition;

    _nPlanes = _siPlanesParameters->getSiPlanesNumber();
    if (_siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT) ++_nPlanes;
    
    
    
    // TGeo manager initialisation
    
	//Pixel Geometry manager creation
	_pixGeoMgr = new EUTelGenericPixGeoMgr();
}

EUTelGeometryTelescopeGeoDescription::~EUTelGeometryTelescopeGeoDescription() {
	delete _geoManager;
	delete _pixGeoMgr;
}

/**
 * Initialise ROOT geometry objects from external .root file
 * @param tgeofilename name of .root file
 */
void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( string tgeofilename ) {
//    #ifdef USE_TGEO
    // get access to ROOT's geometry manager
    
    _geoManager = TGeoManager::Import( tgeofilename.c_str() );
    if( !_geoManager ) {
        streamlog_out( WARNING ) << "Can't read file " << tgeofilename << endl;
    }

    _geoManager->CloseGeometry();
//    #endif //USE_TGEO
}

/**
 * Initialise ROOT geometry objects from GEAR objects
 * 
 * @param geomName name of ROOT geometry object
 * @param dumpRoot dump automatically generated ROOT geometry file for further inspection
 */
void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( std::string& geomName, bool dumpRoot = false ) {
//    #ifdef USE_TGEO
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
   
    
    // Create top world volume containing telescope/DUT geometry
    
    
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
   
 
   
   // Iterate over registered GEAR objects and construct their TGeo representation
   const Double_t PI = 3.141592653589793;
   const Double_t DEG = 180./PI; 
  
   double xc, yc, zc;   // volume center position 
   double alpha, beta, gamma;
   
   IntVec::const_iterator itrPlaneId;
   for ( itrPlaneId = _sensorIDVec.begin(); itrPlaneId != _sensorIDVec.end(); ++itrPlaneId ) {
       
       std::stringstream strId;
       strId << *itrPlaneId;
       
       // Get sensor center position
       xc = siPlaneXPosition( *itrPlaneId );
       yc = siPlaneYPosition( *itrPlaneId );
       zc = siPlaneZPosition( *itrPlaneId );
       
       // Get sensor orientation
       alpha = siPlaneXRotation( *itrPlaneId )*DEG; // [rad] in degrees !
       beta  = siPlaneYRotation( *itrPlaneId )*DEG; // [rad]
       gamma = siPlaneZRotation( *itrPlaneId )*DEG; // [rad]
       
       // Spatial translations of the sensor center
       string stTranslationName = "matrixTranslationSensor";
       stTranslationName.append( strId.str() );
       TGeoTranslation* pMatrixTrans = new TGeoTranslation( stTranslationName.c_str(), xc, yc, zc );
       //ALL clsses deriving from TGeoMatrix are not owned by the ROOT geometry manager, invoking RegisterYourself() transfers
       //ownership and thus ROOT will clean up
       pMatrixTrans->RegisterYourself();      
       
       // Spatial rotation around sensor center
       // TGeoRotation requires Euler angles in degrees
       string stRotationName = "matrixRotationSensorX";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotX = new TGeoRotation( stRotationName.c_str(), 0.,  alpha, 0.);                // around X axis
       stRotationName = "matrixRotationSensorY";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotY = new TGeoRotation( stRotationName.c_str(), 90., beta,  0.);                // around Y axis (combination of rotation around Z axis and new X axis)
       stRotationName = "matrixRotationSensorBackY";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotY1 = new TGeoRotation( stRotationName.c_str(), -90., 0.,  0.);                    // restoration of original orientation (valid in small angle approximataion ~< 15 deg)
       stRotationName = "matrixRotationSensorZ";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotZ = new TGeoRotation( stRotationName.c_str(), 0. , 0.,        gamma);         // around Z axis
       
       // Combined rotation in several steps
       TGeoRotation* pMatrixRot = new TGeoRotation( *pMatrixRotX );
       pMatrixRot->MultiplyBy( pMatrixRotY );
       pMatrixRot->MultiplyBy( pMatrixRotY1 );
       pMatrixRot->MultiplyBy( pMatrixRotZ );
       pMatrixRot->RegisterYourself();      
      
       pMatrixRotX->RegisterYourself();
       pMatrixRotY->RegisterYourself();
       pMatrixRotY1->RegisterYourself(); 
       pMatrixRotZ->RegisterYourself();
 
       // Combined translation and orientation
       TGeoCombiTrans* combi = new TGeoCombiTrans( *pMatrixTrans, *pMatrixRot );
       combi->RegisterYourself();   
 
       // Construction of sensor objects
       
       // Construct object medium. Required for radiation length determination

       // assume SILICON, though all information except of radiation length is ignored
       double a       = 28.085500;     
       double z       = 14.000000;
       double density = 2.330000;
       double radl    = siPlaneMediumRadLen( *itrPlaneId );
       double absl    = 45.753206;
       string stMatName = "materialSensor";
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
       string stMedName = "mediumSensor";
       stMedName.append( strId.str() );
       TGeoMedium* pMed = new TGeoMedium( stMedName.c_str(), numed, pMat, par );
       
       // Construct object shape
       // Shape: Box type: TGeoBBox
       // TGeo requires half-width of box side
       dx = siPlaneXSize( *itrPlaneId ) / 2.;
       dy = siPlaneYSize( *itrPlaneId ) / 2.;
       dz = siPlaneZSize( *itrPlaneId ) / 2.;
       TGeoShape *pBoxSensor = new TGeoBBox( "BoxSensor", dx, dy, dz );
       // Volume: volume_Sensor1
       
       // Geometry navigation package requires following names for objects that have an ID
       // name:ID
       string stVolName = "volume_SensorID:";
       stVolName.append( strId.str() );

       _planePath.insert( std::make_pair(*itrPlaneId, "/volume_World_1/"+stVolName+"_1") );

       TGeoVolume* pvolumeSensor = new TGeoVolume( stVolName.c_str(), pBoxSensor, pMed );
       pvolumeSensor->SetVisLeaves( kTRUE );
       pvolumeWorld->AddNode(pvolumeSensor, 1/*(*itrPlaneId)*/, combi);
	
	//this line tells the pixel geometry manager to load the pixel geometry into the plane			
        _pixGeoMgr->addPlane( *itrPlaneId, geoLibName( *itrPlaneId), stVolName);
   } // loop over sensorID

    _geoManager->CloseGeometry();
    _isGeoInitialized = true;
    // Dump ROOT TGeo object into file
    if ( dumpRoot ) _geoManager->Export( geomName.c_str() );

//    #endif //USE_TGEO
    return;
}

/** Determine id of the sensor in which point is locate
 * 
 * @param globalPos 3D point in global reference frame
 * @return sensorID or -999 if the point in outside of sensor volume
 */
int EUTelGeometryTelescopeGeoDescription::getSensorID( const float globalPos[] ) const {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getSensorID() " << std::endl;
    
    _geoManager->FindNode( globalPos[0], globalPos[1], globalPos[2] );

    const char* volName = const_cast < char* > ( geo::gGeometry( )._geoManager->GetCurrentVolume( )->GetName( ) );

    streamlog_out( DEBUG0 ) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") found in volume: " << volName << std::endl;
    
	/*std::vector< std::string > tokens = Utility::stringSplit(std::string( volName ), ":", false );
    // sensor id must be stored in the last token
    int id = -999;
    std::size_t found = tokens.back().find_first_of("0123456789");
    if ( found != std::string::npos ) id = atoi( tokens.back().c_str() );
    int sensorID = -999;
    if ( std::find( _sensorIDVec.begin(), _sensorIDVec.end(), id ) != _sensorIDVec.end() ) sensorID = id;
    else streamlog_out(DEBUG3) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was not found inside any sensor!" << std::endl; 
    streamlog_out( DEBUG0 ) << "SensorID: " << sensorID << std::endl;
    return sensorID;*/

	std::vector<std::string> split = Utility::stringSplit( std::string( volName ), "/", false);
	//std::string sensor =  split.at(1);

	int sensorID = -999;
	streamlog_out(DEBUG3) << "init sensorID  : " << sensorID  <<  " " << volName << std::endl;

        if( split.size() == 1 && split[0].length() > 16 ) {

          streamlog_out(DEBUG3) << "split[0] " << split[0] << std::endl;
          streamlog_out(DEBUG3) << "split[0].substr(0,16) " << split[0].substr(0,16) << std::endl;
          int strLength = split[0].length(); 
          streamlog_out(DEBUG3) << "split[0].substr(16, strLength ) " << split[0].substr(16, strLength ) << std::endl;

          //since we check bounds, no need for vector.at() but use [], it saves cycles :-)
	  if (  (split[0].substr(0,16) == "volume_SensorID:") )
	  {
                sensorID = strtol( (split[0].substr(16, strLength )).c_str(), NULL, 10 );
		streamlog_out(DEBUG3) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was found at :" << sensorID << std::endl;
          }
	  else
	  {
		streamlog_out(DEBUG3) << "Point (" << globalPos[0] << "," << globalPos[1] << "," << globalPos[2] << ") was not found inside any sensor!" << std::endl;
		//Maybe exception??
	  }
        }
        streamlog_out(DEBUG3) << "sensorID  : " << sensorID  << std::endl;

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
 * @param globalPos (x,y,z) in global coordinate system
 * @param localPos (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::master2Local( const double globalPos[], double localPos[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::master2Local() " << std::endl;

		_geoManager->GetMother(1);
		_geoManager->GetCurrentNode()->MasterToLocal( globalPos, localPos );
    
    streamlog_out(DEBUG0) << std::fixed;
    streamlog_out(DEBUG0) << "Global coordinates:" << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << globalPos[0] << std::setw(10) << std::setprecision(5) << globalPos[1] << std::setw(10) << std::setprecision(5) << globalPos[2] << std::endl;
    streamlog_out(DEBUG0) << "Local coordinates: " << std::endl;
    streamlog_out(DEBUG0) << std::setw(10) << std::setprecision(5) << localPos[0] << std::setw(10) << std::setprecision(5) << localPos[1] << std::setw(10) << std::setprecision(5) << localPos[2] << std::endl;
}


void EUTelGeometryTelescopeGeoDescription::local2masterHit(EVENT::TrackerHit* hit_input, IMPL::TrackerHitImpl* hit_output, LCCollection * hitCollectionOutput){
    streamlog_out(DEBUG2) << "START------------------EUTelGeometryTelescopeGeoDescription::local2MasterHit()-------------------------------------- " << std::endl;
		//Get input sensor ID and properties
		UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
		int sensorID = hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_input))["sensorID"];
		int properties = hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_input))["properties"];
		//Now get ist local position
		const double * localPos =  hit_input->getPosition();
		double globalPos[3];
		//Now determine the position in global coordinates.
		local2Master(sensorID, localPos, globalPos);
		//Fill the new hit_output with information
		hit_output->setPosition(globalPos);
		hit_output->setCovMatrix( hit_input->getCovMatrix());
	  	hit_output->setType( hit_input->getType() );
		UTIL::CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING, hitCollectionOutput);
  		idHitEncoder["sensorID"] =  sensorID;
		//This warns the user if global flag has not been set
		if(properties != kHitInGlobalCoord){
			streamlog_out(WARNING5) << " The properties cell ID is not global as expected!  " << std::endl;
		}
		else{
			streamlog_out(WARNING5) << "The properties cell ID is global. Are you sure this hit is in local coordinates?" << std::endl;
		}
			
		idHitEncoder["properties"] = kHitInGlobalCoord;

  	// This is part were we store the encoded information on the hit
  	idHitEncoder.setCellID( hit_output );

    streamlog_out(DEBUG2) << "END------------------EUTelGeometryTelescopeGeoDescription::local2MasterHit()-------------------------------------------------------- " << std::endl;
}

void EUTelGeometryTelescopeGeoDescription::master2localHit(EVENT::TrackerHit* hit_input, IMPL::TrackerHitImpl* hit_output, LCCollection * hitCollectionOutput){
    streamlog_out(DEBUG2) << "START------------------EUTelGeometryTelescopeGeoDescription::master2localHit()-------------------------------------- " << std::endl;
		//Get information about the input hit
		UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
		int sensorID = hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_input))["sensorID"];
		int properties = hitDecoder(static_cast< IMPL::TrackerHitImpl* >(hit_input))["properties"];
		const double * globalPos =  hit_input->getPosition();
		double localPos[3];
		//Change the coordinate from global to local
		master2Local(globalPos, localPos);
		//Fill information on the new local hit
		hit_output->setPosition(localPos);
		hit_output->setCovMatrix( hit_input->getCovMatrix());
	  	hit_output->setType( hit_input->getType() );
		UTIL::CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING, hitCollectionOutput);

	  	idHitEncoder["sensorID"] =  sensorID;
		///Warn the user if the global coordinate flag has been set correctly
		if(properties == kHitInGlobalCoord){
			streamlog_out(WARNING5) << " The properties cell ID is global as expected!  " << std::endl;
		}
		else{
			streamlog_out(WARNING5) << "The properties cell ID is not global. Are you sure this hit is in local coordinates?" << std::endl;
		}
		
		idHitEncoder["properties"] = 0;

  	// This is part were we store the encoded information on the hit
  	idHitEncoder.setCellID( hit_output );

    streamlog_out(DEBUG2) << "END------------------EUTelGeometryTelescopeGeoDescription::master2localHit()-------------------------------------------------------- " << std::endl;
}

/**
 * Vector coordinate transformation from global reference frame to local reference frame.
 * Corresponding volume is determined automatically.
 * 
 * @param globalVec (x,y,z) in global coordinate system
 * @param localVec (x,y,z) in local coordinate system
 */
void EUTelGeometryTelescopeGeoDescription::local2MasterVec( int sensorID, const double localVec[], double globalVec[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::master2LocalVec() " << std::endl;
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
 * Global-to-local coordinate transformation matrix.
 * Corresponding volume is determined automatically.
 * 
 * @param globalPos (x,y,z) in global coordinate system
 * @return 
 */
const TGeoHMatrix* EUTelGeometryTelescopeGeoDescription::getHMatrix( const double globalPos[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getHMatrix() " << std::endl;
    _geoManager->FindNode( globalPos[0], globalPos[1], globalPos[2] );    
    const TGeoHMatrix* globalH = _geoManager->GetCurrentMatrix();
    return globalH;
}

/**
 * Retrieve magnetic field object.
 * 
 * @return reference to gear::BField object
 */
const gear::BField& EUTelGeometryTelescopeGeoDescription::getMagneticFiled() const {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getMagneticFiled() " << std::endl;
    return marlin::Global::GEAR->getBField();
}



/**
 * @return EUTelGenericPixGeoDescr const * for given plane
 * ID, essential for user interfacing pixel data!
 */
EUTelGenericPixGeoDescr* EUTelGeometryTelescopeGeoDescription::getPixGeoDescr( int planeID ){
    return _pixGeoMgr->getPixGeoDescr(planeID);
}

/**
 * @return path of the plane as a std::string
 */
std::string EUTelGeometryTelescopeGeoDescription::getPlanePath( int planeID ){
    return _planePath.find(planeID)->second;
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
    
    float rad = 0.;        // integral of radiation length in units of X0
    
    const double mm2cm = 0.1;
    
    /* TGeo uses cm and grams as internal units e.g. in radiation length and density. Telescope/LCIO uses mm. Therefore this routine is full of 
     annoying conversion factors */    
    
    const double stepLenght2 = ( globalPosFinish[0] - globalPosStart[0] )*( globalPosFinish[0] - globalPosStart[0] ) +
                               ( globalPosFinish[1] - globalPosStart[1] )*( globalPosFinish[1] - globalPosStart[1] ) +
                               ( globalPosFinish[2] - globalPosStart[2] )*( globalPosFinish[2] - globalPosStart[2] );
    
    const double stepLenght  = TMath::Sqrt( stepLenght2 );

    // don't need conversion factor to for calculation of directions
    const double xp  = ( globalPosFinish[0] - globalPosStart[0] )/stepLenght;
    const double yp  = ( globalPosFinish[1] - globalPosStart[1] )/stepLenght;
    const double zp  = ( globalPosFinish[2] - globalPosStart[2] )/stepLenght;

    streamlog_out(DEBUG0) << "Start point (x,y,z):" << globalPosStart[0] << "," << globalPosStart[1] << "," << globalPosStart[2] << std::endl;
    streamlog_out(DEBUG0) << "Finish point (x,y,z):" << globalPosFinish[0] << "," << globalPosFinish[1] << "," << globalPosFinish[2] << std::endl;
    streamlog_out(DEBUG0) << "Direction (nx,ny,nz):" << xp << "," << yp << "," << zp << std::endl;
    
    double snext;
    double pt[3], loc[3];
    double epsil = 1.E-7;
    double lastrad = 0.;
    int ismall       = 0;
    int nbound       = 0;
    float length     = 0.;
    TGeoMedium *med;
    TGeoShape *shape;
    
    // Get starting node
    gGeoManager->InitTrack( globalPosStart[0]/*mm*/, globalPosStart[1]/*mm*/, globalPosStart[2]/*mm*/, xp, yp, zp );
    TGeoNode *nextnode = gGeoManager->GetCurrentNode( );
    
    double currentStep = stepLenght /*mm*/;
    // Loop over all, encountered during the propagation, volumes 
    while ( nextnode ) {
        med = NULL;
        
	// Check if current point is inside silicon sensor. Radiation length of silicon sensors is accounted in thin scatterers of GBL.
        bool isBoundaryVolume = false;
        if ( gGeoManager->IsSameLocation( globalPosStart[0], globalPosStart[1], globalPosStart[2] ) ||
             gGeoManager->IsSameLocation( globalPosFinish[0], globalPosFinish[1], globalPosFinish[2] ) ) isBoundaryVolume = true;
        
        if ( nextnode ) med = nextnode->GetVolume()->GetMedium();
        else return 0.;
        
        shape = nextnode->GetVolume()->GetShape();
        
        // make a step to the next intersection point
        if ( currentStep > 1.e-9 /*mm*/ ) nextnode = gGeoManager->FindNextBoundaryAndStep( currentStep /*mm*/ );
        else return rad;
        
        snext  = gGeoManager->GetStep() /*mm*/;
        
        // Small steps treatment
        if ( snext < 1.e-8 /*mm*/ ) {
            ismall++;
            
            // Terminate calculation if too many small steps done
            if ( ismall > 3 ) {
                streamlog_out( WARNING1 ) << "ERROR: Small steps in: " << gGeoManager->GetPath() << " shape=" << shape->ClassName() << endl;
                return rad;
            }

            // increase step size (epsilon) and advance along the particle direction
            memcpy( pt, gGeoManager->GetCurrentPoint(), 3 * sizeof (double) );
            const double *dir = gGeoManager->GetCurrentDirection();
            for ( Int_t i = 0; i < 3; i++ ) pt[i] += epsil * dir[i];
            snext = epsil;
            length += snext;
            
            // Ignore start and finish volumes if required
            if ( skipBoundaryPonitsVolumes && isBoundaryVolume ) {
                rad += 0.;
            } else {
                rad += lastrad*snext;
            }
            
            gGeoManager->CdTop( );
            nextnode = gGeoManager->FindNode( pt[0], pt[1], pt[2] );    // Check if particle is crossed the boundary
            if ( gGeoManager->IsOutside() ) return rad;                 // leave if not
            TGeoMatrix *mat = gGeoManager->GetCurrentMatrix();          
            mat->MasterToLocal( pt, loc );
            if ( !gGeoManager->GetCurrentVolume()->Contains( loc ) ) {
                gGeoManager->CdUp();
                nextnode = gGeoManager->GetCurrentNode();               // move to new volume
            }
            continue;
        } else {
            ismall = 0;
        }
        
        // Normal steps case
        nbound++;
        length += snext;
        currentStep -= snext;
        if ( med ) {
            double radlen = med->GetMaterial()->GetRadLen() /*cm*/;
            if ( radlen > 1.e-9 && radlen < 1.e10 ) {
                
                lastrad = 1. / radlen * mm2cm;
                
                // Ignore start and finish volumes if required
                if ( skipBoundaryPonitsVolumes && isBoundaryVolume ) {
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
        }
    }
    
    streamlog_out(DEBUG1) << "--------EUTelGeometryTelescopeGeoDescription::findRadLengthIntegral()--------" << std::endl;
    
    return rad;
}

//
// straight line - shashlyk plane assembler
//
int EUTelGeometryTelescopeGeoDescription::findNextPlane(  double* lpoint,  double* ldir, float* newpoint )
{
// 
   if(newpoint==0)
   {
      streamlog_out ( ERROR0 ) << "::findNextPlane;  newpoint array is void, can not continue..."<<endl;
      return -100;
   }

   double normdir = TMath::Sqrt(ldir[0]*ldir[0]+ldir[1]*ldir[1]+ldir[2]*ldir[2]); 
   streamlog_out ( DEBUG0 ) << "::findNextPlane lpoint: "  << lpoint[0] << " " << lpoint[1] << " "<< lpoint[2] << " " << endl;
   ldir[0] = ldir[0]/normdir; 
   ldir[1] = ldir[1]/normdir; 
   ldir[2] = ldir[2]/normdir;
   streamlog_out ( DEBUG0 ) << "::findNextPlane ldir  : "  << ldir  [0] << " " << ldir  [1] << " "<< ldir  [2] << " " << endl;
 
   for(int ip=0;ip<3;ip++) 
   {
     newpoint[ip] = static_cast<float> (lpoint[ip]);
   }  
   int currentSensorID = getSensorID(newpoint); 
   
   gGeoManager->InitTrack( lpoint, ldir );
   TGeoNode *node = gGeoManager->GetCurrentNode( );

   Int_t inode    = node->GetIndex();
   Int_t i        = 0;

   streamlog_out ( DEBUG0 ) << "::findNextPlane look for next node, starting at node: " << node << " id: " << inode  << " currentSensorID: " << currentSensorID << endl;
 
//   double kStep = 1e-03;
   while(  node = gGeoManager->FindNextBoundaryAndStep(  ) )
   {
       inode = node->GetIndex();
       streamlog_out ( DEBUG0 ) << "::findNextPlane found next node: " << node << " id: " << inode << endl;
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
 
       streamlog_out ( DEBUG0 ) << "::findNextPlane i=" << i  << " " << inode << " " << ipoint[0]  << " " << ipoint[1] << " " << ipoint[2]  << " sensorID:" << sensorID <<  endl;
       if(sensorID >= 0 && sensorID != currentSensorID ) return sensorID;
   }

 
   return -100;
}
//This will take in a global coordinate and direction and output the new global point on the next sensor. Also it outputs the sesnorID if it matches output or -100 if not.
int EUTelGeometryTelescopeGeoDescription::findNextPlaneEntrance(  double* lpoint,  double* ldir, int nextSensorID, float* newpoint )
{
   if(newpoint==0)
   {
      streamlog_out ( ERROR0 ) << "::findNextPlaneEntrance newpoint array is void, can not continue..."<<endl;
      return -100;
   }
   
	//initialise direction and location in global telescope coordinates
   _geoManager->InitTrack( lpoint, ldir );
 
   TGeoNode *node = _geoManager->GetCurrentNode( ); //Return the volume i.e 'node' that contains that point.
   Int_t inode =  node->GetIndex();
   Int_t i=0;

   streamlog_out ( DEBUG0 ) << "::findNextPlaneEntrance node: " << node << " id: " << inode << endl;
 
	//Keep looping until you have left this plane volume and are at another. Note FindNextBoundaryAndStep will only take you to the next volume 'node' it will not enter it.
   while( node = _geoManager->FindNextBoundaryAndStep( ) )
   {
       inode = node->GetIndex();
       const double* point = _geoManager->GetCurrentPoint(); //This will be the new global coordinates after the move
       const double* dir   = _geoManager->GetCurrentDirection(); //This will be the same direction. If there was magnetic field then this would change automatically. However may need Geant4 for this to work. ???
       double ipoint[3] ;
       double idir[3]   ;

       for(int ip=0;ip<3;ip++) 
       {
         ipoint[ip] = point[ip];
         idir[ip]   = dir[ip];
         if(ip==2) ipoint[ip]+=0.001 ; // assumption !!! step by one um into the new volume // new volume is thicker than 1 um
         newpoint[ip] = static_cast<float> (ipoint[ip]);
       }  
       int sensorID = getSensorID(newpoint); 
       i++;     
      
       _geoManager->SetCurrentPoint( ipoint);
       _geoManager->SetCurrentDirection( idir);
 
       streamlog_out ( DEBUG0 ) << "Loop number" << i  << ". Index: " << inode << ". Current global point: " << ipoint[0]  << " " << ipoint[1] << " " << ipoint[2]  << " sensorID: " << sensorID << ". Input of expect next sensor: " << nextSensorID << endl;
       //if( sensorID <0 ) continue;  
       if( sensorID == nextSensorID ) return sensorID;
   }
 
   streamlog_out ( DEBUG0 ) << "::findNextPlaneEntrance node: " << node << " id: " << inode << " sensorID= " << nextSensorID << " not found" << " returning: 0" << endl;
 
   return -100;

}


