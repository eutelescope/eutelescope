/* 
 * File:   EUTelGeometryTelescopeGeoDescription.cpp
 * 
 */
#include "EUTelGeometryTelescopeGeoDescription.h"

// C++
#include <algorithm>
#include <string>
#include <cstring>

// MARLIN
#include "marlin/Global.h"
#include "marlin/VerbosityLevels.h"

//GEAR
#include "GEAR.h" //for GEAR exceptions

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

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;
using namespace geo;
using namespace std;


unsigned EUTelGeometryTelescopeGeoDescription::_counter = 0;

EUTelGeometryTelescopeGeoDescription& EUTelGeometryTelescopeGeoDescription::getInstance( gear::GearMgr* _g ) {
      streamlog_out ( DEBUG0) << "  EUTelGeometryTelescopeGeoDescription::getInstance --- BEGIN ---" << std::endl; 
      static  EUTelGeometryTelescopeGeoDescription instance;

      unsigned i = EUTelGeometryTelescopeGeoDescription::_counter;
      streamlog_out ( DEBUG0) << "  EUTelGeometryTelescopeGeoDescription::getInstance --- iter: " << i << std::endl;  
      if( i < 1 )  { // do it only once !
         instance.setGearManager(_g);
         instance.readGear();
      }
 
      instance.counter();
      streamlog_out ( DEBUG0) << "  EUTelGeometryTelescopeGeoDescription::getInstance --- END --- " << std::endl; 

      return instance;
}

size_t EUTelGeometryTelescopeGeoDescription::nPlanes( ) const {
    return _nPlanes;
}

const EVENT::DoubleVec& EUTelGeometryTelescopeGeoDescription::siPlanesZPositions( ) const {
    return _siPlaneZPosition;
}

/** set methods */ 
void EUTelGeometryTelescopeGeoDescription::setPlaneXPosition( int planeID, double value ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneXPosition[ _sensorIDtoZOrderMap[ planeID ] ] = value ;
}

void EUTelGeometryTelescopeGeoDescription::setPlaneYPosition( int planeID, double value ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneYPosition[ _sensorIDtoZOrderMap[ planeID ] ] = value ;
}

void EUTelGeometryTelescopeGeoDescription::setPlaneZPosition( int planeID, double value ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneZPosition[ _sensorIDtoZOrderMap[ planeID ] ] = value ;
}


void EUTelGeometryTelescopeGeoDescription::setPlaneXRotation( int planeID, double value ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneXRotation[ _sensorIDtoZOrderMap[ planeID ] ] = value ;
}

void EUTelGeometryTelescopeGeoDescription::setPlaneYRotation( int planeID, double value ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneYRotation[ _sensorIDtoZOrderMap[ planeID ] ] = value ;
}

void EUTelGeometryTelescopeGeoDescription::setPlaneZRotation( int planeID, double value ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneZRotation[ _sensorIDtoZOrderMap[ planeID ] ] = value ;
}

void EUTelGeometryTelescopeGeoDescription::setPlaneXRotationRadians( int planeID, double value /* in Radians */ ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneXRotation[ _sensorIDtoZOrderMap[ planeID ] ] = value * DEG;
}

void EUTelGeometryTelescopeGeoDescription::setPlaneYRotationRadians( int planeID, double value /* in Radians */ ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneYRotation[ _sensorIDtoZOrderMap[ planeID ] ] = value * DEG;
}

void EUTelGeometryTelescopeGeoDescription::setPlaneZRotationRadians( int planeID, double value /* in Radians */ ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() )  _siPlaneZRotation[ _sensorIDtoZOrderMap[ planeID ] ] = value * DEG ;
}

/** get methods */
float  EUTelGeometryTelescopeGeoDescription::siPlaneRotation1( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneRotation1[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}
float  EUTelGeometryTelescopeGeoDescription::siPlaneRotation2( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneRotation2[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}
float  EUTelGeometryTelescopeGeoDescription::siPlaneRotation3( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneRotation3[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}
float  EUTelGeometryTelescopeGeoDescription::siPlaneRotation4( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneRotation4[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}




/** get methods */
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

double/* in Radians */ EUTelGeometryTelescopeGeoDescription::siPlaneXRotationRadians( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneXRotation[ _sensorIDtoZOrderMap[ planeID ] ]* RADIAN;
    return -999.;
}

double/* in Radians */ EUTelGeometryTelescopeGeoDescription::siPlaneYRotationRadians( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneYRotation[ _sensorIDtoZOrderMap[ planeID ] ]* RADIAN;
    return -999.;
}

double/* in Radians */ EUTelGeometryTelescopeGeoDescription::siPlaneZRotationRadians( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneZRotation[ _sensorIDtoZOrderMap[ planeID ] ]* RADIAN;
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneXSize( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneXSize[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneYSize( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneYSize[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneZSize( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneZSize[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneRadLength( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneRadLength[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneXNpixels( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneXNpixels[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneYNpixels( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneYNpixels[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}


double EUTelGeometryTelescopeGeoDescription::siPlaneXPitch( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneXPitch[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneYPitch( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneYPitch[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}


double EUTelGeometryTelescopeGeoDescription::siPlaneXResolution( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneXResolution[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}

double EUTelGeometryTelescopeGeoDescription::siPlaneYResolution( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) return _siPlaneYResolution[ _sensorIDtoZOrderMap[ planeID ] ];
    return -999.;
}


std::string EUTelGeometryTelescopeGeoDescription::geoLibName( int planeID ) {

    std::map<int,int>::const_iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);

    int iplane = it->second ;
    std::string& s = _geoLibName[ iplane ];
    streamlog_out(DEBUG1) << s.c_str() << std::endl;
 
    if ( it != _sensorIDtoZOrderMap.end() ) return s;
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

void EUTelGeometryTelescopeGeoDescription::readSiPlanesLayout() {

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
	    std::cout << "No Geometry field found in GEAR file, assuming CAST for all planes" << std::endl;
    	    for(int i = 0; i < _nPlanes; i++)
            {
		    geometryNameParameters.push_back("CAST");
	    }
    }
	    
    setSiPlanesLayoutID( _siPlanesParameters->getSiPlanesID() ) ;
   
    // create an array with the z positions of each layer
    for (int iPlane = 0; iPlane < _nPlanes; iPlane++) {
        _siPlaneXPosition.push_back(_siPlanesLayerLayout->getLayerPositionX(iPlane));
        _siPlaneYPosition.push_back(_siPlanesLayerLayout->getLayerPositionY(iPlane));
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getLayerPositionZ(iPlane));
        _siPlaneXRotation.push_back(_siPlanesLayerLayout->getLayerRotationZY(iPlane));
        _siPlaneYRotation.push_back(_siPlanesLayerLayout->getLayerRotationZX(iPlane));
        _siPlaneZRotation.push_back(_siPlanesLayerLayout->getLayerRotationXY(iPlane));
 
        _siPlaneRotation1.push_back(_siPlanesLayerLayout->getSensitiveRotation1(iPlane));
        _siPlaneRotation2.push_back(_siPlanesLayerLayout->getSensitiveRotation2(iPlane));
        _siPlaneRotation3.push_back(_siPlanesLayerLayout->getSensitiveRotation3(iPlane));
        _siPlaneRotation4.push_back(_siPlanesLayerLayout->getSensitiveRotation4(iPlane));
      
	_siPlaneXSize.push_back(_siPlanesLayerLayout->getSensitiveSizeX(iPlane));
        _siPlaneYSize.push_back(_siPlanesLayerLayout->getSensitiveSizeY(iPlane));
        _siPlaneZSize.push_back(_siPlanesLayerLayout->getSensitiveThickness(iPlane));
 
        _siPlaneXNpixels.push_back(_siPlanesLayerLayout->getSensitiveNpixelX(iPlane));
        _siPlaneYNpixels.push_back(_siPlanesLayerLayout->getSensitiveNpixelY(iPlane));
        _siPlaneXPitch.push_back(_siPlanesLayerLayout->getSensitivePitchX(iPlane));
        _siPlaneYPitch.push_back(_siPlanesLayerLayout->getSensitivePitchY(iPlane));
        _siPlaneXResolution.push_back(_siPlanesLayerLayout->getSensitiveResolution(iPlane)); // should be ResolutionX
        _siPlaneYResolution.push_back(_siPlanesLayerLayout->getSensitiveResolution(iPlane)); // should be ResolutionY
        
        _siPlaneRadLength.push_back(_siPlanesLayerLayout->getSensitiveRadLength(iPlane));
	_geoLibName.push_back(geometryNameParameters[iPlane]);

    }


    // sort the array with increasing z
    std::sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());

    // clear the sensor ID vector
    _sensorIDVec.clear();

    // clear the sensor ID map
    _sensorIDVecMap.clear();
    _sensorIDtoZOrderMap.clear();


    for (int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++) {
        int sensorID = _siPlanesLayerLayout->getID(iPlane);
        _sensorIDVec.push_back(sensorID);
        _sensorIDVecMap.insert(std::make_pair(sensorID, iPlane));

        // count number of the sensors to the left of the current one:
        int sensorsToTheLeft = 0;
        int kposition = _siPlanesLayerLayout->getSensitivePositionZ(iPlane);
        for (int jPlane = 0; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++)
            if (_siPlanesLayerLayout->getSensitivePositionZ(jPlane) + 1e-06 < kposition  ) sensorsToTheLeft++;

        _sensorZOrderToIDMap.insert(std::make_pair(sensorsToTheLeft, sensorID));        
        _sensorIDtoZOrderMap.insert(std::make_pair(sensorID, sensorsToTheLeft));
    }


    _nPlanes = _siPlanesParameters->getSiPlanesNumber();

}

void EUTelGeometryTelescopeGeoDescription::readTrackerPlanesLayout() {

 
    // sensor-planes in geometry navigation:
    _trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&( _gearManager->getTrackerPlanesParameters()));
    _trackerPlanesLayerLayout = const_cast< gear::TrackerPlanesLayerLayout*> (&(_trackerPlanesParameters->getTrackerPlanesLayerLayout()));
    
    setSiPlanesLayoutID( _trackerPlanesParameters->getLayoutID() ) ;

    // clear the sensor ID vector
    _sensorIDVec.clear();
    // clear the sensor ID map
    _sensorIDVecMap.clear();
    _sensorIDtoZOrderMap.clear();

    // data memberS::

    _nPlanes = 0; // should be filed based on the length of the sensor vector.// after the loop

    // create an array with the z positions of each layer
    int nLayers = _trackerPlanesLayerLayout->getNLayers();
    for (int iLayer = 0; iLayer < nLayers; iLayer++) {
        gear::TrackerPlanesLayerImpl* _trackerPlanesLayerImpl = const_cast< gear::TrackerPlanesLayerImpl*>  (_trackerPlanesLayerLayout->getLayer( iLayer) );
        streamlog_out(DEBUG1) << " dlayer : " << iLayer << " " << nLayers  << " at " << _trackerPlanesLayerImpl << " " ;

        int nsensitive = _trackerPlanesLayerImpl->getNSensitiveLayers() ;
        streamlog_out(DEBUG1) << " contains " << nsensitive << " sensitive (sub)layers " << std::endl;

        gear::TrackerPlanesSensitiveLayerImplVec& vector = _trackerPlanesLayerImpl->getSensitiveLayerVec();
       
        for (int iSensLayer = 0; iSensLayer < nsensitive; iSensLayer++) {       

            gear::TrackerPlanesSensitiveLayerImpl& sensitiveLayer = vector.at(iSensLayer);
            int sensorID =   sensitiveLayer.getID();
            streamlog_out(DEBUG1) << " iLayer " << iLayer << " sens: " << iSensLayer << " id :" << sensorID << std::endl;

                   
            _siPlaneXPosition.push_back( sensitiveLayer.getPositionX() );
            _siPlaneYPosition.push_back( sensitiveLayer.getPositionY() );
            _siPlaneZPosition.push_back( sensitiveLayer.getPositionZ() );
            _siPlaneXRotation.push_back( sensitiveLayer.getRotationZY() );
            _siPlaneYRotation.push_back( sensitiveLayer.getRotationZX() );
            _siPlaneZRotation.push_back( sensitiveLayer.getRotationXY() );
        
            _siPlaneRotation1.push_back( 1.0 );
            _siPlaneRotation2.push_back( 0.0 );
            _siPlaneRotation3.push_back( 0.0 );
            _siPlaneRotation4.push_back( 1.0 );

	    _siPlaneXSize.push_back( sensitiveLayer.getSizeX() );
            _siPlaneYSize.push_back( sensitiveLayer.getSizeY() );
            _siPlaneZSize.push_back( sensitiveLayer.getThickness() );
 
            _siPlaneXNpixels.push_back( sensitiveLayer.getNpixelX() );
            _siPlaneYNpixels.push_back( sensitiveLayer.getNpixelY() );
            _siPlaneXPitch.push_back( sensitiveLayer.getPitchX() );
            _siPlaneYPitch.push_back( sensitiveLayer.getPitchY() );
            _siPlaneXResolution.push_back( sensitiveLayer.getResolutionX() );
            _siPlaneYResolution.push_back( sensitiveLayer.getResolutionY() );
        
            _siPlaneRadLength.push_back( sensitiveLayer.getRadLength() );
          
            _geoLibName.push_back( sensitiveLayer.getInfo().c_str() );

            _sensorIDVec.push_back(sensorID);
            _sensorIDVecMap.insert(std::make_pair(sensorID, iLayer)); // what if there are more then 1 sensore per layer?
            streamlog_out(DEBUG1) << " iter: " << _sensorIDVec.at( _sensorIDVec.size()-1 ) << " " << sensorID << " " << sensitiveLayer.getInfo() .c_str() << std::endl; 
        }
    }
 
    _nPlanes =  _sensorIDVec.size(); 
 
    for(int i=0; i< _siPlaneZPosition.size(); i++){ 
      int sensorsToTheLeft = 0;
      int sensorID = _sensorIDVec.at(i);

      for(int j=0; j< _siPlaneZPosition.size(); j++){ 
        if( _siPlaneZPosition.at(j) < _siPlaneZPosition.at(i) - 1e-06 ) sensorsToTheLeft++;
      }
       _sensorZOrderToIDMap.insert(std::make_pair(sensorsToTheLeft, sensorID));        
       _sensorIDtoZOrderMap.insert(std::make_pair(sensorID, sensorsToTheLeft));
    }
}

EUTelGeometryTelescopeGeoDescription::EUTelGeometryTelescopeGeoDescription() :
_siPlanesDefined (false),
_telPlanesDefined (false),
_gearManager( marlin::Global::GEAR ),
_siPlanesParameters(0),
_siPlanesLayerLayout(0),
_trackerPlanesParameters(0),
_trackerPlanesLayerLayout(0),
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
_siPlaneRotation1(),
_siPlaneRotation2(),
_siPlaneRotation3(),
_siPlaneRotation4(),
_siPlaneXSize(),
_siPlaneYSize(),
_siPlaneZSize(),
_siPlaneXPitch(),
_siPlaneYPitch(),
_siPlaneXNpixels(),
_siPlaneYNpixels(),
_siPlaneXResolution(),
_siPlaneYResolution(),
_siPlaneRadLength(),
_geoLibName(),
_nPlanes(0),
_isGeoInitialized(false),
_geoManager(0)
{
	//Set ROOTs verbosity to only display error messages or higher (so info will not be streamed to stderr)
	gErrorIgnoreLevel =  kError;  

	//Pixel Geometry manager creation
	_pixGeoMgr = new EUTelGenericPixGeoMgr();
}

void EUTelGeometryTelescopeGeoDescription::readGear() {
	
    if ( _gearManager == 0 ) {
        streamlog_out(ERROR2) << "The GearMgr is not available, for an unknown reason." << std::endl;
        throw eutelescope::InvalidGeometryException("GEAR manager is not initialised");
    }

    _siPlanesDefined = false;
    _telPlanesDefined = false;

    try{
      _siPlanesParameters = const_cast< gear::SiPlanesParameters*> (&(_gearManager->getSiPlanesParameters()));
      streamlog_out(MESSAGE1)  << "gear::SiPlanes : " << _siPlanesParameters << std::endl;
      _siPlanesDefined = true;
    }catch(...){
      streamlog_out(WARNING)   << "gear::SiPlanes NOT found " << std::endl;
    }

    try{
      _trackerPlanesParameters = const_cast< gear::TrackerPlanesParameters*> (&(_gearManager->getTrackerPlanesParameters()));
      streamlog_out(MESSAGE1)  << "gear::TrackerPlanes : " << _trackerPlanesParameters << std::endl;
      _telPlanesDefined = true;
    }catch(...){
      streamlog_out(WARNING)   << "gear::TrackerPlanes NOT found "  << std::endl;
    }

    if( _siPlanesDefined ){
      readSiPlanesLayout();
    }
    else if( _telPlanesDefined ){
      readTrackerPlanesLayout();
    }

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
    
    _geoManager = TGeoManager::Import( tgeofilename.c_str() );
    if( !_geoManager ) {
        streamlog_out( WARNING ) << "Can't read file " << tgeofilename << endl;
    }

    _geoManager->CloseGeometry();
}

/**
 *
 *
 */
void EUTelGeometryTelescopeGeoDescription::translateSiPlane2TGeo(TGeoVolume* pvolumeWorld, int SensorId ){

  
   double xc, yc, zc;   // volume center position 
   double alpha, beta, gamma;
   double rot1, rot2, rot3, rot4; // for backward compatibility with previous GEAR 

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

       rot1 = siPlaneRotation1( SensorId );
       rot2 = siPlaneRotation2( SensorId );
       rot3 = siPlaneRotation3( SensorId );
       rot4 = siPlaneRotation4( SensorId );

       // Spatial translations of the sensor center
       string stTranslationName = "matrixTranslationSensor";
       stTranslationName.append( strId.str() );
       TGeoTranslation* pMatrixTrans = new TGeoTranslation( stTranslationName.c_str(), xc, yc, zc );
       //ALL clsses deriving from TGeoMatrix are not owned by the ROOT geometry manager, invoking RegisterYourself() transfers
       //ownership and thus ROOT will clean up
       pMatrixTrans->RegisterYourself();      

// initial rot1,2,3,4 rotate and flip: 
       string stRotationName = "matrixRotationSensorRotatenFlip";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotFlip = new TGeoRotation( stRotationName.c_str(), 0.,  0., 0.);                // around X axis
       pMatrixRotFlip->RegisterYourself();
       if(   abs(rot2)<1e-6 &&    abs(rot3)<1e-6  ) {
          // do not rotate it's diagonal

         if(  abs(rot1+1.)<1e-6  ) {
          // flip X axis
          pMatrixRotFlip->ReflectX(0);
         }

         if(  abs(rot4+1.)<1e-6  ) {
          // flip Y axis
          pMatrixRotFlip->ReflectY(0);
         }
       }

       if(  abs(rot1)<1e-6 &&  abs(rot4)<1e-6 ) {
         // do rotate it's OFF diagonal
         pMatrixRotFlip->RotateZ(90.);
         // X->Y     0  1
         // Y->-X   -1  0

         if(  abs(rot2-1.)<1e-6 &&  abs(rot3+1.)<1e-6  ) {
          //         0  1
          //        -1  0
          // do nothing  
         }

         if(  abs(rot2-1.)<1e-6 &&  abs(rot3-1.)<1e-6  ) {
          //         0  1
          //        -1  0
          pMatrixRotFlip->ReflectY(0);
         }

         if(  abs(rot2+1.)<1e-6 &&  abs(rot3+1.)<1e-6  ) {
          //         0  1
          //        -1  0
          pMatrixRotFlip->ReflectX(0);
         }

         if(  abs(rot2+1.)<1e-6 &&  abs(rot3-1.)<1e-6  ) {
          //         0  1
          //        -1  0
          pMatrixRotFlip->ReflectX(0);
          pMatrixRotFlip->ReflectY(0);
         }
       }
      
       // Spatial rotation around sensor center
       // TGeoRotation requires Euler angles in degrees
       stRotationName = "matrixRotationSensorX";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotX = new TGeoRotation( stRotationName.c_str(), 0.,  alpha, 0.);                // around X axis
       pMatrixRotX->RegisterYourself();

       stRotationName = "matrixRotationSensorY";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotY = new TGeoRotation( stRotationName.c_str(), 90., beta,  0.);                // around Y axis (combination of rotation around Z axis and new X axis)
       stRotationName = "matrixRotationSensorBackY";
       pMatrixRotY->RegisterYourself();

       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotY1 = new TGeoRotation( stRotationName.c_str(), -90., 0.,  0.);                    // restoration of original orientation (valid in small angle approximataion ~< 15 deg)
       pMatrixRotY1->RegisterYourself(); 

       stRotationName = "matrixRotationSensorZ";
       stRotationName.append( strId.str() );
       TGeoRotation* pMatrixRotZ = new TGeoRotation( stRotationName.c_str(), 0. , 0.,        gamma);         // around Z axis
       pMatrixRotZ->RegisterYourself();
      
       // Combined rotation in several steps
       TGeoRotation* pMatrixRot = new TGeoRotation( *pMatrixRotX );
 
       pMatrixRot->MultiplyBy( pMatrixRotFlip );
       pMatrixRot->MultiplyBy( pMatrixRotY );
       pMatrixRot->MultiplyBy( pMatrixRotY1 );
       pMatrixRot->MultiplyBy( pMatrixRotZ );
       pMatrixRot->RegisterYourself();      
      

       // Combined translation and orientation
       TGeoCombiTrans* combi = new TGeoCombiTrans( *pMatrixTrans, *pMatrixRot );
       combi->RegisterYourself();   
 
       // Construction of sensor objects
       
       // Construct object medium. Required for radiation length determination

       // assume SILICON, though all information except of radiation length is ignored
       double a       = 28.085500;     
       double z       = 14.000000;
       double density = 2.330000;
       double radl    = siPlaneRadLength( SensorId );
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
       Double_t dx = siPlaneXSize( SensorId ) / 2.;
       Double_t dy = siPlaneYSize( SensorId ) / 2.;
       Double_t dz = siPlaneZSize( SensorId ) / 2.;
       TGeoShape *pBoxSensor = new TGeoBBox( "BoxSensor", dx, dy, dz );
       // Volume: volume_Sensor1
       
       // Geometry navigation package requires following names for objects that have an ID
       // name:ID
       string stVolName = "volume_SensorID:";
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

/** Determine id of the sensor in which point is locate
 * 
 * @param globalPos 3D point in global reference frame
 * @return sensorID or -999 if the point in outside of sensor volume
 */
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


void EUTelGeometryTelescopeGeoDescription::master2Localtwo(int sensorID, const double globalPos[], double localPos[] ) {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::master2Local() " << std::endl;

    _geoManager->cd( _planePath[sensorID].c_str() );
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
 * Retrieve magnetic field object.
 * 
 * @return reference to gear::BField object
 */
const gear::BField& EUTelGeometryTelescopeGeoDescription::getMagneticFiled() const {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getMagneticFiled() " << std::endl;
    return _gearManager->getBField();
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
    double epsil = 1.E-7;
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
                streamlog_out( WARNING1 ) << "ERROR: Small steps in: " << gGeoManager->GetPath() << " shape=" << shape->ClassName() << endl;
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

/**
* Find closest surface intersected by the track and propagate track to that point
* @param input: ts track state
* @param input: The plane you want to find the intersection.
* @param input: Pointer to fill with the new global coordinates   
* @return planeID. If there was a problem return -999.
*/
int EUTelGeometryTelescopeGeoDescription::findIntersectionWithCertainID( float x0, float y0, float z0, float px, float py, float pz, float _beamQ, int nextPlaneID, float* output) {
streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::findIntersection()" << std::endl;
 
	// Set position and momentum vector//////////////////////////////////
  TVector3 trkVec(x0,y0,z0);
	TVector3 pVec(px,py,pz);
	/////////////////////////////////////////////////////////////////////

	streamlog_out(DEBUG5) << "  Global positions: "<< x0 <<"  "<< y0 <<"  "<< z0 << " Momentum: "<< pVec[0]<<","<<pVec[1]<<","<<pVec[2]<<","<< std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////  

 
  // Find magnetic field at that point and then the components/////////////////////////////////// 
  gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction. Why do we need this assumption. Equations of motion do not seem to dictate this.
  const gear::BField&   B = geo::gGeometry().getMagneticFiled();
	const double bx         = B.at( vectorGlobal ).x();
	const double by         = B.at( vectorGlobal ).y();
	const double bz         = B.at( vectorGlobal ).z();
  TVector3 hVec(bx,by,bz);
	const double H = hVec.Mag();
  //////////////////////////////////////////////////////////////////////////////////////////////

  // Calculate track momentum from track parameters and fill some useful variables///////////////////////////////////////////////////////////
  const double p = pVec.Mag();
	const double mm = 1000.;
  const double k = -0.299792458/mm*_beamQ*H;
  const double rho = k/p; 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////       
				      
	//Determine geometry of sensor to be used to determine the point of intersection.//////////////////////////////////////
  TVector3 norm = geo::gGeometry().siPlaneNormal( nextPlaneID  );       
  TVector3 sensorCenter( geo::gGeometry().siPlaneXPosition( nextPlaneID  ), geo::gGeometry().siPlaneYPosition( nextPlaneID  ), geo::gGeometry().siPlaneZPosition( nextPlaneID  ) );
  TVector3 delta = trkVec - sensorCenter;
  TVector3 pVecCrosH = pVec.Cross( hVec.Unit() );
	////////////////////////////////////////////////////////////////////////////////////////////////////////


  if ( streamlog_level(DEBUG5) ) {
		streamlog_out (DEBUG5) << "-------------------------------------------" << std::endl;
	  streamlog_out (DEBUG5) << "Current point (X,Y,Z): " << std::setw(15) << x0  << std::setw(15) << y0 << std::setw(15) << z0 << std::endl;
	  streamlog_out (DEBUG5) << "Next PlaneID : " << nextPlaneID << std::endl;
	  streamlog_out (DEBUG5) << "Normal vector" << std::endl;
	  norm.Print();
	  streamlog_out (DEBUG5) << "P x H vector" << std::endl;
	  pVecCrosH.Print();
	  streamlog_out (DEBUG5) << "Rho: " << rho << std::endl;
	  streamlog_out (DEBUG5) << "P: " << p << std::endl;
  }


	//Solution to the plane equation and the curved line intersection will be an quadratic with the coefficients. The solution is the arc length along the curve
	const double a = -0.5 * rho * ( norm.Dot( pVecCrosH ) ) / p;
  const double b = norm.Dot( pVec ) / p;
  const double c = norm.Dot( delta );
	//////////////////////////////////////////////////////////////////////////////////////////////////////// 
   
  std::vector< double > sol = Utility::solveQuadratic(a,b,c); // solutions are sorted in ascending order. This is a vector of arc length
	double solution = ( sol[0] > 0. ) ? sol[0] : ( ( sol[0] < 0. && sol[1] > 0. ) ? sol[1] : -1. ); //choose solution with minimum arc length
	if(solution < 0){
		streamlog_out ( DEBUG3 ) << "Track intersection was not found" << std::endl;
		return -999;
	}
			
	//Determine the global position from arc length.             
  TVector3 newPos;
	newPos = getXYZfromArcLength(x0, y0,z0,px,py,pz,_beamQ,solution);
	output[0]=newPos[0]; 				output[1]=newPos[1]; 				output[2]=newPos[2];
				
	streamlog_out (DEBUG5) << "Solutions for arc length: " << std::setw(15) << sol[0] << std::setw(15) << sol[1] << std::endl;
	streamlog_out (DEBUG5) << "Final solution (X,Y,Z): " << std::setw(15) << output[0]  << std::setw(15) << output[1]  << std::setw(15) << output[2] << std::endl;

        
  streamlog_out(DEBUG2) << "-------------------------EUTelGeometryTelescopeGeoDescription::findIntersection()--------------------------" << std::endl;
        
  return nextPlaneID;
}
//This function determined the xyz position in global coordinates using the state and arc length of the track s.
TVector3 EUTelGeometryTelescopeGeoDescription::getXYZfromArcLength( float x0, float y0, float z0, float px, float py, float pz, float _beamQ, float s) const {
	streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getXYZfromArcLength()" << std::endl;

  // Fill the postion and momentun into vector
	TVector3 pos( x0, y0, z0 );
	TVector3 pVec(px, py, pz );	
	//////////////////////////////////////////////////////
                
  // Get magnetic field vector
  gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
  const double bx         = B.at( vectorGlobal ).x();
  const double by         = B.at( vectorGlobal ).y();
  const double bz         = B.at( vectorGlobal ).z();
  TVector3 hVec(bx,by,bz);
               
	const double H = hVec.Mag();
  const double p = pVec.Mag();
	const double mm = 1000.;
 	const double k = -0.299792458/mm*_beamQ*H;
  const double rho = k/p;
        
	if ( fabs( k ) > 1.E-6  ) {
		// Non-zero magnetic field case
		TVector3 pCrossH = pVec.Cross(hVec.Unit());
		TVector3 pCrossHCrossH = pCrossH.Cross(hVec.Unit());
		const double pDotH = pVec.Dot(hVec.Unit());
		TVector3 temp1 = pCrossHCrossH;	temp1 *= ( -1./k * sin( rho * s ) );
		TVector3 temp2 = pCrossH;       temp2 *= ( -1./k * ( 1. - cos( rho * s ) ) );
		TVector3 temp3 = hVec;          temp3 *= ( pDotH / p * s );
		pos += temp1;
		pos += temp2;
		pos += temp3;
        } else {
		// Vanishing magnetic field case

		
		const double cosA =  px/p;      // Calculate cos of the angle between Z(beam) and X(solenoid field axis) //NEED TO MAKE SURE THAT TX=PX/P
		const double cosB = py/p ;        // Calculate cos of the angle between Z(beam) and Y
		pos.SetX( x0 + cosA * s );
		pos.SetY( y0 + cosB * s );
		pos.SetZ( z0 + 1./p * pVec.Z() * s );
	}
        
	streamlog_out(DEBUG2) << "---------------------------------EUTelKalmanFilter::getXYZfromArcLength()------------------------------------" << std::endl;
        
	return pos;
}

//This function given position/momentum of a particle. Will give you the approximate jacobian at any point along the track. This effectively relates changes in the particle position/momentum at the original to some distant point. 
//So if I change the initial position by x amount how much will all the other variables position/momentum at the new position change? This is what the Jacobian tells you.

   /** Calculate track parameters propagation jacobian for given track state
     *  and propagation distance. The expressions were derived in parabolic approximation
     *  valid for small values of propagation distance |dz| < 10cm. Can be iterated if necessary.
     * 
     * @param ts track state
     * @param dz propagation distance
     * @return 
     */
  TMatrix EUTelGeometryTelescopeGeoDescription::getPropagationJacobianF( float x0, float y0, float z0, float px, float py, float pz, float _beamQ, float dz ) {
        streamlog_out( DEBUG2 ) << "EUTelGeometryTelescopeGeoDescription::getPropagationJacobianF()" << std::endl;
	// The formulas below are derived from equations of motion of the particle in
        // magnetic field under assumption |dz| small. Must be valid for |dz| < 10 cm

	const double mm = 1000.;
	const double k = 0.299792458/mm;

	TVector3 pVec(px, py, pz );	

	// Get track parameters
	const double invP = _beamQ/pVec.Mag();
        const double tx0 = px/pVec.Mag(); //NEED TO DOUBLE CHECK THAT TX = PX/P
        const double ty0 = py/pVec.Mag();

        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
        const double Bx         = B.at( vectorGlobal ).x();
        const double By         = B.at( vectorGlobal ).y();
        const double Bz         = B.at( vectorGlobal ).z();
        
        const double sqrtFactor = sqrt( 1. + tx0*tx0 + ty0*ty0 );

	const double Ax = sqrtFactor * (  ty0 * ( tx0 * Bx + Bz ) - ( 1. + tx0*tx0 ) * By );
	const double Ay = sqrtFactor * ( -tx0 * ( ty0 * By + Bz ) + ( 1. + ty0*ty0 ) * Bx );

	// Partial derivatives
	//const double dAxdtx0 = tx0 * Ax / (sqrtFactor*sqrtFactor) + sqrtFactor*( ty0*Bx - 2. * tx0 * By );
	const double dAxdty0 = ty0 * Ax / (sqrtFactor*sqrtFactor) + sqrtFactor*( tx0*Bx + Bz );
	const double dAydtx0 = tx0 * Ay / (sqrtFactor*sqrtFactor) + sqrtFactor*( -ty0*By - Bz );
	//const double dAydty0 = ty0 * Ay / (sqrtFactor*sqrtFactor) + sqrtFactor*( -tx0*By + 2. * ty0 * Bx );

	const double dxdtx0 = dz;
	const double dxdty0 = 0.5 * invP * k * dz*dz * dAxdty0;

	const double dydtx0 = 0.5 * invP * k * dz*dz * dAydtx0;
	const double dydty0 = dz;

	const double dtxdty0 = invP * k * dz * dAxdty0;
	const double dtydtx0 = invP * k * dz * dAydtx0;

	const double dxdinvP0 = 0.5 * k * dz*dz * Ax;
	const double dydinvP0 = 0.5 * k * dz*dz * Ay;

	const double dtxdinvP0 = k * dz * Ax;
	const double dtydinvP0 = k * dz * Ay;

	// Fill-in matrix elements
	TMatrix jacobianF(5,5);
	jacobianF.UnitMatrix();
	jacobianF[0][2] = dxdtx0;	jacobianF[0][3] = dxdty0;	jacobianF[0][4] = dxdinvP0;
	jacobianF[1][2] = dydtx0;	jacobianF[1][3] = dydty0;	jacobianF[1][4] = dydinvP0;
	jacobianF[2][3] = dtxdty0;	jacobianF[2][4] = dtxdinvP0;
	jacobianF[3][2] = dtydtx0;	jacobianF[3][4] = dtydinvP0;
        
        if ( streamlog_level(DEBUG0) ){
             streamlog_out( DEBUG0 ) << "Propagation jacobian: " << std::endl;
            jacobianF.Print();
        }
	
        streamlog_out( DEBUG2 ) << "-----------------------------EUTelGeometryTelescopeGeoDescription::getPropagationJacobianF()-------------------------------" << std::endl;

	return jacobianF;
        
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
    for (int iPlane = 0; iPlane < _nPlanes; iPlane++) {
        int sensorID =  _sensorIDVec.at(iPlane);
        std::cout << " iplane " << iPlane << " [" << _nPlanes << "] id :" <<  setw(3) << sensorID  ;
        std::cout << " xPos: " << setw(10) << setprecision(4) <<  _siPlaneXPosition.at( iPlane);
        std::cout << " yPos: " << setw(10) << setprecision(4) <<  _siPlaneYPosition.at( iPlane);
        std::cout << " zPos: " << setw(10) << setprecision(4) <<  _siPlaneZPosition.at( iPlane);
        std::cout << " xRot: " << setw(10) << setprecision(4) <<  _siPlaneXRotation.at( iPlane);
        std::cout << " yRot: " << setw(10) << setprecision(4) <<  _siPlaneYRotation.at( iPlane);
        std::cout << " zRot: " << setw(10) << setprecision(4) <<  _siPlaneZRotation.at( iPlane);
        std::cout << "       "                                   << std::endl;

        siplanesLayerLayout->setLayerPositionX(  iPlane, _siPlaneXPosition.at(iPlane) );
        siplanesLayerLayout->setLayerPositionY(  iPlane, _siPlaneYPosition.at(iPlane) );
        siplanesLayerLayout->setLayerPositionZ(  iPlane, _siPlaneZPosition.at(iPlane) );
        siplanesLayerLayout->setLayerRotationZY( iPlane, _siPlaneXRotation.at(iPlane) );
        siplanesLayerLayout->setLayerRotationZX( iPlane, _siPlaneYRotation.at(iPlane) );
        siplanesLayerLayout->setLayerRotationXY( iPlane, _siPlaneZRotation.at(iPlane) );

    }


    // ------- add to GearMgr ----
    if( _gearManager != 0 ) {
      
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
 
            for( int iplane = 0; iplane < _sensorIDVec.size(); iplane++ ) {
              int sensorID =  _sensorIDVec.at(iplane);
              if( sensitiveLayer.getID() !=  _sensorIDVec.at( iplane) ) continue;  
 
              std::cout << " iLayer " << setw(3)<< iLayer << "["<< setw(3) << nLayers <<"] sens: " << setw(3) << iSensLayer << "[" << setw(3) << nsensitive << "] id :" <<  setw(3) << sensorID  ;
              std::cout << " xPos: " << setw(10) << setprecision(4) << _siPlaneXPosition.at( iplane);
              std::cout << " yPos: " << setw(10) << setprecision(4) <<  _siPlaneYPosition.at( iplane);
              std::cout << " zPos: " << setw(10) << setprecision(4) <<  _siPlaneZPosition.at( iplane);
              std::cout << " xRot: " << setw(10) << setprecision(4) <<  _siPlaneXRotation.at( iplane);
              std::cout << " yRot: " << setw(10) << setprecision(4) <<  _siPlaneYRotation.at( iplane);
              std::cout << " zRot: " << setw(10) << setprecision(4) <<  _siPlaneZRotation.at( iplane);
              std::cout << "       "                                   << std::endl;

              sensitiveLayer.setPositionX( _siPlaneXPosition.at( iplane) );
              sensitiveLayer.setPositionY( _siPlaneYPosition.at( iplane) );
              sensitiveLayer.setPositionZ( _siPlaneZPosition.at( iplane) );

              sensitiveLayer.setRotationZY( _siPlaneXRotation.at( iplane) );
              sensitiveLayer.setRotationZX( _siPlaneYRotation.at( iplane) );
              sensitiveLayer.setRotationXY( _siPlaneZRotation.at( iplane) );

            }
                
   
        }
    }

    // ------- add to GearMgr ----
    if( _gearManager != 0 ) {
      
     _gearManager->setTrackerPlanesParameters( trackerplanesParameters ) ;

    }


    streamlog_out( MESSAGE1 ) << "EUTelGeometryTelescopeGeoDescription::updateTrackerPlanesLayout() --- OVER ---- " << std::endl;
}

void EUTelGeometryTelescopeGeoDescription::updateGearManager() {

 streamlog_out( MESSAGE1 ) << "EUTelGeometryTelescopeGeoDescription::updateGearManager() --- START ---- " << std::endl;

 if( _siPlanesDefined ){
   updateSiPlanesLayout();
 }
 else if( _telPlanesDefined ){
   updateTrackerPlanesLayout();
 }

 streamlog_out( MESSAGE1 ) << "EUTelGeometryTelescopeGeoDescription::updateGearManager() --- OVER ---- " << std::endl;
}

