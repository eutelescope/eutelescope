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
     // streamlog_out ( DEBUG4) << "  EUTelGeometryTelescopeGeoDescription::getInstance --- BEGIN ---" << std::endl; 
      static  EUTelGeometryTelescopeGeoDescription instance;

      unsigned i = EUTelGeometryTelescopeGeoDescription::_counter;
    //  streamlog_out ( DEBUG4) << "  EUTelGeometryTelescopeGeoDescription::getInstance --- iter: " << i << std::endl;  
      if( i < 1 )  { // do it only once !
         instance.setGearManager(_g);
         instance.readGear();
      }
 
      instance.counter();
    //  streamlog_out ( DEBUG4) << "  EUTelGeometryTelescopeGeoDescription::getInstance --- END --- " << std::endl; 

      return instance;
}

size_t EUTelGeometryTelescopeGeoDescription::nPlanes( ) const {
    return _nPlanes;
}

const EVENT::DoubleVec& EUTelGeometryTelescopeGeoDescription::siPlanesZPositions( ) const {
    return _siPlaneZPosition;
}

/** set methods */ 
void EUTelGeometryTelescopeGeoDescription::setInitialDisplacementToFirstPlane(float initialDisplacement){
	_initialDisplacement = initialDisplacement;
}
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
float EUTelGeometryTelescopeGeoDescription::getInitialDisplacementToFirstPlane() const {
	return _initialDisplacement;
}
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
//Note  that to determine these axis we MUST use the geometry class after initialisation. By this I mean directly from the root file create.
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneNormal( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) {
				const double zAxisLocal[3]  = {0,0,1};
				double zAxisGlobal[3]; 
				local2MasterVec(planeID,zAxisLocal , zAxisGlobal ); 
        TVector3 normVec(zAxisGlobal[0],zAxisGlobal[1],zAxisGlobal[2]);
        return normVec;
    }else{
			TVector3 defaultZ(0,0,1);
			return defaultZ;
		}
}
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneXAxis( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) {
			const double xAxisLocal[3]  = {1,0,0};
			double xAxisGlobal[3]; 
			local2MasterVec(planeID,xAxisLocal , xAxisGlobal ); 
			TVector3 xVec(xAxisGlobal[0],xAxisGlobal[1],xAxisGlobal[2]);
			return xVec;
    }else{
			TVector3 defaultX(1,0,0);
			return defaultX;
		}
}
TVector3 EUTelGeometryTelescopeGeoDescription::siPlaneYAxis( int planeID ) {
    std::map<int,int>::iterator it;
    it = _sensorIDtoZOrderMap.find(planeID);
    if ( it != _sensorIDtoZOrderMap.end() ) {
			const double yAxisLocal[3]  = {0,1,0};
			double yAxisGlobal[3]; 
			local2MasterVec(planeID,yAxisLocal , yAxisGlobal ); 
			TVector3 yVec(yAxisGlobal[0],yAxisGlobal[1],yAxisGlobal[2]);
			return yVec;
    }else{
			TVector3 defaultY (0,1,0);
			return defaultY;
	  }
}
const std::map<int, int>& EUTelGeometryTelescopeGeoDescription::sensorIDstoZOrder( ) const {
    return _sensorIDtoZOrderMap;
}
std::map<int, int>& EUTelGeometryTelescopeGeoDescription::sensorZOrderToIDWithoutExcludedPlanes() {
	return _sensorZOrderToIDWithoutExcludedPlanes;
}
std::map<int, int>& EUTelGeometryTelescopeGeoDescription::sensorIDToZOrderWithoutExcludedPlanes() {
	return _sensorIDToZOrderWithoutExcludedPlanes;
}
void EUTelGeometryTelescopeGeoDescription::initialisePlanesToExcluded(FloatVec planeIDs ){
	int counter=0;
	for(int i = 0 ; i <_sensorZOrderToIDMap.size(); ++i){
		bool excluded=false;
		for(int  j =0; j< planeIDs.size(); ++j){
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
		throw(lcio::Exception( Utility::outputColourString("The number of Planes-Excluded is not correct. This could be a problem with geometry", "RED")));
	}else{
		streamlog_out(DEBUG0) << Utility::outputColourString("The correct number of planes have been excluded","GREEN") << std::endl;
	}
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
       			//Note we set rotations 1,2,3,4 automatically here since we don't want to add more gear information that uses this format. Only the old Si planes should use that.  
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
	double rot1, rot3; // for backward compatibility with previous GEAR. We only need 2 entries from gear file for fast z rotation function in TGeoRotations. 

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
	rot3 = siPlaneRotation3( SensorId );
	//rot1=sin and rot=cos
	//We must check that the input is correct. Since there is only 1 degree of freedom.
	float	possibleUnit = sqrt(pow(rot1,2)+pow(rot3,2));
	if(possibleUnit>1.1 or possibleUnit<0.9){//Check that the squares equal 1 and give some rounding error room. 
		streamlog_out(ERROR5) << "SensorID: " << SensorId << ". Rot Squares sum =  " <<possibleUnit << std::endl;   
		throw(lcio::Exception(Utility::outputColourString("The integer rotations do not represent a z rotation since there squares don't sum to 1. Note we have four variables in gear file as input BUT A SINGLE DEGREE OF FREEDOM!", "RED"))); 	
	}
	double integerRotations[2]={rot3,rot1};
	//Create spatial TGeoTranslation object.
	string stTranslationName = "matrixTranslationSensor";
	stTranslationName.append( strId.str() );
	TGeoTranslation* pMatrixTrans = new TGeoTranslation( stTranslationName.c_str(), xc, yc, zc );
	//ALL clsses deriving from TGeoMatrix are not owned by the ROOT geometry manager, invoking RegisterYourself() transfers
	//ownership and thus ROOT will clean up
	pMatrixTrans->RegisterYourself();      

	//Create TGeoRotation object. 
	//Translations are of course just positional changes in the global frame.
	//Note that each subsequent rotation is using the new coordinate system of the last transformation all the way back to the global frame.
	//The order is:
	//Integer Z rotation. This only uses the top left element. This for future gear types should be removed.
	//Z rotations specified by in degrees.
	//X rotations 
	//Y rotations
	TGeoRotation * pMatrixRotCombined = new TGeoRotation();
	pMatrixRotCombined->FastRotZ(integerRotations);//Z Rotation (Integer). This will rotate a vector around the z axis using the right hand rule
	pMatrixRotCombined->RotateZ(gamma);//Z Rotation (degrees)//This will again rotate a vector around z axis usign the right hand rule.  
	pMatrixRotCombined->RotateX(alpha);//X Rotations (degrees)//This will rotate a vector usign the right hand rule round the x-axis
	pMatrixRotCombined->RotateY(beta);//Y Rotations (degrees)//Same again for Y axis 
	pMatrixRotCombined->RegisterYourself();//We must allow the matrix to be used by the TGeo manager.
	// Combined translation and orientation
	TGeoCombiTrans* combi = new TGeoCombiTrans( *pMatrixTrans, *pMatrixRotCombined );
	//This is to print to screen the rotation and translation matrices used to transform from local to global frame.
	streamlog_out(MESSAGE9) << "THESE MATRICES ARE USED TO TAKE A POINT IN THE LOCAL FRAME AND MOVE IT TO THE GLOBAL FRAME."  << std::endl;   
	streamlog_out(MESSAGE9) << "SensorID: " << SensorId << " Rotation matrix for this object."  << std::endl;   
	const double* rotationMatrix =  combi->GetRotationMatrix();	
	streamlog_out(MESSAGE9) << setw(10) <<rotationMatrix[0]<<"  "<<rotationMatrix[1]<<"   "<<rotationMatrix[2]<<endl;
	streamlog_out(MESSAGE9) << setw(10) <<rotationMatrix[3]<<"  "<<rotationMatrix[4]<<"   "<<rotationMatrix[5]<<endl;
	streamlog_out(MESSAGE9) << setw(10) <<rotationMatrix[6]<<"  "<<rotationMatrix[7]<<"   "<<rotationMatrix[8]<<endl;

	//streamlog_out(MESSAGE9) << setw(10) <<rotationMatrix[0] << setw(10) <<rotationMatrix[1]<< setw(10) <<rotationMatrix[2]<< setw(10)<<endl<<endl; 
	//streamlog_out(MESSAGE9) << setw(10) <<rotationMatrix[3] << setw(10) <<rotationMatrix[4]<< setw(10) <<rotationMatrix[5]<< setw(10)<<endl<<endl; 
	//streamlog_out(MESSAGE9) << setw(10) <<rotationMatrix[6] << setw(10) <<rotationMatrix[7]<< setw(10) <<rotationMatrix[8]<< setw(10)<<endl<<endl; 
	const double* translationMatrix =  combi->GetTranslation();	
	streamlog_out(MESSAGE9) << "SensorID: " << SensorId << " Translation vector for this object."  << std::endl;   
	streamlog_out(MESSAGE9) << setw(10) <<translationMatrix[0] << setw(10) <<translationMatrix[1]<< setw(10) <<translationMatrix[2]<< setw(10)<<endl; 

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
 * Retrieve magnetic field object.
 * 
 * @return reference to gear::BField object
 */
const gear::BField& EUTelGeometryTelescopeGeoDescription::getMagneticFiled() const {
    streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getMagneticFiled()----BEGIN " << std::endl;
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
* Find closest surface intersected by the track and return the position
*/
int EUTelGeometryTelescopeGeoDescription::findIntersectionWithCertainID( float x0, float y0, float z0, float px, float py, float pz, float beamQ, int nextPlaneID, float outputPosition[],TVector3& outputMomentum, float& arcLength) {
streamlog_out(DEBUG5) << "EUTelGeometryTelescopeGeoDescription::findIntersectionWithCertainID()------BEGIN" << std::endl;
	//positions are in mm
  TVector3 trkVec(x0,y0,z0);
	TVector3 pVec(px,py,pz);
	streamlog_out(DEBUG5) << "  Global positions (Input): "<< x0 <<"  "<< y0 <<"  "<< z0 << " Momentum: "<< px<<","<<py<<","<<pz<< std::endl;
 
  // Find magnetic field at that point and then the components/////////////////////////////////// 
  gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction. Why do we need this assumption. Equations of motion do not seem to dictate this.
  const gear::BField&   B = geo::gGeometry().getMagneticFiled();
	const double bx         = B.at( vectorGlobal ).x();
	const double by         = B.at( vectorGlobal ).y();
	const double bz         = B.at( vectorGlobal ).z();
	streamlog_out (DEBUG5) << "The magnetic field vector (x,y,z): "<<bx<<","<<by<<","<<bz << std::endl;
	streamlog_out (DEBUG5) << "Beam charge: "<<beamQ<< std::endl;

	//B field is in units of Tesla
  TVector3 hVec(bx,by,bz);
	const double H = hVec.Mag();

  // Calculate track momentum from track parameters and fill some useful variables///////////////////////////////////////////////////////////
  const double p = pVec.Mag();
	const double constant =  -0.299792458; //This is a constant used in the derivation of this equation. This is the distance light travels in a nano second    
	const double mm = 1000;
	const double combineConstantsAndMagneticField = (constant*beamQ*H)/mm;
  const double rho = combineConstantsAndMagneticField/p; 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////       
				      
	//Determine geometry of sensor to be used to determine the point of intersection.//////////////////////////////////////
  TVector3 norm = geo::gGeometry().siPlaneNormal( nextPlaneID  );       
	streamlog_out (DEBUG5) << "The normal of the plane is (x,y,z): "<<norm[0]<<","<<norm[1]<<","<<norm[2]<< std::endl;
  TVector3 sensorCenter( geo::gGeometry().siPlaneXPosition( nextPlaneID  ), geo::gGeometry().siPlaneYPosition( nextPlaneID  ), geo::gGeometry().siPlaneZPosition( nextPlaneID  ) );
	streamlog_out (DEBUG5) << "The sensor centre (x,y,z): "<<sensorCenter[0]<<","<<sensorCenter[1]<<","<<sensorCenter[2] << std::endl;
  TVector3 delta = (trkVec - sensorCenter);//Must be in mm since momentum is.
  TVector3 pVecCrosH = pVec.Cross(hVec.Unit());
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
	  streamlog_out (DEBUG5) << "Delta "<< std::endl;
		delta.Print();
  }
	//Solution to the plane equation and the curved line intersection will be an quadratic with the coefficients. The solution is the arc length along the curve
	const double a = -0.5 * rho * ( norm.Dot( pVecCrosH ) ) / p;
  const double b = norm.Dot( pVec ) / p;
  const double c = norm.Dot( delta );
	//////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //Solution must be in femto metres 
  std::vector< double > sol = Utility::solveQuadratic(a,b,c); // solutions are sorted in ascending order. This is a vector of arc length
	double solution = ( sol[0] > 0. ) ? sol[0] : ( ( sol[0] < 0. && sol[1] > 0. ) ? sol[1] : -1. ); //choose solution with minimum arc length
	if(solution < 0){
		streamlog_out ( DEBUG3 ) << "Track intersection was not found" << std::endl;
		return -999;
	}
			
	//Determine the global position from arc length.             
  TVector3 newPos;
	TVector3 newMomentum;
	newPos = getXYZfromArcLength(trkVec,pVec,beamQ,solution);
	newMomentum = getXYZMomentumfromArcLength(pVec, trkVec, beamQ, solution);
	outputPosition[0]=newPos[0]; 				outputPosition[1]=newPos[1]; 				outputPosition[2]=newPos[2];
	outputMomentum[0]=newMomentum[0]; 				outputMomentum[1]=newMomentum[1]; 				outputMomentum[2]=newMomentum[2];
	arcLength = solution;		
	streamlog_out (DEBUG5) << "Solutions for arc length: " << std::setw(15) << sol[0] << std::setw(15) << sol[1] << " Final output arc length: " << arcLength <<std::endl;
	streamlog_out (DEBUG5) << "Final Solution momentum(X,Y,Z): " << std::setw(15) << outputMomentum[0]  << std::setw(15) << outputMomentum[1]  << std::setw(15) << outputMomentum[2] << std::endl;
	streamlog_out (DEBUG5) << "Final solution (X,Y,Z): " << std::setw(15) << outputPosition[0]  << std::setw(15) << outputPosition[1]  << std::setw(15) << outputPosition[2] << std::endl;

        
  streamlog_out(DEBUG2) << "-------------------------EUTelGeometryTelescopeGeoDescription::findIntersection()--------------------------" << std::endl;
        
  return nextPlaneID;
}
//This will calculate the momentum at a arc length away given initial parameters.
TVector3 EUTelGeometryTelescopeGeoDescription::getXYZMomentumfromArcLength(TVector3 momentum, TVector3 globalPositionStart, float charge, float arcLength ){
	float mm= 1000;
	TVector3 T = momentum.Unit();//This is one coordinate axis of curvilinear coordinate system.	
	const gear::BField&   Bfield = geo::gGeometry().getMagneticFiled();
	gear::Vector3D vectorGlobal(globalPositionStart[0],globalPositionStart[1],globalPositionStart[2]);//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
	const double Bx = (Bfield.at( vectorGlobal ).x());//We times bu 0.3 due to units of other variables. See paper. Must be Tesla
	const double By = (Bfield.at( vectorGlobal ).y());
	const double Bz = (Bfield.at( vectorGlobal ).z());
	TVector3 B;
	B[0]=Bx*0.3; B[1]=By*0.3; B[2]=Bz*0.3;
	TVector3 H = (B.Unit());
	const float alpha = (H.Cross(T)).Mag();
	const float gamma = H.Dot(T);
	const float Q = -(B.Mag())*(charge/(momentum.Mag()));//You could use end momentum since it must be constant
	float theta = (Q*arcLength)/mm;//divide by 1000 to convert to metres 
	TVector3 N = (H.Cross(T)).Unit();
	const float cosTheta = cos(theta);
	const float sinTheta = sin(theta);
	const float oneMinusCosTheta = (1-cos(theta));
	TVector3 momentumEndUnit = gamma*oneMinusCosTheta*H+cosTheta*T+alpha*sinTheta*N;
	streamlog_out ( DEBUG0 ) << "Momentum direction (Unit Vector): " <<  momentumEndUnit[0]<<" , "<<momentumEndUnit[1] <<" , "<<momentumEndUnit[2]<<std::endl;
	TVector3 momentumEnd = momentumEndUnit*(momentum.Mag());
	streamlog_out ( DEBUG0 ) << "Momentum: " <<  momentumEnd[0]<<" , "<<momentumEnd[1] <<" , "<<momentumEnd[2]<<std::endl;

	return momentumEnd;
}
//This function determined the xyz position in global coordinates using the state and arc length of the track s.
TVector3 EUTelGeometryTelescopeGeoDescription::getXYZfromArcLength(TVector3 pos ,TVector3 pVec , float beamQ, double s) const {
	streamlog_out(DEBUG2) << "EUTelGeometryTelescopeGeoDescription::getXYZfromArcLength()" << std::endl;
                
  // Get magnetic field vector
  gear::Vector3D vectorGlobal( pos[0], pos[1], pos[2] );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
  const double bx         = B.at( vectorGlobal ).x();
  const double by         = B.at( vectorGlobal ).y();
  const double bz         = B.at( vectorGlobal ).z();
  TVector3 hVec(bx,by,bz);
               
	const double H = hVec.Mag();
  const double p = pVec.Mag();
	const double constant = -0.299792458; 
	const double mm = 1000;
	const double combineConstantsAndMagneticField = (constant*beamQ*H)/mm;
	const double k = combineConstantsAndMagneticField;
  const double rho = combineConstantsAndMagneticField/p; 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////       

	if ( fabs( k ) > 0  ){
		// Non-zero magnetic field case
		TVector3 pCrossH = pVec.Cross(hVec.Unit());
		TVector3 pCrossHCrossH = pCrossH.Cross(hVec.Unit());
		const double pDotH = pVec.Dot(hVec.Unit());
		TVector3 temp1 = pCrossHCrossH;	temp1 *= ( (-1./k) * sin( rho * s ) );
		TVector3 temp2 = pCrossH;       temp2 *= ( (-1./k) * ( 1. - cos( rho * s ) ) );
		TVector3 temp3 = hVec;          temp3 *= ( (pDotH / p) * s );
		pos += temp1;
		pos += temp2;
		pos += temp3;
	}else{
		// Vanishing magnetic field case. //Here you just determine the fraction of P in each direction and this must be the fraction of s that this direction gets
		const double cosA =  pVec[0]/p;      // Calculate cos of the angle between Z(beam) and X(solenoid field axis) //NEED TO MAKE SURE THAT TX=PX/P
		const double cosB = pVec[1]/p ;        // Calculate cos of the angle between Z(beam) and Y
		pos.SetX( pos[0] + cosA * s );
		pos.SetY( pos[1] + cosB * s );
		pos.SetZ( pos[2] + 1./p * pVec.Z() * s );
	}
        
	streamlog_out(DEBUG2) << "---------------------------------EUTelKalmanFilter::getXYZfromArcLength()------------------------------------" << std::endl;
        
	return pos;
}
//Here we define the transformation between the curvilinear and the local frame. Note the local frame we have is defined as the local frame of the telescope. 
//We do this since the local frame is arbitrary.
//However we must describe the local frame in the curvilinear frame. Since we connect the two frames via the same global frame. 
//The transformations are the same as described below.
//This is unlike the curvilinear frame which can not have the particles moving in z due to construction. 
//Therefore we follow the paper: Derivation of Jacobians for the propagation of covariance matrices of track parameters in homogeneous magnetic fields. 
//Which implies that anytime we are describing the curvilinear system using our global momentums then we need to transform the from our system to theirs.
//Within the function this will be clearly labelled
//This is a simple transform our x becomes their(curvilinear y), our y becomes their z and z becomes x
//However this is ok since we never directly access the curvilinear system. It is only a bridge between two local systems. 
TMatrixD EUTelGeometryTelescopeGeoDescription::getLocalToCurvilinearTransformMatrix(TVector3 globalMomentum, int  planeID, float charge){
	const gear::BField&   Bfield = geo::gGeometry().getMagneticFiled();
	gear::Vector3D vectorGlobal(0.1,0.1,0.1);//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
	//Magnetic field must be changed to curvilinear coordinate system. Since this is used in the curvilinear jacobian/////////////////////////////////////////////////////////////////////////////////////////
	const double Bx = (Bfield.at( vectorGlobal ).z())*0.3;//We times bu 0.3 due to units of other variables. See paper. Must be Tesla
	const double By = (Bfield.at( vectorGlobal ).y())*0.3;
	const double Bz = (Bfield.at( vectorGlobal ).x())*0.3;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TVector3 B;
	B[0]=Bx; B[1]=By; B[2]=Bz;
	TVector3 H = B.Unit();
	////////////////////////////////////We transform the momentum to curvilinear frame. Since this is what is used to describe the curvilinear frame///////////////////// 
	TVector3 curvilinearGlobalMomentum;
	curvilinearGlobalMomentum[0]=globalMomentum[2];curvilinearGlobalMomentum[1]=globalMomentum[1];curvilinearGlobalMomentum[2]=globalMomentum[0];
	TVector3 T = curvilinearGlobalMomentum.Unit();//With no magnetic field this will point in z direction.	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	const float cosLambda = sqrt(T[0] * T[0] + (T[1]) * (T[1]));
	//////////////////////////////////////////We use the global curvilinear frame z direction to create the other unit axis/////////////// 
	TVector3 zGlobalNormal;
	zGlobalNormal[0]=0;zGlobalNormal[1]=0;zGlobalNormal[2]=1;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TVector3 U = (zGlobalNormal.Cross(T)).Unit(); 
	TVector3 V = (T.Cross(U)); 	
	streamlog_out(DEBUG0)<<"The Z(V) axis of the curvilinear system in the global system (Note not the same as global telescope frame)."<< std::endl; 
	streamlog_message( DEBUG0, V.Print();, std::endl; );
	streamlog_out(DEBUG0)<<"The Y(U) axis of the curvilinear system in the global system (Note not the same as global telescope frame)"<< std::endl; 
	streamlog_message( DEBUG0, U.Print();, std::endl; );
	streamlog_out(DEBUG0)<<"The X(T) axis of the curvilinear system in the global system  (Note not the same as global telescope frame). Should be beam direction."<< std::endl; 
	streamlog_message( DEBUG0, T.Print();, std::endl; );
	TVector3 ITelescopeFrame;
	TVector3 KTelescopeFrame;
	TVector3 JTelescopeFrame;
	TVector3 I;
	TVector3 K;
	TVector3 J;

	///This is the EUTelescope local z direction.//////////////////////////// 
	streamlog_out(DEBUG0)<<"Set Z(I) with correct sensor rotation. Plane: "<< planeID  << std::endl; 
	ITelescopeFrame = geo::gGeometry().siPlaneNormal( planeID  );       
	I[0]=ITelescopeFrame[2];I[1]=ITelescopeFrame[1];I[2]=ITelescopeFrame[0];
	//////////////////////////////////////////////////////////////////////////////
	//This is the EUTelescope local Y direction////////////////////////////////////////
	streamlog_out(DEBUG0)<<"Set Y(K) with correct sensor rotation. Plane: "<< planeID  << std::endl; 
	KTelescopeFrame = geo::gGeometry().siPlaneYAxis( planeID  ); //This is the y direction of the local frame in global coordinates.      
	K[0]=KTelescopeFrame[2];K[1]=KTelescopeFrame[1];K[2]=KTelescopeFrame[0];
	//This is the EUTelescope local x direction//////////////////////////////////
	streamlog_out(DEBUG0)<<"Set X(J) with correct sensor rotation. Plane: "<< planeID  << std::endl; 
	JTelescopeFrame  = geo::gGeometry().siPlaneXAxis( planeID  ); //X direction      
	J[0]=JTelescopeFrame[2];J[1]=JTelescopeFrame[1];J[2]=JTelescopeFrame[0];
	///////////////////////////////////////////////////////////////////
	streamlog_out(DEBUG0)<<"The Z(J) axis of local system in the global curvilinear system"<< std::endl; 
	streamlog_message( DEBUG0, J.Print();, std::endl; );
	streamlog_out(DEBUG0)<<"The Y(K) axis of local system in the global curvilinear system"<< std::endl; 
	streamlog_message( DEBUG0, K.Print();, std::endl; );
	streamlog_out(DEBUG0)<<"The X(I) axis of local system in the global curvilinear system"<< std::endl; 
	streamlog_message( DEBUG0, I.Print();, std::endl; );

	TVector3 N = (H.Cross(T)).Unit();
	//
	const double alpha = (H.Cross(T)).Mag();
	const double Q = -(B.Mag())*(charge/(curvilinearGlobalMomentum.Mag()));//You could use end momentum since it must be constant
	//
	const double TDotI = T.Dot(I);
	const double TDotJ = T.Dot(J);
	const double TDotK = T.Dot(K);
	const double VDotJ = V.Dot(J);
	const double VDotK = V.Dot(K);
	const double VDotN = V.Dot(N);
	const double UDotJ = U.Dot(J);
	const double UDotK = U.Dot(K);
	const double UDotN = U.Dot(N);
	TMatrix jacobian(5,5);
	jacobian.Zero();
	jacobian[0][0]=1; 
	                 jacobian[1][1]=TDotI*VDotJ;             jacobian[1][2]=TDotI*VDotK;             jacobian[1][3]=-alpha*Q*TDotJ*VDotN;             jacobian[1][4]=-alpha*Q*TDotK*VDotN;
	                 jacobian[2][1]=(TDotI*UDotJ)/cosLambda; jacobian[2][2]=(TDotI*UDotK)/cosLambda; jacobian[2][3]=(-alpha*Q*TDotJ*UDotN)/cosLambda; jacobian[2][4]=(-alpha*Q*TDotK*UDotN)/cosLambda;
																																																   jacobian[3][3]=UDotJ;													  jacobian[3][4]=UDotK;
																																																   jacobian[4][3]=VDotJ;													  jacobian[4][4]=VDotK;
	return jacobian;
}
//This is described in Derivations of Jacobians for the propagation of covariance matrices of track parameters in homogeneous magnetic fields. A satrandie, W Wittek
//This papaer describes the one letter variables. 
//s must be in metres
//momentum must be GeV/c
/*
TMatrix EUTelGeometryTelescopeGeoDescription::getPropagationJacobianCurvilinear(TVector3 globalMomentumStart, TVector3 globalMomentumEnd, TVector3 globalPositionStart, TVector3 globalPositionEnd, float charge, float arcLength){
	const float lambda0 = asin(globalMomentumStart[2]/(globalMomentumStart.Mag()));//This will be in radians.
	const float lambda = asin(globalMomentumEnd[2]/(globalMomentumEnd.Mag()));//This will be in radians.
	const float phi0 = asin(globalMomentumStart[1]/(globalMomentumStart.Mag())*cos(lambda0));
	const float phi = asin(globalMomentumEnd[1]/(globalMomentumEnd.Mag())*cos(lambda));
	const float cosLambda0 = cos(lambda0);
	const float cosLambda = cos(lambda);

	//The global normal will alway be in this direction.
	TVector3 zGlobalNormal;
	zGlobalNormal[0]=0;zGlobalNormal[1]=0;zGlobalNormal[2]=1;
	TVector3 M0;
	M0[0] = globalPositionStart[0]; M0[1] = globalPositionStart[1]; M0[2] = globalPositionStart[2];
	TVector3 M;
	M[0] = globalPositionEnd[0]; M[1] = globalPositionEnd[1]; M[2] = globalPositionEnd[2];
	TVector3 M0MinusM = M0-M;
	TVector3 T0 = globalMomentumStart.Unit();//This is one coordinate axis of curvilinear coordinate system.	
	TVector3 U0 = (zGlobalNormal.Cross(T0)).Unit();//This is the next coordinate axis
	TVector3 V0 = (T0.Cross(U0));	 
	TVector3 T = globalMomentumEnd.Unit();//This is one coordinate axis of curvilinear coordinate system.	
	TVector3 U = (zGlobalNormal.Cross(T)).Unit();//This is the next coordinate axis
	TVector3 V = (T.Cross(U));	
	const gear::BField&   Bfield = geo::gGeometry().getMagneticFiled();
	gear::Vector3D vectorGlobal(globalPositionStart[0],globalPositionStart[1],globalPositionStart[1]);//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
	const double Bx = (Bfield.at( vectorGlobal ).x())*0.3;//We times bu 0.3 due to units of other variables. See paper. Must be Tesla
	const double By = (Bfield.at( vectorGlobal ).y())*0.3;
	const double Bz = (Bfield.at( vectorGlobal ).z())*0.3;
	TVector3 B;
	B[0]=Bx; B[1]=By; B[2]=Bz;
	TVector3 H = B.Unit();
	const float alpha = (H.Cross(T)).Mag();
	const float gamma = H.Dot(T);
	const float Q = -(B.Mag())*(charge/(globalMomentumStart.Mag()));//You could use end momentum since it must be constant
	const float qpInv = pow((charge/(globalMomentumStart.Mag())),-1);
	float theta = Q*arcLength;
	TVector3 N = (H.Cross(T)).Unit();
	//Some terms used in the deriviatives.
	const float NDotV = N.Dot(V);
	const float NDotU = N.Dot(U);
	const float VDotM0MinusM = V.Dot(M0MinusM);
	const float V0DotV = V0.Dot(V); 
	const float V0DotU = V0.Dot(U);
	const float V0DotT = V0.Dot(T);
	const float V0DotN = V0.Dot(N);
	const float HDotV0 = H.Dot(V0);
	const float HDotV = H.Dot(V);
	const float HDotT = H.Dot(T);
	const float HDotU = H.Dot(U);
	const float HDotU0 = H.Dot(U0);
	const float HCrossV0DotV = (H.Cross(V0)).Dot(V);
	const float HCrossU0DotV = (H.Cross(U0)).Dot(V);
	const float HCrossV0DotU = (H.Cross(V0)).Dot(U);
	const float HCrossU0DotU = (H.Cross(U0)).Dot(U);
	const float U0DotV = U0.Dot(V);
	const float U0DotU = U0.Dot(U);
	const float U0DotT = U0.Dot(T);
	const float U0DotN = U0.Dot(N);
	const float UDotM0MinusM = U.Dot(M0MinusM);
	const float cosTheta = cos(theta);
	const float oneMinusCosTheta = (1-cos(theta));
	const float sinTheta =sin(theta);
	const float thetaMinusSinTheta = (theta - sin(theta));
	//Here is the derivatives for the jacobain
	//mometum
	const float dQPdQP0 = 1; //A.1  //A.1 is the equation number in the paper
	//lambda
	const float dLambdadQP0= -alpha*Q*qpInv*(NDotV)*(T.Dot(M0MinusM));//A.2
	const float dLambdadLambda0 =  cosTheta*V0DotV + sinTheta*HCrossV0DotV + oneMinusCosTheta*HDotV0*HDotV + alpha*NDotV*(-sinTheta*V0DotT) + alpha*(oneMinusCosTheta)*V0DotN - thetaMinusSinTheta*HDotT*HDotV0;//A.3 
	const float dLambdadPhi0 = cosLambda0*(cosTheta*U0DotV+sinTheta*HCrossU0DotV + oneMinusCosTheta*HDotU0*HDotV +alpha*NDotV*(-sinTheta*U0DotT + alpha*oneMinusCosTheta*U0DotN - thetaMinusSinTheta*HDotT*HDotU0));//A.4  
	const float dLambdadx0 = -alpha*Q*NDotV*U0DotT;//A.5
	const float dLambdady0 = -alpha*Q*NDotV*V0DotT;//A.6
	//phi
	const float dPhidQP0 = -((alpha*Q)/cosLambda)*qpInv*NDotU*(T.Dot(M0MinusM));//A.7 
	const float dPhidLambda0 = (1/cosLambda)*(cosTheta*V0DotU+ sinTheta*HCrossV0DotU+oneMinusCosTheta*HDotV0*HDotU +alpha*NDotU*(-sinTheta*V0DotT + alpha*oneMinusCosTheta*V0DotN -thetaMinusSinTheta*HDotT*HDotV0));//A.8
	const float dPhidPhi0 = (cosLambda0/cosLambda)*(cosTheta*U0DotU + sinTheta*HCrossU0DotU+oneMinusCosTheta*HDotU0*HDotU +alpha*NDotU*(-sinTheta*U0DotT+alpha*oneMinusCosTheta*U0DotN-thetaMinusSinTheta*HDotT*HDotU0));//A.9 
	const float dPhidx0 = -(alpha*Q/cosLambda)*NDotU*U0DotT;//A.10
	const float dPhidy0 = -(alpha*Q/cosLambda)*NDotU*V0DotT;//A.11
	//x
	const float dxdQP0 = qpInv*UDotM0MinusM;//A.12
	const float dxdLambda0 = (sinTheta/Q)*V0DotU + (oneMinusCosTheta/Q)*HCrossV0DotU+(thetaMinusSinTheta/Q)*HDotV0*HDotU;//A.13
	const float dxdPhi0 = cosLambda0*((sinTheta/Q)*U0DotU + (oneMinusCosTheta/Q)*HCrossU0DotU+(thetaMinusSinTheta/Q)*HDotU0*HDotU);//A.14
	const float dxdx0 = U0DotU;//A.15
	const float dxdy0 = V0DotU;//A.16
	//y
	const float dydQP0 = qpInv*VDotM0MinusM;//A.17
	const float dydLambda0 = (sinTheta/Q)*V0DotV + (oneMinusCosTheta/Q)*HCrossV0DotV+(thetaMinusSinTheta/Q)*HDotV0*HDotV;//A.18 
	const float dydPhi0 = cosLambda0*((sinTheta/Q)*U0DotV +(oneMinusCosTheta/Q)*HCrossU0DotV+(thetaMinusSinTheta/Q)*HDotU0*HDotV);//A.19
	const float dydx0 = U0DotV;//A.20
	const float dydy0 = V0DotV;//A.21
	//Create matrix from these terms
	TMatrix jacobian(5,5);
	jacobian.Zero();
	jacobian[0][0] = dQPdQP0; 
	jacobian[1][0] = dLambdadQP0; jacobian[1][1]=dLambdadLambda0; jacobian[1][2]=dLambdadPhi0; jacobian[1][3]=dLambdadx0; jacobian[1][4]=dLambdady0;
	jacobian[2][0] = dPhidQP0; jacobian[2][1]=dPhidLambda0; jacobian[2][2]=dPhidPhi0; jacobian[2][3]=dPhidx0; jacobian[2][4]=dPhidy0;
	jacobian[3][0] = dxdQP0; jacobian[3][1]=dxdx0; jacobian[3][2]=dxdPhi0; jacobian[3][3]=dxdx0; jacobian[3][4]=dxdy0;
	jacobian[4][0] = dydQP0; jacobian[4][1]=dydy0; jacobian[4][2]=dydPhi0; jacobian[4][3]=dydx0; jacobian[4][4]=dydy0;

	return jacobian;

}
*/
/*
 * \param [in] ds   (3D) arc length to end point
 *  * \param [in] qbyp q/p (signed inverse momentum)
 *   * \param [in] t1   track direction at start point
 *    * \param [in] t2   track direction at end point
 *     * \param [in] b   B (magnetic field)
 *      * \return (5*5) propagation matrix
 *       */
//Note that the curvilinear frame that this jacobian has been derived in does not work for particles moving in z-direction. 
//This is due to a construction that assumes tha beam pipe is in the z-direction. Since it is used for collider experiements. 
//We therefore have to change to the coordinate system used in paper before we apply this jacobian.
//This is ok since we never access the curvilinear system directly but always through the local system which is defined in the local frame of the telescope
//I.e Telescope x becomes y, y becomes z and z becomes x.
TMatrixD EUTelGeometryTelescopeGeoDescription::getPropagationJacobianCurvilinear(float ds , float  qbyp,  TVector3 t1w, TVector3 t2w) {
	TVector3 t1(t1w[2],t1w[1],t1w[0]);//This is need to change to claus's coordinate system
	TVector3 t2(t2w[2],t2w[1],t2w[0]);
	t1.Unit();
	t2.Unit();
//	if(t1.Mag() != 1.0 or t2.Mag() != 1.0){//TO DO: This statement does not work for some reason.
//		streamlog_out(MESSAGE9) << "The magnitude of the two vectors is:  "<< t1.Mag()<< " , "<<t2.Mag() << std::endl;
//		throw(lcio::Exception(Utility::outputColourString("The calculation of the jacobian must be performed by unit vectors!.", "RED"))); 	
//	}
	const gear::BField&   Bfield = getMagneticFiled();
	//Must transform the B field to be in the frame used in the curvilinear frame 
	gear::Vector3D vectorGlobal(0.1,0.1,0.1);//Since field is homogeneous this seems silly but we need to specify a position to geometry to get B-field.
	/////////////////////////////////////////////////////////Must also change the magnetic field to be in the correct coordinate system///////////
	const double Bx = (Bfield.at( vectorGlobal ).z())*0.3*pow(10,-3); //We times but 0.3 due to units of other variables. See paper. Must be Tesla
	const double By = (Bfield.at( vectorGlobal ).y())*0.3*pow(10,-3);
	const double Bz = (Bfield.at( vectorGlobal ).x())*0.3*pow(10,-3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TVector3 b;
	b[0]=Bx; b[1]=By; b[2]=Bz;
	streamlog_out( DEBUG2 ) << "EUTelGeometryTelescopeGeoDescription::getPropagationJacobianCurvilinear()------BEGIN" << std::endl;
	streamlog_out( DEBUG2 ) <<"This is the input to the jacobian"<<endl;  
	streamlog_out( DEBUG2 ) <<"The arc length: " <<ds <<endl;
	streamlog_out( DEBUG2 ) <<"The curvature: "<< qbyp <<endl; 
	streamlog_out(DEBUG0)<<"The unit momentum start "<< std::endl; 
	streamlog_message( DEBUG0, t1.Print();, std::endl; );
	streamlog_out(DEBUG0)<<"The unit momentum end "<< std::endl; 
	streamlog_message( DEBUG0, t2.Print();, std::endl; );
	streamlog_out(DEBUG0)<<"The unit Magnetic field  "<< std::endl; 
	streamlog_message( DEBUG0, b.Print();, std::endl; );
	TMatrixD ajac(5, 5);
	TVector3  bc  = b;//This is b*c. speed of light in 1 nanosecond
	ajac.UnitMatrix(); 
	const double qp = -bc.Mag(); // -|B*c|
	const double q = qp * qbyp; // Q
	if (q == 0.) {
		// line
 		ajac[3][2] = ds * sqrt(t1[0] * t1[0] + t1[1] * t1[1]);
		ajac[4][1] = ds;
	} else {
		// helix
		// at start
		const double cosl1 = sqrt(t1[0] * t1[0] + t1[1] * t1[1]);
		// at end
		const double cosl2 = sqrt(t2[0] * t2[0] + t2[1] * t2[1]);
		const double cosl2Inv = 1. / cosl2;
		// magnetic field direction
		TVector3 hn(bc.Unit());
		// (signed) momentum
		const double pav = 1.0 / qbyp;
		//
		const double theta = q * ds;
		const double sint = sin(theta);
		const double cost = cos(theta);
		const double gamma = hn.Dot(t2); // H*T
		TVector3 an1 = hn.Cross(t1); // HxT0
		TVector3 an2 = hn.Cross(t2); // HxT
		// U0, V0
		const double au1 = 1. / sqrt(t1[0]*t1[0]+t1[1]*t1[1]);
		TVector3 u1(-au1 * t1[1], au1 * t1[0], 0.);
		TVector3 v1(-t1[2] * u1[1], t1[2] * u1[0], t1[0] * u1[1] - t1[1] * u1[0]);
		// U, V
		const double au2 = 1. /sqrt(t2[0]*t2[0]+t2[1]*t2[1]);
		TVector3 u2(-au2 * t2[1], au2 * t2[0], 0.);
		TVector3 v2(-t2[2] * u2[1], t2[2] * u2[0], t2[0] * u2[1] - t2[1] * u2[0]);
		//
		const double anv = -hn.Dot(u2); // N*V=-H*U
		const double anu = hn.Dot(v2);  // N*U= H*V
		const double omcost = 1. - cost;
		const double tmsint = theta - sint;
		// M0-M
		TVector3 dx(-(gamma * tmsint * hn[0] + sint * t1[0] + omcost * an1[0]) / q,
		-(gamma * tmsint * hn[1] + sint * t1[1] + omcost * an1[1]) / q,
		-(gamma * tmsint * hn[2] + sint * t1[2] + omcost * an1[2]) / q);
		// HxU0
		TVector3 hu1 = hn.Cross(u1);
		// HxV0
		TVector3 hv1 = hn.Cross(v1);
		// some.Dot products
		const double u1u2 = u1.Dot(u2), u1v2 = u1.Dot(v2), v1u2 = v1.Dot(u2), v1v2 = v1.Dot(v2);
		const double hu1u2 = hu1.Dot(u2), hu1v2 = hu1.Dot(v2), hv1u2 = hv1.Dot(u2), hv1v2 = hv1.Dot(v2);
		const double hnu1 = hn.Dot(u1), hnv1 = hn.Dot(v1), hnu2 = hn.Dot(u2), hnv2 = hn.Dot(v2);
		const double t2u1 = t2.Dot(u1), t2v1 = t2.Dot(v1);
		const double t2dx = t2.Dot(dx), u2dx = u2.Dot(dx), v2dx = v2.Dot(dx);
		const double an2u1 = an2.Dot(u1), an2v1 = an2.Dot(v1);
		// jacobian
		// 1/P
		ajac[0][0] = 1.;
		// Lambda
		ajac[1][0] = -qp * anv * t2dx;
		ajac[1][1] = cost * v1v2 + sint * hv1v2 + omcost * hnv1 * hnv2 + anv * (-sint * t2v1 + omcost * an2v1 - gamma * tmsint * hnv1);
		ajac[1][2] = cosl1
		* (cost * u1v2 + sint * hu1v2 + omcost * hnu1 * hnv2 + anv * (-sint * t2u1 + omcost * an2u1 - gamma * tmsint * hnu1));
		ajac[1][3] = -q * anv * t2u1;
		ajac[1][4] = -q * anv * t2v1;
		// Phi
		ajac[2][0] = -qp * anu * t2dx * cosl2Inv;
		ajac[2][1] = cosl2Inv
		* (cost * v1u2 + sint * hv1u2 + omcost * hnv1 * hnu2 + anu * (-sint * t2v1 + omcost * an2v1 - gamma * tmsint * hnv1));
		ajac[2][2] = cosl2Inv * cosl1
		* (cost * u1u2 + sint * hu1u2 + omcost * hnu1 * hnu2 + anu * (-sint * t2u1 + omcost * an2u1 - gamma * tmsint * hnu1));
		ajac[2][3] = -q * anu * t2u1 * cosl2Inv;
		ajac[2][4] = -q * anu * t2v1 * cosl2Inv;
		// Xt
		ajac[3][0] = pav * u2dx;
		ajac[3][1] = (sint * v1u2 + omcost * hv1u2 + tmsint * hnu2 * hnv1) / q;
		ajac[3][2] = (sint * u1u2 + omcost * hu1u2 + tmsint * hnu2 * hnu1) * cosl1 / q;
		ajac[3][3] = u1u2;
		ajac[3][4] = v1u2;
		// Yt
		ajac[4][0] = pav * v2dx;
		ajac[4][1] = (sint * v1v2 + omcost * hv1v2 + tmsint * hnv2 * hnv1) / q;
		ajac[4][2] = (sint * u1v2 + omcost * hu1v2 + tmsint * hnv2 * hnu1) * cosl1 / q;
		ajac[4][3] = u1v2;
		ajac[4][4] = v1v2;
	}
	streamlog_out( DEBUG2 ) << "EUTelGeometryTelescopeGeoDescription::getPropagationJacobianCurvilinear()------END" << std::endl;

		return ajac;
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
  TMatrix EUTelGeometryTelescopeGeoDescription::getPropagationJacobianF( float x0, float y0, float z0, float px, float py, float pz, float beamQ, float dz ) {
        streamlog_out( DEBUG2 ) << "EUTelGeometryTelescopeGeoDescription::getPropagationJacobianF()" << std::endl;
	// The formulas below are derived from equations of motion of the particle in
	// magnetic field under assumption |dz| small. Must be valid for |dz| < 10 cm
	const double k = 0.299792458*pow(10,-4);//Here is the conversion from k=(GeV/c) KG mm^-1
	const double pxMomentum = px;//change momentum to GeV/c 
	const double pyMomentum = py;//change momentum to GeV/c 
	const double pzMomentum = pz;//change momentum to GeV/c 
	TVector3 pVec(pxMomentum, pyMomentum, pzMomentum);	
	// Get track parameters
	const double invP = beamQ/pVec.Mag();//This should be in 1/(GeV/c)
	const double tx0 = px/pVec.Mag();//Unitless  
	const double ty0 = py/pVec.Mag();

        // Get magnetic field vector
	gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
	const double Bx         = (B.at( vectorGlobal ).x())*10;//times but 10 to convert from Tesla to KiloGauss. 1 T = 10^4 Gauss.
	const double By         = (B.at( vectorGlobal ).y())*10;
	const double Bz         = (B.at( vectorGlobal ).z())*10;
	
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
	jacobianF[3][1] = dxdtx0;	jacobianF[3][2] = dxdty0;	jacobianF[3][0] = dxdinvP0;
	jacobianF[4][1] = dydtx0;	jacobianF[4][2] = dydty0;	jacobianF[4][0] = dydinvP0;
	jacobianF[1][2] = dtxdty0;	jacobianF[1][0] = dtxdinvP0;
	jacobianF[2][1] = dtydtx0;	jacobianF[2][0] = dtydinvP0;
        
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

