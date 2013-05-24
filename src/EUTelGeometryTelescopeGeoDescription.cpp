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

// ROOT
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoBBox.h"
#include "TVectorD.h"
#include "TMath.h"

using namespace eutelescope;

using namespace std;

EUTelGeometryTelescopeGeoDescription& EUTelGeometryTelescopeGeoDescription::getInstance() {
    static EUTelGeometryTelescopeGeoDescription instance;
    return instance;
}

EUTelGeometryTelescopeGeoDescription::EUTelGeometryTelescopeGeoDescription() {

    // -------          Сopy-paste from another class           ----------- //

    // Check if the GEAR manager is not corrupted, otherwise stop

    if (!marlin::Global::GEAR) {
        streamlog_out(ERROR2) << "The GearMgr is not available, for an unknown reason." << std::endl;
        throw eutelescope::InvalidGeometryException("GEAR manager is not initialised");
    }

    // sensor-planes in geometry navigation:
    _siPlanesParameters = const_cast<gear::SiPlanesParameters*> (&(marlin::Global::GEAR->getSiPlanesParameters()));
    _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));

    // create an array with the z positions of each layer
    for (int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++) {
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getLayerPositionZ(iPlane));
        _sensorIDMap.insert(std::make_pair(_siPlanesLayerLayout->getLayerPositionZ(iPlane), this->_siPlanesLayerLayout->getID(iPlane)));
    }

    if (_siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT) {
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getDUTPositionZ());
        _sensorIDMap.insert(std::make_pair(_siPlanesLayerLayout->getDUTPositionZ(), _siPlanesLayerLayout->getDUTID()));
    }

    // sort the array with increasing z
    std::sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());

    // clear the sensor ID vector
    _sensorIDVec.clear();

    // clear the sensor ID map
    _sensorIDVecMap.clear();
    _sensorIDtoZOrderMap.clear();

    // clear the sensor ID vector (z-axis order)
    _sensorIDVecZOrder.clear();

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

        _sensorIDVecZOrder.push_back(sensorsToTheLeft);
        _sensorIDtoZOrderMap.insert(std::make_pair(sensorID, sensorsToTheLeft));
    }

    delete [] keepZPosition;

    _nPlanes = _siPlanesParameters->getSiPlanesNumber();
    if (_siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT) ++_nPlanes;
    
    
    
    // TGeo manager initialisation
    
}

EUTelGeometryTelescopeGeoDescription::EUTelGeometryTelescopeGeoDescription(const EUTelGeometryTelescopeGeoDescription& orig) {
}

EUTelGeometryTelescopeGeoDescription::~EUTelGeometryTelescopeGeoDescription() {
    delete _geoManager;
}

void EUTelGeometryTelescopeGeoDescription::initializeTGeoDescription( string tgeofilename ) {
//    #ifdef USE_TGEO
    // get access to ROOT's geometry manager
    
    _geoManager = TGeoManager::Import( tgeofilename.c_str() );
    if( !_geoManager ) {
        streamlog_out( WARNING ) << "Can't read file " << tgeofilename << endl;
    }
//    _geoManager->CloseGeometry();
//    #endif //USE_TGEO
}

/** From ROOT's geometry stress test */
//void EUTelGeometryTelescopeGeoDescription::findRad(Double_t x, Double_t y, Double_t z,
//        Double_t theta, Double_t phi, Int_t &nbound, Float_t &length, Float_t &safe, Float_t &rad, Bool_t verbose) {
//    
//   Double_t xp  = TMath::Sin(theta)*TMath::Cos(phi);
//   Double_t yp  = TMath::Sin(theta)*TMath::Sin(phi);
//   Double_t zp  = TMath::Cos(theta);
//   Double_t snext;
//   char path[256];
//   Double_t pt[3];
//   Double_t loc[3];
//   Double_t epsil = 1.E-2;
//   Double_t lastrad = 0.;
//   Int_t ismall = 0;
//   nbound = 0;
//   length = 0.;
//   safe   = 0.;
//   rad    = 0.;
//   TGeoMedium *med;
//   TGeoShape *shape;
//   gGeoManager->InitTrack(x,y,z,xp,yp,zp);
//             
//   TGeoNode *nextnode = gGeoManager->GetCurrentNode();
//   safe = gGeoManager->Safety();
//   while (nextnode) {
//      med = 0;
//      if (nextnode) med = nextnode->GetVolume()->GetMedium();
//      else return;      
//      shape = nextnode->GetVolume()->GetShape();
//      nextnode = gGeoManager->FindNextBoundaryAndStep();
//      snext  = gGeoManager->GetStep();
//      if (snext<1.e-8) {
//         ismall++;
//         if (ismall > 3) {
//            streamlog_out << "ERROR: Small steps in: " << gGeoManager->GetPath() << " shape=" << shape->ClassName() << endl;
//            return;
//         }   
//         memcpy(pt,gGeoManager->GetCurrentPoint(),3*sizeof(Double_t));
//         const Double_t *dir = gGeoManager->GetCurrentDirection();
//         for (Int_t i=0;i<3;i++) pt[i] += epsil*dir[i];
//         snext = epsil;
//         length += snext;
//         rad += lastrad*snext;
//         gGeoManager->CdTop();
//         nextnode = gGeoManager->FindNode(pt[0],pt[1],pt[2]);
//         if (gGeoManager->IsOutside()) return;
//         TGeoMatrix *mat = gGeoManager->GetCurrentMatrix();
//         mat->MasterToLocal(pt,loc);
//         if (!gGeoManager->GetCurrentVolume()->Contains(loc)) {
//            gGeoManager->CdUp();
//            nextnode = gGeoManager->GetCurrentNode();
//         }   
//         continue;
//      } else {
//         ismall = 0;
//      }      
//      nbound++;
//      length += snext;
//      if (med) {
//         Double_t radlen = med->GetMaterial()->GetRadLen();
//         if (radlen>1.e-5 && radlen<1.e10) {
//            lastrad = med->GetMaterial()->GetDensity()/radlen;
//            rad += lastrad*snext;
//         } else {
//            lastrad = 0.;
//         }      
//             streamlog_out(DEBUG0) << "STEP #" << nbound << " " << path << endl;
//             streamlog_out(DEBUG0) << "   step=" << snext << " length=" << length 
//                           << " rad=" << med->GetMaterial()->GetDensity()*snext/med->GetMaterial()->GetRadLen()
//                           << " " << med->GetName() << endl;
//             streamlog_out(DEBUG0) << gGeoManager->GetPath() << endl;
//      }
//   }   
//}