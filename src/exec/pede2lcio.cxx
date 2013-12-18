//
#include "marlin/VerbosityLevels.h"

// eutelescope includes
#include "EUTelAlignmentConstant.h"
#include "anyoption.h"

// lcio includes
#include <IO/LCWriter.h>
#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCCollectionVec.h>

// GEAR
#include "gear/GearMgr.h"
#include "gear/gearxml/GearXML.h"
#include "gearimpl/Util.h"
#include "gear/SiPlanesLayerLayout.h"
#include "gear/SiPlanesParameters.h"

// system includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <functional>
#include <algorithm>
#include <memory>

#include "TRotation.h"
#include "TVector3.h"

using namespace std;

struct CollectionWriter {
    lcio::LCCollectionVec* _constantsCollection;
        void operator()( std::pair<const int, eutelescope::EUTelAlignmentConstant*>& pair ) {
          _constantsCollection->push_back( pair.second );
                cout << (*pair.second) << endl;
        }
} colWriter;

void prepareGEAR( const string& oldGearfileName, const string& newGearfileName, const map< int, eutelescope::EUTelAlignmentConstant* >& alignmentConstants ) {
    
    streamlog_out(MESSAGE4) << "Reading " << oldGearfileName << std::endl;
    streamlog_out(MESSAGE4) << "GEAR file " << newGearfileName << " will be generated." << std::endl;
    
    gear::GearXML gearXML( oldGearfileName ) ;

    gear::GearMgr* gearManager = gearXML.createGearMgr() ;
    
    if (!gearManager) {
        cerr << "Cannot instantiate GEAR manager" << std::endl;
        return;
    }

    // sensor-planes in geometry navigation:
    gear::SiPlanesParameters* siPlanesParameters = const_cast<gear::SiPlanesParameters*> (&(gearManager->getSiPlanesParameters()));
    gear::SiPlanesLayerLayout* siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> (&(siPlanesParameters->getSiPlanesLayerLayout()));

    // update positions and orientations of the planes
    // TODO: set appropriate new values for new GEAR file

    streamlog_out(MESSAGE4) << "Combined alignment (current GEAR + MILLE corrections below):" << std::endl;    

    map< int, eutelescope::EUTelAlignmentConstant* >::const_iterator itrAlignmentConstant;

    streamlog_out(MESSAGE4) << "Plane ID" << setw(20) << "X shift" << setw(20) << "Y shift" << setw(20) << "Z shift" << setw(20)
                     << "X ROTATION" << setw(20) << "Y ROTATION" << setw(20) << "Z ROTATION" << std::endl;
    
    for (int iPlane = 0; iPlane < siPlanesLayerLayout->getNLayers(); iPlane++) {
        int sensorID = siPlanesLayerLayout->getSensitiveID(iPlane);
        if( ( itrAlignmentConstant = alignmentConstants.find( sensorID ) ) != alignmentConstants.end() ) {


 	    const double alpha = siPlanesLayerLayout->getLayerRotationZY(iPlane);
	    const double beta  = siPlanesLayerLayout->getLayerRotationZX(iPlane);
	    const double gamma = siPlanesLayerLayout->getLayerRotationXY(iPlane);

            streamlog_out(MESSAGE4) << "former " << sensorID << setw(20) << siPlanesLayerLayout->getLayerPositionX(iPlane)  << 
			 	 setw(20) << siPlanesLayerLayout->getLayerPositionY(iPlane)  <<
				 setw(20) << siPlanesLayerLayout->getLayerPositionZ(iPlane)  <<
                 		 setw(20) << alpha  <<
                                 setw(20) << beta   <<
                                 setw(20) << gamma  << std::endl;

	    TRotation invR;
	    invR.RotateX(alpha);
	    invR.RotateY(beta);
	    invR.RotateZ(gamma);
	    invR.Invert();
	
	    const double dalpha = (*itrAlignmentConstant).second->getAlpha();
            const double dbeta  = (*itrAlignmentConstant).second->getBeta();
            const double dgamma = (*itrAlignmentConstant).second->getGamma();

	    TRotation invDeltaR;
	    invDeltaR.RotateX((*itrAlignmentConstant).second->getAlpha());
            invDeltaR.RotateY((*itrAlignmentConstant).second->getBeta());
            invDeltaR.RotateZ((*itrAlignmentConstant).second->getGamma());
	    invDeltaR.Invert();

	    const double dr0x = (*itrAlignmentConstant).second->getXOffset();
	    const double dr0y = (*itrAlignmentConstant).second->getYOffset();
	    const double dr0z = (*itrAlignmentConstant).second->getZOffset();

	    TVector3 delta_r0( dr0x, dr0y, dr0z );
//	    delta_r0 *= invDeltaR;
//	    delta_r0 *= invR;
//	    delta_r0 = invR*(invDeltaR*delta_r0);
            delta_r0 = invR*delta_r0;
 
             streamlog_out(MESSAGE4) << "invR:"<< std::endl;
             streamlog_out(MESSAGE4) << " X: " << setw(20) << invR[0][0] << " " << invR[0][1] << " " << invR[0][2]  << std::endl;
             streamlog_out(MESSAGE4) << " Y: " << setw(20) << invR[1][0] << " " << invR[1][1] << " " << invR[1][2]  << std::endl;
             streamlog_out(MESSAGE4) << " Z: " << setw(20) << invR[2][0] << " " << invR[2][1] << " " << invR[2][2]  << std::endl;
             streamlog_out(MESSAGE4) << ""<< std::endl;
 
           
//	    delta_r0 *= invR;

#ifdef GEAR_MAJOR_VERSION 
#if GEAR_VERSION_GE( 17,4)  
// ZY and ZX rotations are calculated wrongly yet, do not implement:
// XYZ shifts and XY rotation seems to be correct
//
            siPlanesLayerLayout-> setLayerPositionX(iPlane, siPlanesLayerLayout->getLayerPositionX(iPlane) +  delta_r0.X() ) ;
            siPlanesLayerLayout-> setLayerPositionY(iPlane, siPlanesLayerLayout->getLayerPositionY(iPlane) +  delta_r0.Y() ) ;
            siPlanesLayerLayout-> setLayerPositionZ(iPlane, siPlanesLayerLayout->getLayerPositionZ(iPlane) +  delta_r0.Z() ) ;
            siPlanesLayerLayout->setLayerRotationZY(iPlane, alpha - dalpha );
            siPlanesLayerLayout->setLayerRotationZX(iPlane, beta  - dbeta  );
            siPlanesLayerLayout->setLayerRotationXY(iPlane, gamma - dgamma );
#endif
#endif       
            streamlog_out(MESSAGE4) << "align by shifts (in local frame) " << std::endl;
            streamlog_out(MESSAGE4) << " by: X' " << setw(20) << dr0x;
            streamlog_out(MESSAGE4) << " by: Y' " << setw(20) << dr0y;
            streamlog_out(MESSAGE4) << " by: Z' " << setw(20) << dr0z << std::endl;

            streamlog_out(MESSAGE4) << "align by rotations (in local frame) " << std::endl;
            streamlog_out(MESSAGE4) << " by: al " << setw(20) << dalpha;
            streamlog_out(MESSAGE4) << " by: be " << setw(20) << dbeta;
            streamlog_out(MESSAGE4) << " by: ga " << setw(20) << dgamma << std::endl;

            streamlog_out(MESSAGE4) << "rotated from local to global: " << std::endl;
            streamlog_out(MESSAGE4) << " by: X" << setw(20) <<  delta_r0.X() ;
            streamlog_out(MESSAGE4) << " by: Y" << setw(20) <<  delta_r0.Y();
            streamlog_out(MESSAGE4) << " by: Z" << setw(20) <<  delta_r0.Z() << std::endl;

            streamlog_out(MESSAGE4) << "new : " << sensorID << setw(20) << siPlanesLayerLayout->getLayerPositionX(iPlane)  << 
				setw(20) << siPlanesLayerLayout->getLayerPositionY(iPlane)  <<
				setw(20) << siPlanesLayerLayout->getLayerPositionZ(iPlane)  <<
		               	setw(20) << siPlanesLayerLayout->getLayerRotationZY(iPlane)  <<
                 	  	setw(20) << siPlanesLayerLayout->getLayerRotationZX(iPlane)  <<
                  		setw(20) << siPlanesLayerLayout->getLayerRotationXY(iPlane)  << std::endl;
        }
    }

    
    streamlog_out(MESSAGE4) << "Not implemented" << std::endl;
    gear::GearXML::createXMLFile( gearManager, newGearfileName );
    
    return;
}

int main( int argc, char ** argv ) {

    
  auto_ptr< AnyOption > option( new AnyOption );

  string usageString = 
    "\n"
    "Eutelescope alignment constant converter\n"
    "\n"
    "This program is used to convert the output of pede into\n"
    "a LCIO and/or GEAR file containing alignment constants.\n"
    "\n"
    "Usage: pede2lcio [options] conversion_file lciofile.slcio [old GEAR file] [new GEAR file]\n\n"
    "-h --help       To print this help\n"
    "-g --gear       To generate GEAR file"
    "\n";

  option->addUsage( usageString.c_str() );
  option->setFlag( "help", 'h');
  option->setFlag( "gear", 'g');

  // process the command line arguments
  option->processCommandArgs( argc, argv );

  if ( option->getFlag('h') || option->getFlag( "help" ) ) {
    option->printUsage();
    return 0;
  }

  if ( !( ( option->getArgc() == 2 ) || ( option->getArgc() == 4 ) ) ) {
    option->printUsage();
    return 0;
  }

  string pedeFileName, lcioFileName;
  pedeFileName = option->getArgv(0);
  lcioFileName = option->getArgv(1);

  // check if the lcio file has the extension
  if ( lcioFileName.rfind( ".slcio", string::npos ) == string::npos ) {
    lcioFileName.append( ".slcio" );
  }

  // check GEAR flag
  string oldGearFileName, newGearFileName;
  bool wantGEAR = false;
  if ( option->getFlag('g') || option->getFlag( "gear" ) ) {
      
      streamlog_out(MESSAGE4) << "Requested generation of GEAR" << std::endl;
      
    wantGEAR = true;
    oldGearFileName = option->getArgv(2);
    if ( oldGearFileName.rfind( ".xml", string::npos ) == string::npos ) {
         oldGearFileName.append( ".xml" );
    }
    newGearFileName = option->getArgv(3);
    if ( lcioFileName.rfind( ".xml", string::npos ) == string::npos ) {
         lcioFileName.append( ".xml" );
    }
    
    streamlog_out(MESSAGE4) << oldGearFileName << " " << lcioFileName << std::endl;
  }
  
  streamlog_out(MESSAGE4) << "Converting " << pedeFileName << " in " << lcioFileName << std::endl;

  // try to open the input file. This should be a text file
  ifstream pedeFile( pedeFileName.c_str(), ios::in );

  map< int, eutelescope::EUTelAlignmentConstant* > constants_map;
  
  if ( pedeFile.fail() ) {

    cerr << "Error opening the " << pedeFileName << std::endl;
    return -1;

  } else {

    // open the LCIO output file
    lcio::LCWriter * lcWriter = lcio::LCFactory::getInstance()->createLCWriter();

    try {
      lcWriter->open( lcioFileName.c_str() , lcio::LCIO::WRITE_NEW );
    } catch ( lcio::IOException& e ) {
      cerr << e.what() << std::endl;
      return -2;
    }

    // write an almost empty run header
    lcio::LCRunHeaderImpl * lcHeader  = new lcio::LCRunHeaderImpl;
    lcHeader->setRunNumber( 0 );
    lcWriter->writeRunHeader(lcHeader);

    delete lcHeader;

    lcio::LCEventImpl * event = new lcio::LCEventImpl;
    event->setRunNumber( 0 );
    event->setEventNumber( 0 );

    lcio::LCTime * now = new lcio::LCTime;
    event->setTimeStamp( now->timeStamp() );
    delete now;

    lcio::LCCollectionVec * constantsCollection = new lcio::LCCollectionVec( lcio::LCIO::LCGENERICOBJECT );

    vector< string > tokens;
    istringstream tokenizer;
    string line;
    string buffer;
    double value = 0.;
    double err = 0.;

    int sensorID = 0;

    while ( getline( pedeFile, line ) ) {

      bool goodLine = false;

        if ( !line.empty() ) goodLine = true;
     
        if( !goodLine ) continue;

        tokens.clear();
        tokenizer.clear();
        tokenizer.str( line );

        while ( tokenizer >> skipws >> buffer ) {
          tokens.push_back( buffer ) ;
        }

        if ( ( tokens.size() == 5 ) || ( tokens.size() == 6 ) ) goodLine = true;

        if( !goodLine ) continue;

        value =  atof(tokens[4].c_str());
        sensorID = atoi(tokens[3].c_str());
        if( tokens.size() == 6 ) err = atof(tokens[5].c_str());
        if( tokens.size() == 5 ) err = 0.;
        
        if( constants_map.find( sensorID ) == constants_map.end() ) {
            eutelescope::EUTelAlignmentConstant * constant = new eutelescope::EUTelAlignmentConstant;
            constants_map[sensorID] = constant;
        }
        
        if( tokens[2].compare("shift") == 0 ) {
            if( tokens[1].compare("X") == 0 ) {
                constants_map[sensorID]->setXOffset( value );
                constants_map[sensorID]->setXOffsetError( err ) ;
            }
            if( tokens[1].compare("Y") == 0 ) {
                constants_map[sensorID]->setYOffset( value );
                constants_map[sensorID]->setYOffsetError( err ) ;
            }
            if( tokens[1].compare("Z") == 0 ) {
                constants_map[sensorID]->setZOffset( value );
                constants_map[sensorID]->setZOffsetError( err ) ;
            }
        }
        if( tokens[2].compare("rotation") == 0 ) {
            if( tokens[1].compare("YZ") == 0 ) {
                constants_map[sensorID]->setAlpha( value );
                constants_map[sensorID]->setAlphaError( err ) ;
            }
            if( tokens[1].compare("XZ") == 0 ) {
                constants_map[sensorID]->setBeta( value );
                constants_map[sensorID]->setBetaError( err ) ;
            }
            if( tokens[1].compare("XY") == 0 ) {
                constants_map[sensorID]->setGamma( value );
                constants_map[sensorID]->setGammaError( err ) ;
            }
        }

        constants_map[sensorID]->setSensorID( sensorID );

    }

    // Process GEAR output if requested
    if ( wantGEAR ) prepareGEAR(oldGearFileName, newGearFileName, constants_map);
    
    streamlog_out(MESSAGE4) << "Alignment corrections:" << std::endl;
    
    colWriter._constantsCollection = constantsCollection;
    for_each( constants_map.begin(), constants_map.end(), colWriter );

    event->addCollection( constantsCollection, "alignment" );
    lcWriter->writeEvent( event );
    delete event;

    lcWriter->close();
  }


  pedeFile.close();


  return 0;
}
