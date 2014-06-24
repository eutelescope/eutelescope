//
#include "marlin/VerbosityLevels.h"

// eutelescope includes
#include "EUTelAlignmentConstant.h"
#include "anyoption.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// lcio includes
#include <IO/LCWriter.h>
#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCCollectionVec.h>

// GEAR
#include "gearimpl/Util.h"
#include "gearxml/GearXML.h"
#include "gear/GearMgr.h"
#include "gear/GEAR.h"


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
using namespace eutelescope;

 

struct CollectionWriter {
    lcio::LCCollectionVec* _constantsCollection;
        void operator()( std::pair<const int, EUTelAlignmentConstant*>& pair ) {
          _constantsCollection->push_back( pair.second );
                cout << (*pair.second) << endl;
        }
} colWriter;

void prepareGEAR( const string& oldGearfileName, const string& newGearfileName, const map< int, EUTelAlignmentConstant* >& alignmentConstants ) {
    
    streamlog_out(MESSAGE4) << "Reading " << oldGearfileName << std::endl;
    streamlog_out(MESSAGE4) << "GEAR file " << newGearfileName << " will be generated." << std::endl;
    
    gear::GearXML gearXML( oldGearfileName ) ;
    gear::GearMgr* gearManager = gearXML.createGearMgr() ;
    
    if (!gearManager) {
        cerr << "Cannot instantiate GEAR manager" << std::endl;
        return;
    }


    // Getting access to geometry description
    std::string name("test.root");
    geo::gGeometry( gearManager ).initializeTGeoDescription(name,false);



    // update positions and orientations of the planes
    // TODO: set appropriate new values for new GEAR file

    streamlog_out(MESSAGE4) << "Combined alignment (current GEAR + MILLE corrections below):" << std::endl;    

    map< int, EUTelAlignmentConstant* >::const_iterator itrAlignmentConstant;

    streamlog_out(MESSAGE4) <<  setw(10) << "" << setw(8) << "Plane ID" << setw(13) << setprecision(4) << "X "       << setw(13) << setprecision(4) << "Y "       << setw(13) << setprecision(4) << "Z "
                                                                        << setw(13) << setprecision(4) << "ZY(Xrot)" << setw(13) << setprecision(4) << "ZX(Yrot)" << setw(13) << setprecision(4) << "XY(Zrot)" 
                                                                        << setw(13) << setprecision(4) << "     in global coordinates, the shifts (X,Y,Z)  " << std::endl;

    double xplane = 0.;
    double yplane = 0.;
    double zplane = 0.;
    double xrot   = 0.;
    double yrot   = 0.;
    double zrot   = 0.;



    std::vector <int> sensorIDsVec = geo::gGeometry().sensorIDsVec();    
    for (int iPlane = 0; iPlane < sensorIDsVec.size(); iPlane++) {

        int sensorID = sensorIDsVec.at(iPlane);
        if( ( itrAlignmentConstant = alignmentConstants.find( sensorID ) ) != alignmentConstants.end() ) {


            xplane = geo::gGeometry().siPlaneXPosition(sensorID) ; 
            yplane = geo::gGeometry().siPlaneYPosition(sensorID) ; 
            zplane = geo::gGeometry().siPlaneZPosition(sensorID) ;
 	    xrot   = geo::gGeometry().siPlaneXRotation(sensorID) ;
 	    yrot   = geo::gGeometry().siPlaneYRotation(sensorID) ;
	    zrot   = geo::gGeometry().siPlaneZRotation(sensorID) ;

            streamlog_out(MESSAGE4) << std::endl << 
                                 setw(10) << "original " << std::fixed <<
                                 setw( 8) << sensorID << 
                                 setw(13) << setprecision(4) << xplane   << 
                                 setw(13) << setprecision(4) << yplane   << 
                                 setw(13) << setprecision(4) << zplane   << 
                 		 setw(13) << setprecision(4) << xrot     <<
                                 setw(13) << setprecision(4) << yrot     <<
                                 setw(13) << setprecision(4) << zrot     << std::endl;

	    TRotation invR;
	    invR.RotateX( xrot );
	    invR.RotateY( yrot );
	    invR.RotateZ( zrot );
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
 
             streamlog_out(MESSAGE2) << "invR:"<< std::endl;
             streamlog_out(MESSAGE2) << " X: " << setw(20) << invR[0][0] << " " << invR[0][1] << " " << invR[0][2]  << std::endl;
             streamlog_out(MESSAGE2) << " Y: " << setw(20) << invR[1][0] << " " << invR[1][1] << " " << invR[1][2]  << std::endl;
             streamlog_out(MESSAGE2) << " Z: " << setw(20) << invR[2][0] << " " << invR[2][1] << " " << invR[2][2]  << std::endl;
             streamlog_out(MESSAGE2) << ""<< std::endl;
 
           
//	    delta_r0 *= invR;

//#ifdef GEAR_MAJOR_VERSION 
//#if GEAR_VERSION_GE( 17,4)  
// ZY and ZX rotations are calculated wrongly yet, do not implement:
// XYZ shifts and XY rotation seems to be correct
//
            geo::gGeometry().setPlaneXPosition(sensorID,  xplane  +  delta_r0.X() ) ;
            geo::gGeometry().setPlaneYPosition(sensorID,  yplane  +  delta_r0.Y() ) ;
            geo::gGeometry().setPlaneZPosition(sensorID,  zplane  +  delta_r0.Z() ) ;
            geo::gGeometry().setPlaneXRotation(sensorID, (xrot  - dalpha)  ) ;
            geo::gGeometry().setPlaneYRotation(sensorID, (yrot  - dbeta )  ) ;
            geo::gGeometry().setPlaneZRotation(sensorID, (zrot  - dgamma)  ) ;
//#endif
//#endif       
            streamlog_out(MESSAGE4) << setw(10) << "align  " << setw( 8) << " " ;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) << dr0x;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) << dr0y;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) << dr0z; 

            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) << dalpha;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) << dbeta;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) << dgamma ;

            streamlog_out(MESSAGE4) << setw( 3) << " " ;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<  delta_r0.X() ;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<  delta_r0.Y();
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<  delta_r0.Z() << std::endl;

            xplane = geo::gGeometry().siPlaneXPosition(sensorID) ; 
            yplane = geo::gGeometry().siPlaneYPosition(sensorID) ; 
            zplane = geo::gGeometry().siPlaneZPosition(sensorID) ;
 	    xrot   = geo::gGeometry().siPlaneXRotation(sensorID) ;
 	    yrot   = geo::gGeometry().siPlaneYRotation(sensorID) ;
	    zrot   = geo::gGeometry().siPlaneZRotation(sensorID) ;

            streamlog_out(MESSAGE4) <<
                                 setw(10) << "new : " << 
                                 setw( 8) << " "  << 
                                 setw(13) << setprecision(4) << xplane   << 
                                 setw(13) << setprecision(4) << yplane   << 
                                 setw(13) << setprecision(4) << zplane   << 
                 		 setw(13) << setprecision(4) << xrot     <<
                                 setw(13) << setprecision(4) << yrot     <<
                                 setw(13) << setprecision(4) << zrot     << std::endl;


        }
    }

  
    geo::gGeometry().updateGearManager();
  
    gear::GearXML::createXMLFile( gearManager, newGearfileName ) ;

    streamlog_out(MESSAGE4) << "Not implemented" << std::endl;
    
    return;
}

int main( int argc, char ** argv ) {

  streamlog::out.init( std::cout , "pede2lcio output stream") ;

  streamlog::logscope scope(streamlog::out) ;

  scope.setLevel<streamlog::MESSAGE3>() ;



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
    if ( newGearFileName.rfind( ".xml", string::npos ) == string::npos ) {
         newGearFileName.append( ".xml" );
    }
    
    streamlog_out(MESSAGE4) << " oldGear: " << oldGearFileName << " newGear: " << newGearFileName << std::endl;
  }
  
  streamlog_out(MESSAGE4) << "Converting " << pedeFileName << " in " << lcioFileName << std::endl;

  // try to open the input file. This should be a text file
  ifstream pedeFile( pedeFileName.c_str(), ios::in );

  map< int, EUTelAlignmentConstant* > constants_map;
  
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
            EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
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
