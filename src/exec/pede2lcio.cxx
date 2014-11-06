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

    map< int, EUTelAlignmentConstant* >::const_iterator itrAlignmentConstant;

   double xplane = 0.;
    double yplane = 0.;
    double zplane = 0.;
    double xrot   = 0.;
    double yrot   = 0.;
    double zrot   = 0.;

		streamlog_out(MESSAGE9) << "Local translation come from Millepede." << std::endl;
		streamlog_out(MESSAGE9) << "Global will be added to gear geometry parameters. " << std::endl;


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

	    TRotation invR;
	    invR.RotateX( xrot );
	    invR.RotateY( yrot );
	    invR.RotateZ( zrot );
	    invR.Invert();
	
	    const double dalpha = ((*itrAlignmentConstant).second->getAlpha())*(180/M_PI);
            const double dbeta  = ((*itrAlignmentConstant).second->getBeta())*(180/M_PI);
            const double dgamma = ((*itrAlignmentConstant).second->getGamma())*(180/M_PI);

	    TRotation invDeltaR;
	    invDeltaR.RotateX((*itrAlignmentConstant).second->getAlpha());
			invDeltaR.RotateY((*itrAlignmentConstant).second->getBeta());
			invDeltaR.RotateZ((*itrAlignmentConstant).second->getGamma());
	    invDeltaR.Invert();

	    const double dr0x = (*itrAlignmentConstant).second->getXOffset();
	    const double dr0y = (*itrAlignmentConstant).second->getYOffset();
	    const double dr0z = (*itrAlignmentConstant).second->getZOffset();
	
			const double posLocalDiff[3] = {dr0x,dr0y,dr0z};
			double delta_r0[3];
			geo::gGeometry().local2MasterVec(sensorID,posLocalDiff, delta_r0);//Here we transform the local alignment position offsets to global position offsets.
			const double posTest[3]={1,0,0};
			double posTestOutput[3];
			geo::gGeometry().local2Master(sensorID,posTest, posTestOutput);
			streamlog_out(MESSAGE9)<<"Here we have the test for sensor " << sensorID <<std::endl;
			streamlog_out(MESSAGE9)<<posTestOutput[0]<<"  "<<posTestOutput[1]<<"  "<<posTestOutput[2]<<endl;	
			const double angleLocalDiff[3]={dalpha,dbeta,dgamma};
			double delta_angle[3];
			//IMPORTANT:Note the transformation of the angles assumes that they transform like a vector. This is not true unless the angles are small.   
			geo::gGeometry().local2MasterVec(sensorID,angleLocalDiff, delta_angle);//Here we transform the local alignment angle offsets to global angle offsets.

 
           
//	    delta_r0 *= invR;

//#ifdef GEAR_MAJOR_VERSION 
//#if GEAR_VERSION_GE( 17,4)  
// ZY and ZX rotations are calculated wrongly yet, do not implement:
// XYZ shifts and XY rotation seems to be correct
//
            geo::gGeometry().setPlaneXPosition(sensorID,  xplane  + delta_r0[0]  ) ;
            geo::gGeometry().setPlaneYPosition(sensorID,  yplane  + delta_r0[1]  ) ;
            geo::gGeometry().setPlaneZPosition(sensorID,  zplane  + delta_r0[2]  ) ;
            geo::gGeometry().setPlaneXRotation(sensorID, (xrot  - angleLocalDiff[0])  ) ;
            geo::gGeometry().setPlaneYRotation(sensorID, (yrot  - angleLocalDiff[1] )  ) ;
            geo::gGeometry().setPlaneZRotation(sensorID, (zrot  - angleLocalDiff[2])  ) ;
//#endif
//#endif       
					 streamlog_out(MESSAGE9) << "Input and output alignment shift (translations) for sensor: "<<sensorID << std::endl;
            streamlog_out(MESSAGE4) << setw(10) << "Align Translations (Local) x,y,z  " << setw( 8) << " " ;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<posLocalDiff[0] ;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<posLocalDiff[1];
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<posLocalDiff[2]<<endl; 

            streamlog_out(MESSAGE4) << setw(10) << "Align Translations (Global) x,y,z " << setw( 8) << " " ;
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<delta_r0[0];
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<delta_r0[1];
            streamlog_out(MESSAGE4) << setw(13) << setprecision(4) <<delta_r0[2]<<endl; 

            xplane = geo::gGeometry().siPlaneXPosition(sensorID) ; 
            yplane = geo::gGeometry().siPlaneYPosition(sensorID) ; 
            zplane = geo::gGeometry().siPlaneZPosition(sensorID) ;
 	    xrot   = geo::gGeometry().siPlaneXRotationRadians(sensorID) ;
 	    yrot   = geo::gGeometry().siPlaneYRotationRadians(sensorID) ;
	    zrot   = geo::gGeometry().siPlaneZRotationRadians(sensorID) ;

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
