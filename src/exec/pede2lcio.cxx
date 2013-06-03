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

// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <functional>
#include <algorithm>
#include <memory>

using namespace std;

struct CollectionWriter {
    lcio::LCCollectionVec* _constantsCollection;
        void operator()( std::pair<const int, eutelescope::EUTelAlignmentConstant*>& pair ) {
          _constantsCollection->push_back( pair.second );
                cout << (*pair.second) << endl;
        }
} colWriter;

int main( int argc, char ** argv ) {

  auto_ptr< AnyOption > option( new AnyOption );

  string usageString = 
    "\n"
    "Eutelescope alignment constant converter\n"
    "\n"
    "This program is used to convert the output of pede into\n"
    "a LCIO file containing alignment constants.\n"
    "\n"
    "Usage: pede2lcio [options] conversion_file lciofile.slcio\n\n"
    "-h --help       To print this help\n"
    "\n";

  option->addUsage( usageString.c_str() );
  option->setFlag( "help", 'h');

  // process the command line arguments
  option->processCommandArgs( argc, argv );

  if ( option->getFlag('h') || option->getFlag( "help" ) ) {
    option->printUsage();
    return 0;
  }

  if ( ( option->getArgc() == 0 ) || 
       ( option->getArgc() != 2 ) ) {
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

  cout << "Converting " << pedeFileName << " in " << lcioFileName << endl;

  // try to open the input file. This should be a text file
  ifstream pedeFile( pedeFileName.c_str(), ios::in );

  if ( pedeFile.fail() ) {

    cerr << "Error opening the " << pedeFileName << endl;
    return -1;

  } else {

    // open the LCIO output file
    lcio::LCWriter * lcWriter = lcio::LCFactory::getInstance()->createLCWriter();

    try {
      lcWriter->open( lcioFileName.c_str() , lcio::LCIO::WRITE_NEW );
    } catch ( lcio::IOException& e ) {
      cerr << e.what() << endl;
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

    map< int, eutelescope::EUTelAlignmentConstant* > constants_map;

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
