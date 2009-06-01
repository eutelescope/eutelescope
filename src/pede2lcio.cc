// eutelescope includes
#include "EUTelAlignmentConstant.h"

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

using namespace std;

void usage(const char * programName) {
  cout << "Eutelescope alignment constant converter" << endl << endl
       << "This program is used to convert the output of pede " << endl
       << "into a LCIO file containing alignment constants " << endl << endl
       << "Usage: " << endl
       << programName << " pedefile.txt lciofile.slcio " << endl;
}


int main( int argc, char ** argv ) {

  if ( argc != 3 ) {
    usage(argv[0]);
    return 0;
  }

  string pedeFileName, lcioFileName;
  for ( int iArg = 1 ; iArg < 2; ++iArg ) {
    if ( ( string(argv[iArg]) == "-h" ) ||
         ( string(argv[iArg]) == "--help") ) {
      usage(argv[0]);
      return 0;
    }
  }
  pedeFileName = string(argv[1]);
  lcioFileName = string(argv[2]);

  // check if the lcio file has the extension
  if ( lcioFileName.rfind( ".slcio", string::npos ) == string::npos ) {
    lcioFileName.append( ".slcio" );
  }

  cout << "Converting " << pedeFileName << " in " << lcioFileName << endl;

  // try to open the input file. This should be a text file
  ifstream pedeFile( pedeFileName.c_str(), ios::in );

  if ( pedeFile.fail() ) {

    cerr << "Error opening the " << pedeFileName << endl;
    return 1;

  } else {

    // open the LCIO output file
    lcio::LCWriter * lcWriter = lcio::LCFactory::getInstance()->createLCWriter();

    try {
      lcWriter->open( lcioFileName.c_str() , lcio::LCIO::WRITE_NEW );
    } catch ( lcio::IOException& e ) {
      cerr << e.what() << endl;
      return 2;
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


    vector<double > tokens;
    stringstream tokenizer;
    string line;
    double buffer;

    // get the first line and throw it away since it is a
    // comment!
    getline( pedeFile, line );

    int sensorID = 0;

    while ( true  ) {

      eutelescope::EUTelAlignmentConstant * constant = new eutelescope::EUTelAlignmentConstant;

      constant->setSensorID( sensorID );
      ++sensorID;

      bool goodLine = true;

      for ( unsigned int iParam = 0 ; iParam < 3 ; ++iParam ) {

        getline( pedeFile, line );

        if ( line.empty() ) {
          goodLine = false;
        }

        tokens.clear();
        tokenizer.clear();
        tokenizer.str( line );

        while ( tokenizer >> buffer ) {
          tokens.push_back( buffer ) ;
        }

        if ( ( tokens.size() == 3 ) || ( tokens.size() == 6 ) ) {
          goodLine = true;
        } else goodLine = false;

        bool isFixed = ( tokens.size() == 3 );
//         if ( isFixed ) {
//           cout  << "Parameter " << tokens[0] << " is at " << ( tokens[1] / 1000 )
//                << " (fixed)"  << endl;
//         } else {
//           cout << "Parameter " << tokens[0] << " is at " << (tokens[1] / 1000 )
//                << " +/- " << ( tokens[4] / 1000 )  << endl;
//         }

        if ( iParam == 0 ) {
          constant->setXOffset( tokens[1] / 1000 );
          if ( ! isFixed ) {
            double err  = tokens[4] / 1000;
            constant->setXOffsetError( err ) ;
          }
        }
        if ( iParam == 1 ) {
          constant->setYOffset( tokens[1] / 1000 ) ;
          if ( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000 ) ;
        }
        if ( iParam == 2 ) {
          constant->setGamma( tokens[1]  ) ;
          if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
        }

      }

      // right place to add the constant to the collection
      if ( goodLine ) {
        constantsCollection->push_back( constant );
        cout << (*constant) << endl;
      }
      else delete constant;
    
      if ( !pedeFile ) break;
    }


    event->addCollection( constantsCollection, "alignment" );
    lcWriter->writeEvent( event );
    delete event;

    lcWriter->close();
  }


  pedeFile.close();


  return 0;
}
