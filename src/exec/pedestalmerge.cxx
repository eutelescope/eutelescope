// eutelescope includes ""
#include "anyoption.h"
#include "EUTELESCOPE.h"

// lcio includes <>
#include <IO/LCWriter.h>
#include <IO/LCReader.h>
#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>

//system includes <>
#include <glob.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include <map>

using namespace std;
using namespace IMPL;



int main( int argc, char ** argv ) {

  auto_ptr< AnyOption > option( new AnyOption );

  string usageString =
    "\n"
    "This program can be used to merge together more than one pedestal run \n"
    "in particular when the user wants to collect in one file one the telescope \n"
    "reference plane pedestal info with the DUT ones\n"
    "\n"
    "pedestalmerge [option] -o outputfile.slcio file1.slcio file2.slcio [fileN.slcio]\n"
    "\n"
    "-h --help         Print this help\n";

  option->addUsage( usageString.c_str() );
  option->setFlag( "help", 'h');
  option->setOption( "output", 'o' );

  option->processCommandArgs( argc,  argv );

  if ( option->getFlag('h') || option->getFlag( "help" ) ) {
    option->printUsage();
    return 0;
  }


  if ( option->getValue( "output" ) == NULL ) {
    cerr << "Please provide an output file name using -o option" << endl;
    return 2;
  }

  string outputFileName = option->getValue( "output" );
  // check if the output lcio file has the extension
  if ( outputFileName.rfind( ".slcio", string::npos ) == string::npos ) {
    outputFileName.append( ".slcio" );
  }

  // the input files may be using wildcards
  glob_t globbuf;
  for ( size_t iArg = 0 ; iArg < static_cast<size_t>(option->getArgc()); ++iArg ) {
    if ( iArg == 0 ) glob( option->getArgv( iArg ), 0, NULL, &globbuf);
    else  glob( option->getArgv( iArg ), GLOB_APPEND, NULL, &globbuf);
  }

  // moving to a vector of string because it's easier
  vector< string > inputFileNames( &globbuf.gl_pathv[0], &globbuf.gl_pathv[ globbuf.gl_pathc ] );

  // check how many good input file we got. If it is less than 1, no need to go further.
  if ( inputFileNames.size() <= 1 ) {
    cerr << "Please provide at least two valid input files" << endl;
    return 1;
  }

  // print some information
  cout << "Target file: " << outputFileName << endl;
  for ( size_t iFile = 0; iFile < inputFileNames.size() ; ++iFile ) {
    cout << "Input file: " << inputFileNames.at( iFile ) << endl;
  }

  // open the LCIO output file
  lcio::LCWriter * lcWriter = lcio::LCFactory::getInstance()->createLCWriter();

  try {
      lcWriter->open( outputFileName.c_str() , lcio::LCIO::WRITE_NEW );
  } catch ( lcio::IOException& e ) {
    cerr << e.what() << endl;
    return 3;
  }


  // write an almost empty run header
  lcio::LCRunHeaderImpl * lcHeader  = new lcio::LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );
  lcWriter->writeRunHeader(lcHeader);
  delete lcHeader;

  // prepare an event
  lcio::LCEventImpl * event = new lcio::LCEventImpl;
  event->setRunNumber( 0 );
  event->setEventNumber( 0 );

  lcio::LCTime * now = new lcio::LCTime;
  event->setTimeStamp( now->timeStamp() );
  delete now;

  // not all the collections are merged, but only the following types:
  string acceptedType = "";
  acceptedType.append( lcio::LCIO::TRACKERRAWDATA );
  acceptedType.append( lcio::LCIO::TRACKERDATA );

  // we need to loop over input files twice. The first time we need to read all the names of the 
  // input collections.
  map< string , string > collectionNameTypeMap;

  // the reading factory
  lcio::LCReader * lcReader = lcio::LCFactory::getInstance()->createLCReader();

  for ( size_t iFile = 0 ; iFile < inputFileNames.size(); ++iFile ) {

    try {
      lcReader->open( inputFileNames.at( iFile ).c_str() );

      lcio::LCEventImpl * inputEvent = dynamic_cast< lcio::LCEventImpl* > ( lcReader->readNextEvent() ) ;
      vector< string > inputCollectionNames = *inputEvent->getCollectionNames();
      for ( size_t iCol = 0; iCol < inputCollectionNames.size(); ++iCol ) {
        lcio::LCCollectionVec * inputCollection = dynamic_cast< lcio::LCCollectionVec* > ( inputEvent->getCollection( inputCollectionNames.at( iCol ).c_str() ) );
        string type = inputCollection->getTypeName();
        if ( acceptedType.find( type ) != string::npos ) {
          map< string, string >::iterator iter = collectionNameTypeMap.find( inputCollectionNames.at( iCol ) );
          if ( iter != collectionNameTypeMap.end() ) {
            if ( iter->second != type ) {
              cerr << "Error! Collection " << inputCollectionNames.at( iCol ) << " is found to be both " << type << " and " << iter->second << endl;
              return 4;
            }
          }
          collectionNameTypeMap.insert( make_pair ( inputCollectionNames.at( iCol ), type ) );
        }
      }

      lcReader->close();

    } catch ( lcio::IOException& e ) {
      cerr << e.what() << endl;
    }

  }

  // print the names of the collections that we will be merging and in the same time prepare a map to store the physical output collections. 
  map< string , lcio::LCCollectionVec *> collectionMap;
  map< string , string >::iterator collectionIterator = collectionNameTypeMap.begin();

  cout << "Input collections found: " << endl;
  while ( collectionIterator != collectionNameTypeMap.end() ) {
    string name = collectionIterator->first;
    string type = collectionIterator->second ;
    cout << "-->\t" << name <<  " (" << type << ")" << endl;
    collectionMap[ name ] = new lcio::LCCollectionVec( type );

    if ( type == lcio::LCIO::TRACKERRAWDATA )  {
      lcio::CellIDEncoder<TrackerRawDataImpl>    encoderEncoder( eutelescope::EUTELESCOPE::MATRIXDEFAULTENCODING, collectionMap[ name ]);
    } else if ( type == lcio::LCIO::TRACKERDATA )  {
      lcio::CellIDEncoder<TrackerDataImpl>       encoderEncoder( eutelescope::EUTELESCOPE::MATRIXDEFAULTENCODING, collectionMap[ name ]);
    }


    ++collectionIterator;
  }

  // reloop once again, but this time do the real job! 
  for ( size_t iFile = 0 ; iFile < inputFileNames.size(); ++iFile ) {

    try {

      lcReader->open( inputFileNames.at( iFile ).c_str() );

      lcio::LCEventImpl * inputEvent = dynamic_cast< lcio::LCEventImpl* > ( lcReader->readNextEvent() ) ;

      // loop over all the expected collection names
      map< string, string >::iterator mapIter = collectionNameTypeMap.begin() ;

      while ( mapIter != collectionNameTypeMap.end() ) {

        string type = mapIter->second;

        try {
          lcio::LCCollectionVec * inputCollection = dynamic_cast< lcio::LCCollectionVec * > ( inputEvent->getCollection( mapIter->first ) );

          for ( size_t iElement = 0 ; iElement < inputCollection->size() ; ++iElement ) {

            if ( type == lcio::LCIO::TRACKERRAWDATA ) {
              lcio::TrackerRawDataImpl * input = dynamic_cast< lcio::TrackerRawDataImpl * > ( inputCollection->getElementAt( iElement ) ) ;
              lcio::TrackerRawDataImpl * output = new lcio::TrackerRawDataImpl;
              output->setADCValues( input->getADCValues()  ) ;
              output->setCellID0  ( input->getCellID0() ) ;
              output->setCellID1  ( input->getCellID1() ) ;
              output->setTime     ( input->getTime()    ) ;
              collectionMap[ mapIter->first]->addElement( output );
            } else if ( type == lcio::LCIO::TRACKERDATA ) {
              lcio::TrackerDataImpl * input  = dynamic_cast< lcio::TrackerDataImpl * > ( inputCollection->getElementAt( iElement ) ) ;
              lcio::TrackerDataImpl * output = new lcio::TrackerDataImpl;
              output->setChargeValues( input->getChargeValues()  ) ;
              output->setCellID0  ( input->getCellID0() ) ;
              output->setCellID1  ( input->getCellID1() ) ;
              output->setTime     ( input->getTime()    ) ;
              collectionMap[ mapIter->first]->addElement( output );
            }

          }


        } catch ( lcio::DataNotAvailableException& e ) {
        }
        ++mapIter;
      }

      lcReader->close();

    } catch ( lcio::IOException& e ) {
      cerr << e.what() << endl;
    }
  }

  map<string, lcio::LCCollectionVec * >::iterator collectionIter = collectionMap.begin();
  while ( collectionIter != collectionMap.end() ) {
    event->addCollection( collectionIter->second, collectionIter->first );
    ++collectionIter;
  }

  lcWriter->writeEvent( event );
  delete event;

  lcWriter->close();

  return 0;

}
