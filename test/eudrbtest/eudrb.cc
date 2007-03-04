// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: eudrb.cc,v 1.1 2007-03-04 18:39:35 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */


// personal include ".h"
#include "EUTelEUDRBReader.h"

#ifdef USE_GSL
// gsl include <.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

// system include <>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;
using namespace eutelescope;

const int nEvent     = 100;
const int nDetector  = 1;
const int nChannel   = 4;
const int xNPixel    = 66;
const int yNPixel    = 256;
const int nFrame      = 3;
const string fileName = "eudrb.dat";
const double noise    = 5;



int main(int argc, char ** argv ) {

  ofstream file( fileName.c_str(), ios::out | ios::binary );

#ifdef USE_GSL
  gsl_rng * randomGenerator = gsl_rng_alloc( gsl_rng_taus );
#endif  

  EUDRBFileHeader fileHeader;
  fileHeader.numberOfEvent    = nEvent;
  fileHeader.numberOfDetector = nDetector;
  fileHeader.nXPixel          = xNPixel;
  fileHeader.nYPixel          = yNPixel;
  fileHeader.eventSize        = sizeof(fileHeader) + sizeof(EUDRBTrailer) + 
    ( xNPixel * yNPixel * nChannel * nFrame ) / 2 * sizeof(int);
  fileHeader.dataSize         = ( xNPixel * yNPixel * nChannel * nFrame ) / 2 * sizeof(int);
  fileHeader.chACBitMask      = 0x0FFF0000;
  fileHeader.chACRightShift   = 16;
  fileHeader.chBDBitMask      = 0x00000FFF;
  fileHeader.chBDRightShift   = 0;

  file.write( reinterpret_cast< char * > (&fileHeader), sizeof( fileHeader ) );

  int * buffer = new int[ ( xNPixel * yNPixel * nChannel * nFrame ) / 2 ];
  
  for (int iEvent = 0; iEvent < nEvent; iEvent++ ) {
    EUDRBEventHeader eventHeader;
    eventHeader.eventNumber   = iEvent;
    eventHeader.triggerNumber = iEvent;
    file.write(reinterpret_cast< char * > (&eventHeader), sizeof( eventHeader ) );

    for ( int iRecord = 0; iRecord < ( xNPixel * yNPixel * nChannel * nFrame ) / 2; iRecord++ ) {

#ifdef USE_GSL      
      short pixel0     = ( 2048 +  static_cast<short> (gsl_ran_gaussian( randomGenerator , noise)) ) & 0xFFF;
      short pixel1     = ( 2048 +  static_cast<short> (gsl_ran_gaussian( randomGenerator , noise)) ) & 0xFFF; 
#else
      short pixel0     = ( 2048 - static_cast<short> (noise) + (rand() / (RAND_MAX / static_cast<short> (2 * noise) )) );
      short pixel1     = ( 2048 - static_cast<short> (noise) + (rand() / (RAND_MAX / static_cast<short> (2 * noise) )) );
#endif
      buffer[iRecord]  = ( pixel0 << fileHeader.chACRightShift ) + pixel1;
      
    }

    file.write(reinterpret_cast< char * > (&buffer[0]), fileHeader.dataSize ) ;
    
    EUDRBTrailer eventTrailer;
    eventTrailer.trailer = 0x89abcdef;
    file.write(reinterpret_cast<char*>(&eventTrailer), sizeof (eventTrailer) );


  }
  
  delete [] buffer;
#ifdef USE_GSL
  gsl_rng_free( randomGenerator );
#endif
  file.close();

  return 0;

}
