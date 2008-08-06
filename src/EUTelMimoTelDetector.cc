// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Author Loretta Negrini, Univ. Insubria <mailto:loryneg@gmail.com>
// Author Silvia Bonfanti, Univ. Insubria <mailto:silviafisica@gmail.com>
// Version $Id: EUTelMimoTelDetector.cc,v 1.1 2008-08-06 20:37:00 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelMimoTelDetector.h" 

// system includes <>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;
using namespace eutelescope;

EUTelMimoTelDetector::EUTelMimoTelDetector() : EUTelBaseDetector() {

  _xMin = 0;
  _xMax = 263;
  
  _yMin = 0;
  _yMax = 255;
  
  _markerPos.push_back( 0 );
  _markerPos.push_back( 1 );
  _markerPos.push_back( 66 );
  _markerPos.push_back( 67 );
  _markerPos.push_back( 132 );
  _markerPos.push_back( 133 );
  _markerPos.push_back( 198 );
  _markerPos.push_back( 199 );

  _signalPolarity = -1;

  _name = "MimoTel";

  _xPitch = 0.03;
  _yPitch = 0.03;

}

void EUTelMimoTelDetector::setMode( string mode ) {
  
  _mode = mode;

}

void EUTelMimoTelDetector::print( ostream& os ) const {

  size_t w = 35;

  string pol = "negative";
  if ( _signalPolarity > 0 ) pol = "positive"; 

  os << setw( w ) << "Detector name" << _name << endl
     << setw( w ) << "Mode" << _mode << endl
     << setw( w ) << "Pixel along x from" << _xMin << " to " << _xMax << endl
     << setw( w ) << "Pixel along y from" << _yMin << " to " << _yMax << endl
     << setw( w ) << "Pixel pitches " << _xPitch << ", " << _yPitch << endl
     << setw( w ) << "Signal polarity" << pol << endl;
   
  if ( hasMarker() ) {

    os << "Detector has the following colomns used as markers: "<< endl;

    vector< size_t >::const_iterator iter = getMarkerPosition().begin();
    while ( iter != getMarkerPosition().end() ) {
      os << "x = " << (*iter) << endl;
      
      ++iter;
    }

  }
}
