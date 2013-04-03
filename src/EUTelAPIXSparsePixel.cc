// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id$
 

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
#include "EUTelAPIXSparsePixel.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace eutelescope;


EUTelAPIXSparsePixel::EUTelAPIXSparsePixel()
: EUTelBaseSparsePixel(),
  _xCoord(0),
  _yCoord(0),
  _signal(0),
  _chip(0),
  _time(0)
{
  _noOfElements = 5,
  _type = kEUTelAPIXSparsePixel;
}

EUTelAPIXSparsePixel::EUTelAPIXSparsePixel(short xCoord, short yCoord, short signal, short chip, short time)
: EUTelBaseSparsePixel(),
  _xCoord(xCoord),
  _yCoord(yCoord),
  _signal(signal),
  _chip(chip),
  _time(time)
{
  _noOfElements = 5,
  _type = kEUTelAPIXSparsePixel;
}

EUTelAPIXSparsePixel::EUTelAPIXSparsePixel(const EUTelAPIXSparsePixel &orig)
: EUTelBaseSparsePixel(),
  _xCoord(orig.getXCoord()),
  _yCoord(orig.getYCoord()),
  _signal(static_cast< short >(orig.getSignal())),
  _chip(orig.getChip()),
  _time(static_cast< short >(orig.getTime()))
{
  _noOfElements = orig.getNoOfElements();
  _type = orig.getSparsePixelType();
}


unsigned int EUTelAPIXSparsePixel::getNoOfElements() const {
  return _noOfElements;
}

SparsePixelType EUTelAPIXSparsePixel::getSparsePixelType() const {
  return _type;
}

void EUTelAPIXSparsePixel::print(std::ostream& os) const {
  int bigWidth = 50;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
  os << endl;
  int width = 20;
  os << setw(width) << setiosflags(ios::left) << "Type: "     << _type << endl
     << setw(width) << "Elements: " << _noOfElements << endl
     << setw(width) << "x coord: "  << _xCoord << endl
     << setw(width) << "y coord: "  << _yCoord << endl
     << setw(width) << "signal: "  << _signal << endl
     << setw(width) << "chip: "  << _chip << endl
     << setw(width) << "time: "  << _time << endl;
  for ( int i = 0 ; i < bigWidth ; ++i ) {
    os << "-";
  }
}
