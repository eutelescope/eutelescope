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
#include "EUTelGenericSparsePixel.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace eutelescope;

//Constructor without parameters, all values are assigned zero
EUTelGenericSparsePixel::EUTelGenericSparsePixel(): 
	_xCoord(0),
	_yCoord(0),
	_signal(0),
	_time(0)
{
	_noOfElements = 4;
	_type = kEUTelGenericSparsePixel;
}



//Constructor with all four parameters
EUTelGenericSparsePixel::EUTelGenericSparsePixel(short xCoord, short yCoord, float signal, short time):
	_xCoord(xCoord),
	_yCoord(yCoord),
	_signal(signal),
	_time(time)
{
	_noOfElements = 4;
	_type = kEUTelGenericSparsePixel;
}

//Constructor with only X,Y,signal. The time is thus set to zero
EUTelGenericSparsePixel::EUTelGenericSparsePixel(short xCoord, short yCoord, float signal):
	_xCoord(xCoord),
	_yCoord(yCoord),
	_signal(signal),
	_time(0)
{
	_noOfElements = 4;
	_type = kEUTelGenericSparsePixel;
}


unsigned int EUTelGenericSparsePixel::getNoOfElements() const
{
	return _noOfElements;
}


SparsePixelType EUTelGenericSparsePixel::getSparsePixelType() const
{
	return _type;
}

void EUTelGenericSparsePixel::print(std::ostream& os) const
{
	int bigWidth = 50;
	for ( int i = 0 ; i < bigWidth ; ++i ) 
	{
		os << "-";
  	}
  	os << std::endl;
  	int width = 20;
  	os 	<< std::setw(width) << std::setiosflags(std::ios::left) << "Type: "     << _type << std::endl
     	<< std::setw(width) << "Elements: " << _noOfElements << std::endl
     	<< std::setw(width) << "x coord: "  << _xCoord << std::endl
     	<< std::setw(width) << "y coord: "  << _yCoord << std::endl
     	<< std::setw(width) << "signal: "  << _signal << std::endl
 	<< std::setw(width) << "time: "  << _time << std::endl;
 	for ( int i = 0 ; i < bigWidth ; ++i ) 
	{
		os << "-";
	}
}
