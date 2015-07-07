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
#include "EUTelBinaryPixel.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace eutelescope;

//Default constructor, returns pixel with all fields set to zero
EUTelBinaryPixel::EUTelBinaryPixel():
	EUTelGenericSparsePixel(),
	_hitTime(0),
	_frameTime(0)
{
	_noOfElementsDerived = 6;
	_typeDerived = kEUTelBinaryPixel;
}

//Constructor taking all possible six arguments
// First four: same as GenericSparsePixel
// Last two: additional arguments to store precise time information 
EUTelBinaryPixel::EUTelBinaryPixel(short xCoord, short yCoord, float signal, short time, double hitTime, double frameTime):
	EUTelGenericSparsePixel(xCoord, yCoord, signal, time),
	_hitTime(hitTime),
	_frameTime(frameTime)
{
	_noOfElementsDerived = 6;
	_typeDerived = kEUTelBinaryPixel;
}

//Constructor taking a EUTelGenericSparsePixel, all geometry related entries are set to zero
EUTelBinaryPixel::EUTelBinaryPixel(EUTelGenericSparsePixel& genericPixel):
	EUTelGenericSparsePixel(genericPixel),
	_hitTime(0),
	_frameTime(0)
{
	_noOfElementsDerived = 6;
	_typeDerived = kEUTelBinaryPixel;
}

//Constructor taking a EUTelGenericSparsePixel and the two time stamps
EUTelBinaryPixel::EUTelBinaryPixel(EUTelGenericSparsePixel& genericPixel, double hitTime, double frameTime):
	EUTelGenericSparsePixel(genericPixel),
	_hitTime(hitTime),
	_frameTime(frameTime)
{
	_noOfElementsDerived = 6;
	_typeDerived = kEUTelBinaryPixel;
}


unsigned int EUTelBinaryPixel::getNoOfElements() const
{
	return _noOfElementsDerived;
}

SparsePixelType EUTelBinaryPixel::getSparsePixelType() const 
{
	return _typeDerived;
}

void EUTelBinaryPixel::print(std::ostream& os) const 
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
	 	<< std::setw(width) << "time: "  << _time << std::endl
	 	<< std::setw(width) << "hit time: "  << _hitTime << std::endl
	 	<< std::setw(width) << "frame time: "  << _frameTime << std::endl;
	  for( int i = 0 ; i < bigWidth ; ++i )
	{
		os << "-";
	}
}
