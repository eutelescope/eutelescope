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
#include "EUTelGeometricPixel.h"

// system includes <>
#include <iostream>
#include <iomanip>

using namespace eutelescope;

//Default constructor, returns pixel with all fields set to zero
EUTelGeometricPixel::EUTelGeometricPixel():
	EUTelGenericSparsePixel(),
	_posX(0),
	_posY(0),
	_boundX(0),
	_boundY(0)
{
	_noOfElementsDerived = 8;
	_typeDerived = kEUTelGeometricPixel;
}

//Constructor taking all possible eight arguments
EUTelGeometricPixel::EUTelGeometricPixel(short xCoord, short yCoord, float signal, short time, float posX, float posY, float boundX, float boundY):
	EUTelGenericSparsePixel(xCoord, yCoord, signal, time),
	_posX(posX),
	_posY(posY),
	_boundX(boundX),
	_boundY(boundY)
{
	_noOfElementsDerived = 8;
	_typeDerived = kEUTelGeometricPixel;
}

//Constructor taking a EUTelGenericSparsePixel, all geometry related entries are set to zero
EUTelGeometricPixel::EUTelGeometricPixel(EUTelGenericSparsePixel& genericPixel):
	EUTelGenericSparsePixel(genericPixel),
	_posX(0),
	_posY(0),
	_boundX(0),
	_boundY(0)
{
	_noOfElementsDerived = 8;
	_typeDerived = kEUTelGeometricPixel;
}

//Constructor taking a EUTelGenericSparsePixel and all the four geometry related parameters
EUTelGeometricPixel::EUTelGeometricPixel(EUTelGenericSparsePixel& genericPixel, float posX, float posY, float boundX, float boundY):
	EUTelGenericSparsePixel(genericPixel),
	_posX(posX),
	_posY(posY),
	_boundX(boundX),
	_boundY(boundY)
{
	_noOfElementsDerived = 8;
	_typeDerived = kEUTelGeometricPixel;
}


unsigned int EUTelGeometricPixel::getNoOfElements() const
{
	return _noOfElementsDerived;
}

SparsePixelType EUTelGeometricPixel::getSparsePixelType() const 
{
	return _typeDerived;
}

void EUTelGeometricPixel::print(std::ostream& os) const 
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
	 	<< std::setw(width) << "posX: "  << _posX << std::endl
	 	<< std::setw(width) << "posY: "  << _posY << std::endl
	 	<< std::setw(width) << "boundX: "  << _boundX << std::endl
	 	<< std::setw(width) << "boundY: "  << _boundY << std::endl;
	for( int i = 0 ; i < bigWidth ; ++i )
	{
		os << "-";
	}
}
