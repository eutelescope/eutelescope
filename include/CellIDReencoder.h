#ifndef CELLIDREENCODER_H
#define CELLIDREENCODER_H

#include "UTIL/CellIDEncoder.h"

namespace lcio{
namespace UTIL{

template <class T>
class CellIDReencoder : public CellIDEncoder<T>{

public:
	CellIDReencoder( std::string encoding, EVENT::LCCollection* col ): CellIDEncoder<T>( encoding, col )
	{
	}
	
	void readValues(T* hit)
	{
		 long64 val = ( long64(hit->getCellID0()& 0xffffffff) ) | ( long64(hit->getCellID1()) << 32 );
		 this->setValue( val );
	}
};

}} //namespaces

#endif
