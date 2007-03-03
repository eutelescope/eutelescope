// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelEtaFunctionImpl.cc,v 1.2 2007-03-03 08:56:26 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelEtaFunctionImpl.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

// system includes <>
#include <vector>
#include <algorithm>

using namespace lcio;
using namespace eutelescope;
using namespace std;

EUTelEtaFunctionImpl::EUTelEtaFunctionImpl(int nBin) : IMPL::LCGenericObjectImpl(0,0,2 * nBin) {
  _typeName        = "Eta function";
  _dataDescription = "The first nBin double values are the bin centers, the second nBin doubles values"
    " the corresponding eta function value";
  _isFixedSize     = true;
}

EUTelEtaFunctionImpl::EUTelEtaFunctionImpl(int nBin, vector<double > centerVec, vector<double > valueVec) : 
  IMPL::LCGenericObjectImpl(0,0,2 * nBin) {
    _typeName        = "Eta function";
    _dataDescription = "The first nBin double values are the bin centers, the second nBin doubles values"
      " the corresponding eta function value";
    _isFixedSize     = true;
   
 
    unsigned int indexCenter;
    unsigned int indexValue ;
   
    for (indexCenter = 0, indexValue = indexCenter + getNDouble() / 2; 
	 indexCenter < centerVec.size(); indexCenter++, indexValue++ ) {
      setDoubleVal(indexCenter, centerVec[indexCenter]);
      setDoubleVal(indexValue,  valueVec[indexCenter]);
    }
      
}

void EUTelEtaFunctionImpl::setBinCenterVector(vector<double > center) {

  unsigned int indexCenter;
  for ( indexCenter = 0; indexCenter < center.size(); indexCenter++ ) {
    setDoubleVal( indexCenter, center[indexCenter] );
  }
}

void EUTelEtaFunctionImpl::setEtaValueVector(vector<double > value) {

  unsigned int index;
  int shift = getNDouble() / 2;
  for ( index = 0; index < value.size(); index++ ) {
    setDoubleVal( index + shift, value[index] );
  }
}
  

const vector<double > EUTelEtaFunctionImpl::getBinCenterVector() const {

  vector<double > center;
  int index;
  for (index = 0; index < getNDouble() / 2; index++) {
    center.push_back(getDoubleVal(index));
  }

  return center;

}

const vector<double > EUTelEtaFunctionImpl::getEtaValueVector() const {

  vector<double > value;
  int index;
  for (index = getNDouble() / 2; index < getNDouble(); index++) {
    value.push_back(getDoubleVal(index));
  }

  return value;

}

double EUTelEtaFunctionImpl::getEtaFromCoG(double x) const {

  typedef vector<double >::const_iterator DoubleIter;
  
  DoubleIter xLeft    = lower_bound(getCoGBeginConstIterator(), getCoGEndConstIterator(), x);
  DoubleIter xRight   = xLeft + 1;
  DoubleIter etaLeft  = getEtaBeginConstIterator() +  ( xLeft - getCoGBeginConstIterator());
  DoubleIter etaRight = etaLeft + 1;

  return *etaLeft + ( *etaLeft - *etaRight ) / ( *xLeft - *xRight ) * ( x - *xLeft) ;

}

vector<double >::const_iterator EUTelEtaFunctionImpl::getCoGBeginConstIterator() const {
  return _doubleVec.begin();
}

vector<double >::const_iterator EUTelEtaFunctionImpl::getCoGEndConstIterator() const {
  return _doubleVec.begin() + (getNDouble() / 2 );
}

vector<double >::const_iterator EUTelEtaFunctionImpl::getEtaBeginConstIterator() const {
  return _doubleVec.begin() + (getNDouble() / 2);
} 

vector<double >::const_iterator EUTelEtaFunctionImpl::getEtaEndConstIterator() const {
  return _doubleVec.end();
} 
