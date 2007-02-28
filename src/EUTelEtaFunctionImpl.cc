#include "EUTelEtaFunctionImpl.h"

#include <lcio.h>
#include <IMPL/LCGenericObjectImpl.h>

#include <vector>

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
  

vector<double > EUTelEtaFunctionImpl::getBinCenterVector() const {

  vector<double > center;
  int index;
  for (index = 0; index < getNDouble() / 2; index++) {
    center.push_back(getDoubleVal(index));
  }

  return center;

}

vector<double > EUTelEtaFunctionImpl::getEtaValueVector() const {

  vector<double > value;
  int index;
  for (index = getNDouble() / 2; index < getNDouble(); index++) {
    value.push_back(getDoubleVal(index));
  }

  return value;

}
