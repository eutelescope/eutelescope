// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#include "EUTelPseudo1DHistogram.h"


// system includes
#include <cstdlib>

using namespace eutelescope;


//=============================================================================

EUTelPseudo1DHistogram::EUTelPseudo1DHistogram(int NOfBins, double min, double max) {

  _FullNumberOfBins = NOfBins + 2;
  _NumberOfBins = NOfBins;
  _MinValue = min;
  _MaxValue = max;
  _NOfEntries = new int[_FullNumberOfBins];
  _Content = new double[_FullNumberOfBins];
  _UpperIntervalLimit = new double[_FullNumberOfBins-1];
  for (int i = 0; i < _FullNumberOfBins; ++i) {
    _NOfEntries[i] = 0;
    _Content[i] = 0.0;
  }

  _BinWidth = (fabs(_MaxValue - _MinValue)) / _NumberOfBins;
  for (int i = 0; i < _FullNumberOfBins-1; ++i) {
    _UpperIntervalLimit[i] = _MinValue + i*_BinWidth;
  }

}

//=============================================================================

EUTelPseudo1DHistogram::~EUTelPseudo1DHistogram() {

  _FullNumberOfBins = 0;
  delete[] _NOfEntries;
  _NOfEntries = 0;
  delete[] _Content;
  _Content = 0;
  delete[] _UpperIntervalLimit;
  _UpperIntervalLimit = 0;

}

//=============================================================================

void EUTelPseudo1DHistogram::clearContent() {

  for (int i = 0; i < _FullNumberOfBins; ++i) {
    _NOfEntries[i] = 0;
    _Content[i] = 0.0;
  }

}

//=============================================================================

void EUTelPseudo1DHistogram::fill(double x, double w) {

  if (x < _MinValue) {
    ++_NOfEntries[0];
    _Content[0] += w;  
  }
  else if (x > _MaxValue) {
    ++_NOfEntries[_FullNumberOfBins-1];
    _Content[_FullNumberOfBins-1] += w;  
  }
  else if (x == _MaxValue) {
    ++_NOfEntries[_FullNumberOfBins-2];
    _Content[_FullNumberOfBins-2] += w;  
  }
  else {
    int index = static_cast< int >(floor(((static_cast< double >(_NumberOfBins))/(_MaxValue - _MinValue))*(x - _MinValue))) + 1;
    ++_NOfEntries[index];
    _Content[index] += w;
  }

}

//=============================================================================

int EUTelPseudo1DHistogram::findBin(double x) {

  if (x < _MinValue) {
    return 0;
  }
  else if (x > _MaxValue) {
    return _FullNumberOfBins-1;
  }
  else if (x == _MaxValue) {
    return _FullNumberOfBins-2;
  }
  else {
    return static_cast< int >(floor(((static_cast< double >(_NumberOfBins))/(_MaxValue - _MinValue))*(x - _MinValue))) + 1;
  }

}

//=============================================================================

double EUTelPseudo1DHistogram::getBinContent(int bin) {

  if ( isInRange(bin) ) {
    return _Content[bin];
  }
  else {
    std::cout << "Requested bin not in range of the 'EUTelPseudo1DHistogram'." << std::endl;
    return -1.0;
  }
  
}

//=============================================================================

int EUTelPseudo1DHistogram::getNumberOfEntries(int bin) {

  if ( isInRange(bin) ) { 
    return _NOfEntries[bin];
  }
  else {
    std::cout << "Requested bin not in range of the 'EUTelPseudo1DHistogram'." << std::endl;
    return -1;
  }

}

//=============================================================================

bool EUTelPseudo1DHistogram::isInRange(int bin) {

  if ( (bin >= 0) && (bin < _FullNumberOfBins) ) { 
    return true;
  }
  else return false;

}

//=============================================================================

double EUTelPseudo1DHistogram::integral(int startbin, int endbin) {

  if ( isInRange(startbin) && isInRange(endbin) ) {    
    double result = 0.0;
    for (int i = startbin; i < abs(endbin - startbin)+startbin+1; ++i) {
      result += getBinContent(i);
    }
    return result;
  }
  else {
    std::cout << "At least one requested bin is not in range of the 'EUTelPseudo1DHistogram'."
	      << std::endl;
    return -1.0;
  }
  
}

//=============================================================================

void EUTelPseudo1DHistogram::printContent() {

  std::cout << "bin" << "     " << "content" << "     " << "entries" << "     " 
	    << "interval" << std::endl;

  std::cout.precision(3);
  for (int i = 0; i < _FullNumberOfBins; ++i) {
    std::cout << i << "          " << _Content[i] << "          " << _NOfEntries[i] 
	      << "          ";
    if (i==0) {
      std::cout << "(-inf," << _UpperIntervalLimit[i] << "]" << std::endl;
    }
    else if (i==_FullNumberOfBins-1) {
      std::cout << "[" << _UpperIntervalLimit[i-1] << ",inf)" << std::endl;
    }
    else {
      std::cout << "[" << _UpperIntervalLimit[i-1] << "," <<  _UpperIntervalLimit[i]
		<< "]" << std::endl;
    }
  }
}

//=============================================================================

int EUTelPseudo1DHistogram::getNumberOfBins() {
  return _NumberOfBins; 
}

//=============================================================================

double EUTelPseudo1DHistogram::getBinCenter(int index){
  if (isInRange(index)) {
    return _UpperIntervalLimit[index] - (0.5 * _BinWidth) ;
  } else {
    std::cout << "At least one requested bin is not in range of the 'EUTelPseudo1DHistogram'."
	      << std::endl;
    return -1.0;
  }
}
