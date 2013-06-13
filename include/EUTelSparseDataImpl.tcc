// Version: $Id$
// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELSPARSEDATAIMPL_TCC
#define EUTELSPARSEDATAIMPL_TCC

namespace eutelescope {

  template<class PixelType> 
  EUTelSparseDataImpl<PixelType>::EUTelSparseDataImpl(IMPL::TrackerDataImpl * data) :
    _trackerData(),
    _nElement(0),
    _type(kUnknownPixelType)
  {

    std::auto_ptr<PixelType> pixel ( new PixelType );
    _nElement     = pixel->getNoOfElements();
    _type         = pixel->getSparsePixelType();
    _trackerData  = data;
  } 
  
  template<class PixelType>
  unsigned int EUTelSparseDataImpl<PixelType>::size() const {
    return _trackerData->getChargeValues().size() / _nElement;
  }

  template<class PixelType> 
  EUTelSparseDataImpl<PixelType>::EUTelSparseDataImpl(const EUTelSparseDataImpl &z) : _trackerData(NULL), _nElement(0), _type(0) {
    _trackerData->setCellID0(z->trackerData()->getCellID0());
    _trackerData->setCellID1(z->trackerData()->getCellID1());
    _trackerData->setTime(z->trackerData()->getTime());
    _trackerData->setChargeValues(z->trackerData()->chargeValues());
    _nElement = z->getNElement();
    _type = z->getType();
  }
  
  template<class PixelType> 
  EUTelSparseDataImpl<PixelType>& EUTelSparseDataImpl<PixelType>::operator = (const EUTelSparseDataImpl &z){
    if (this == &z) return *this;  //This handles self assignment
    _trackerData->setCellID0(z->trackerData()->getCellID0());
    _trackerData->setCellID1(z->trackerData()->getCellID1());
    _trackerData->setTime(z->trackerData()->getTime());
    _trackerData->setChargeValues(z->trackerData()->chargeValues());
    _nElement = z->getNElement();
    _type = z->getType();
    return *this;
  }

  template<class PixelType>
  std::list<std::list< unsigned int> > EUTelSparseDataImpl<PixelType>::findNeighborPixels(double minDistance) const {

    PixelType * pixel      = new PixelType;
    PixelType * otherPixel = new PixelType;

    std::list< std::list< unsigned int> > listOfList;

    std::vector<int> status(size(), 0);
    for ( unsigned int iPixel = 0 ; iPixel < size() ; iPixel++ ) {
      
      // grouped pixel is a vector containing the pixels that are
      // fulfilling the minDistance condition
      std::list<unsigned int> groupedPixel;
      
      getSparsePixelAt( iPixel, pixel);

      if ( status[iPixel] == 0 ) {
	groupedPixel.push_back(iPixel);
	status[iPixel] = iPixel + 1;
	
	// first of all, find all the pixels close to the current one
	// and add them to the groupedPixel vector.
	for ( unsigned int iPixel2 = iPixel + 1 ; iPixel2 < size() ; iPixel2++ ) {
	  if ( status[iPixel2] == 0 ) {
	    getSparsePixelAt(iPixel2, otherPixel);
	    if ( distance(pixel, otherPixel) <= minDistance ) {
	      groupedPixel.push_back(iPixel2);
	      status[iPixel2] = iPixel + 1;
	    }
	  }
	}
	
	if ( groupedPixel.size() > 1 ) {
	  // the current group of pixel has more than one pixel in
	  // it. A part from the first, we have to check if the other
	  // are close enough to other pixels
	  std::list<unsigned int>::iterator firstIter = groupedPixel.begin();
	  ++firstIter;
	  while ( firstIter != groupedPixel.end() ) {
	    for ( unsigned int iPixel2 = (*firstIter) + 1 ; iPixel2 < size(); iPixel2++ ) {
	      if ( status[iPixel2] == 0 ) {
		if ( distance( getSparsePixelAt(*firstIter, pixel) , getSparsePixelAt(iPixel2, otherPixel) ) <= minDistance ) {
		  groupedPixel.push_back(iPixel2);
		  status[iPixel2] = iPixel + 1;
		}
	      }
	    }
	    ++firstIter;
	  }
	}
      }
      
      if ( !groupedPixel.empty() )   listOfList.push_back( groupedPixel );
      
    }
    

    delete pixel;
    delete otherPixel;
    
    return listOfList;

  }




}


#endif
