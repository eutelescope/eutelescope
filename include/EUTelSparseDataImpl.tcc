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
  EUTelSparseDataImpl<PixelType>::EUTelSparseDataImpl(IMPL::TrackerDataImpl * data) {
    // to work properly the template class has to be or at least
    // inherit from EUTelBaseSparsePixel
    EUTelBaseSparsePixel      * goodPixel;
    std::auto_ptr<PixelType>    pixel(new PixelType);
    goodPixel = dynamic_cast<EUTelBaseSparsePixel *>( pixel.get() );
    if ( goodPixel != 0x0 ) {
      // the template class should be a good sparse pixel so we can
      // continue...
      _nElement     = goodPixel->getNoOfElements();
      _type         = goodPixel->getSparsePixelType();
      _trackerData  = data;
    } else {
    // the template class is not inheriting from
    // EUTelBaseSparsePixel, so it cannot be used, throw an
    // exception...
      throw InvalidParameterException(std::string("The template parameter is not valid"));
    }
  } 
  
  template<class PixelType>
  void EUTelSparseDataImpl<PixelType>::addSparsePixel(PixelType * pixel) {
    
    if ( _type == kEUTelSimpleSparsePixel ) {
      EUTelSimpleSparsePixel * simplePixelType = dynamic_cast<EUTelSimpleSparsePixel*> ( pixel ) ;
      _trackerData->chargeValues().push_back( static_cast<float> (simplePixelType->getXCoord()) );
      _trackerData->chargeValues().push_back( static_cast<float> (simplePixelType->getYCoord()) );
      _trackerData->chargeValues().push_back( static_cast<float> (pixel->getSignal()) );
    } else if ( _type == kUnknownPixelType ) {
      throw UnknownDataTypeException("Unknown sparse pixel type");
    }
    
  }

  template<class PixelType>
  unsigned int EUTelSparseDataImpl<PixelType>::size() const {
    return _trackerData->chargeValues().size() / _nElement;
  }

  
//   template<class PixelType> 
//   PixelType * EUTelSparseDataImpl<PixelType>::getSparsePixelAt(unsigned int index, PixelType *) { return 0x0; }	

  template<class PixelType>
  std::list<std::list< unsigned int> > EUTelSparseDataImpl<PixelType>::findClusters(double minDistance) const {

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
	status[iPixel] = iPixel;
	
	// first of all, find all the pixels close to the current one
	// and add them to the groupedPixel vector.
	for ( unsigned int iPixel2 = iPixel ; iPixel2 < size() ; iPixel2++ ) {
	  if ( status[iPixel2] == 0 ) {
	    getSparsePixelAt(iPixel2, otherPixel);
	    if ( distance(pixel, otherPixel) <= minDistance ) {
	      groupedPixel.push_back(iPixel2);
	      status[iPixel2] = iPixel;
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
		  status[iPixel2] = iPixel;
		}
	      }
	    }
	    ++firstIter;
	  }
	}
      }
      
      listOfList.push_back( groupedPixel );
      
    }
    

    delete pixel;
    delete otherPixel;
    
    return listOfList;

  }




}


#endif
