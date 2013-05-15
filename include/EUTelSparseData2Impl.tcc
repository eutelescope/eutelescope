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
#ifndef EUTELSPARSEDATA2IMPL_TCC
#define EUTELSPARSEDATA2IMPL_TCC

namespace eutelescope {

  template<class PixelType>
  EUTelSparseData2Impl<PixelType>::EUTelSparseData2Impl(IMPL::TrackerDataImpl * data)
  : _trackerData(data),
    _nElement(),
    _type(),
    _pixelVec(),
    _isPositionSorted(false),
    _isSignalSorted(false),
    _isOriginalOrder(true)
  {

    std::auto_ptr<PixelType> pixel ( new PixelType );
    _nElement     = pixel->getNoOfElements();
    _type         = pixel->getSparsePixelType();
    _trackerData  = data;

    _pixelVec.clear();
    _isPositionSorted = false;
    _isSignalSorted   = false;
    _isOriginalOrder  = true;

	fillPixelVec();
    /*for ( unsigned int index = 0 ; index <  _trackerData->getChargeValues().size() ; index += 3 ) {
      _pixelVec.push_back( PixelType( (short) _trackerData->getChargeValues()[ index ],
				      (short) _trackerData->getChargeValues()[ index + 1 ],
				      (short) _trackerData->getChargeValues()[ index + 2 ] ) );
    }*/


  } 
 /* 
  template<class PixelType>
  EUTelSparseData2Impl<PixelType>::EUTelSparseData2Impl( const EUTelSparseData2Impl & z ) {

  }

  template<class PixelType>
  EUTelSparseData2Impl<PixelType> & EUTelSparseData2Impl<PixelType>::operator=( const EUTelSparseData2Impl & z ) {
    return *this;
  }
*/
  template<class PixelType>
  unsigned int EUTelSparseData2Impl<PixelType>::size() const {
    return _pixelVec.size();
  }

  
  template<class PixelType>
  std::list<std::list< unsigned int> > 
  EUTelSparseData2Impl<PixelType>::findNeighborPixels(double  minSignal )   const {


    typedef typename std::vector<PixelType > PixelVector;
    typedef typename PixelVector::iterator PixelVectorIterator;
    PixelVectorIterator foundPixel;
    PixelVectorIterator firstYPixel;
    PixelVectorIterator lastYPixel;
    PixelVectorIterator pixelBegin   = _pixelVec.begin();
    PixelVectorIterator pixelEnd     = _pixelVec.end();
    PixelVectorIterator currentPixel;

    // as a first thing sort by position 
    if ( ! _isPositionSorted ) sortByPosition() ;
    
    // prepare the return listOfList
    std::list< std::list< unsigned int> > listOfList;    

    // prepare a status vector to avoid double counting. All the
    // vector components are initialized to 0 that is to say good
    // pixel not belonging to any cluster candidate
    std::vector<int > status( size(), 0 );
    
    for ( unsigned int iPixel = 0 ; iPixel < size() ; iPixel++ ) {

      streamlog_out_T ( DEBUG1 ) << "Starting from pixel " << iPixel << std::endl;

      if ( status[iPixel] == 0 ) {
	
	streamlog_out_T ( DEBUG1 ) << "--> Status good " << std::endl;
	
	std::list<unsigned int > groupedPixel;
	groupedPixel.push_back( iPixel );
	
	// mark the current pixel as already belonging to group of
	// neighbour pixels
	status[iPixel] = 1;

	// prepare an iterator for this list
	std::list<unsigned int >::iterator indexIter = groupedPixel.begin();

	// 
	bool isFirstPixelOfTheGroup = true;
	
	// start a loop over all the pixels already in the list
	while ( indexIter != groupedPixel.end() ) {

	  // this is the current pixel ( x, y )
	  int xCoord = _pixelVec[ (*indexIter) ].getXCoord();
	  int yCoord = _pixelVec[ (*indexIter) ].getYCoord();
	  currentPixel = pixelBegin + (*indexIter);
	  int xTest, yTest, indexTest;
	  bool firstRowFound = false;

	  streamlog_out_T ( DEBUG1 ) << "--> X, Y " << xCoord << ", " << yCoord << std::endl;

	  if ( ! isFirstPixelOfTheGroup ) {
	    
	    // for every pixel in the list (but the first) we have to
	    // check if the 8 pixels around are good and present in the
	    // original list. 
	    // 
	    // Those are:
	    // *1* --> ( x - 1, y - 1 )
	    // *2* --> ( x    , y - 1 )
	    // *3* --> ( x + 1, y - 1 )
	    //
	    // *4* --> ( x - 1, y     )
	    // *5* --> ( x + 1, y     )
	    //
	    // *6* --> ( x - 1, y + 1 )
	    // *7* --> ( x    , y + 1 )
	    // *8* --> ( x + 1, y + 1 )
	    //
	    // Those pixels are naturally divided into 3 groups having
	    // the same y coordinate
	    

	    // this is the first group
	    // let's start from *1* ( x - 1, y - 1 )
	    xTest = xCoord - 1;
	    yTest = yCoord - 1;
	    
	    // find the first pixel having y = yCoord - 1.
	    // This pixel has to be found before the current pixel 
	    firstYPixel = find_if ( pixelBegin, currentPixel, EUTelBaseSparsePixel::HasYCoord<PixelType> ( yTest ) );
	    if ( firstYPixel != currentPixel ) {
	      firstRowFound = true;
	      streamlog_out_T ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				       << (*firstYPixel ) << std::endl;
	      
	      // now start scanning the matrix until one of the
	      // following condition is verified
	      // 
	      // a. a pixel with a different y is found
	      // b. all pixels *1*, *2* and *3* have been found
	      // c. both *2* and *3* are found
	      // d. pixel *3* is found
	      // e. the last pixel in the matrix is found
	      
	      lastYPixel = firstYPixel;
	      while ( true ) {

		// checking for the end of the matrix and in case just
		// break the loop
		if ( lastYPixel == pixelEnd ) 
		  break;
		
		// checking if we are at the end of the row (so it is
		// changing y), in case just break the loop
		if ( (*lastYPixel).getYCoord() != yTest )  {
		  break;
		}
		
		// checking if the current iterator is pointing to
		// pixel *1*, continue incrementing the iterator
		if ( (*lastYPixel).getXCoord() == xCoord - 1 ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *1* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() >= minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *1* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *1* status is bad" << std::endl;
		  }
		  ++lastYPixel;

		  // checking for the end of the matrix and in case just
		  // break the loop
		  if ( lastYPixel == pixelEnd ) 
		    break;
		}

		// checking if the current iterator is pointing to
		// pixel *2*, in case continue incrementing the
		// iterator
		if ( (*lastYPixel).getXCoord() == xCoord ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *2* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() >= minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *2* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *2* status is bad" << std::endl;
		  }
		    
		  ++lastYPixel;

		  // checking for the end of the matrix and in case just
		  // break the loop
		  if ( lastYPixel == pixelEnd ) 
		    break;
		}

		// checking if the current iterator is pointing to
		// pixel *3*, in case break the loop since we found
		// already everything we need!
		if ( (*lastYPixel).getXCoord() == xCoord + 1 ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *3* " << std::endl
					     << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() > minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *3* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *3* status is bad" << std::endl;
		  }
		  ++lastYPixel;
		  break;
		}

		// checking if we passed all searched pixels, in case
		// break the loop...
		if ( (*lastYPixel).getXCoord() > xCoord + 1 ) {
		  break;
		}

		++lastYPixel;
	      }
	    } else {
	      firstRowFound = false;
	    }


	    // let's do the second group now!
	    xTest = xCoord - 1;
	    yTest = yCoord ;

	    // if some pixels from the first group have been found,
	    // then the firstYPixel is the lastYPixel of the previous
	    // group, other we have to look for it again
	    if ( firstRowFound )  firstYPixel = lastYPixel;
	    else {
	      firstYPixel = find_if ( _pixelVec.begin(), currentPixel + 1 , EUTelBaseSparsePixel::HasYCoord<PixelType> ( yTest ) );
	    }

	    streamlog_out_T ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				     << (*firstYPixel ) << std::endl;

	    
	    // now start scanning the matrix until one of the
	    // following condition is verified
	    // 
	    // a. a pixel with a different y is found
	    // b. all pixels *4*, *5* have been found
	    // c. pixel *5* is found
	    // d. the last pixel in the matrix is found
	    lastYPixel = firstYPixel;
	    while ( true ) {

	      // checking for the end of the matrix and in case just
	      // break the loop
	      if ( lastYPixel == pixelEnd ) 
		break;
	      
	      // checking if we are at the end of the row (so it is
	      // changing y), in case just break the loop
	      if ( (*lastYPixel).getYCoord() != yTest )  {
		break;
	      }

	      // checking if the current iterator is pointing to
	      // pixel *4*, continue incrementing the iterator
	      if ( (*lastYPixel).getXCoord() == xCoord - 1 ) {
		streamlog_out_T ( DEBUG1 ) << "--> Found pixel *4* " << std::endl
					     << (*lastYPixel ) << std::endl;
		indexTest = lastYPixel - pixelBegin;
		if ( status[indexTest] == 0 ) {
		  if ( (*lastYPixel).getSignal() > minSignal ) {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *4* status is good" << std::endl;
		    status[ indexTest ] = 1;
		    groupedPixel.push_back( indexTest );
		  }
		} else {
		  streamlog_out_T ( DEBUG1 ) << "--> Pixel *4* status is bad" << std::endl;
		}
		++lastYPixel;
		
		// checking for the end of the matrix and in case just
		// break the loop
		if ( lastYPixel == pixelEnd ) 
		  break;
		
	      }

	      // checking if the current iterator is pointing to
	      // pixel *5*, in case break the loop since we found
	      // already everything we need!
	      if ( (*lastYPixel).getXCoord() == xCoord + 1 ) {
		streamlog_out_T ( DEBUG1 ) << "--> Found pixel *5* " << std::endl
					     << (*lastYPixel ) << std::endl;
		indexTest = lastYPixel - pixelBegin;
		if ( status[indexTest] == 0 ) {
		  if ( (*lastYPixel).getSignal() > minSignal ) {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *5* status is good" << std::endl;
		    status[ indexTest ] = 1;
		    groupedPixel.push_back( indexTest );
		  }
		} else {
		  streamlog_out_T ( DEBUG1 ) << "--> Pixel *5* status is bad" << std::endl;
		}
		++lastYPixel;
		break;
	      }
	      
	      // checking if we passed all searched pixels, in case
	      // break the loop...
	      if ( (*lastYPixel).getXCoord() > xCoord + 1 ) {
		break;
	      }	    

	      ++lastYPixel;
	    }

	    
	    // let's move on with the third group of pixels
	    // 
	    // now considering pixel *6*
	    xTest = xCoord - 1;
	    yTest = yCoord + 1;
	    
	    // as first pixel with y = yCoord + 1 we can consider the
	    // lastYPixel of the previous row. This is found for sure
	    // because at least one pixel with y = yCoord exists
	    firstYPixel = lastYPixel;
	    
	    // this may be not the first, so we can improve...
	    firstYPixel = find_if ( firstYPixel, pixelEnd, EUTelBaseSparsePixel::HasYCoord<PixelType> ( yTest ) );

	    if ( firstYPixel != pixelEnd ) {
	      streamlog_out_T ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				       << (*firstYPixel ) << std::endl;



	      // now start scanning the matrix until one of the
	      // following condition is verified
	      // 
	      // a. a pixel with a different y is found
	      // b. all pixels *6*, *7* and *8* have been found
	      // c. both *7* and *8* are found
	      // d. pixel *8* is found
	      // e. the last pixel in the matrix is found
	      
	      lastYPixel = firstYPixel;
	      while ( true ) {
		
		// checking for the end of the matrix and in case just
		// break the loop
		if ( lastYPixel == pixelEnd ) 
		  break;

		// checking if we are at the end of the row (so it is
		// changing y), in case just break the loop
		if ( (*lastYPixel).getYCoord() != yTest )  {
		  break;
		}

		// checking if the current iterator is pointing to
		// pixel *6*, continue incrementing the iterator
		if ( (*lastYPixel).getXCoord() == xCoord - 1 ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *6* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() > minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *6* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *6* status is bad" << std::endl;
		  }
		  ++lastYPixel;

		  // checking for the end of the matrix and in case just
		  // break the loop
		  if ( lastYPixel == pixelEnd ) 
		    break;
		}


		// checking if the current iterator is pointing to
		// pixel *7*, in case continue incrementing the
		// iterator
		if ( (*lastYPixel).getXCoord() == xCoord ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *7* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() > minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *7* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *7* status is bad" << std::endl;
		  }
		  ++lastYPixel;

		  // checking for the end of the matrix and in case just
		  // break the loop
		  if ( lastYPixel == pixelEnd ) 
		    break;
		}

		// checking if the current iterator is pointing to
		// pixel *8*, in case break the loop since we found
		// already everything we need!
		if ( (*lastYPixel).getXCoord() == xCoord + 1 ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *8* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() > minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		  }
		  ++lastYPixel;
		  break;
		}
		
		// checking if we passed all searched pixels, in case
		// break the loop...
		if ( (*lastYPixel).getXCoord() > xCoord + 1 )   break;
		

		++lastYPixel;
		
	      }
	    }	      

	  } else {
	    
	    // for the first pixel of the group we just have to check
	    // four neighbor because the other should have been
	    // already found.
	    // 
	    // Those are:
	    // *5* --> ( x + 1, y     )
	    //
	    // *6* --> ( x - 1, y + 1 )
	    // *7* --> ( x    , y + 1 )
	    // *8* --> ( x + 1, y + 1 )
	    //
	    
	    // looking for pixel *5* is very easy: either this is
	    // currentPixel + 1 or it doens't exist.
	    xTest = xCoord + 1;
	    yTest = yCoord;
	    streamlog_out_T ( DEBUG1 ) << "--> Checking the presence of pixel *5* (" << xTest << ", " << yTest << ")" << std::endl;
	    foundPixel = currentPixel + 1;
	    
	    if ( foundPixel != pixelEnd ) {
	      
	      if ( ( (*foundPixel).getXCoord() == xTest ) && 
		   ( (*foundPixel).getYCoord() == yTest ) ) {
		streamlog_out_T ( DEBUG1 ) << "--> Found pixel *5* " << std::endl
					   << (*foundPixel ) << std::endl;
		indexTest = foundPixel - pixelBegin;
		if ( status[ indexTest ] == 0 ) {
		  if ( (*foundPixel).getSignal() > minSignal ) {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *5* status is good" << std::endl;
		    status[ indexTest ] = 1;
		    groupedPixel.push_back( indexTest );
		  }
		} else {
		  streamlog_out_T ( DEBUG1 ) << "--> Pixel *5* status is bad" << std::endl;
		}
	      } else {
		streamlog_out_T ( DEBUG1 ) << "--> NOT Found pixel *5* " << std::endl;
	      }
	    } else {
	      streamlog_out_T ( DEBUG1 ) << "--> NOT Found pixel *5* because end of list! " << std::endl;
	    }
	    
	    // now move to the next row
	    xTest = xCoord - 1;
	    yTest = yCoord + 1;
	    
	    // find the first pixel having y = yCoord + 1
	    // this pixel has to be found after the current pixel
	    firstYPixel = find_if ( currentPixel, pixelEnd, EUTelBaseSparsePixel::HasYCoord<PixelType> ( yTest ) ) ;
	    if ( firstYPixel != pixelEnd ) {
	      streamlog_out_T ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				       << (*firstYPixel ) << std::endl;
	      
	      
	      // now start scanning the matrix until one of the
	      // following condition is verified
	      // 
	      // a. a pixel with a different y is found
	      // b. all pixels *6*, *7* and *8* have been found
	      // c. both *7* and *8* are found
	      // d. pixel *8* is found
	      // e. the last pixel in the matrix is found
	      
	      lastYPixel = firstYPixel;
	      while ( true ) {
		// checking for the end of the matrix and in case just
		// break the loop
		if ( lastYPixel == pixelEnd ) 
		  break;
		
		// checking if we are at the end of the row (so it is
		// changing y), in case just break the loop
		if ( (*lastYPixel).getYCoord() != yTest )   break;
		
		
		// checking if the current iterator is pointing to
		// pixel *6*, continue incrementing the iterator
		if ( (*lastYPixel).getXCoord() == xCoord - 1 ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *6* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() > minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *6* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *6* status is bad" << std::endl;
		  }
		  ++lastYPixel;
		  
		  // checking for the end of the matrix and in case just
		  // break the loop
		  if ( lastYPixel == pixelEnd ) 
		    break;
		}
		
		// checking if the current iterator is pointing to
		// pixel *7*, in case continue incrementing the
		// iterator
		if ( (*lastYPixel).getXCoord() == xCoord ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *7* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() > minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *7* status is good" << std::endl
						 << (*lastYPixel ) << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *7* status is bad" << std::endl;
		  }
		  ++lastYPixel;
		  
		  // checking for the end of the matrix and in case just
		  // break the loop
		  if ( lastYPixel == pixelEnd ) 
		    break;
		}

		// checking if the current iterator is pointing to
		// pixel *8*, in case break the loop since we found
		// already everything we need!
		if ( (*lastYPixel).getXCoord() == xCoord + 1 ) {
		  streamlog_out_T ( DEBUG1 ) << "--> Found pixel *8* " << std::endl
					   << (*lastYPixel ) << std::endl;
		  indexTest = lastYPixel - pixelBegin;
		  if ( status[indexTest] == 0 ) {
		    if ( (*lastYPixel).getSignal() > minSignal ) {
		      streamlog_out_T ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl
						 << (*lastYPixel ) << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    }
		  } else {
		    streamlog_out_T ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		  }
		  ++lastYPixel;
		  break;
		}
		
		// checking if we passed all searched pixels, in case
		// break the loop...
		if ( (*lastYPixel).getXCoord() > xCoord + 1 )  break;
		
		++lastYPixel;
		
	      }
	      
	      
	    } else {
	      // no pixels with y = yCoord + 1, nothing else to do! 
	      streamlog_out_T ( DEBUG1 ) << "--> No pixels found with Y = " << yTest << std::endl;
	    }
	    
	    isFirstPixelOfTheGroup = false;
	    
	  }
	    
	  ++indexIter;
	}
	  
	listOfList.push_back( groupedPixel );
      }
    }
      
    return listOfList;
    
  }
    
    


}
  

#endif
