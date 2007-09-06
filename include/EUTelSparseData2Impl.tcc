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
  EUTelSparseData2Impl<PixelType>::EUTelSparseData2Impl(IMPL::TrackerDataImpl * data) {

    std::auto_ptr<PixelType> pixel ( new PixelType );
    _nElement     = pixel->getNoOfElements();
    _type         = pixel->getSparsePixelType();
    _trackerData  = data;

    _pixelVec.clear();
    _isPositionSorted = false;
    _isSignalSorted   = false;
    _isOriginalOrder  = true;

    for ( unsigned int index = 0 ; index <  _trackerData->chargeValues().size() ; index += 3 ) {
      _pixelVec.push_back( PixelType( (short) _trackerData->chargeValues()[ index ],
				      (short) _trackerData->chargeValues()[ index + 1 ],
				      (short) _trackerData->chargeValues()[ index + 2 ] ) );
    }


  } 
  
  template<class PixelType>
  unsigned int EUTelSparseData2Impl<PixelType>::size() const {
    return _pixelVec.size();
  }

  
  template<class PixelType>
  std::list<std::list< unsigned int> > 
  EUTelSparseData2Impl<PixelType>::findNeighborPixels(double /* minDistance  */)   const {


    typedef typename std::vector<PixelType> PixelVector;
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

      streamlog_out ( DEBUG1 ) << "Starting from pixel " << iPixel << std::endl;

      if ( status[iPixel] == 0 ) {
	
	streamlog_out ( DEBUG1 ) << "--> Status good " << std::endl;
	
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

	  streamlog_out ( DEBUG1 ) << "--> X, Y " << xCoord << ", " << yCoord << std::endl;

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
	      streamlog_out ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				      << (*firstYPixel ) << std::endl;
	      
	      // now find the last pixel having y = yCoord - 1
	      lastYPixel = firstYPixel;
	      while ( lastYPixel != currentPixel ) {
		if ( (*lastYPixel).getYCoord() != yTest ) {
		  break;
		}
		++lastYPixel;
	      }
	      streamlog_out ( DEBUG1 ) << "--> Last pixel with Y = " << yTest << " is " << std::endl
				       << *(lastYPixel - 1 ) << std::endl;

	      // now the pixel under test should be found somewhere in
	      // between firstYPixel and lastYPixel
	      streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *1* (" << xTest << ", " << yTest << ")" << std::endl;
	      foundPixel = find_if ( firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
	      
	      if ( foundPixel != lastYPixel ) {
		streamlog_out ( DEBUG1 ) << "--> Found pixel *1* " << std::endl;
		indexTest = foundPixel - pixelBegin;
		if ( status[ indexTest ] == 0 ) {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *1* status is good" << std::endl;
		  status[ indexTest ] = 1;
		  groupedPixel.push_back( indexTest );
		} else {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *1* status is bad" << std::endl;
		}
		
		// since pixel *1* has been found, pixel *2* (if
		// exists) has to be the next one in the ordered
		// matrix. Just increment the iterator and check if it
		// is it!
		++xTest;
		++foundPixel;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *2* (" << xTest << ", " << yTest << ")" << std::endl;
		if ( foundPixel != pixelEnd ) {

		  if ( ( (*foundPixel).getXCoord() == xTest ) && 
		       ( (*foundPixel).getYCoord() == yTest ) ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *2* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *2* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *2* status is bad" << std::endl;
		    }
		    
		    // since pixel *2* has been found, pixel *3* (if
		    // exists) has to be the next one in the ordered
		    // matrix. Just increment the iterator and check if
		    // it is it!

		    ++xTest; 
		    ++foundPixel;
		    streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *3* (" << xTest << ", " << yTest << ")" << std::endl;
		    if ( foundPixel != pixelEnd ) {

		      if ( ( (*foundPixel).getXCoord() == xTest ) && 
			   ( (*foundPixel).getYCoord() == yTest ) ) {
			streamlog_out ( DEBUG1 ) << "--> Found pixel *3* " << std::endl;
			indexTest = foundPixel - pixelBegin;
			if ( status[ indexTest ] == 0 ) {
			  streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is good" << std::endl;
			  status[ indexTest ] = 1;
			  groupedPixel.push_back( indexTest );
			} else {
			  streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is bad" << std::endl;
			}
		      } else {
			streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *3* " << std::endl;
		      }

		    } else {
		      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *3* because end of list!" << std::endl;
		    }

		  } else {
		    // it is not pixel *2*, but it can be pixel
		    // *3*. let's check!
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *2* " << std::endl;

		    ++xTest;
		    streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *3* (" << xTest << ", " << yTest << ")" << std::endl;
		    if ( ( (*foundPixel).getXCoord() == xTest ) && 
			 ( (*foundPixel).getYCoord() == yTest ) ) {
		      streamlog_out ( DEBUG1 ) << "--> Found pixel *3* " << std::endl;
		      indexTest = foundPixel - pixelBegin;
		      if ( status[ indexTest ] == 0 ) {
			streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is good" << std::endl;
			status[ indexTest ] = 1;
			groupedPixel.push_back( indexTest );
		      } else {
			streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is bad" << std::endl;
		      }
		    } 
		  }
		} else {
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *2* because end of list!" << std::endl;
		}
	      } else {
		streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *1* " << std::endl;
		// pixel *1* not found, but pixel *2* can still
		// exist. let's look for it
		++xTest;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *2* (" << xTest << ", " << yTest << ")" << std::endl;
		foundPixel = find_if ( firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
		if ( foundPixel != lastYPixel ) {
		  streamlog_out ( DEBUG1 ) << "--> Found pixel *2* " << std::endl;
		  indexTest = foundPixel - pixelBegin;
		  if ( status[ indexTest ] == 0 ) {
		    streamlog_out ( DEBUG1 ) << "--> Pixel *2* status is good" << std::endl;
		    status[ indexTest ] = 1;
		    groupedPixel.push_back( indexTest );
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> Pixel *2* status is bad" << std::endl;
		  }
		  
		  // since pixel *2* has been found, pixel *3* (if
		  // exists) has to be the next one in the ordered
		  // matrix. Just increment the iterator and check if
		  // it is it! 
		  ++xTest;
		  ++foundPixel;
		  streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *3* (" << xTest << ", " << yTest << ")" << std::endl;

		  if ( foundPixel != pixelEnd ) {

		    if ( ( (*foundPixel).getXCoord() == xTest ) && 
			 ( (*foundPixel).getYCoord() == yTest ) ) {
		      streamlog_out ( DEBUG1 ) << "--> Found pixel *3* " << std::endl;
		      indexTest = foundPixel - pixelBegin;
		      if ( status[ indexTest ] == 0 ) {
			streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is good" << std::endl;
			status[ indexTest ] = 1;
			groupedPixel.push_back( indexTest );
		      } else {
			streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is bad" << std::endl;
		      }
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *3* because end of list!" << std::endl;
		  }

		} else { 
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *2* " << std::endl;
		  // pixel *2* not found, but pixel *3* can still be
		  // there somewhere. let's look for it
		  ++xTest;
		  streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *3* (" << xTest << ", " << yTest << ")" << std::endl;
		  foundPixel = find_if (firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
		  if ( foundPixel != lastYPixel ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *3* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *3* status is bad" << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *3* " << std::endl;
		  }
		}
		
	      }
	      
	      firstRowFound = true;
	    } else {
	      // no pixels with y = yCoord - 1, so the first group is
	      // done! 
	      firstRowFound = false;
	      streamlog_out ( DEBUG1 ) << "--> No pixels found with Y = " << yTest << std::endl;
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

	    streamlog_out ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				     << (*firstYPixel ) << std::endl;

	    if ( firstYPixel == currentPixel ) {
	      // this means that the first pixel of this row is the
	      // starting pixel, so pixel *4* doesn't exist.
	      // look directly for pixel *5* incrementing
	      // firstYPixel by one and also xTest
	      //
	      lastYPixel == currentPixel;
	      ++xTest ;
	      foundPixel = firstYPixel + 1;
	      streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *5* (" << xTest << ", " << yTest << ")" << std::endl;
	      
	      if ( foundPixel != pixelEnd ) {

		if ( ( (*foundPixel).getXCoord() == xTest ) && 
		     ( (*foundPixel).getYCoord() == yTest ) ) {
		  streamlog_out ( DEBUG1 ) << "--> Found pixel *5* " << std::endl;
		  indexTest = foundPixel - pixelBegin;
		  if ( status[ indexTest ] == 0 ) {
		    streamlog_out ( DEBUG1 ) << "--> Pixel *5* status is good" << std::endl;
		    status[ indexTest ] = 1;
		    groupedPixel.push_back( indexTest );
		  }
		  lastYPixel == currentPixel;
		} else {
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* " << std::endl;
		}
	      } else {
		streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* because end of list!" << std::endl;
	      }
	      
	    } else if ( ( firstYPixel != currentPixel ) &&
			( firstYPixel != currentPixel + 1 ) ) {
	      streamlog_out ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				       << (*firstYPixel ) << std::endl;
	      // this meas that there are other pixels on the same
	      // row on the current one but with lower x.	      
	      // we also need to find the last pixel 
	      lastYPixel = firstYPixel;
	      while ( lastYPixel != pixelEnd ) {
		if ( (*lastYPixel).getYCoord() != yTest ) {
		  break;
		}
		++lastYPixel;
	      }
	      streamlog_out ( DEBUG1 ) << "--> Last pixel with Y = " << yTest << " is " << std::endl
				       << *(lastYPixel - 1 ) << std::endl;

	      // Check if *4* exist. Pixel *4* has to be found between
	      // firstYPixel and currentPixel
	      streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *4* (" << xTest << ", " << yTest << ")" << std::endl;
	      foundPixel = find_if ( firstYPixel, currentPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ) ) ;
	      if ( foundPixel != currentPixel ) {
		streamlog_out ( DEBUG1 ) << "--> Found pixel *4* " << std::endl;
		indexTest = foundPixel - pixelBegin;
		if ( status[ indexTest ] == 0 ) {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *4* status is good" << std::endl;
		  status[ indexTest ] = 1;
		  groupedPixel.push_back( indexTest );
		} else {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *4* status is bad" << std::endl;
		}


		// if *4* is found, then the *5* has to be two
		// positions ahead this one.
		foundPixel += 2;
		xTest += 2;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *5* (" << xTest << ", " << yTest << ")" << std::endl;
		if ( foundPixel < pixelEnd ) {
		
		  if ( ( (*foundPixel).getXCoord() == xTest ) && 
		       ( (*foundPixel).getYCoord() == yTest ) ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *5* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *5* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *5* status is bad" << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* " << std::endl;
		  }
		} else {
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* because end of list" << std::endl;
		}

	      } else {
		streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *4* " << std::endl;
		// pixel *4* not found. foundPixel is now pointing to
		// currentPixel. So if pixel *5* exists has to be
		// found one position ahead! 
		xTest += 2;
		++foundPixel;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *5* (" << xTest << ", " << yTest << ")" << std::endl;
		if ( foundPixel != pixelEnd ) {

		  if ( ( (*foundPixel).getXCoord() == xTest ) && 
		       ( (*foundPixel).getYCoord() == yTest ) ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *5* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *5* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *5* status is bad" << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* " << std::endl;
		  }

		} else {
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* because end of list" << std::endl;
		}
	      }
	    } else {
	      // this is a very unlikely case that it is not expected
	      // to happen because at least one pixel (the current
	      // one) has to be found on the second row. So complain
	      // and exit!!!!
	      streamlog_out ( ERROR4 ) << "Fatal error with sparse data re-clustering. Exiting" << std::endl;
	      exit(-1);
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
	      streamlog_out ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				       << (*firstYPixel ) << std::endl;

	      // look for the last pixel having the y = yCoord + 1
	      lastYPixel = firstYPixel;
	      while ( lastYPixel != pixelEnd ) {
		if ( (*lastYPixel).getYCoord() != yTest ) {
		  break;
		}
		++lastYPixel;
	      }
	      streamlog_out ( DEBUG1 ) << "--> Last pixel with Y = " << yTest << " is " << std::endl
				       << *(lastYPixel - 1 ) << std::endl;	      
	      
	      // look for pixel *6*
	      streamlog_out ( DEBUG1 ) <<  "--> Checking the presence of pixel *6* (" << xTest << ", " << yTest << ")" << std::endl;
	      foundPixel = find_if ( firstYPixel, lastYPixel,  EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ) ) ;
	      if ( foundPixel != lastYPixel ) {
		streamlog_out ( DEBUG1 ) << "--> Found pixel *6* " << std::endl;
		indexTest = foundPixel - pixelBegin;
		if ( status[ indexTest ] == 0 ) {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *6* status is good" << std::endl;
		  status[ indexTest ] = 1;
		  groupedPixel.push_back( indexTest );
		} else {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *6* status is bad" << std::endl;
		}
		
		// if *6* is found, then the *7* has to be the next
		// one in the ordered matrix
		++foundPixel;
		++xTest;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *7* (" << xTest << ", " << yTest << ")" << std::endl;
		if ( foundPixel != pixelEnd ) {
		
		  if ( ( (*foundPixel).getXCoord() == xTest ) && 
		       ( (*foundPixel).getYCoord() == yTest ) ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *7* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is bad" << std::endl;
		    }
		    
		    // since pixel *7* has been found, pixel *8* (if
		    // exists) has to be the next one in the ordered
		    // matrix. Just increment the iterator and check it
		    // it is it! 
		    ++xTest;
		    ++foundPixel;
		    streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		    if ( foundPixel != pixelEnd ) {
		      if ( ( (*foundPixel).getXCoord() == xTest ) && 
			   ( (*foundPixel).getYCoord() == yTest ) ) {
			streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
			indexTest = foundPixel - pixelBegin;
			if ( status[ indexTest ] == 0 ) {
			  streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
			  status[ indexTest ] = 1;
			  groupedPixel.push_back( indexTest );
			} else {
			  streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
			}
		      } else {
			streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		      }
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* because end of list! " << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *7* " << std::endl;
		    // it is not pixel *7*, then it might be pixel
		    // *8*. Increment only xTest....
		    ++xTest;
		    streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		    if ( ( (*foundPixel).getXCoord() == xTest ) && 
			 ( (*foundPixel).getYCoord() == yTest ) ) {
		      streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
		      indexTest = foundPixel - pixelBegin;
		      if ( status[ indexTest ] == 0 ) {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
			status[ indexTest ] = 1;
			groupedPixel.push_back( indexTest );
		      } else {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		      }
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		    }
		  } 
		} else {
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *7* because end of list! " << std::endl;
		}
	      } else {
		// pixel *6* not found but pixel *7* can still
		// exist. Let's look for it!
		++xTest;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *7* (" << xTest << ", " << yTest << ")" << std::endl;
		foundPixel = find_if ( firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
		if ( foundPixel != lastYPixel ) {
		  streamlog_out ( DEBUG1 ) << "--> Found pixel *7* " << std::endl;
		  indexTest = foundPixel - pixelBegin;
		  if ( status[ indexTest ] == 0 ) {
		    streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is good" << std::endl;
		    status[ indexTest ] = 1;
		    groupedPixel.push_back( indexTest );
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is bad" << std::endl;
		  }
		  
		  // pixel *7* is found, so pixel *8* either is the
		  // next one or it doesn't exist! 
		  ++xTest;
		  ++foundPixel;
		  streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		  if ( foundPixel != pixelEnd ) {

		    if ( ( (*foundPixel).getXCoord() == xTest ) && 
			 ( (*foundPixel).getYCoord() == yTest ) ) {
		      streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
		      indexTest = foundPixel - pixelBegin;
		      if ( status[ indexTest ] == 0 ) {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
			status[ indexTest ] = 1;
			groupedPixel.push_back( indexTest );
		      } else {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		      }
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* because end of list!" << std::endl;
		  }
		} else {
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *7* " << std::endl;		    
		  // pixel *7* not found, what about *8*?
		  ++xTest;
		  streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		  foundPixel = find_if (firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
		  if ( foundPixel != lastYPixel ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		  }
		}
	      }
	      
	    } else {
	      streamlog_out ( DEBUG1 ) << "No pixels found with Y = " << yTest << std::endl;
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
	    streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *5* (" << xTest << ", " << yTest << ")" << std::endl;
	    foundPixel = currentPixel + 1;
	    
	    if ( foundPixel != pixelEnd ) {
	      
	      if ( ( (*foundPixel).getXCoord() == xTest ) && 
		   ( (*foundPixel).getYCoord() == yTest ) ) {
		streamlog_out ( DEBUG1 ) << "--> Found pixel *5* " << std::endl;
		indexTest = foundPixel - pixelBegin;
		if ( status[ indexTest ] == 0 ) {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *5* status is good" << std::endl;
		  status[ indexTest ] = 1;
		  groupedPixel.push_back( indexTest );
		} else {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *5* status is bad" << std::endl;
		}
	      } else {
		streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* " << std::endl;
	      }
	    } else {
	      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *5* because end of list! " << std::endl;
	    }

	    // now move to the next row
	    xTest = xCoord - 1;
	    yTest = yCoord + 1;
	    
	    // find the first pixel having y = yCoord + 1
	    // this pixel has to be found after the current pixel
	    firstYPixel = find_if ( currentPixel, pixelEnd, EUTelBaseSparsePixel::HasYCoord<PixelType> ( yTest ) ) ;
	    if ( firstYPixel != pixelEnd ) {
	      streamlog_out ( DEBUG1 ) << "--> First pixel with Y = " << yTest << " is " << std::endl
				       << (*firstYPixel ) << std::endl;
	      
	      // now find the last pixel having y = yCoord + 1
	      lastYPixel = firstYPixel;
	      while ( lastYPixel != pixelEnd ) {
		if ( (*lastYPixel).getYCoord() != yTest ) {
		  break;
		}
		++lastYPixel;
	      }
	      streamlog_out ( DEBUG1 ) << "--> Last pixel with Y = " << yTest << " is " << std::endl
				       << *(lastYPixel - 1 ) << std::endl;

	      // now the pixel under test should be found somewhere in
	      // between firstYPixel and lastYPixel
	      streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *6* (" << xTest << ", " << yTest << ")" << std::endl;
	      foundPixel = find_if ( firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
	      if ( foundPixel != lastYPixel ) {
		streamlog_out ( DEBUG1 ) << "--> Found pixel *6* " << std::endl;
		indexTest = foundPixel - pixelBegin;
		if ( status[ indexTest ] == 0 ) {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *6* status is good" << std::endl;
		  status[ indexTest ] = 1;
		  groupedPixel.push_back( indexTest );
		} else {
		  streamlog_out ( DEBUG1 ) << "--> Pixel *6* status is bad" << std::endl;
		}

		// since pixel *6* has been found, pixel *7* if exists is
		// the following one
		++xTest;
		++foundPixel;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *7* (" << xTest << ", " << yTest << ")" << std::endl;
		
		if ( foundPixel != pixelEnd ) {

		  if ( ( (*foundPixel).getXCoord() == xTest ) && 
		       ( (*foundPixel).getYCoord() == yTest ) ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *7* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is bad" << std::endl;
		    }	

		    // since pixel *7* has been found, pixel *8* if
		    // exists is the following one
		    ++xTest;
		    ++foundPixel;
		    streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		    if ( foundPixel != pixelEnd ) {

		      if ( ( (*foundPixel).getXCoord() == xTest ) && 
			   ( (*foundPixel).getYCoord() == yTest ) ) {
			streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
			indexTest = foundPixel - pixelBegin;
			if ( status[ indexTest ] == 0 ) {
			  streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
			  status[ indexTest ] = 1;
			  groupedPixel.push_back( indexTest );
			} else {
			  streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
			}
		      } else {
			streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		      }
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* because end of list" << std::endl;
		    }
		  } else {
		    // it is not pixel *7*, but it can be pixel *8*
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *7* " << std::endl;
		    ++xTest;
		    streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		    if ( ( (*foundPixel).getXCoord() == xTest ) && 
			 ( (*foundPixel).getYCoord() == yTest ) ) {
		      streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
		      indexTest = foundPixel - pixelBegin;
		      if ( status[ indexTest ] == 0 ) {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
			status[ indexTest ] = 1;
			groupedPixel.push_back( indexTest );
		      } else {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		      }
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		    }
		  }
		} else {
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *7* because end of list" << std::endl;
		}
	      } else {
		streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *6* " << std::endl;
		// pixel *6* not found, but pixel *7* can still exist.
		++xTest;
		streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *7* (" << xTest << ", " << yTest << ")" << std::endl;
		foundPixel = find_if ( firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
		if ( foundPixel != lastYPixel ) {
		  streamlog_out ( DEBUG1 ) << "--> Found pixel *7* " << std::endl;
		  indexTest = foundPixel - pixelBegin;
		  if ( status[ indexTest ] == 0 ) {
		    streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is good" << std::endl;
		    status[ indexTest ] = 1;
		    groupedPixel.push_back( indexTest );
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> Pixel *7* status is bad" << std::endl;
		  }

		  // pixel *7* is found, so pixel *8* is the next one
		  ++xTest;
		  ++foundPixel;
		  streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		  if ( foundPixel != pixelEnd ) {

		    if ( ( (*foundPixel).getXCoord() == xTest ) && 
			 ( (*foundPixel).getYCoord() == yTest ) ) {
		      streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
		      indexTest = foundPixel - pixelBegin;
		      if ( status[ indexTest ] == 0 ) {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
			status[ indexTest ] = 1;
			groupedPixel.push_back( indexTest );
		      } else {
			streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		      }
		    }  else {
		      streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* because end of list" << std::endl;
		  }
		} else {
		  // pixel *7* is not found, but pixel *8* can still
		  // be there somewhere.
		  streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *7* " << std::endl;
		  ++xTest;
		  streamlog_out ( DEBUG1 ) << "--> Checking the presence of pixel *8* (" << xTest << ", " << yTest << ")" << std::endl;
		  foundPixel = find_if (firstYPixel, lastYPixel, EUTelBaseSparsePixel::HasXCoord<PixelType> ( xTest ));
		  if ( foundPixel != lastYPixel ) {
		    streamlog_out ( DEBUG1 ) << "--> Found pixel *8* " << std::endl;
		    indexTest = foundPixel - pixelBegin;
		    if ( status[ indexTest ] == 0 ) {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is good" << std::endl;
		      status[ indexTest ] = 1;
		      groupedPixel.push_back( indexTest );
		    } else {
		      streamlog_out ( DEBUG1 ) << "--> Pixel *8* status is bad" << std::endl;
		    }
		  } else {
		    streamlog_out ( DEBUG1 ) << "--> NOT Found pixel *8* " << std::endl;
		  }
		}
	      }
	    } else {
	      // no pixels with y = yCoord + 1, nothing else to do! 
	      streamlog_out ( DEBUG1 ) << "--> No pixels found with Y = " << yTest << std::endl;
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
