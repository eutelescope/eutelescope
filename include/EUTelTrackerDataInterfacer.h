/*
 * This source code is part of the Eutelescope package of Marlin.
 * You are free to use this source files for your own development as
 * long as it stays in a public research context. You are not
 * allowed to use it for commercial purpose. You must put this
 * header with author names in all development based on this file.
 *
 */
#ifndef EUTELTRACKERDATAINTERFACER_HCC
#define EUTELTRACKERDATAINTERFACER_HCC

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelBaseSparsePixel.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/TrackerDataImpl.h>

#include <functional>

namespace eutelescope {

//! Implementation of the EUTelescope raw tracker data structure.
/*! This is the interface class for the raw data access class.
 *  It is used to interface the data stored in a LCIO::TrackerData
 *  via the pixel classes introduced in EUTelescope.
 *  For the implementation see: EUTelTrackerDataInterfacerImpl
 *  Given that different pixel types store a different amount of 
 *  information, the implemantations have to be templated in the pixel type.
 *  Different specializations for various pixel typed have to be 
 *  introduced. EUTelescope cliuster classes are intended to make use of 
 *  this interfacing class, thus it being the centrally managed piece
 *  of code, determined to interfacer LCIO::TrackerData.
 *  This abstract class was introduced to allow for polymorphic usage 
 *  of EUTelTrackerDataInterfacerImpl. It allows to retrieve 
 *  EUTelBaseSparsePixel without actually knowing the pixel type 
 *  (obviously the correct EUTelTrackerDataImpl has to be 
 *  instantiated somehow)
 */
class EUTelTrackerDataInterfacer{
  protected:
	//!	The vector holding std::reference_wrapper to allow to loop over EUTelBaseSparsePixel const &
	/*!	In order to write generic processors, we want to instantiate a EUTelTrackerDataInterfacer 
	 *	with a factory, not caring about the concrete implementation. Such instances can only loop
	 *  over the very basic pixel type - EUTelBaseSparsePixel.
	 *  Since teh concrete implementation will hold a std::vector of the actual pixel types, we
	 *  use a reference to those to loop over them. Since (certain) operations on the actual pixel 
	 *  vector (in the implementing class) may cause a reallocate, we need to carefully track those
	 *  operations and invalidate this vector. 
	 */
	mutable std::vector<std::reference_wrapper<EUTelBaseSparsePixel const>> _refVec;
	//! Flag to store the validity state of the _refVec
	mutable bool _refVecValid = false;
	//! Method to validate _refVec
	/*! This method will be called whenever an operation on the actual pixel vector might cause a
	 *	reallocation. It clears the _refVec and will loop over the pixel vector and insert
	 *	std::reference_wrapper for each element
	 */
	virtual void validateRefVec() const = 0;
  public:
	//! Destructor
	virtual ~EUTelTrackerDataInterfacer() = default;

	//! Returns a const & to the underlying pixel vector
	std::vector<std::reference_wrapper<EUTelBaseSparsePixel const>> const & getPixels() const {
		if(!_refVecValid) this->validateRefVec();
		return _refVec;
	}


	//! Push back a pixel, very similar to STL containers
	/*! Note that the function takes a EUTelBaseSparsePixel reference, since EUTelBaseSparsePixel
	 *	is purely virtual, a reference to the actual pixel type must be used. A 
	 *	dynamic_cast<PixelType const &> will ensure this - and, given it fails, throw an exception. 
	 *	Using the actual implementing class and 
	 *	EUTelTrackerDataInterfacerImpl<PixelType>::push_back(PixelType const & pixel) will not have
	 *	this cast overhead (but is not as generic).
	 */ 
	virtual void push_back(EUTelBaseSparsePixel const & pixel) = 0; //throws std::bad_cast
    
	//! begin() for EUTelTrackerDataInterfacer
	/*! Loops over the _refVec, note that this function is not virtual. The return type of this is a
	 *	std::const_iterator<std::reference_wrapper<EUTelBaseSparsePixel const>>>
	 */
	auto begin() const -> decltype(this->_refVec.cbegin()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.cbegin();
	}

	//! end() for EUTelTrackerDataInterfacer
	/*! See comments in begin()
	 */
	auto end() const -> decltype(this->_refVec.cend()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.cend();
	}

	//! at() behaving very similar to STL container's::at()
	/*! Contrary to operator[], at() has a range check and might throw an exception. This function is
	 *	virtual and will be overriden by the actual implementation. The return type is a 
	 *	EUTelBaseSparsePixel const & and the overriden implementation will return a PixelType const &
	 *	which is allowed as C++ supports covariance in return types
	 */
	virtual auto at(size_t i) const -> decltype(_refVec.at(i).get()) { //throws std::out_of_range
		if(!_refVecValid) this->validateRefVec();
		return _refVec.at(i).get();
	}
   
	//! operator[] behaving very similar to STL container's::operator[]
	/*! See comments in at()
	 */
	virtual auto operator[](size_t i) const -> decltype(_refVec.operator[](i).get()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.operator[](i).get();
	}
 
	//! Get the number of sparse pixels in the collection
	virtual auto size() const -> decltype(_refVec.size()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.size();
	}

	//! Check if no pixels are present in the collection
	virtual auto empty() const -> decltype(_refVec.empty()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.empty();
	}
};
} //namespace
#endif
