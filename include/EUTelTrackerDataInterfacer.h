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
 *  (obviously the correct EUTelTrackerDataImpl hast to be 
 *  instantiated somehow)
 */
class EUTelTrackerDataInterfacer{
  protected:
	mutable std::vector<std::reference_wrapper<EUTelBaseSparsePixel const>> _refVec;
	mutable bool _refVecValid = false;
	
	virtual void validateRefVec() const = 0;
  public:
    //! Destructor
    virtual ~EUTelTrackerDataInterfacer() {}

    //! Get one of the sparse pixel
    /*! This method is used to get one of the sparse pixel contained
     * into the TrackerData.
     *
     * @param index Index of the sparse pixel within the collection
     * @param pixel A pointer to the output pixel
     * @return A pointer to the output pixel same as @c pixel
     */
    virtual void getSparsePixelAt(size_t index, std::unique_ptr<EUTelBaseSparsePixel> & pixel) const = 0;

    virtual void push_back(EUTelBaseSparsePixel const & pixel) = 0;
    
	auto begin() const -> decltype(this->_refVec.cbegin()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.cbegin();
    }

	auto end() const -> decltype(this->_refVec.cend()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.cend();
    }

	virtual auto at(size_t i) const -> decltype(_refVec.at(i).get()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.at(i).get();
    }
   
	virtual auto operator[](size_t i) const -> decltype(_refVec.operator[](i).get()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.operator[](i).get();
    }
 
    //! Get the number of sparse pixels in the collection
    /*! This utility can be used to know how many pixels are contained
     * in the TrackerData.
     *
     * @return the size of TrackerData measured in sparse
     * pixels.
     */
    virtual auto size() const -> decltype(_refVec.size()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.size();
    }

    virtual auto empty() const -> decltype(_refVec.empty()) {
		if(!_refVecValid) this->validateRefVec();
		return _refVec.empty();
    }

  };
} //namespace
#endif
