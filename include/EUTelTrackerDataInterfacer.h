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
    virtual EUTelBaseSparsePixel* getSparsePixelAt(unsigned int index, EUTelBaseSparsePixel* pixel) const = 0;

    //! Get the number of sparse pixels in the collection
    /*! This utility can be used to know how many pixels are contained
     * in the TrackerData.
     *
     * @return the size of TrackerData measured in sparse
     * pixels.
     */
    virtual unsigned int size() const = 0;

    //! Expose the TrackerDataImpl to the public
    /*! This method is used to allow a direct and public access to the
     * TrackerDataImpl used to collect all the sparse data
     * information.
     *
     * @return The TrackerDataImpl with all the sparse data.
     */
    //IMPL::TrackerDataImpl* trackerData() = 0;
  };
} //namespace
#endif
