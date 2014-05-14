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


  //! TODO: Implementation of the EUTelescope sparse data structure.
  /*!  TODO: Within the EUTelescope framework input data can be provided both
   *  TODO: in non zero suppressed mode, i.e. one ADC signal for each pixel
   *  TODO: in the detector irrespectively of the signal amplitude; or in
   *  TODO: This sparse data container differs from EUTelSparseDataImpl
   *  TODO: because it stores a local copy of all the pixels in the data
   *  TODO: sample. This allows sorting and vector operations resulting in a
   */
  class EUTelTrackerDataInterfacer{

  public:
    //! Default constructor
   // EUTelTrackerDataInterfacer(IMPL::TrackerDataImpl* data);
	
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


  private:
    //! This is the TrackerDataImpl
    /*! This is the object where the sparse data information are
     * collected all together.
     */
   // IMPL::TrackerDataImpl* _trackerData;

    //! Number of elements in the sparse pixel
    /*! This value is initialized in the constructor and taken from
     * the template class.
     */
    //unsigned int _nElement;

    //! Sparse pixel type
    /*! This enumerator value is set in the constructor and taken from
     * the template class.
     */
   // SparsePixelType _type;

    //! Local copy of the pixels
   // mutable std::vector<PixelType > _pixelVec; 
  };

} //namespace
#endif
