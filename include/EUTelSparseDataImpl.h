// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef  EXPERIMENTAL
#ifndef EUTELSPARSEDATAIMPL_H
#define EUTELSPARSEDATAIMPL_H

// personal includes ".h"

// marlin includes ".h"

// lcio includes <.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <vector>


namespace eutelescope {

  // forward declaration
  class EUTelSparsePixel;
  class EUTelClusterImpl;

  //! Implementation of the EUTelescope sparse data structure
  /*! Within the EUTelescope framework input data can be provided both
   *  in non zero suppressed mode, i.e. one ADC signal for each pixel
   *  in the detector irrespectively of the signal amplitude; or in
   *  the so called zero suppress data format. In this case, the DAQ
   *  output to the PC only the information concerning pixels having a
   *  signal above a certain (user defined) threshold. While the
   *  handling of NZS data can be easily done via a
   *  LCIO::TrackerRawDataImpl class, this is no more true for sparse
   *  data for the following reasons:
   *
   *  @li The output for pixels above threshold will be a float number
   *  because it will be already corrected for pedestal and common
   *  noise.
   *
   *  @li Not all pixels are present in the output, so it is
   *  meaningless to use a structure like the LCIO::TrackerRawDataImpl
   *  where the geometrical position of each pixel can be
   *  reconstructed using the pixel readout order.
   *
   *  @li For each pixel storing the analog signal is not enough. At
   *  least also the pixel address has to associated to the
   *  signal. Advanced studies would benefit from having also other
   *  information of the pixel, like the actual pedestal and noise
   *  values, the threshold used for the selection and, last but not
   *  least, the initial raw value.
   *
   *  For all these reasons the use of a standard TrackerRawData is
   *  impossible and another way out has to be found. A possibility
   *  would be to use a LCGenericObject containing all the needed
   *  information but this is adding a performance penalty. Another
   *  possibility is to use a TrackerData as container but within the
   *  adcValues vector instead of storing only the pixel signal also
   *  the pixel coordinates and other things can be stored in the
   *  following manner:
   *
   *  @code
   *  adcValues.push_back(pixelXCoord);
   *  adcValues.push_back(pixelYCoord);
   *  adcValues.push_back(pixelSignal);
   *  adcValues.push_back(threshold);
   *  adcValues.push_back(pedestal);
   *  adcValues.push_back(noise);
   *  adcValues.push_back(rawValue);
   *  adcValues.push_back(spare); // to be properly encoded if needed.
   *  @encode
   *
   *  For the time being the second possibility has been chosen mainly
   *  because the usage of a LCGenericObject is introducing some
   *  performance penalties.
   *
   *  @todo Ask Frank which kind of penalties are coming from the use
   *  of a LCGenericObject, I mean, access to disk, compression
   *  factor or whatever. This is important to understand if those
   *  penalties are really relevant for our purpose. 
   *
   *  To handle the complicated data structure behind this
   *  TrackerData, a specific helper class has been designed on
   *  purpose (eutelescope::EUTelSparsePixel).
   *
   *  An important feature added to this class is the possibility to
   *  return a vector of EUTelClusterImpl as a result of a clusterization
   *  process run on this sparse data. 
   *  
   *  @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @Version $Id: EUTelSparseDataImpl.h,v 1.1 2007-06-11 22:18:04 bulgheroni Exp $
   */ 

  class EUTelSparseDataImpl {

  public:
    //! Default constructor
    EUTelSparseDataImpl(IMPL::TrackerDataImpl * data);
   
    //! Destructor
    virtual ~EUTelSparseDataImpl() { /* NOOP */ ; }
    
    //! Add a sparse pixel
    /*! This method is used to add to the current TrackerDataImpl a
     *  new sparse pixel with all the pieces of information.
     */
    void addSparsePixel(EUTelSparsePixel * pixel);
    
    //! Perform cluster search
    /*! This is a very important method for this class since it is
     *  used to group together nearby pixels present in this
     *  TrackerDataImpl and return a collection of EUTelClusterImpl.
     * 
     *  For the time being only one clustering strategy has been
     *  envisaged for sparse data. This is base on the distance of two
     *  pixels. If this is below the threshold distance set by the
     *  user, than the two pixels are assumed to belong to the same
     *  cluster. Remember that this distance is measure in pixel units
     *  and so two pixels sharing one side are 1 pixel far apart,
     *  while two pixels with a touching corner are
     *  <code>sqrt(2)</code> far apart.
     *
     *  @param minDistance The minimum distance between two pixels in
     *  order to consider the two making a cluster.
     *
     *  @return A vector of pointers to EUTelClusterImpl containing
     *  the cluster search result.
     */
    virtual std::vector<EUTelClusterImpl *> findClusters(double minDistance);

    //! Get one of the sparse pixel
    /*! This method is used to get one of the sparse pixel contained
     *  into the TrackerData. 
     *
     *  @param index Index of the sparse pixel within the collection
     *
     *  @return A EUTelSparsePixel with the information concerning the
     *  pixel number @a index in the collection
     */ 
    EUTelSparsePixel * getSparsePixelAt(int index);
    
    //! Get the number of sparse pixels in the collection
    /*! This utility can be used to know how many pixels are contained
     *  in the TrackerData.
     *
     *  @return the size of TrackerData measured in sparse
     *  pixels. Returns 0x0 in case @a index is out of range.
     */ 
    unsigned int size() const ;

    //! Expose the TrackerDataImpl to the public
    /*! This method is used to allow a direct and public access to the
     *  TrackerDataImpl used to collect all the sparse data
     *  information.
     *  
     *  @return The TrackerDataImpl with all the sparse data.
     */
    IMPL::TrackerDataImpl * trackerData(); 

  private:
    //! This is the TrackerDataImpl
    /*! This is the object where the sparse data information are
     *  collected all together.
     */
    IMPL::TrackerDataImpl * _trackerData;

    //! Static constant.
    /*! This static constant initialized to 8 is used to set the
     *  number of entries of the TrackerData each EUTelSparsePixel is using.
     */ 
    static unsigned const int _nElement;

  };
 
}
#endif
#endif
