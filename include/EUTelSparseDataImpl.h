// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELSPARSEDATAIMPL_H
#define EUTELSPARSEDATAIMPL_H

// personal includes ".h"
#include "EUTelVirtualCluster.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelExceptions.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <vector>
#include <string>
#include <memory>

namespace eutelescope {

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
   *  @Version $Id: EUTelSparseDataImpl.h,v 1.2 2007-08-18 21:49:40 bulgheroni Exp $
   */ 

  template<class PixelType>
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
    void addSparsePixel(PixelType * pixel);
    
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
     *  <code>sqrt(2)</code> far apart. So, setting the minDistance to
     *  sqrt(2) is equivalent to ask "touching pixels".
     *
     *  @param minDistance The minimum distance between two pixels in
     *  order to consider the two making a cluster.
     *
     *  @return A vector of pointers to EUTelClusterImpl containing
     *  the cluster search result.
     */
    virtual std::vector<EUTelVirtualCluster *> findClusters(double minDistance);

    //! Get one of the sparse pixel
    /*! This method is used to get one of the sparse pixel contained
     *  into the TrackerData. 
     *
     *  @param index Index of the sparse pixel within the collection
     *
     *  @return A EUTelSparsePixel with the information concerning the
     *  pixel number @a index in the collection
     */ 
    PixelType * getSparsePixelAt(unsigned int index);
    
    //! Get the number of sparse pixels in the collection
    /*! This utility can be used to know how many pixels are contained
     *  in the TrackerData.
     *
     *  @return the size of TrackerData measured in sparse
     *  pixels. 
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

    //! Number of elements in the sparse pixel
    /*! This value is initialized in the constructor and taken from
     *  the template class.
     */ 
    unsigned int _nElement;

    //! Sparse pixel type
    /*! This enumerator value is set in the constructor and taken from
     *  the template class.
     */
    SparsePixelType _type;

  };
 
  template<class PixelType> 
  EUTelSparseDataImpl<PixelType>::EUTelSparseDataImpl(IMPL::TrackerDataImpl * data) {
    // to work properly the template class has to be or at least
    // inherit from EUTelBaseSparsePixel
    EUTelBaseSparsePixel      * goodPixel;
    std::auto_ptr<PixelType>    pixel(new PixelType);
    if ( goodPixel = dynamic_cast<PixelType *>( pixel.get() ) ) {
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
  std::vector<EUTelVirtualCluster *> EUTelSparseDataImpl<PixelType>::findClusters(double minDistance) {
    
    std::vector<int> status(size, 0);
    for ( unsigned int iPixel = 0 ; iPixel < size() ; iPixel++ ) {

      // grouped pixel is a vector containing the pixels that are
      // fulfilling the minDistance condition
      std::vector<int> groupedPixel;

      // get the pixel under investigation
      EUTelBaseSparsePixel * pixel = getSparsePixelAt(iPixel);

      if ( status[iPixel] == 0 ) {
	groupedPixel.push_back(iPixel);
	status[iPixel] = iPixel;
	
	// first of all, find all the pixels close to the current one
	// and add them to the groupedPixel vector.
	for ( unsigned int iPixel2 = iPixel ; iPixel2 < size() ; iPixel2++ ) {
	  if ( status[iPixel2] == 0 ) {
	    EUTelBaseSparsePixel * otherPixel = getSparsePixelAt(iPixel2);
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
	  std::vector<int>::iterator firstIter = groupedPixel.begin() + 1;
	  while ( firstIter != groupedPixel.end() ) {
	    for ( unsigned int iPixel2 = (*firstIter) + 1 ; iPixel2 < size(); iPixel2++ ) {
	      if ( status[iPixel2] == 0 ) {
		if ( distance( getSparsePixelAt(*firstIter) , getSparsePixelAt(iPixel2) ) <= minDistance ) {
		  groupedPixel.push_back(iPixel2);
		  status[iPixel2] = iPixel;
		}
	      }
	    }
	    ++firstIter;
	  }
	}
      }

      

      // that's the right place to build the cluster

    }

    return 0x0;
  }

  template<>
  EUTelSimpleSparsePixel * EUTelSparseDataImpl<EUTelSimpleSparsePixel>::getSparsePixelAt(unsigned int index) {
    EUTelSimpleSparsePixel * pixel = new EUTelSimpleSparsePixel;
    pixel->setXCoord( static_cast<int> ( _trackerData->chargeValues()[index * _nElement]     ) );
    pixel->setYCoord( static_cast<int> ( _trackerData->chargeValues()[index * _nElement + 1] ) );
    pixel->setSignal( static_cast<int> ( _trackerData->chargeValues()[index * _nElement + 2] ) );
    return pixel;
  }

  template<class PixelType>
  unsigned int EUTelSparseDataImpl<PixelType>::size() const {
    return _trackerData->chargeValues().size() / _nElement;
  }
					 

}
#endif
