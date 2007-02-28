// -*- C++ -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELFFCLUSTERIMPL_H
#define EUTELFFCLUSTERIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>


namespace eutelescope {

  //! Implementation of the a fixed frame cluster for the EUDET telescope.
  /*! Within the Eutelescope environment, a cluster is a group of
   *  neighboring pixels having a signal or signal to noise such that
   *  all the thresholds are passed.
   *  
   *  This is a special implementation of a cluster having a fixed
   *  frame, it means that the cluster shape is known a-priori to be
   *  rectangular, having the seed pixel (the one with the highest
   *  signal) exactly in the cluster center. This requirement is
   *  setting some limitations:
   *
   *  \li To have the seed pixel in the center the sizes of the cluster
   *  edge have to be some odd number of pixels.
   *  
   *  \li To avoid cluster merging the minimum distance between two
   *  seed pixels has to exceed <code>sqrt( pow(xSize,2) +
   *  pow(ySize,2) )</code>
   *
   *  \li Since the cluster shape is a-priori known, to reconstruct
   *  its position only the seed pixel coordinates are needed together
   *  with the sizes along the two axes. This is requiring that all
   *  clusters have to have the same number of pixels (<code>xSize *
   *  ySize</code>), so in the case of seed pixels founds on a
   *  detector edge or close to it, the corresponding cluster has to
   *  be padded with fictitious null signal pixels not belonging to the
   *  detector.
   *
   *  To store such an object a TrackerData object can be used. The
   *  geometrical information are stored into the two CellID fields
   *  using the CellIDEncoder with the
   *  EUTELESCOPE::CLUSTERDEFAULTENCODING encoding. With this encoding
   *  the following information are recorded into each cluster:
   *
   *  \li <b>sensorID</b>: the identification number of the
   *  detector. All the analysis steps are performed simultaneously on
   *  all detectors of the telescope. So for each identified cluster
   *  we need to know at which detector it is actually belonging.
   *
   *  \li <b>clusterID</b>: the identification number of the
   *  cluster. At this analysis/reconstruction stage there is actually
   *  no limitation on the maximum number of clusters found on a
   *  plane. 
   *
   *  \li <b>xSeed</b> and <b>ySeed</b>: are the coordinates (in
   *  pixels) of the seed pixel.
   *
   *  \li <b>xCluSize</b> and <b>yCluSize</b>: are the sizes of the
   *  cluster edges in pixel unit. Those are bound to be odd numbers
   *  to allow the seed pixel to be the cluster center.
   *
   *  This reimplementation is allowing a quick access to the encoded
   *  information along with some useful methods to have direct access
   *  to cluster related quantities like the total charge, the charge
   *  center of gravity and others.
   *
   *  @see EUTELESCOPE::CLUSTERDEFAULTENCODING for the standard
   *  encoding of fixed frame clusters
   *
   *  @todo Test the charge center of mass method.
   *  
   *  @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @Version $Id: EUTelFFClusterImpl.h,v 1.4 2007-02-28 08:14:19 bulgheroni Exp $
   */ 

  class EUTelFFClusterImpl : public IMPL::TrackerDataImpl {

  public:
    //! Default constructor
    EUTelFFClusterImpl();

    //! Destructor
    virtual ~EUTelFFClusterImpl() { /* NOOP */ ; }
    
    //! Get the detector ID
    /*! This method is used to get from the CellID the detector
     *  identification number
     *
     *  @return the detector ID 
     */
    inline int getDetectorID()  const {
      
      int   rhs = 0;
      lcio::long64 mask  = 0x1F;
      lcio::long64 cell0 = static_cast<lcio::long64> (getCellID0());
      return static_cast<int> ( ( cell0 & mask ) >> rhs );
      
    }

    //! Get the cluster identification number
    /*! This method is used to return the cluster identification number
     *  @return the clusterID
     */
    inline int getClusterID() const {
      
      int rhs = 5;
      lcio::long64 mask = 0x1FE0;
      lcio::long64 cell0 = static_cast<lcio::long64> (getCellID0());
      return static_cast<int> ( ( cell0  & mask ) >> rhs );
    }


    //! Get the seed pixel coordinates
    /*! This method is used to obtain the seed pixel coordinate from
     *  the current cluster.
     *
     *  @param xSeed reference to the x coordinate of the seed pixel
     *  @param ySeed reference to the y coordinate of the seed pixel
     */
    inline void getSeedCoord(int& xSeed, int& ySeed) const {

      lcio::long64 cell0 = static_cast<lcio::long64> (getCellID0());
      lcio::long64 cell1 = static_cast<lcio::long64> (getCellID1());

      {
	// first parameter block
	int rhs = 13;
	lcio::long64 mask = (0xFFF << rhs);
	xSeed = static_cast<int>  ( ( cell0 & mask ) >> rhs ) ;
      }

      { 
	// second parameter block
	// being on the field edge, two masks are required

	int          rhs0  =  25;
	int          lhs1  =   7;
	lcio::long64 mask0 = 0xFE000000;
	lcio::long64 mask1 = 0x1F;

	ySeed = static_cast<int> (( (cell0 & mask0) >> rhs0 ) | ( (cell1 & mask1)  << lhs1 ));
      }
    }

    //! Get the cluster size along the two directions
    /*! This method is used to obtain the cluster size along the two
     *  directions.
     *
     *  @param xSize reference to the cluster size along x
     *  @param ySize reference to the cluster size along y
     */
    inline void getClusterSize(int& xSize, int& ySize) const {
      
      lcio::long64 cell1 = static_cast<lcio::long64> (getCellID1());
      
      { 
	// first parameter block
	int rhs = 5;
	lcio::long64 mask = ( 0x1F << rhs );
	
	xSize = static_cast<int> ( ( cell1 & mask ) >> rhs );
      }

      {
	// second parameter block
	int rhs = 10;
	lcio::long64 mask = ( 0x1F << rhs );
	
	ySize = static_cast<int> ( ( cell1 & mask ) >> rhs );
      }
    }

    //! Get cluster quality 
    /*! This method is used to check the quality of the current
     *  clusters. Possible values of cluster qualities are given by
     *  the ClusterQuality enum.
     *
     *  @return the current cluster quality
     */
    inline ClusterQuality getClusterQuality() const {
      lcio::long64 cell1 = static_cast<lcio::long64> (getCellID1());

      int rhs = 15;
      lcio::long64 mask = ( 0x1F << rhs );
      
      return static_cast<ClusterQuality> ( (cell1 & mask) >> rhs ) ;
 
    }

    //! Set the cluster quality flag
    /*! This method is used to apply the cluster quality flag to the
     *  current cluster. It modifies directly the CellID1 content
     *
     *  @param quality It is the cluster quality bit mask using a
     *  eutelescope::ClusterQuality enum
     */ 
    void setClusterQuality(ClusterQuality quality);

    //! Get distance from another cluster
    /*! This method is used to calculate the distance between to
     *  clusters belonging to the same detector plane (having the same
     *  detectorID). It can be used by a processor aming to separate
     *  merging clusters or just flagging those with kClusterMergin
     *  quality attributes. The returned distance is in pixel units
     *  and it is calculated using the seed position and neither the
     *  CoG nor the Eta corrected cluster center. Those corrections,
     *  in fact, are too small to severily affect the cluster
     *  separation.
     * 
     *  @param otherCluster The other cluster for distance calculation
     *  @return The distance in pixel unit from @a this to @a otherCluster
     */
    float getDistance(EUTelFFClusterImpl * otherCluster) const;

    //! Get the cluster external radius
    /*! This method is used to calculate the fixed frame cluster
     *  external circle radius.
     *
     *  @return The external circle radius in pixel unit
     */ 
    float getExternalRadius() const;

    //! Get the total cluster signal
    /*! This method is used to calculate the total charge contained
     *  inside the current cluster
     *  
     *  @return the total cluster charge
     */
    float getTotalCharge() const;

    //! Get the center of gravity shift
    /*! This method is used to calculate the signed distance of the
     *  charge center of gravity from the seed coordinates. With this
     *  method all pixels in the clusters are considered. This is a
     *  very useful method for the estimation of the eta function.
     *  
     *  @param xCoG reference to the CoG shift along x
     *  @param yCoG reference to the CoG shift along y
     */
    void getCenterOfGravityShift(float& xCoG, float& yCoG) const;

    //! Get the center of gravity shift using only a n x m cluster.
    /*! This method is used to calculate the CoG shift, but only
     *  considering pixels belonging to a rectangular <code> n x
     *  m</code> cluster centered around the seed pixel. If this
     *  sub-cluster is greater or equal to the full cluster then
     *  getCenterOfGravityShift(float&, float&) is used.
     *
     *  @param xCoG reference to the CoG shift along x
     *  @param yCoG reference to the CoG shift along y
     *  @param xSize maximum number of pixels along x
     *  @param ySize maximum number of pixels along y
     */ 
    void getCenterOfGravityShift(float& xCoG, float& yCoG, int xSize, int ySize) const ;

    //! Get the center of gravity shift using only @a n pixels with the
    //highest signal.
    /*! This is another method to get the CoG shift from the current
     *  cluster. With this method only @a n pixels with the highest
     *  signal are used. In case, the selected number of pixels is
     *  greater than the cluster size, then
     *  getCenterOfGravityShift(float&, float&) is used
     *
     *  @param xCoG reference to the CoG shift along x
     *  @param yCoG reference to the CoG shift along y
     *  @param n number of pixels with the highest signal
     */ 
    void getCenterOfGravityShift(float& xCoG, float& yCoG, int n) const ;
    

    //! Get the cluster charge CoG
    /*! This method is used to calculate the current cluster charge
     *  center of gravity. Using those number to identify the cluster
     *  center is inducing a certain bias since the probability to
     *  find the cluster center has to be flat within a pixel.
     *
     *  @param xCoG reference to the charge center of gravity along x
     *  @param yCoG reference to the charge center of gravity along y
     */
    void getCenterOfGravity(float& xCoG, float& yCoG) const;


  };
 
}
#endif
