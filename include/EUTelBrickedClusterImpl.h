// Version: $Id$
/*
*   This source code is part of the Eutelescope package of Marlin.
*   You are free to use this source files for your own development as
*   long as it stays in a public research context. You are not
*   allowed to use it for commercial purpose. You must put this
*   header with author names in all development based on this file.
*
*/

#ifndef EUTELBRICKEDCLUSTERIMPL_H
#define EUTELBRICKEDCLUSTERIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"

// marlin includes ".h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <iostream>


namespace eutelescope {


    class EUTelBrickedClusterImpl : public EUTelVirtualCluster {

    public:
    //! Default constructor
    EUTelBrickedClusterImpl(TrackerDataImpl * data);

    //! Destructor
    virtual ~EUTelBrickedClusterImpl() { /* NOOP */ ; }

    //! Set the pixel noise values
    /*! This method is used to set the noise values. The Fixed Frame
    *  cluster implementation does not have in the TrackerData the
    *  noise value of the pixels. This method allows to attach at the
    *  current cluster a vector of floats with one entry for each
    *  pixel representing the noise value. The noise values in the
    *  vector must be ordered with the same order the pixel signals
    *  have been loaded into the TrackerData.
    *
    *  Note that those noise values cannot be saved on disk and that
    *  they will be lost when the EUTelBrickedClusterImpl object will be
    *  destroyed.
    *
    *  When this method is called, the _noiseSetSwitch is set to true
    *  and all the other methods calling noise features will be made
    *  available.
    *
    *  @param noiseValues The vector containing the pixel noise.
    */
    void setNoiseValues(std::vector<float > noiseValues) ;

    //! Get the detector ID
    /*! This method is used to get from the CellID the detector
    *  identification number
    *
    *  @return the detector ID
    */
   inline int getDetectorID()  const {
        UTIL::CellIDDecoder<TrackerDataImpl>cellDecoder( EUTELESCOPE::CLUSTERDEFAULTENCODING );
        return cellDecoder(_trackerData)["sensorID"];
        }

    //! Get the cluster central pixel
    /*! Due to the definition of EUTelBrickedClusterImpl the central pixel
    *  is the seed pixel.
    *
    *  @param xCenter reference to the x coordinate of the central
    *  pixel
    *  @param yCenter reference to the y coordinate of the central
    *  pixel
    */
    inline void getCenterCoord(int& xCenter, int& yCenter) const {
        getSeedCoord(xCenter, yCenter);
    }

    //! Get the seed pixel coordinates
    /*! This method is used to obtain the seed pixel coordinate from
    *  the current cluster.
    *
    *  @param xSeed reference to the x coordinate of the seed pixel
    *  @param ySeed reference to the y coordinate of the seed pixel
    */
    inline void getSeedCoord(int& xSeed, int& ySeed) const {
        UTIL::CellIDDecoder<TrackerDataImpl>cellDecoder( EUTELESCOPE::CLUSTERDEFAULTENCODING );
	xSeed = cellDecoder(_trackerData)["xSeed"];
	ySeed = cellDecoder(_trackerData)["ySeed"];
    }

    //! Get the cluster size along the two directions
    /*! This method is used to obtain the cluster size along the two
    *  directions.
    *
    *  @param xSize reference to the cluster size along x
    *  @param ySize reference to the cluster size along y
    */
    inline void getClusterSize(int& xSize, int& ySize) const {
     UTIL::CellIDDecoder<TrackerDataImpl>cellDecoder( EUTELESCOPE::CLUSTERDEFAULTENCODING );
	xSize = cellDecoder(_trackerData)["xCluSize"];
	ySize = cellDecoder(_trackerData)["yCluSize"];
    }

    //! Get cluster quality
    /*! This method is used to check the quality of the current
    *  clusters. Possible values of cluster qualities are given by
    *  the ClusterQuality enum.
    *
    *  @return the current cluster quality
    */
    inline ClusterQuality getClusterQuality() const {
        UTIL::CellIDDecoder<TrackerDataImpl>cellDecoder( EUTELESCOPE::CLUSTERDEFAULTENCODING );
        return static_cast<ClusterQuality>(static_cast<lcio::long64>(cellDecoder(_trackerData)["quality"]));
        }

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
    float getDistance(EUTelVirtualCluster * otherCluster) const;

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

    //! Get seed charge
    /*! This method is used to get the seed pixel charge, mainly for
    *  histogram filling.
    *
    *  @return the seed pixel charge
    */
    float getSeedCharge() const;

    //! Get the cluster charge with N pixels
    /*! This method is used to obtain the cluster charge if
    *  considering only a certain number of pixels with the highest
    *  signal.
    *
    *  @param nPixel The number of pixels to consider
    *  @return The cluster charge using only nPixel
    */
    float getClusterCharge(int nPixel) const;

    //! Calculate the cluster charge with different number of pixels
    /*! This method is a better and faster replacement of the
    *  getClusterCharge(int) method. This one is actually avoiding to
    *  re-sort the signal vector all the times it is called.
    *
    *  This method is getting the signal array and sort it once in
    *  descending order. Then the sorted array is integrated
    *  according to the number of pixels in the @c nPixels vector.
    *
    *  @param nPixels The list of number of pixels
    *  @return The charges for each number of pixels
    */
    std::vector<float> getClusterCharge(std::vector<int > nPixels) const;

    //! Get the cluster charge within a subframe
    /*! This method is used to obtain the cluster charge considering
    *  only pixels belonging to a subframe of the full cluster.
    *
    *  @param  xSize Odd number to define the x size of the subframe
    *  @param  ySize Odd number to define the y size of the subframe
    *  @return The charge of the cluster subframe
    */
    float getClusterCharge(int xSize, int ySize) const;

    //! Get the center of gravity shift
    /*! This method is used to calculate the signed distance of the
    *  charge center of gravity from the seed coordinates. With this
    *  method all pixels in the clusters are considered. This is a
    *  very useful method for the estimation of the eta function. (SEE BELOW)
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

    //! Get the center of gravity shift using only @a n pixels
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


    //! Get the center of gravity shift FOR ETA:
    /*! For eta we need the cog offset from the original pixel.
    *   -> But we do NOT the correction of the x coordinate!
    *   (seed_XCoordinateCorrection is: in case the seed pixel
    *    lies in an even row, the cog shift will reflect an
    *    additional 0.5 offset correcting the different bricked
    *    pixel geometry compared to a normal pixel matrix)
    *
    *   So this does:
    *    Let the normal COG Shift (with(!) seed coordinate correction)
    *    be computed and adjust the result:
    *    If the seed pixel row is even, then add the 0.5 that was substracted before!
    *
    *  @param xCoG reference to the CoG shift along x
    *  @param yCoG reference to the CoG shift along y
    */
    void getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(float& xCoG, float& yCoG) const ;

    //! Get the center of gravity shift FOR ETA (n most significant pixels):
    /*! Same as above, but only for the n most significant pixels.
    *
    *  @param xCoG reference to the CoG shift along x
    *  @param yCoG reference to the CoG shift along y
    *  @param n number of pixels with the highest signal
    */
    void getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(float& xCoG, float& yCoG, int n) const ;

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

    //! Get the pixel noise values
    /*! This method returns the pixel noise value vector if it was
    *  properly set. Otherwise it throws a DataNotAvailableException.
    *
    *  @return a vector of float with the pixel noise values.
    */
    std::vector<float > getNoiseValues() const;


    //! Get the cluster noise
    /*! This method returns the full cluster noise. This is evaluated
    *  according to the following formula:
    *
    *  @code
    *  N = sqrt( N_1^2 + N_2^2 + ... + N_m^2 )
    *  @endcode
    */
    float getClusterNoise() const;

    //! Get the cluster SNR
    /*! This method is used to calculate the cluster signal to noise
    *  ratio.
    *
    *  @return The cluster SNR for the current cluster
    */
    float getClusterSNR() const ;

    //! Get seed pixel SNR
    /*! This method is used to calculate the seed pixel signal to
    *  noise ratio. In this implementation, the seed pixel is both
    *  the one with the highest signal and the central one.
    *
    *  @return The seed pixel SNR
    */
    float getSeedSNR() const ;

    //! Get the cluster N SNR
    /*! This method returns the SNR of the cluster considering only
    *  the N most significant pixels. The pixel significance is based
    *  on a signal (and not SNR) basis.
    *
    *  @param nPixel The number of pixel to consider in the cluster
    *  @return The cluster N SNR.
    */
    float getClusterSNR(int nPixel) const ;

    //! Calculate the cluster SNR with different number of pixels
    /*! This method is a better and faster replacement of the
    *  getClusterSNR(int) method. This one is actually avoiding to
    *  re-sort the signal vector all the times it is called.
    *
    *  @param nPixels The list of number of pixels
    *  @return The SNRs for each number of pixels
    */
    std::vector<float > getClusterSNR(std::vector<int > nPixels) const ;

    //! Get the cluster N x M SNR
    /*! This method returns the SNR when considering only a
    *  rectangular subframe of N x M pixel centered around the seed
    *  pixel.
    *
    *  @param xSize Odd number to define the x size of the subframe
    *  @param ySize Odd number to define the y size of the subframe
    *  @return The SNR of the cluster subframe
    */
    float getClusterSNR(int xSize, int ySize) const ;

    //! Return a pointer to the TrackerDataImpl
    /*! This method is used to expose to the public the
    *  TrackerDataImpl member.
    *
    *  @return The pointer of _trackerData
    */
    IMPL::TrackerDataImpl * trackerData() { return _trackerData; }

    //! Print method
    /*! This method is used to print out the content of the clusters
    *
    *  @param os The input output stream
    */
    void print(std::ostream& os)  const;

    void debugOutput() const;

    private:

    //! Noise values vector
    std::vector<float > _noiseValues;

    //! Noise set switch
    /*! By default this switch is set to false and it is disabling all
    *  methods referring to the _noiseValue. This bool is set to true
    *  only when the setNoiseValues() is called and all the
    *  crosschecks have been passed.
    */
    bool _noiseSetSwitch;

    //! Removes unwanted signals in a full frame interpreted as bricked!
    /*!
    *  This removes unwanted signals in the linearized matrices
    *  (represented by a vector).
    *
    *  This could be seen as the conversion from a normal
    *  FixedFrame to a Bricked Pixel Frame.
    *
    *  Noise as well as charges should be 'cleaned' by this function
    *  after they were set for the first time! They will be filled with
    *  all the signals of a 3x3 frame, but in a bricked frame two
    *  pixels are not needed! Leaving them in would cause wrong results
    *  especially when computing the noise for a cluster where there are
    *  two signals missing but the noise is available.
    *  That way we can re-use all the functions of a Fixed Frame and
    *  obtain the results of a Cluster that is missing a few entries
    *  because of the special bricked pixel structure assumed.
    */
    void setOutsiderValuesInVectorInterpretedAsBrickedMatrix(std::vector<float>& v, float val) const;

    void outputBricked3x3MatrixFromVector(const std::vector<float>& v) const;

    };
/*
    void getCenterOfGravityShiftBrickedPixelStructure(float& x, float& y)
    {
        //!DUMMY
        printf("fu");
        return;
    }
*/
}
#endif
