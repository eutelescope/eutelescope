/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELHITMAKER_H
#define EUTELHITMAKER_H

#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelUtility.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesLayerLayout.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRunHeader.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace eutelescope {

  //! Hit maker processor
  /*! Beyond this cryptic name there is a simple as important
   *  processor. This is the place were clusters found in the
   *  clustering processor and saved into TrackerPulse objects become
   *  hits having a certain position in the outern space, a covariance
   *  matrix and all the other attributes needed to use this point
   *  into a track fitting algorithm.
   *
   *  This is the place where the cluster center in the local detector
   *  frame of reference is translated into a geometrical position
   *  according to the telescope description provided in the GEAR
   *  environment. For this reason, a cross check on the availability
   *  of a GEAR instance and of a GEAR geometry description XML file
   *  is done during the init phase.
   *
   *  The center of the cluster can be calculate both using a standard
   *  linear procedure, or using the Eta function. In this case the
   *  user must activate the "Apply Eta correction" function and
   *  provide the name of a suitable collection containing Eta
   *  information for each of the detector.
   *
   *  <h4>Frame of reference conversion</h4>
   *
   *  Data coming from the DAQ and from the previous analysis steps
   *  are organized into Tracker(Raw)Data objects with pixel signals
   *  saved into a 1D array in a row by row arrangement. Which is the
   *  first pixel to be readout depends on the hardware. For the
   *  MimoTel, the first read pixel is the one occupying the top left
   *  most position as shown in the picture here below. This means
   *  that the last pixel will be the one on the bottom right
   *  corner. Generally speaking this local frame of reference has to
   *  be converted in the telescope frame of reference that is
   *  right-handed with the z-axis parallel to the beam directions,
   *  the y-axis going from the bottom to the top side and the x-axis
   *  from the right to the left hand side.
   *
   *  \image html frameOfReference.png "Frame of references: the detector on the
   * left and the telescope one on the right."
   *
   *  Moreover, the detector FoR can be rotated with respect to the
   *  telescope one. This rotation in the x-y plane is expressed by
   *  the two 2D vectors  xPointing[2], yPointing[2].
   *
   *  Some examples:
   *
   *  \li The two FoRs are perfectly aligned:
   *     xPointing = {  1 ,  0 }
   *     yPointing = {  0 ,  1 }
   *
   *  \li The detector is oriented as in the figure above with the x
   *  axis horizontal going left to right and the y axis vertical
   *  going top to bottom.
   *     xPointing = { -1 ,  0 }
   *     yPointing = {  0 , -1 }
   *
   *  Among all possible rotations, four are particularly important
   *  because of the way the sensor support can be inserted into the
   *  telescope mechanics. Those are described into the figure here
   *  below:
   *  \image html orientation.png "Four important orientations"
   *
   *  <b>Orientation 1</b> represents the case in which the sensor is
   *  inserted up right in the telescope setup and the beam is hitting
   *  the front surface of the sensor. This configuration corresponds
   *  to:<br>
   *      xPointing = { -1 ,  0 }
   *      yPointing = {  0 , -1 }
   *
   *  <b>Orientation 2</b> represents the case the sensor is inserted
   *  up side down facing the beam. This is coded as:<br>
   *      xPointing = {  1 ,  0 }
   *      yPointing = {  0 ,  1 }
   *
   *  <b>Orientation 3</b> represents the case of the sensor in up
   *  right position but back illuminated. This is coded as:<br>
   *      xPointing = {  1 ,  0 }
   *      yPointing = {  0 , -1 }
   *
   *  <b>Orientation 4</b> represents the sensor up side down and back
   *  illuminated. This is obtained with:<br>
   *      xPointing = { -1 ,  0 }
   *      yPointing = {  0 ,  1 }
   *
   *  Other possible rotations are not currently available in the
   *  geometry description.
   *
   *  <h4>Control histograms</h4>
   *  If MARLIN_USE_AIDA is defined and a AIDAProcessor is activated,
   *  then some control histrograms are filled. There are mainly three
   *  families of histograms:
   *
   *  @li <b>Local hit map</b>: one histogram per plane, it contains
   *  the hit position in the detector frame of reference. The x and y
   *  axes have been already converted in millimeter.
   *
   *  @li <b>Telescope hit map</b>: one histogram per plane, it
   *  contains the hit position in the telescope frame of
   *  reference. This histo is actually the same of the previous one,
   *  but all translation and rotation to move the local frame of
   *  reference to the telescope one have been performed already.
   *
   *  @li <b>Density plot</b>: one global histogram, it contains into
   *  a 3D frame all detected hit. If the run contains enough hits and
   *  they are uniformly distributed on the sensor surface, this plot
   *  should resemble the telescope geometry.
   *
   *  @li <b>Cluster center</b>: a 2D histogram containing the center
   *  of gravity within a pixel.
   *
   *  @li <b>Cluster center Eta</b>: as above but eta corrected.
   *
   *  @li <b>Cluster center X/Y</b>: projection along X and Y of the
   *  cluster center 2D histogram
   *
   *  @li <b>Cluster center X/Y eta</b>: as above but eta corrected.
   *
   *  <h4>GEAR Geometry file</h4>
   *  The geometry information about plane positions and orientations
   *  is
   *
   *  <h4>Good programming note</h4>
   *
   *  This processor is useful if and only if it has access to
   *  geometry information through a GEAR implementation. Not all
   *  EUTelescope users will build Marlin using GEAR, especially who
   *  is interested to perform only single detector
   *  characterization. So to avoid compilation error and complains,
   *  this class is enclosed within a
   *  \code
   *  \#ifdef USE_GEAR
   *     ....
   *  \#endif
   *  \endcode  block.
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Tracker pulse</b>: A collection with cluster
   *  information. Along with that also the original TrackerData
   *  (suitably reimplemented as a EUTelVirtualCluster) should be
   *  available in the event.
   *
   *  <b>Eta functions</b>: Optionally if the user wants to apply the
   *  eta correction, a collection of EtaFunctionImpl should be loaded
   *  as condition file before this processor.
   *
   *  <h4>Output collections</h4>
   *  <b>Tracker hit</b>: A collection of TrackerHit, ready to be used
   *  by a track fitting algorithm
   *
   *  @param TrackerPulseCollectionName The name of the tracker pulse
   *  collection.
   *
   *  @param HitCollectionName The name of the output Tracker Hit
   *  collection.
   *
   *  @param EtaCollectionName A vector of strings with the name of
   *  the eta collection along x and y.
   *
   *  @param EtaSwitch A boolean to switch on and off the eta
   *  corrections.
   *
   *  @param CoGAlgorithm The center of gravity is calculated always
   *  using the same algorithm but different results can be obtained
   *  changing the number of pixels considered in the calculation. The
   *  user via the steering file can choose to use the "Full" cluster
   *  or a cluster made by only the first N most significant pixels
   *  ("NPixel") or using a submatrix made by N times M pixels
   *  centered around the seed ("NxMPixel").
   *
   *  @param NPixel Parameter used only when CoGAlgorithm is
   *  "NPixel". This is the number of most significant pixels to be
   *  used for the CoG calculation.
   *
   *  @param NxMPixel Parameter used only when CoGAlgorithm is
   *  "NxMPixel". This vector contains respectively the number of
   *  pixel along x and y to be used.
   *
   *  @param Enable3DHisto This parameter can be use to enable /
   *  disable the filling of the density plot mentioned above. The
   *  reason for this is that such an histogram may require a huge
   *  amount of memory and consequently slowing down the full
   *  processing.
   *
   */

  class EUTelHitMaker : public marlin::Processor {

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelHitMaker)

  public:
    //! Returns a new instance of EUTelHitMaker
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelHitMaker.
     */
    virtual Processor *newProcessor() { return new EUTelHitMaker; }

    //! Default constructor
    EUTelHitMaker();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     *
     *  @param run the LCRunHeader of the current run
     */
    virtual void processRunHeader(LCRunHeader *run);

    //! Called every event
    /*! This is called for each event in the file. Each element of the
     *  pulse collection is scanned and the center of the cluster is
     *  translated into the external frame of reference thanks to the
     *  GEAR geometry description.
     *
     *  The cluster center might be calculated using a standard linear
     *  charge center of gravity algorithm or applying a more
     *  sophisticated non-linear eta function. This behaviour is
     *  steered by the user from the steering file.
     *
     *  @throw UnknownDataTypeException if the cluster type is unknown
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent(LCEvent *evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();

    //! Histogram booking
    /*! Some control histograms are filled during this procedure in
     *  order to be able to perform easy check on the quality of the
     *  output hits and also to understand if the frame of reference
     *  conversion has been properly done. Of course this method is
     *  effectively doing something only in the case MARLIN_USE_AIDA.
     *
     *  @since v00-00-09 the histogram booking is done on-the-fly,
     *  in the meaning the each time a new sensorID is encountered
     *  all the corresponding histograms are booked.
     */
    void bookHistos(int sensorID);

  protected:
    //! TrackerPulse collection name
    /*! This is the name of the collection holding the pulse
     *  information. The other collection containing the original
     *  TrackerData is accessed through the references in the
     *  TrackerPulse object. So no need to get it directly.
     */
    std::string _pulseCollectionName;

    //! TrackerHit collection name
    /*! This is the name the user wants to give to the output hit
     *  collection.
     */
    std::string _hitCollectionName;

    //! Coordinates reference frame switch
    bool _switchLocalCoordinates;
    
    //! Fill histogram switch
    bool _histogramSwitch;

  private:
    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

    //! Set of booked histogram
    /*  This helper set is used by the on-the-fly histogram booking
     *  procedure.
     */
    std::set<int> _alreadyBookedSensorID;

	//! Histogram maps 
    #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    std::map<int, AIDA::IBaseHistogram *> _hitLocalHistos;
    std::map<int, AIDA::IBaseHistogram *> _hitTelescopeHistos;
    #endif
  };

  //! A global instance of the processor
  EUTelHitMaker gEUTelHitMaker;
}
#endif
#endif
