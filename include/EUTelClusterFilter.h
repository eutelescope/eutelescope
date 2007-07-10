// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCLUSTERFILTER_H
#define EUTELCLUSTERFILTER_H 1

// eutelescope includes ".h" 
#include "EUTelROI.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 

// system includes <>
#include <map>
#include <string>
#include <vector>
#include <sstream>

namespace eutelescope {

  class EUTelVirtualCluster;

  //! Cluster filter
  /*! This processor is used during the analysis chain to perform a
   *  further filtering on the available clusters. The main advantage
   *  of having this processor is that, especially when the user is
   *  getting used with a new detector, the cluster search
   *  (EUTelClusteringProcessor) should be executed with loose cuts,
   *  and then when the expected positions of peak signal is known,
   *  the previous results can be reprocessed very quickly just
   *  running this cluster filter.
   *
   *  Being born exactly during a beam test, this processor has to be
   *  considered always work in progress and it will take sometimes
   *  before it will reach its final and stable configuration.
   *
   *  The idea is to let the user select from as much as wide as
   *  possible list of cuts the selection criteria. For the time being
   *  it is particular focused on clusters (be careful still not hits,
   *  they are still in the detector frame of reference). 
   * 
   *  All cuts are declared on a detector basis. This means that if in
   *  the telescope there are 4 sensors and the user wants to cut on
   *  clusters having a total charge in excess than 100 ADC counts, so
   *  she/he has to provide the following in the steering file:
   *
   *  @code
   *  <parameter name="ClusterMinTotalCharge" type="FloatVec"> 100 100 100 100 </parameter>
   *  @endcode
   *
   *  If instead she/he wants to cuts on the charge of a cluster but
   *  considering only the first N most significant pixels, the
   *  steering file should something like this:
   *   
   *  @code
   *  <parameter name="ClusterNMinCharge" type="FloatVec"> 9 75 43 55 87 </parameter>
   *  @endcode
   *
   *  This last parameter will set a cut on the cluster charge
   *  (considering only 9 pixels) with a threshold of 75 ADC on the
   *  first sensor, 43 on the second and so on and so forth. Note that
   *  the first digit represents the number of significant pixels
   *
   *  In a similar way, a selection criterion can be based on the
   *  charge collected by a smaller squared sub cluster made, for
   *  example, by 3 x 3 pixels only.
   *
   *  @code
   *  <parameter name="ClusterNxNMinCharge" type"FloatVec"> 3 90 95 120 </parameter>
   *  @encode
   *
   *  As above the first number if the subcluster size.
   *
   *  The seed pixel can also be used as a cut, again in ADC value and
   *  one per sensor.
   *
   *  The user can select also upon the minimum and the maximum number
   *  of clusters per plane. So for example the following code will
   *  selects only events in each for each sensors there are between 1
   *  and 5 clusters.
   *  
   *  @code
   *  <parameter name="MinNoOfCluster" type="IntVec"> 1 1 1 1  </parameter>
   *  <parameter name="MaxNoOfCluster" type="IntVec"> 5 5 5 5  </parameter>
   *  @endcode
   *
   *  To switch off all those cuts the user can set them to zero or to
   *  a negative value. The system will then disable their
   *  functionality.
   *
   *  Another possible selection criterion is based on the cluster
   *  quality. In this case it is worth to remember that the cluster
   *  quality is defined via a specific enum ClusterQuality reported
   *  below:
   *
   *  @code
   *    enum ClusterQuality {
   *      kGoodCluster       = 0,
   *      kIncompleteCluster = 1L << 0,
   *      kBorderCluster     = 1L << 1,
   *      kMergedCluster     = 1L << 2
   *    };
   *  @endcode
   *
   *  This means that a good cluster corresponds to 0, an incomplete
   *  one to 1, etc. It is possible to sum up different cluster
   *  qualities, for example a merged cluster (4) that is also on the
   *  border (2) will have quality equal to 6.
   *  The idea behind is that the processor will select all clusters
   *  having the quality the users specify.

   *  Another selection can be done of the cluster position. For
   *  example the user may wants to keep only clusters having the
   *  cluster center within a certain range of pixels according to the
   *  following parametrization:
   *
   *  @code
   *  <parameter name="InsideRegion" type="FloatVec"> detID xBotLeft  yBotLeft xTopRight yTopRight  </parameter> 
   *  <parameter name="OutsideRegion" type="FloatVec"> detID xBotLeft  yBotLeft xTopRight yTopRight  </parameter> 
   *  @endcode
   *
   *  <h4> Input collection </h4>
   *  A TrackerPulse collection containing the clusters to be filtered.
   *  
   *  <h4> Output collection </h4>
   *  A TrackerPulse collection containing only the clusters having
   *  passed all the selection criteria.
   *
   *  @param PulseCollectionName Name of the input cluster
   *  collection to be filtered
   *
   *  @param FilteredPulseCollectionName Name of the filtered output
   *  collection. 
   *
   *  @param ClusterMinTotalCharge This is the minimum value allowed
   *  for the cluster total charge in ADC counts. For this parameter,
   *  one value for each detector should be provided. To switch it
   *  off, set all values to 0 or to a negative value.
   *
   *  @param ClusterNMinCharge This is a selection criterion very
   *  similar to ClusterMinTotalCharge but it acts not on the total
   *  charge but on the charge collected by the first N most
   *  significant pixels in the cluster. To specify the selection, the
   *  user has to specify ( noOfDetector + 1 ) parameters. In fact,
   *  the first parameter is the number of significant pixels. To
   *  switch it off is enough to set to zero the first value.
   *
   *  @param ClusterNxNMinCharge This is a selection criterion similar
   *  to the previous one but acting on the charge collected by a
   *  smaller squared cluster. The first number in the vector is the
   *  cluster size. The same rules about the number of vector
   *  components and switching off of ClusterNMinCharge apply also here.
   *
   *  @param SeedMinCharge This is the minimum allowed charge
   *  collected by the seed pixel (i. e. the one with the highest
   *  signal). The user has to specify one floating value for each
   *  detector. Set everything to zero, or to a negative value to
   *  disable the cut.
   *
   *  @param ClusterQuality This selection is based on the cluster
   *  quality as defined by the eutelescope::ClusterQuality
   *  enumeration. The user has to specify one integer number for each
   *  detector. To disable the selection put all negative numbers.
   *
   *  @param InsideRegion This is a sort of geographical cut. Only
   *  clusters belonging to the defined ROI will be accepted. To
   *  specify the ROI the user has to provide 5 numbers, the first one
   *  being the sensor identification, the other four the coordinates
   *  of the bottom left and top right corner. Multiple ROIs also on
   *  the same sensor can be applied. To switch it off it is enough to
   *  set the sensor ID to a negative number.
   *
   *  @param OutsideRegion This is very similar to InsideRegion, but
   *  it works the other way round. Only cluster outside the ROI are
   *  accepted. 
   *
   *  @param MinClusterPerPlane This selection is working on a sensor
   *  level and no more on a cluster level. Only events where a
   *  certain minimum number of clusters per plane has been found will
   *  be accepted. The user has to specify one value for each plane,
   *  of course setting a 0 is disabling the cut.
   *
   *  @param MaxClusterPerPlane As MinClusterPerPlane above but
   *  working in the other way round. Only events where the total
   *  number of clusters per plane was not exceeding a certain limit
   *  are accepted. The user has to specify one value for each
   *  plane. Setting a negative number is disabling the cut.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelClusterFilter.h,v 1.6 2007-07-10 14:31:15 bulgheroni Exp $
   *
   *
   */

  class EUTelClusterFilter : public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTeleClusterSeparationProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTeleClusterSeparationProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelClusterFilter;
    }

    //! Default constructor 
    /*! This is the place where all the processor parameters are
     *  assigned to their local variable in the class. Those variables
     *  will be crosscheck for consistency in the init() method.
     */ 
    EUTelClusterFilter ();

    //! Called at the job beginning.
    /*! This is printing out all the parameters and also making a
     *  brief summary of which selection criteria are really used,
     *  which are deactivated and their values. 
     *  This is also the place where depending all the on/off switches
     *  are set and cuts are checked for consistency.
     *
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. 
     * 
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is the method. For each event, the input tracker pulse
     *  collection is obtained and a loop on all clusters is
     *  started. There are a few different kinds of selection
     *  criteria:
     *  
     *  @li  <b>Single cluster based</b>: These cuts can be applied on
     *  a single cut and they don't require the knowledge of anything
     *  else. An example is the cluster pulse cut.
     *
     *  @li  <b>Detector based</b>: These cuts have to be applied
     *  taking into account all the hits found on the same
     *  detector. An example is the rejection of events having a too
     *  high number of clusters. This means that all the clusters of a
     *  sensor should be processed first.
     *
     *  @li  <b>Telescope based</b>: These are even more complicated
     *  because they require a global knowledge of the event. There
     *  are no selction criteria of this kind implemented yet.
     *
     *  All active selection criteria are tested against each
     *  cluster. This is because at the end of the job a complete
     *  rejection summary is displayed and the user can evaluate the
     *  strength of each single item.
     * 
     *  @param evt The input LCEvent
     *
     *  @throw UnknownDataTypeException if the cluster type stored in
     *  the TrackerPulse is unknown.
     */
    virtual void processEvent (LCEvent * evt);


    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. It can be used to fill check
     *  plots. For the time being there is nothing to check and do in
     *  this slot.
     * 
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check (LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is finished. A
     *  part from printing a good bye message, it also show the
     *  selection summary.
     *  
     */
    virtual void end();

    //! Check if the total cluster charge is above a certain value
    /*! This is used to select clusters having a total integrated
     *  charge above a certain value. This threshold value is given on
     *  a per detector basis and stored into the
     *  _clusterMinTotalChargeVec.
     *
     *  @param cluster The cluster under test.
     *  @return True if the @c cluster has a charge below its own threshold.
     * 
     */ 
    bool isAboveMinTotalCharge(EUTelVirtualCluster * cluster) const ;

    //! Check if the total cluster charge is below a certain value
    /*! This is used to select clusters having a total integrated
     *  charge below a certain value. This threshold value is given on
     *  a per detector basis and stored into the
     *  _clusterMaxTotalChargeVec.    .
     *
     *  @todo Implement it! 
     *  @param cluster The cluster under test.
     *  @return True if the @c cluster has a charge below its own threshold.
     * 
     */
    bool isBelowMaxTotalCharge(EUTelVirtualCluster * cluster) const { return true; }


    //! Check against the charge collected by N pixels
    /*! This is working in a similar way to the isAboveMinTotalCharge
     *  but is comparing not the total charge but the charge collected
     *  by the first N most significant pixels
     *  
     *  The thresholds are stored into a vector on a detector
     *  basis. The first number is the number of pixels to be
     *  considered.
     * 
     *  @return True if the charge is above threshold
     *  @param cluster The cluster under test.
     */
    bool isAboveNMinCharge(EUTelVirtualCluster * cluster) const;

    //! Check against the charge collected by N x N pixels
    /*! This cut is working on the charge collected by a subframe N x
     *  N pixels wide centered around the seed.
     *
     *  @param cluster The cluster under test.
     *  @return True if the charge is above threshold.
     */
    bool isAboveNxNMinCharge(EUTelVirtualCluster * cluster) const;


    //! Seed pixel cut
    /*! This is used to select clusters having a seed pixel charge
     *  above the specified threshold
     *  
     *  @return True if the seed pixel charge is above threshold
     *  @param cluster The cluster under test.
     */
    bool isAboveMinSeedCharge(EUTelVirtualCluster * cluster) const;

    //! Quality cut
    /*! This is a selection cut based on the cluster quality. Only
     *  clusters having the given clusters are accepted.
     *
     *  The quality is provided into a vector of integer by the
     *  user. To disable this cut put a negative value into the
     *  quality vector.
     *
     *  @return True if the quality is correct
     *  @param cluster The cluster under test.
     */ 
    bool hasQuality(EUTelVirtualCluster * cluster) const;

    //! Minimum cluster number
    /*! This selection criterion can be used to require a minimum
     *  number of cluster for each detector. To switch the cut off,
     *  just put its corresponding value to 0 or any negative
     *  number. 
     *
     *  @return True if there are more than the required cluster per plane
     *  @param clusterVec A vector with the number of clusters per plane
     */ 
    bool areClusterEnough(std::vector<int > clusterVec) const;

    //! Maximum cluster number
    /*! This selection criterion can be used to limit the maximum
     *  number of cluster reconstructed for each plane. To switch the
     *  cut off, set the corresponding value to any negative
     *  number. Be careful this selection criteria has an opposite
     *  logic with respect to the other, i. e. the cut is passed when
     *  the there are <b>NO</b> too many clusters
     * 
     *  @return True if there are TOO MANY clusters
     *  @param clusterVec A vector with the number of clusters per plane
     */ 
    bool areClusterTooMany(std::vector<int > clusterVec) const;

    //! Inside the ROI
    /*! This selection criterion can be used to get only clusters
     *  having the center within a certain ROI.
     *
     *  @return True if the cluster center is inside the ROI
     *  @param cluster The cluster under test.
     *
     */ 
    bool isInsideROI(EUTelVirtualCluster * cluster) const;

    //! Outside the ROI
    /*! This selection criterion can be used to get only clusters
     *  having the center outside a certain ROI.
     *
     *  @return True if the cluster center is outside the ROI
     *  @param cluster The cluster under test.
     *
     */ 
    bool isOutsideROI(EUTelVirtualCluster * cluster) const;

    //! Print the rejection summary
    /*! To better understand which cut is more important, a rejection
     *  counter is kept updated during the processing and at the end
     *  it is printed out;
     *
     *  @return an output stream object to be printed out
     */
    std::string printSummary() const;
 
  protected:

    //! Input pulse collection name.
    /*! This is the name of the input pulse collection name
     */
    std::string _inputPulseCollectionName;

    //! Output pulse collection name.
    /*! This is the name of the output pulse collection name
     */
    std::string _outputPulseCollectionName;

    //! Threshold for the minimum total cluster charge
    /*! This is a vector of the same size as the number of detectors
     *  in the telescope and for each detector there is a float number
     *  representing the minimum allowed total cluster charge.
     * 
     *  This selection is switched off for a sensor when this value is
     *  lesser equal to zero.
     * 
     */ 
    std::vector<float > _minTotalChargeVec;

    //! Thresholds for the N pixel cluster
    /*! This vector contains the thresholds for the minimum allowed
     *  charge collected by a cluster considering only the first N
     *  most significant pixels.
     *
     *  The number of components of this vector should be a integer
     *  multiple of @c _noOfDetectors + 1. This is because the first
     *  digit for each set is the number of pixels in the cluster
     *  while the other @c _noOfDetectors are the different thresholds
     *
     *  To switch it off, just put N = 0.
     *
     */
    std::vector<float > _minNChargeVec;

    //! Thresholds for the N x N pixel cluster.
    /*! This vector contains the thresholds for the minimum allowed
     *  charge collected by a cluster made by the N x N pixels around
     *  the seed. 
     * 
     *  The number of components of this vector should be a integer
     *  multiple of @c _noOfDetectors + 1. This is because the first
     *  digit for each set is the number of pixels in the cluster
     *  while the other @c _noOfDetectors are the different thresholds
     *
     *  To switch it off, just put N = 0.
     */
    std::vector<float > _minNxNChargeVec;

    //! Thresholds for the seed pixel charge
    /*! This vector contains the thresholds for the minimum allowed
     *  seed charge in the cluster.
     * 
     *  The number of components in this vector should be equal to the
     *  number of detectors in the telescope.
     *
     *  To switch it off set the components to 0
     *
     */
    std::vector<float > _minSeedChargeVec;

    //! Qualities for the clusters
    /*! This vector contains the quality information for the
     *  cluster. To switch it off set a negative number.
     *  Those numbers corresponds to the eutelescope::ClusterQuality
     *  enum type.
     */ 
    std::vector<int > _clusterQualityVec;

    //! Minimum number of cluster per plane
    /*! This vector of contains the minimum number of allowed clusters
     *  per plane. To switch it off set it to zero or to any negative
     *  numbers.
     *
     */
    std::vector<int > _minClusterNoVec;

    //! Maximum number of cluster per plane
    /*! This vector of contains the maximum number of allowed clusters
     *  per plane. To switch it off set a negative number.
     *
     */
    std::vector<int > _maxClusterNoVec;    
    
    //! A vector of region of interest
    /*! This is a vector of region of interests and it is built using
     *  the information provided by the user in the steering file. All
     *  clusters should be centered inside those ROI's
     */ 
    std::vector<EUTelROI > _insideROIVec;

    //! A vector of region of interest 
    /*! This is a vector of region of interests and it is built using
     *  the information provided by the user in the steering file. All
     *  clusters should be centered outside those ROI's.
     */ 
    std::vector<EUTelROI > _outsideROIVec;

  private:

    //! A temporary vector for the inside ROI
    /*! The reason for this temporary array is that from the steering
     *  file the user can specify only values and not directly a
     *  EUTelROI. So the user specified values are stored in the
     *  _tempInsideROI and then moved to the _insideROIVec soon after.
     *
     *  @todo Ask Frank why Marlin is crashing if in the
     *  registerProcessorParameter we ask to put the values in
     *  something that is not a class member
     */
    std::vector<float > _tempInsideROI;

    //! A temporary vector for the outside ROI
    /*! The reason for this temporary array is that from the steering
     *  file the user can specify only values and not directly a
     *  EUTelROI. So the user specified values are stored in the
     *  _tempOutsideROI and then moved to the _outsideROIVec soon after.
     *
     *  @todo Ask Frank why Marlin is crashing if in the
     *  registerProcessorParameter we ask to put the values in
     *  something that is not a class member
     */
    std::vector<float > _tempOutsideROI;

    //! Switch for the minimum total cluster charge
    bool _minTotalChargeSwitch;

    //! Switch for the minimum N pixel cluster charge
    bool _minNChargeSwitch;

    //! Switch for the minimum N x N pixel cluster charge
    bool _minNxNChargeSwitch;

    //! Switch for the minimum seed charge
    bool _minSeedChargeSwitch;

    //! Switch for the cluster quality
    bool _clusterQualitySwitch;

    //! Switch for the minimum cluster number
    bool _minClusterNoSwitch;

    //! Switch for the maximum cluster number
    bool _maxClusterNoSwitch;

    //! Switch for the insideROI selection
    bool _insideROISwitch;

    //! Switch for the outsideROI selection
    bool _outsideROISwitch;
    
    //! The number of detectors
    int _noOfDetectors;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Rejection summary map
    mutable std::map<std::string, std::vector<unsigned int > > _rejectionMap;

  public:

    //! Helper predicate class
    /*! The main use of this class is to find a ROI having the same
     *  detector ID as the current cluster. It is used as a predicate
     *  in find_if.
     *
     *  @author Antonio Bulgheroni, INFN  <mailto:antonio.bulgheroni@gmail.com>
     *  @version $Id: EUTelClusterFilter.h,v 1.6 2007-07-10 14:31:15 bulgheroni Exp $
     */
    class HasSameID {
    public:
      //! Default constructor.
      /*! 
       *  @param id is the detector ID we want to check against
       */ 
      HasSameID(int id) : _id(id) {;}

      //! Overload operator()
      /*! This is the operator used by the predicate function in the
       *  find_if algorithm. 
       *
       *  The idea is that I have a list or ROI's and I want to find
       *  the first one having the same detector ID of the cluster I'm
       *  testing. 
       *  
       *  @param roi The ROI under test
       *  @return True if the ROI is located on the same sensor
       *  identified by @c _id
       */
      bool operator()(EUTelROI roi) const { return (roi.getDetectorID() == _id) ; }
      
    private:
      //! The detector ID
      int _id;
    };
  



  };

  //! A global instance of the processor
  EUTelClusterFilter gEUTelClusterFilter;      

}
#endif

