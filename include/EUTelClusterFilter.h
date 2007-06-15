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
   *  central pixel within a certain range of pixels according to the
   *  following parametrization:
   *
   *  @code
   *  <parameter name="InsideRegion" type="IntVec"> detID xBotLeft  yBotLeft xTopRight yTopRight  </parameter> 
   *  <parameter name="OutsideRegion" type="IntVec"> detID xBotLeft  yBotLeft xTopRight yTopRight  </parameter> 
   *  @endcode
   *
   *  @param PulseCollectionName Name of the input cluster
   *  collection to be filtered
   *  
   *  @param FilteredPulseCollectionName Name of the filtered output
   *  collection. 
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelClusterFilter.h,v 1.2 2007-06-15 15:04:07 bulgheroni Exp $
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
    /*  This is the method. For each event, the input tracker pulse
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


    //! Seed pixel cut
    /*! This is used to select clusters having a seed pixel charge
     *  above the specified threshold
     *  
     *  @return True if the seed pixel charge is above threshold
     *  @param cluster The cluster under test.
     */
    bool isAboveMinSeedCharge(EUTelVirtualCluster * cluster) const;

    //! Print the rejection summary
    /*! To better understand which cut is more important, a rejection
     *  counter is kept updated during the processing and at the end
     *  it is printed out;
     *
     *  @return an output stream object to be printed out
     */
    std::stringstream& printSummary() const;
 
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

  private:

    //! Switch for the minimum total cluster charge
    bool _minTotalChargeSwitch;

    //! Switch for the minimum N pixel cluster charge
    bool _minNChargeSwitch;

    //! Switch for the minimum seed charge
    bool _minSeedChargeSwitch;

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

  };

  //! A global instance of the processor
  EUTelClusterFilter gEUTelClusterFilter;      

}
#endif

