// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
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

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#ifdef USE_GEAR
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>
#endif

// lcio includes <.h> 
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

// system includes <>
#include <string>
#include <vector>
#include <map>

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
   *  Good programming note: If you have a look at the code, you will
   *  see that all the parts of code related to GEAR are placed within
   *  a @code #ifdef USE_GEAR ...  #endif @endcode even if the use of
   *  such a processor is completly meaning less without the use of
   *  GEAR. This is true, but if someone wants to use EUTelescope just
   *  to calculate pedestals and to find clusters then GEAR is not
   *  needed. Those #ifdef/#endif blocks are assuring no compilation
   *  errors also for those kind of users.
   *
   *  <h4>Input</h4> 
   *  
   *  <b>TrackerPulse</b>: A collection with cluster
   *  information. Along with that also the original TrackerData
   *  (suitably reimplemented as a EUTelVirtualCluster) should be
   *  available in the event.
   *
   *  <b>GEAR geometry description</b>: During the init() phase this
   *  processor will check if a GEAR geometry description is available
   *  and if Marlin is actually configured to use GEAR. During the
   *  processRunHeader the consistency of the given geometry with the
   *  one stored on the cluster input file is check in order not to
   *  apply the wrong geometry description.
   *
   *  <b>EtaFunction</b>: Optionally if the user wants to apply the
   *  eta correction, a collection of EtaFunctionImpl should be loaded
   *  as condition file before this processor.
   *
   *  <h4>Output</h4>
   *  <b>TrackerHit</b>: A collection of TrackerHit, ready to be used
   *  by a track fitting algorithm
   *
   *  @param TrackerPulseCollectionName The name of the tracker pulse colelction
   *  @param EtaCollectionNames A vector of strings with the name of
   *  the eta collection along x and y.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelHitMaker.h,v 1.1 2007-05-25 05:16:55 bulgheroni Exp $
   *
   *
   */

  class EUTelHitMaker : public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTelHitMaker
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelHitMaker.
     */
    virtual Processor * newProcessor() {
      return new EUTelHitMaker;
    }

    //! Default constructor 
    EUTelHitMaker ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     * 
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. Each element of the
     *  pulse collection is scanned and the center of the cluster is
     *  translated into the external frame of reference thanks to the
     *  GEAR geometry description. 
     *
     *  The cluster center might be calculate using a standard linear
     *  charge center of gravity algortihm or applying a more
     *  sophisticated non linear eta function. This behaviour is
     *  regulated by the user from the steering file.
     * 
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished. 
     */
    virtual void end();


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

    //! Eta function collection names
    /*! This is a vector containg the name of the eta collection along
     *  the two directions (say x and y). AS for pedestal information
     *  the Eta function has to be calculated beforehand and then
     *  loaded as a condition collection by a suitable
     *  ConditionProcessor.
     */
    std::vector< std::string > _etaCollectionNames;

    //! Switch to apply eta correction
    /*! Eta correction is very important to obtain the maximum
     *  possibile spatial resolution, but having the possibility to
     *  switch off the correction might be usefull for debug
     *  purposes. 
     *
     *  _etaCorrection == 1 means that the correction will be applied
     */ 
    int _etaCorrection ;

  private:
    
    //! Run number 
    int _iRun;
    
    //! Event number
    int _iEvt;

    //! Conversion ID map. 
    /*! In the data file, each cluster is tagged with a detector ID
     *  identify the sensor it belongs to. In the geometry
     *  description, there are along with the sensors also "passive"
     *  layers and other stuff. Those are identify by a layerindex. So
     *  we need a conversion table to go from the detectorID to the
     *  layerindex.
     */ 
    std::map< int, int > _conversionIdMap;


#ifdef USE_GEAR
    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */ 
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *  
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */ 
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;
#endif
    

  };

  //! A global instance of the processor
  EUTelHitMaker gEUTelHitMaker;      

}
#endif
