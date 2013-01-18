/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCOPYPEDESTALPROCESSOR_H
#define EUTELCOPYPEDESTALPROCESSOR_H 1

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 
#include <IMPL/LCCollectionVec.h>

// system includes <>



namespace eutelescope {

  //! Copy pedestal processor 
  /*! This Marlin processor is used to copy pedestal / noise / status
   *  collection from a condition collection to the current
   *  event. This is an un-avoidable step when the condition
   *  collections are read via SimpleFileHandler because it makes them
   *  read-only. This is, in fact, a limitation all the time we would
   *  like to update the pedestal and noise during the process and
   *  when we would like to use the status to mark hit pixels.
   *
   *  <h4>Input conditions</h4> 
   *  
   *  <b>Pedestal condition</b>: Input pedestal condition.
   *  <br> <b>Noise condition</b>: Input noise condition. <br>
   *  <b>Status condition</b>: Input status condition. <br>
   *
   *  <h4>Output collections</h4>
   *
   *  <b>Pedestal collection</b>: Output pedestal collection.<br>
   *  <b>Noise collection</b>: Output noise collection. <br> <b>Status
   *  Collection</b>: Output status collection.
   *
   *  @param PedestalConditionName This is the name of the pedestal
   *  condition
   *
   *  @param NoiseConditionName This ts the name of the noise
   *  condition
   * 
   *  @param StatusConditionName This is the name of the status
   *  condition
   *
   *  @param PedestalCollectionName This is the name of the pedestal
   *  collection
   *
   *  @param NoiseCollectionName This is the name of the noise
   *  collection
   *
   *  @param StatusCollectionName This is the name of the status collection.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   *
   */

  class EUTelCopyPedestalProcessor : public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTelCopyPedestalProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelCopyPedestalProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelCopyPedestalProcessor;
    }

    //! Default constructor 
    EUTelCopyPedestalProcessor ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters. 
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
    /*! This is called for each event in the file. During the first
     *  event the input conditions are copied to the local
     *  collections. Then the local collections are added to the event
     *  but their ownership is moved to this processor.
     * 
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished. Since the owner of the local copied collections is
     *  this processor, it is its own responsibility to delete
     *  them. This is done in the end() call back.
     */
    virtual void end();


  protected:
    
    //! Pedestal condition name
    /*! Input pedestal condition name. Default value is pedestalDB.
     */
    std::string _pedestalConditionName;

    //! Noise condition name
    /*! Input noise condition name. Default value is noiseDB.
     */
    std::string _noiseConditionName;

    //! Status condition name
    /*! Input status condition name. Default value is statusDB.
     */
    std::string _statusConditionName;

    //! Pedestal collection name
    /*! Input pedestal collection name. Default value is pedestal.
     */
    std::string _pedestalCollectionName;

    //! Noise collection name
    /*! Input noise collection name. Default value is noise.
     */
    std::string _noiseCollectionName;

    //! Status collection name
    /*! Input status collection name. Default value is status.
     */
    std::string _statusCollectionName;
    
    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

  private:

    //! Pedestal collection
    IMPL::LCCollectionVec * _pedestalCollectionVec;

    //! Noise collection
    IMPL::LCCollectionVec * _noiseCollectionVec;

    //! Status collection
    IMPL::LCCollectionVec * _statusCollectionVec;
    
  };

  //! A global instance of the processor
  EUTelCopyPedestalProcessor gEUTelCopyPedestalProcessor;      

}
#endif
