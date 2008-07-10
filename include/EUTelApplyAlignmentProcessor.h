// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELAPPLYALIGNMENTPROCESSOR_H
#define EUTELAPPLYALIGNMENTPROCESSOR_H 1

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 

// system includes <>
#include <iostream>
#include <map>

namespace eutelescope {

  //! Apply the alignment correction to the input hit collection
  /*! The main goal of this processor is to apply the alignment
   *  correction to the input hit collection. 
   *
   *  The alignment constants are, as well, contained into collections
   *  of EUTelAlignmentConstant objects and can be loaded from an
   *  external LCIO file or from a DB query. 
   *
   *  <h4>Input collections</h4>
   *  <br><b>InputHit</b>. 
   *  This is a collection of TrackerHit containing all the hits that
   *  need to be corrected for misalignment.
   *  
   *  <br><b>AlignmentConstant</b>. 
   *  This is a collection of EUTelAlignmentConstant. This is a
   *  LCGenericObject  containing all the needed alignment constants
   *  calculated by previous processors.
   *  
   *  <h4>Output</h4>
   *  
   *  <br><b>OutputHit</b>. 
   *  This is a collection of TrackerHit with the correct hits.
   *
   *  @param CorrectionMethod There are actually several different
   *  methods to apply the alignment constants. Here below a list of
   *  available methods:
   *     0. <b>Shifts only</b> No rotations angle will be taken into
   *  account.
   *     1. <b>Rotation first</b> The rotational matrix is applied
   *  before the translation.
   *     2. <b>Translation first</b> The hits are first shifted and
   *  then rotated. 
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelApplyAlignmentProcessor.h,v 1.2 2008-07-10 17:00:02 bulgheroni Exp $
   *
   *
   */

  class EUTelApplyAlignmentProcessor:public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTelApplyAlignmentProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelApplyAlignmentProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelApplyAlignmentProcessor;
    }

    //! Default constructor 
    EUTelApplyAlignmentProcessor ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and reset all needed data
     *  members. 
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
    /*! This is called for each event in the file. A few consistency
     *  checks are done on the first event to see if the alignment
     *  constants and the input collection are compatible.
     * 
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. For the time being there is
     *  nothing to check and do in this slot.
     * 
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check (LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is finished. It
     *  prints only a goodbye message
     */
    virtual void end();


  protected:

    //! Input collection name.
    /*! This is the name of the input hit collection.
     */
    std::string _inputHitCollectionName;

    //! Alignment constant collection name
    /*! This is the name of the collection containing the alignment
     *  constants. This should be the results of the execution of a
     *  EUTelMille and pede.
     */
    std::string _alignmentCollectionName;

    //! Output collection name.
    /*! This is the name of the output hit collection.
     */
    std::string _outputHitCollectionName;

    //! Correction method
    /*! There are actually several different
     *  methods to apply the alignment constants. Here below a list of
     *  available methods:
     *     0. <b>Shifts only</b> No rotations angle will be taken into
     *  account.
     *     1. <b>Rotation first</b> The rotational matrix is applied
     *  before the translation.
     *     2. <b>Translation first</b> The hits are first shifted and
     *  then rotated. 
     */
    int _correctionMethod;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Look Up Table for the sensor ID
    std::map< int, int > _lookUpTable;

  };

  //! A global instance of the processor
  EUTelApplyAlignmentProcessor gEUTelApplyAlignmentProcessor;      

}
#endif
