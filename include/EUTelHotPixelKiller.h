// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELHOTPIXELKILLER
#define EUTELHOTPIXELKILLER 1

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/LCCollectionVec.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <map>


namespace eutelescope {

  //! Processor to mask hot pixels
  /*! This processor is used to keep hot matrix out from the analysis
   *  procedure. This processor is based on the idea that if a pixel
   *  is found to be a part of a cluster to often, probably it is a noisy
   *  pixel and should be removed from the game.
   *
   *  The input status collection has to be writable, so if it is read
   *  in from a file, it has to be copied locally (EUTelCopyPedestalProcessor).
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Status collection </b> the current status collection
   *
   *  <h4>Output collections</h4>
   *
   *  @param NoOfEventPerCycle The number of event to take before
   *  proceeding with the masking
   *  @param MaxAllowedFiringFreq This number [0,1] represents the
   *  maximum allowed firing frequency. Set it to a suitable value
   *  depending on the occupancy.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelHotPixelKiller.h,v 1.3 2008-08-23 12:30:51 bulgheroni Exp $
   *
   */

  class EUTelHotPixelKiller : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelHotPixelKiller
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelHotPixelKiller.
     */
    virtual Processor * newProcessor() {
      return new EUTelHotPixelKiller;
    }

    //! Default constructor
    EUTelHotPixelKiller ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and performs some asserts about
     *  the value of the provided parameters
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run LCRunHeader of the this current run
     *
     *  @throw InvalidParameterException if a paramter is wrongly set
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. If the current @c
     *  evt is flagged to be used for update, then the selected
     *  algorithm wrapper is called
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     *
     *  @throw InvalidParameterException if information in the cellID
     *  are inconsistence
     */
    virtual void processEvent (LCEvent * evt);

    //! Check call back
    /*! This method is called every event just after the processEvent
     *  one. For the time being it is just calling the pixel
     *  monitoring protected method
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void check(LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished. Just printing a good bye message
     */
    virtual void end();


  protected:


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! Histogram with the firing frequency 2D distribution
    static std::string _firing2DHistoName;

    //! Histogram with the firing cumulative 1D distribution
    static std::string _firing1DHistoName;

    //! book histogram method
    void bookAndFillHistos();


#endif

    //! Print the summary
    std::string printSummary() const ;


    //! Status collection name
    /*! Input status collection name. Default value is status.
     */
    std::string _statusCollectionName;

    //! Number of events for update cycle
    int _noOfEventPerCycle;

    //! Maximum allowed firing frequency
    float _maxAllowedFiringFreq;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Number of detector planes in the run
    /*! This is the total number of detector saved into this input
     *  file <br>NOTE: Pedestal, noise and especially common mode
     *  suppression is * done on a detector base. It is retrieved from
     *  a runHeader * parameter.
     */
    int _noOfDetectors;

    //! First pixel along X
    /*! This array of int is used to store the number of the first
     *  pixel along the X direction
     */
    IntVec _minX;

    //! Last pixel along X
    /*! This array of int is used to store the number of the last
     *  pixel along the X direction
     */
    IntVec _maxX;

    //! First pixel along Y
    /*! This array of int is used to store the number of the first
     *  pixel along the Y direction
     */
    IntVec _minY;

    //! Last pixel along Y
    /*! This array of int is used to store the number of the last
     *  pixel along the Y direction
     */
    IntVec _maxY;


    //! Total number of cycle
    int _totalNoOfCycle;

  private:

    //! The current updating cycle number
    unsigned short _iCycle;

    //! A vector with the number of killed pixel
    std::vector< std::vector< unsigned short > > _killedPixelVec;

    //! A vector with the firing frequency value
    std::vector< std::vector< unsigned short > > _firingFreqVec;


  };

  //! A global instance of the processor
  EUTelHotPixelKiller gEUTelHotPixelKiller;

}
#endif
