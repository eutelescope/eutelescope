// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELAUTOPEDESTALNOISEPROCESSOR
#define EUTELAUTOPEDESTALNOISEPROCESSOR

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 
#include <IMPL/LCCollectionVec.h>

// system includes <>



namespace eutelescope {

  //! Automatic pedestal and noise processor
  /*! Especially when dealing with self-biased MAPS, having a specific
   *  pedestal run can be considered useless since all pixels have a
   *  mean value very close to 0 and the typical noise distribution is
   *  peaked around a single value.
   *  
   *  The general idea behind this process is that, we can start the
   *  analysis procedure with fictitious pedestal and noise values
   *  valid for a whole sensor that are going then updated during the
   *  procedure itself (@see PedestalUpdateProcessor)
   * 
   *  The user has to specify the initial pedestal and noise values
   *  for each detector as a steering parameters via the _initPedestal
   *  and _initNoise float vector.
   *
   *  <h4>Input</h4> 
   *  
   *  <b>InitPedestal</b>: Float vector containing the initial value
   *  for pixel pedestal. The size of this vector should be equal to
   *  the number of detectors in the telescope <br> <b>InitNoise</b>:
   *  Float vector containing the initial value for the pixel
   *  noise. The size of this vector should be equal to the number of
   *  detector in the telescope <br> <b>PedestalCollectionName</b>:
   *  this is the name of the output pedestal collection.<br>
   *  <b>NoiseCollectionName</b>: this is the name of the output noise
   *  collection. <br> <b>StatusCollectionName</b>: this is the name
   *  of the output status collection.
   *
   *  <h4>Output</h4>
   *  <b>PedestalCollection</b><br><b>NoiseCollection</b><br><b>StatusCollection</b>
   *
   *  @param InitPedestal the initial pedestal value (vector of float)
   *  @param InitNoise the initial noise value (vector of float)
   *  @param PedestalCollectionName
   *  @param NoiseCollectionName
   *  @param StatusCollectionName
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelAutoPedestalNoiseProcessor.h,v 1.4 2007-05-21 11:37:33 bulgheroni Exp $
   *
   *
   */

  class EUTelAutoPedestalNoiseProcessor : public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTelAutoPedestalNoiseProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelAutoPedestalNoiseProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelAutoPedestalNoiseProcessor;
    }

    //! Default constructor 
    EUTelAutoPedestalNoiseProcessor ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters. 
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. From the run header the number of detectors is
     *  retrieved and compared with the size of the init values
     *  floating vectors. In the case the input vectors are shorter
     *  than the detector number, they are padded with the last
     *  entry. In the opposite case, the last exceeding elements are
     *  removed.
     * 
     *  The vectors containing the min and max pixels along the two
     *  directions are retrieved as well because they are used to set
     *  the CellID in the CollectionVec.
     * 
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. During the first
     *  event the pedestal / noise / status collections are prepared
     *  with the initialization values provided by the user. Those
     *  collections are taken off from the current event, in order to
     *  have them available to all the events.
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

    //! Initial pedestal value
    /*! This vector of floats is used to store the initial values of
     *  pedestal. The size of this vector should match the number of
     *  detectors in the telescope setup. This vector is provided by
     *  the user from steering file.
     */
    FloatVec _initPedestal;

    //! Initial noise value
    /*! This vector of floats is used to store the initial values of
     *  noise. The size of this vector should match the number of
     *  detectors in the telescope setup. This vector is provided by
     *  the user from steering file.
     */
    FloatVec _initNoise;
    
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
    
  };

  //! A global instance of the processor
  EUTelAutoPedestalNoiseProcessor gEUTelAutoPedestalNoiseProcessor;      

}
#endif
