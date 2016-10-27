/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELRAWDATASPARSIFIER_H
#define EUTELRAWDATASPARSIFIER_H 1

// eutelescope includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>

// system includes <>
#include <vector>


namespace eutelescope {

  //! Zero suppression processor
  /*! This processor is emulating the behaviour of the EUDRB when it
   *  is working in ZS mode. It uses as input the same information the
   *  EUTelCalibrationProcessor uses, i.e. the raw data full frame
   *  collections (one for each sensor), and from the pedestal
   *  database, the pedestal, noise and status collection.
   *
   *  The output of this processor is again a TrackerData collection
   *  but filled using sparsified pixels exactly in the same format
   *  the EUDRB would have produced.
   *
   *  The number of parameters is very limited, in fact the user can
   *  only select which kind of sparsified pixel has to be used using
   *  the SparsePixelType enumerator and the sigma cut.
   *
   *
   *  @see EUTelBaseSparsePixel
   *  @see SparsePixelType
   *
   *
   *  the calibration (mainly pedestal value) as previously saved into
   *  a condition file or database.
   *
   *  The main things done by this processors are the following:
   *
   *
   *  <h4>Input collections</h4>
   *  <br><b>RawDataCollection</b>. This is a collection of
   *  TrackerRawData containing all pixel signals as they are readout
   *  by the DAQ system. This collection is coming from the input
   *  slcio file.
   *
   *  <br><b>PedestalCollection</b>. This is a collection of
   *  TrackerData containing all pixel pedestal values as they are
   *  calculated by EUTelPedestalNoiseProcessor. This collection is
   *  coming in from the condition file.
   *
   *  <br><b>NoiseCollection</b>. This is a collection of TrackerData
   *  contaning all pixel noise values as they are calculated by
   *  EUTelPedestalNoiseProcessor. This collection is coming from
   *  the condition file and is used only if the user switched on the
   *  common mode calculation.
   *
   *  <br><b>StatusCollection</b>. This is a collection of
   *  TrackerRawData containing the status of each pixels. This
   *  collection is added to the event by the condition processor and
   *  is used only if common mode flag is active.
   *
   *  <h4>Output</h4>
   *
   *  <br><b>SparseDataCollection</b>. This is a collection of TrackerData
   *  containing sparsified pixel of the type the user selected with
   *  only pixels having signal in excess the sigma value multiplied
   *  by the noise.
   *
   *  @param RawDataCollectionName Name of the input raw data collection
   *
   *  @param PedestalCollectionName Name of the input (condition)
   *  pedestal collection
   *
   *  @param NoiseCollectionName Name of the input (condition) noise
   *  collection
   *
   *  @param StatusCollectionName Name of the input (condition) status
   *  collection
   *
   *  @param SparsePixelType This is the type of pixel the user wants
   *  to use for sparsified pixel storage. The SparsePixelType enum
   *  should be for this.
   *
   *  @param SigmaCutVector This is a vector of float with one
   *  component per plane. Each floating number is used as a
   *  multiplicative factor in the definition of the selection
   *  threshold according to
   *
   *  @code
   *  Threshold[iPixel][iPlane] = Sigma[iPlane] * Noise[iPixel][iPlane]
   *  @endcode
   *
   *  @param SparseDataCollectionName The name of the output
   *  zero suppressed data collection
   *
   *  @warning This processor assumes that the input rawData, the
   *  pedestal, the noise and the status collection are aligned
   *  according to the sensorID. In other words, the i-th element of
   *  all the input collections has the same sensorID. If, this is not
   *  the case, <b>DO NOT USE</b> this processor and ask to the
   *  developers for a better implementation.
   *
   *  //todo Whenever we will have some few minutes to spare, the
   *  reason for the above warning message can be easily removed
   *  preparing a map for each ancillary collection (pede, noise,
   *  status) linking the position of each element and its sensorID.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */

  class EUTelRawDataSparsifier:public marlin::Processor {

  public:


    //! Returns a new instance of EUTelRawDataSparsifier
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelRawDataSparsifier.
     */
    virtual Processor * newProcessor() {
      return new EUTelRawDataSparsifier;
    }

    //! Default constructor
    EUTelRawDataSparsifier ();

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
     *  @since v00-00-09 the total number of detector is no more taken
     *  from the run header but it is guessed from the number of
     *  elements in the input collection.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. In the case this is
     *  the first event, then a cross-check is done on the
     *  compatibility of input raw data and pedestal collections. They
     *  have to contain the same number of detector planes and each of
     *  them should have exactly the same number of pixels. Of course
     *  this check is not exhaustive, but at least should avoid
     *  segmentation fault.
     *
     *  @throw IncompatibleDataSetException in the case the two
     *  collections are found to be incompatible
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
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
    /*! This method is called when the loop on events is finished. It
     *  prints only a goodbye message
     */
    virtual void end();


  protected:

    //! Input collection name.
    /*! For the time being we have just one collection that can be used as input
     */
    std::string _rawDataCollectionName;

    //! Pedestal collection name.
    /*! We have three different output collections that are going to be
     *  saved in the file within the first event possibly so that the
     *  next processor will have access to the pedestal, noise and
     *  status immediately.  We have a collection for the pedestal, one
     *  for the noise and another for the status.  The first two are
     *  self explaining, while few more words are needed for the status
     *  one.  This one is used, for example, to set a bad pixel mask or
     *  to flag pixel as already belonging to an already reconstructed
     *  cluster.
     */
    std::string _pedestalCollectionName;

    //! Noise collection name.
    /*! See EUTelRawDataSparsifier::_pedestalCollectionName for the detailed description
     */
    std::string _noiseCollectionName;

    //! Status collection name.
    /*! See EUTelRawDataSparsifier::_pedestalCollectionName for the detailed description
     */
    std::string _statusCollectionName;

    //! Sparsified data collection name.
    /*! The name of the output sparsified data collection. Usually
     *  simply "data"
     */
    std::string _sparsifiedDataCollectionName;

    //! Type of sparsified pixel
    /*! Which information of the pixel passing the zero suppression
     *  can be stored into different data structure. The user can
     *  select which one in the steering file using the
     *  SparsePixelType enumerator
     */
    int _pixelType;

    //! Sigma cut vector
    /*! One component for plane, it is used for the threshold
     *  calculation.
     */
    std::vector<float > _sigmaCutVec;

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

    //! Number of detector planes in the run
    /*! This is the total number of detector saved into this input
     *  file
     */
    size_t _noOfDetector;
  };

  //! A global instance of the processor
  EUTelRawDataSparsifier gEUTelRawDataSparsifier;

}
#endif
