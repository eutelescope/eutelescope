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
   *  procedure itself.
   *
   *  @see eutelescope::EUTelPedestalUpdateProcessor
   *
   *  The user has to specify the initial pedestal and noise values
   *  for each detector as a steering parameters via the _initPedestal
   *  and _initNoise float vector.
   *
   *  <h4>Input collection:</h4>
   *  None
   *
   *  <h4>Output collection:</h4>
   *  <b>PedestalCollection</b><br>
   *  <b>NoiseCollection</b><br>
   *  <b>StatusCollection</b>
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
   *  @param InitPedestal The initial pedestal value. This is a vector
   *  of float, one value for each detector in the telescope.
   *  @param InitNoise The initial noise value. This is a vector of
   *  float, one value for each detector in the telescope.
   *  @param PedestalCollectionName This is the name of the pedestal
   *  collection.
   *  @param NoiseCollectionName This is the name of the noise collection.
   *  @param StatusCollectionName This is the name of the status collection.
   *
   *  <h4> Steering file example</h4>
   *  Here below is the piece of XML file automatically generated and
   *  concerning this processor.
   *  @verbatim
     <processor name="MyEUTelAutoPedestalNoiseProcessor" type="EUTelAutoPedestalNoiseProcessor">
      <!--EUTelAutoPedestalNoiseProcessor produces initial pedestal / noise / status with user provided values-->
      <!--The initial value of noise (one value for detector)-->
      <!--parameter name="InitNoiseValue" type="FloatVec">1 1 1 1 1 1  </parameter-->
      <!--The initial value of pedestal (one value for detector)-->
      <!--parameter name="InitPedestalValue" type="FloatVec">0 0 0 0 0 0  </parameter-->
      <!--The maximum pixel along x (default 255, one value per detector)-->
      <!--parameter name="MaxXVector" type="IntVec">255 255 255 255 255 255  </parameter-->
      <!--The maximum pixel along y (default 255, one value per detector)-->
      <!--parameter name="MaxYVector" type="IntVec">255 255 255 255 255 255  </parameter-->
      <!--The minimum pixel along x (default 0, one value per detector)-->
      <!--parameter name="MinXVector" type="IntVec">0 0 0 0 0 0  </parameter-->
      <!--The minimum pixel along y (default 0, one value per detector)-->
      <!--parameter name="MinYVector" type="IntVec">0 0 0 0 0 0  </parameter-->
      <!--Noise local collection-->
      <parameter name="NoiseCollectionName" type="string" lcioOutType="TrackerData">noise </parameter>
      <!--Pedestal local collection-->
      <parameter name="PedestalCollectionName" type="string" lcioOutType="TrackerData">pedestal </parameter>
      <!--The sensorID for the generated collection (one per detector)-->
      <!--parameter name="SensorIDVec" type="IntVec">0 1 2 3 4 5  </parameter-->
      <!--Pixel status collection-->
      <parameter name="StatusCollectionName" type="string" lcioOutType="TrackerRawData">status </parameter>
     </processor>
     @endverbatim
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   *
   *
   */

  class EUTelAutoPedestalNoiseProcessor : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelAutoPedestalNoiseProcessor
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
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
     *  is incremented. 
     *
     *  @since Version v00-00-09 the detector number is no more taken from
     *  the run header. We assume that the number of initial values for
     *  the pedestal collection properly represents the number of
     *  detectors in the telescope. Anyway a compatibility check is done
     *  when the pedestal are subtracted. The detector boundaries are
     *  as well taken from parameters and no more from the run
     *  header. In other words, the user has to specify which are the
     *  detector minimum and maximum pixel along both directions. By
     *  default the values for a mimotel sensor are used. In this new
     *  version, the user can assign to each element of the output
     *  collections a sensorID different from the position of the
     *  element in the collection. 
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

    //! The sensorID vector
    /*! This vector contains the sensorID for the generated pede,
     *  noise and status collection.
     *
     *  @since This parameter has been introduced since version
     *  v00-00-09 because previously it was assumed that the output
     *  collection sensorID was numbered from 0 to the number of
     *  sensors. To get rid from all the geometrical information
     *  contained into the run header, now the number of detector is
     *  guessed from the number of elements in the initial pedestal
     *  value, the detector boundaries (min, max, X and Y) are
     *  explicitly declared by the user and the sensorID can as well
     *  be chosen.
     */
    IntVec _sensorIDVec;

  };

  //! A global instance of the processor
  EUTelAutoPedestalNoiseProcessor gEUTelAutoPedestalNoiseProcessor;

}
#endif
