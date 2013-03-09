/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMAPSDIGI_H
#define EUTELMAPSDIGI_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// Include for Track Detailed Simulation (temporarily in Eutelescope)
#include "TDSPixelsChargeMap.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/TrackerDataImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif


// system includes <>
#include <string>
#include <vector>
#include <map>

#if defined ( USE_ALLPIX )
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"

#include "AllPixEventAction.hh"
#include "G4SDManager.hh"
#include "G4PrimaryVertex.hh"
#include "AllPixMedipix2Digitizer.hh"
#include "AllPixFEI3StandardDigitizer.hh"
#include "AllPixMimosa26Digitizer.hh"
#include "AllPixTimepixDigitizer.hh"
#include "AllPixMCTruthDigitizer.hh"
#include "AllPixLETCalculatorDigitizer.hh"
#endif

namespace eutelescope {

  //! MAPS digitization processor
  /*! Realistic model based on the beam test results is used to
   *  calculate the expected response of MAPS sensors for the Geant4
   *  simulated hits, as returned by Mokka.
   *
   *  First step: conversion from global reference frame to local
   *  reference frame of individual sensor (based on the GEAR geometry
   *  description) was adopted from EUTelHitMaker.
   *
   *  Each Mokka hit is divided into a number of smaller steps (based
   *  on energy deposit and path length). For each step corresponding
   *  charge deposit in sensor pixels is calculated by 2D integration
   *  of the expected charge density over the pixel surface. Charge
   *  capture in the silicon (signal attenuation) and charge
   *  reflection from the epitaxial layer boundary are taken into
   *  account.
   *
   *  The core of the algorithm is implemented in TDS (Track Detailed
   *  Simulation) package by Piotr Niezurawski (pniez@fuw.edu.pl)
   *
   *  @author Aleksander Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
   *  @version $Id$
   *
   */

  class EUTelMAPSdigi : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelMAPSdigi
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelMAPSdigi.
     */
    virtual Processor * newProcessor() {
      return new EUTelMAPSdigi;
    }

    //! Default constructor
    EUTelMAPSdigi ();

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
    /*! This is called for each event in the file.
     *
     *  Each simulated hit is translated from the global to the local
     *  frame of reference  thanks to the GEAR geometry description.
     *
     *  The track segment is then divided into smaller fragments and
     *  the expected charge sharing between pixels is calculated based
     *  on the realistic parametrization.
     *
     *  Finally, all pixels fired are put into the output collection
     *  in the format corresponding to sparcified raw data.
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);

    //! Called if the detector is Mimosa26
    /*! 
     */
    virtual void packMimosa26( int digiIndex );

    //! Called if the detector is FEI4    
    /*! 
     */
    virtual void packFEI4( int digiIndex );


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
     */
    void bookHistos();

    void bookHist1D(std::string HistoName, std::string HistoTitle, int iDet,
                    int xNBin, double xMin, double xMax);

    void bookProf1D(std::string HistoName, std::string HistoTitle, int iDet,
                    int xNBin, double xMin, double xMax, double vMin, double vMax);

    void bookHist2D(std::string HistoName, std::string HistoTitle, int iDet,
                    int xNBin, double xMin, double xMax, int yNBin, double yMin, double yMax);

    //! Histogram filling

    void fillHist1D(std::string HistName, int layerIndex, double xVal, double wVal=1.);

    void fillProf1D(std::string HistName, int layerIndex, double xVal, double wVal);

    void fillHist2D(std::string HistName, int layerIndex, double xVal, double yVal, double wVal=1.);


  protected:

    //
    //  Definitions of the parameters needed by the Track Detailed
    //  Simulation (TDS) algorithm  start here
    //


    //! Maxmimu charge per simulation step [e]
    /*! Maximum charge for single integration step.
     *  If more charge is deposited in a single Mokka step, the step
     *  is divided into smaller ones.
     */
     double  _maxStepInCharge;

    //! Maximum length of single simulation step
    /*! Mokka hits with longer path have to be divided into smaller steps
     */

     double _maxStepInLength;

     //! Maximum range in pixels for charge diffusion calculations
     /*! Charge sharing is calculated for all pixels closer
      * in length and width to the seed pixel than the maximum given
      */

     int _integMaxNumberPixelsAlongL;
     int _integMaxNumberPixelsAlongW;

     //! Attenuation length in silicon [mm]
     /*! Lambda parameter describing charge attenuation in diffusion
      */

     double _chargeAttenuationLength;

    //! Charge reflection coefficient
    /*! Given fraction of the charge, spreading isotropically from the
     *   ionization point is reflected from epitaxial layer boundary.
     */

     double _chargeReflectedContribution;

    //! Number of GSL integration steps

    int _gslFunctionCalls;

     //! Numbers of bins in table storing integration results
     /*! To speed up charge integration results are stored in a 5D
      * table. The more bins the more precise results, but more time
      * is needed to fill the table. Once the table is filled
      * simulation is very fast as no numerical integration has to be
      * made any more. Different numbers of bins can be used in lenght
      * and width of the sensor, and in depth.
      */

    int _integPixelSegmentsAlongL;
    int _integPixelSegmentsAlongW;
    int _integPixelSegmentsAlongH;

    //! Flag for using one integrations storage
    /*! Integration can speed up significantly when using same
     *  integration storage for all sensors. However, this makes sense
     *  only if all sensors have same geometry.
     */

    bool _useCommonIntegrationStorage;

    //! Ionization energy in silicon [eV]
    /*! To account for Poisson charge fluctuations energy deposit
     *   returned by Mokka has to be converted to units of elementary
     *   charge
     */

    double _ionizationEnergy;

    //! Digitized detector ID list
    /*! Different digitization parameters (charge scaling, poisson
     *  fluctuations, gain, noise) can be used for each detector
     *  layer. This vector gives order of detector IDs for following
     *  parameter vectors. If this vector is empty, same values are
     *  used for all planes (first elements of parameter vectors)
     */

    std::vector< int >   _DigiLayerIDs;


    //! Scaling of the deposited charge (collection efficiency)
    /*! Because of additional charge losses, charge collected by the
     *   detector can be smaller than expected from uniform
     *   diffusion. This scaling factor can also be used to decrease or
     *   increase effects of Poisson fluctuations.
     */

    std::vector< float > _depositedChargeScaling;

    //! Poisson smearing flag
    /*! If this flag is set, charge collected on the pixel (number of
     *  elementary charges) is smeared according to the Poisson
     *  distribution. 
     */

    std::vector< int >  _applyPoissonSmearing;


    //! ADC gain in ADC counts per unit charge
    /*! Charge collected on the pixel is amplified and converted to
     *   ADC counts.
     */

    std::vector< float > _adcGain;

    //! ADC gain variation
    /*! ADC gain variation can be used to describe possible readout
     *  fluctuations and other effects which result in effective
     *  "noise" proportional to the signal
     */

    std::vector< float > _adcGainVariation;


    //! ADC noise in ADC counts
    /*! Gaussian readout noise, independent on the signal level, which
     *  is added after charge conversion. 
     */

    std::vector< float > _adcNoise;

    //! ADC offeset in ADC counts
    /*! Constant pedestal value, which is added to all pixels
     *  after charge conversion and before zero suppression.
     */

    std::vector< float > _adcOffset;

    //! Zero Suppression Threshold
    /*! Threshold (in ADC counts) for removing empty pixels from pixel
     *  collection. If set to 0 (default) no zero suppression is applied.
     */

    std::vector< float > _zeroSuppressionThreshold;


    //! Sensor readout type flag
    /*! This flag decides about the final conversion of the calculated
     *  signal (in ADC counts) to the pixel signal stored in an output
     *  collection. Considered readout types are: 
     *    \li  1 - digital: final signal resulting from the applied gain,
     *             noise and offset values is stored as a short int
     *             number (if passed the threshold cut)
     *    \li  2 - binary: value 1 is stored for all fired pixels
     */

    std::vector< int >  _pixelReadoutType;

    //! ADC range
    /*! Maximum value, which can be stored as a pixel signal.
     *  If set, signal is also constrained to be positive or 0.
     *  Used only for digital readout (readout type set to 1).
     */

    std::vector< int > _adcRange;

    //
    // Other processor parameter definitions start here
    //

    //! SimTrackerHit collection name
    /*! This is the name of the collection holding the simulated hits
     *  from Mokka.
     */
    std::string _simhitCollectionName;


    //! Output pixel collection name
    /*! This is the name of the collection used to store the results
     *  of digitization.
     */
    std::string _pixelCollectionName;


    //! Debug Event Count
    /*! Print out debug and information
     *  messages only for one out of given number of events. If
     *  zero, no debug information is printed.
     */
   int _debugCount ;

    //! Flag for filling charge profiles
    /*! Charge profiles present average pixel charge, as a function of
     *  pixel number in cluster. Pixels in cluster are sorted either
     *  by |charge| or by |charge|/distance_to_seed (weighted
     *  profile). Filling profiles can be time consuming as it
     *  requires sorting pixel collection few times.
     */

    bool _fillChargeProfiles;

  private:

    // Local functions used in the algorithm

    //! Checks track segment start and end point
    /*! Track segment start and end points are compared with sensor
     *  dimensions. If track is going outside the segment, length of
     *  the track and position of its center are corrected.
     */

    double CheckPathLimits();

    //! Apply sensor rotation as described in GEAR
    /*! Rotation using Euler angles, as defined in new GEAR geometry
     *  description. Part of transformation from telescope frame to
     *  local sensor frame.
     */

    void InvEulerRotation(double* _telPos, double* _gRotation); 

  
    //
    // Local variables
    //

    TrackerDataImpl * zsFrame;
 
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

    //! Digitization layer ID map.
    /*! This map is used to match input digitization parameters read
     *  from steering file to detector ID. 
     */
    std::map< int, int > _digiIdMap;

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

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    //! Name of the local hit map histo
    /*! The histogram pointed by this name is a 2D histo. The x and y
     *  axes correspond to the pixel detector axes in its own local
     *  frame of reference already converted in millimeter. There is
     *  a histo like this for each detector in the geometry.
     */
    static std::string _hitHistoLocalName;

    //! Name of the hit map histo
    /*! The histogram pointed by this name is a 2D histo. The x and y
     *  axes correspond to the pixel detector axes already in the
     *  telescope frame of reference. There is a histo like this for
     *  each detector in the geometry.
     */
    static std::string _hitHistoTelescopeName;

    //! Name of the pixel map histo
    /*! The histogram pointed by this name is a 2D histo. The x and y
     *  axes correspond to the pixel number in its own local frame. There is
     *  a histo like this for each detector in the geometry.
     */
    static std::string _pixelHistoName;

    //! Name for the collected charge distribution

    static std::string _chargeHistoName;

    //! Name for the total ADC signal distribution

    static std::string _signalHistoName;

    //! Name for the pixel multiplicity distribution

    static std::string _multipHistoName;

    //! Name for the cluster charge profile (ordered in charge)

    static std::string _chargeProfileName;

    //! Name for the cluster charge profile (ordered in
    //! charge/seed_distance ratio)

    static std::string _weightedProfileName;

#endif

    //! Fill histogram switch
    /*! This boolean switch was initially introduced for debug reason
     *  but then we realized that it could stay there and protect
     *  against missing AIDA::Processor.
     *
     */
    bool _histogramSwitch;

    /*!
     *  Following variables contain information about Mokka hit
     *  transformed to the local (sensor) coordinate frame
     */

    //! Mokka hit position in local frame

    double _localPosition[3];

    //! Particle momentum (as stored by Mokka) in local frame

    double _localMomentum[3];

    //! Particle direction in local frame (normalized to 1)

    double _localDirection[3];

    //! Path length corresponding to Mokka hit

    double _mokkaPath;

    //! Energy deposit corresponding to Mokka hit

    double _mokkaDeposit;

    /*!
     *  Following variables contain information about sensor
     *  dimensions in the local coordinate frame
     */

    //! Sensor size in local reference frame
    /*! Sensor dimensions can be different in local reference frame
     *  as X can swap with Y due to rotation.
     */

    double _localSize[3];

    //! Sensor pitch in XY in local reference frame
    /*! Sensor pitch can be different in local reference frame
     *  X can swap with Y due to rotation
     */

    double _localPitch[3];

    //! Pixel Charge Map from TDS
    /*! Main structure used by Track Detailed Simulation (TDS)
     *  package. Energy deposit given by Mokka is converted to pixel
     *  charges and added to this map. Methods for adding
     *  fluctuations, noise, taking into account gain and threshold
     *  cut are also available.
     */

    TDS::TDSPixelsChargeMap  *_pixelChargeMap;

    std::map< int,  TDS::TDSPixelsChargeMap *>  _pixelChargeMapCollection;

    //! Vector of TDS pixels
    /*! Vector of TDS pixels is used as output contained for Track
     *   Detailed Simulation
     */

    std::vector<TDS::TDSPixel > _vectorOfPixels;
    std::vector<TDS::TDSPixel >::iterator _pixelIterator;

    //! Integration storage pointer for TDS
    TDS::TDSIntegrationStorage * _integrationStorage;


    //! Map for the TrackerData output collection
    std::map<int , lcio::TrackerDataImpl * > _trackerDataMap;

#if defined( USE_ALLPIX )
    //! digitiser interface to ALLPIX library
    vector<G4String> digitizerModulesNames;
    G4String GetNewName(G4String oldName, G4String toErase, G4String toAppend)
    {

	string oldName_s = oldName.data();
	string toErase_s = toErase.data();
	string toAppend_s = toAppend.data();
	string fixedString;

	size_t pos = oldName_s.find_last_of(toErase) - toErase_s.size();
	if(pos == string::npos) // couldn't find it
		return oldName;


	fixedString = oldName_s.substr(0,pos);
	fixedString.append(toAppend_s);

	return G4String(fixedString);
    }

    void SetupDetectors ();
#endif
  };

  //! A global instance of the processor
  EUTelMAPSdigi gEUTelMAPSdigi;

}
#endif
#endif
