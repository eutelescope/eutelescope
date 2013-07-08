// Author Aleksander Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 */

// built only if GEAR is used
#if defined( USE_GEAR )

// ROOT includes:
#include "TVector3.h"

// eutelescope includes ".h"
#include "EUTelMAPSdigi.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelSimpleSparsePixel.h"

// Include for Track Detailed Simulation (temporarily in Eutelescope)

// #include "TDSPixelsChargeMap.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogram3D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <iomanip>


// __endofheader__

using namespace std;
using namespace marlin;
using namespace gear;
using namespace TDS;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelMAPSdigi::_hitHistoLocalName          = "HitHistoLocal";
std::string EUTelMAPSdigi::_hitHistoTelescopeName      = "HitHistoTelescope";
std::string EUTelMAPSdigi::_pixelHistoName             = "PixelHisto";
std::string EUTelMAPSdigi::_chargeHistoName             = "ChargeHisto";
std::string EUTelMAPSdigi::_signalHistoName             = "SignalHisto";
std::string EUTelMAPSdigi::_multipHistoName             = "MultiplicityHisto";
std::string EUTelMAPSdigi::_chargeProfileName             = "ChargeProfile";
std::string EUTelMAPSdigi::_weightedProfileName             = "WeightedCharge";
#endif


EUTelMAPSdigi::EUTelMAPSdigi () 
: Processor("EUTelMAPSdigi"),
  _maxStepInCharge(0.0),
  _maxStepInLength(0.0),
  _integMaxNumberPixelsAlongL(0),
  _integMaxNumberPixelsAlongW(0),
  _chargeAttenuationLength(0.0),
  _chargeReflectedContribution(0.0),
  _gslFunctionCalls(0),
  _integPixelSegmentsAlongL(0),
  _integPixelSegmentsAlongW(0),
  _integPixelSegmentsAlongH(0),
  _useCommonIntegrationStorage(false),
  _ionizationEnergy(0.0),
  _DigiLayerIDs(),
  _depositedChargeScaling(),
  _applyPoissonSmearing(),
  _adcGain(),
  _adcGainVariation(),
  _adcNoise(),
  _adcOffset(),
  _zeroSuppressionThreshold(),
  _pixelReadoutType(),
  _adcRange(),
  _simhitCollectionName(""),
  _pixelCollectionName(""),
  _debugCount(0),
  _fillChargeProfiles(false),
  zsFrame(NULL),
  _iRun(0),
  _iEvt(0),
  _conversionIdMap(),
  _digiIdMap(),
  _siPlanesParameters(NULL),
  _siPlanesLayerLayout(NULL),
  _aidaHistoMap(),
  _histogramSwitch(false),
  _mokkaPath(0.0),
  _mokkaDeposit(0.0),
  _pixelChargeMap(NULL),
  _pixelChargeMapCollection(),
  _vectorOfPixels(),
  _pixelIterator(),
  _integrationStorage(NULL),
  _trackerDataMap()
   {

  // modify processor description
  _description =
    "EUTelMAPSdigi is used to calculate expected charge sharing between pixels"
    "for simulated Mokka hits";

  registerInputCollection(LCIO::SIMTRACKERHIT,"SimHitCollectionName",
                          "Simulated (Mokka) hit collection name",
                          _simhitCollectionName, string ( "telEUTelescopeCollection" ));

  registerOutputCollection(LCIO::TRACKERDATA,"PixelCollectionName",
                           "Collection name for simulated raw data",
                           _pixelCollectionName, string ("simEUTelescopeData"));

  // Track Detailes Simulation parameters

  registerProcessorParameter ("MaxStepInCharge",
                              "Maximum charge per simulation step [e]",
                              _maxStepInCharge,  static_cast < double > (100.));

  registerProcessorParameter ("MaxStepInLength",
                              "Maximum step length for single hit",
                              _maxStepInLength,  static_cast < double > (0.002));

  registerProcessorParameter ("IntegMaxNumberPixelsAlongL",
                              "Maximum diffusion range in pixels along sensor length ",
                              _integMaxNumberPixelsAlongL,  static_cast < int > (5));

  registerProcessorParameter ("IntegMaxNumberPixelsAlongW",
                              "Maximum diffusion range in pixels along sensor width ",
                              _integMaxNumberPixelsAlongW,  static_cast < int > (5));

  registerProcessorParameter ("ChargeAttenuationLength",
                              "Charge attenuation length in diffusion",
                              _chargeAttenuationLength,  static_cast < double > (0.055));

  registerProcessorParameter ("ChargeReflectedContribution",
                              "Charge reflection coefficient",
                              _chargeReflectedContribution,  static_cast < double > (1.0));

  registerProcessorParameter ("GSL function calls",
                              "Number of function calls in one GSL integration",
                              _gslFunctionCalls,  static_cast < int > (500));

  registerProcessorParameter ("IntegPixelSegmentsAlongL",
                              "Number of bins along sensor lenght for storing integration results",
                              _integPixelSegmentsAlongL,  static_cast < int > (16));


  registerProcessorParameter ("IntegPixelSegmentsAlongW",
                              "Number of bins along sensor width for storing integration results",
                              _integPixelSegmentsAlongW,  static_cast < int > (16));


  registerProcessorParameter ("IntegPixelSegmentsAlongH",
                              "Number of bins along sensor depth for storing integration results",
                              _integPixelSegmentsAlongH,  static_cast < int > (16));


  registerProcessorParameter ("Use common integration storage",
                              "Common integration storage can be used for all sensors if they have same geometry",
                              _useCommonIntegrationStorage,  static_cast < bool > (true));

  registerProcessorParameter ("IonizationEnergy",
                              "Ionization energy in silicon [eV]",
                              _ionizationEnergy,  static_cast < double > (3.6));

  // Digitization parameters

  std::vector< int > initLayerIDs;

  std::vector< float > initLayerScaling;
  initLayerScaling.push_back(1.0);

  std::vector< int > initLayerFlag;
  initLayerFlag.push_back(0);

  std::vector< float > initLayerGain;
  initLayerGain.push_back(1.);

  std::vector< float > initLayerGainVar;
  initLayerGainVar.push_back(0.);

  std::vector< float > initLayerNoise;
  initLayerNoise.push_back(0.);

  std::vector< float > initLayerOffset;
  initLayerOffset.push_back(0.);

  std::vector< int > initLayerRange;
  initLayerRange.push_back(0);

  std::vector< float > initLayerThreshold;
  initLayerThreshold.push_back(0.);

  std::vector< int > initLayerType;
  initLayerType.push_back(1);


  registerOptionalParameter ("DigiLayerIDs",
                             "Ids of layers which are considered in digitization",
                             _DigiLayerIDs, initLayerIDs);


  registerProcessorParameter ("DepositedChargeScaling",
                              "Scaling of the deposited charge (collection efficiency)",
                              _depositedChargeScaling , initLayerScaling );


  registerProcessorParameter ("ApplyPoissonSmearing",
                              "Flag for applying Poisson smearing of the collected charge",
                              _applyPoissonSmearing,  initLayerFlag );


  registerProcessorParameter ("AdcGain",
                              "ADC gain  in ADC counts per unit charge",
                              _adcGain,  initLayerGain );



  registerProcessorParameter ("AdcGainVariation",
                              "ADC gain variation ",
                              _adcGainVariation ,  initLayerGainVar );



  registerProcessorParameter ("AdcNoise",
                              "ADC noise in ADC counts",
                              _adcNoise ,  initLayerNoise );



  registerProcessorParameter ("AdcOffset",
                              "Constant pedestal offset in ADC counts",
                              _adcOffset ,  initLayerOffset );


  registerProcessorParameter ("AdcRange",
                              "Maximum value which can be stored as a pixel signal",
                              _adcRange ,  initLayerRange );


  registerProcessorParameter ("ZeroSuppressionThreshold",
                              "Threshold (in ADC counts) for removing empty pixels",
                              _zeroSuppressionThreshold ,  initLayerThreshold );


  registerProcessorParameter ("SensorReadoutType",
            "Flag for setting sensor readout type: 1 for digital, 2 for binary",
                              _pixelReadoutType,  initLayerType );



  // Other processor parameters


  registerProcessorParameter ("DebugEventCount",
                              "Print out every DebugEnevtCount event",
                              _debugCount,  static_cast < int > (10));


  registerProcessorParameter ("Fill charge distribution profiles",
                              "Charge profiles present average pixel charge, as a function of pixel number in cluster",
                              _fillChargeProfiles,  static_cast < bool > (false));


}


void EUTelMAPSdigi::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;


  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  streamlog_out ( ERROR4 ) <<  "Marlin was not built with GEAR support." << endl
                           <<  "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }
  
  if ( _DigiLayerIDs.empty() )
  {
    streamlog_out ( ERROR4 ) <<  "The lsit of sensors is not defined in the steering file. Pls go back to your template and specify _DigiLayerIDs." << endl;
    exit(-1);
  }
 
  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _histogramSwitch = true;

#endif

  // Check if digitization parameters are not missing (!)

  if(_depositedChargeScaling.empty() ||
     _applyPoissonSmearing.empty() ||
     _adcGain.empty() ||
     _adcGainVariation.empty() ||
     _adcNoise.empty() ||
     _adcRange.empty() ||
     _pixelReadoutType.empty() ||
     _zeroSuppressionThreshold.empty() ||
     _adcOffset.empty() )
    {
  streamlog_out ( ERROR4 ) <<  "Digitization parameters not given!" << endl
                    <<  "You have to specify all digitization parameters" << endl;
  exit(-1);
    }


  // Check if correct number of digitization parameters is given

  /*  unsigned int nTelPlanes = _siPlanesParameters->getSiPlanesNumber();

  if(_DigiLayerIDs.size()>0 && _DigiLayerIDs.size()!=nTelPlanes)
    {
  streamlog_out ( ERROR4 ) <<  "Wrong number of Layer IDs for digitization" << endl
                  <<  "Size of DigiLayerIDs have to match GEAR description" << endl;
  exit(-1);
    }

  if(_DigiLayerIDs.size()==0)nTelPlanes=1;

  if(_depositedChargeScaling.size() != nTelPlanes ||
     _applyPoissonSmearing.size() != nTelPlanes ||
     _adcGain.size() != nTelPlanes ||
     _adcGainVariation.size() != nTelPlanes ||
     _adcNoise.size() != nTelPlanes ||
     _adcRange.size() != nTelPlanes ||
     _pixelReadoutType.size() != nTelPlanes ||
     _zeroSuppressionThreshold.size() != nTelPlanes ||
     _adcOffset.size() != nTelPlanes )
    {
 streamlog_out ( ERROR4 ) <<  "Wrong size of digitization parameter vector" << endl;
  exit(-1);
    }
*/
  // prepare digitization parameter index map

  if( !_DigiLayerIDs.empty() )
    for(unsigned int id=0; id<_DigiLayerIDs.size(); id++)
      _digiIdMap.insert( make_pair(_DigiLayerIDs.at(id), id ) );


  // Book common integration storage for all sensors, if this was
  // requested

  if (_useCommonIntegrationStorage)
    {
      _integrationStorage = new TDSIntegrationStorage(_integPixelSegmentsAlongL,
                                                      _integPixelSegmentsAlongW, _integPixelSegmentsAlongH);
      streamlog_out( MESSAGE4 )  << " Common TDS integration storage initialized  " << endl;
    }

#if defined( USE_ALLPIX)
  digitizerModulesNames.clear();
#endif


}

void EUTelMAPSdigi::processRunHeader (LCRunHeader * rdr) {

// This processor should only be used for simulated data
// Check if the input was generated with Mokka


  streamlog_out( MESSAGE4 )  << "  Run : " << rdr->getRunNumber()
                             << " - "      << rdr->getDetectorName()
                             << ":  "      << rdr->getDescription()  << endl ;

  string simulator = rdr->parameters().getStringVal("SimulatorName");

  if(simulator != "")
    {
      streamlog_out( MESSAGE4 )  << simulator << " input file recognized, version  "
                                 << rdr->parameters().getStringVal("SimulatorVersion")  << endl;
    }
  else
    {
      streamlog_out ( ERROR4 ) << "Error during consistency check: " << endl
                               << "EUTelMAPSdigi processor can only run on simulated data"  << endl
                               << "but SimulatorName not set !?" << endl;
      exit(-1);
    }


  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() );

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file.  But it is not always set for simulation output


  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() )
    streamlog_out ( WARNING0 ) <<  "Warning during the geometry consistency check: " << endl
                               << "The run header says the GeoID is " << header->getGeoID() << endl
                               << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID() << endl;



  // now book histograms plz...
  if ( isFirstEvent() )  bookHistos();

  // increment the run counter
  ++_iRun;
}

#if defined( USE_ALLPIX )
void EUTelMAPSdigi::SetupDetectors () 
{
        int itr        = 0;
	int detectorId = 0;
	digitizerModulesNames.clear();

	// Digit manager
	G4DigiManager * fDM = G4DigiManager::GetDMpointer();

		G4String hcName = "";

		// Now build the digit Collection name
		G4String digitColectionName = _simhitCollectionName;
         	G4String digitizerName = "Mimos26";

		// Creating an instance of the actual digitizer, and keep pointer through the interface
		AllPixDigitizerInterface * dmPtr;

		G4String digitSuffix = "_";
		digitSuffix += digitizerName;
		digitSuffix += "Digitizer";
		G4String digitizerModName =  GetNewName(hcName, "HitsCollection", digitSuffix);


		digitizerModulesNames.push_back(digitizerModName);

                    
	        // __beginofdigitlist__
		if (digitizerName == "FEI3Standard") {
			AllPixFEI3StandardDigitizer * dp = new AllPixFEI3StandardDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		} else if (digitizerName == "Medipix2") {
			AllPixMedipix2Digitizer * dp = new AllPixMedipix2Digitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		} else if (digitizerName == "Mimosa26") {
			AllPixMimosa26Digitizer * dp = new AllPixMimosa26Digitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		} else if (digitizerName == "Timepix") {
			AllPixTimepixDigitizer * dp = new AllPixTimepixDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		}else if (digitizerName == "MCTruth") {
			AllPixMCTruthDigitizer * dp = new AllPixMCTruthDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
			dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
			cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
		}// Included by newdigitizer.sh script --> LETCalculator
		else if (digitizerName == "LETCalculator") {
					AllPixLETCalculatorDigitizer * dp = new AllPixLETCalculatorDigitizer(digitizerModulesNames[itr] , hcName, digitColectionName);
					dmPtr = static_cast<AllPixDigitizerInterface *> (dp);
					cout << "    Setting up a " << digitizerName << " digitizer for det : " << detectorId << endl;
				}
 
        	// __endofdigitlist__
		else {
			G4cout << "    can't find digitizer with name : " << digitizerName << G4endl;
			exit(1);
		}

		// push back the digitizer
//		m_digiPtrs.push_back( dmPtr );
		fDM->AddNewModule( dmPtr );
}
#endif


void EUTelMAPSdigi::processEvent (LCEvent * event) {

  int  debug = ( _debugCount>0 && _iEvt%_debugCount == 0);

  debug = 0;

  if (_iEvt % _debugCount == 0  || debug)
    streamlog_out( MESSAGE4 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
                              << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;

  ++_iEvt;

  event->parameters().setValue( "EventType",2);

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

// Event type not set in the Mokka output
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  }


  LCCollectionVec * simhitCollection   = 0;

  try {

    simhitCollection = static_cast<LCCollectionVec*> (event->getCollection( _simhitCollectionName ));

  } catch (DataNotAvailableException& e  ) {

    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }


  // prepare also an output collection
  auto_ptr< lcio::LCCollectionVec > zsDataCollection( new LCCollectionVec ( LCIO::TRACKERDATA ) ) ;

  // prepare also the corresponding encoder
  CellIDEncoder< TrackerDataImpl > zsDataEncoder ( EUTELESCOPE::ZSDATADEFAULTENCODING, zsDataCollection.get() ) ;

    int detectorID    = -99; // it's a non sense
    int oldDetectorID = -100;
    int mapDetectorID = -100;

    int    layerIndex = -99;
    double xZero = 0., yZero = 0., zZero = 0. ;
    double xSize = 0., ySize = 0.;
    double zThickness = 0.;
    double xPitch = 0., yPitch = 0.;
    double xPointing[2] = { 1., 0. }, yPointing[2] = { 1., 0. };

    double gRotation[3] = { 0., 0., 0.}; // not rotated



    for ( int iHit = 0; iHit < simhitCollection->getNumberOfElements(); iHit++ ) 
    {
      if(debug>1) streamlog_out ( MESSAGE5 ) << "iHit = " << iHit << " :of: " << simhitCollection->getNumberOfElements() << endl;
      
      SimTrackerHitImpl * simhit   = static_cast<SimTrackerHitImpl*> ( simhitCollection->getElementAt(iHit) );
     
      // there could be several hits belonging to the same
      // detector. So update the geometry information only if this new
      // cluster belongs to a different detector.
      detectorID = simhit->getCellID0();
      if( detectorID > 99 ) continue;

      if( !_DigiLayerIDs.empty() && _DigiLayerIDs[0] <  10 ) detectorID -=  1; 

      //
 
      bool detector_type_choice = false;
      for(unsigned int iPlane = 0; iPlane < _DigiLayerIDs.size(); iPlane++)
      {
        if(detectorID == _DigiLayerIDs[iPlane] ) detector_type_choice = true;
      }

      if( detector_type_choice == false ) 
      {
        // clear stored pixel map
        //if( _pixelChargeMap != 0) _pixelChargeMap->clear();
        //if( _vectorOfPixels.size() != 0) _vectorOfPixels.clear();

        //break;
      }
 
      if(debug>1) streamlog_out ( MESSAGE5 ) << "Sensitive layer ID = " << detectorID << " found in SimTracker collection" << endl;
 
      if ( detectorID != oldDetectorID )
      {

        oldDetectorID = detectorID;

        if ( _conversionIdMap.size() != static_cast< unsigned >( _siPlanesParameters->getSiPlanesNumber()) ) 
        {
          // first of all try to see if this detectorID already belong to
          if ( _conversionIdMap.find( detectorID ) == _conversionIdMap.end() ) 
          {
            // this means that this detector ID was not already inserted,
            // so this is the right place to do that
            for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) 
            {
              if ( _siPlanesLayerLayout->getSensitiveID(iLayer) == detectorID ) 
              {
                _conversionIdMap.insert( make_pair( detectorID, iLayer ) );
                if(debug>1) streamlog_out ( MESSAGE5 ) << "Sensitive layer ID = " << detectorID << " found in GEAR" << endl;
                break;
              }
            }
          }
        }

        // perfect! The full geometry description is now coming from the
        // GEAR interface. Let's keep the finger xed!
        layerIndex   = _conversionIdMap[detectorID];

        xZero        = _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
        yZero        = _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
        zZero        = _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
        zThickness   = _siPlanesLayerLayout->getSensitiveThickness(layerIndex); // mm
        xPitch       = _siPlanesLayerLayout->getSensitivePitchX(layerIndex);    // mm
        yPitch       = _siPlanesLayerLayout->getSensitivePitchY(layerIndex);    // mm
        xSize        = _siPlanesLayerLayout->getSensitiveSizeX(layerIndex);     // mm
        ySize        = _siPlanesLayerLayout->getSensitiveSizeY(layerIndex);     // mm
        xPointing[0] = _siPlanesLayerLayout->getSensitiveRotation1(layerIndex); // was -1 ;
        xPointing[1] = _siPlanesLayerLayout->getSensitiveRotation2(layerIndex); // was  0 ;
        yPointing[0] = _siPlanesLayerLayout->getSensitiveRotation3(layerIndex); // was  0 ;
        yPointing[1] = _siPlanesLayerLayout->getSensitiveRotation4(layerIndex); // was -1 ;
        int npixelX = _siPlanesLayerLayout->getSensitiveNpixelX( layerIndex );
        int npixelY = _siPlanesLayerLayout->getSensitiveNpixelY( layerIndex );
        xSize = xPitch*npixelX;
        ySize = yPitch*npixelY;


        try
          {
          gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(layerIndex); // Euler alpha ;
          gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(layerIndex); // Euler alpha ;
          gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(layerIndex); // Euler alpha ;

        // input angles are in DEGREEs !!!
        // translate into radians
          gRotation[0] =  gRotation[0]*3.1415926/180.; // 
          gRotation[1] =  gRotation[1]*3.1415926/180.; //
          gRotation[2] =  gRotation[2]*3.1415926/180.; //
          }
          catch(...)
          {
              streamlog_out( DEBUG5 ) << "No sensor rotation is given in the GEAR steering file, assume NONE." << endl;
          }


        if (  ( xPointing[0] == xPointing[1] ) && ( xPointing[0] == 0 ) ) 
        {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;
        }

        if (  ( yPointing[0] == yPointing[1] ) && ( yPointing[0] == 0 ) ) 
        {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;

        }

      }  // End of if ( detectorID != oldDetectorID )


      // get the position and momentum of particle generating energy deposit

      const double* hitpos = simhit->getPosition();
      const float*  hitmom = simhit->getMomentum();

      if (debug>1)
        streamlog_out( MESSAGE5 ) << "SimHit in global frame  at X = " <<  hitpos[0]
                                <<  " Y = " <<  hitpos[1]
                                <<  " Z = " <<  hitpos[2]
                                <<  "  detector ID = " << detectorID << endl;

      // 2D histograms of simulated hits in global (telescope) reference frame

      fillHist2D(_hitHistoTelescopeName,layerIndex,hitpos[0], hitpos[1] );



      // transform position and momentum to local (sensor) reference frame
      // sensor rotation is described by two vectors (from implementation in HitMaker):
      // (xPointing[0],yPointing[0],0) gives direction of sensor X axis in telescope reference frame
      // (xPointing[1],yPointing[1],0) gives direction of sensor Y axis in telescope reference frame
      // direction of Z axis is unchanged, only pointing can flip


      // first the translation to the frame with origin at the sensor center

      double senspos[3];

      senspos[0] = hitpos[0] -  xZero;
      senspos[1] = hitpos[1] -  yZero;
      senspos[2] = hitpos[2] -  zZero;


      // NEW: apply sensor rotation (if defined in GEAR)

      InvEulerRotation(senspos, gRotation);

      double sensmom[3];

      sensmom[0] = static_cast< double >( hitmom[0]);
      sensmom[1] = static_cast< double >( hitmom[1]);
      sensmom[2] = static_cast< double >( hitmom[2]);

      InvEulerRotation(sensmom, gRotation);

      // now move to the sensor corner. Sign determines which way to move

      // X axis
      double sign = 1;

      if      ( xPointing[0] < -0.7 )       sign = -1 ;
      else if ( xPointing[0] >  0.7 )       sign =  1 ;
      else {
        if       ( xPointing[1] < -0.7 )    sign = -1 ;
        else if  ( xPointing[1] >  0.7 )    sign =  1 ;
      }

      senspos[0] +=  sign * xSize/2;

      // Y axis

      if      ( yPointing[0] < -0.7 )       sign = -1 ;
      else if ( yPointing[0] >  0.7 )       sign =  1 ;
      else {
        if       ( yPointing[1] < -0.7 )    sign = -1 ;
        else if  ( yPointing[1] >  0.7 )    sign =  1 ;
      }

      senspos[1] +=  sign * ySize/2;

      // Z axis direction from determinant of the XY rotation matrix

      sign = xPointing[0] * yPointing[1] - xPointing[1] * yPointing[0] ;

      senspos[2] += sign * zThickness/2;

      // Last step: rotation taking into account sensor orientation.
      // Reverse rotation has to be applied!
      // "sign" is the determinant of the rotation matrix


      _localPosition[0] =  yPointing[1]/sign * senspos[0] - xPointing[1]/sign * senspos[1] ;
      _localPosition[1] = -yPointing[0]/sign * senspos[0] + xPointing[0]/sign * senspos[1] ;
      _localPosition[2] = sign * senspos[2] ;

      if (debug)
        streamlog_out( DEBUG4 ) << "Position in local frame at X = " << _localPosition[0]
                                <<  " Y = " << _localPosition[1]
                                <<  " Z = " << _localPosition[2] << endl;


      // Fill local position histogram

      fillHist2D(_hitHistoLocalName,layerIndex,_localPosition[0], _localPosition[1]);


      // momentum transformation (only rotation)


      _localMomentum[0] =  yPointing[1]/sign * sensmom[0] - xPointing[1]/sign * sensmom[1] ;
      _localMomentum[1] = -yPointing[0]/sign * sensmom[0] + xPointing[0]/sign * sensmom[1] ;
      _localMomentum[2] = sign * sensmom[2] ;

      // Normalized direction vector

      double ptot;
      ptot=sqrt(_localMomentum[0]*_localMomentum[0]+_localMomentum[1]*_localMomentum[1]+_localMomentum[2]*_localMomentum[2]);


      for(int idir=0; idir<3; idir++)
        _localDirection[idir]=_localMomentum[idir]/ptot;


      // Sensor size in the local (sensor) frame (X can swap with Y due to rotation!)


      _localSize[0] =  yPointing[1]/sign * xSize - xPointing[1]/sign * ySize ;
      if(_localSize[0] < 0)_localSize[0] = -_localSize[0];

      _localSize[1] = -yPointing[0]/sign * xSize + xPointing[0]/sign * ySize ;
      if(_localSize[1] < 0)_localSize[1] = -_localSize[1];



      _localSize[2] = zThickness;

      if(debug>1) streamlog_out( MESSAGE5 ) <<
                               " _localSize[0] = " << _localSize[0] << 
                               " _localSize[1] = " << _localSize[1] << 
                               " _localSize[2] = " << _localSize[2] << endl;


      // From EUTelHitMaker code it looks like sensor pitches given in
      // GEAR are already in the local frame of reference:

      _localPitch[0] = xPitch;
      _localPitch[1] = yPitch;
      _localPitch[2] = zThickness;



      // Remaining information on the simulated hit:

//      _mokkaDeposit =simhit->getdEdx(); // OBSOLETE member function
      _mokkaDeposit =simhit->getEDep(); 
      _mokkaPath=simhit->getPathLength();

      ///////////////////////////////////////////////////////////////////////////
      //
      // Now we are in the local coordinate frame, and we should have
      // all information needed for digitization
      //
      //  Hit position in the local frame of reference:  _localPosition[3]
      //  Particle momentum [MeV] in local frame:  _localMomentum[3]
      //  Particle direction in local frame:  _localDirection[3]
      //  Total energy deposit:  _mokkaDeposit
      //  Total particle path lenght:  _mokkaPath
      //
      //  Sensor boundaries in the local frame are:
      //    (0,0,0) to _localSize[3]
      //
      //  Sensor pitch is given by _localPitch[2] (X and Y only)
      //
      //////////////////////////////////////////////////////////////////////////



      // We have to initialize TDS pixel charge map if not done yet

      if ( detectorID != mapDetectorID && detectorID < 100 ) 
      {
        mapDetectorID = detectorID;

        if ( _pixelChargeMapCollection.size() != static_cast< unsigned >( _DigiLayerIDs.size()) )
        {
          // first of all try to see if this detectorID already belong to
          if ( _pixelChargeMapCollection.find( detectorID ) == _pixelChargeMapCollection.end()  ) 
          {
            // this means that the map for was not defined yet,
            // so this is the right place to do that
            std::cout << " Initialising a TDS Charge Map " << std::endl;

            // current TDS requires that negative sensor thickness is given!
            _pixelChargeMap = new TDSPixelsChargeMap(_localSize[0],_localSize[1],-_localSize[2]);

            std::cout << "iHit: " << iHit << " " << _pixelChargeMap << endl;

            _pixelChargeMapCollection.insert( make_pair( detectorID,_pixelChargeMap) );

            // set detector type (Default MAPS)
            if( detectorID >=20 ) _pixelChargeMap->setDetectorType("FEI4");
 
            // Set additional geometry parameters
            _pixelChargeMap->setPixelLength(_localPitch[0]);
            _pixelChargeMap->setPixelWidth(_localPitch[1]);

            // Set parameters of the charge diffusion model (input parameters)

            _pixelChargeMap->setLambda(_chargeAttenuationLength);
            _pixelChargeMap->setReflectedContribution(_chargeReflectedContribution);

            // Initialize integration parameters (all taken from input
            // processor parameters)

            _pixelChargeMap->initializeIntegration(_maxStepInLength,_maxStepInCharge,
                                                   _integMaxNumberPixelsAlongL, _integMaxNumberPixelsAlongW, _gslFunctionCalls);


            streamlog_out( MESSAGE4 )  << " Pixel map initialized  for detectorID "
                                       << detectorID << endl;


            // Initialize integration storage

            if (!_useCommonIntegrationStorage)
              {
                _integrationStorage = new TDSIntegrationStorage(_integPixelSegmentsAlongL,
                                                                _integPixelSegmentsAlongW, _integPixelSegmentsAlongH);

                streamlog_out( MESSAGE4 )  << " TDS integration storage initialized  for detectorID "
                                           << detectorID << endl;

              }

            _pixelChargeMap->setPointerToIntegrationStorage(_integrationStorage);
          }
        }

        // Get pointer to the TDS pixel charge map
        _pixelChargeMap = _pixelChargeMapCollection[detectorID];

      }


      // First step is to calculate track segment starting point and end point
      // One has to check if not going outside the sensor when moving by +/- path/2

      double pathFrac = CheckPathLimits();

      if (debug && pathFrac < 0.9999)
        streamlog_out( DEBUG4 ) << "Track path inside sensor reduced by factor " << pathFrac << endl;

      if (pathFrac <= 0.)
        streamlog_out( WARNING4 ) << "Track path outside sensor !? " << endl
                                  << "Position in local frame at X = " << _localPosition[0]
                                  <<  " Y = " << _localPosition[1]
                                  <<  " Z = " << _localPosition[2] << endl;

      //
      //  Here comes the core of the algorithm
      //
      if(debug>1) 
           streamlog_out( DEBUG4 ) << "going to convert mokka deposit [GeV] to charge units of elementary charge " << endl;
      if(debug>1) 
           streamlog_out( DEBUG4 ) << "_mokkaPath " << _mokkaPath << endl;
      if(debug>1) 
           streamlog_out( DEBUG4 ) << "_mokkaDeposit"<< _mokkaDeposit << endl;



      if(_mokkaPath > 0 && _mokkaDeposit > 0.)
        {
        // Convert Mokka deposit [GeV] to charge in units of elementary
        // charge

        _mokkaDeposit *= 1000000000./_ionizationEnergy;

        // Set the Geant/Mokka step to be considered
        // TDS requires position in depth to be negative ! Shift it by
        // thickness to keep direction unchanged


        if (debug>1 ||  ( _localPosition[2]-_localSize[2] > 0)  || ( _localPosition[2] > 0.) )
          streamlog_out( MESSAGE5 ) << " sensor = " << detectorID << 
                                     " _localPosition[0] = " << _localPosition[0] <<
                                     " _localPosition[1] = " << _localPosition[1] <<
                                     " _localPosition[2] = " << _localPosition[2] << endl << 
                                     " _localSize[0]     = " << _localSize[0]     <<
                                     " _localSize[1]     = " << _localSize[1]     <<
                                     " _localSize[2]     = " << _localSize[2]     << endl <<
                                     " _ localDirection[0]= " << _localDirection[0]<<
                                     " _localDirection[1]= " << _localDirection[1]<<
                                     " _localDirection[2]= " << _localDirection[2]<< endl <<
                                     " _mokkaPath        = " << _mokkaPath  <<
                                     " _mokkaDeposit     = " << _mokkaDeposit << endl;

        if( _localPosition[2] > 0.) continue; 
        TDSStep step(_localPosition[0], _localPosition[1], _localPosition[2]-_localSize[2],
                     _localDirection[0],_localDirection[1],_localDirection[2],
                     _mokkaPath, _mokkaDeposit);

        // distribute charge among pixels (here all the work is done!)

        _pixelChargeMap->update(step);

        if (debug>0)
          streamlog_out( DEBUG4 ) << "Adding Mokka deposit of " << _mokkaDeposit <<  endl;
        }
    }
    // end of loop over SimTrackerHit collection


//====================================================================
//
// Process collected charges and write all pixels to output collection
//
//====================================================================

    // prepare the output TrackerData to host the SparsePixel

for(unsigned int idet = 0; idet < _DigiLayerIDs.size(); idet++)
{
    int idetectorID = _DigiLayerIDs[idet];

    // check if the we already have a tracker data for this
    // detector id
    if( idetectorID< 10)
    {
      for ( int iDetector = 0 ; iDetector < _siPlanesParameters->getSiPlanesNumber(); iDetector++ ) 
      {
        if( _siPlanesLayerLayout->getSensitiveID(iDetector) > 10 ) continue;
        TrackerDataImpl * tmp_zsFrame      = new TrackerDataImpl;
        zsDataEncoder[ "sensorID" ]        = _siPlanesLayerLayout->getSensitiveID(iDetector);
        zsDataEncoder[ "sparsePixelType" ] = kEUTelSimpleSparsePixel;
        zsDataEncoder.setCellID( tmp_zsFrame );
        _trackerDataMap[ _siPlanesLayerLayout->getSensitiveID(iDetector) ] = tmp_zsFrame;
        if(debug>0)streamlog_out( MESSAGE5 ) << " map : ("<< iDetector <<") "  << _siPlanesLayerLayout->getSensitiveID(iDetector) << endl;
      }
      break;
    }
    else
        if( idetectorID>=20 && idetectorID < 100 )
    {
      for ( int iDetector = 0 ; iDetector < _siPlanesParameters->getSiPlanesNumber(); iDetector++ ) 
      {
        if( _siPlanesLayerLayout->getSensitiveID(iDetector) <  20 ) continue;
        TrackerDataImpl * tmp_zsFrame      = new TrackerDataImpl;
        zsDataEncoder[ "sensorID" ]        = _siPlanesLayerLayout->getSensitiveID(iDetector);
 	zsDataEncoder["sparsePixelType"]   = kEUTelAPIXSparsePixel;
        zsDataEncoder.setCellID( tmp_zsFrame );
        _trackerDataMap[ _siPlanesLayerLayout->getSensitiveID(iDetector) ] = tmp_zsFrame;
        if(debug>0)streamlog_out( MESSAGE5 ) << " map : ("<< iDetector <<") " << _siPlanesLayerLayout->getSensitiveID(iDetector) << endl;
      }
      break; 
    }
    else
    {
      streamlog_out( ERROR5 ) << "Detector ID ["<< idetectorID << "] is not specific to the telescope planes (0-9, 20-29)... exit" << endl;
      exit(-1); 
    }
}

    // Loop over defined detectors

    std::map< int,  TDSPixelsChargeMap *>::iterator mapIterator;

    for(mapIterator = _pixelChargeMapCollection.begin();
        mapIterator != _pixelChargeMapCollection.end(); ++mapIterator)
      {
        detectorID = mapIterator->first;


        layerIndex   = _conversionIdMap[detectorID];

        int digiIndex = 0;

        if( !_DigiLayerIDs.empty() ) digiIndex = _digiIdMap[detectorID];

        _pixelChargeMap = mapIterator->second;


        // Charge processing
        if (debug)
          streamlog_out( MESSAGE5 ) << " _pixelChargeMap " << _pixelChargeMap << " det " << detectorID << " lay " << layerIndex << endl;
       
        double totalCharge=_pixelChargeMap->getTotalCharge();

        if (debug)
          streamlog_out( DEBUG4 ) << "Total charge collected in detector " << detectorID <<
            " : " <<   totalCharge <<  endl;

        // Collected charge histogram

        fillHist1D(_chargeHistoName, layerIndex, totalCharge);

        // Scaling of charge deposited at each pixel

        if(_depositedChargeScaling[digiIndex]!=1.)
          _pixelChargeMap->scaleCharge(static_cast< double >(_depositedChargeScaling[digiIndex]));


        // Poisson smearing (if requested)

        if(_applyPoissonSmearing[digiIndex])
          {
            _pixelChargeMap->applyPoissonFluctuations(false);

            if (debug)
              streamlog_out( DEBUG4 ) << "Charge after Poisson fluctuations: "
                                      << _pixelChargeMap->getTotalCharge() << endl;

          }

        // ADC gain, noise, pedestal

        _pixelChargeMap->applyGain(static_cast< double >(_adcGain[digiIndex]), static_cast< double >(_adcGainVariation[digiIndex]), static_cast< double >(_adcNoise[digiIndex]), static_cast< double >(_adcOffset[digiIndex]));

        // TDS allow for negative charges and negative thresholds.
        // However we assume here that threshold has to be
        // positive. If zero or negative threshold value is set, not
        // threshold correction is applied.

        if(_zeroSuppressionThreshold[digiIndex]>0)
          _pixelChargeMap->applyThresholdCut(static_cast< double >(_zeroSuppressionThreshold[digiIndex]));

        totalCharge=_pixelChargeMap->getTotalCharge();

        if (debug)
          streamlog_out( DEBUG4 ) << "Total signal after gain and ZS: "
                                  << totalCharge << endl;


        // Signal histogram (charge after procession and ADC conversions)

        fillHist1D(_signalHistoName, layerIndex, totalCharge);


        // Get vector of fired pixels from map

        _vectorOfPixels = _pixelChargeMap->getVectorOfPixels();


        // Pixel multiplicity histogram

        fillHist1D(_multipHistoName, layerIndex, static_cast< double >(_vectorOfPixels.size()));


        // Data copying to output stream starts here
        // =========================================

        // check if the we already have a tracker data for this
        // detector id
        zsFrame = NULL;
        if ( _trackerDataMap.find( detectorID ) == _trackerDataMap.end() ) {
          // this is the first time we have to deal with such a
          // detector. so first of all let's create the
          // corresponding TrackerData
          //
          // clear stored pixel map
          _pixelChargeMap->clear();
          _vectorOfPixels.clear();

          // unable to find the proper zsFrame
          if(debug) streamlog_out ( ERROR4 ) << "Unable to find the TrackerData corresponding to detectorID " << detectorID << endl;
          continue; //throw SkipEventException(this);
        } else {
          zsFrame = _trackerDataMap[ detectorID ];
        }


        if( detectorID< 10 )
        {
           //streamlog_out ( MESSAGE5 )<< "magic. Mimosa type"<< endl;
           packMimosa26( digiIndex );
        }
        else
            if( detectorID >= 10 && detectorID < 20 )
        {
           streamlog_out ( DEBUG5 ) << "magic. FEI3   type - Currently no function packFEI3"<< endl;
        }
        else
            if( detectorID >= 20 && detectorID < 30 )
        {
           //streamlog_out ( MESSAGE5 )<< "magic. FEI4   type"<< endl;
           packFEI4( digiIndex );
        }
        else
        {

        }
///////////////////////////aaaaaaaaaaaaaaa


        // List all fired pixels

        int nPixel = _vectorOfPixels.size() ;

        if (debug>1)
        {
            streamlog_out( MESSAGE4 ) <<  "Detector ID = " << detectorID
                                      << ", " << nPixel << " pixels fired, "
                                      << "total charge deposited: " << totalCharge << endl;

            int ipixel = 0;
            for(_pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); ++_pixelIterator)
            {
              streamlog_out( MESSAGE5 ) <<  " Pixel [" <<  ipixel << "] at  (" << _pixelIterator->getIndexAlongL()
                                      <<  "," <<  _pixelIterator->getIndexAlongW()
                                      <<  ") with Q = " << _pixelIterator->getCharge() << endl;
              ipixel++;
            }


// Store pixel charges in the 2D histogram (charge map)
// Only for events with debug output!

            for(_pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); ++_pixelIterator)
            {
              double xpixel = _pixelIterator->getIndexAlongL() + 0.5;
              double ypixel = _pixelIterator->getIndexAlongW() + 0.5;
              fillHist2D(_pixelHistoName, layerIndex, xpixel,ypixel,_pixelIterator->getCharge());
            }

        }// end of if (debug)


        // charge profile filling (for algorithm tuning)
        // consider only pixels around the highest pixel
        // to make this selection we have to sort pixels again
        // skip cases when there is no or one pixel only (unlikely)

        if(_fillChargeProfiles && nPixel > 1 )
          {
            // Maps for pixel sorting in abs(charge) and in
            // abs(charge)/distance_to_seed^2

            std::multimap<double,double, std::greater<double> > pixelAbsMap;
            std::multimap<double,double, std::greater<double> > pixelWeightedMap;


            std::multimap<double,double, std::greater<double> >::iterator pixelAbsIterator;
            std::multimap<double,double, std::greater<double> >::iterator pixelWeightedIterator;


            int dLmax = _integMaxNumberPixelsAlongL/2;
            int dWmax = _integMaxNumberPixelsAlongW/2;

            // In the new TDS version
            // TDSPixelsChargeMap::getVectorOfPixels() returns  sorted
            // vector of pixels, starting from the highest deposit -
            // it can be used as a seed pixel (assuming we consider
            // one cluster only)

            int seedL =   _vectorOfPixels[0].getIndexAlongL();
            int seedW =   _vectorOfPixels[0].getIndexAlongW();

            double clusterCharge=0.;
            int intClusterCharge=0;

            // store cluster pixels in maps for sorting (maps are
            // sorted by the key value)

            for(_pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); )
              {
                int pixelL = _pixelIterator->getIndexAlongL();
                int pixelW = _pixelIterator->getIndexAlongW();

                if(pixelL > seedL+dLmax || pixelL < seedL-dLmax ||
                   pixelW > seedW+dWmax || pixelW < seedW-dWmax )
                  ++_pixelIterator;
                else
                  {
                   double dist = (pixelL-seedL)*(pixelL-seedL)+(pixelW-seedW)*(pixelW-seedW);
                   double charge = _pixelIterator->getCharge();

                   // ADC digitization simulation

                   int icharge = static_cast< int >(charge);
                   if(charge<0.) icharge--;

                   // double absCharge = fabs(charge);
                   double absCharge = abs(icharge);

                   double weight = absCharge;
                   if(dist>0.5)weight/=sqrt(dist);

                   pixelAbsMap.insert( make_pair( absCharge, static_cast< double >( icharge)));
                   pixelWeightedMap.insert( make_pair( weight, static_cast< double >( icharge)));

                   clusterCharge+=charge;
                   intClusterCharge+=icharge;
                   ++_pixelIterator;
                  }
              }


            // fill charge profiles (maps are sorted by definition!)

           int iPixel=0;
           for(pixelAbsIterator=pixelAbsMap.begin();
               pixelAbsIterator != pixelAbsMap.end(); ++pixelAbsIterator)
              {
              fillProf1D(_chargeProfileName,layerIndex,iPixel+0.5, pixelAbsIterator->second);
              iPixel++;
              }


           iPixel=0;
           for(pixelWeightedIterator=pixelWeightedMap.begin();
               pixelWeightedIterator != pixelWeightedMap.end(); ++pixelWeightedIterator)
              {
              fillProf1D(_weightedProfileName,layerIndex,iPixel+0.5, pixelWeightedIterator->second);
              iPixel++;
              }


            // Test output at the end of processing
            //
            if(debug>1)
                streamlog_out( DEBUG4 ) <<  "Seed at L = " << seedL << " W = " << seedW
                               << ", " << pixelAbsMap.size() << " pixels in cluster " << endl
                               << "Total cluster charge: " << clusterCharge
                               << ", integer cluster charge: " << intClusterCharge << endl;

          }
        // end of if(_fillChargeProfiles...)

        // clear stored pixel map
        _pixelChargeMap->clear();
        _vectorOfPixels.clear();

      } // end of loop over defined pixel maps

        // before clearing the content of the trackerDataMap, I need
        // to push back to the output collection all the trackerdata I
        // have
    map< int, TrackerDataImpl * >::iterator trackerDataIter = _trackerDataMap.begin();
    while ( trackerDataIter != _trackerDataMap.end() ) 
    {

      zsDataCollection->push_back( trackerDataIter->second ) ;

      ++trackerDataIter;
    }


    _trackerDataMap.clear();

    // now, dulcis in fundo, add this collection to the event!
    event->addCollection( zsDataCollection.release(), _pixelCollectionName );



    if ( isFirstEvent() ) _isFirstEvent = false;

}

void EUTelMAPSdigi::packMimosa26( int digiIndex ) 
{
  int debug = 0;

//  streamlog_out ( MESSAGE4 )  << "packing Mimosa26" << endl;

  if( zsFrame == 0 )
  {
     streamlog_out( ERROR5 ) << " zsFrame is absent. emergency exit" << endl;
     exit(-1);
  } 
        // now I have to prepare a temporary storage for the
        // sparse pixel
        auto_ptr< EUTelSparseDataImpl< EUTelSimpleSparsePixel > > sparseFrame( new EUTelSparseDataImpl< EUTelSimpleSparsePixel > ( zsFrame ) );
        auto_ptr< EUTelSimpleSparsePixel > sparsePixel( new EUTelSimpleSparsePixel );
 

        // and also a temporary storage for the sparse pixel
					
        //sparseFrame->addSparsePixel( thisHit );
	//tmphits.push_back( thisHit );
        for ( _pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); ++_pixelIterator ) 
        {

          // move the pixel information from the pixelIterator to
          // the sparsePixel
          sparsePixel->setXCoord( _pixelIterator->getIndexAlongL() );
          sparsePixel->setYCoord( _pixelIterator->getIndexAlongW() );

	  // Take readout type and ADC range into account

          int pixelValue = 1;        // For binary readout all pixels
				     // above threshold are set to 1

          if( _pixelReadoutType[digiIndex] == 1 )  // digital readout
	    {
            double pixelCharge = _pixelIterator->getCharge();

            pixelValue = static_cast< int >(pixelCharge);

            if(pixelCharge<0.) pixelValue--; 

 
           if(_adcRange[digiIndex] > 0)
	      {
              if(pixelValue > _adcRange[digiIndex])pixelValue=_adcRange[digiIndex];
              if(pixelValue < 0 )pixelValue=0;
	      }
	    }

          if (debug>1)
                     streamlog_out( MESSAGE4 ) <<  "pixel charge " << pixelValue << " X " <<  _pixelIterator->getIndexAlongL() << " " <<  _pixelIterator->getIndexAlongW() << endl;
          sparsePixel->setSignal( pixelValue );

          // add the sparse pixel to the sparse frame
          sparseFrame->addSparsePixel( sparsePixel.get() );

        }
        // end of data copying


}

void EUTelMAPSdigi::packFEI4( int digiIndex ) 
{

  int debug = 0;

//  streamlog_out ( MESSAGE4 )  << "packing FEI4" << endl;

  if( zsFrame == 0 )
  {
     streamlog_out( ERROR5 ) << " zsFrame is absent. emergency exit" << endl;
     exit(-1);
  } 
        // now I have to prepare a temporary storage for the
        // sparse pixel
 	auto_ptr< EUTelSparseDataImpl< EUTelAPIXSparsePixel   > > sparseFrame( new EUTelSparseDataImpl< EUTelAPIXSparsePixel   > ( zsFrame ) );
        auto_ptr< EUTelAPIXSparsePixel> sparsePixel( new EUTelAPIXSparsePixel );

        for ( _pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); ++_pixelIterator ) 
        {
          if(debug>1) streamlog_out ( MESSAGE5 ) << " pixelIterator:" <<  _pixelIterator->getIndexAlongL() << " " <<  _pixelIterator->getIndexAlongW()  << endl;  
          // move the pixel information from the pixelIterator to
          // the sparsePixel
          sparsePixel->setXCoord( _pixelIterator->getIndexAlongL() );
          sparsePixel->setYCoord( _pixelIterator->getIndexAlongW() );

	  // Take readout type and ADC range into account

          int pixelValue = 1;        // For binary readout all pixels
				     // above threshold are set to 1

          if( _pixelReadoutType[digiIndex] == 1 )  // digital readout
	    {
            double pixelCharge = _pixelIterator->getCharge();

            pixelValue = static_cast< int >(pixelCharge);

            if(pixelCharge<0.) pixelValue--; 

 
           if(_adcRange[digiIndex] > 0)
	      {
              if(pixelValue > _adcRange[digiIndex])pixelValue=_adcRange[digiIndex];
              if(pixelValue < 0 )pixelValue=0;
	      }
	    }

          if (debug>1)
                     streamlog_out( MESSAGE4 ) <<  "pixel charge " << pixelValue << " X " <<  _pixelIterator->getIndexAlongL() << " " <<  _pixelIterator->getIndexAlongW() << endl;
          sparsePixel->setSignal( pixelValue );

          // add the sparse pixel to the sparse frame
          sparseFrame->addSparsePixel( sparsePixel.get() );

        }
        // end of data copying

}


void EUTelMAPSdigi::end() {

  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelMAPSdigi::fillHist1D(string HistName, int layerIndex, double xVal, double wVal)
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  string tempHistoName;
  if ( _histogramSwitch ) {
    {
      stringstream ss;
      ss << HistName << "_" << layerIndex ;
      tempHistoName = ss.str();
    }
    if ( AIDA::IHistogram1D* histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[ tempHistoName ]) )
      histo->fill(xVal,wVal);
    else {
      streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                << ".\nDisabling histogramming from now on " << endl;
      _histogramSwitch = false;
    }

  }
#endif

}

void EUTelMAPSdigi::fillHist2D(string HistName, int layerIndex, double xVal, double yVal, double wVal)
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  string tempHistoName;
  if ( _histogramSwitch ) {
    {
      stringstream ss;
      ss << HistName << "_" << layerIndex ;
      tempHistoName = ss.str();
    }
    if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[ tempHistoName ]) )
      histo->fill(xVal,yVal,wVal);
    else {
      streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                << ".\nDisabling histogramming from now on " << endl;
      _histogramSwitch = false;
    }

  }
#endif

}

void EUTelMAPSdigi::fillProf1D(string HistName, int layerIndex, double xVal, double wVal)
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  string tempHistoName;
  if ( _histogramSwitch ) {
    {
      stringstream ss;
      ss << HistName << "_" << layerIndex ;
      tempHistoName = ss.str();
    }
    if ( AIDA::IProfile1D* histo = dynamic_cast<AIDA::IProfile1D*>(_aidaHistoMap[ tempHistoName ]) )
      histo->fill(xVal,wVal);
    else {
      streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                << ".\nDisabling histogramming from now on " << endl;
      _histogramSwitch = false;
    }

  }
#endif

}

void EUTelMAPSdigi::bookHist1D(string HistoName, string HistoTitle, int iDet,
                               int xNBin, double xMin, double xMax)
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  if ( _histogramSwitch ) {

    try {
      string tempHistoName;


      {
        stringstream ss ;
        ss <<  HistoName << "_" << iDet ;
        tempHistoName = ss.str();
      }

      AIDA::IHistogram1D * Histo = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(),xNBin, xMin, xMax);
      if ( Histo ) {
        Histo->setTitle(HistoTitle);
        _aidaHistoMap.insert( make_pair( tempHistoName, Histo ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << tempHistoName << ".\n"
                                  << "Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

    } catch (lcio::Exception& e ) {

      streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
      string answer;
      while ( true ) {
        streamlog_out ( ERROR1 ) <<  "[q]/[c]" << endl;
        cin >> answer;
        transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
        if ( answer == "q" ) {
          exit(-1);
        } else if ( answer == "c" )
          _histogramSwitch = false;
        break;
      }
    }

  }

#endif

}


void EUTelMAPSdigi::bookHist2D(string HistoName, string HistoTitle, int iDet,
                               int xNBin, double xMin, double xMax, int yNBin, double yMin, double yMax)
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  if ( _histogramSwitch ) {

    try {
      string tempHistoName;


      {
        stringstream ss ;
        ss <<  HistoName << "_" << iDet ;
        tempHistoName = ss.str();
      }

      AIDA::IHistogram2D * Histo = AIDAProcessor::histogramFactory(this)->createHistogram2D( tempHistoName.c_str(),
                                                                                             xNBin, xMin, xMax, yNBin, yMin, yMax );
      if ( Histo ) {
        Histo->setTitle(HistoTitle);
        _aidaHistoMap.insert( make_pair( tempHistoName, Histo ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << tempHistoName << ".\n"
                                  << "Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

    } catch (lcio::Exception& e ) {

      streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
      string answer;
      while ( true ) {
        streamlog_out ( ERROR1 ) <<  "[q]/[c]" << endl;
        cin >> answer;
        transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
        if ( answer == "q" ) {
          exit(-1);
        } else if ( answer == "c" )
          _histogramSwitch = false;
        break;
      }
    }

  }

#endif

}


void EUTelMAPSdigi::bookProf1D(string HistoName, string HistoTitle, int iDet,
                               int xNBin, double xMin, double xMax, double vMin, double vMax)
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  if ( _histogramSwitch ) {

    try {
      string tempHistoName;


      {
        stringstream ss ;
        ss <<  HistoName << "_" << iDet ;
        tempHistoName = ss.str();
      }

      AIDA::IProfile1D * Histo = AIDAProcessor::histogramFactory(this)->createProfile1D( tempHistoName.c_str(),xNBin, xMin, xMax, vMin, vMax);
      if ( Histo ) {
        Histo->setTitle(HistoTitle);
        _aidaHistoMap.insert( make_pair( tempHistoName, Histo ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << tempHistoName << ".\n"
                                  << "Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

    } catch (lcio::Exception& e ) {

      streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
      string answer;
      while ( true ) {
        streamlog_out ( ERROR1 ) <<  "[q]/[c]" << endl;
        cin >> answer;
        transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
        if ( answer == "q" ) {
          exit(-1);
        } else if ( answer == "c" )
          _histogramSwitch = false;
        break;
      }
    }

  }

#endif

}


void EUTelMAPSdigi::bookHistos() {

  streamlog_out ( MESSAGE4 ) <<  "Booking histograms " << endl;

  for (int iDet = 0 ; iDet < _siPlanesParameters->getSiPlanesNumber(); iDet++)
    {

      // 2D position in local frame

      double xMin =  0;
      double xMax =  _siPlanesLayerLayout->getSensitiveSizeX ( iDet );

      double yMin =  0;
      double yMax =  _siPlanesLayerLayout->getSensitiveSizeY ( iDet );

      int xNBin =  200;
      int yNBin =  200;

      bookHist2D(_hitHistoLocalName,"Hit map in the detector local frame of reference",
                 iDet,xNBin, xMin, xMax, yNBin, yMin, yMax );


      // 2D position in telescope frame

      // 2 should be enough because it
      // means that the sensor is wrong
      // by all its size.
      double safetyFactor = 2.0;

      xMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) -
                              ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet ) ));
      xMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) +
                              ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet )));

      yMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) -
                              ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )));
      yMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) +
                              ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )) );

      xNBin = 400;
      yNBin = 400;

      bookHist2D(_hitHistoTelescopeName,"Hit map in the telescope frame of reference",
                 iDet,xNBin, xMin, xMax, yNBin, yMin, yMax );


      // Collected charge deposit

      xNBin =  500;
      xMin =  0;
      xMax = 5000.;

      bookHist1D(_chargeHistoName, "Collected charge", iDet,xNBin, xMin, xMax);

      // Final sensor signal

      int detectorID=_siPlanesLayerLayout->getSensitiveID(iDet);
      int digiIndex = 0;

      if(!_DigiLayerIDs.empty() )digiIndex = _digiIdMap[detectorID];
      xMax*=_depositedChargeScaling[digiIndex]*_adcGain[digiIndex];

      bookHist1D(_signalHistoName, "Sensor signal", iDet,xNBin, xMin, xMax);

      // Pixel multiplicity

      xNBin = 200;

      xMin = -0.5;
      xMax = 199.5;

      bookHist1D(_multipHistoName, "Pixel multiplicity", iDet,xNBin, xMin, xMax);

      // Charge profiles

      if(_fillChargeProfiles)
        {

          xNBin = _integMaxNumberPixelsAlongL*_integMaxNumberPixelsAlongL;
          xMin = 0;
          xMax = xNBin;
          yMin = -1000000;
          yMax = +1000000;

          bookProf1D(_chargeProfileName, "Charge profile, sorted by |Q|", iDet,xNBin, xMin, xMax, yMin, yMax);

          bookProf1D(_weightedProfileName, "Charge profile, sorted by |Q|/r", iDet,xNBin, xMin, xMax, yMin, yMax);
        }


      // Pixel map to show generated cluster shapes

      xNBin =  _siPlanesLayerLayout->getSensitiveNpixelX(iDet);
      yNBin =  _siPlanesLayerLayout->getSensitiveNpixelY(iDet);

      xMin =  0;
      xMax =  xNBin;

      yMin =  0;
      yMax =  yNBin;

      bookHist2D(_pixelHistoName,"Pixel map",
                 iDet,xNBin, xMin, xMax, yNBin, yMin, yMax );

    } // end of loop over layers

}
//
// ===============================================================================
//
//  Private function members
//


double EUTelMAPSdigi::CheckPathLimits()
{
  // For most tracks the whole path should fit inside a sensor

  double pathStart=0.;
  double pathEnd=1.;

  // Check all sensor dimensions !
  // (Probably checking Z would be enough :-)
/*
  for(int idim=0; idim<3 ; idim++)
    {
      double side = _localPitch[idim]/2.;
 
std::cout << "aidim : " << _localSize[idim] << std::endl;
std::cout << "bidim : " << idim << " pathStart : " << pathStart << " localPosition:" << _localPosition[idim] << " localDirection:"<< _localDirection[idim]<<std::endl;
     if( pathStart < side - _localPosition[idim]/(_mokkaPath*_localDirection[idim]) )
       {
        pathStart = side - _localPosition[idim]/(_mokkaPath*_localDirection[idim]);
std::cout << "bidim : " << idim << " pathStart : " << pathStart << " localPosition:" << _localPosition[idim] << " localDirection:"<< _localDirection[idim]<<std::endl;
       }
 
std::cout << "cidim : " << idim << " pathStart : " << pathStart << " localPosition:" << _localPosition[idim] << " localDirection:"<< _localDirection[idim]<<std::endl;
    if(_localPosition[idim]+(pathStart- side   )*_mokkaPath*_localDirection[idim] > _localSize[idim])
       {
        pathStart= side+ (_localSize[idim]- _localPosition[idim])/(_mokkaPath*_localDirection[idim]);
std::cout << "didim : " << idim << " pathStart : " << pathStart << " localPosition:" << _localPosition[idim] << " localDirection:"<< _localDirection[idim]<<std::endl;
       }

//std::cout << " idim : " << idim << " pathEnd   : " << pathEnd   << " localPosition:" << _localPosition[idim] << " localDirection:"<< _localDirection[idim]<<std::endl;
    if(_localPosition[idim]+(pathEnd-side)*_mokkaPath*_localDirection[idim] > _localSize[idim])
       {
        pathEnd= side + (_localSize[idim]- _localPosition[idim])/(_mokkaPath*_localDirection[idim]);
       }
    if(_localPosition[idim]+(pathEnd-side)*_mokkaPath*_localDirection[idim] < 0)
       {
        pathEnd= side - _localPosition[idim]/(_mokkaPath*_localDirection[idim]);
       }  
    }
*/
  //

  if(  pathStart < 0. || pathStart > 1. || pathEnd < 0. || pathEnd > 1.)
    streamlog_out ( WARNING0 ) << "Warning in checking path limits: out of range " << endl
                               << "Start point = " << pathStart << "  End point = " << pathEnd << endl;

  if(  pathStart < 0. ) pathStart = 0.;
  if(  pathStart > 1. ) pathStart = 1.;
  if(  pathEnd < 0. ) pathEnd = 0.;
  if(  pathEnd > 1. ) pathEnd = 1.;

  // To avoid numerical problems

  if(pathEnd < pathStart)
    {
      _mokkaPath=0.;
      return 0.;
    }

  // Apply corrections to track length and track position

  double midShift=(pathEnd+pathStart-1.)/2.;

  for(int idim=0; idim<3 ; idim++)
    _localPosition[idim]+=_mokkaPath*_localDirection[idim]*midShift;

  _mokkaPath*=(pathEnd - pathStart);

  return pathEnd - pathStart;

}


//
// Implementation of new gear parameters: rotation along Euler angles
//  (modified from EUTelHitMaker)

void EUTelMAPSdigi::InvEulerRotation(double* _telPos, double* _gRotation) {
   
 try{
        double doesNothing = _telPos[2];
        doesNothing++;  //Gets rid of compiler warnings
    }
    catch(...)
    {
        throw InvalidParameterException("_telPos[] array can not be accessed \n");
    }

 try{
        double doesNothing = _gRotation[2];
        doesNothing++;  //Gets rid of compiler warnings
    }
    catch(...)
    {
        throw InvalidParameterException("_gRotation[] array can not be accessed \n"); 
    }


    TVector3 _RotatedSensorHit( _telPos[0], _telPos[1], _telPos[2] );


    // Inverse rotation
    // ----------------
    // Order reversed and minus signed added, with respect to
    // EUTelHitMaker implementation. Inverse transformation is
    // approximate only (assuming rotations along X and Y are small). 

    if( TMath::Abs(_gRotation[0]) > 1e-6 )    _RotatedSensorHit.RotateZ( -_gRotation[0] ); // in XY
    if( TMath::Abs(_gRotation[1]) > 1e-6 )    _RotatedSensorHit.RotateY( -_gRotation[1] ); // in ZX 
    if( TMath::Abs(_gRotation[2]) > 1e-6 )    _RotatedSensorHit.RotateX( -_gRotation[2] ); // in ZY

    _telPos[0] = _RotatedSensorHit.X();
    _telPos[1] = _RotatedSensorHit.Y();
    _telPos[2] = _RotatedSensorHit.Z();
 
}


#endif
