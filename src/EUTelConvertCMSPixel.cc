// Version: $Id: CMSPixelReader.cc 2557 2013-04-19 13:16:02Z spanns $

#ifdef USE_GEAR

// CMSPixel includes
#include "EUTelConvertCMSPixel.h"

// EUTelescope includes
#include "EUTELESCOPE.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

// GEAR includes
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// Marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// Marlin AIDA include
#include "marlin/AIDAProcessor.h"
// AIDA includes
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// LCIO includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <Exceptions.h>

// System includes
#include <vector>
#include <map>
#include <string.h>
#include <memory>

using namespace std;
using namespace marlin;
using namespace CMSPixel;
using namespace eutelescope;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelConvertCMSPixel::_hitMapHistoName              = "hitMap";
std::string EUTelConvertCMSPixel::_hitMapTrigHistoName              = "hitMapTrig";
std::string EUTelConvertCMSPixel::_hitMapCutHistoName              = "hitMapCut";
std::string EUTelConvertCMSPixel::_pulseHeightHistoName        	= "pulseHeight";
std::string EUTelConvertCMSPixel::_triggerPhaseHistoName        = "triggerPhase";
std::string EUTelConvertCMSPixel::_triggerPhaseHitHistoName        = "triggerPhaseHit";
std::string EUTelConvertCMSPixel::_triggerPhaseHitCutHistoName        = "triggerPhaseHitCut";
std::string EUTelConvertCMSPixel::_dcolMonitorHistoName        = "dcolMonitor";
std::string EUTelConvertCMSPixel::_dcolMonitorEvtHistoName        = "dcolMonitorEvt";
#endif


EUTelConvertCMSPixel::EUTelConvertCMSPixel ():DataSourceProcessor  ("EUTelConvertCMSPixel") {

  _description =
    "Reads PSI46 testboard data files and creates LCEvents (zero suppressed data as sparsePixel).\n"
    "Make sure to not specify any LCIOInputFiles in the steering in order to read CMSPixel files.";

  registerProcessorParameter ("FileName", "Input file",
                              _fileName, std::string ("mtb.bin"));

  registerProcessorParameter ("addressLevelsFile", "Address levels calibration file for the TBM and ROC address encoding levels.",
                              _levelsFile, std::string ("addressParameters.dat"));

  registerProcessorParameter ("runNumber", "RunNumber",
                              _srunNumber, std::string("000001"));

  registerProcessorParameter ("writeEmptyEvents", "Enable or disable the writing of events with no hit in all sensor planes.",
                              _writeEmptyEvents, static_cast < bool >(false));

  registerOutputCollection (LCIO::TRACKERDATA, "sparseDataCollectionName", "Name of the output sparsified data collection",
                            _sparseDataCollectionName, string("sparse"));

  registerProcessorParameter ("ROC_type", "Choose the ROC type. Can be: psi46v2, psi46xdb, psi46dig_trig, psi46dig, psi46digv2_b, psi46digv2.",
                              _ROC_type, static_cast < int >(0));

  registerProcessorParameter ("TB_type", "Choose the Testboard type. Can be: PSI_ATB, PSI_DTB, RAL",
                              _TB_type, static_cast < std::string >("RAL"));

  IntVec std_planes;
  for ( size_t iDetector = 0 ; iDetector < 16; ++iDetector ) {
    std_planes.push_back( iDetector );
  }
  
  registerOptionalParameter ("shufflePlanes", "int vector to hold the telescope plane IDs in the order in which they get the readout token.",
			     _shufflePlanes, std_planes);

  IntVec std_cut;
    std_cut.push_back(0);
    std_cut.push_back(0);
    std_cut.push_back(51);
    std_cut.push_back(0);
    std_cut.push_back(79);

  registerOptionalParameter ("cutHitmap", "Gives the possibilty to cut our a rectangle from the hitmap and evaluate the trigger phase only for events with pixel hits in this region. Format is: ROC COL_MIN COL_MAX ROW_MIN ROW_MAX. Example: 5 34 52 33 54 gives a rectange on ROC5 between 34 <= x <= 52; 33 <= y <= 54.",
			     _cutHitmap, std_cut);
                              
  registerOptionalParameter ("haveTBMheaders", "Switch TBM mode on and off. This gives the possibility to read data without real or emulated TBM headers as seen from soem testboard FPGAs. TRUE will look for TBM headers and trailers.",
			   _haveTBM, static_cast < bool >(true));

  //  registerOptionalParameter ("denseData", "Switch between 4/12bit readout (PSI ATB/DTB) and IPBus 16bit readout (RAL board).", _useIPBus, static_cast < bool >(false));

  registerOptionalParameter("debugDecoder","Set decoder verbosity level: QUIET, SUMMARY, ERROR, WARNING, INFO, DEBUG, DEBUG1-4", _debugSwitch, static_cast< std::string > ( "SUMMARY" ) );
  
  registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling", _fillHistos, static_cast< bool > ( true ) );
 
}


void EUTelConvertCMSPixel::init () {
 
  // Print the processor parameters:
  printParameters ();
    
  // Set processor back, initialize the runnumber:
  _runNumber = atoi(_srunNumber.c_str());
  _isFirstEvent = true;
  eventNumber = 0;
  timestamp_event1 = 0;
    
  if(_haveTBM) flags += FLAG_HAVETBM;
  if(_useIPBus) flags += FLAG_16BITS_PER_WORD;

  // Set Decoder Logging level:
  Log::ReportingLevel() = Log::FromString(_debugSwitch);
    
}


void EUTelConvertCMSPixel::initializeGeometry() 
{
  streamlog_out( MESSAGE5 ) << "Initializing geometry" << endl;

  _noOfROC = 0;
  _noOfXPixel = 0;
  _noOfYPixel = 0;	

  _siPlanesParameters  = const_cast< gear::SiPlanesParameters*  > ( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast< gear::SiPlanesLayerLayout* > ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _layerIndexMap.clear();
  for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); ++iLayer ) {
    _layerIndexMap.insert( make_pair( _siPlanesLayerLayout->getID( iLayer ), iLayer ) );
  }

  _noOfROC = _siPlanesLayerLayout->getNLayers();
	
  // We only use identical telescope planes, so reading the parameters from the first should be fine:
  _noOfXPixel = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[0] );
  _noOfYPixel = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[0] );

  if ( _noOfROC == 0 || _noOfXPixel == 0 || _noOfYPixel == 0 ) {
    streamlog_out( WARNING ) << "Unable to initialize the geometry. Please check GEAR file." << endl;
    _isGeometryReady = false;
  } else {
    _isGeometryReady = true;
  }
  streamlog_out( MESSAGE5 ) << "Active SensorPlanes: " << _noOfROC << endl;
  streamlog_out( MESSAGE5 ) << "Pixels in X: " << _noOfXPixel << endl;
  streamlog_out( MESSAGE5 ) << "Pixels in Y: " << _noOfYPixel << endl;

}


void EUTelConvertCMSPixel::readDataSource (int Ntrig) 
{
   
  EUTelEventImpl *evt = NULL;
  CMSPixelFileDecoder *readout;

  // Initialize geometry:
  initializeGeometry();
  if(!_isGeometryReady) throw InvalidGeometryException ("Wrong geometry file?");

  for(unsigned int i = 0; i < _noOfROC; i++) {
    streamlog_out( DEBUG5 ) << "ROC " << i << " -> TelescopePlane " << _shufflePlanes[i] << endl;
  }
       
  // Construct new decoder:
  if(strcmp(_TB_type.c_str(),"RAL") == 0) {
    streamlog_out( DEBUG5 ) << "Constructing RAL testboard decoder..." << endl;
    streamlog_out( DEBUG5 ) << "Parameters: file: " << _fileName << " nROCs: " << _noOfROC 
			    << " flags: " << flags << " rocType: " << _ROC_type << endl;
    readout = new CMSPixelFileDecoderRAL(_fileName.c_str(),_noOfROC,flags,_ROC_type);
  }
  else if(strcmp(_TB_type.c_str(),"PSI_DTB") == 0) {
    streamlog_out( DEBUG5 ) << "Constructing PSI_DTB testboard decoder..." << endl;
    streamlog_out( DEBUG5 ) << "Parameters: file: " << _fileName << " nROCs: " << _noOfROC 
			    << " flags: " << flags << " rocType: " << _ROC_type 
			    << " levels: " << _levelsFile << endl;
    readout = new CMSPixelFileDecoderPSI_DTB(_fileName.c_str(),_noOfROC,flags,_ROC_type,_levelsFile.c_str());
  }
  else if(strcmp(_TB_type.c_str(),"PSI_ATB") == 0) {
    streamlog_out( DEBUG5 ) << "Constructing PSI_ATB testboard decoder..." << endl;
    streamlog_out( DEBUG5 ) << "Parameters: file: " << _fileName << " nROCs: " << _noOfROC 
			    << " flags: " << flags << " rocType: " << _ROC_type
			    << " levels: " << _levelsFile << endl;
    readout = new CMSPixelFileDecoderPSI_ATB(_fileName.c_str(),_noOfROC,flags,_ROC_type,_levelsFile.c_str());
  }
  else {
    throw DataNotAvailableException("Could not determine correct testboard type. Check configuration.");
  }


  // Loop while we have input data, break points set in Decoder call:
  while (true) {

      // Check if it's the first event:
      if(_isFirstEvent) {
	// We are in the first event, so type BORE. Write the run header.
	auto_ptr<IMPL::LCRunHeaderImpl> lcHeader  ( new IMPL::LCRunHeaderImpl );
	auto_ptr<EUTelRunHeaderImpl>    runHeader ( new EUTelRunHeaderImpl (lcHeader.get()) );
	runHeader->addProcessor( type() );
	runHeader->lcRunHeader()->setDescription(" Events read from CMSPixel input file: " + _fileName);
	runHeader->lcRunHeader()->setRunNumber (_runNumber);
	runHeader->setHeaderVersion (0.0001);
	runHeader->setDataType (EUTELESCOPE::CONVDATA);
	runHeader->setDateTime ();
	runHeader->addIntermediateFile (_fileName);
	runHeader->addProcessor (_processorName);
	runHeader->setNoOfDetector(_noOfROC);
	runHeader->setMinX(IntVec(_noOfROC, 0));
	runHeader->setMaxX(IntVec(_noOfROC, _noOfXPixel - 1));
	runHeader->setMinY(IntVec(_noOfROC, 0));
	runHeader->setMaxY(IntVec(_noOfROC, _noOfYPixel - 1));
	runHeader->lcRunHeader()->setDetectorName("CMSPixelTelescope");

	// Process the run header:
	ProcessorMgr::instance()->processRunHeader(static_cast<lcio::LCRunHeader*>(lcHeader.release()));
            
	// Book histogramms:
	if(_fillHistos) bookHistos();

	_isFirstEvent = false;
      }

      // Trigger counter:
      if(eventNumber >= Ntrig) {
	streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;
	streamlog_out ( MESSAGE5 ) << "  End of processing: reached MaxRecordNumber (" << Ntrig << ")" << endl;
	streamlog_out ( MESSAGE5 ) << "  If you want to process more events check your steerfile." << endl;                                    
	break;
      }
      eventNumber++;

    
      // Initialize pixel vector:
      std::vector< pixel > event_data;
      CMSPixel::timing evt_timing;

      // Read next event from file, containing all ROCs / pixels for one trigger:
      status = readout->get_event(&event_data,evt_timing);

      // Fill the trigger phase histogram for all events, even empty ones:
      (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[_triggerPhaseHistoName]))->fill((int)evt_timing.trigger_phase);

      // Get timestamp from first event:
      if(timestamp_event1 == 0) timestamp_event1 = evt_timing.timestamp;

      if(status <= DEC_ERROR_NO_MORE_DATA) {
	streamlog_out (ERROR) << "Decoder returned error " << status << std::endl;
	// We didn't write single event - it was just impossible to open/read the data file:
	if(eventNumber == 1) {
	  streamlog_out ( WARNING ) << "The data file contained no valid event." << endl;
	  readout->statistics.print();
	  throw DataNotAvailableException("Failed to read from data file.");
	}
	// Else: we just reached EOF.
	else break;
      }
      else if(status <= DEC_ERROR_NO_TBM_HEADER ) {
	streamlog_out (WARNING5) << "There was an exception while processing event " << eventNumber << ". Will continue with next event." << std::endl;
	continue;
      }
      else if(status <= DEC_ERROR_INVALID_ROC_HEADER) {
	streamlog_out (DEBUG5) << "Issue with ROC header or pixel address in event " << eventNumber << ". Event will be used anyway." << std::endl;
      }
      else if(!_writeEmptyEvents && status == DEC_ERROR_EMPTY_EVENT) {
	streamlog_out (DEBUG5) << "Event " << eventNumber << " is empty. Continuing with next." << std::endl;
	continue;
      }
      
      streamlog_out(DEBUG5) << "Event read: " << event_data.size() << " hits" << std::endl;
      for(std::vector<pixel>::const_iterator it = event_data.begin(); it < event_data.end(); ++it) {
	streamlog_out(DEBUG5) << "ROC" << (*it).roc << " x" << (*it).col << " y" << (*it).row 
			      << " ph" << (*it).raw << std::endl;
      }

      LCCollectionVec * sparseDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);


      // Fill the trigger phase histograms of events containing a hit:
      if(!event_data.empty()) (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[_triggerPhaseHitHistoName]))->fill((int)evt_timing.trigger_phase);
      // Initialize bool to write trigger phase for events within the cut boundaries:
      bool cut_done = false;

      // Initialize iterator ROC counter:
      std::vector<pixel>::const_iterator it = event_data.begin();
        
      // Now loop over all ROC chips to be read out:
      for(uint16_t iROC = 0; iROC < _noOfROC; iROC++) {
	//      while(it != event_data.end()) {

	//streamlog_out(DEBUG5) << "Processing ROC " << (*it).roc << std::endl;
	streamlog_out(DEBUG5) << "Processing ROC " << iROC << std::endl;
	//iROC = (*it).roc;

	// Prepare the sensor's header:
	TrackerDataImpl * sparse = new TrackerDataImpl();
	CellIDEncoder<TrackerDataImpl> sparseDataEncoder(EUTELESCOPE::ZSDATADEFAULTENCODING, sparseDataCollection);
	sparseDataEncoder["sensorID"]        = iROC;//(*it).roc;
	sparseDataEncoder["sparsePixelType"] = static_cast<int>(1);
	sparseDataEncoder.setCellID(sparse);
	EUTelSparseDataImpl<EUTelSimpleSparsePixel> sparseData(sparse) ;
        
	// Now add all the pixel hits to that sensor:
	while(it != event_data.end() && iROC == (*it).roc) {
	  streamlog_out(DEBUG5) << "At ROC " << (*it).roc << ", still having hit data..." << std::endl;
	  // Create a new pixel to be stored:
	  auto_ptr<EUTelSimpleSparsePixel> sparsePixel(new EUTelSimpleSparsePixel);
	  sparsePixel->setXCoord( static_cast<int>((*it).col ));
	  sparsePixel->setYCoord( static_cast<int>((*it).row ));
	  sparsePixel->setSignal( static_cast<short>((*it).raw));
	  streamlog_out(DEBUG0) << (*sparsePixel.get()) << endl;
	  
	  // Fill histogramms if necessary:
	  if(_fillHistos) fillHistos((*it).col, (*it).row, (*it).raw, (*it).roc, evt_timing.timestamp, evt_timing.trigger_phase, eventNumber);

	  // Fill the trigger phase histogram (hit/cut) if we have pixel hits within the cut boundaries:
	  if((*it).roc == _cutHitmap[0] 
	     && (*it).col >= _cutHitmap[1] && (*it).col <= _cutHitmap[2]
	     && (*it).row >= _cutHitmap[3] && (*it).row <= _cutHitmap[4]) {
	      string tempHistoName = _hitMapCutHistoName + "_d" + to_string((*it).roc);
	      (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >((*it).col), static_cast<double >((*it).row), 1.);

	      if(!cut_done) {
		(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[_triggerPhaseHitCutHistoName]))->fill((int)evt_timing.trigger_phase);
		cut_done = true;
	      }
	  }

	  sparseData.addSparsePixel(sparsePixel.get());
	  
	  // Move on to next pixel hit:
	  ++it;
	}
	
	streamlog_out(DEBUG5) << sparseData.size() << " pixel hits stored for this ROC." << std::endl;
	sparseDataCollection->push_back( sparse );
            
      }

      streamlog_out(DEBUG5) << "We stored all pixel hits on " << sparseDataCollection->size() 
			    << " ROCs, finalising event now..." << std::endl;

      // Start constructing current event:
      evt = new EUTelEventImpl();
      evt->setDetectorName("CMSPixelTelescope");
      evt->setEventType(kDE);
      evt->setRunNumber (_runNumber);
      evt->setEventNumber (eventNumber);
      LCTime * now = new LCTime();
      evt->setTimeStamp(evt_timing.timestamp);
      delete now;
      // ...and write it out:
      evt->addCollection (sparseDataCollection, _sparseDataCollectionName);
      ProcessorMgr::instance()->processEvent(static_cast<LCEventImpl*> (evt));
     
      delete evt;

    }; // end of while (true)


  // Write last event with type EORE
  eventNumber++;    
  evt = new EUTelEventImpl();
  evt->setDetectorName("CMSPixelTelescope");
  evt->setEventType(kEORE);
  evt->setRunNumber (_runNumber);
  evt->setEventNumber (eventNumber);
  LCTime * now = new LCTime();
  evt->setTimeStamp(now->timeStamp());
  delete now;
  ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (evt));
  streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    
  streamlog_out ( MESSAGE5 ) << "  Write EORE as event " << evt->getEventNumber() << endl;
        
  // Print the readout statistics, invoked by the destructor:
  streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    
  delete readout;
  streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    

  // Delete the EORE event:    
  delete evt;
}


void EUTelConvertCMSPixel::end () {
  message<MESSAGE5> ("Successfully finished") ;
}


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void EUTelConvertCMSPixel::fillHistos (int xCoord, int yCoord, int value, int sensorID, int64_t timestamp, int triggerphase, int evt) {

  string tempHistoName;
			
  tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
  (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >(xCoord), static_cast<double >(yCoord), 1.);

  tempHistoName = _hitMapTrigHistoName + "_d" + to_string( sensorID ) + "_trph" + to_string( triggerphase );
  (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >(xCoord), static_cast<double >(yCoord), 1.);
			
  tempHistoName = _pulseHeightHistoName + "_d" + to_string( sensorID );
  (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(value);

  // RAL IPBus boards report 1us ticks as timestamp. Here we display 10 seconds;
  double time = (timestamp - timestamp_event1)/1E4;
  tempHistoName = _dcolMonitorHistoName + "_d" + to_string( sensorID );
  (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double>(time), static_cast<double >(xCoord), 1.);

  tempHistoName = _dcolMonitorEvtHistoName + "_d" + to_string( sensorID );
  (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double>(evt), static_cast<double >(xCoord), 1.);

}


void EUTelConvertCMSPixel::bookHistos() {
	
  streamlog_out ( MESSAGE5 )  << "Booking histograms " << endl;

  string tempHistoName;
  string basePath;
  for (unsigned int iDetector = 0; iDetector < _noOfROC; iDetector++) {

    basePath = "detector_" + to_string( iDetector );
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath.append("/");

    tempHistoName = _hitMapHistoName + "_d" + to_string( iDetector );
    AIDA::IHistogram2D * hitMapHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), _noOfXPixel, 0, _noOfXPixel, _noOfYPixel, 0, _noOfYPixel);
    _aidaHistoMap.insert(make_pair(tempHistoName, hitMapHisto));
    hitMapHisto->setTitle("Hit map;pixels X;pixels Y");

    for(unsigned int iTriggerPhase = 0; iTriggerPhase < 8; iTriggerPhase++) {
      tempHistoName = _hitMapTrigHistoName + "_d" + to_string( iDetector ) + "_trph" + to_string( iTriggerPhase );
      AIDA::IHistogram2D * hitMapTrigHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), _noOfXPixel, 0, _noOfXPixel, _noOfYPixel, 0, _noOfYPixel);
      _aidaHistoMap.insert(make_pair(tempHistoName, hitMapTrigHisto));
      string temptitle = "Hit map for trigger phase " + to_string( iTriggerPhase ) + "; pixels X; pixels Y";
      hitMapTrigHisto->setTitle(temptitle);
    }

    tempHistoName = _hitMapCutHistoName + "_d" + to_string( iDetector );
    AIDA::IHistogram2D * hitMapCutHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), _noOfXPixel, 0, _noOfXPixel, _noOfYPixel, 0, _noOfYPixel);
    _aidaHistoMap.insert(make_pair(tempHistoName, hitMapCutHisto));
    hitMapCutHisto->setTitle("Hit map (cut on pixels at scintillator position);pixels X;pixels Y");

    string pulseHeightTitle = "pulse height ROC" + to_string( iDetector ) + ";ADC counts;# events";
    tempHistoName = _pulseHeightHistoName + "_d" + to_string( iDetector );
    AIDA::IHistogram1D * pulseHeightHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 525,-1050,1050);
    _aidaHistoMap.insert(make_pair(tempHistoName, pulseHeightHisto));
    pulseHeightHisto->setTitle(pulseHeightTitle.c_str());

    string dcolMonitorTitle = "DCOL hits over time, ROC" + to_string( iDetector ) + ";time / 10ms;column ID;hits per 10ms";
    tempHistoName = _dcolMonitorHistoName + "_d" + to_string( iDetector );
    AIDA::IHistogram2D * dcolMonitorHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), 1000, 0, 1000, 52, -0.5, 51.5);
    _aidaHistoMap.insert(make_pair(tempHistoName, dcolMonitorHisto));
    dcolMonitorHisto->setTitle(dcolMonitorTitle.c_str());

    string dcolMonitorEvtTitle = "DCOL hits over event #, ROC" + to_string( iDetector ) + ";event # / 100;column ID;hits per 100 events";
    tempHistoName = _dcolMonitorEvtHistoName + "_d" + to_string( iDetector );
    AIDA::IHistogram2D * dcolMonitorEvtHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), 5000, 0, 500000, 52, -0.5, 51.5);
    _aidaHistoMap.insert(make_pair(tempHistoName, dcolMonitorEvtHisto));
    dcolMonitorEvtHisto->setTitle(dcolMonitorEvtTitle.c_str());
				
  }

  tempHistoName = _triggerPhaseHistoName;
  AIDA::IHistogram1D * triggerPhaseHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName.c_str(), 10,-1,8);
  _aidaHistoMap.insert(make_pair(tempHistoName, triggerPhaseHisto));
  triggerPhaseHisto->setTitle("Trigger Phase;phase bits; # events");

  tempHistoName = _triggerPhaseHitHistoName;
  AIDA::IHistogram1D * triggerPhaseHitHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName.c_str(), 10,-1,8);
  _aidaHistoMap.insert(make_pair(tempHistoName, triggerPhaseHitHisto));
  triggerPhaseHitHisto->setTitle("Trigger Phase of events w/ pixel hit;phase bits; # events w/ pixel hits");

  tempHistoName = _triggerPhaseHitCutHistoName;
  AIDA::IHistogram1D * triggerPhaseHitCutHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D(tempHistoName.c_str(), 10,-1,8);
  _aidaHistoMap.insert(make_pair(tempHistoName, triggerPhaseHitCutHisto));
  triggerPhaseHitCutHisto->setTitle("Trigger Phase of events w/ pixel hit at the scintillator position;phase bits; # events w/ pixel hits");

  streamlog_out ( MESSAGE5 )  << "end of Booking histograms " << endl;
}
#endif // USE_AIDA || MARLIN_USE_AIDA

#endif // USE_GEAR
