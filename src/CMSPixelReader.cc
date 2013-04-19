// Version: $Id$
/*========================================================================*/
/*          CMSPixel file converter (RAW->LCIO)                           */
/*          Author: Simon Spannagel (s.spannagel@cern.ch)                 */
/*          Created       23 feb 2012                                     */
/*          Last modified 19 jul 2012                                     */
/*========================================================================*/


#ifdef USE_GEAR

// CMSPixel includes
#include "CMSPixelReader.h"
#include "CMSPixelDecoder.h"

// EUTelescope includes
#include "EUTELESCOPE.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSimpleSparsePixel.h"
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

using namespace std;
using namespace marlin;
using namespace CMSPixel;
using namespace eutelescope;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string CMSPixelReader::_hitMapHistoName             	= "hitMap";
std::string CMSPixelReader::_pulseHeightHistoName          	= "pulseHeight";
#endif


CMSPixelReader::CMSPixelReader ():DataSourceProcessor  ("CMSPixelReader") {

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
  registerOutputCollection (LCIO::TRACKERDATA, "sparseDataCollectionName",
                            "Name of the output sparsified data collection",
                            _sparseDataCollectionName, string("sparse"));
  registerOptionalParameter ("digitalROC", "Choose whether you have digital ROC data to decode or analogue.",
                              _digitalROC, static_cast < bool >(false));

  IntVec std_planes;
  for ( size_t iDetector = 0 ; iDetector < 16; ++iDetector ) {
    std_planes.push_back( iDetector );
  }
  
  registerOptionalParameter ("shufflePlanes", "int vector to hold the telescope plane IDs in the order in which they get the readout token.",
                              _shufflePlanes, std_planes);
                              
  registerOptionalParameter ("lazyDecoding", "Switch between strict and lazy decoding. Lazy decoding will not check for correct TBM trailers and a valid data length. DO NOT USE unless you have a binary file that can't be decoded using strict mode (e.g. corrupted TBM trailers).",
                              _lazyDecoding, static_cast < bool >(false));
  registerOptionalParameter("eventSelection","Select the events to process: 0 - all, 1 - only with correct No. of ROC headers, 2 - only with corr. ROC headers and without bit errors in them.",
                            _event_selection, static_cast< int > ( 0 ) );
  registerOptionalParameter ("haveTBMheaders", "Switch TBM mode on and off. This gives the possibility to read data without real or emulated TBM headers as seen from soem testboard FPGAs. TRUE will look for TBM headers and trailers.",
                              _haveTBM, static_cast < bool >(true));
  registerOptionalParameter ("useIPBus", "Switch from USB readout format (Altera board) to the IPBus readout format (Xilinx board).",
                              _useIPBus, static_cast < bool >(false));
  registerOptionalParameter("debugDecoder","Decoder DEBUG mode level: 0 - off, 1 - debug, 2 - deep debug.",
                            _debugSwitch, static_cast< int > ( 1 ) );
	registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling", _fillHistos, static_cast< bool > ( true ) );
                              
}


void CMSPixelReader::init () {
 
    // Print the processor parameters:
    printParameters ();
    if( _lazyDecoding ) streamlog_out ( WARNING ) << "YOU ARE USING CMSPixelReader WITH lazyDecoding! ONLY DO THIS IF STRICT DECODING FAILS FOR A GIVEN BINARY!" << endl;
    
    // Set processor back to write the header:
    _runNumber = atoi(_srunNumber.c_str());
    _isFirstEvent = true;
    eventNumber = 0;
    iROC = 0;
    
    if(_writeEmptyEvents) flags += FLAG_EMPTYEVENTS;
    if(_lazyDecoding) flags += FLAG_LAZYDECODING;
    if(_haveTBM) flags += FLAG_HAVETBM;
    if(_useIPBus) flags += FLAG_IPBUS;
    
}


void CMSPixelReader::initializeGeometry() {

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


void CMSPixelReader::readDataSource (int Ntrig) 
   {
   
    EUTelEventImpl *event = NULL;
    CMSPixelDecoder *readout;
    // Initialize event vector:
    vector< vector< CMS_event > > event_data;

    // Initialize geometry:
    initializeGeometry();
    if(!_isGeometryReady) throw InvalidGeometryException ("Wrong geometry file?");

   for(unsigned int i = 0; i < _noOfROC; i++) {
        streamlog_out( MESSAGE5 ) << "ROC " << i << " -> TelescopePlane " << _shufflePlanes[i] << endl;
   }
       
    // Open data file
    if(_digitalROC) 
        readout = new CMSPixelDecoderDigital(_fileName.c_str(),&status,_noOfROC,flags,_event_selection,_debugSwitch);
    else
        readout = new CMSPixelDecoderAnalogue(_fileName.c_str(),&status,_noOfROC,flags,_levelsFile.c_str(),_event_selection,_debugSwitch);


    if (status==-1) {   
        streamlog_out ( ERROR5 ) << "Problem opening file " << filename << "." << endl;
        throw DataNotAvailableException("No data file available.");
    }
    else if (status < -1) {
        streamlog_out ( ERROR5 ) << "Problem with address levels file " << levelsFile << "." << endl;
        throw InvalidGeometryException("Addresslevels file does not match (check Number of ROCs in GEAR geometry).");
    }

    
    
    // Loop while we have input data
    while (true)
    {
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
            ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> ( lcHeader.release()) );
            
            // Book histogramms:
    		if ( _fillHistos ) bookHistos();

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

    
        // Clear vector and read next event from file, containing all ROCs / pixels for one trigger:
        event_data.clear();
        
        // Read next event from data source:
        status = readout->get_event(&event_data);
        
        if(status<0) {
            sprintf(exception,"There was a major problem processing event #%i", eventNumber);
            throw EventException(exception);
        }
        else if(status>0) {
            // We didn't write single event - it was just impossible to open/read the data file:
            if(eventNumber == 1) {
                streamlog_out ( WARNING ) << "The data file contained no valid event." << endl;
                throw DataNotAvailableException("Failed to read from data file.");
            }
            // Else: we just reached EOF.
            else break;
        }


        LCCollectionVec * sparseDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);


        // Shuffle the planes according to the given processor parameter:
        std::vector<int> _shuffled = _shufflePlanes;
        for(unsigned int i = 0; i < _noOfROC; i++) {
            if(static_cast< int >(i) != _shuffled[i]) {
                // Swap the two vectors:
                event_data[i].swap(event_data[_shuffled[i]]);
                // Prevent the re-swapping of the elements:
                _shuffled[_shuffled[i]] = _shuffled[i];
            }
        }

        // Initialize ROC counter:
        iROC = 0;
        
        // Now loop over all ROC chips to be read out:
        do { 
            
            TrackerDataImpl * sparse = new TrackerDataImpl();
            CellIDEncoder<TrackerDataImpl> sparseDataEncoder(EUTELESCOPE::ZSDATADEFAULTENCODING, sparseDataCollection);
            sparseDataEncoder["sensorID"]        = iROC;
            sparseDataEncoder["sparsePixelType"] = static_cast<int> ( 1 );
            sparseDataEncoder.setCellID(sparse);

            EUTelSparseDataImpl<EUTelSimpleSparsePixel>  sparseData( sparse ) ;
                    
            // Check if we have an empty event:
            if(!event_data.empty()) {
                for (unsigned int ipx = 0; ipx < event_data[iROC].size(); ipx++) {
                    auto_ptr<EUTelSimpleSparsePixel> sparsePixel( new EUTelSimpleSparsePixel );
                    
                    sparsePixel->setXCoord( static_cast<int> (event_data[iROC][ipx].col ));
                    sparsePixel->setYCoord( static_cast<int> (event_data[iROC][ipx].row ));
                    sparsePixel->setSignal( static_cast<short> (event_data[iROC][ipx].raw));
                    streamlog_out ( DEBUG0 ) << (*sparsePixel.get()) << endl;
            		
            		//fill histogramms if necessary:
            		if ( _fillHistos ) fillHistos ( event_data[iROC][ipx].col, event_data[iROC][ipx].row, event_data[iROC][ipx].raw, iROC );
                    sparseData.addSparsePixel( sparsePixel.get() );
                } 
            }
         
            sparseDataCollection->push_back( sparse );
            
            iROC++;            
        } while (iROC<event_data.size());


        // Start constructing current event:
        event = new EUTelEventImpl();
        event->setDetectorName("CMSPixelTelescope");
        event->setEventType(kDE);
        event->setRunNumber (_runNumber);
        event->setEventNumber (eventNumber);
        LCTime * now = new LCTime();
        event->setTimeStamp(now->timeStamp());
        delete now;
        // ...and write it out:
        event->addCollection (sparseDataCollection, _sparseDataCollectionName);
        ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (event));
        
        delete event;

    }; // end of while (true)


    // Write last event with type EORE
    eventNumber++;    
    event = new EUTelEventImpl();
    event->setDetectorName("CMSPixelTelescope");
    event->setEventType(kEORE);
    event->setRunNumber (_runNumber);
    event->setEventNumber (eventNumber);
    LCTime * now = new LCTime();
    event->setTimeStamp(now->timeStamp());
    delete now;
    ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (event));
    streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    
    streamlog_out ( MESSAGE5 ) << "  Write EORE as event " << event->getEventNumber() << endl;
        
    // Print the readout statistics, invoked by the destructor:
    streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    
    delete readout;
    streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    

    // Delete the EORE event:    
    delete event;
 }


void CMSPixelReader::end () {
   message<MESSAGE5> ("Successfully finished") ;
 }


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void CMSPixelReader::fillHistos (int xCoord, int yCoord, int value, int sensorID) {

			string tempHistoName;
			
			tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
			(dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >(xCoord), static_cast<double >(yCoord), 1.);
			
			tempHistoName = _pulseHeightHistoName + "_d" + to_string( sensorID );
			(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(value);
}


void CMSPixelReader::bookHistos() {
	
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
		hitMapHisto->setTitle("Hit map");

		string pulseHeightTitle = "pulse height ROC" + to_string( iDetector );
		tempHistoName = _pulseHeightHistoName + "_d" + to_string( iDetector );
		AIDA::IHistogram1D * pulseHeightHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 525,-1050,1050);
		_aidaHistoMap.insert(make_pair(tempHistoName, pulseHeightHisto));
		pulseHeightHisto->setTitle(pulseHeightTitle.c_str());
				
	}
	streamlog_out ( MESSAGE5 )  << "end of Booking histograms " << endl;
}
#endif // USE_AIDA || MARLIN_USE_AIDA

#endif // USE_GEAR
