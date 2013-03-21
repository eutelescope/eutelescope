// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTELESCOPE.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelExceptions.h"
#include "EUTelMimosa26Generator.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h> 
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// system includes <>
#include <vector>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelMimosa26Generator::EUTelMimosa26Generator () :Processor("EUTelMimosa26Generator") {

  // modify processor description
  _description =
    "EUTelMimosa26Generator transform full frame raw data in ZS calibrated data mimicking the EUDRB behavior. ";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "ZsDataCollectionName",
			   "Input zs data collection",
			   _zsDataCollectionName, string ("zsdata"));

  registerInputCollection (LCIO::TRACKERDATA, "PedestalCollectionName",
			   "Pedestal from the condition file",
			   _pedestalCollectionName, string ("pedestal"));

  registerInputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
			   "Noise from the condition file",
			   _noiseCollectionName, string("noise"));

  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
			   "Pixel status from the condition file",
			   _statusCollectionName, string("status"));

  registerOutputCollection (LCIO::TRACKERDATA, "SparsifiedMimosa26DataCollectionName",
			    "Name of the output sparsified data collection",
			    _sparsifiedMimosa26DataCollectionName, string("data"));

  registerProcessorParameter("SparsePixelType", "Type of sparsified pixel data structure (use SparsePixelType enum)",
			     _pixelType , static_cast<int> ( 1 ) );
  
  vector<float > sigmaCutVecExample;
  sigmaCutVecExample.push_back(2.5);
  sigmaCutVecExample.push_back(2.5);
  sigmaCutVecExample.push_back(2.5);
  sigmaCutVecExample.push_back(2.5);
  sigmaCutVecExample.push_back(2.5);
  sigmaCutVecExample.push_back(2.5);

  registerProcessorParameter("SigmaCut","A vector of float containing for each plane the multiplication factor for the noise",
			     _sigmaCutVec, sigmaCutVecExample);

  
}


void EUTelMimosa26Generator::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

}

void EUTelMimosa26Generator::processRunHeader (LCRunHeader * rdr) {

 
  auto_ptr<EUTelRunHeaderImpl> runHeader (new EUTelRunHeaderImpl(rdr) );
  runHeader->addProcessor( type() );

  _noOfDetector = runHeader->getNoOfDetector();

  // let's check if the number of sigma cut components is the same of
  // the detector number.
  if ( static_cast< unsigned >(_noOfDetector) != _sigmaCutVec.size() ) {
    streamlog_out( WARNING2 ) << "The number of values in the sigma cut does not match the number of detectors\n"
			      << "Changing SigmaCutVec consequently." << endl;
    _sigmaCutVec.resize(_noOfDetector, _sigmaCutVec.back());
  }
 
  // increment the run counter
  ++_iRun;

}


void EUTelMimosa26Generator::processEvent (LCEvent * event) {


  

  if (_iEvt % 10 == 0) 
    streamlog_out( MESSAGE4 ) << "Processing event " 
			      << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
			      << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
			      << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;

  // right place to increment the event counter
  ++_iEvt;   

 
  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {
    streamlog_out( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  }  
 
  try {
    
    LCCollectionVec * inputCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_zsDataCollectionName));
    LCCollectionVec * pedestalCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_pedestalCollectionName));
    LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_noiseCollectionName));
    LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
    CellIDDecoder<TrackerRawDataImpl> cellDecoder(statusCollectionVec);
 
    if (isFirstEvent()) {
      
      // this is the right place to cross check wheter the pedestal and
      // the input data are at least compatible. I mean the same number
      // of detectors and the same number of pixels in each place.
      
      if ( (inputCollectionVec->getNumberOfElements() != pedestalCollectionVec->getNumberOfElements()) ||
	   (inputCollectionVec->getNumberOfElements() != _noOfDetector) ) {
	stringstream ss;
	ss << "Input data and pedestal are incompatible" << endl
	   << "Input collection has    " << inputCollectionVec->getNumberOfElements()    << " detectors," << endl
	   << "Pedestal collection has " << pedestalCollectionVec->getNumberOfElements() << " detectors," << endl
	   << "The expected number is  " << _noOfDetector << endl;
	throw IncompatibleDataSetException(ss.str());
      }

      
      _isFirstEvent = false;
    }
   
    LCCollectionVec * sparsifiedMimosa26DataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
   
   
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
   
      TrackerDataImpl  * zsData   = dynamic_cast < TrackerDataImpl * >(inputCollectionVec->getElementAt(iDetector));
      // TrackerDataImpl     * pedestal  = dynamic_cast < TrackerDataImpl * >   (pedestalCollectionVec->getElementAt(iDetector));
      TrackerDataImpl     * noise     = dynamic_cast < TrackerDataImpl * >   (noiseCollectionVec->getElementAt(iDetector));
      TrackerRawDataImpl     * noise2     = dynamic_cast < TrackerRawDataImpl * >   (noiseCollectionVec->getElementAt(iDetector));
      TrackerRawDataImpl  * status    = dynamic_cast < TrackerRawDataImpl * >(statusCollectionVec->getElementAt(iDetector));
      TrackerDataImpl     * sparsified = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl> sparseDataEncoder(EUTELESCOPE::ZSDATADEFAULTENCODING, sparsifiedMimosa26DataCollection);
      int sensorID = static_cast<int> ( cellDecoder(status)["sensorID"] );
      sparseDataEncoder["sensorID"]        = sensorID;
      sparseDataEncoder["sparsePixelType"] = static_cast<int> ( _pixelType );
      sparseDataEncoder.setCellID(sparsified);
      EUTelMatrixDecoder matrixDecoder(cellDecoder, noise2);
      
      float sigmaCut = _sigmaCutVec[sensorID];

      if ( _pixelType == kEUTelSimpleSparsePixel ) {

	EUTelSparseDataImpl<EUTelSimpleSparsePixel>  sparseData( sparsified ) ;
    
        auto_ptr<EUTelSparseDataImpl<EUTelSimpleSparsePixel > >
          sparseinputData(new EUTelSparseDataImpl<EUTelSimpleSparsePixel> ( zsData ));
        auto_ptr<EUTelSimpleSparsePixel > sparseinputPixel( new EUTelSimpleSparsePixel );
        
        for ( unsigned int iPixel = 0; iPixel < sparseinputData->size(); iPixel++ ) {
        sparseinputData->getSparsePixelAt( iPixel, sparseinputPixel.get() );
        const int index  = matrixDecoder.getIndexFromXY( sparseinputPixel->getXCoord(), sparseinputPixel->getYCoord() );
        const float signal = sparseinputPixel->getSignal();
        const float n = noise->getChargeValues()[ index ];
        const float threshold = sigmaCut * n;
        if ( 
            ( status->getADCValues()[ index ] == EUTELESCOPE::GOODPIXEL ) &&
            signal > threshold
            ) {
          auto_ptr<EUTelSimpleSparsePixel> sparsePixel( new EUTelSimpleSparsePixel );
          sparsePixel->setXCoord( sparseinputPixel->getXCoord() );
          sparsePixel->setYCoord( sparseinputPixel->getYCoord() );
          sparsePixel->setSignal(1);
         
          sparseData.addSparsePixel( sparsePixel.get() );
         
        }

      }
        } else if ( _pixelType == kUnknownPixelType ) {
	throw UnknownDataTypeException("Unknown or not valid sparse pixel type");
      }

      sparsifiedMimosa26DataCollection->push_back( sparsified );
    }


    evt->addCollection(sparsifiedMimosa26DataCollection, _sparsifiedMimosa26DataCollectionName);
 
  } catch (DataNotAvailableException& e) {
    streamlog_out ( ERROR2 ) <<  e.what() << "\n" << "Skipping this event " << endl;
    throw SkipEventException(this);
  }

}

void EUTelMimosa26Generator::check (LCEvent * /*evt*/) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}

void EUTelMimosa26Generator::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;

}

