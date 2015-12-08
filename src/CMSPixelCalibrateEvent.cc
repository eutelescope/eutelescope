/*========================================================================*/
/*          CMSPixel CalibrateEvent (Calibration of raw detector data)    */
/*          Author: Simon Spannagel                                       */
/*                (simon.spannagel@student.kit.edu or s.spannagel@cern.ch)*/
/*          Created       23 feb 2012                                     */
/*          Last modified 04 apr 2012                                     */
/*========================================================================*/

// based on EUTelCalibrateEvent:
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

#ifdef USE_GEAR

#include "CMSPixelCalibrateEvent.h"

// EUTelescope includes
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelGenericSparsePixel.h"

// Marlin includes
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// LCIO includes
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCEventImpl.h>

// GEAR includes
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

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

// system includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <memory>
#include <TMath.h>
#include <vector>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string CMSPixelCalibrateEventProcessor::_pulseVcalHistoName             = "pulseVcal";
std::string CMSPixelCalibrateEventProcessor::_pulseHeightHistoName           = "pulseHeight";
std::string CMSPixelCalibrateEventProcessor::_NANMapHistoName                = "NANmap";
#endif


CMSPixelCalibrateEventProcessor::CMSPixelCalibrateEventProcessor () :Processor("CMSPixelCalibrateEventProcessor") {

    _description =
    "CMSPixelCalibrateEventProcessor calibrates the input data according to the tanh fit performed by a module full test with psi46expert.";

    registerInputCollection (LCIO::TRACKERDATA, "sparseDataCollectionName",
                           "Input zero suppressed data collection",
                           _sparseDataCollectionName, string ("sparse"));

    registerOutputCollection (LCIO::TRACKERDATA, "CalibratedDataCollectionName",
                            "Name of the output calibrated data collection",
                            _calibratedDataCollectionName, string("data"));

    registerProcessorParameter("SparsePixelType", "Type of sparsified pixel data structure (use SparsePixelType enum)",
                             _pixelType , static_cast<int> ( 1 ) );

    registerProcessorParameter ("calibrationFile", "Calibration file prefix containing the p0-p3 parameters for the Tanh fit. Use %i for the ROC number.",
                              _calibrationFile, std::string ("phCalibrationFit60_C.dat"));
    registerProcessorParameter ("calibrationType", "Switch between calibration input data types phCalibration (0)  and Gaintanh calibration (1).",
                              _phCalibration, static_cast< bool > ( 1 ) );
	registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling", _fillHistos, static_cast< bool > ( false ) );
	
}


void CMSPixelCalibrateEventProcessor::initializeGeometry() {

	streamlog_out( MESSAGE4 ) << "Initializing geometry" << endl;

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


void CMSPixelCalibrateEventProcessor::initializeCalibration() throw ( marlin::StopProcessingException ) {

    streamlog_out( MESSAGE5 ) << "Read Calibration Data for ROC ";
    
    for(unsigned int i = 0; i < _noOfROC; i++) {

        std::vector< cal_param > cal_roc;
        char calibrationFile[100];
        std::stringstream cf;
        cf << _calibrationFile << i;
        strcpy(calibrationFile,cf.str().c_str());
        streamlog_out( DEBUG5 ) << "ROC" << i << " File: " << calibrationFile << endl;

        std::ifstream* file = new std::ifstream(calibrationFile);

        if ( !file->is_open() ){
            streamlog_out( WARNING ) << "Unable to initialize calibration for ROC" << i << " - could not open file!" << endl;
            throw StopProcessingException( this ) ;
        }

        // Skip reading labels
        char dummyString[100];
        for ( int iskip = 0; iskip < 15; iskip++ ) *file >> dummyString;

        // Write calibration parameters into the struct:
        while(!file->eof()) {
            cal_param dummycal;
            *file >> dummycal.par0 >> dummycal.par1 >> dummycal.par2 >> dummycal.par3 >> dummyString >> dummyString >> dummyString;
            // Push back into ROC vector:
            cal_roc.push_back(dummycal);            
        }

        streamlog_out( MESSAGE5 ) << i << " ";
        streamlog_out( DEBUG5 ) << endl << "First pixel: " << cal_roc[0].par0 << " " << cal_roc[0].par1 << " " << cal_roc[0].par2 << " " << cal_roc[0].par3 << endl;

        delete file;  

        // Check size of vector:
        if(cal_roc.size() < _noOfXPixel * _noOfYPixel) {
            streamlog_out( WARNING ) << "Unable to initialize calibration for ROC" << i << " - wrong number of pixel calibration data" << endl;
            throw StopProcessingException( this ) ;
        }
        
        calibration.push_back(cal_roc);
        
    } // end looping over noOfROC
    
    streamlog_out( MESSAGE5 ) << endl;
  
}

void CMSPixelCalibrateEventProcessor::initializeGaintanhCalibration() throw ( marlin::StopProcessingException ) {

    streamlog_out( MESSAGE5 ) << "Read calibration data for ROC ";
    
    int icol;
    int irow;
    double am;
    double ho;
    double ga;
    double vo;
    double aa = 1.0;
    
    for(unsigned int i = 0; i < _noOfROC; i++) {

        std::vector< cal_param > cal_roc;
        char calibrationFile[100];
        std::stringstream cf;
	cf << _calibrationFile << i;
	strcpy(calibrationFile,cf.str().c_str());
	streamlog_out( DEBUG5 ) << "ROC" << i << " File: " << calibrationFile << endl;

        std::ifstream* file = new std::ifstream(calibrationFile);

        if ( !file->is_open() ){
            streamlog_out( WARNING ) << "Unable to initialize calibration for ROC" << i << " - could not open file!" << endl;
            throw StopProcessingException( this ) ;
        }

        char dummyString[100];
        // Write calibration parameters into the struct:
        while( *file >> dummyString ){
            cal_param dummycal;
            *file >> icol;
            *file >> irow;
            *file >> am;    //Amax
            *file >> ho;    //horz offset
            *file >> ga;    //gain [ADC/large Vcal]
            *file >> vo;    //vert offset

            dummycal.par0 = am*aa;  //amax
            dummycal.par1 = ga;     //gain
            dummycal.par2 = ho;     //horz
            dummycal.par3 = vo*aa;  //vert
            cal_roc.push_back(dummycal);
        }


        streamlog_out( MESSAGE5 ) << i << " ";
        streamlog_out( DEBUG5 ) << endl << "First pixel: " << cal_roc[0].par0 << " " << cal_roc[0].par2 << " " << cal_roc[0].par1 << " " << cal_roc[0].par3 << endl;

        delete file;  

        // Check size of vector:
        if(cal_roc.size() < _noOfXPixel * _noOfYPixel) {
            streamlog_out( WARNING ) << "Unable to initialize calibration for ROC" << i << " - wrong number of pixel calibration data" << endl;
            throw StopProcessingException( this ) ;
        }
        
        calibration.push_back(cal_roc);
        
    } // end looping over noOfROC
    
    streamlog_out( MESSAGE5 ) << endl;
}

void CMSPixelCalibrateEventProcessor::init () {

    printParameters ();

    // Set the run and event counters to zero:
    _iRun = 0;
    
    // Initialize geometry:
    initializeGeometry();
    if(!_isGeometryReady) throw InvalidGeometryException ("Wrong geometry file?");
        
    // Initialize Calibration:
    if(_phCalibration) initializeCalibration();
    else initializeGaintanhCalibration();
    
    // Book histogramms:
    if ( _fillHistos ) bookHistos();

}


void CMSPixelCalibrateEventProcessor::processRunHeader (LCRunHeader * rdr) {

    auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );
    runHeader->addProcessor( type() );

    // Increment the run counter
    ++_iRun;

}



void CMSPixelCalibrateEventProcessor::processEvent (LCEvent * event) {

    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

    if ( evt->getEventType() == kEORE ) {
        streamlog_out ( DEBUG5 ) << "EORE found: nothing else to do." << endl;
        return;
    } else if ( evt->getEventType() == kUNKNOWN ) {
        streamlog_out ( WARNING ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                                   << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }

    try {

        LCCollectionVec * inputCollectionVec    = dynamic_cast<LCCollectionVec*> (evt->getCollection(_sparseDataCollectionName));
        CellIDDecoder<TrackerDataImpl> cellDecoder(inputCollectionVec);

        LCCollectionVec * correctedDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
            
        for (unsigned int iDetector = 0; iDetector < inputCollectionVec->size(); iDetector++) {

            TrackerDataImpl  * sparseData   = dynamic_cast < TrackerDataImpl * >(inputCollectionVec->getElementAt(iDetector));

            TrackerDataImpl     * corrected = new TrackerDataImpl;
            CellIDEncoder<TrackerDataImpl> idDataEncoder(EUTELESCOPE::ZSDATADEFAULTENCODING, correctedDataCollection);
            idDataEncoder["sensorID"] = static_cast<int> (cellDecoder(sparseData)["sensorID"]);
            idDataEncoder["sparsePixelType"] = static_cast<int> (cellDecoder( sparseData )["sparsePixelType"]);
            idDataEncoder.setCellID(corrected);

            EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>  correctedData( corrected ) ;

            auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> > pixelData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>( sparseData ));

            // Loop over all pixels in the sparseData object.
            EUTelGenericSparsePixel Pixel;
            auto_ptr<EUTelGenericSparsePixel> correctedPixel( new EUTelGenericSparsePixel );

            for ( unsigned int iPixel = 0; iPixel < pixelData->size(); iPixel++ ) {
                pixelData->getSparsePixelAt( iPixel, &Pixel);

                correctedPixel->setXCoord( Pixel.getXCoord() );
                correctedPixel->setYCoord( Pixel.getYCoord() );
                
                int iPix = Pixel.getXCoord()*_noOfYPixel + Pixel.getYCoord();

                double corrected;
		bool rangecheck = true;
                if(_phCalibration) rangecheck = calTanH(corrected,Pixel.getSignal(),calibration[iDetector][iPix].par0,calibration[iDetector][iPix].par1,calibration[iDetector][iPix].par2,calibration[iDetector][iPix].par3);
		else rangecheck = calWeibull(corrected,Pixel.getSignal());
                
	        if(rangecheck) {
		  correctedPixel->setSignal( static_cast< short int >(corrected));
                    
                    // Filling histogramms if needed:
               		if ( _fillHistos ) fillHistos ( static_cast< int >(Pixel.getSignal()), static_cast< int >(correctedPixel->getSignal()), iDetector );
               		
               		// Debug output:
                    streamlog_out ( DEBUG5 ) << "evt" << evt->getEventNumber() << " ROC" << iDetector << " Pixel " << Pixel.getXCoord() << " " << Pixel.getYCoord() << ": " << Pixel.getSignal() << " -> " << correctedPixel->getSignal() << endl;
                }
                else {
                    streamlog_out ( WARNING ) << "evt" << evt->getEventNumber() << " ROC" << iDetector << " Pixel " << Pixel.getXCoord() << " " << Pixel.getYCoord() << ": failed to calibrate! Skipping." << endl;
                    // Fill the NAN histogram:
                    
                    string tempHistoName = _NANMapHistoName + "_d" + to_string( iDetector );
			(dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(Pixel.getXCoord(), Pixel.getYCoord(), 1.);
			
                    // This event contains pixel that cannot be calibrated and has to be skipped.
                    throw SkipEventException(this);  
                }
                correctedData.addSparsePixel( correctedPixel.get() );			
                
            }
            correctedDataCollection->push_back( corrected );

        }

        evt->addCollection(correctedDataCollection, _calibratedDataCollectionName);

    } catch (DataNotAvailableException& e) {
          streamlog_out  ( WARNING ) <<  "No input collection found on event " << event->getEventNumber()
                                      << " in run " << event->getRunNumber() << endl;
    }

}



void CMSPixelCalibrateEventProcessor::check (LCEvent * /* evt */ ) {

    // nothing to check here.
}


void CMSPixelCalibrateEventProcessor::end() {

    streamlog_out ( MESSAGE5 ) <<  "Successfully finished" << endl;
}


bool CMSPixelCalibrateEventProcessor::calTanH(double &corr, double y, double p0, double p1, double p2, double p3) {
  // Check for ATanh boundaries, values should be in  (-1,1)
  if(-1 < (y-p3)/p2 && (y-p3)/p2 < 1) {
    corr = (TMath::ATanH((y-p3)/p2) + p1)/p0;
    return true;
  }
  else return false;
}

bool CMSPixelCalibrateEventProcessor::calWeibull(double &corr, double y) {
  //corr = (pow( -log( 1.0 - Ared / ma9 ), 1.0/expo[col][row]) * Gain[col][row] + horz[col][row] ) * keV;
  corr = y;
  streamlog_out( ERROR5 ) << "Calibration mode not yet implemented. Choose different one." << endl;
  return false;
}

bool CMSPixelCalibrateEventProcessor::calLinear(double &corr, double y) {
  corr = y;
  streamlog_out( ERROR5 ) << "Calibration mode not yet implemented. Choose different one." << endl;
  return false;
}


  //bool CMSPixelCalibrateEventProcessor::checkBoundaries(double &corr, double y, double p0, double p1, double p2, double p3) {
  //
  //    if(_phCalibration) {
  //      // Check for ATanh boundaries, values should be in  (-1,1)
  //       if(-1 < (y-p3)/p2 && (y-p3)/p2 < 1) {
  //        corr = (TMath::ATanH((y-p3)/p2) + p1)/p0;
  //        return true;
  //    }
  //    else return false;
  //}
  //else {
  //      double Aout = 1.0; // Aout in no-TBM units!
  //    double Ared = Aout - p3;//sub vert off, see gaintanh2ps.C
  //    double ma9 = p0;
  //
  //    if( Ared >  ma9-1 ) ma9 = Ared + 2;
  //    if( Ared <  -ma9+1 ) ma9 = Ared - 2;
  //
  //    // calibrate into ke units:
  //
  //    //cout << "atanh(" << Ared / ma9 << ") = " << TMath::ATanH( Ared / ma9 ) << endl;
  //
  //    corr = ( TMath::ATanH( Ared / ma9 ) * p1 + p2 ) * 0.45  ; // [ke]
  //    return true;
  ////}
  //
  //}


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void CMSPixelCalibrateEventProcessor::fillHistos (int raw, int calibrated, int sensorID) {

			string tempHistoName;

			tempHistoName = _pulseHeightHistoName + "_d" + to_string( sensorID );
			(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(raw);
			
			tempHistoName = _pulseVcalHistoName + "_d" + to_string( sensorID );
			(dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(calibrated);
}


void CMSPixelCalibrateEventProcessor::bookHistos() {
	
	streamlog_out ( MESSAGE5 )  << "Booking histograms " << endl;

	string tempHistoName;
	string basePath;

	int    signalNBin  = 1999;
	double signalMin   = 0.;
	double signalMax   = 2000.;
		
	for (unsigned int iDetector = 0; iDetector < _noOfROC; iDetector++) {

		basePath = "detector_" + to_string( iDetector );
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		string pulseVcalTitle = "pulse height in Vcal units ROC" + to_string( iDetector );
		tempHistoName = _pulseVcalHistoName + "_d" + to_string( iDetector );
		AIDA::IHistogram1D * pulseVcalHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin,signalMin,signalMax);
		_aidaHistoMap.insert(make_pair(tempHistoName, pulseVcalHisto));
		pulseVcalHisto->setTitle(pulseVcalTitle.c_str());

		string pulseHeightTitle = "pulse height ROC" + to_string( iDetector );
		tempHistoName = _pulseHeightHistoName + "_d" + to_string( iDetector );
		AIDA::IHistogram1D * pulseHeightHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), signalNBin/4,-1000,1000);
		_aidaHistoMap.insert(make_pair(tempHistoName, pulseHeightHisto));
		pulseHeightHisto->setTitle(pulseHeightTitle.c_str());
		
		tempHistoName = _NANMapHistoName + "_d" + to_string( iDetector );
		AIDA::IHistogram2D * NANMapHisto =
		AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), _noOfXPixel, 0, _noOfXPixel, _noOfYPixel, 0, _noOfYPixel);
		_aidaHistoMap.insert(make_pair(tempHistoName, NANMapHisto));
		NANMapHisto->setTitle("Map of pixels with NAN response");
				
	}
	streamlog_out ( MESSAGE5 )  << "end of Booking histograms " << endl;
}
#endif // USE_AIDA || MARLIN_USE_AIDA

#endif
