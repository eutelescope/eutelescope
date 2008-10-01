// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelApplyAlignmentProcessor.cc,v 1.8 2008-10-01 14:54:06 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelApplyAlignmentProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes
#include <AIDA/IHistogram3D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
#endif

#ifdef USE_GEAR
// gear includes <.h>
#include "marlin/Global.h"
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#endif


// system includes <>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

#if defined(USE_GEAR)
using namespace gear;
#endif

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelApplyAlignmentProcessor::_densityPlotBeforeAlignName = "DensityPlotBeforeAlign";
std::string EUTelApplyAlignmentProcessor::_densityPlotAfterAlignName  = "DensityPloAfterAlign";
std::string EUTelApplyAlignmentProcessor::_hitHistoBeforeAlignName    = "HitHistoBeforeAlign";
std::string EUTelApplyAlignmentProcessor::_hitHistoAfterAlignName     = "HitHistoAfterAlign";
#endif

EUTelApplyAlignmentProcessor::EUTelApplyAlignmentProcessor () :Processor("EUTelApplyAlignmentProcessor") {

  // modify processor description
  _description =
    "Apply to the input hit the alignment corrections";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName",
                           "The name of the input hit collection",
                           _inputHitCollectionName, string ("hit"));

  registerInputCollection (LCIO::LCGENERICOBJECT, "AlignmentConstantName",
                           "Alignment constant from the condition file",
                           _alignmentCollectionName, string ("alignment"));

  registerOutputCollection (LCIO::TRACKERHIT, "OutputHitCollectionName",
                            "The name of the output hit collection",
                            _outputHitCollectionName, string("correctedHit"));


  // now the optional parameters
  registerProcessorParameter ("CorrectionMethod",
                              "Available methods are:\n"
                              " 0 --> shift only \n"
                              " 1 --> rotation first \n"
                              " 2 --> shift first ",
                              _correctionMethod, static_cast<int > (1));

  // the histogram on / off switch
  registerOptionalParameter("HistogramSwitch","Enable or disable histograms",
                            _histogramSwitch, static_cast< bool > ( 1 ) );


}


void EUTelApplyAlignmentProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

#if defined(USE_GEAR)
  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));
#endif

#if defined(MARLIN_USE_AIDA) || defined(USE_AIDA)
  _histogramSwitch = true;
#endif

}

void EUTelApplyAlignmentProcessor::processRunHeader (LCRunHeader * rdr) {

  // convert the run header into something mode easy to digest
  auto_ptr<EUTelRunHeaderImpl> runHeader( new EUTelRunHeaderImpl( rdr ) );

  string eudrbGlobalMode = runHeader->getEUDRBMode() ;

  transform( eudrbGlobalMode.begin(), eudrbGlobalMode.end(), eudrbGlobalMode.begin(), ::tolower);

  _hasNZSData = false;

  if ( ( eudrbGlobalMode == "raw2" ) ||
       ( eudrbGlobalMode == "raw3" ) ||
       ( eudrbGlobalMode == "mixed" ) ) {
    _hasNZSData = true;
  }

  _hasZSData = false;

  if ( ( eudrbGlobalMode == "zs" ) ||
       ( eudrbGlobalMode == "mixed" ) ) {
    _hasZSData = true;
  }

  // increment the run counter
  ++_iRun;

}


void EUTelApplyAlignmentProcessor::processEvent (LCEvent * event) {

  if ( _iEvt % 10 == 0 )
    streamlog_out ( MESSAGE4 ) << "Processing event "
                               << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                               << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                               << setfill(' ')
                               << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
  ++_iEvt;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }



  // the cell decoder to get the sensor ID
  CellIDDecoder<TrackerDataImpl> * clusterCellDecoder = NULL; //( originalZSDataCollectionVec );
  CellIDDecoder<TrackerDataImpl> * clusterZSCellDecoder = NULL;


  try {

    LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));

    // I also need the original data collection for ZS and NZS data
    LCCollectionVec * originalDataCollectionVec = NULL;
    LCCollectionVec * originalZSDataCollectionVec = NULL;


    if ( _hasNZSData ) {
      originalDataCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection( "original_data" ) );
      clusterCellDecoder = new CellIDDecoder<TrackerDataImpl>(  originalDataCollectionVec );
    }
    if ( _hasZSData ) {
      originalZSDataCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection( "original_zsdata" ) );
      clusterZSCellDecoder = new CellIDDecoder<TrackerDataImpl>(  originalZSDataCollectionVec );
    }




    if (isFirstEvent()) {

      bookHistos();

      streamlog_out ( MESSAGE ) << "The alignment collection contains: " <<  alignmentCollectionVec->size()
                                << " planes " << endl;

      for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) {

        EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
        _lookUpTable[ alignment->getSensorID() ] = iPos;

      }

#ifndef NDEBUG
      // print out the lookup table
      map< int , int >::iterator mapIter = _lookUpTable.begin();
      while ( mapIter != _lookUpTable.end() ) {
        streamlog_out ( DEBUG ) << "Sensor ID = " << mapIter->first
                                << " is in position " << mapIter->second << endl;
        ++mapIter;
      }
#endif

      _isFirstEvent = false;
    }


    LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);

    for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {

      TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;

      // now we have to understand which layer this hit belongs to.
      LCObjectVec        clusterVec = inputHit->getRawHits();
      TrackerDataImpl  * cluster    = dynamic_cast< TrackerDataImpl *>  ( clusterVec[0] );

      int sensorID;

//       if ( _hasZSData && _hasNZSData ) {
//      // it means that this is a MIXED run still I don't know what
//      // to do
//      streamlog_out ( ERROR ) << "This processor is unable to deal with MIXED data. Sorry for quitting..." << endl;
//      exit(-01);
//       }
      if ( _hasNZSData ) sensorID = (*clusterCellDecoder)( cluster ) ["sensorID"] ;
      if ( _hasZSData  ) sensorID = (*clusterZSCellDecoder)( cluster ) ["sensorID"]   ;


      TrackerHitImpl   * outputHit  = new TrackerHitImpl;
      outputHit->setType( inputHit->getType() );
      outputHit->rawHits() = clusterVec;

      // now that we know at which sensor the hit belongs to, we can
      // get the corresponding alignment constants
      map< int , int >::iterator  positionIter = _lookUpTable.find( sensorID );

      double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
      double   outputPosition[3]  = { 0., 0., 0. };

      if ( positionIter != _lookUpTable.end() ) {

#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) ) && defined(USE_GEAR)
        string tempHistoName;
        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss  << _hitHistoBeforeAlignName << "_" << sensorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( inputPosition[0], inputPosition[1] );
          }
        }


        AIDA::IHistogram3D * histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotBeforeAlignName ] );
        if ( histo3D ) histo3D->fill( inputPosition[0], inputPosition[1], inputPosition[2] );
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }
#endif

        EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * >
          ( alignmentCollectionVec->getElementAt( positionIter->second ) );

        if ( _correctionMethod == 0 ) {

          // this is the shift only case

          outputPosition[0] = inputPosition[0] - alignment->getXOffset();
          outputPosition[1] = inputPosition[1] - alignment->getYOffset();
          outputPosition[2] = inputPosition[2] - alignment->getZOffset();

        } else if ( _correctionMethod == 1 ) {

          // this is the rotation first

          // first the rotation
          outputPosition[0] = inputPosition[0] + alignment->getGamma() * inputPosition[1] + alignment->getBeta() * inputPosition[2] ;
          outputPosition[1] = -1 * alignment->getGamma() * inputPosition[0] + inputPosition[1] + alignment->getAlpha() * inputPosition[2];
          outputPosition[2] = -1 * alignment->getBeta()  * inputPosition[0] + alignment->getAlpha() * inputPosition[1] + inputPosition[2];

          // second the shift
          outputPosition[0] -= alignment->getXOffset();
          outputPosition[1] -= alignment->getYOffset();
          outputPosition[2] -= alignment->getZOffset();

        } else if ( _correctionMethod == 2 ) {

          // this is the translation first

          // first the shifts
          inputPosition[0] -= alignment->getXOffset();
          inputPosition[1] -= alignment->getYOffset();
          inputPosition[2] -= alignment->getZOffset();

          // second the rotation
          outputPosition[0] = inputPosition[0] + alignment->getGamma() * inputPosition[1] + alignment->getBeta() * inputPosition[2] ;
          outputPosition[1] = -1 * alignment->getGamma() * inputPosition[0] + inputPosition[1] + alignment->getAlpha() * inputPosition[2];
          outputPosition[2] = -1 * alignment->getBeta()  * inputPosition[0] + alignment->getAlpha() * inputPosition[1] + inputPosition[2];

        }

#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) ) && defined(USE_GEAR)
        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss  << _hitHistoAfterAlignName << "_" << sensorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( inputPosition[0], inputPosition[1] );
          }
        }


        histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotAfterAlignName ] );
        if ( histo3D ) histo3D->fill( inputPosition[0], inputPosition[1], inputPosition[2] );
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }
#endif


      } else {

        // this hit belongs to a plane whose sensorID is not in the
        // alignment constants. So the idea is to eventually advice
        // the users if running in DEBUG and copy the not aligned hit
        // in the new collection.
        streamlog_out ( DEBUG ) << "Sensor ID " << sensorID << " not found. Skipping alignment for hit "
                                << iHit << endl;

        for ( size_t i = 0; i < 3; ++i )
          outputPosition[i] = inputPosition[i];

      }

      outputHit->setPosition( outputPosition ) ;
      outputCollectionVec->push_back( outputHit );


    }

    evt->addCollection( outputCollectionVec, _outputHitCollectionName );

    delete clusterCellDecoder;
    delete clusterZSCellDecoder;

  } catch (DataNotAvailableException& e) {
    if ( clusterCellDecoder ) delete clusterCellDecoder;
    if ( clusterZSCellDecoder ) delete clusterZSCellDecoder;
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

}

void EUTelApplyAlignmentProcessor::bookHistos() {

#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) ) && defined(USE_GEAR)

  try {
    streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;

    string tempHistoName;

    // histograms are grouped into folders named after the
    // detector. This requires to loop on detector now.
    for (int iDet = 0 ; iDet < _siPlanesParameters->getSiPlanesNumber(); iDet++) {
      
      int sensorID = _siPlanesLayerLayout->getID( iDet ) ;

      string basePath;
      {
        stringstream ss ;
        ss << "plane-" << sensorID;
        basePath = ss.str();
      }
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());
      basePath = basePath + "/";

      // 2 should be enough because it
      // means that the sensor is wrong
      // by all its size.
      double safetyFactor = 2.0;

      double xMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) -
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet ) ));
      double xMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) +
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet )));

      double yMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) -
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )));
      double yMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) +
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )) );

      int xNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelX( iDet );
      int yNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelY( iDet );

      {
        stringstream ss ;
        ss <<  _hitHistoBeforeAlignName << "_" << sensorID ;
        tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * hitHistoBeforeAlign =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( ( basePath + tempHistoName ).c_str(),
                                                                  xNBin, xMin, xMax, yNBin, yMin, yMax );

      if ( hitHistoBeforeAlign ) {
        hitHistoBeforeAlign->setTitle("Hit map in the telescope frame of reference before align");
        _aidaHistoMap.insert( make_pair ( tempHistoName, hitHistoBeforeAlign ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }


      {
        stringstream ss ;
        ss <<  _hitHistoAfterAlignName << "_" << sensorID  ;
        tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * hitHistoAfterAlign =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( ( basePath + tempHistoName ).c_str(),
                                                                  xNBin, xMin, xMax, yNBin, yMin, yMax );

      if ( hitHistoAfterAlign ) {
        hitHistoAfterAlign->setTitle("Hit map in the telescope frame of reference after align");
        _aidaHistoMap.insert( make_pair ( tempHistoName, hitHistoAfterAlign ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }
    }

    // we have to found the boundaries of this histograms. Let's take
    // the outer positions in all directions
    double xMin  =      numeric_limits< double >::max();
    double xMax  = -1 * numeric_limits< double >::max();
    int    xNBin = numeric_limits< int >::min();

    double yMin  =      numeric_limits< double >::max();
    double yMax  = -1 * numeric_limits< double >::max();
    int    yNBin = numeric_limits< int >::min();

    for ( int iPlane = 0 ; iPlane < _siPlanesParameters->getSiPlanesNumber(); ++iPlane ) {

      // x axis
      xMin  = min( _siPlanesLayerLayout->getSensitivePositionX( iPlane ) - ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeX( iPlane )), xMin);
      xMax  = max( _siPlanesLayerLayout->getSensitivePositionX( iPlane ) + ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeX( iPlane )), xMax);
      xNBin = max( _siPlanesLayerLayout->getSensitiveNpixelX( iPlane ), xNBin );

      // y axis
      yMin  = min( _siPlanesLayerLayout->getSensitivePositionY( iPlane ) - ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeY( iPlane )), yMin);
      yMax  = max( _siPlanesLayerLayout->getSensitivePositionY( iPlane ) + ( 0.5  * _siPlanesLayerLayout->getSensitiveSizeY( iPlane )), yMax);
      yNBin = max( _siPlanesLayerLayout->getSensitiveNpixelY( iPlane ), yNBin );

    }


    // since we may still have alignment problem, we have to take a
    // safety factor on the x and y direction especially.
    // here I take something less than 2 because otherwise I will have
    // a 200MB histogram.
    double safetyFactor = 1.2;

    double xDistance = std::abs( xMax - xMin ) ;
    double xCenter   = ( xMax + xMin ) / 2 ;
    xMin  = xCenter - safetyFactor * ( xDistance / 2 );
    xMax  = xCenter + safetyFactor * ( xDistance / 2 );
    xNBin = static_cast< int > ( xNBin * safetyFactor );

    // generate the x axis binning
    vector< double > xAxis;
    double step = xDistance / xNBin;
    for ( int i = 0 ; i < xNBin ; ++i ) {
      xAxis.push_back ( xMin + i * step );
    }

    double yDistance = std::abs( yMax - yMin ) ;
    double yCenter   = ( yMax + yMin ) / 2 ;
    yMin  = yCenter - safetyFactor * ( yDistance / 2 );
    yMax  = yCenter + safetyFactor * ( yDistance / 2 );
    yNBin = static_cast< int > ( yNBin * safetyFactor );

    // generate the y axis binning
    vector< double > yAxis;
    step = yDistance / yNBin;
    for ( int i = 0 ; i < yNBin ; ++i ) {
      yAxis.push_back( yMin + i * step ) ;
    }


    // generate the z axis but not equally spaced!
    double safetyMargin = 10; // this is mm
    vector< double > zAxis;
    for ( int i = 0 ; i < 2 * _siPlanesParameters->getSiPlanesNumber(); ++i ) {
      double zPos =  _siPlanesLayerLayout->getSensitivePositionZ( i/2 );
      zAxis.push_back( zPos - safetyMargin) ;
      ++i;
      zAxis.push_back( zPos + safetyMargin );
    }


    AIDA::IHistogram3D * densityPlot = AIDAProcessor::histogramFactory(this)->createHistogram3D( _densityPlotBeforeAlignName ,
                                                                                                 "Hit position in the telescope frame of reference before align",
                                                                                                 xAxis, yAxis, zAxis, "");

    if ( densityPlot ) {
      _aidaHistoMap.insert( make_pair ( _densityPlotBeforeAlignName, densityPlot ) ) ;
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (_densityPlotBeforeAlignName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
      _histogramSwitch = false;
    }


    densityPlot = AIDAProcessor::histogramFactory(this)->createHistogram3D( _densityPlotAfterAlignName ,
                                                                            "Hit position in the telescope frame of reference after align",
                                                                            xAxis, yAxis, zAxis, "");

    if ( densityPlot ) {
      _aidaHistoMap.insert( make_pair ( _densityPlotAfterAlignName, densityPlot ) ) ;
    } else {
      streamlog_out ( ERROR1 )  << "Problem booking the " << (_densityPlotAfterAlignName) << ".\n"
                                << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
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
#endif // AIDA && GEAR
}


void EUTelApplyAlignmentProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelApplyAlignmentProcessor::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;

}

