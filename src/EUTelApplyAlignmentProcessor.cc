// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Authors
// Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Joerg Behr, Hamburg Uni/DESY  <joerg.behr@desy.de> 
// Slava Libov, DESY <mailto:vladyslav.libov@desy.de>
// Igor Rubinskiy, DESY <mailto:igorrubinsky@gmail.com>
// 
// Version $Id: EUTelApplyAlignmentProcessor.cc,v 1.17 2009-07-30 17:19:19 jbehr Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *   
 *
 */

#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelApplyAlignmentProcessor.h"
#include "EUTelAlignmentConstant.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelAPIXSparsePixel.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelAPIXSparseClusterImpl.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"
#
// ROOT includes:
#include "TVector3.h"
#include "TVector2.h"
#include "TMatrix.h"



#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes ".h"
#include <AIDA/IHistogram3D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
#endif

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
//#include <TrackerHitImpl2.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// ROOT includes ".h"
#include <TVectorD.h>
#include <TMatrixD.h>


// system includes <>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <memory>
#include <string>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

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

  registerOptionalParameter("ReferenceCollection","This is the name of the reference it collection (init at 0,0,0)",
                            _referenceHitCollectionName, static_cast< string > ( "reference" ) );
   registerOptionalParameter("ApplyToReferenceCollection","Do you want the reference hit collection to be corrected by the shifts and tilts from the alignment collection? (default - false )",
                            _applyToReferenceHitCollection, static_cast< bool   > ( false ));
 


  // now the optional parameters
  registerProcessorParameter ("CorrectionMethod",
                              "Available methods are:\n"
                              " 0 --> shift only \n"
                              " 1 --> rotation first \n"
                              " 2 --> shift first ",
                              _correctionMethod, static_cast<int > (1));

  // now the optional parameters
  registerProcessorParameter ("ApplyAlignmentDirection",
                              "Available directinos are:\n"
                              " 0 --> direct  \n"
                              " 1 --> reverse ",
                              _applyAlignmentDirection, static_cast<int > (0));


  // the histogram on / off switch
  registerOptionalParameter("HistogramSwitch","Enable or disable histograms",
                            _histogramSwitch, static_cast< bool > ( 0 ) );


  EVENT::StringVec	_alignmentCollectionSuffixExamples;
  _alignmentCollectionSuffixExamples.push_back("alignment");
  
  registerProcessorParameter ("alignmentCollectionNames",
                            "List of alignment collections that were applied to the DUT",
                            _alignmentCollectionSuffixes, _alignmentCollectionSuffixExamples);

  registerOptionalParameter("DoAlignCollection","Implement geometry shifts and rotations as described in alignmentCollectionName ",
                            _doAlignCollection, static_cast< bool > ( 0 ) );

  registerOptionalParameter("DoGear","Implement geometry shifts and rotations as described in GEAR steering file ",
                            _doGear, static_cast< bool > ( 0 ) );

  // DEBUG parameters :
  // turn ON/OFF debug features 
  registerOptionalParameter("DEBUG","Enable or disable DEBUG mode ",
                            _debugSwitch, static_cast< bool > ( 0 ) );
  registerOptionalParameter("Alpha","Rotation Angle around X axis",
                            _alpha, static_cast< double > ( 0.00 ) );
  registerOptionalParameter("Beta","Rotation Angle around Y axis",
                            _beta, static_cast< double > ( 0.00 ) );
  registerOptionalParameter("Gamma","Rotation Angle around Z axis",
                            _gamma, static_cast< double > ( 0.00 ) );
  registerOptionalParameter("PrintEvents", "Events number to have DEBUG1 printed outs (default=10)",
                            _printEvents, static_cast<int> (10) );

}


void EUTelApplyAlignmentProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
  }

#if defined(MARLIN_USE_AIDA) || defined(USE_AIDA)
  _histogramSwitch = true;
#endif
  _isFirstEvent = true;
  fevent = true;

  _referenceHitVec = 0;
 
}

void EUTelApplyAlignmentProcessor::processRunHeader (LCRunHeader * runHeader) {

//  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl( rdr ) ) ;
//  runHeader->addProcessor( type() );
//
//  // increment the run counter
//  ++_iRun;

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();
  // convert to string
  char buf[256];
  sprintf(buf, "%i", runNr);
  std::string runNr_str(buf);

  message<MESSAGE> ( log() << "Processing run header " << _nRun
                     << ", run nr " << runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  //  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  message<MESSAGE> ( log() << detectorName << " : " << detectorDescription ) ;

  // pick up correct alignment collection
	_alignmentCollectionNames.clear();
	for (unsigned i = 0; i < _alignmentCollectionSuffixes.size(); i++) {
		std::string	temp = _alignmentCollectionSuffixes[i];
		_alignmentCollectionNames.push_back(temp);
		cout << _alignmentCollectionNames[i] << endl;
	}

}


void EUTelApplyAlignmentProcessor::processEvent (LCEvent * event) {

    if( _alignmentCollectionNames.size() <= 0 )
    {
        streamlog_out ( ERROR ) << "Alignment collections are UNDEFINED, the processor can not continue. EXIT " << endl;
        throw StopProcessingException(this);       
    }
    else
    {    
 
        if ( fevent )
            streamlog_out ( MESSAGE4 ) << "Processing run  "  
                << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                << " Number of defined alignment collection is  " << _alignmentCollectionNames.size()
                << endl;
        
        if ( fevent && _applyToReferenceHitCollection ) 
            _referenceHitVec     = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
 
/*
        if ( fevent && _applyToReferenceHitCollection ) 
            { 
               _referenceHitVecAligned = new LCCollectionVec( LCIO::LCGENERICOBJECT );
               event->addCollection( _referenceHitVecAligned, "refhitaligned"                        );
            }
 */
       
        for ( int i = _alignmentCollectionNames.size() - 1; i >= 0; i--) 
        {
            // read the first available alignment collection
            // CAUTION 1: it might be important to keep the order of alignment collections (if many) given in the opposite direction
            // CAUTION 2: to be controled via steering files
            //
            _alignmentCollectionName = _alignmentCollectionNames.at(i);
            
            if( GetApplyAlignmentDirection() == 0 )
            {
                if( _doGear )
                {
                    ApplyGear6D(event);
                }
                if( _doAlignCollection )
                {
                    Direct(event);
                }
            }
            else
                if( GetApplyAlignmentDirection() == 1 )
                {
                    if( _doAlignCollection )
                    {
                        Reverse(event);     
                    }
                    if( _doGear )
                    {
                        RevertGear6D(event);
                    }
                }
                else
                {
                    throw StopProcessingException(this); 
                }
            
        }
 
        if(fevent)
        {
            _isFirstEvent = false;
            fevent = false;
        }
        _iEvt++; 
    }
/* 
    if ( _applyToReferenceHitCollection ) 
    {
      LCCollectionVec * ref    = static_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName ));
      for(size_t ii = 0 ; ii <  ref->getNumberOfElements(); ii++)
      {
        EUTelReferenceHit * refhit = static_cast< EUTelReferenceHit*> ( ref->getElementAt(ii) ) ;
        printf("FIN sensorID: %5d dx:%5.3f dy:%5.3f dz:%5.3f  alfa:%5.3f beta:%5.3f gamma:%5.3f \n",
        refhit->getSensorID(   ),                      
        refhit->getXOffset(    ),
        refhit->getYOffset(    ),
        refhit->getZOffset(    ),
        refhit->getAlpha(),
        refhit->getBeta(),
        refhit->getGamma()    );
      }
    }
*/
}

void EUTelApplyAlignmentProcessor::ApplyGear6D( LCEvent *event) 
{
 

  if ( _iEvt % 100 == 0 )
    streamlog_out ( MESSAGE4 ) << "Processing event  (ApplyGear6D) "
                               << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                               << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                               << setfill(' ')
                               << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
//++_iEvt;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) 
  {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  }
  else if ( evt->getEventType() == kUNKNOWN ) 
  {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


  try 
  {

    LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    
    if (fevent) 
    {

     
#ifndef NDEBUG
        // print out the lookup table
        map< int , int >::iterator mapIter = _lookUpTable.begin();
        while ( mapIter != _lookUpTable.end() ) 
        {
            streamlog_out ( DEBUG ) << "Sensor ID = " << mapIter->first
                                    << " is in position " << mapIter->second << endl;
            ++mapIter;
        }
#endif
    }
    
    
    LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);

    for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) 
    {

      TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;

      // now we have to understand which layer this hit belongs to.
      int sensorID = guessSensorID( inputHit );

      if ( _conversionIdMap.size() != (unsigned) _siPlanesParameters->getSiPlanesNumber() ) 
      {
          // first of all try to see if this sensorID already belong to
          if ( _conversionIdMap.find( sensorID ) == _conversionIdMap.end() ) 
          {
              // this means that this detector ID was not already inserted,
              // so this is the right place to do that
          
              for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) 
              {
                  if ( _siPlanesLayerLayout->getID(iLayer) == sensorID ) 
                  {
                      _conversionIdMap.insert( make_pair( sensorID, iLayer ) );
                      break;
                  }
              }
          }
      }

      int layerIndex   = _conversionIdMap[sensorID];

      // determine z position of the plane
	  // 20 December 2010 @libov
      float	z_sensor = 0;
	  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) 
      {
          if (sensorID == _siPlanesLayerLayout->getID( iPlane ) ) 
          {
              z_sensor = _siPlanesLayerLayout -> getSensitivePositionZ( iPlane ) + 0.5 * _siPlanesLayerLayout->getSensitiveThickness( iPlane );
              break;
          }
      }

      // copy the input to the output, at least for the common part
      TrackerHitImpl   * outputHit  = new TrackerHitImpl;
      outputHit->setType( inputHit->getType() );
      outputHit->rawHits() = inputHit->getRawHits();


      double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
      double   outputPosition[3]  = { 0., 0., 0. };

      if ( 1==1 ) 
      {

          double telPos[3]    = {0., 0., 0.};
          double gRotation[3] = { 0., 0., 0.}; // not rotated

          if ( _debugSwitch )
          {
              telPos[0] = 0.;
              telPos[1] = 0.;
              telPos[2] = 0.; 
              gRotation[0] = _alpha;
              gRotation[1] = _beta ;
              gRotation[2] = _gamma;
          }
          else
          {
              gRotation[0]    = _siPlanesLayerLayout->getLayerRotationXY(layerIndex); // Euler alpha ;
              gRotation[1]    = _siPlanesLayerLayout->getLayerRotationZX(layerIndex); // Euler alpha ;
              gRotation[2]    = _siPlanesLayerLayout->getLayerRotationZY(layerIndex); // Euler alpha ;

              // input angles are in DEGREEs !!!
              // translate into radians
    
              gRotation[0]  =   gRotation[0] *3.1415926/180.; // 
              gRotation[1]  =   gRotation[1] *3.1415926/180.; //
              gRotation[2]  =   gRotation[2] *3.1415926/180.; //

              telPos[0]  =  _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
              telPos[1]  =  _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
              telPos[2]  =  _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
              
          }
 
      
          if( _iEvt < _printEvents )
          {
                if ( _debugSwitch ) 
                {
                    streamlog_out ( MESSAGE )  << "Debugmode ON " << endl;                                   
                }
                
                streamlog_out ( MESSAGE )  << "_applyGear6D " << endl;
                streamlog_out ( MESSAGE )  << " telPos[0] = " << telPos[0]  << endl;
                streamlog_out ( MESSAGE )  << " telPos[1] = " << telPos[1]  << endl;
                streamlog_out ( MESSAGE )  << " telPos[2] = " << telPos[2]  << endl;
                streamlog_out ( MESSAGE )  << " gRotation[0] = " << gRotation[0]  << endl;
                streamlog_out ( MESSAGE )  << " gRotation[1] = " << gRotation[1]  << endl;
                streamlog_out ( MESSAGE )  << " gRotation[2] = " << gRotation[2]  << endl;
          }

          // rotations first
 
          outputPosition[0] = inputPosition[0];
          outputPosition[1] = inputPosition[1];
          outputPosition[2] = inputPosition[2] - z_sensor;
          
          _EulerRotation( sensorID,outputPosition, gRotation);

          // then the shifts
//          outputPosition[0]  = telPos[0];
//          outputPosition[1]  = telPos[1];
//          outputPosition[2]  = telPos[2];

          outputPosition[2] += z_sensor;
      }
      else
      {
          // this hit belongs to a plane whose sensorID is not in the
          // alignment constants. So the idea is to eventually advice
          // the users if running in DEBUG and copy the not aligned hit
          // in the new collection.
//          streamlog_out ( DEBUG3 ) << "Sensor ID " << sensorID << " not found. Skipping alignment for hit "
//                                    << iHit << endl;

          for ( size_t i = 0; i < 3; ++i )
          {
              outputPosition[i] = inputPosition[i];
          }
          
      }

      if ( _iEvt < _printEvents )
      {
         streamlog_out ( MESSAGE ) << "ApplyGear: INPUT: Sensor ID " << sensorID << " " << inputPosition[0] << " " << inputPosition[1] << " " << inputPosition[2] << " " << z_sensor << endl;                
         streamlog_out ( MESSAGE ) << "ApplyGear: OUTPUT:Sensor ID " << sensorID << " " << outputPosition[0] << " " << outputPosition[1] << " " << outputPosition[2] << " " << z_sensor << endl;                
      }

      outputHit->setPosition( outputPosition ) ;
      outputCollectionVec->push_back( outputHit );
    }

    evt->addCollection( outputCollectionVec, _outputHitCollectionName );

  } catch (DataNotAvailableException& e) {
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

}


void EUTelApplyAlignmentProcessor::RevertGear6D( LCEvent *event) 
{
 
  if ( _iEvt % 100 == 0 )
    streamlog_out ( MESSAGE4 ) << "Processing event  (RevertGear6D) "
                               << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                               << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                               << setfill(' ')
                               << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
//  ++_iEvt;


  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) 
  {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  }
  else if ( evt->getEventType() == kUNKNOWN ) 
  {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }


  try 
  {

    LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    
    if (fevent) 
    {

     
#ifndef NDEBUG
        // print out the lookup table
        map< int , int >::iterator mapIter = _lookUpTable.begin();
        while ( mapIter != _lookUpTable.end() ) 
        {
            streamlog_out ( DEBUG ) << "Sensor ID = " << mapIter->first
                                    << " is in position " << mapIter->second << endl;
            ++mapIter;
        }
#endif
    }
    
    
    LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);

    for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) 
    {

      TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;

      // now we have to understand which layer this hit belongs to.
      int sensorID = guessSensorID( inputHit );

      if ( _conversionIdMap.size() != (unsigned) _siPlanesParameters->getSiPlanesNumber() ) 
      {
          // first of all try to see if this sensorID already belong to
          if ( _conversionIdMap.find( sensorID ) == _conversionIdMap.end() ) 
          {
              // this means that this detector ID was not already inserted,
              // so this is the right place to do that
          
              for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) 
              {
                  if ( _siPlanesLayerLayout->getID(iLayer) == sensorID ) 
                  {
                      _conversionIdMap.insert( make_pair( sensorID, iLayer ) );
                      break;
                  }
              }
          }
      }
      
      int layerIndex   = _conversionIdMap[sensorID];

      // determine z position of the plane
      // 20 December 2010 @libov
      float	z_sensor = 0;
      for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) 
      {
          if (sensorID == _siPlanesLayerLayout->getID( iPlane ) ) 
          {
              z_sensor = _siPlanesLayerLayout -> getSensitivePositionZ( iPlane ) + 0.5 * _siPlanesLayerLayout->getSensitiveThickness( iPlane );
              break;
          }
      }

// retrieve the refhit cooridantes (eventual offset of the sensor) 
      double x_refhit = 0.; 
      double y_refhit = 0.; 
      double z_refhit = 0.; 
      
      if( _referenceHitVec == 0 )
      {
        streamlog_out(MESSAGE) << "_referenceHitVec is empty" << endl;
      }
      else
      {
//        streamlog_out(MESSAGE) << "reference Hit collection name : " << _referenceHitCollectionName << endl;
 
       for(size_t ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
       {
        EUTelReferenceHit * refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
        if( sensorID != refhit->getSensorID() )
        {
       //   streamlog_out(MESSAGE) << "Looping through a varity of sensor IDs" << endl;
          continue;
        }
        else
        {
          x_refhit =  refhit->getXOffset();
          y_refhit =  refhit->getYOffset();
          z_refhit =  refhit->getZOffset();
          if( _iEvt < _printEvents )
          {
            streamlog_out(MESSAGE) << "Sensor ID and Alignment plane ID match!" << endl;
            streamlog_out(MESSAGE) << "x_refhit " << x_refhit  << endl; 
            streamlog_out(MESSAGE) << "y_refhit " << y_refhit  << endl; 
            streamlog_out(MESSAGE) << "z_refhit " << z_refhit  << endl; 
          }
          break;
        }
       }
      }

      // copy the input to the output, at least for the common part
      TrackerHitImpl   * outputHit  = new TrackerHitImpl;
      outputHit->setType( inputHit->getType() );
      outputHit->rawHits() = inputHit->getRawHits();


      double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
      double   outputPosition[3]  = { 0., 0., 0. };

   
      if ( 1==1  ) 
      {

          double telPos[3]    = {0., 0., 0.};
          double gRotation[3] = { 0., 0., 0.}; // not rotated


          if ( _debugSwitch )
          {
              telPos[0] = 0.;
              telPos[1] = 0.;
              telPos[2] = 0.; 
              gRotation[0] = _alpha;
              gRotation[1] = _beta ;
              gRotation[2] = _gamma;
          }
          else
          {
              gRotation[0]    = _siPlanesLayerLayout->getLayerRotationXY(layerIndex); // Euler alpha ;
              gRotation[1]    = _siPlanesLayerLayout->getLayerRotationZX(layerIndex); // Euler alpha ;
              gRotation[2]    = _siPlanesLayerLayout->getLayerRotationZY(layerIndex); // Euler alpha ;

              // input angles are in DEGREEs !!!
              // translate into radians
    
              gRotation[0]  =   gRotation[0] *3.1415926/180.; // 
              gRotation[1]  =   gRotation[1] *3.1415926/180.; //
              gRotation[2]  =   gRotation[2] *3.1415926/180.; //

              telPos[0]  =  _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
              telPos[1]  =  _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
              telPos[2]  =  _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
              
          }
      
          if( _iEvt < _printEvents )
          {
                if ( _debugSwitch ) 
                {
                    streamlog_out ( MESSAGE )  << "Debugmode ON " << endl;                                   
                }
                
                streamlog_out ( MESSAGE )  << "_revertGear6D " << endl;
                streamlog_out ( MESSAGE )  << " telPos[0] = " << telPos[0]  << endl;
                streamlog_out ( MESSAGE )  << " telPos[1] = " << telPos[1]  << endl;
                streamlog_out ( MESSAGE )  << " telPos[2] = " << telPos[2]  << endl;
                streamlog_out ( MESSAGE )  << " gRotation[0] = " << gRotation[0]  << endl;
                streamlog_out ( MESSAGE )  << " gRotation[1] = " << gRotation[1]  << endl;
                streamlog_out ( MESSAGE )  << " gRotation[2] = " << gRotation[2]  << endl;
          }
// new implmentation // Rubinskiy 11.11.11

// undo the shifts = go back to the center of the sensor frame (rotations unchanged)
      outputPosition[0] = inputPosition[0] - x_refhit;
      outputPosition[1] = inputPosition[1] - y_refhit;
      outputPosition[2] = inputPosition[2] - z_refhit;

      if ( _iEvt < _printEvents )
      {
         streamlog_out ( MESSAGE ) << "RevertGear: intermediate: "<<outputPosition[0] << " " <<outputPosition[1] << " " <<outputPosition[2] <<  endl;
      } 
      TVector3 iCenterOfSensorFrame( outputPosition[0], outputPosition[1], outputPosition[2] );

      iCenterOfSensorFrame.RotateZ( -gRotation[0] );
      iCenterOfSensorFrame.RotateY( -gRotation[1] );
      iCenterOfSensorFrame.RotateX( -gRotation[2] ); // first rotaiton in ZY plane -> around X axis (gamma)

      outputPosition[0] = iCenterOfSensorFrame(0);
      outputPosition[1] = iCenterOfSensorFrame(1);
      outputPosition[2] = iCenterOfSensorFrame(2);

         
          
/* old implementation
          // rotations first
 
          outputPosition[0] = inputPosition[0];
          outputPosition[1] = inputPosition[1];
          outputPosition[2] = inputPosition[2] - z_sensor;
          
          _EulerRotationInverse( sensorID,outputPosition, gRotation);

          // then the shifts
//          outputPosition[0] -= telPos[0];
//          outputPosition[1] -= telPos[1];
//          outputPosition[2] -= telPos[2];

          outputPosition[2] += z_sensor;
*/
      }
      else
      {
          // this hit belongs to a plane whose sensorID is not in the
          // alignment constants. So the idea is to eventually advice
          // the users if running in DEBUG and copy the not aligned hit
          // in the new collection.
          streamlog_out ( DEBUG ) << "Sensor ID " << sensorID << " not found. Skipping alignment for hit "
                                    << iHit << endl;

          for ( size_t i = 0; i < 3; ++i )
          {
              outputPosition[i] = inputPosition[i];
          }
          
      }

      if ( _iEvt < _printEvents )
      {
        streamlog_out ( MESSAGE ) << "RevertGear: INPUT: Sensor ID " << sensorID << " " << inputPosition[0] << " " << inputPosition[1] << " " << inputPosition[2] <<  endl;
        streamlog_out ( MESSAGE ) << "RevertGear: OUTPUT:Sensor ID " << sensorID << " " << outputPosition[0] << " " << outputPosition[1] << " " << outputPosition[2] << endl;
      }

      outputHit->setPosition( outputPosition ) ;
      outputCollectionVec->push_back( outputHit );
 
      TransformToLocalFrame(outputHit,evt); 
    }

    evt->addCollection( outputCollectionVec, _outputHitCollectionName );

  } catch (DataNotAvailableException& e) {
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }


}

void EUTelApplyAlignmentProcessor::Direct(LCEvent *event) {

  if ( _iEvt % 100 == 0 )
    streamlog_out ( MESSAGE4 ) << "Processing event  (ApplyAlignment Direct) "
                               << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                               << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                               << setfill(' ')
                               << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
//  ++_iEvt;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) 
  {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  }
  else if ( evt->getEventType() == kUNKNOWN ) 
  {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  try 
  {

    LCCollectionVec* inputCollectionVec         = dynamic_cast < LCCollectionVec* > (evt->getCollection(_inputHitCollectionName));
    LCCollectionVec* alignmentCollectionVec     = dynamic_cast < LCCollectionVec* > (evt->getCollection(_alignmentCollectionName));
    
    if (fevent) 
    {
//        LCCollection* alignmentCollectionVec     = static_cast < LCCollection* > (evt->getCollection(_alignmentCollectionName));          
        streamlog_out ( MESSAGE ) << "The alignment collection ["<< _alignmentCollectionName.c_str() <<"] contains: " <<  alignmentCollectionVec->getNumberOfElements() << " planes " << endl;    
    
        if(alignmentCollectionVec->size() > 0 )
        {
            streamlog_out ( MESSAGE ) << "alignment sensorID: " ;
 
            for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) 
            {
                EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
                _lookUpTable[ alignment->getSensorID() ] = iPos;
                streamlog_out ( MESSAGE ) << iPos << " " ;
                if ( _applyToReferenceHitCollection ) 
                {
                 AlignReferenceHit( evt,  alignment); 
                }
            }
            streamlog_out ( MESSAGE ) << endl;
            if ( _applyToReferenceHitCollection ) 
            {
            }
        }

        if ( _histogramSwitch )
        {          
            bookHistos();
        }
          
        streamlog_out ( MESSAGE ) << "The alignment collection ["<< _alignmentCollectionName.c_str() <<"] contains: " <<  alignmentCollectionVec->getNumberOfElements()
                                  << " planes " << endl;   

#ifndef NDEBUG
        // print out the lookup table
        map< int , int >::iterator mapIter = _lookUpTable.begin();
        while ( mapIter != _lookUpTable.end() ) 
        {
            streamlog_out ( DEBUG ) << "Sensor ID = " << mapIter->first
                                    << " is in position " << mapIter->second << endl;
            ++mapIter;
        }
#endif
    }


    LCCollectionVec * outputCollectionVec = 0;
    try{
        outputCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_outputHitCollectionName));
    }catch(...){
        outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);
    }

    for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {

      TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
//if(inputHit==0)continue;
      // now we have to understand which layer this hit belongs to.
      int sensorID = guessSensorID( inputHit );

      //find proper alignment colleciton:
            double alpha = 0.;
            double beta  = 0.;
            double gamma = 0.;
            double offsetX = 0.;
            double offsetY = 0.;
            double offsetZ = 0.;
 
      // now that we know at which sensor the hit belongs to, we can
      // get the corresponding alignment constants
      map< int , int >::iterator  positionIter = _lookUpTable.find( sensorID );

//      printf(" positionIter %5d at %5d \n", sensorID, positionIter->second );
      if ( positionIter == _lookUpTable.end() )
         {
//..           printf("wrong sensorId %5d ??\n", sensorID );
//           continue; do nothing as if alignment == 0.
         }
      else
         {
           EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * >  ( alignmentCollectionVec->getElementAt( positionIter->second ) );

           alpha = -alignment->getAlpha();
           beta  = -alignment->getBeta();
           gamma = -alignment->getGamma();
           offsetX = alignment->getXOffset();
           offsetY = alignment->getYOffset();
           offsetZ = alignment->getZOffset();
         }  



//        printf("alignment %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n" 
//                ,alignment->getXOffset()
//                ,alignment->getYOffset()
//                ,alignment->getZOffset()
//                ,alignment->getAlpha()
//                ,alignment->getBeta() 
//                ,alignment->getGamma()
//                );


//printf("%5d sensorID: %5d at %p \n", iHit, sensorID, inputHit );

/* obsolete, since we have refit collection now: Rubinskiy, 04.11.2011
      // determine z position of the plane
	  // 20 December 2010 @libov
      float	z_sensor = 0;
	  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) 
      {
          if (sensorID == _siPlanesLayerLayout->getID( iPlane ) ) 
          {
              z_sensor = _siPlanesLayerLayout -> getSensitivePositionZ( iPlane ) + 0.5 * _siPlanesLayerLayout->getSensitiveThickness( iPlane );
              break;
          }
      }
*/
      // refhit = center-of-the-sensor coordinates:
      double x_refhit = 0.; 
      double y_refhit = 0.; 
      double z_refhit = 0.; 
      
      if( _referenceHitVec == 0 )
      {
        streamlog_out(MESSAGE) << "_referenceHitVec is empty" << endl;
      }
      else
      {
//        streamlog_out(MESSAGE) << "reference Hit collection name : " << _referenceHitCollectionName << endl;
 
       for(size_t ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
       {
        EUTelReferenceHit * refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
        if( sensorID != refhit->getSensorID() )
        {
       //   streamlog_out(MESSAGE) << "Looping through a varity of sensor IDs" << endl;
          continue;
        }
        else
        {
//          streamlog_out(MESSAGE) << "Sensor ID and Alignment plane ID match!" << endl;
          x_refhit =  refhit->getXOffset();
          y_refhit =  refhit->getYOffset();
          z_refhit =  refhit->getZOffset();

          x_refhit += offsetX;
          y_refhit += offsetY;
          z_refhit += offsetZ;

          break;
        } 
/*        printf("PRE sensorID: %5d dx:%5.3f dy:%5.3f dz:%5.3f  [%5.2f %5.2f %5.2f]\n",
        refhit->getSensorID(   ),                      
        refhit->getXOffset(    ),
        refhit->getYOffset(    ),
        refhit->getZOffset(    ), x_refhit, y_refhit, z_refhit
        );
*/
       }
      }
      // copy the input to the output, at least for the common part
      TrackerHitImpl   * outputHit  = new TrackerHitImpl;
      outputHit->setType( inputHit->getType() );
      outputHit->rawHits() = inputHit->getRawHits();

      // hit coordinates in the center-of-the sensor frame (axis coincide with the global frame)
      double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
      inputPosition[0] = inputPosition[0] - x_refhit;
      inputPosition[1] = inputPosition[1] - y_refhit;
      inputPosition[2] = inputPosition[2] - z_refhit;

      // initial setting of the sensor = center-of-the sensor coordinates in the global frame
      double   outputPosition[3]  = { 0., 0., 0. };
      outputPosition[0] = x_refhit;
      outputPosition[1] = y_refhit;
      outputPosition[2] = z_refhit;
/*
        printf("PRE sensorID: %5d   [%5.2f %5.2f %5.2f]\n",
                sensorID, x_refhit, y_refhit, z_refhit
        );
*/


//      if ( positionIter != _lookUpTable.end() ) {


#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )
        string tempHistoName;
        AIDA::IHistogram3D *histo3D; 
        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss  << _hitHistoBeforeAlignName << "_" << sensorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( inputPosition[0], inputPosition[1] );
          }
          else
            {
              streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                        << ".\nDisabling histogramming from now on " << endl;
              _histogramSwitch = false;
            }
          histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotBeforeAlignName ] );
          if ( histo3D ) histo3D->fill( inputPosition[0], inputPosition[1], inputPosition[2] );
          else {
            streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << _densityPlotBeforeAlignName
                                      << ".\nDisabling histogramming from now on " << endl;
            _histogramSwitch = false;
          }
        }
#endif

        
        if ( _correctionMethod == 0 ) 
        {

            // this is the shift only case
            outputPosition[0] += inputPosition[0] - offsetX;                
            outputPosition[1] += inputPosition[1] - offsetY; 
            outputPosition[2] += inputPosition[2] - offsetZ;

        } 
        else if ( _correctionMethod == 1 ) 
        {

/*            double alpha = 0.;
            double beta  = 0.;
            double gamma = 0.;
            double offsetX = 0.;
            double offsetY = 0.;
            double offsetZ = 0.;
*/
            if ( _debugSwitch )
            {
                alpha = _alpha;
                beta  = _beta;
                gamma = _gamma; 
                offsetX = 0.;
                offsetY = 0.;
                offsetZ = 0.;
            }
            else
            {
  //              alpha = alignment->getAlpha();
  //              beta  = alignment->getBeta();
  //              gamma = alignment->getGamma();
  //              offsetX = alignment->getXOffset();
  //              offsetY = alignment->getYOffset();
  //              offsetZ = alignment->getZOffset();
            }
            if( _iEvt < _printEvents )
            {
                if ( _debugSwitch ) 
                {
                    streamlog_out ( MESSAGE )  << "Debugmode ON " << endl;                                   
                }
                
                streamlog_out ( MESSAGE )  << "_correctionMethod == rotation first " << endl;
                streamlog_out ( MESSAGE )  << " alignment->getAlpha() = " << alpha  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getBeta()  = " <<  beta  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getGamma() = " << gamma << endl;
                streamlog_out ( MESSAGE )  << " alignment->getXOffest() = " << offsetX  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getYOffest() = " << offsetY  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getZOffest() = " << offsetZ  << endl;
            }
        
            // this is the rotation first
            TVector3 iCenterOfSensorFrame(  inputPosition[0],  inputPosition[1],  inputPosition[2] );

            iCenterOfSensorFrame.RotateX( alpha );
            iCenterOfSensorFrame.RotateY( beta  );
            iCenterOfSensorFrame.RotateZ( gamma );

            outputPosition[0] += iCenterOfSensorFrame(0);
            outputPosition[1] += iCenterOfSensorFrame(1);
            outputPosition[2] += iCenterOfSensorFrame(2);

 /* orig
            // first the rotation (matrix layout)
            outputPosition[0] +=                inputPosition[0] + (-1)*gamma * inputPosition[1] +      beta  * inputPosition[2] ;
            outputPosition[1] +=        gamma * inputPosition[0] +              inputPosition[1] + (-1)*alpha * inputPosition[2] ;
            outputPosition[2] += (-1) * beta  * inputPosition[0] +      alpha * inputPosition[1] +              inputPosition[2] ;
  */

//	x = x_temp * (1 + alpha * alpha ) + ( (-1) * gamma - alpha * beta) * y_temp + ( (-1) * beta + alpha * gamma) * z_temp;
//	y = x_temp * (gamma - alpha * beta ) + (1 + beta * beta) * y_temp + ((-1) * alpha - beta * gamma) * z_temp;
//	z = x_temp * (beta + alpha * gamma ) + (alpha - gamma * beta) * y_temp + ( 1 + gamma * gamma) * z_temp;
//            outputPosition[0] +=  (1+alpha*alpha)*inputPosition[0] + (-1)*gamma * inputPosition[1] +      beta  * inputPosition[2] ;
//            outputPosition[1] +=        gamma * inputPosition[0] +              inputPosition[1] + (-1)*alpha * inputPosition[2] ;
//            outputPosition[2] += (-1) * beta  * inputPosition[0] +      alpha * inputPosition[1] +              inputPosition[2] ;



            // second the shift
            outputPosition[0] -= offsetX; 
            outputPosition[1] -= offsetY; 
            outputPosition[2] -= offsetZ; 
//            printf(" final %5d  %5.2f %5.2f %5.2f :: %5.2f %5.2f %5.2f \n", sensorID, outputPosition[0], outputPosition[1], outputPosition[2],
//                                  offsetX, offsetY, offsetZ );
 
//            outputPosition[2] += z_sensor ;
          
        }
        else if ( _correctionMethod == 2 ) 
        {
/*            double alpha = 0.;
            double beta  = 0.;
            double gamma = 0.;
            double offsetX = 0.;
            double offsetY = 0.;
            double offsetZ = 0.;
*/
            if ( _debugSwitch )
            {
                alpha = _alpha;
                beta  = _beta;
                gamma = _gamma; 
                offsetX = 0.;
                offsetY = 0.;
                offsetZ = 0.;
            }
            else
            {
   //             alpha = alignment->getAlpha();
   //             beta  = alignment->getBeta();
   //             gamma = alignment->getGamma();
   //             offsetX = alignment->getXOffset();
   //             offsetY = alignment->getYOffset();
   //             offsetZ = alignment->getZOffset();
            }
            if( _iEvt < _printEvents )
            {
                if ( _debugSwitch ) 
                {
                    streamlog_out ( MESSAGE )  << "Debugmode ON " << endl;                                   
                }
                
                streamlog_out ( MESSAGE )  << "_correctionMethod == rotation first " << endl;
                streamlog_out ( MESSAGE )  << " alignment->getAlpha() = " << alpha  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getBeta()  = " <<  beta  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getGamma() = " << gamma << endl;
                streamlog_out ( MESSAGE )  << " alignment->getXOffest() = " << offsetX  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getYOffest() = " << offsetY  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getZOffest() = " << offsetZ  << endl;
            }

          // this is the translation first
          // first the shifts
          inputPosition[0] -= offsetX;
          inputPosition[1] -= offsetY; 
          inputPosition[2] -= offsetZ;

          // second the rotation (matrix layout)
          outputPosition[0] +=                inputPosition[0] + (-1)* gamma * inputPosition[1] +      beta  * inputPosition[2]  ;
          outputPosition[1] +=        gamma * inputPosition[0] +               inputPosition[1] + (-1)*alpha * inputPosition[2]  ;
          outputPosition[2] += (-1) * beta  * inputPosition[0] +       alpha * inputPosition[1] +              inputPosition[2]  ;

//          outputPosition[2] += z_sensor;
        }

#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) ) 
        if ( _histogramSwitch ) {
          {
            stringstream ss;
            ss  << _hitHistoAfterAlignName << "_" << sensorID ;
            tempHistoName = ss.str();
          }
          if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) {
            histo->fill( outputPosition[0], outputPosition[1] );
          }
          else {
            streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                      << ".\nDisabling histogramming from now on " << endl;
            _histogramSwitch = false;
          }
          histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotAfterAlignName ] );
          if ( histo3D ) histo3D->fill( outputPosition[0], outputPosition[1], outputPosition[2] );
          else {
            streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << _densityPlotAfterAlignName
                                      << ".\nDisabling histogramming from now on " << endl;
            _histogramSwitch = false;
          }
        }
#endif

/*
      } 
      else 
      {

        // this hit belongs to a plane whose sensorID is not in the
        // alignment constants. So the idea is to eventually advice
        // the users if running in DEBUG and copy the not aligned hit
        // in the new collection.
//        streamlog_out ( DEBUG ) << "Sensor ID " << sensorID << " not found. Skipping alignment for hit "
//                                << iHit << endl;

        for ( size_t i = 0; i < 3; ++i )
        {
            outputPosition[i] = inputPosition[i];
        }

      }
*/

      if ( _iEvt < _printEvents )
      {
         streamlog_out ( MESSAGE ) << "DIRECT: INPUT: Sensor ID " << sensorID << " " << inputPosition[0] << " " << inputPosition[1] << " " << inputPosition[2] <<  endl;   
         streamlog_out ( MESSAGE ) << "DIRECT: OUTPUT:Sensor ID " << sensorID << " " << outputPosition[0] << " " << outputPosition[1] << " " << outputPosition[2]  << endl;
      }

      outputHit->setPosition( outputPosition ) ;
      outputCollectionVec->push_back( outputHit );
    }

    try{
       evt->addCollection( outputCollectionVec, _outputHitCollectionName );
    }catch(...){
//       printf("can not add Collection %s \n", _outputHitCollectionName.c_str());
    }

  } catch (DataNotAvailableException& e) {
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

}

void EUTelApplyAlignmentProcessor::Reverse(LCEvent *event) {

    if ( _iEvt % 100 == 0 )
        streamlog_out ( MESSAGE4 ) << "Processing event (ApplyAlignment Reverse) "
                               << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                               << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                               << setfill(' ')
                               << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;
//    ++_iEvt;


    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

    if ( evt->getEventType() == kEORE ) 
    {
        streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
        return;
    }
    else if ( evt->getEventType() == kUNKNOWN ) 
    {
        streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
            << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }


    try 
    {
        LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
        LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));
    
        if (fevent) 
        {
            LCCollectionVec * alignmentCollectionVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));          
            streamlog_out ( MESSAGE ) << "The alignment collection ["<< _alignmentCollectionName.c_str() <<"] contains: " <<  alignmentCollectionVec->size() << " planes " << endl;    
    
            if(alignmentCollectionVec->size() > 0 )
            {
                streamlog_out ( MESSAGE ) << "alignment sensorID: " ;
                for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) 
                {
                    EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
                    _lookUpTable[ alignment->getSensorID() ] = iPos;
                    streamlog_out ( MESSAGE ) << iPos << " " ;
                }
                streamlog_out ( MESSAGE ) << endl;
            }
            
            if ( _histogramSwitch ) 
            {
//                crashes if commented out , echeck why ??
//                bookHistos();
            }
            
            streamlog_out ( MESSAGE ) << "The alignment collection ["<< _alignmentCollectionName.c_str() <<"] contains: " <<  alignmentCollectionVec->size()
                                      << " planes " << endl;   

#ifndef NDEBUG
            // print out the lookup table
            map< int , int >::iterator mapIter = _lookUpTable.begin();
            while ( mapIter != _lookUpTable.end() ) 
            {
                streamlog_out ( DEBUG ) << "Sensor ID = " << mapIter->first
                                        << " is in position " << mapIter->second << endl;
                ++mapIter;
            }
#endif
        }


        LCCollectionVec * outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);

        for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) 
        {

            TrackerHitImpl   * inputHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;

            // now we have to understand which layer this hit belongs to.
            int sensorID = guessSensorID( inputHit );

            // determine z position of the plane
            // 20 December 2010 @libov

            float	z_sensor = 0;
            for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) 
            {
                if (sensorID == _siPlanesLayerLayout->getID( iPlane ) ) 
                {
                    z_sensor = _siPlanesLayerLayout -> getSensitivePositionZ( iPlane ) + 0.5 * _siPlanesLayerLayout->getSensitiveThickness( iPlane );
                    break;
                }
            }

            // copy the input to the output, at least for the common part
            TrackerHitImpl   * outputHit  = new TrackerHitImpl;
            outputHit->setType( inputHit->getType() );
            outputHit->rawHits() = inputHit->getRawHits();

            // now that we know at which sensor the hit belongs to, we can
            // get the corresponding alignment constants
            map< int , int >::iterator  positionIter = _lookUpTable.find( sensorID );

            double * inputPosition      = const_cast< double * > ( inputHit->getPosition() ) ;
            double   outputPosition[3]  = { 0., 0., 0. };

            if ( positionIter != _lookUpTable.end() ) 
            {

#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )
                string tempHistoName;
                AIDA::IHistogram3D *histo3D; 
                if ( _histogramSwitch ) 
                {
                    {
                        stringstream ss;
                        ss  << _hitHistoBeforeAlignName << "_" << sensorID ;
                        tempHistoName = ss.str();
                    }
                    if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) 
                    {
                        histo->fill( inputPosition[0], inputPosition[1] );
                    }
                    else
                    {
                        streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                                  << ".\nDisabling histogramming from now on " << endl;
                        _histogramSwitch = false;
                    }
                    
                    histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotBeforeAlignName ] );
                    if ( histo3D ) 
                    {
                        histo3D->fill( inputPosition[0], inputPosition[1], inputPosition[2] );
                    }
                    else 
                    {
                        streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << _densityPlotBeforeAlignName
                                                  << ".\nDisabling histogramming from now on " << endl;
                        _histogramSwitch = false;
                    }
                }
#endif

                EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( positionIter->second ) );

//        printf("alignment %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n" 
//                ,alignment->getXOffset()
//                ,alignment->getYOffset()
//                ,alignment->getZOffset()
//                ,alignment->getAlpha()
//                ,alignment->getBeta() 
//                ,alignment->getGamma()
//                );
        
                if ( _correctionMethod == 0 ) 
                {                   

                    // this is the shift only case

                    outputPosition[0] = inputPosition[0] + alignment->getXOffset();
                    outputPosition[1] = inputPosition[1] + alignment->getYOffset();
                    outputPosition[2] = inputPosition[2] + alignment->getZOffset();

                }
                else
                    if ( _correctionMethod == 1 ) 
                    {
                        // this is the rotation first
 
                        double alpha = 0.;
                        double beta  = 0.;
                        double gamma = 0.;
            double offsetX = 0.;
            double offsetY = 0.;
            double offsetZ = 0.;

            if ( _debugSwitch )
            {
                alpha = _alpha;
                beta  = _beta;
                gamma = _gamma; 
                offsetX = 0.;
                offsetY = 0.;
                offsetZ = 0.;
            }
            else
            {
                alpha = alignment->getAlpha();
                beta  = alignment->getBeta();
                gamma = alignment->getGamma();
                offsetX = alignment->getXOffset();
                offsetY = alignment->getYOffset();
                offsetZ = alignment->getZOffset();
            }
            if( _iEvt < _printEvents )
            {
                if ( _debugSwitch ) 
                {
                    streamlog_out ( MESSAGE )  << "Debugmode ON " << endl;                                   
                }
                
                streamlog_out ( MESSAGE )  << "_correctionMethod == rotation first " << endl;
                streamlog_out ( MESSAGE )  << " alignment->getAlpha() = " << alpha  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getBeta()  = " <<  beta  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getGamma() = " << gamma << endl;
                streamlog_out ( MESSAGE )  << " alignment->getXOffest() = " << offsetX  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getYOffest() = " << offsetY  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getZOffest() = " << offsetZ  << endl;
            }


       
                        // second the rotation
                        // for the inverse matrix derivation see paper log book 19/01/2011
                        // libov@mail.desy.de

                        // local variables (X,Y,Z) position of a hit
                        // x_temp (original) -> x (rotated)
                        // (x_temp,y_temp,z_temp) -> (x,y,z) 
                        // rotation matrix is:
                        //      1 + a^2    -g - a*b  -b + a*g
                        //      g - a*b     1 + b^2  -a - b*g 
                        //      b + a*g     a - g*b   1 + g^2  
                        // 
                        double	x_temp = inputPosition[0];
                        double	y_temp = inputPosition[1];
                        double	z_temp = inputPosition[2] - z_sensor;
// no correct any more due to changes in convention of the angle signs!
// further more a more clear way would be to multiply an input vector by an inverse rotation matrix  
// Igor Rubinsky 09-10-2011

                        // rotation first
                        double x = x_temp * (1 + alpha * alpha )    + ( (-1) * gamma - alpha * beta) * y_temp    + ( (-1) * beta + alpha * gamma) * z_temp;
                        double y = x_temp * (gamma - alpha * beta ) + (1 + beta * beta) * y_temp                 + ((-1) * alpha - beta * gamma) * z_temp;
                        double z = x_temp * (beta + alpha * gamma ) + (alpha - gamma * beta) * y_temp            + ( 1 + gamma * gamma) * z_temp;

                        double det = 1 + alpha * alpha + beta * beta + gamma * gamma;

                        x = x / det;
                        y = y / det;
                        z = z / det;

                        x += alignment->getXOffset();
                        y += alignment->getYOffset();
                        z += alignment->getZOffset();

                        // now final coordinates
                        outputPosition[0] = x ;
                        outputPosition[1] = y;
                        outputPosition[2] = z + z_sensor ;
          
                    }
                    else
                        if ( _correctionMethod == 2 ) 
                        {
                            // this is the translation first

                            double alpha = 0.;
                            double beta  = 0.;
                            double gamma = 0.;
            double offsetX = 0.;
            double offsetY = 0.;
            double offsetZ = 0.;

            if ( _debugSwitch )
            {
                alpha = _alpha;
                beta  = _beta;
                gamma = _gamma; 
                offsetX = 0.;
                offsetY = 0.;
                offsetZ = 0.;
            }
            else
            {
                alpha = alignment->getAlpha();
                beta  = alignment->getBeta();
                gamma = alignment->getGamma();
                offsetX = alignment->getXOffset();
                offsetY = alignment->getYOffset();
                offsetZ = alignment->getZOffset();
            }
            if( _iEvt < _printEvents )
            {
                if ( _debugSwitch ) 
                {
                    streamlog_out ( MESSAGE )  << "Debugmode ON " << endl;                                   
                }
                
                streamlog_out ( MESSAGE )  << "_correctionMethod == rotation first " << endl;
                streamlog_out ( MESSAGE )  << " alignment->getAlpha() = " << alpha  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getBeta()  = " <<  beta  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getGamma() = " << gamma << endl;
                streamlog_out ( MESSAGE )  << " alignment->getXOffest() = " << offsetX  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getYOffest() = " << offsetY  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getZOffest() = " << offsetZ  << endl;
            }


                            if ( _debugSwitch )
                            {
                               
                                alpha = _alpha;
                                beta  = _beta;
                                gamma = _gamma;
                            }
                            else
                            {
                                alpha = alignment->getAlpha();
                                beta  = alignment->getBeta();
                                gamma = alignment->getGamma();
                            }
                            if( _iEvt < _printEvents )
                            {
                                if ( _debugSwitch ) 
                                {
                                    streamlog_out ( MESSAGE )  << "Debugmode ON " << endl;                                   
                                }
                                streamlog_out ( MESSAGE )  << "_correctionMethod == translation  first " << endl;
                                streamlog_out ( MESSAGE )  << " alignment->getAlpha() = " << alpha  << endl;
                                streamlog_out ( MESSAGE )  << " alignment->getBeta() = " <<  beta  << endl;
                                streamlog_out ( MESSAGE )  << " alignment->getGamma() = " << gamma << endl;
                                streamlog_out ( MESSAGE )  << " alignment->getXOffest() = " << alignment->getXOffset() << endl;
                                streamlog_out ( MESSAGE )  << " alignment->getYOffest() = " << alignment->getYOffset() << endl;
                                streamlog_out ( MESSAGE )  << " alignment->getZOffest() = " << alignment->getZOffset() << endl;
                            }
 
                            
                            double	x_temp = inputPosition[0];
                            double	y_temp = inputPosition[1];
                            double	z_temp = inputPosition[2] - z_sensor;

                            // first the shifts
                            x_temp += alignment->getXOffset();
                            y_temp += alignment->getYOffset();
                            z_temp += alignment->getZOffset();


                            // second the rotation (matrix layout)
// no correct! see comment above
// Igor Rubinsky 09-10-2011
                            double x = x_temp * (1 + alpha * alpha ) + ( (-1) * gamma - alpha * beta) * y_temp + ( (-1) * beta + alpha * gamma) * z_temp;
                            double y = x_temp * (gamma - alpha * beta ) + (1 + beta * beta) * y_temp + ((-1) * alpha - beta * gamma) * z_temp;
                            double z = x_temp * (beta + alpha * gamma ) + (alpha - gamma * beta) * y_temp + ( 1 + gamma * gamma) * z_temp;

                            double det = 1 + alpha * alpha + beta * beta + gamma * gamma;

                            x = x / det;
                            y = y / det;
                            z = z / det;


                            // now final coordinates
                            outputPosition[0] = x ;
                            outputPosition[1] = y;
                            outputPosition[2] = z + z_sensor;
                        }

#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) ) 

                if ( _histogramSwitch ) 
                {
                    {
                        stringstream ss;
                        ss  << _hitHistoAfterAlignName << "_" << sensorID ;
                        tempHistoName = ss.str();
                    }
                    if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*> ( _aidaHistoMap[ tempHistoName ] )) 
                    {
                        histo->fill( outputPosition[0], outputPosition[1] );
                    }
                    else 
                    {
                        streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                                  << ".\nDisabling histogramming from now on " << endl;
                        _histogramSwitch = false;
                    }
                    
                    histo3D = dynamic_cast<AIDA::IHistogram3D*> (_aidaHistoMap[ _densityPlotAfterAlignName ] );
                    if ( histo3D )
                    {
                        histo3D->fill( outputPosition[0], outputPosition[1], outputPosition[2] );
                    }
                    else 
                    {
                        streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << _densityPlotAfterAlignName
                                                  << ".\nDisabling histogramming from now on " << endl;
                        _histogramSwitch = false;
                    }
                }
#endif

            }
            else 
            {

                // this hit belongs to a plane whose sensorID is not in the
                // alignment constants. So the idea is to eventually advice
                // the users if running in DEBUG and copy the not aligned hit
                // in the new collection.

                streamlog_out ( DEBUG   ) << "Sensor ID " << sensorID << " not found. Skipping alignment for hit "
                                          << iHit << endl;
                
                for ( size_t i = 0; i < 3; ++i )
                {
                    outputPosition[i] = inputPosition[i];
                }
            }

            if ( _iEvt < _printEvents )
            {
               streamlog_out ( MESSAGE ) << "Reverse: INPUT: Sensor ID " << sensorID << " " << inputPosition[0] << " " << inputPosition[1] << " " << inputPosition[2] << " " << endl;                
               streamlog_out ( MESSAGE ) << "Reverse: OUTPUT:Sensor ID " << sensorID << " " << outputPosition[0] << " " << outputPosition[1] << " " << outputPosition[2] << " " << endl;                
            }
            outputHit->setPosition( outputPosition ) ;
            outputCollectionVec->push_back( outputHit );
        }

        evt->addCollection( outputCollectionVec, _outputHitCollectionName );

    }
    catch (DataNotAvailableException& e) 
    {
        streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                    << " in run " << event->getRunNumber() << endl;
    }
    
}

void EUTelApplyAlignmentProcessor::bookHistos() {

#if ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

  try {
    streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;

    string tempHistoName;

    // histograms are grouped into folders named after the
    // detector. This requires to loop on detector now.
    for (int iDet = 0 ; iDet < _siPlanesParameters->getSiPlanesNumber(); iDet++) 
    {
      int sensorID = _siPlanesLayerLayout->getID( iDet ) ;
 
      streamlog_out ( MESSAGE4 ) <<  "Booking histograms for sensorID: " << sensorID << endl;

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
      double safetyFactor = 1.0;

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
    double safetyFactor = 1.0;

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

    streamlog_out ( MESSAGE4 ) <<  "Booking histograms DONE" << endl;

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


void EUTelApplyAlignmentProcessor::check (LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelApplyAlignmentProcessor::end() {
  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;
 

  delete [] _siPlaneZPosition;

}

int EUTelApplyAlignmentProcessor::guessSensorID( TrackerHitImpl * hit ) {

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;
  double * hitPosition = const_cast<double * > (hit->getPosition());

//  LCCollectionVec * referenceHitVec     = dynamic_cast < LCCollectionVec * > (evt->getCollection( _referenceHitCollectionName));
  if( _referenceHitVec == 0)
  {
    streamlog_out(MESSAGE) << "_referenceHitVec is empty" << endl;
    return 0;
  }

        try
        {
            LCObjectVec clusterVector = hit->getRawHits();

            EUTelVirtualCluster *cluster = 0;
// printf("cluster type %5d \n", hit->getType() );
            if ( hit->getType() == kEUTelBrickedClusterImpl ) {

               // fixed cluster implementation. Remember it
               //  can come from
               //  both RAW and ZS data
   
               cluster = new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
                
            } else if ( hit->getType() == kEUTelDFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } else if ( hit->getType() == kEUTelFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } 
            else if ( hit->getType() == kEUTelAPIXClusterImpl ) 
            {
              
//              cluster = new EUTelSparseClusterImpl< EUTelAPIXSparsePixel >
//                 ( static_cast<TrackerDataImpl *> ( clusterVector[ 0 ]  ) );

                // streamlog_out(MESSAGE4) << "Type is kEUTelAPIXClusterImpl" << endl;
                TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );
                // streamlog_out(MESSAGE4) << "Vector size is: " << clusterVector.size() << endl;

                cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);
	      
	        // CellIDDecoder<TrackerDataImpl> cellDecoder(clusterFrame);
                eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);
                
            }
            else if ( hit->getType() == kEUTelSparseClusterImpl ) 
            {
               cluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > ( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            }

            if(cluster != 0)
            {
            int sensorID = cluster->getDetectorID();
//            printf("ApplyAlignement :: actual sensorID: %5d \n", sensorID );
            return sensorID;
            }
          }
          catch(...)
          {
            printf("guess SensorID crashed \n");
          }


/*
      for(size_t ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
      {
        EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
//        printf(" _referenceHitVec %p refhit %p \n", _referenceHitVec, refhit);
        
        TVector3 hit3d( hitPosition[0], hitPosition[1], hitPosition[2] );
        TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
 
        double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
//        printf("iPlane %5d   hitPos:  [%8.3f;%8.3f%8.3f]  distance: %8.3f \n", refhit->getSensorID(), hitPosition[0],hitPosition[1],hitPosition[2], distance  );
        if ( distance < minDistance ) 
        {
           minDistance = distance;
           sensorID = refhit->getSensorID();
//           printf("sensorID: %5d \n", sensorID );
        }    

      }
 */
/*
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) 
  {
//      printf("iPlane %5d   hitPos:  %8.3f  siZpos: %8.3f \n", iPlane, hitPosition[2] , _siPlaneZPosition[ iPlane ] );
      double distance = std::abs( hitPosition[2] - _siPlaneZPosition[ iPlane ] );
      if ( distance < minDistance ) 
      {
          minDistance = distance;
          sensorID = _siPlanesLayerLayout->getID( iPlane );
      }
  }
  if  ( _siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT ) 
  {
      double distance = std::abs( hitPosition[2] - _siPlanesLayerLayout->getDUTPositionZ() );
      if( distance < minDistance )
      {
        minDistance = distance;
        sensorID = _siPlanesLayerLayout->getDUTID();
      }
  }
  if ( minDistance > 10  ) //mm 
  {
    // advice the user that the guessing wasn't successful 
    streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
      "Please check the consistency of the data with the GEAR file " << endl;
    throw SkipEventException(this);
  }
*/
  return sensorID;
}


void EUTelApplyAlignmentProcessor::TransformToLocalFrame(TrackerHitImpl* outputHit, LCEvent * ev) 
{
// revert alignment, in an inverse order... this part should be fixed // perhaps it's obsolete already// Rubinskiy 12.11.2011
//	for ( int i = _alignmentCollectionNames.size() - 1; i >= 0; i--) 
//      {
//          revertAlignment (x, y, z, _alignmentCollectionNames[i], ev );
//	}

        double *outputPosition = const_cast< double * > ( outputHit->getPosition() ) ;

        // now we have to understand which layer this hit belongs to.
        int sensorID = guessSensorID( outputHit );

        if ( _conversionIdMap.size() != (unsigned) _siPlanesParameters->getSiPlanesNumber() ) 
        {
          // first of all try to see if this sensorID already belong to
          if ( _conversionIdMap.find( sensorID ) == _conversionIdMap.end() ) 
          {
              // this means that this detector ID was not already inserted,
              // so this is the right place to do that
          
              for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) 
              {
                  if ( _siPlanesLayerLayout->getID(iLayer) == sensorID ) 
                  {
                      _conversionIdMap.insert( make_pair( sensorID, iLayer ) );
                      break;
                  }
              }
          }
        }
      
        int layerIndex   = _conversionIdMap[sensorID];


        double xSize = 0., ySize = 0.;
        double zThickness = 0.;
        double xPitch = 0., yPitch = 0.;
        double xPointing[2] = { 1., 0. }, yPointing[2] = { 1., 0. };

        xPitch       = _siPlanesLayerLayout->getSensitivePitchX(layerIndex);    // mm
        yPitch       = _siPlanesLayerLayout->getSensitivePitchY(layerIndex);    // mm
        xSize        = _siPlanesLayerLayout->getSensitiveSizeX(layerIndex);     // mm
        ySize        = _siPlanesLayerLayout->getSensitiveSizeY(layerIndex);     // mm
        xPointing[0] = _siPlanesLayerLayout->getSensitiveRotation1(layerIndex); // was -1 ;
        xPointing[1] = _siPlanesLayerLayout->getSensitiveRotation2(layerIndex); // was  0 ;
        yPointing[0] = _siPlanesLayerLayout->getSensitiveRotation3(layerIndex); // was  0 ;
        yPointing[1] = _siPlanesLayerLayout->getSensitiveRotation4(layerIndex); // was -1 ;

        double sign = 0;
          if      ( xPointing[0] < -0.7 )     sign = -1 ;
          else if ( xPointing[0] > 0.7 )      sign =  1 ;
          else 
          {
           if       ( xPointing[1] < -0.7 )   sign = -1 ;
           else if  ( xPointing[1] > 0.7 )    sign =  1 ;
          }
          outputPosition[0] -=  (-1)* sign * xSize/2;   //apply shifts few lines later

          if      ( yPointing[0] < -0.7 )     sign = -1 ;
          else if ( yPointing[0] > 0.7 )      sign =  1 ;
          else 
          {
           if       ( yPointing[1] < -0.7 )   sign = -1 ;
           else if  ( yPointing[1] > 0.7 )    sign =  1 ;
          }
          outputPosition[1] -= (-1)* sign * ySize/2;   // apply shifts few lines later

        TMatrix flip0(2,1);
          flip0(0,0) = outputPosition[0];
          flip0(1,0) = outputPosition[1];

        TMatrix flip1(2,1);
          flip1(0,0) = 0.;
          flip1(1,0) = 0.;

        TMatrix flip(2,2);
          flip(0,0) = xPointing[0];  
          flip(1,0) = yPointing[0];  
          flip(0,1) = xPointing[1];  
          flip(1,1) = yPointing[1];  

        TMatrix antiflip(2,2);
          antiflip = flip.Invert();
         
        flip1.Mult(antiflip, flip0);

        TVector2 clusterCenter( flip1(0,0)/xPitch-0.5, flip1(1,0)/yPitch-0.5 ); 
        if ( _iEvt < _printEvents )
        {
          streamlog_out ( MESSAGE ) << "RevertGear: matrix: " << flip1(0,0) << ":" << flip1(1,0) << endl;
          streamlog_out ( MESSAGE ) << "RevertGear: cluster coordinates [pitch:"<<xPitch<<":"<<yPitch<<"]:" << clusterCenter.X() << " " << clusterCenter.Y() << endl;
        }

        EUTelVirtualCluster  * cluster;
        ClusterType type = static_cast<ClusterType>(static_cast<int>(outputHit->getType()));
        if ( _iEvt < _printEvents )
        {
          streamlog_out ( MESSAGE ) << "RevertGear: hit/cluster type:" << type << endl;
        }

        // prepare a LCObjectVec to store the current cluster
        LCObjectVec clusterVec = outputHit->getRawHits();
        for(int i=0;i<clusterVec.size();i++)
        {      
          if ( _iEvt < _printEvents )
          {
            streamlog_out ( MESSAGE ) << "[" << i << "]" ;
          }
          if ( type == kEUTelSparseClusterImpl ) 
          {
            cluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > ( static_cast<TrackerDataImpl *> ( clusterVec[i]  ) );
          } 
          else if ( type == kEUTelAPIXClusterImpl ) 
          {
        
            cluster = new EUTelSparseClusterImpl< EUTelAPIXSparsePixel > ( static_cast<TrackerDataImpl *> ( clusterVec[i]  ) );
//            printf("ID: %5d:  x:%5d y:%5d tot:%5.1f\n", cluster->getDetectorID(), xsize0, ysize0,  cluster->getTotalCharge() );
          }
          else 
          {
            streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
            throw UnknownDataTypeException("Cluster type unknown");
          }
          float xCoG(0.0f), yCoG(0.0f);
          cluster->getCenterOfGravity(xCoG, yCoG);
          double xDet = (xCoG + 0.5) * xPitch;
          double yDet = (yCoG + 0.5) * yPitch;

          double telPos[3];
          telPos[0] = xPointing[0] * xDet + xPointing[1] * yDet;
          telPos[1] = yPointing[0] * xDet + yPointing[1] * yDet;

          // now the translation
          // not sure about the sign. At least it is working for the current
          // configuration but we need to double check it
          double sign = 0;
          if      ( xPointing[0] < -0.7 )       sign = -1 ;
          else if ( xPointing[0] > 0.7 )       sign =  1 ;
          else 
          {
           if       ( xPointing[1] < -0.7 )    sign = -1 ;
           else if  ( xPointing[1] > 0.7 )    sign =  1 ;
          }
          telPos[0] +=  (-1)* sign * xSize/2;   //apply shifts few lines later

          if      ( yPointing[0] < -0.7 )       sign = -1 ;
          else if ( yPointing[0] > 0.7 )       sign =  1 ;
          else 
          {
           if       ( yPointing[1] < -0.7 )    sign = -1 ;
           else if  ( yPointing[1] > 0.7 )    sign =  1 ;
          }
          telPos[1] += (-1)* sign * ySize/2;   // apply shifts few lines later

          if ( _iEvt < _printEvents )
          {
            streamlog_out ( MESSAGE ) << "[CoG:" << xCoG << ":" << yCoG << " -> Det:" << xDet << ":" << yDet << "] == " << telPos[0] << ":" << telPos[1] ;
          }
        }
 
        if ( _iEvt < _printEvents )
        { 
          streamlog_out ( MESSAGE ) << endl;
        }

}

void EUTelApplyAlignmentProcessor::revertAlignment(double & x, double & y, double & z, std::string	collectionName, LCEvent * lcevent) 
{

	// in this function, some parts of the EUTelApplyAlignmentProcessor are used

	// get the alignment constant object
	// first, get the alignment collection
	LCCollectionVec * alignmentCollectionVec = dynamic_cast < LCCollectionVec * > ( lcevent->getCollection(collectionName) );

    // next, find the alignment constant corresponding to the DUT
	EUTelAlignmentConstant * c = NULL;
    
    int _manualDUTid = 10;

	for ( size_t iPos = 0; iPos < alignmentCollectionVec->size(); ++iPos ) {

		c = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iPos ) );
		if (c -> getSensorID() == _manualDUTid ) break;	// this means we found the alignment constant corresponding
    													// to the DUT; the pointer to it is now stored in c and can
														// be furhter used
	}

	if ( c == NULL ) {
		cout << "Was not possible to found alignment constant, terminating" << endl;
		abort();
	}

	// now apply the constants

	// the way we apply constants is correctionMethod1 (first the rotations, second the shifts)
	// to revert the alignment transformation properly, the constants have to be applied in
	// a reverted way, i.e. first the shifts, second the rotations

	// not that the sign is different for all the constants wrt to what is done in EUTelApplyAlignmentProcessor -
	// the transformation is reverted

	// first the shifts
	x += c->getXOffset();
	y += c->getYOffset();
	z += c->getZOffset();

	double	x_temp = x;
	double	y_temp = y;
	double	z_temp = z;

	double	alpha = c -> getAlpha();
	double	beta  = c -> getBeta();
	double	gamma = c -> getGamma();

	// second the rotation
	// for the inverse matrix derivation see paper log book 19/01/2011
	// libov@mail.desy.de
// no correct any more due to changes in convention of the angle signs!
// further more a more clear way would be to multiply an input vector by an inverse rotation matrix  
// Igor Rubinsky 09-10-2011


	x = x_temp * (1 + alpha * alpha ) + ( (-1) * gamma - alpha * beta) * y_temp + ( (-1) * beta + alpha * gamma) * z_temp;
	y = x_temp * (gamma - alpha * beta ) + (1 + beta * beta) * y_temp + ((-1) * alpha - beta * gamma) * z_temp;
	z = x_temp * (beta + alpha * gamma ) + (alpha - gamma * beta) * y_temp + ( 1 + gamma * gamma) * z_temp;

	double det = 1 + alpha * alpha + beta * beta + gamma * gamma;

	x = x / det;
	y = y / det;
	z = z / det;
}



void EUTelApplyAlignmentProcessor::_EulerRotation(int sensorID, double* _telPos, double* _gRotation) {
   
    try{
        double t = _telPos[2];
    }
    catch(...)
    {
        throw InvalidParameterException("_telPos[] array can not be accessed \n");
    }

    TVector3 _RotatedSensorHit( _telPos[0], _telPos[1], _telPos[2] );

    if( TMath::Abs(_gRotation[2]) > 1e-6 )    _RotatedSensorHit.RotateX( _gRotation[2] ); // in ZY
    if( TMath::Abs(_gRotation[1]) > 1e-6 )    _RotatedSensorHit.RotateY( _gRotation[1] ); // in ZX 
    if( TMath::Abs(_gRotation[0]) > 1e-6 )    _RotatedSensorHit.RotateZ( _gRotation[0] ); // in XY

    _telPos[0] = _RotatedSensorHit.X();
    _telPos[1] = _RotatedSensorHit.Y();
    _telPos[2] = _RotatedSensorHit.Z();
 
}


void EUTelApplyAlignmentProcessor::_EulerRotationInverse(int sensorID, double* _telPos, double* _gRotation) {
   
    try{
        double t = _telPos[2];
    }
    catch(...)
    {
        throw InvalidParameterException("_telPos[] array can not be accessed \n");
    }

    TVector3 _RotatedSensorHit( _telPos[0], _telPos[1], _telPos[2] );

    if( TMath::Abs(_gRotation[0]) > 1e-6 )    _RotatedSensorHit.RotateZ( -1.*_gRotation[0] ); // in XY
    if( TMath::Abs(_gRotation[1]) > 1e-6 )    _RotatedSensorHit.RotateY( -1.*_gRotation[1] ); // in ZX 
    if( TMath::Abs(_gRotation[2]) > 1e-6 )    _RotatedSensorHit.RotateX( -1.*_gRotation[2] ); // in ZY
 

    _telPos[0] = _RotatedSensorHit.X();
    _telPos[1] = _RotatedSensorHit.Y();
    _telPos[2] = _RotatedSensorHit.Z();
 
}

void EUTelApplyAlignmentProcessor::AlignReferenceHit(EUTelEventImpl * evt, EUTelAlignmentConstant * alignment )
{
    int iPlane = alignment->getSensorID();
    streamlog_out(MESSAGE)<< "AlignReferenceHit for " <<  alignment->getSensorID() << endl;
    double            alpha = alignment->getAlpha();
    double            beta  = alignment->getBeta();
    double            gamma = alignment->getGamma();
    double            offsetX = alignment->getXOffset();
    double            offsetY = alignment->getYOffset();
    double            offsetZ = alignment->getZOffset();
                
                streamlog_out ( MESSAGE )  << " alignment->getAlpha() = " << alpha  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getBeta()  = " <<  beta  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getGamma() = " << gamma << endl;
                streamlog_out ( MESSAGE )  << " alignment->getXOffest() = " << offsetX  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getYOffest() = " << offsetY  << endl;
                streamlog_out ( MESSAGE )  << " alignment->getZOffest() = " << offsetZ  << endl;

    if( evt == 0 )
    {
      streamlog_out(ERROR) << "EMPTY event!" << endl;
      return;
    } 
    
    if( _referenceHitVec == 0 )
    {
      streamlog_out(MESSAGE) << "_referenceHitVec is empty" << endl;
    }
   
    streamlog_out(MESSAGE) << "reference Hit collection name : " << _referenceHitCollectionName << endl;
 
      for(size_t ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
      {
        EUTelReferenceHit * refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
        if( iPlane != refhit->getSensorID() )
        {
       //   streamlog_out(MESSAGE) << "Looping through a varity of sensor IDs" << endl;
          continue;
        }
        else
        {
          streamlog_out(MESSAGE) << "Sensor ID and Alignment plane ID match!" << endl;
        }
        printf("PRE sensorID: %5d dx:%5.3f dy:%5.3f dz:%5.3f  alfa:%5.3f beta:%5.3f gamma:%5.3f \n",
        refhit->getSensorID(   ),                      
        refhit->getXOffset(    ),
        refhit->getYOffset(    ),
        refhit->getZOffset(    ),
        refhit->getAlpha(),
        refhit->getBeta(),
        refhit->getGamma()    );

        refhit->setXOffset( refhit->getXOffset() - offsetX );
        refhit->setYOffset( refhit->getYOffset() - offsetY );
        refhit->setZOffset( refhit->getZOffset() - offsetZ );
        TVector3 _RotatedVector( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );

//orig        _RotatedVector.RotateZ(  -alpha        ); // in ZY
//orig        _RotatedVector.RotateX(  -beta         ); // in ZY
//orig        _RotatedVector.RotateY(  -gamma        ); // in XY 
        _RotatedVector.RotateX(  -alpha        ); // in ZY
        _RotatedVector.RotateY(  -beta         ); // in ZY
        _RotatedVector.RotateZ(  -gamma        ); // in XY 
 
        refhit->setAlpha( _RotatedVector[0] );
        refhit->setBeta( _RotatedVector[1] );
        refhit->setGamma( _RotatedVector[2] );
//    referenceHitCollection->push_back( refhit );
 
        printf("AFT sensorID: %5d dx:%5.3f dy:%5.3f dz:%5.3f  alfa:%5.3f beta:%5.3f gamma:%5.3f \n",
        refhit->getSensorID(   ),                      
        refhit->getXOffset(    ),
        refhit->getYOffset(    ),
        refhit->getZOffset(    ),
        refhit->getAlpha(),
        refhit->getBeta(),
        refhit->getGamma()    );
   
      //  _referenceHitVecAligned->push_back( refhit );
     }

}

#endif
