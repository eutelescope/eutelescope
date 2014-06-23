
// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version: $Id$
// Date 2007.09.10

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined(USE_GEAR) && ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

// eutelescope inlcudes
#include "EUTelFitTuple.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"

// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>


#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static members mainly used to name histograms
std::string EUTelFitTuple::_FitTupleName  = "EUFit";


EUTelFitTuple::EUTelFitTuple() : Processor("EUTelFitTuple") {

  // modify processor description
  _description = "Prepare n-tuple with track fit results" ;


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACK,
                           "InputCollectionName" ,
                           "Name of the input Track collection"  ,
                           _inputColName ,
                           std::string("telescopetracks") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                           "InputDUTCollectionName" ,
                           "Name of the input DUT hit collection"  ,
                           _inputDUTColName ,
                           std::string("duthit") ) ;

  // other processor parameters:

  registerProcessorParameter ("MissingValue",
                              "Value used for missing measurements",
                              _missingValue,  static_cast < double > (-100.));


  registerProcessorParameter ("UseManualDUT",
                              "Flag for manual DUT selection",
                              _useManualDUT,  static_cast < bool > (false));

  registerProcessorParameter ("ManualDUTid",
                              "Id of sensor layer which should be used as DUT",
                              _manualDUTid,  static_cast < int > (0));

  registerProcessorParameter ("DistMax",
                              "Maximum allowed distance between fit and matched DUT hit",
                              _distMax,  static_cast < double > (0.1));

  registerProcessorParameter ("MaxZDistance",
                              "Maximum allowed distance between hit Z coordinate and the plane's Z coordinate",
                              _maxZDistance,  static_cast < double > (1.));


  std::vector<float > initAlign;
  initAlign.push_back(0.);
  initAlign.push_back(0.);
  initAlign.push_back(0.);

  registerProcessorParameter ("DUTalignment",
                              "Alignment corrections for DUT: shift in X, Y and rotation around Z",
                              _DUTalign, initAlign);

}


void EUTelFitTuple::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _tluTimeStamp = 0 ;


  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR5> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  // Read geometry information from GEAR

  message<MESSAGE5> ( log() << "Reading telescope geometry description from GEAR ") ;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));


// Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

// Check for DUT

  if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )
    {
      _iDUT = _nTelPlanes ;
      _nTelPlanes++;
    }
  else
    _iDUT = -1 ;

// Read position in Z (for sorting)

  _planeSort = new int[_nTelPlanes];
  _planePosition   = new double[_nTelPlanes];

  for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
    {
      _planePosition[ipl]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
      _planeSort[ipl]=ipl;
    }

  if(_iDUT>0)
    {
      _planePosition[_iDUT]=_siPlanesLayerLayout->getDUTPositionZ();
      _planeSort[_iDUT]=_iDUT;
    }

  // Binary sorting

  bool sorted;
  do{
    sorted=false;
    for(int iz=0; iz<_nTelPlanes-1 ; iz++)
      if(_planePosition[iz]>_planePosition[iz+1])
        {
          double _posZ = _planePosition[iz];
          _planePosition[iz] = _planePosition[iz+1];
          _planePosition[iz+1] = _posZ;

          int _idZ = _planeSort[iz];
          _planeSort[iz] = _planeSort[iz+1];
          _planeSort[iz+1] = _idZ;

          sorted=true;
        }

  }while(sorted);

// Book local geometry arrays

  _planeID         = new int[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];

// Fill remaining layer parameters

  for(int iz=0; iz < _nTelPlanes ; iz++)
    {
      int ipl=_planeSort[iz];

      double resolution;

      if(ipl != _iDUT )
        {
          _planeID[iz]=_siPlanesLayerLayout->getID(ipl);
          resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
        }
      else
        {
          _planeID[iz]=_siPlanesLayerLayout->getDUTID();
          resolution = _siPlanesLayerLayout->getDUTSensitiveResolution();
        }

      _isActive[iz] = (resolution > 0);


    }

  // Get new DUT position (after sorting)

  for(int iz=0;iz< _nTelPlanes ; iz++)
    if(_planeSort[iz]==_iDUT)
      {
        _iDUT=iz;
        break;
      }

  // DUT position can be changed by processor parameter

  if(_useManualDUT)
    {
      bool _manualOK=false;

      for(int iz=0; iz < _nTelPlanes ; iz++)
        if(_planeID[iz]==_manualDUTid)
          {
            _iDUT=iz;
            _manualOK=true;
          }

      if(!_manualOK)
        {
          message<ERROR5> ( log() << "Manual DUT flag set, layer not found ID = "
                           << _manualDUTid
                           << "\n Program will terminate! Correct geometry description!");
          exit(-1);
        }
    }

  _zDUT=_planePosition[_iDUT];

  // Print out geometry information

  message<MESSAGE5> ( log() << "Telescope configuration with " << _nTelPlanes << " planes" );


  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      stringstream ss ;
      if(ipl == _iDUT)
        ss << "D.U.T.  plane" ;
      else
        if(_isActive[ipl])
          ss << "Active  plane" ;
        else
          ss << "Passive plane" ;

      ss << "  ID = " << _planeID[ipl]
         << "  at Z [mm] = " << _planePosition[ipl];

      message<MESSAGE5> ( log() << ss.str() );
    }


  // Allocate arrays for track fitting

  _isMeasured      = new bool[_nTelPlanes];
  _isFitted        = new bool[_nTelPlanes];

  _measuredX     = new double[_nTelPlanes];
  _measuredY     = new double[_nTelPlanes];
  _measuredZ     = new double[_nTelPlanes];
  _measuredQ     = new double[_nTelPlanes];
  _fittedX       = new double[_nTelPlanes];
  _fittedY       = new double[_nTelPlanes];
  _fittedEX       = new double[_nTelPlanes];
  _fittedEY       = new double[_nTelPlanes];
  _nhits         = new int[_nTelPlanes];
  _ClusterSize         = new double[_nTelPlanes];



// Book histograms
  bookHistos();


}

void EUTelFitTuple::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  _runNr = runHeader->getRunNumber();

  message<MESSAGE5> ( log() << "Processing run header " << _nRun
                     << ", run nr " << _runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  message<MESSAGE5> ( log() << detectorName << " : " << detectorDescription ) ;

  int nDet = subDets->size();

  if(nDet)message<MESSAGE5> ( log() << nDet << " subdetectors defined :" );
  stringstream ss;
  for(int idet=0;idet<nDet;idet++)  message<MESSAGE5> (log()  << idet+1 << " : " << subDets->at(idet) );


}

void EUTelFitTuple::processEvent( LCEvent * event ) {

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
  if ( euEvent->getEventType() == kEORE ) {
    message<DEBUG5> ( "EORE found: nothing else to do." );
    return;
  }

  _nEvt ++ ;
  _evtNr        = event->getEventNumber();
  _tluTimeStamp = static_cast<long int> (event->getTimeStamp());

  LCCollection* col;
  try {
    col = event->getCollection( _inputColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputColName << "from event " << event->getEventNumber() << " in run " << event->getRunNumber() << endl;
    return;
  }


  LCCollection* hitcol = NULL;
  bool _DUTok=true;

  try {
    hitcol = event->getCollection( _inputDUTColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
//     message<ERROR5> ( log() << "Not able to get collection "
//                      << _inputDUTColName
//                      << "\nfrom event " << event->getEventNumber()
//                      << " in run " << event->getRunNumber()  );
    _DUTok=false;
  }


  // Loop over tracks in input collections

  int nTrack = col->getNumberOfElements()  ;

  message<DEBUG5> ( log() << "Total of " << nTrack << " tracks in input collection " );

  int nDUT = 0;
  if(_DUTok) nDUT = hitcol->getNumberOfElements()  ;

  message<DEBUG5> ( log() << "Total of " << nDUT << " hits in input collection " );

  // nhits for each plane
  for(int i=0;i<_nTelPlanes;i++) _nhits[i] = 0;
  for(int ihit=0; ihit< nDUT ; ihit++)
  {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihit) ) ;
      // Hit position
      const double * pos = meshit->getPosition();
      // We find plane number of the hit by looking at the Z position
      double distMin = 1.; int hitPlane = -1 ;
      for(int ipl=0;ipl<_nTelPlanes;ipl++)
      {
          double dist =  pos[2] - _planePosition[ipl] ;
          if(dist*dist < distMin*distMin)
          {
              hitPlane=ipl;
              distMin=dist;
              if(hitPlane!=-1) _nhits[hitPlane]++;
          }
      }
  }

  for(int itrack=0; itrack< nTrack ; itrack++)
    {
      Track * fittrack = dynamic_cast<Track*>( col->getElementAt(itrack) ) ;

      // Hit list assigned to track

      std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

      // Copy hits assign to the track to local table
      // Assign hits to sensor planes


      int nHit =   trackhits.size();

      message<DEBUG5> ( log() << "Track " << itrack << " with " << nHit << " hits, Chi2 = "
                                << fittrack->getChi2() << "/" << fittrack->getNdf());


      // Clear plane tables

      for(int ipl=0;ipl<_nTelPlanes;ipl++)
        {
          _isMeasured[ipl]=false;
          _isFitted[ipl]=false;

          _measuredX[ipl]=_missingValue;
          _measuredY[ipl]=_missingValue;
          _measuredZ[ipl]=_missingValue;
          _measuredQ[ipl]=_missingValue;

          _fittedX[ipl]=_missingValue;
          _fittedY[ipl]=_missingValue;
          _fittedEX[ipl]=_missingValue;
          _fittedEY[ipl]=_missingValue;

          _ClusterSize[ipl]=0;


        }

      // setup cellIdDecoder to decode the hit properties
      CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

      // Loop over hits and fill hit tables

      for(int ihit=0; ihit< nHit ; ihit++)
        {
          TrackerHit * meshit = trackhits.at(ihit);

          // Hit position

          const double * pos = meshit->getPosition();

          // We find plane number of the hit
          // by looking at the Z position

          int hitPlane = -1 ;
	  // Initialise variable for comparing distance to planes:
	  double distMin = _maxZDistance;

          for(int ipl=0;ipl<_nTelPlanes;ipl++)
            {
              double dist =  pos[2] - _planePosition[ipl] ;

              if(dist*dist < distMin*distMin)
                {
                  hitPlane=ipl;
                  distMin=dist;
                }
            }

          // Ignore hits not matched to any plane

          if(hitPlane<0)
            {
              message<ERROR5> ( log() << "Hit outside telescope plane at z [mm] = "  << pos[2] );
              continue;
            }


          if( (hitCellDecoder(meshit)["properties"] & kFittedHit) == 0 )
            {
              // Measured hits

              _isMeasured[hitPlane]=true;

              _measuredX[hitPlane]=pos[0];
              _measuredY[hitPlane]=pos[1];
              _measuredZ[hitPlane]=pos[2];

              // Get cluster charge

              _measuredQ[hitPlane]=0.;

              EVENT::LCObjectVec rawdata =  meshit->getRawHits();

              if(rawdata.size()>0 && rawdata.at(0)!=NULL ) {
		EUTelSparseClusterImpl< EUTelSimpleSparsePixel > * cluster =
		  new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > ( static_cast<TrackerDataImpl*> (rawdata.at(0))) ;
		_measuredQ[hitPlane]=cluster->getTotalCharge();

		// Get cluster size
		_ClusterSize[hitPlane]=(double)cluster->size();
	      }

              message<DEBUG5> ( log() << "Measured hit in plane " << hitPlane << " at  X = "
                                        << pos[0] << ", Y = " << pos[1] << ", Q = " << _measuredQ[hitPlane] );

            }
          else
            {
              // Fitted hits
              const EVENT::FloatVec cov = meshit->getCovMatrix();

              _isFitted[hitPlane]=true;

              _fittedX[hitPlane]=pos[0];
              _fittedY[hitPlane]=pos[1];
              if(cov.at(0)>0.) _fittedEX[hitPlane]=sqrt(cov.at(0));
              if(cov.at(0)>0.) message<DEBUG> ( log()<<"Test : "<<_fittedEX[hitPlane]);
              if(cov.at(2)>0.) _fittedEY[hitPlane]=sqrt(cov.at(2));

              message<DEBUG5> ( log() << "Fitted  hit  in plane " << hitPlane << " at  X = "
                                        << pos[0] << ", Y = " << pos[1] );

            }

        }

      // Fill n-tuple

      int icol=0;
      _FitTuple->fill(icol++,_nEvt);
      _FitTuple->fill(icol++,_runNr);
      _FitTuple->fill(icol++,_evtNr);
      _FitTuple->fill(icol++,_tluTimeStamp); // new! TLU timestamp
      _FitTuple->fill(icol++,nTrack); // new! TLU timestamp
      _FitTuple->fill(icol++,fittrack->getNdf());
      _FitTuple->fill(icol++,fittrack->getChi2());

      for(int ipl=0; ipl<_nTelPlanes;ipl++)
        {
          _FitTuple->fill(icol++,_measuredX[ipl]);
          _FitTuple->fill(icol++,_measuredY[ipl]);
          _FitTuple->fill(icol++,_measuredZ[ipl]);
          _FitTuple->fill(icol++,_measuredQ[ipl]);
          _FitTuple->fill(icol++,_fittedX[ipl]);
          _FitTuple->fill(icol++,_fittedY[ipl]);
          _FitTuple->fill(icol++,_fittedEX[ipl]);
          _FitTuple->fill(icol++,_fittedEY[ipl]);
          _FitTuple->fill(icol++,(double)_nhits[ipl]);
          _FitTuple->fill(icol++,(double)_ClusterSize[ipl]);

        }

      //  Look for closest DUT hit


      double dutX=_missingValue;
      double dutY=_missingValue;
      double dutR=_missingValue;
      double dutQ=_missingValue;
      double dutClusterSize=_missingValue;

      if(_DUTok)
        {
          double distmin=_distMax*_distMax;
          int imin=-1;

          for(int ihit=0; ihit< nDUT ; ihit++)
            {
              TrackerHit * meshit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihit) ) ;

              // Hit position

              const double * pos = meshit->getPosition();

              double dist = pos[2] - _zDUT;

              if(dist*dist < 1)
                {
                  // Apply alignment corrections

                  double corrX  =   pos[0]*cos(_DUTalign.at(2))
                    +pos[1]*sin(_DUTalign.at(2))
                    +_DUTalign.at(0);

                  double corrY  =  pos[1]*cos(_DUTalign.at(2))
                    -pos[0]*sin(_DUTalign.at(2))
                    +_DUTalign.at(1);

                  double distXY=
                    (corrX-_fittedX[_iDUT])*(corrX-_fittedX[_iDUT])
                    + (corrY-_fittedY[_iDUT])*(corrY-_fittedY[_iDUT]);

                  if(distXY<distmin)
                    {
                      imin=ihit;
                      distmin=distXY;
                      dutX=corrX;
                      dutY=corrY;
                    }

                }

            }

          // Try to get DUT cluster charge

          if(imin>=0)
            {
              TrackerHit * meshit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(imin) ) ;

              EVENT::LCObjectVec rawdata =  meshit->getRawHits();

              if(rawdata.size()>0 && rawdata.at(0)!=NULL )
                {
                  //EUTelVirtualCluster * cluster = new EUTelFFClusterImpl ( static_cast<TrackerDataImpl*> (rawdata.at(0))) ;
                  EUTelSparseClusterImpl< EUTelSimpleSparsePixel > * cluster =
                      new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > ( static_cast<TrackerDataImpl*> (rawdata.at(0))) ;
                  dutQ=cluster->getTotalCharge();
                  // Get cluster size
                  dutClusterSize=(double)cluster->size();
                }

              dutR=sqrt(distmin);

              message<DEBUG5> ( log() << "Matched DUT hit at X = " << dutX << "   Y = " << dutY
                                        << "   Dxy = " << dutR << "   Q = " << dutQ );
            }
          else
            message<DEBUG5> ( log() << "DUT hit not matched !" );


          // End of if(_DUTok)
        }


      _FitTuple->fill(icol++,dutX);
      _FitTuple->fill(icol++,dutY);
      _FitTuple->fill(icol++,dutR);
      _FitTuple->fill(icol++,dutQ);
      _FitTuple->fill(icol++,dutClusterSize);
      _FitTuple->fill(icol++,nTrack);

      _FitTuple->addRow();

      // End of loop over tracks
    }


  return;
}



void EUTelFitTuple::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelFitTuple::end(){

  //   std::cout << "EUTelFitTuple::end()  " << name()
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;


  message<MESSAGE5> ( log() << "N-tuple with "
                     << _FitTuple->rows() << " rows created" );


  // Clean memory

  delete [] _planeSort ;
  delete [] _planePosition ;
  delete [] _planeID ;
  delete [] _isActive ;

  delete [] _isMeasured ;
  delete [] _isFitted ;
  delete [] _measuredX  ;
  delete [] _measuredY  ;
  delete [] _measuredZ  ;
  delete [] _measuredQ  ;
  delete [] _fittedX ;
  delete [] _fittedY ;
  delete [] _fittedEX ;
  delete [] _fittedEY ;
  delete [] _nhits ;
  delete [] _ClusterSize ;



}



void EUTelFitTuple::bookHistos()
{


  message<MESSAGE5> ( log() << "Booking fit n-tuple \n" );

  std::vector<std::string> _columnNames;
  std::vector<std::string> _columnType;

  _columnNames.push_back("Event");
  _columnType.push_back("int");

  _columnNames.push_back("RunNr");
  _columnType.push_back("int");

  _columnNames.push_back("EvtNr");
  _columnType.push_back("int");

  _columnNames.push_back("TLUtime");
  _columnType.push_back("long int");

  _columnNames.push_back("Track");
  _columnType.push_back("int");

  _columnNames.push_back("Ndf");
  _columnType.push_back("int");

  _columnNames.push_back("Chi2");
  _columnType.push_back("float");

  const char * _varName[] = { "measX", "measY" , "measZ", "measQ", "fitX", "fitY","fitEX","fitEY","nhits","ClusterSize" };

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    for(int ivar=0; ivar<10;ivar++)
      {
        stringstream ss;
        ss << _varName[ivar] << "_" << ipl;
        _columnNames.push_back(ss.str());
        _columnType.push_back("double");
      }

  // DUT variables

  _columnNames.push_back("dutX");
  _columnType.push_back("double");

  _columnNames.push_back("dutY");
  _columnType.push_back("double");

  _columnNames.push_back("dutR");
  _columnType.push_back("double");

  _columnNames.push_back("dutQ");
  _columnType.push_back("double");

  _columnNames.push_back("dutClusterSize");
  _columnType.push_back("double");

  _columnNames.push_back("nTrack");
  _columnType.push_back("int");

  _FitTuple=AIDAProcessor::tupleFactory(this)->create(_FitTupleName, _FitTupleName, _columnNames, _columnType, "");


  message<DEBUG5> ( log() << "Booking completed \n\n");

  return;
}

#endif // GEAR && AIDA
