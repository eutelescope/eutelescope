// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 */

// eutelescope includes ".h"
#include "EUTelCalculateEtaProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelEtaFunctionImpl.h"
#include "EUTelPseudo1DHistogram.h"
#include "EUTelExceptions.h"

// marlin includes ".h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include "marlin/AIDAProcessor.h"
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#include <AIDA/IAxis.h>
#endif

// #if defiend(USE_ROOT) || defined(MARLIN_USE_ROOT)
// #include "ROOTProcessor.h"
// #include <TH1.h>
// #endif

// lcio includes <.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

// system includes <>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <cstdlib>
#include <map>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

#ifdef MARLIN_USE_HISTOGRAM
string EUTelCalculateEtaProcessor::_cogHistogramXName = "CoG_X";
string EUTelCalculateEtaProcessor::_cogHistogramYName = "CoG_Y";
string EUTelCalculateEtaProcessor::_cogIntegralXName  = "Integral_CoG_X";
string EUTelCalculateEtaProcessor::_cogIntegralYName  = "Integral_CoG_Y";
string EUTelCalculateEtaProcessor::_etaHistoXName     = "EtaProfile_X";
string EUTelCalculateEtaProcessor::_etaHistoYName     = "EtaProfile_Y";
string EUTelCalculateEtaProcessor::_cogHisto2DName    = "CoG_Histo2D";
#endif

const double EUTelCalculateEtaProcessor::_min = -0.5;
const double EUTelCalculateEtaProcessor::_max =  0.5;

EUTelCalculateEtaProcessor::EUTelCalculateEtaProcessor () : Processor("EUTelCalculateEtaProcessor") {

  // modify processor description
  _description =
    "EUTelCalculateEtaProcessor calculates the eta function for a given set of clusters";

  _isEtaCalculationFinished = false;

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERPULSE, "ClusterCollectionName",
                           "Input cluster collection",
                           _clusterCollectionName, string ("cluster"));

  registerProcessorParameter("EventNumber",
                             "Write here how many events you want to use for eta calculation (-1 for all)",
                             _nEvent, static_cast<int> ( -1 ));

  IntVec noOfBinExample;
  noOfBinExample.push_back(1000);
  noOfBinExample.push_back(1000);
  registerProcessorParameter("NumberOfBins",
                             "Write here in how many bins the seed pixel should be divided (x and y)",
                             _noOfBin, noOfBinExample, noOfBinExample.size());

  registerProcessorParameter("ClusterQualitySelection",
                             "To use only kGoodQuality write 0 here",
                             _clusterQuality, static_cast<int> ( 0 ));
  registerProcessorParameter("ClusterTypeSelection",
                             "Write FULL: full cluster, NxMPixel: for a NxM sub-cluster, NPixel: to use only N pixel",
                             _clusterTypeSelection, string("FULL"));

  IntVec xyCluSizeExample;
  xyCluSizeExample.push_back(3);
  xyCluSizeExample.push_back(3);
  registerProcessorParameter("NxMPixelClusterSize",
                             "The size along x and y of the subcluster (only for NxMPixel)",
                             _xyCluSize, xyCluSizeExample, xyCluSizeExample.size());

  registerProcessorParameter("NPixelSize",
                             "The number of pixel with the highest signal (only for NPixel)",
                             _nPixel, static_cast<int> ( 5 ));

  registerProcessorParameter("EtaXCollectionName",
                             "Set the name of the Eta collection along x",
                             _etaXCollectionName, string("xEtaCondition"));

  registerProcessorParameter("EtaYCollectionName",
                             "Set the name of the Eta collection along y",
                             _etaYCollectionName, string("yEtaCondition"));

  registerProcessorParameter("OutputEtaFileName",
                             "This is the name of the output condition file",
                             _outputEtaFileName, string("etafile"));

  registerOptionalParameter("RejectSinglePixelCluster","reject single pixel cluster. 1=reject, 0=keep, 2=reject clusters with two pixels, where the second pixel is not diagonal to the seed. ",_rejectsingplepixelcluster, static_cast <int> (0));



}


void EUTelCalculateEtaProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // reset stl vectors
  _cogHistogramX.clear();
  _cogHistogramY.clear();
  _integralHistoX.clear();
  _integralHistoY.clear();

  if(_rejectsingplepixelcluster != 0 && _rejectsingplepixelcluster != 1 && _rejectsingplepixelcluster != 2)
    {
      streamlog_out ( ERROR4 ) << "the parameter RejectSinglePixelCluster must set to 0,1 or 2, but it is "<<  _rejectsingplepixelcluster<< endl;
      exit(-1);
    }
}

void EUTelCalculateEtaProcessor::processRunHeader (LCRunHeader * rdr) {
  std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
  runHeader->addProcessor( type() );

  _noOfDetector = runHeader->getNoOfDetector();
  _detectorName = runHeader->lcRunHeader()->getDetectorName();

  int tempEvent;
  if ( Global::parameters->getIntVal("MaxRecordNumber") == 0 ) {
    tempEvent = runHeader->getNoOfEvent();
  } else {
    tempEvent = min( runHeader->getNoOfEvent(), //seems that
                                                //runHeader->getNoOfEvent()
                                                //is always eq 0???
                     Global::parameters->getIntVal("MaxRecordNumber") ) - 1;
  }

  if ( ( _nEvent == -1 ) || ( _nEvent >= tempEvent  && tempEvent > 0 ) ) {
    _nEvent = tempEvent;
  }


  if ( !_isEtaCalculationFinished ) {

    // prepare the output file
    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

    try {
      lcWriter->open( _outputEtaFileName, LCIO::WRITE_NEW );
    } catch (IOException& e) {
      streamlog_out ( ERROR4 ) << e.what() << endl
                               << "Sorry for quitting." << endl;
      exit(-1);
    }

    LCRunHeaderImpl    * lcHeader  = new LCRunHeaderImpl;
    EUTelRunHeaderImpl * newHeader = new EUTelRunHeaderImpl(lcHeader);
    newHeader->lcRunHeader()->setRunNumber(runHeader->lcRunHeader()->getRunNumber());
    newHeader->lcRunHeader()->setDetectorName(runHeader->lcRunHeader()->getDetectorName());
    newHeader->setHeaderVersion(runHeader->getHeaderVersion());
    newHeader->setDataType(runHeader->getDataType());
    newHeader->setDateTime();
    newHeader->setDAQHWName(runHeader->getDAQHWName());
    newHeader->setDAQHWVersion(runHeader->getDAQHWVersion());
    newHeader->setDAQSWName(runHeader->getDAQSWName());
    newHeader->setDAQSWVersion(runHeader->getDAQSWVersion());
    newHeader->setNoOfEvent(runHeader->getNoOfEvent());
    newHeader->setNoOfDetector(runHeader->getNoOfDetector());
    newHeader->setMinX(runHeader->getMinY());
    newHeader->setMaxX(runHeader->getMaxX());
    newHeader->setMinY(runHeader->getMinY());
    newHeader->setMaxY(runHeader->getMaxY());
    newHeader->addProcessor( type());
    newHeader->setGeoID(runHeader->getGeoID());

    lcWriter->writeRunHeader(lcHeader);
    delete newHeader;
    delete lcHeader;
    lcWriter->close();
  }


  // increment the run counter
  ++_iRun;

}


void EUTelCalculateEtaProcessor::processEvent (LCEvent * event) {

  ++_iEvt;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) {

    //streamlog_out ( DEBUG4 ) << "Found kEORE, calling finalizeProcessor()" << endl;
    //finalizeProcessor();
    //return;


    streamlog_out ( DEBUG4 ) << "Found kEORE, calling finalizeProcessor() now!" << endl;
    finalizeProcessor();
    return;

  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  if ( !_isEtaCalculationFinished ) {

    try {
      LCCollectionVec * clusterCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_clusterCollectionName));
      CellIDDecoder<TrackerPulseImpl> cellDecoder(clusterCollectionVec);

      for (int iCluster = 0; iCluster < clusterCollectionVec->getNumberOfElements() ; iCluster++) {

        TrackerPulseImpl   * pulse = dynamic_cast<TrackerPulseImpl *>  ( clusterCollectionVec->getElementAt(iCluster) );
        int temp = cellDecoder(pulse)["type"];
        ClusterType type = static_cast<ClusterType>( temp );

        // all clusters have to inherit from the virtual cluster (that is
        // a TrackerDataImpl with some utility methods).
        EUTelVirtualCluster    * cluster;
        if ( type == kEUTelDFFClusterImpl ) {

          // digital fixed cluster implementation. Remember it can come from
          // both RAW and ZS data
          cluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );

        } else if ( type == kEUTelFFClusterImpl ) {

          // fixed cluster implementation. Remember it can come from
          // both RAW and ZS data
          cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) );
          //streamlog_out ( MESSAGE2 ) <<  "seen a kEUTelFFClusterImpl" << endl; //!HACK TAKI

        } else if ( type == kEUTelBrickedClusterImpl ) {

          // bricked cluster implementation
          // Remember it can come from both RAW and ZS data
          cluster = new EUTelBrickedClusterImpl( static_cast<TrackerDataImpl*> (pulse->getTrackerData()) ); //!HACK TAKI
          //streamlog_out ( MESSAGE2 ) <<  "seen a kEUTelBrickedClusterImpl" << endl; //!HACK TAKI

        } else if ( type == kEUTelSparseClusterImpl ) {

          // ok the cluster is of sparse type, but we also need to know
          // the kind of pixel description used. This information is
          // stored in the corresponding original data collection.

          LCCollectionVec * sparseClusterCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection("original_zsdata"));
          TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt( 0 ));
          CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
          SparsePixelType pixelType = static_cast<SparsePixelType> ( static_cast<int> ( anotherDecoder( oneCluster )["sparsePixelType"] ));

          // now we know the pixel type. So we can properly create a new
          // instance of the sparse cluster
          if ( pixelType == kEUTelGenericSparsePixel ) {
            cluster = new EUTelSparseClusterImpl< EUTelGenericSparsePixel >
              ( static_cast<TrackerDataImpl *> ( pulse->getTrackerData()  ) );
          } else {
            streamlog_out ( ERROR4 ) << "Unknown pixel type.  Sorry for quitting." << endl;
            throw UnknownDataTypeException("Pixel type unknown");
          }
          //streamlog_out ( MESSAGE2 ) <<  "seen a kEUTelSparseClusterImpl" << endl; //!HACK TAKI

        } else {
          streamlog_out ( ERROR4 ) <<  "Unknown cluster type. Sorry for quitting" << endl;
          throw UnknownDataTypeException("Cluster type unknown");
        }

        int detectorID = cluster->getDetectorID();
        float xShift, yShift;



        if ( cluster->getClusterQuality() == static_cast<ClusterQuality> (_clusterQuality) )
          {

            EUTelBrickedClusterImpl* p_tmpBrickedCluster = NULL;
            if ( type == kEUTelBrickedClusterImpl )
              {
                p_tmpBrickedCluster = dynamic_cast< EUTelBrickedClusterImpl* >(cluster);
                //Static of cluster to EUTelBrickedClusterImpl* was done for sure in the case of
                //( type == kEUTelBrickedClusterImpl ).
                //So this cast must work as well!
                //This is just a (different) pointer to the same memory as "cluster". So no additional delete needed.
              }

            if ( _clusterTypeSelection == "FULL" ) {

              if (p_tmpBrickedCluster)
                {
                  //streamlog_out ( MESSAGE2 ) <<  "DEBUG: doing eta FULL on a bricked cluster!" << endl;
                  p_tmpBrickedCluster->getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(xShift, yShift);
                }
              else
                {
                  cluster->getCenterOfGravityShift(xShift, yShift);
                }

            } else if ( _clusterTypeSelection == "NxMPixel" ) {

              if (p_tmpBrickedCluster)
                {
                  streamlog_out ( WARNING4 ) <<  "NxM not applicable for a bricked cluster!! Doing FULL!" << endl;
                  p_tmpBrickedCluster->getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(xShift, yShift);
                }
              else
                {
                  cluster->getCenterOfGravityShift(xShift, yShift, _xyCluSize[0], _xyCluSize[1]);
                }

            } else if ( _clusterTypeSelection == "NPixel" ) {

              if (p_tmpBrickedCluster)
                {
                  //streamlog_out ( MESSAGE2 ) <<  "DEBUG: doing eta NPixel on a bricked cluster!" << endl;
                  p_tmpBrickedCluster->getCenterOfGravityShiftWithOutGlobalSeedCoordinateCorrection(xShift, yShift, _nPixel);
                }
              else
                {
                  cluster->getCenterOfGravityShift(xShift, yShift, _nPixel);
                }

            }


            //#define TAKI_DEBUG_ETA 1
#ifdef TAKI_DEBUG_ETA
            if (p_tmpBrickedCluster)
              {
                streamlog_out ( MESSAGE2 ) << endl;
                streamlog_out ( MESSAGE2 ) <<  "Just done ETA on a BrickedCluster!" << endl;
                p_tmpBrickedCluster->debugOutput();
              }
#endif //TAKI_DEBUG_ETA


            // look for the proper pseudo histogram before filling
            // it. In case the corresponding pseudo histogram is not yet
            // available, book it on the fly!
            if ( _cogHistogramX.find( detectorID ) == _cogHistogramX.end() ) {
              _cogHistogramX[ detectorID ]  = new EUTelPseudo1DHistogram(_noOfBin[0], _min, _max);
              _integralHistoX[ detectorID ] = new EUTelPseudo1DHistogram(_noOfBin[0], _min, _max);
            }
            if ( _cogHistogramY.find( detectorID ) == _cogHistogramY.end() ) {
              _cogHistogramY[ detectorID ] = new EUTelPseudo1DHistogram(_noOfBin[1], _min, _max);
              _integralHistoY[ detectorID ] = new EUTelPseudo1DHistogram(_noOfBin[0], _min, _max);
            }
            //is this a single pixel cluster? 
       
            bool spc_cut = false;
            if( _rejectsingplepixelcluster == 2)
              {
                spc_cut = type != kEUTelDFFClusterImpl &&
                  (abs(static_cast<double>(xShift)) < numeric_limits< double >::min()  ||
                   abs(static_cast<double>(yShift)) <  numeric_limits< double >::min());
              }
            else if(_rejectsingplepixelcluster == 1)
              {
                spc_cut = type != kEUTelDFFClusterImpl &&
                  abs(static_cast<double>(xShift)) < numeric_limits< double >::min()  &&
                  abs(static_cast<double>(yShift)) <  numeric_limits< double >::min(); 
              }

            if(!spc_cut)
              {
                _cogHistogramX[detectorID]->fill(static_cast<double>(xShift), 1.0);
                _cogHistogramY[detectorID]->fill(static_cast<double>(yShift), 1.0);
              }
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            {

              if ( _alreadyBookedSensorID.find( detectorID ) == _alreadyBookedSensorID.end()) {
                // need booking!

                // the path first!
                string path = "detector_" + to_string( detectorID ) ;
                AIDAProcessor::tree( this )->mkdir( path.c_str() );
                path.append( "/" );

                // Center of gravity along X
                string name  = _cogHistogramXName + "_" + to_string(detectorID);
                string title =  "CoG shift along X on detector " + to_string( detectorID );
                AIDA::IHistogram1D * cogHistoX =
                  AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), _noOfBin[0], _min, _max);
                cogHistoX->setTitle(title.c_str());
                _aidaHistoMap.insert( make_pair(name, cogHistoX) );

                // Center of gravity along y
                name  = _cogHistogramYName + "_" + to_string(detectorID);
                title =  "CoG shift along Y on detector " + to_string( detectorID );
                AIDA::IHistogram1D * cogHistoY =
                  AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), _noOfBin[1], _min, _max);
                cogHistoY->setTitle(title.c_str());
                _aidaHistoMap.insert( make_pair(name, cogHistoY) );

                // Integral along x
                name  =  _cogIntegralXName + "_" + to_string( detectorID ) ;
                title = "Integral CoG (x) shift histogram on " + to_string( detectorID );
                AIDA::IHistogram1D * cogIntegralHistoX =
                  AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), _noOfBin[0], _min, _max);
                cogIntegralHistoX->setTitle(title.c_str());
                _aidaHistoMap.insert( make_pair(name, cogIntegralHistoX) );

                // Integral along y
                name  = _cogIntegralYName + "_" + to_string( detectorID ) ;
                title = "Integral CoG (y) shift histogram on " + to_string( detectorID );
                AIDA::IHistogram1D * cogIntegralHistoY =
                  AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), _noOfBin[1], _min, _max);
                cogIntegralHistoY->setTitle(title.c_str());
                _aidaHistoMap.insert( make_pair(name, cogIntegralHistoY) );

                // Eta histogram along x
                name  = _etaHistoXName + "_" + to_string( detectorID ) ;
                title = "Eta profile x for detector " + to_string( detectorID ) ;
                AIDA::IProfile1D * etaHistoX =
                  AIDAProcessor::histogramFactory(this)->createProfile1D( (path + name).c_str(), _noOfBin[0], _min, _max, _min, _max);
                etaHistoX->setTitle(title.c_str());
                _aidaHistoMap.insert( make_pair(name, etaHistoX) );

                // Eta histogram along y
                name  = _etaHistoYName + "_" + to_string( detectorID ) ;
                title = "Eta profile y for detector " + to_string( detectorID ) ;
                AIDA::IProfile1D * etaHistoY =
                  AIDAProcessor::histogramFactory(this)->createProfile1D( (path + name).c_str(), _noOfBin[1], _min, _max, _min, _max);
                etaHistoY->setTitle(title.c_str());
                _aidaHistoMap.insert( make_pair(name, etaHistoY) );

                // 2D histo with CoG
                name  = _cogHisto2DName + "_" + to_string( detectorID ) ;
                title = "2D Histo with the CoG within the seed pixel for detector " + to_string( detectorID ) ;
                AIDA::IHistogram2D * cogHisto2D = AIDAProcessor::histogramFactory(this)
                  ->createHistogram2D( (path + name).c_str(), _noOfBin[0], _min, _max, _noOfBin[1], _min, _max);
                cogHisto2D->setTitle(title.c_str());
                _aidaHistoMap.insert( make_pair(name, cogHisto2D) );

                _alreadyBookedSensorID.insert( detectorID ) ;
              }
           

              if(!spc_cut)
                {
                  string name = _cogHistogramXName + "_" + to_string( detectorID );
                  (dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]))->fill(xShift);

                  name = _cogHistogramYName + "_" + to_string( detectorID );
                  (dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]))->fill(yShift);

                  name = _cogHisto2DName + "_" + to_string( detectorID );
                  (dynamic_cast< AIDA::IHistogram2D* > (_aidaHistoMap[name]))->fill(xShift, yShift);
                }
            }
#endif
          }
        else
          {
            //streamlog_out ( MESSAGE2 ) <<  "CLUSTER QUALITY NOT GOOD ENOUGH!!!" << endl;
          }

        delete cluster;
      }

    } catch ( lcio::DataNotAvailableException & e) {
      return ;
    }
  }

  if(_iEvt>=_nEvent && _nEvent != -1 && !_isEtaCalculationFinished ){
    streamlog_out (  DEBUG4 ) << _nEvent << " events done, calling finalizeProcessor()" << endl;
    finalizeProcessor();
    return;
  }

  if ( isFirstEvent() ) _isFirstEvent = false;

  setReturnValue( "isEtaCalculationFinished" , _isEtaCalculationFinished);

}



void EUTelCalculateEtaProcessor::check (LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelCalculateEtaProcessor::finalizeProcessor() {

  if ( _isEtaCalculationFinished ) return;

  double integral = 0;

  streamlog_out ( MESSAGE4 ) << "Writing the output eta file " << _outputEtaFileName << endl;

  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

  try {
    lcWriter->open( _outputEtaFileName, LCIO::WRITE_APPEND);
  } catch (IOException& e ) {
    streamlog_out ( ERROR4 ) <<  e.what() << endl
                             << "Sorry for quitting." << endl;
    exit(-1);
  }

  LCEventImpl * event = new LCEventImpl();
  event->setDetectorName(_detectorName);
  event->setRunNumber(_iRun);

  LCTime * now = new LCTime;
  event->setTimeStamp(now->timeStamp());
  delete now;

  LCCollectionVec * etaXCollection = new LCCollectionVec(LCIO::LCGENERICOBJECT);
  LCCollectionVec * etaYCollection = new LCCollectionVec(LCIO::LCGENERICOBJECT);

  map< int, EUTelPseudo1DHistogram * >::iterator iter = _cogHistogramX.begin();
  while ( iter != _cogHistogramX.end() ) {

    int iDetector = iter->first;

    for (int iBin = 1; iBin <= _cogHistogramX[iDetector]->getNumberOfBins(); iBin++ ) {
      double x = _cogHistogramX[iDetector]->getBinCenter(iBin);
      integral = _cogHistogramX[iDetector]->integral(1, iBin);
      _integralHistoX[iDetector]->fill(x, integral);

    }

    vector<double > etaBinCenter;
    vector<double > etaBinValue;

    for (int iBin = 1; iBin <= _integralHistoX[iDetector]->getNumberOfBins(); iBin++) {
      etaBinCenter.push_back( _integralHistoX[iDetector]->getBinCenter(iBin) );
      etaBinValue.push_back(  _integralHistoX[iDetector]->getBinContent(iBin) / integral - 0.5 ); // - 0.5);
    }


    EUTelEtaFunctionImpl * etaX = new EUTelEtaFunctionImpl(etaBinCenter.size(), etaBinCenter, etaBinValue);
#if ETA_VERSION >= 2
    etaX->setSensorID( iDetector ) ;
#endif
    etaXCollection->push_back(etaX);

    etaBinCenter.clear();
    etaBinValue.clear();

    for (int iBin = 1; iBin <= _cogHistogramY[iDetector]->getNumberOfBins(); iBin++ ) {
      double y = _cogHistogramY[iDetector]->getBinCenter(iBin);
      integral = _cogHistogramY[iDetector]->integral(1, iBin);
      _integralHistoY[iDetector]->fill(y, integral);
    }

    for (int iBin = 1; iBin <= _integralHistoY[iDetector]->getNumberOfBins(); iBin++) {
      etaBinCenter.push_back( _integralHistoY[iDetector]->getBinCenter(iBin) );
      etaBinValue.push_back(  _integralHistoY[iDetector]->getBinContent(iBin) / integral - 0.5 ); // 0.5);
    }

    EUTelEtaFunctionImpl * etaY = new EUTelEtaFunctionImpl(etaBinCenter.size(), etaBinCenter, etaBinValue);
#if ETA_VERSION >= 2
    etaY->setSensorID( iDetector ) ;
#endif
    etaYCollection->push_back(etaY);


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    integral = 0;
    string name = _cogIntegralXName + "_" + to_string( iDetector );
    AIDA::IHistogram1D * integralHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    name = _cogHistogramXName + "_" + to_string( iDetector );
    AIDA::IHistogram1D * cogHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    integral = 0;
    for (int iBin = 0; iBin < cogHisto->axis().bins(); iBin++) {
      double x = cogHisto->axis().binLowerEdge(iBin) + 0.5 * cogHisto->axis().binWidth(iBin);
      integral += cogHisto->binHeight(iBin);
      integralHisto->fill(x, integral);
    }

    name = _etaHistoXName + "_" + to_string( iDetector );
    AIDA::IProfile1D * etaXHisto = dynamic_cast< AIDA::IProfile1D* > (_aidaHistoMap[name]);

    for (int iBin = 0; iBin < integralHisto->axis().bins(); iBin++) {
      double x = integralHisto->axis().binLowerEdge(iBin) + 0.5 * integralHisto->axis().binWidth(iBin);
      double v = integralHisto->binHeight(iBin);
      etaXHisto->fill(x, v/integral - 0.5 );
    }

    name = _cogIntegralYName + "_" + to_string( iDetector );
    integralHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    name = _cogHistogramYName  + "_" + to_string( iDetector );
    cogHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    integral = 0;
    for (int iBin = 0; iBin < cogHisto->axis().bins(); iBin++) {
      double y = cogHisto->axis().binLowerEdge(iBin) + 0.5 * cogHisto->axis().binWidth(iBin);
      integral += cogHisto->binHeight(iBin);
      integralHisto->fill(y, integral);
    }

    name =  _etaHistoYName  + "_" + to_string( iDetector );
    AIDA::IProfile1D * etaYHisto = dynamic_cast< AIDA::IProfile1D* > (_aidaHistoMap[name]);

    for (int iBin = 0; iBin < integralHisto->axis().bins(); iBin++) {
      double y = integralHisto->axis().binLowerEdge(iBin) + 0.5 * integralHisto->axis().binWidth(iBin);
      double v = integralHisto->binHeight(iBin);
      etaYHisto->fill(y, v/integral - 0.5 );
    }

#endif

    ++iter;
  }

  event->addCollection(etaXCollection, _etaXCollectionName);
  event->addCollection(etaYCollection, _etaYCollectionName);

  lcWriter->writeEvent(event);
  delete event;

  lcWriter->close();


  _isEtaCalculationFinished = true;
  setReturnValue( "isEtaCalculationFinished" , _isEtaCalculationFinished);
  //  throw RewindDataFilesException(this);

}


void EUTelCalculateEtaProcessor::end() {

  if(!_isEtaCalculationFinished ){
    streamlog_out (  DEBUG4 ) << "calling finalizeProcessor()" << endl;
    finalizeProcessor();
  }

  streamlog_out ( MESSAGE2 ) <<  "Successfully finished" << endl;
}

