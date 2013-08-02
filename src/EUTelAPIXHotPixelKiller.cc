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
#include "EUTelAPIXHotPixelKiller.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseData2Impl.h"
#include "EUTelSparseCluster2Impl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include "marlin/AIDAProcessor.h"
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#endif

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IO/LCWriter.h>

#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>
// lccd
//#include <lccd/DBInterface.hh>

// system includes <>
#include <map>
#include <memory>

using namespace std;
using namespace marlin;
using namespace eutelescope;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
string EUTelAPIXHotPixelKiller::_firing2DHistoName = "Firing2D";
string EUTelAPIXHotPixelKiller::_firing1DHistoName = "Firing1D";
#endif


EUTelAPIXHotPixelKiller::EUTelAPIXHotPixelKiller () : Processor("EUTelAPIXHotPixelKiller") 
{

  // modify processor description
  _description =
    "EUTelAPIXHotPixelKiller periodically check for pixel singing loud too often and remove them from the analysis";


  registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                           "Input of Zero Suppressed data",
                           _zsDataCollectionName, string ("zsdata") );

  registerInputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
                           "Noise (input) collection name",
                           _noiseCollectionName, string("noise"));

  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                           "Pixel status (input) collection name",
                           _statusCollectionName, string("status"));

 
  registerProcessorParameter("BuildHotPixelDatabase",
                             "This flag is used to initialise simple data decoding and hot pixel finder (0-no, 1-yes)",
                             _flagBuildHotPixelDatabase, static_cast< int>( 0 ) );

  registerProcessorParameter("NoOfEventPerCycle",
                             "The number of events to be considered for each update cycle",
                             _noOfEventPerCycle, static_cast<int>( 100 ) );

  registerProcessorParameter("MaxAllowedFiringFreq",
                             "This float number [0,1] represents the maximum allowed firing frequency\n"
                             "within the selected number of event per cycle",
                             _maxAllowedFiringFreq, static_cast<float> (0.2) );

  registerProcessorParameter("TotalNoOfCycle",
                             "The total number of hot pixel cycle",
                             _totalNoOfCycle, static_cast<int> ( 10 ) );

  registerOptionalParameter("HotPixelDBFile","This is the name of the LCIO file name with the output hotpixel db (add .slcio)",
                             _hotpixelDBFile, static_cast< string > ( "hotpixel.slcio" ) );

  registerOptionalParameter("ExcludedPlanes", "The list of sensor ids that have to be excluded from the clustering.",
                             _ExcludedPlanes, std::vector<int> () );

  registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection to be saved into the output slcio file",
                             _hotPixelCollectionName, static_cast< string > ( "hotpixel" ));


}



void EUTelAPIXHotPixelKiller::init () 
{
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // reset the cycle number
  _iCycle = 0;

  // reset the vector with killed pixels
  _killedPixelVec.clear();

  // reset the vector with the firing frequency
  _firingFreqVec.clear();

  // reset hotpixel map vectors
  _hitIndexMapVec.clear();
  _inverse_hitIndexMapVec.clear();
  _pixelMapVec.clear();

}

void EUTelAPIXHotPixelKiller::modifyEvent( LCEvent * /*event*/ )
{
  //This function does nothing, but must be included due to being inherited from an abstract base class
  return;
}

void EUTelAPIXHotPixelKiller::processRunHeader (LCRunHeader * rdr ) {

  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl ( rdr ) );
  runHeader->addProcessor( type() );

  // increment the run counter
  ++_iRun;

  // reset the event counter
  _iEvt = 0;

}


void EUTelAPIXHotPixelKiller::initializeGeometry( LCEvent * event ) 
{

    LCCollectionVec * collection = dynamic_cast< LCCollectionVec *> ( event->getCollection( _statusCollectionName ) );
    if ( collection == 0 )
    { 
        streamlog_out( WARNING2 ) <<  " statusCollectionVec is not found " <<  endl;       
    }

    _noOfDetectors = collection->size();

    CellIDDecoder<TrackerRawDataImpl > decoder( collection );

    for ( size_t iDetector = 0 ; iDetector < collection->size() ; ++iDetector ) 
    {
        TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( collection->getElementAt( iDetector ) ) ;
        int sensorID = decoder( status ) [ "sensorID" ] ;

        _minX[ sensorID ] = decoder( status ) [ "xMin" ];
        _minY[ sensorID ] = decoder( status ) [ "yMin" ];
        _maxX[ sensorID ] = decoder( status ) [ "xMax" ];
        _maxY[ sensorID ] = decoder( status ) [ "yMax" ];

        _sensorIDVec.push_back( sensorID );
    }  


    // 
    // now another map relating the position in the ancillary
    // collections (noise, pedestal and status) with the sensorID
    // 
    _ancillaryIndexMap.clear();

    // 
    // if HotPixelKiller is intended to run standalone 
    // open the noise collection and create the _ancillaryIndexMap
    //     
    if( getBuildHotPixelDatabase() != 0 )    
    {
        try
        {
            // this is the exemplary ancillary collection
            LCCollectionVec * noiseCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _noiseCollectionName ) );

            // prepare also a cell decoder
            CellIDDecoder< TrackerDataImpl > noiseDecoder( noiseCollectionVec );

            for ( size_t iDetector = 0 ; iDetector < noiseCollectionVec->size(); ++iDetector ) 
            {
                TrackerDataImpl * noise = dynamic_cast< TrackerDataImpl * > ( noiseCollectionVec->getElementAt ( iDetector ) );
                if( noise == 0 )
                {
                    streamlog_out( WARNING2 ) <<  " noise TrackerDataImpl is not found " <<  endl;
                }
                _ancillaryIndexMap.insert( make_pair( noiseDecoder( noise ) ["sensorID"], iDetector ) );
            }
        } catch (  lcio::DataNotAvailableException ) {
            streamlog_out( WARNING2 ) << "Unable to initialize the geometry. Trying with the following event" << endl;
            throw SkipEventException( this ) ;
        }
    }
}


// 
// this method will run only standalone
// 
void EUTelAPIXHotPixelKiller::HotPixelFinder(EUTelEventImpl  *evt)
{
    if (evt == 0 )
    {
        exit(-1);
    }

    // get the collections of interest from the event.
    LCCollectionVec * zsInputCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _zsDataCollectionName ));
    LCCollectionVec * statusCollectionVec   = dynamic_cast < LCCollectionVec * > (evt->getCollection( _statusCollectionName ));
    LCCollectionVec * noiseCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection( _noiseCollectionName ));

    // prepare some decoders
    CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
    CellIDDecoder<TrackerDataImpl> statusDecoder( statusCollectionVec );
    CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );


    for ( unsigned int iDetector = 0 ; iDetector < zsInputCollectionVec->size(); iDetector++ ) 
    {
        
        // get the TrackerData and guess which kind of sparsified data it
        // contains.

        TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( iDetector ) );
        SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );

        if (type != kEUTelAPIXSparsePixel  ) 
        {
          std::cout << " pixel is not of APIX type " << std::endl ;
        }

        int _sensorID            = static_cast<int > ( cellDecoder( zsData )["sensorID"] );
        int  sensorID            = _sensorID;


        //if this is an excluded sensor go to the next element

        bool foundexcludedsensor = false;
        for(size_t j = 0; j < _ExcludedPlanes.size(); ++j)
        {
            if(_ExcludedPlanes[j] == _sensorID)
            {
                foundexcludedsensor = true;
            }
        }
        if(foundexcludedsensor)  continue;

        // get the noise and the status matrix with the right detectorID
        TrackerRawDataImpl * status = 0;

        // the noise map. we only need this map for decoding issues.
        TrackerDataImpl    * noise  = 0;
   
        // get the noise and the status matrix with the right detectorID
        status = dynamic_cast<TrackerRawDataImpl*>(statusCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));        
        vector< short > statusVec = status->adcValues();

        //the noise map. we only need this map for decoding issues.
        noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));

        // prepare the matrix decoder
        EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );

        // now prepare the EUTelescope interface to sparsified data.  
//old        auto_ptr<EUTelSparseDataImpl<EUTelSimpleSparsePixel > >  sparseData(new EUTelSparseDataImpl<EUTelSimpleSparsePixel> ( zsData ));
        auto_ptr<EUTelSparseDataImpl<EUTelAPIXSparsePixel > > apixData(new EUTelSparseDataImpl<EUTelAPIXSparsePixel> ( zsData ));

        streamlog_out ( DEBUG1 ) << "Processing sparse data on detector " << _sensorID << " with "
                                 << apixData->size() << " pixels " << endl;
        
        for ( unsigned int iPixel = 0; iPixel < apixData->size(); iPixel++ ) 
        {
            // loop over all pixels in the apixData object.      
            EUTelAPIXSparsePixel *sparsePixel =  new EUTelAPIXSparsePixel() ;

            apixData->getSparsePixelAt( iPixel, sparsePixel );
            int decoded_XY_index = matrixDecoder.getIndexFromXY( sparsePixel->getXCoord(), sparsePixel->getYCoord() ); // unique pixel index !!

            if( _hitIndexMapVec[iDetector].find( decoded_XY_index ) == _hitIndexMapVec[iDetector].end() )
            {
                int old_size = status->adcValues().size();
                // increment
              
                int new_size     = old_size + 1 ;
                int last_element = old_size;
              
                status->adcValues().resize( new_size  );
                _hitIndexMapVec[iDetector].insert ( make_pair ( decoded_XY_index, last_element ) );
                _inverse_hitIndexMapVec[iDetector].insert ( make_pair ( last_element, decoded_XY_index ) );

              
                status->adcValues()[ last_element ] = EUTELESCOPE::HITPIXEL ;  // adcValues is a vector, there fore must address the elements incrementally
                
                _pixelMapVec[iDetector].insert ( make_pair( decoded_XY_index, sparsePixel ) );                     // one more map, get the pixel point bny its unique index
//                printf("--last_element:%7d;  pixel %7d, index %7d, pointer %7d %7d \n",
//                        last_element, iPixel, decoded_XY_index, _pixelMapVec[iDetector][ decoded_XY_index]->getXCoord(), _pixelMapVec[iDetector][ decoded_XY_index]->getYCoord()  );
            }
            else
            {
//                printf("idet: %5d, status: %p, addressing known pixel decoded: %7d  orig: %7d  status->size:%7d\n", 
//                        iDetector, status, decoded_XY_index, _hitIndexMapVec[iDetector][ decoded_XY_index ], status->adcValues().size()  );
                status->adcValues()[ _hitIndexMapVec[iDetector][ decoded_XY_index]  ] = EUTELESCOPE::HITPIXEL ;
            }
        }
    }    

}



void EUTelAPIXHotPixelKiller::processEvent (LCEvent * event) 
{

    if( event == 0 )
    {
        streamlog_out ( WARNING2 ) <<  "event does not exist!. skip " <<  endl;       
        return;
    }
    
    if ( _iCycle > static_cast< unsigned short >( _totalNoOfCycle ) ) return;

    if (_iEvt % 1000 == 0)
        streamlog_out( MESSAGE4 ) << "Processing event "
            << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
            << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
            << " (Total = " << setw(10) << (_iCycle * _noOfEventPerCycle) + _iEvt << ")"
            << resetiosflags(ios::left) << endl;
    
    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
    if ( evt->getEventType() == kEORE ) 
    {
        streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  endl;
        return;
    }
    else if ( evt->getEventType() == kUNKNOWN ) 
    {
        streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber()
            << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }

    try 
    {
        LCCollectionVec * statusCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _statusCollectionName ) );
        CellIDDecoder<TrackerRawDataImpl>      statusCellDecoder( statusCollectionVec );

    // if first event initilize 
    // sensorIDVec and sensors boundaries [min,max] for [X, Y]
    //  
    if ( isFirstEvent() ) 
    {
        streamlog_out ( MESSAGE5 ) << "in initializeGeometry : " << event << "; will prepare map of pixels too " << endl;;
        initializeGeometry( event );

        _firingFreqVec.clear();
        if( getBuildHotPixelDatabase() != 0 )
        {
            _hitIndexMapVec.clear();
            _inverse_hitIndexMapVec.clear();
            _pixelMapVec.clear();
        }
        
        for ( int iDetector = 0; iDetector < statusCollectionVec->getNumberOfElements() ; iDetector++) 
        {  
           streamlog_out ( MESSAGE5 ) << " First event :: adding Detector " << iDetector << endl;          
           _firingFreqVec.resize( iDetector+1 );
           
           if( getBuildHotPixelDatabase() != 0 )
            {
                if(static_cast< int >( _hitIndexMapVec.size()) < iDetector+1 )            _hitIndexMapVec.resize(iDetector+1);
                if( static_cast< int >(_inverse_hitIndexMapVec.size()) < iDetector+1 )    _inverse_hitIndexMapVec.resize(iDetector+1);
                if( static_cast< int >(_pixelMapVec.size()) < iDetector+1 )               _pixelMapVec.resize(iDetector+1);
            }
        }
        
        _isFirstEvent = false;
    }
    
    if( getBuildHotPixelDatabase() != 0 )
    {
        HotPixelFinder(evt);
    }

     
    
    for ( int iDetector = 0; iDetector < statusCollectionVec->getNumberOfElements() ; iDetector++) 
    {
        TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( statusCollectionVec->getElementAt( iDetector ) );
        if( _firingFreqVec[iDetector].size() < status->getADCValues().size() )
        {
            _firingFreqVec[iDetector].resize( status->getADCValues().size() )  ;
        }
    }
    

    for ( int iDetector = 0; iDetector < statusCollectionVec->getNumberOfElements() ; iDetector++) 
    {
        TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( statusCollectionVec->getElementAt( iDetector ) );
        vector< short > statusVec = status->adcValues();
        streamlog_out ( DEBUG3 ) << 
            " freq loop: idet=" << iDetector << 
            " statusVec.size=" <<  statusVec.size() << endl;
        
        // some incremental index value to status of each pixel 
        for ( unsigned int index = 0; index < statusVec.size(); index++ ) 
        {
            streamlog_out ( DEBUG3 )  << 
                " index "<< index << 
                " all " << statusVec.size() << 
                " statusVec[index]: " <<  statusVec[index] << 
                " hit:" << EUTELESCOPE::HITPIXEL << endl;
            streamlog_out ( DEBUG3)  << 
                " get index " << index << 
                " (tot " << statusVec.size() << ")" <<
                " (orig decoded_XY_index = " << _inverse_hitIndexMapVec[iDetector][index] <<
                " freq " << _firingFreqVec[iDetector][index] <<
                " status " <<  statusVec[ index ]  << endl;
            
           if( statusVec[ index ] == EUTELESCOPE::HITPIXEL ) 
            {
                _firingFreqVec[ iDetector ][ index ] += 1;
                status->adcValues()[ index ] = EUTELESCOPE::GOODPIXEL;
                if( _firingFreqVec[ iDetector ][ index ]  > 1   )
                {
                    streamlog_out ( DEBUG3) << 
                        " idet: " << iDetector <<
                        " get index " << index <<
                        " (tot " << statusVec.size() << ")" <<
                        " (orig decoded_XY_index = " <<  _inverse_hitIndexMapVec[iDetector][index] << 
                        " freq: " <<  _firingFreqVec[ iDetector ][ index ] <<
                        " statusVec[index]: " <<  statusVec[ index ]  << endl;
                }
            }
        }
    }
    
    ++_iEvt;
  } 
  catch (lcio::DataNotAvailableException& e ) 
  {
    streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << endl;
    std::cout <<  "Input collection not found in the current event. Skipping..." << std::endl;
    return;    
  } 
  catch ( ParseException& e ) 
  {
      std::cout <<  e.what() << "\n" << std::endl;
      return;     
  }

}


void EUTelAPIXHotPixelKiller::end() {

  streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
  streamlog_out ( MESSAGE4 ) << printSummary() << endl;

}


string EUTelAPIXHotPixelKiller::printSummary() const {

  stringstream ss;
  int bigSpacer   = 10;
  int tinySpacer  =  5;
  int smallSpacer = 12;
  int fullLineWidth = bigSpacer + tinySpacer + _noOfDetectors * smallSpacer + 1;
  stringstream doubleLine;
  stringstream singleLine;
  for ( int i = 0; i < fullLineWidth ; i++ ) {
    doubleLine << "=";
    singleLine << "-";
  }

  if ( _killedPixelVec.empty() ) {
    return "" ;
  }

  ss << doubleLine.str() << endl
     << " Hot pixel killer summary " << endl
     << doubleLine.str() << endl;


  for ( unsigned int iCycle = 0 ; iCycle < _killedPixelVec[0].size(); iCycle++ ) {
    ss << " " << setiosflags( ios::left ) << setw(bigSpacer) << "Cycle num:"
       << setiosflags( ios::right ) << setw(tinySpacer) << iCycle << resetiosflags( ios::left) ;
    for ( unsigned int iDetector = 0; iDetector < _killedPixelVec.size(); iDetector++ ) {
      ss << setw(smallSpacer) << _killedPixelVec[iDetector][iCycle] ;
    }
    ss << "\n" << singleLine.str() << endl;
  }
  ss << "\n" << doubleLine.str() << endl;

  return ss.str();
}



void EUTelAPIXHotPixelKiller::check( LCEvent * event ) 
{

    if ( _iCycle > static_cast< unsigned short >( _totalNoOfCycle ) )
    {
        return;
    }

    
    if ( _iEvt == _noOfEventPerCycle -1 ) 
    {
        try 
        {            
            LCCollectionVec * statusCollectionVec = dynamic_cast< LCCollectionVec * > ( event->getCollection( _statusCollectionName ) );
            CellIDDecoder<TrackerRawDataImpl>      statusCellDecoder( statusCollectionVec );
            
            for ( unsigned int iDetector = 0; iDetector < _firingFreqVec.size(); iDetector++ ) 
            {
                if ( _iCycle == 0 ) 
                {
                    _killedPixelVec.push_back( vector< unsigned short > (0, 0) );
                }

                TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( statusCollectionVec->getElementAt( iDetector ) );
                vector< short > statusVec = status->adcValues();
                unsigned short killerCounter = 0;
                
                for ( unsigned int iPixel = 0; iPixel < _firingFreqVec[iDetector].size(); iPixel++ ) 
                {
                    if ( _firingFreqVec[iDetector][ iPixel ] / ( static_cast< double >( _iEvt ) ) > _maxAllowedFiringFreq ) 
                    {
                        streamlog_out ( DEBUG5 ) << " Pixel " << iPixel << " on detector " << _sensorIDVec.at( iDetector )
                            << " is firing too often (" << _firingFreqVec[iDetector][iPixel] / (static_cast< double >( _iEvt ) )
                            << "). Masking it now on! " << endl;
                        status->adcValues()[ iPixel ] = EUTELESCOPE::FIRINGPIXEL;
                        ++killerCounter;
                   }
                }
                
                _killedPixelVec[ iDetector ].push_back( killerCounter );
            }


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      bookAndFillHistos();
#endif

      // increment the cycle number
      ++_iCycle;

      if( getBuildHotPixelDatabase() != 0 )
      {
          HotPixelDBWriter(event );
      }

      // reset the _iEvt counter
      _iEvt = 0;

    } catch (lcio::DataNotAvailableException& e ) {
      streamlog_out ( WARNING2 )  << "Input collection not found in the current event. Skipping..." << endl;
      return;
    }
  }
}


void EUTelAPIXHotPixelKiller::HotPixelDBWriter(LCEvent *input_event)
{    
//    lccd::DBInterface dbinterface =  lccd::DBInterface("m26", false);
//    lccd::LCCDTimeStamp timestampe= lccd::LCCDTimeStamp();
//    dbinterface.createSimpleFile 	( timestampe, 	"tag1 ");
//return;

    streamlog_out ( MESSAGE5 ) << "EUTelAPIXHotPixelKiller::HotPixelDBWriter " << endl;
    streamlog_out ( MESSAGE5 ) << "writing out hot pixel db into " << _hotpixelDBFile.c_str() << endl;

    // create data decoder
    LCCollectionVec * statusCollectionVec = dynamic_cast< LCCollectionVec * > ( input_event->getCollection( _statusCollectionName ) );
    CellIDDecoder<TrackerRawDataImpl > decoder( statusCollectionVec );

    // reopen the LCIO file this time in append mode
    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
    LCReader * lcReader = LCFactory::getInstance()->createLCReader();
    LCRunHeaderImpl  *lcHeader = 0;
    LCEventImpl      *event    = 0;
    LCTime           *now      = 0;

    try 
    {
       lcWriter->open( _hotpixelDBFile, LCIO::WRITE_APPEND );
       streamlog_out ( MESSAGE5 ) << _hotpixelDBFile << " was opened for writing" << endl;
    } 
    catch ( IOException& e ) 
    {
      streamlog_out ( ERROR4 ) << e.what() << endl << "Sorry, was not abel to APPEND to a hotpixel file, try open new " << endl;
      try 
      {
          lcWriter->open( _hotpixelDBFile, LCIO::WRITE_NEW );
          // create.write new stuff: 
          // write an almost empty run header
          lcHeader  = new LCRunHeaderImpl;
          lcHeader->setRunNumber( 0 );

          lcWriter->writeRunHeader(lcHeader);
          delete lcHeader;

          event = new LCEventImpl;
          event->setRunNumber( 0 );
          event->setEventNumber( 0 );
          event->setDetectorName("APIX");
          streamlog_out ( MESSAGE5 ) << "event created ok"  << endl;       

          now   = new LCTime;
          event->setTimeStamp( now->timeStamp() );
          delete now;
      }
      catch ( IOException& e ) 
      {
          streamlog_out ( ERROR4 ) << e.what() << endl << "Sorry, was not able to create new file, will quit now. " << endl;
          exit(-1);
      }
    }

    if( event == 0 )
    {
       streamlog_out ( MESSAGE5 ) << "creating new event"  << endl;      
       try
       {
         lcReader->open( _hotpixelDBFile );
         event =  static_cast<LCEventImpl *>  (lcReader->readNextEvent( ));
         streamlog_out ( MESSAGE5 ) << "event read ok"  << endl;
         if( event == 0 )
         {
           lcHeader  = new LCRunHeaderImpl;
           lcHeader->setRunNumber( 0 );
           lcWriter->writeRunHeader(lcHeader);
           delete lcHeader;

           event = new LCEventImpl;
           event->setRunNumber( 0 );
           event->setEventNumber( 0 );
           event->setDetectorName("APIX");
           streamlog_out ( MESSAGE5 ) << "event recreated ok"  << endl;       

           now   = new LCTime;
           event->setTimeStamp( now->timeStamp() );
           delete now; 
         }
      }
      catch(...)
      {
          streamlog_out ( MESSAGE5 ) << "could not read anything"  << endl;       
          exit(-1);
      }
    }
   
    if(event==0)
    {
      streamlog_out ( ERROR5 ) << " event == 0 " << endl;
      return;  
    }

    // create main collection to be saved into the db file 
    LCCollectionVec * hotPixelCollection;

    // create new or open existing hotpixel collection
    try 
    {
        hotPixelCollection = static_cast< LCCollectionVec* > ( event->getCollection( _hotPixelCollectionName  ) );
        streamlog_out ( MESSAGE5 ) << " hotPixelCollection: " << _hotPixelCollectionName << 
                                   " found found with " << hotPixelCollection->getNumberOfElements() << 
                                   " elements " <<  endl; 
        hotPixelCollection->clear();
        streamlog_out ( MESSAGE5 ) << " hotPixelCollection: " << _hotPixelCollectionName << 
                                   " cleared: now " << hotPixelCollection->getNumberOfElements() << 
                                   " elements " <<  endl; 
    }
    catch ( lcio::DataNotAvailableException& e ) 
    {
        hotPixelCollection = new LCCollectionVec( lcio::LCIO::TRACKERDATA );
        event->addCollection( hotPixelCollection, _hotPixelCollectionName );
        streamlog_out ( MESSAGE5 ) << " hotPixelCollection: " << _hotPixelCollectionName << 
                                   " created" <<  endl; 
    }



    for ( unsigned int iDetector = 0; iDetector < _firingFreqVec.size(); iDetector++ ) 
    {
        TrackerRawDataImpl * status = dynamic_cast< TrackerRawDataImpl * > ( statusCollectionVec->getElementAt( iDetector ) );
        int sensorID = decoder( status ) [ "sensorID" ] ;
        streamlog_out ( MESSAGE5 ) <<
                                 " iDetector : " << iDetector  <<
                                 " sensorID : " << sensorID  <<
                                  endl;

        vector< short > statusVec  = status->adcValues();
         
        CellIDEncoder< TrackerDataImpl > hotPixelEncoder  ( eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, hotPixelCollection  );
        hotPixelEncoder["sensorID"]        = sensorID;
        hotPixelEncoder["sparsePixelType"] = eutelescope::kEUTelAPIXSparsePixel;

        // prepare a new TrackerData for the hot Pixel data
        std::auto_ptr<lcio::TrackerDataImpl > currentFrame( new lcio::TrackerDataImpl );
        hotPixelEncoder.setCellID( currentFrame.get() );

        // this is the structure that will host the sparse pixel  
        std::auto_ptr< eutelescope::EUTelSparseDataImpl< eutelescope::EUTelAPIXSparsePixel > >
            sparseFrame( new eutelescope::EUTelSparseDataImpl< eutelescope::EUTelAPIXSparsePixel > ( currentFrame.get() ) );

        for ( unsigned int iPixel = 0; iPixel < _firingFreqVec[iDetector].size(); iPixel++ ) 
        {
            int decoded_XY_index = _inverse_hitIndexMapVec[iDetector][iPixel];

            if ( decoded_XY_index > 0 &&  statusVec[ iPixel ] == EUTELESCOPE::FIRINGPIXEL )                
            {
                 streamlog_out ( MESSAGE5 ) <<
                     " writing out idet: " << iDetector <<
                     " ipixel: " << iPixel <<
                     " statusVec: " << statusVec[ iPixel ]  <<
                     " decoded_XY_index " << decoded_XY_index <<
                     " fired " <<   _firingFreqVec[iDetector][ iPixel ]  / ( static_cast< double >( _iEvt ) ) <<
                     " allowed = " << _maxAllowedFiringFreq << 
                     endl; 
                 sparseFrame->addSparsePixel( _pixelMapVec[iDetector][ decoded_XY_index] );                
            }
        }

        hotPixelCollection->push_back( currentFrame.release() );
    }
 

//    lcReader->close();
// modify to skip if APPENDing an event 
    lcWriter->writeEvent( event );
    
    lcWriter->close();    

    try 
    {
       // open the smae file once again as NEW (overwrite mode!) 
       lcWriter->open( _hotpixelDBFile, LCIO::WRITE_NEW );

       // write out the same event, will be the first event now 
       lcWriter->writeEvent( event );
   
       // close it
       lcWriter->close();    
    }
    catch ( IOException& e ) 
    {
       streamlog_out ( ERROR4 ) << e.what() << endl << "Sorry, was not able to create new file, will quit now. " << endl;
       exit(-1);
    } 

    delete event;
}

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

void EUTelAPIXHotPixelKiller::bookAndFillHistos() 
{

  streamlog_out ( MESSAGE0 ) << "Booking and filling histograms for cycle " << _iCycle << endl;

  string tempHistoName, basePath;
  for ( int iDetector = 0; iDetector < _noOfDetectors; iDetector++ ) 
  { 
    streamlog_out ( MESSAGE0 ) << "-- histograms for iDetector = " << iDetector << " of " << _noOfDetectors << endl;

       
    basePath = "detector_" + to_string( _sensorIDVec.at( iDetector ) ) ;
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());

    basePath += "/cycle_" + to_string( _iCycle );
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath.append("/");


    tempHistoName = _firing2DHistoName + "_d" + to_string( _sensorIDVec.at(iDetector)) + "_c" + to_string( _iCycle ) ;
    int     xBin = _maxX[_sensorIDVec.at( iDetector )] - _minX[_sensorIDVec.at( iDetector )] + 1;
    double  xMin = static_cast<double >(_minX[_sensorIDVec.at( iDetector )]) - 0.5;
    double  xMax = static_cast<double >(_maxX[_sensorIDVec.at( iDetector )]) + 0.5;
    int     yBin = _maxY[_sensorIDVec.at( iDetector )] - _minY[_sensorIDVec.at( iDetector )] + 1;
    double  yMin = static_cast<double >(_minY[_sensorIDVec.at( iDetector )]) - 0.5;
    double  yMax = static_cast<double >(_maxY[_sensorIDVec.at( iDetector )]) + 0.5;


    
    AIDA::IHistogram2D * firing2DHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                xBin, xMin, xMax,yBin, yMin, yMax);

    if( firing2DHisto == 0 )
    {
        streamlog_out ( ERROR5 ) << "CreateHistogram2D for " <<  (basePath + tempHistoName).c_str()  << " failed " << endl;
        streamlog_out ( ERROR5 ) << "Execution stopped, check that your path (" << basePath.c_str() << ")exists  " << endl;
        return;
    }
    
    firing2DHisto->setTitle("Firing frequency map");

    tempHistoName = _firing1DHistoName + "_d" + to_string( _sensorIDVec.at( iDetector) ) + "_c" + to_string( _iCycle );

    int nBin = 100;
    double min = 0;
    double max = 1;
    AIDA::IHistogram1D * firing1DHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
                                                                nBin, min, max );
    firing1DHisto->setTitle("Firing frequency distribution");

    int iPixel = 0;

  //  streamlog_out( MESSAGE5 ) << " iDetector " << iDetector << "  _firingFreqVec.size = " <<  _firingFreqVec.size() << endl;
    
    for (int yPixel = _minY[_sensorIDVec.at( iDetector)]; yPixel <= _maxY[_sensorIDVec.at( iDetector)]; yPixel++) 
    {
      for (int xPixel = _minX[_sensorIDVec.at( iDetector)]; xPixel <= _maxX[_sensorIDVec.at( iDetector)]; xPixel++) 
      {
          //streamlog_out( MESSAGE5 ) << " iPixel = " << iPixel << " _firingFreqVec[ iDetector ].size()= " << _firingFreqVec[ iDetector ].size() <<  endl; 
          if(  iDetector < static_cast< int >(_firingFreqVec.size()) &&  static_cast< int >(_firingFreqVec[ iDetector ].size()) > 0  &&  _firingFreqVec[ iDetector ][ iPixel ] > 0 )
          {
              firing2DHisto->fill(xPixel, yPixel, _firingFreqVec[ iDetector ][ iPixel ] );
              firing1DHisto->fill( _firingFreqVec[ iDetector ][ iPixel ] / ( static_cast< double >( _noOfEventPerCycle ) ) );
              ++iPixel;
          }
      }
    }
  }
  
}



#endif
