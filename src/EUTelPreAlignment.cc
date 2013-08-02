// Version: $Id$
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelPreAlignment.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
//#include <TrackerHitImpl2.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif


// gear includes <.h>
#include "marlin/Global.h"
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// system includes <>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <memory>
#include <cstdio>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

EUTelPreAlign::EUTelPreAlign () :Processor("EUTelPreAlign") {
  _description = "Apply alignment constants to hit collection";

  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName",
                           "The name of the input hit collection",
                           _inputHitCollectionName, string ("hit"));

  registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedID, 0);

  registerOptionalParameter("AlignmentConstantLCIOFile","Name of LCIO db file where alignment constantds will be stored", 
			    _alignmentConstantLCIOFile, std::string( "alignment.slcio" ) );

  registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection that clusters should be checked against (optional). ",
                             _hotPixelCollectionName, static_cast< string > ( "" ));

  registerProcessorParameter ("Events",
                              "How many events should be used for an approximation to the X,Y shifts (pre-alignment)? (default=50000)",
                              _events, static_cast <int> (50000) );

  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("referenceHit") );
 
  registerOptionalParameter("UseReferenceCollection","Do you want the reference hit collection to be used for coordinate transformations?",  _useReferenceHitCollection, static_cast< bool   > ( true ));
 
  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax, std::vector<float > (6,  10.) );

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax, std::vector<float > (6,  10.) );

  registerOptionalParameter ("MinNumberOfCorrelatedHits",
                              "If there are more then this number of correlated hits (planes->track candidate) (default=5)",
                              _minNumberOfCorrelatedHits, static_cast <int> (5) );

  registerOptionalParameter("HistogramFilling","Switch on or off the histogram filling",_fillHistos, static_cast< bool > ( true ) );

}


void EUTelPreAlign::init () {
  // this method is called only once even when the rewind is active
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;  _iEvt = 0;

  _UsefullHotPixelCollectionFound = 0; 

  _referenceHitVec = 0;

  // clear the sensor ID vector
  _sensorIDVec.clear();

  // clear the sensor ID map
  _sensorIDVecMap.clear();
  _sensorIDtoZOrderMap.clear();


  // clear the sensor ID vector (z-axis order)
  _sensorIDVecZOrder.clear();

//
  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    if(_siPlanesLayerLayout->getID(iPlane) == _fixedID){ 
      //Get Zpos of ref plane
      _fixedZ = _siPlanesLayerLayout->getSensitivePositionZ(iPlane); 
    } else {
      //Get 
      _preAligners.push_back( PreAligner(  _siPlanesLayerLayout->getSensitivePitchX(iPlane),
					   _siPlanesLayerLayout->getSensitivePitchY(iPlane),
					   _siPlanesLayerLayout->getSensitivePositionZ(iPlane),
					   _siPlanesLayerLayout->getID(iPlane)) );
    }
  }


   _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
   for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
   {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
    int sensorID = _siPlanesLayerLayout->getID( iPlane );

    _sensorIDVec.push_back( sensorID );
    _sensorIDVecMap.insert( make_pair( sensorID, iPlane ) );

    // count number of the sensors to the left of the current one:
    int _sensors_to_the_left = 0;
    for ( int jPlane = 0 ; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++ ) 
    {
        if( _siPlanesLayerLayout->getLayerPositionZ(jPlane) + 1e-06 < _siPlaneZPosition[ iPlane ] )
        {
            _sensors_to_the_left++;
        }
    }
    streamlog_out ( DEBUG5 ) << " iPlane " << iPlane << " sensor_#_along_Z_axis " << _sensors_to_the_left << "[z= " << setprecision (3) << _siPlaneZPosition[iPlane] << " ] [sensorID " << sensorID << " ]  " << endl;

    _sensorIDVecZOrder.push_back( _sensors_to_the_left );
    _sensorIDtoZOrderMap.insert(make_pair( sensorID, _sensors_to_the_left));
  }

  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
  {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
    int sensorID = _siPlanesLayerLayout->getID( iPlane );
    _sensorIDinZordered.insert(make_pair( _sensorIDtoZOrderMap[ sensorID ], sensorID ) );
  }

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  string tempHistoName = "";
  string basePath; 

  if(_fillHistos) {
    for(unsigned int i = 1; i < _sensorIDVecZOrder.size(); i++)
    {
      int sensorID = _sensorIDinZordered[i];
      basePath = "plane_" + to_string( sensorID );
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());
      basePath.append("/");
 
      tempHistoName = "hitXCorr_0_to_" + to_string( sensorID ) ;
      AIDA::IHistogram1D * histo1Da = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.);
      _hitXCorr.insert( make_pair( sensorID, histo1Da) );
 
      tempHistoName = "hitYCorr_0_to_" + to_string( sensorID) ;
      AIDA::IHistogram1D * histo1Db = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.) ;
      _hitYCorr.insert( make_pair( sensorID, histo1Db) );
    }
  }
#endif

}

void EUTelPreAlign::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl( rdr ) ) ;
  runHeader->addProcessor( type() );
  ++_iRun;
}


void  EUTelPreAlign::FillHotPixelMap(LCEvent *event)
{

  if (_hotPixelCollectionName.empty()) return;

    LCCollectionVec *hotPixelCollectionVec = 0;
    try 
    {
      hotPixelCollectionVec = static_cast< LCCollectionVec* > ( event->getCollection( _hotPixelCollectionName  ) );
      streamlog_out ( DEBUG5 ) << "Hotpixel database " << _hotPixelCollectionName.c_str() << " found" << endl; 
    }
    catch (...)
    {
      streamlog_out ( WARNING5 ) << "Hotpixel database " << _hotPixelCollectionName.c_str() << " not found" << endl; 
      return;
    }

        CellIDDecoder<TrackerDataImpl> cellDecoder( hotPixelCollectionVec );
	
        for(int i=0; i<  hotPixelCollectionVec->getNumberOfElements(); i++)
        {
           TrackerDataImpl* hotPixelData = dynamic_cast< TrackerDataImpl *> ( hotPixelCollectionVec->getElementAt( i ) );
	   SparsePixelType  type         = static_cast<SparsePixelType> (static_cast<int> (cellDecoder( hotPixelData )["sparsePixelType"]));

	   int sensorID              = static_cast<int > ( cellDecoder( hotPixelData )["sensorID"] );

           if( type  == kEUTelAPIXSparsePixel)
           {  
             auto_ptr<EUTelSparseDataImpl<EUTelAPIXSparsePixel > > apixData(new EUTelSparseDataImpl<EUTelAPIXSparsePixel> ( hotPixelData ));
	     EUTelAPIXSparsePixel apixPixel;
             //Push all single Pixels of one plane in the apixPixelVec
             for ( unsigned int iPixel = 0; iPixel < apixData->size(); iPixel++ ) 
             {
               std::vector<int> apixColVec();
               apixData->getSparsePixelAt( iPixel, &apixPixel);

               try
               {
		 _hotPixelMap[sensorID].push_back(std::make_pair(apixPixel.getXCoord(), apixPixel.getYCoord()));
               }
               catch(...)
               {
		 streamlog_out ( ERROR5 ) << " cannot add pixel to hotpixel map!"  << endl; 
               }
             }
           }
           else if( type  ==  kEUTelSimpleSparsePixel )
           {  
              auto_ptr<EUTelSparseClusterImpl< EUTelSimpleSparsePixel > > m26Data( new EUTelSparseClusterImpl< EUTelSimpleSparsePixel >   ( hotPixelData ) );

              std::vector<EUTelSimpleSparsePixel*> m26PixelVec;
	      EUTelSimpleSparsePixel m26Pixel;
  	      //Push all single Pixels of one plane in the m26PixelVec

             for ( unsigned int iPixel = 0; iPixel < m26Data->size(); iPixel++ ) 
             {
              std::vector<int> m26ColVec();
              m26Data->getSparsePixelAt( iPixel, &m26Pixel);

              try
              {
		 _hotPixelMap[sensorID].push_back(std::make_pair(m26Pixel.getXCoord(), m26Pixel.getYCoord()));
              }
              catch(...)
              {
		streamlog_out ( ERROR5 ) << " cannot add pixel to hotpixel map! SensorID: "  << sensorID << ", X:" << m26Pixel.getXCoord() << ", Y:" << m26Pixel.getYCoord() << endl; 
		abort();
              }
             }
           }          
            else
           {
              _UsefullHotPixelCollectionFound = 0;
           }  	
       }
}

void EUTelPreAlign::processEvent (LCEvent * event) {

  if ( isFirstEvent() )
  {

    FillHotPixelMap(event);

   if ( _useReferenceHitCollection ) 
    {
       try{

   _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
       }
       catch(...)
       {
         _referenceHitVec = 0;
         _useReferenceHitCollection = false;
       }
    }
  }

  ++_iEvt;

  if(_iEvt > _events) return;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }



  try {
    LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));


    std::vector<float>   hitX;
    std::vector<float>   hitY;
    std::vector<PreAligner*>   prealign;


    //Loop over hits in fixed plane
    for (size_t ref = 0; ref < inputCollectionVec->size(); ref++) 
    {

      TrackerHitImpl   * refHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( ref ) ) ;
      const double * refPos = refHit->getPosition();
      // identify fixed plane
      if( guessSensorID(refPos) != _fixedID ) continue;
      hitX.clear();
      hitY.clear();
      prealign.clear();

      for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) 
      {
	TrackerHitImpl   * hit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
        if( hitContainsHotPixels(hit) ) continue;
        
        const double * pos = hit->getPosition();
        int iHitID = guessSensorID(pos);  
        if( iHitID == _fixedID ) continue;
        bool gotIt(false);
	for(size_t ii = 0; ii < _preAligners.size(); ii++)
        {
	  PreAligner& pa = _preAligners.at(ii);
          if( pa.getIden() != iHitID  ) { continue; }
	  gotIt = true;

          double correlationX =  refPos[0] - pos[0] ;
          double correlationY =  refPos[1] - pos[1] ;

          int idZ = _sensorIDtoZOrderMap[ iHitID ];
          if( 
              (_residualsXMin[idZ] < correlationX ) && ( correlationX < _residualsXMax[idZ]) &&
              (_residualsYMin[idZ] < correlationY ) && ( correlationY < _residualsYMax[idZ]) 
            )
          {
            hitX.push_back( correlationX );
            hitY.push_back( correlationY );
            prealign.push_back(&pa);
          }
	  break;
	}
	if(not gotIt) 
        {
	  streamlog_out ( ERROR5 ) << "Mismatched hit at " << pos[2] << endl;
	}
      }
        
      if( prealign.size() > static_cast< unsigned int >(_minNumberOfCorrelatedHits) && hitX.size() == hitY.size() )
      {
         for(unsigned int ii = 0 ;ii < prealign.size();ii++)
         {  
              prealign[ii]->addPoint( hitX[ii], hitY[ii] );
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
              if(_fillHistos) {
                (dynamic_cast<AIDA::IHistogram1D*> (_hitXCorr[prealign[ii]->getIden()]))->fill( hitX[ii] );
                (dynamic_cast<AIDA::IHistogram1D*> (_hitYCorr[prealign[ii]->getIden()]))->fill( hitY[ii] );
              }
#endif
         } 
      }
        
    }
  }
  catch (DataNotAvailableException& e) { 
    streamlog_out  ( WARNING2 ) <<  "No input collection " << _inputHitCollectionName << " found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

  if ( isFirstEvent() ) _isFirstEvent = false;
        
}

bool EUTelPreAlign::hitContainsHotPixels( TrackerHitImpl   * hit) 
{

        bool skipHit = false;

	// if no hot pixel map was loaded, just return here
	if (_hotPixelMap.size() == 0) return 0;
 
        try
        {
            LCObjectVec clusterVector = hit->getRawHits();


           if ( hit->getType() == kEUTelSparseClusterImpl ) 
           {
              TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );
              eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > *cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel >(clusterFrame);

              int sensorID = cluster->getDetectorID();
 
              for ( unsigned int iPixel = 0; iPixel < cluster->size(); iPixel++ ) 
              {
                EUTelSimpleSparsePixel m26Pixel;
                cluster->getSparsePixelAt( iPixel, &m26Pixel);
		{
		  try{
		    if (std::find(_hotPixelMap.at(sensorID).begin(), 
				  _hotPixelMap.at(sensorID).end(),
				  std::make_pair(m26Pixel.getXCoord(),m26Pixel.getYCoord()))
			!= _hotPixelMap.at(sensorID).end()){ 
		      skipHit = true; 	      
		      delete cluster;                        			  
		      return true; // if TRUE  this hit will be skipped
		    }
		  }
		  catch(const std::out_of_range& oor){
		    streamlog_out(DEBUG0) << " Could not find hot pixel map for sensor ID " 
					  << sensorID << ": " << oor.what() << endl;
		  }
		}
              }

              delete cluster;
              return 0;
           } 
            else if ( hit->getType() == kEUTelBrickedClusterImpl ) 
            {

              // fixed cluster implementation. Remember it
              //  can come from
              //  both RAW and ZS data
   
              streamlog_out ( WARNING5 ) << " Hit type kEUTelBrickedClusterImpl is not implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;

           } 
            else if ( hit->getType() == kEUTelDFFClusterImpl ) 
            {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              streamlog_out ( WARNING5 ) << " Hit type kEUTelDFFClusterImpl is not implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;
           } 
            else if ( hit->getType() == kEUTelFFClusterImpl ) 
            {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              streamlog_out ( WARNING5 ) << " Hit type kEUTelFFClusterImpl is not implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;
           } 
            else if ( hit->getType() == kEUTelAPIXClusterImpl ) 
            {
             
                TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );
                eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);
                
                int sensorID = apixCluster->getDetectorID();
                for (unsigned int iPixel = 0; iPixel < apixCluster->size(); ++iPixel) 
                {
                    EUTelAPIXSparsePixel apixPixel;
                    apixCluster->getSparsePixelAt(iPixel, &apixPixel);
                    {
		      try{
			if (std::find(_hotPixelMap.at(sensorID).begin(), 
				      _hotPixelMap.at(sensorID).end(),
				      std::make_pair(apixPixel.getXCoord(), apixPixel.getYCoord()))
			    != _hotPixelMap.at(sensorID).end()){ 
                          skipHit = true; 	      
                          delete apixCluster;                        
                          return true; // if TRUE  this hit will be skipped
			}
			else
			  { 
			    skipHit = false; 	      
			  } 
		      }
		      catch(const std::out_of_range& oor){
			streamlog_out(DEBUG0) << " Could not find hot pixel map for sensor ID " 
					      << sensorID << ": " << oor.what() << endl;
		      }
		      
		    }
                }                
                delete apixCluster;
                return skipHit; // if TRUE  this hit will be skipped
           } 
            else
            {
              streamlog_out ( WARNING5 ) << " Hit type is not known and is not implemented in hotPixel finder method, all pixels are considered for PreAlignment." <<  endl;
           }
 
       }
	catch (exception& e)
	  {
	     streamlog_out(ERROR4) << "something went wrong in EUTelPreAlign::hitContainsHotPixels: " << e.what() << endl;
	  }
	catch(...)
	  { 
	    // if anything went wrong in the above return FALSE, meaning do not skip this hit
	    
	    streamlog_out(ERROR4) << "something went wrong in EUTelPreAlign::hitContainsHotPixels " << endl;
	    return 0;
	  }
	
       // if none of the above worked return FALSE, meaning do not skip this hit
       return 0;
}


int EUTelPreAlign::guessSensorID(const double * hit ) 
{

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;

  if( _referenceHitVec == 0 || _useReferenceHitCollection == false )
  {
    // use z information of planes instead of reference vector
    for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) {
      double distance = std::abs( hit[2] - _siPlaneZPosition[ iPlane ] );
      if ( distance < minDistance ) {
	minDistance = distance;
	sensorID = _siPlanesLayerLayout->getID( iPlane );
      }
    }
    if ( minDistance > 30  ) {
      // advice the user that the guessing wasn't successful 
      streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
	"Please check the consistency of the data with the GEAR file: hitPosition[2]=" << hit[2] <<       endl;
    }
    
    return sensorID;
  }

      for(size_t ii = 0 ; ii < static_cast< unsigned int >(_referenceHitVec->getNumberOfElements()); ii++)
      {
        EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
        
        TVector3 hit3d( hit[0], hit[1], hit[2] );
        TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
 
        double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
        if ( distance < minDistance ) 
        {
           minDistance = distance;
           sensorID = refhit->getSensorID();
        }    

      }

  return sensorID;
}



      
void EUTelPreAlign::end() {
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW    );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
    exit(-1);
  }

  streamlog_out ( MESSAGE5 ) << "Writing to " << _alignmentConstantLCIOFile << endl;

  LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );
  lcWriter->writeRunHeader(lcHeader);
  delete lcHeader;
  LCEventImpl * event = new LCEventImpl;
  event->setRunNumber( 0 );
  event->setEventNumber( 0 );
  LCTime * now = new LCTime;
  event->setTimeStamp( now->timeStamp() );
  delete now;

  LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );


  for(size_t ii=0; ii<_sensorIDVec.size(); ii++)
  {
    bool ifound = false;
    for(size_t jj=0; jj< _preAligners.size(); jj++)
    {
      int sensorID = _preAligners.at(jj).getIden();
      if( _sensorIDVec[ii] == sensorID ) { ifound = true; break; }
    }
    if( ifound == false)
    {
      EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
      constant->setXOffset( 0.0 );
      constant->setYOffset( 0.0 );
      constant->setSensorID( _sensorIDVec[ii] );
      constantsCollection->push_back( constant );
      streamlog_out ( MESSAGE5 ) << (*constant) << endl;
      continue; 
    }
  }


  for(size_t ii = 0 ; ii < _preAligners.size(); ii++){
   EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
   if( abs( _preAligners.at(ii).getPeakX() ) <1000. )
      constant->setXOffset( -1.0 * _preAligners.at(ii).getPeakX());
   else
      constant->setXOffset( -1.0 * 0.0                           );
 
   if( abs( _preAligners.at(ii).getPeakY() ) <1000. )
      constant->setYOffset( -1.0 * _preAligners.at(ii).getPeakY());
   else
      constant->setYOffset( -1.0 * 0.0                           );
 
    int sensorID = _preAligners.at(ii).getIden();
    constant->setSensorID( sensorID );
    constantsCollection->push_back( constant );
    streamlog_out ( MESSAGE5 ) << (*constant) << endl;
  }

  streamlog_out( DEBUG5 ) << " adding Collection " << "alignment " << endl;
 
  event->addCollection( constantsCollection, "alignment" );
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
}
#endif
