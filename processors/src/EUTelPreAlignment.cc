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
#include "EUTelGeometryTelescopeGeoDescription.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
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

EUTelPreAlign::EUTelPreAlign(): Processor("EUTelPreAlign")
{
  _description = "Apply alignment constants to hit collection";

  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName", "The name of the input hit collection", _inputHitCollectionName, std::string("hit"));

  registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedID, 0);

  registerOptionalParameter("AlignmentConstantLCIOFile","Name of LCIO db file where alignment constantds will be stored", _alignmentConstantLCIOFile, std::string("alignment.slcio") );

  registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection that clusters should be checked against (optional).", _hotPixelCollectionName, std::string(""));

  registerProcessorParameter ("Events", "How many events should be used for an approximation to the X,Y shifts (pre-alignment)? (default=50000)", _events, 50000 );
 
  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin, std::vector<float > (6, -10.) );

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax, std::vector<float > (6,  10.) );

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax, std::vector<float > (6,  10.) );

  registerOptionalParameter ("MinNumberOfCorrelatedHits", "If there are more then this number of correlated hits (planes->track candidate) (default=5)",
			     _minNumberOfCorrelatedHits, static_cast <int> (5) );

  registerOptionalParameter("HistogramFilling", "Switch on or off the histogram filling", _fillHistos, bool(true) );
  
  registerOptionalParameter("DumpGEAR", "Dump alignment into GEAR file instead of prealignment database", _dumpGEAR, bool(false) );
  
  registerOptionalParameter("NewGEARSuffix", "Suffix for the new GEAR file, set to empty string (this is not default!) to overwrite old GEAR file", _GEARFileSuffix, std::string("_pre") );

  registerOptionalParameter("ExcludedPlanes", "The list of sensor IDs that shall be excluded.", _ExcludedPlanes, std::vector<int>() );

  registerOptionalParameter("ExcludedPlanesXCoord", "The list of sensor IDs for which the X coordinate shall be excluded.", _ExcludedPlanesXCoord, std::vector<int>() );

  registerOptionalParameter("ExcludedPlanesYCoord", "The list of sensor IDs for which the Y coordinate  shall be excluded.", _ExcludedPlanesYCoord, std::vector<int>() );
}

void EUTelPreAlign::init () {
	// this method is called only once even when the rewind is active
	printParameters ();

	_iRun = 0;  _iEvt = 0;

	_sensorIDVec = geo::gGeometry().sensorIDsVec();
	_sensorIDtoZOrderMap.clear();
	for(size_t index = 0; index < _sensorIDVec.size(); index++) {
		_sensorIDtoZOrderMap.insert( std::make_pair(_sensorIDVec.at(index), (int)index) );
	}

	for( std::vector<int>::iterator it = _sensorIDVec.begin(); it != _sensorIDVec.end(); it++) {
		int sensorID = *it;
		if(sensorID == _fixedID) { 
			_fixedZ = geo::gGeometry().siPlaneZPosition(sensorID); 
		} else {
			_preAligners.push_back( PreAligner(	geo::gGeometry().siPlaneXPitch(sensorID)/10.,
					       			geo::gGeometry().siPlaneYPitch(sensorID)/10.,
								geo::gGeometry().siPlaneZPosition(sensorID),
								sensorID ) );	
		}
	}

	#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	std::string tempHistoName = "";
	std::string basePath; 

	if( _fillHistos ) {
		// Allow any plane to be the fixed reference:
		for(size_t i = 0; i < _sensorIDVec.size(); i++) {
			int sensorID = _sensorIDVec.at(i);

			basePath = "plane_" + to_string( sensorID );
			AIDAProcessor::tree(this)->mkdir(basePath.c_str());
			basePath.append("/");

			tempHistoName = "hitXCorr_fixed_to_" + to_string( sensorID ) ;
			AIDA::IHistogram1D * histo1Da = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.);
			_hitXCorr.insert( make_pair( sensorID, histo1Da) );

			tempHistoName = "hitYCorr_fixed_to_" + to_string( sensorID) ;
			AIDA::IHistogram1D * histo1Db = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 100 , -10., 10.) ;
			_hitYCorr.insert( make_pair( sensorID, histo1Db) );
		}
	}
	#endif
}

void EUTelPreAlign::processRunHeader(LCRunHeader* rdr) {
	std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
    runHeader->addProcessor( type() );
	++_iRun;
}


void  EUTelPreAlign::FillHotPixelMap(LCEvent *event)
{

  if( _hotPixelCollectionName.empty()) return;

  LCCollectionVec* hotPixelCollectionVec = nullptr;
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

      int sensorID = static_cast<int>( cellDecoder( hotPixelData )["sensorID"] );

      if( type  ==  kEUTelGenericSparsePixel )
	{  
	  auto m26Data = std::make_unique<EUTelSparseClusterImpl<EUTelGenericSparsePixel>>(hotPixelData);
	  auto& pixelVec = m26Data->getPixels();

	  for( auto& m26Pixel: pixelVec ) {
              try {
		  _hotPixelMap[sensorID].push_back(std::make_pair(m26Pixel.getXCoord(), m26Pixel.getYCoord()));
		} catch(...) {
		  streamlog_out ( ERROR5 ) << " cannot add pixel to hotpixel map! SensorID: "  << sensorID << ", X:" << m26Pixel.getXCoord() << ", Y:" << m26Pixel.getYCoord() << endl; 
		  abort();
		}
	    }
	}          
    }
}

void EUTelPreAlign::processEvent(LCEvent* event)
{
		if( isFirstEvent()) FillHotPixelMap(event);

		++_iEvt;

		if(_iEvt > _events) return;

		EUTelEventImpl* evt = static_cast<EUTelEventImpl*> (event);

		if(  evt->getEventType() == kEORE ) {
				streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
				return;
		} else if(  evt->getEventType() == kUNKNOWN ) {
				streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
						<< " is of unknown type. Continue considering it as a normal Data Event." << endl;
		}

		try
		{
				LCCollectionVec * inputCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
				UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );

				std::vector<float> residX;
				std::vector<float> residY;
				std::vector<PreAligner*> prealign;

				//Loop over hits in fixed plane:
				for( size_t ref = 0; ref < inputCollectionVec->size(); ref++ )
				{

						TrackerHitImpl* refHit = dynamic_cast<TrackerHitImpl*>( inputCollectionVec->getElementAt(ref) );
						const double* refPos = refHit->getPosition();

						int sensorID = hitDecoder(refHit)["sensorID"];

						// identify fixed plane
						if( sensorID != _fixedID ) continue;

						residX.clear();
						residY.clear();
						prealign.clear();

						for( size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++)
						{

								TrackerHitImpl* hit = dynamic_cast<TrackerHitImpl*>( inputCollectionVec->getElementAt(iHit) );
								//Hits with a hot pixel are ignored
								if( hitContainsHotPixels(hit) ) continue;

								const double * pos = hit->getPosition();
								int iHitID = hitDecoder(hit)["sensorID"]; 

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
										  ) {
												residX.push_back( correlationX );
												residY.push_back( correlationY );
												prealign.push_back(&pa);
										}
										break;
								}
								if( not gotIt ) 
								{
										streamlog_out ( ERROR5 ) << "Mismatched hit at " << pos[2] << endl;
								}
						}

						if( prealign.size() > static_cast< unsigned int >(_minNumberOfCorrelatedHits) && residX.size() == residY.size() ) {
								for( unsigned int ii = 0 ;ii < prealign.size(); ii++ ) {

										prealign[ii]->addPoint( residX[ii], residY[ii] );

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
										if( _fillHistos ) {
												( dynamic_cast<AIDA::IHistogram1D*> (_hitXCorr[ prealign[ii]->getIden() ] ) )->fill( residX[ii] );
												( dynamic_cast<AIDA::IHistogram1D*> (_hitYCorr[ prealign[ii]->getIden() ] ) )->fill( residY[ii] );
										}
#endif
								}
						}
				}
		}
		catch( DataNotAvailableException& e) { 
				streamlog_out  ( WARNING2 ) <<  "No input collection " << _inputHitCollectionName << " found on event " << event->getEventNumber()
						<< " in run " << event->getRunNumber() << endl;
		}

		if( isFirstEvent() ) _isFirstEvent = false;

}

bool EUTelPreAlign::hitContainsHotPixels( TrackerHitImpl   * hit) 
{

  // if no hot pixel map was loaded, just return here
  if( _hotPixelMap.size() == 0) return false;

  try
    {
      LCObjectVec clusterVector = hit->getRawHits();

      if ( hit->getType() == kEUTelSparseClusterImpl ) 
	{
	  TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );

	auto cluster = std::make_unique<EUTelSparseClusterImpl<EUTelGenericSparsePixel>>(clusterFrame);
	  int sensorID = cluster->getDetectorID();
	  auto& pixelVec = cluster->getPixels();

	  for( auto& m26Pixel: pixelVec ) {
	      {
		try{
		  if( std::find(_hotPixelMap.at(sensorID).begin(), 
				_hotPixelMap.at(sensorID).end(),
				std::make_pair(m26Pixel.getXCoord(),m26Pixel.getYCoord()))
		      != _hotPixelMap.at(sensorID).end()){ 
		    return true; // if TRUE  this hit will be skipped
		  }
		}
		catch(const std::out_of_range& oor){
		  streamlog_out(DEBUG0) << " Could not find hot pixel map for sensor ID " 
					<< sensorID << ": " << oor.what() << endl;
		}
	      }
	    }

	  return false;
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
      
void EUTelPreAlign::end()
{
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
				int sensorID = _preAligners.at(ii).getIden();
				std::vector<int>::iterator it = find(_ExcludedPlanes.begin(),_ExcludedPlanes.end(),sensorID);
				std::vector<int>::iterator itXCoord = find(_ExcludedPlanesXCoord.begin(),_ExcludedPlanesXCoord.end(),sensorID);
				std::vector<int>::iterator itYCoord = find(_ExcludedPlanesYCoord.begin(),_ExcludedPlanesYCoord.end(),sensorID);

				EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();

				if(it == _ExcludedPlanes.end())
				{
						if( itXCoord == _ExcludedPlanesXCoord.end() && abs( _preAligners.at(ii).getPeakX() ) < 1000 )
								constant->setXOffset( -1.0* _preAligners.at(ii).getPeakX() );
						else
								constant->setXOffset( 0.0 );

						if(  itYCoord == _ExcludedPlanesYCoord.end() && abs( _preAligners.at(ii).getPeakY() ) < 1000. )
								constant->setYOffset( -1.0 * _preAligners.at(ii).getPeakY() );
						else
								constant->setYOffset( 0.0 );
				}
				else
				{
						constant->setXOffset(0.0);
						constant->setYOffset(0.0);
				}

				constant->setSensorID( sensorID );
				constantsCollection->push_back( constant );

				//Also update the EUTelGeometry descr.
				double updatedXOff = geo::gGeometry().siPlaneXPosition(sensorID) + _preAligners.at(ii).getPeakX();
				double updatedYOff = geo::gGeometry().siPlaneYPosition(sensorID) + _preAligners.at(ii).getPeakY();

				geo::gGeometry().setPlaneXPosition(sensorID, updatedXOff);
				geo::gGeometry().setPlaneYPosition(sensorID, updatedYOff);

				streamlog_out ( MESSAGE5 ) << (*constant) << endl;
		}

		//if we dont dump into the new gear file, we write out the old database
		if(!_dumpGEAR)
		{
				LCWriter* lcWriter = LCFactory::getInstance()->createLCWriter();
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


				streamlog_out( DEBUG5 ) << " adding Collection " << "alignment " << endl;

				event->addCollection( constantsCollection, "alignment" );
				lcWriter->writeEvent( event );
				delete event;
				lcWriter->close();
		}
		else
		{
				//Write updated GEAR file
				marlin::StringParameters* MarlinStringParams = marlin::Global::parameters;
				std::string outputFilename = (MarlinStringParams->getStringVal("GearXMLFile")).substr(0, (MarlinStringParams->getStringVal("GearXMLFile")).size()-4);
				streamlog_out(MESSAGE5) << "Writing updated GEAR file with filename: " << outputFilename+"_pre.xml" << std::endl;
				geo::gGeometry().writeGEARFile(outputFilename+_GEARFileSuffix+".xml");

				//in case we don't write out the collection, we need to delete ourself as we don't pass the collection to LCIO
				delete constantsCollection;
		}
}
