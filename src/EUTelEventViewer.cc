/*
 *  Rewritten and alibava features added by Thomas Eichhorn
 *  (2016 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 *
 *  Author Antonio Bulgheroni
 *  (INFN)
 *
 *  email:antonio.bulgheroni@gmail.com
 */

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// this processor is built only if GEAR and MARLINUTIL are available, otherwise
// it is simply skipped
#if defined(USE_GEAR) && defined(USE_MARLINUTIL)

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelEventViewer.h"
#include "EUTelAlignmentConstant.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
# include "marlin/Global.h"

// MarlinUtil
#include "MarlinCED.h"
#include "ClusterShapes.h"
#include "ced_cli.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <vector>
#include <string>
#include <memory>
#include <ctime>
#include <unistd.h>

using namespace std ;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

EUTelEventViewer::EUTelEventViewer() : Processor("EUTelEventViewer"),
	_trackerHitCollectionNameVec(),
	_trackCollectionNameVec(),
	_alignmentCollectionName(""),
	_layerTrackerHit(-1),
	_layerTrack(-1),
	_waitForKeyboard(true),
	_autoForwardDelay(4.0),
	_detModel(99999)
{

	_description = "Event display" ;

	vector<string > trackerHitCollectionNameVecExample;
	trackerHitCollectionNameVecExample.push_back("hit");

	registerInputCollections( LCIO::TRACKERHIT, "TrackerHitCollections","Tracker hit collection names",  _trackerHitCollectionNameVec,trackerHitCollectionNameVecExample );

	vector<string > trackCollectionNameVecExample;
	trackCollectionNameVecExample.push_back("testfittrack");

	registerInputCollections( LCIO::TRACK, "TrackCollections","Track collection names", _trackCollectionNameVec,trackCollectionNameVecExample );

	registerProcessorParameter("LayerTrackerHit","Layer for Tracker Hits",_layerTrackerHit,static_cast< int >(-1));

	registerProcessorParameter("LayerTracks","Layer for Tracks",_layerTrack,static_cast< int >(-1));

	registerProcessorParameter("DetectorModel","Detector Model: 99999 to draw the ideal model taken from the GEAR description (with no rotations), -1 to draw the model described by GEAR and corrected for alignment, -2 to draw the model from GEAR with rotations.",_detModel,static_cast< int >(99999));

	registerProcessorParameter("WaitForKeyboard","Wait for keyboard input before proceeding to next event",_waitForKeyboard,static_cast< bool > (true) );

	registerProcessorParameter("AutoForwardDelay","This is the time in seconds between two following event displays. Enabled only when WaitForKeybord is off", _autoForwardDelay, 4.0);

	registerOptionalParameter("AlignmentConstantName","Alignment constants from the condition file, used in detector model -1",_alignmentCollectionName, string ("alignment"));

	registerOptionalParameter("HitSize","Size of drawn hits",_hitsize,int(3));

	registerOptionalParameter("LineSize","Size of drawn hits",_linesize,int(3));

	registerOptionalParameter("SensorThickness","Additional drawn sensor thickness for visibility in detector models -1 and -2",_sensorthickness,double(1.5));

	registerOptionalParameter("XScale","Scaling factor for all x dimensions",_xscale,double(1.0));

	registerOptionalParameter("YScale","Scaling factor for all y dimensions",_yscale,double(1.0));

	registerOptionalParameter("ZScale","Scaling factor for all z dimensions",_zscale,double(1.0));

	registerOptionalParameter("DUTHighlight","To highlight measured DUT hits, set it's sensor id here. -1 to switch this off",_duthighlight,int(-1));

}


void EUTelEventViewer::init()
{
	printParameters();
	MarlinCED::init(this) ;
}


void EUTelEventViewer::processRunHeader( LCRunHeader * rdr )
{
	std::unique_ptr<EUTelRunHeaderImpl> runHeader = std::make_unique<EUTelRunHeaderImpl>(rdr);
	runHeader->addProcessor(type());
}


void EUTelEventViewer::processEvent( LCEvent * evt )
{

	EUTelEventImpl * event = static_cast<EUTelEventImpl *> ( evt );

	if ( event->getEventType() == kEORE )
	{
		streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
		return;
	} else if ( event->getEventType() == kUNKNOWN ) {
		streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber() << " is of unknown type." << endl;
		streamlog_out ( WARNING2 ) << "Continue considering it as a normal data event." << endl;
	}

	// for convenience
	int tempDetModel = _detModel;

	// try to load the alignment collections
	LCCollectionVec * alignmentCollectionVec = NULL;
	try
	{
		alignmentCollectionVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_alignmentCollectionName));
	} catch ( DataNotAvailableException& e ) {

		if ( _detModel == -1 )
		{

			// the user wants to draw the alignment corrected geometry, but
			// no alignment collection found
			streamlog_out ( WARNING2 ) << "No alignment condition found, using default ideal geometry model (99999)" << endl;
			streamlog_out ( WARNING2 ) << "Scaling will not work, reverting to 1!" << endl;
			tempDetModel = 99999;
			_xscale = 1.0;
			_yscale = 1.0;
			_zscale = 1.0;
		}
	}

	//check
	if (tempDetModel == 99999)
	{
		if ((_xscale != 1.0) || (_yscale != 1.0) || (_zscale != 1.0) )
		{
			streamlog_out ( WARNING2 ) << "Scaling will not work with geometry model 99999, reverting scaling to 1!" << endl;
			_xscale = 1.0;
			_yscale = 1.0;
			_zscale = 1.0;
		}
	}

	// Drawing Geometry
	MarlinCED::newEvent( this, tempDetModel ) ;

	// draw according to gear file and correct for alignment
	if ( tempDetModel == -1 )
	{

		streamlog_out(DEBUG0) << "Using detector model " << tempDetModel << endl;
		const gear::SiPlanesParameters&  siPlanesParameters  = Global::GEAR->getSiPlanesParameters();
		const gear::SiPlanesLayerLayout& siPlanesLayerLayout = siPlanesParameters.getSiPlanesLayerLayout();

		double * rotation = new double[3];
		double * sizes  = new double[3];
		double * center = new double[3];
		unsigned int color = 0xFFFFFF;

		for ( int iLayer = 0 ; iLayer < siPlanesLayerLayout.getNLayers() ; iLayer++ )
		{
			streamlog_out(DEBUG0) << "Drawing layer " << iLayer << " of " << siPlanesLayerLayout.getNLayers() << endl;

			EUTelAlignmentConstant * alignment = static_cast< EUTelAlignmentConstant * > ( alignmentCollectionVec->getElementAt( iLayer ) );

			center[0] = (siPlanesLayerLayout.getSensitivePositionX(iLayer) - alignment->getXOffset())*_xscale;
			center[1] = (siPlanesLayerLayout.getSensitivePositionY(iLayer) - alignment->getYOffset())*_yscale;
			center[2] = (siPlanesLayerLayout.getSensitivePositionZ(iLayer) - alignment->getZOffset())*_zscale;

			sizes[0]  = (siPlanesLayerLayout.getSensitiveSizeX(iLayer)*_xscale);
			sizes[1]  = (siPlanesLayerLayout.getSensitiveSizeY(iLayer)*_yscale);
			sizes[2]  = (siPlanesLayerLayout.getSensitiveThickness(iLayer)+_sensorthickness)*_zscale ;

			rotation[0] = -1.0*siPlanesLayerLayout.getLayerRotationZY(iLayer)- alignment->getAlpha();
			rotation[1] = -1.0*siPlanesLayerLayout.getLayerRotationZX(iLayer)- alignment->getAlpha();
			rotation[2] = -1.0*siPlanesLayerLayout.getLayerRotationXY(iLayer)- alignment->getAlpha();

			ced_geobox_r( sizes, center, rotation, color, 1);
		}
		delete [] center;
		delete [] sizes;
	}

	// draw according to gear file only
	if ( tempDetModel == -2 )
	{

		streamlog_out(DEBUG0) << "Using detector model " << tempDetModel << endl;
		const gear::SiPlanesParameters&  siPlanesParameters  = Global::GEAR->getSiPlanesParameters();
		const gear::SiPlanesLayerLayout& siPlanesLayerLayout = siPlanesParameters.getSiPlanesLayerLayout();

		double * rotation = new double[3];
		double * sizes  = new double[3];
		double * center = new double[3];
		unsigned int color = 0xFFFFFF;

		for ( int iLayer = 0 ; iLayer < siPlanesLayerLayout.getNLayers() ; iLayer++ )
		{
			streamlog_out(DEBUG0) << "Drawing layer " << iLayer << " of " << siPlanesLayerLayout.getNLayers() << endl;

			center[0] = siPlanesLayerLayout.getSensitivePositionX(iLayer)*_xscale;
			center[1] = siPlanesLayerLayout.getSensitivePositionY(iLayer)*_yscale;
			center[2] = siPlanesLayerLayout.getSensitivePositionZ(iLayer)*_zscale;

			sizes[0]  = siPlanesLayerLayout.getSensitiveSizeX(iLayer)*_xscale;
			sizes[1]  = siPlanesLayerLayout.getSensitiveSizeY(iLayer)*_yscale;
			sizes[2]  = (siPlanesLayerLayout.getSensitiveThickness(iLayer)+_sensorthickness)*_zscale ;

			rotation[0] = siPlanesLayerLayout.getLayerRotationZY(iLayer);
			rotation[1] = siPlanesLayerLayout.getLayerRotationZX(iLayer);
			rotation[2] = siPlanesLayerLayout.getLayerRotationXY(iLayer);

			ced_geobox_r( sizes, center, rotation, color, 1);

		}
		delete [] center;
		delete [] sizes;
	}

	// Drawing hit collections
	if ( _layerTrackerHit >= 0 )
	{
		for ( unsigned int iCollection = 0; iCollection < _trackerHitCollectionNameVec.size(); iCollection++ )
		{
			try
			{
				LCCollection * collection = evt->getCollection( _trackerHitCollectionNameVec[iCollection].c_str() );
				streamlog_out(DEBUG0) << "Drawing hit collection " << _trackerHitCollectionNameVec[iCollection].c_str() << endl;

				for ( int iHit = 0; iHit < collection->getNumberOfElements(); iHit++ )
				{
					streamlog_out(DEBUG0) << "Drawing hit " << iHit << " of " << collection->getNumberOfElements() << endl;
					TrackerHitImpl * hit = dynamic_cast<TrackerHitImpl *> ( collection->getElementAt(iHit) ) ;
					float x = static_cast<float > ( hit->getPosition()[0]*_xscale );
					float y = static_cast<float > ( hit->getPosition()[1]*_yscale );
					float z = static_cast<float > ( hit->getPosition()[2]*_zscale );
					unsigned int color =  returnColor(iCollection+1);
					unsigned int type  = CED_HIT_STAR;
					int size = _hitsize;

					UTIL::CellIDDecoder<TrackerHitImpl> hitDecoder ( EUTELESCOPE::HITENCODING );
					int sensorid = hitDecoder(hit)["sensorID"];
					if (sensorid == _duthighlight)
					{
						streamlog_out(DEBUG1) << "Highlighting DUT hit!" << endl;
						color = returnColor(iCollection+2);
						size  = 2*_hitsize;
					}

					ced_hit(x,y,z,type ,size , color);
				}
			}  catch (DataNotAvailableException& e) {
				streamlog_out ( WARNING2 ) << "No input collection (" << _trackerHitCollectionNameVec[iCollection] << " found on " <<  event->getEventNumber() << " in run " << event->getRunNumber() << endl;
			}
		}
	}

	// setup cellIdDecoder to decode the hit properties
	CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

	// Drawing  Tracks
	if (_layerTrack >= 0)
	{
		for ( unsigned int iCollection = 0; iCollection < _trackCollectionNameVec.size(); iCollection++ )
		{
			try
			{
				LCCollection * collection =    evt->getCollection(_trackCollectionNameVec[iCollection].c_str());
				streamlog_out(DEBUG0) << "Drawing track collection " << _trackCollectionNameVec[iCollection].c_str() << endl;

				for ( int iTrack = 0; iTrack < collection->getNumberOfElements(); iTrack++ )
				{

					streamlog_out(DEBUG0) << "Drawing track " << iTrack << " of " << collection->getNumberOfElements() << endl;
					TrackImpl * track = dynamic_cast<TrackImpl*> ( collection->getElementAt(iTrack) );
					TrackerHitVec hitvec = track->getTrackerHits();

					unsigned int color =  returnColor(iCollection);

					// ok now I have everything I'm interested in. The idea
					// behing is that for each fitted hit I will draw a hit and
					// a line to the next plane.
					float xPrev = 0, yPrev = 0, zPrev = 0;
					float xNext, yNext, zNext;

					bool first = true;

					for ( size_t iHit = 0; iHit < hitvec.size()  ; ++iHit )
					{

						streamlog_out(DEBUG0) << "Drawing track hit " << iHit << " of " << hitvec.size() << endl;
						TrackerHitImpl * hit;
						if ( (hit = dynamic_cast<TrackerHitImpl*> ( hitvec[ iHit ] )) != 0x0 )
						{

							// show only hit resulting from fit
							if ( (hitCellDecoder(hit)["properties"] & kFittedHit) > 0 )
							{
								if ( first )
								{
									first = false;
									xPrev = static_cast<float> ( hit->getPosition()[0]*_xscale );
									yPrev = static_cast<float> ( hit->getPosition()[1]*_yscale );
									zPrev = static_cast<float> ( hit->getPosition()[2]*_zscale );
									ced_hit(xPrev, yPrev, zPrev, CED_HIT_CROSS,_hitsize,color);
								} else {
									xNext = static_cast<float> ( hit->getPosition()[0]*_xscale );
									yNext = static_cast<float> ( hit->getPosition()[1]*_yscale );
									zNext = static_cast<float> ( hit->getPosition()[2]*_zscale );
									ced_hit(xNext, yNext, zNext, CED_HIT_CROSS,_hitsize,color);

									ced_line( xPrev, yPrev, zPrev, xNext, yNext, zNext, _layerTrack << CED_LAYER_SHIFT, _linesize, color );

									xPrev = xNext;
									yPrev = yNext;
									zPrev = zNext;

								}
							}
						}
					}
				}
			}  catch(DataNotAvailableException &e) {
				streamlog_out ( WARNING2 ) << "No input collection (" << _trackerHitCollectionNameVec[iCollection] << " found on " <<  event->getEventNumber() << " in run " << event->getRunNumber() << endl;
			}
		}
	}

	// draw the current event!
	MarlinCED::draw( this, _waitForKeyboard ) ;

	// in case we don't have to wait for the kepyboard, so wait
	// _autoForwardDelay second before continue
	if ( ! _waitForKeyboard ){
		usleep( static_cast< int > (_autoForwardDelay * 1000000) );
	}
}


void EUTelEventViewer::check( LCEvent *)
{

}


void EUTelEventViewer::end()
{

}


int EUTelEventViewer::returnColor(int counter)
{

	int icol =  counter % 16;
	int kcol =  0x00ff00;
	if (icol==1)  kcol = 0xAA00ff;
	if (icol==2)  kcol = 0xff0000;
	if (icol==3)  kcol = 0x00ffff;
	if (icol==4)  kcol = 0xffff00;
	if (icol==5)  kcol = 0xff00ff;
	if (icol==6)  kcol = 0xffffff;
	if (icol==7)  kcol = 0x0fffff;
	if (icol==8)  kcol = 0x000088;
	if (icol==9)  kcol = 0x008800;
	if (icol==10) kcol = 0x880000;
	if (icol==11) kcol = 0x008888;
	if (icol==12) kcol = 0x888800;
	if (icol==13) kcol = 0x880088;
	if (icol==14) kcol = 0x888888;
	if (icol==15) kcol = 0x00A888;

	return kcol;
}


#endif // USE_GEAR && USE_CED
