/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IO/LCWriter.h>
#include <IO/LCReader.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <UTIL/LCTOOLS.h>

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <glob.h>
#include <vector>
#include <set>
#include <map>

// eutelescope includes ""
#include "EUTELESCOPE.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "anyoption.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"

#include "CBCClustering.h"


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace IMPL;
using namespace eutelescope;


CBCClustering::CBCClustering ( ) : Processor ( "CBCClustering" ),
_aidaHistoMap ( )
{

    _description = "CBCClustering clusters the CBC data stream.";

    registerProcessorParameter ( "CBCInputCollectionName", "The name of the CBC collection we want to read", _cbcInputCollectionName, string ( "cbc_input" ) );

    registerProcessorParameter ( "CBCDataOutputCollectionName", "The name of the CBC data collection we want to write", _cbcDataOutputCollectionName, string ( "cbc_data_output" ) );

    registerProcessorParameter ( "CBCPulseOutputCollectionName", "The name of the CBC pulse collection we want to write", _cbcPulseOutputCollectionName, string ( "cbc_pulse_output" ) );

    registerProcessorParameter ( "ChannelCount", "The total number of channels in the sensor", _chancount, 1016 );

    registerProcessorParameter ( "MaxClusterCountPerEvent", "The maximum allowed number of clusters in an event, events with more clusters will be discarded", _maxclusters, 4 );

    registerProcessorParameter ( "MaxClusterSize", "The maximum allowed cluster size, larger clusters will be discarded as noise", _maxclustersize, 4 );

    registerProcessorParameter ( "NonSensitiveAxis", "The unsensitive axis of the CBC", _nonsensitiveaxis, string ( "x" ) );

    registerProcessorParameter ( "OutputSensorID", "The sensor id to write", _outputSensorID, 6 );

    registerProcessorParameter ( "ZSMode", "Zero Suppression Mode? 0 for off, 1 for first, 2 for second sensor", _zsmode, 0 );

}


void CBCClustering::init ( )
{
    streamlog_out ( MESSAGE4 ) << "Running init" << endl;
    printParameters ( );

    if ( _nonsensitiveaxis == "x" || _nonsensitiveaxis == "X" )
    {
	_nonsensitiveaxis = "x";
	streamlog_out ( MESSAGE4 ) << "Non-sensitive axis is x!" << endl;
    }
    else if ( _nonsensitiveaxis == "y" || _nonsensitiveaxis == "Y" )
    {
	_nonsensitiveaxis = "y";
	streamlog_out ( MESSAGE4 ) << "Non-sensitive axis is y!" << endl;
    }
    else
    {
	streamlog_out ( ERROR5 ) << "Illegal setting for NonSensitiveAxis! Valid coordinates are x and y!" << endl;
	exit ( -1 );
    }

}


void CBCClustering::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    auto arunHeader = std::make_unique < EUTelRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    bookHistos ( );

}


void CBCClustering::processEvent ( LCEvent * anEvent )
{

    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

	// the collection we read
	LCCollectionVec * inputCollectionVec;

	// the collection we output
	LCCollectionVec * clusterCollection;

	// the sparse collection we output
	LCCollectionVec * sparseClusterCollectionVec = nullptr;
	sparseClusterCollectionVec =  new LCCollectionVec ( LCIO::TRACKERDATA );

	try
	{
	    clusterCollection = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcDataOutputCollectionName ) );
	}
	catch ( lcio::DataNotAvailableException& e )
	{
	    clusterCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
	}

	// find the seed clusters on our data
	try
        {
		// give the collection vec its data
		inputCollectionVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcInputCollectionName ) );

		// loop over collection sizes, just in case
		int noOfDetector = inputCollectionVec -> getNumberOfElements ( );
                for ( int i = 0; i < noOfDetector; ++i ) 
                {
		    TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( inputCollectionVec -> getElementAt ( i ) );

		    FloatVec datavec;
		    datavec = trkdata -> getChargeValues ( );

		    /*
		    if ( _zsmode != 0 && ( _zsmode != ( i + 1 ) ) )
		    {
			size_t temp = datavec.size ( );
			datavec.assign ( temp, 0 );
		    }
		    */

		    int nClusters = 0;

		    //FloatVec outputvec;

		    /*
		    // highly inefficient way...
		    if ( _zsmode == ( i + 1 ) )
		    {

			FloatVec pixel_x;
			FloatVec pixel_y;
			FloatVec pixel_charge;
			FloatVec pixel_time;

			if ( datavec.size ( ) > 0 )
			{
			    // data is x,y,q,t
			    for ( size_t ix = 0; ix <= ( datavec.size ( ) - 4 ); ix = ix + 4 )
			    {
				streamlog_out ( DEBUG1 ) << "chan x " << datavec[ix] << endl;
				streamlog_out ( DEBUG1 ) << "chan y " << datavec[ix+1] << endl;
				streamlog_out ( DEBUG1 ) << "chan q " << datavec[ix+2] << endl;
				streamlog_out ( DEBUG1 ) << "chan t " << datavec[ix+3] << endl;
				pixel_x.push_back ( datavec[ix] );
				pixel_y.push_back ( datavec[ix+1] );
				pixel_charge.push_back ( datavec[ix+2] );
				pixel_time.push_back ( datavec[ix+3] );
			    }
			}

			// fill the channels with a signal if it's there...
			int point = 0;

			for ( int j = 0; j < _chancount; j++ )
			{

			    double background = 0.0;

			    // add the charge, if any
			    if ( _nonsensitiveaxis == "x" )
			    {
				if ( pixel_y.size ( ) > 0 )
				{
				    if ( pixel_y[point] == ( j ) )
				    {
					double signal = 0.0;
					signal = pixel_charge[point];
					background += signal;
					point++;
					streamlog_out ( DEBUG2 ) << "Output chanel "<< j << " signal is " << signal << endl;
				    }
				}
			    }
			    if ( _nonsensitiveaxis == "y" )
			    {
				if ( pixel_x.size ( ) > 0 )
				{
				    if ( pixel_x[point] == ( j ) )
				    {
					double signal = 0.0;
					signal = pixel_charge[point];
					background += signal;
					point++;
					streamlog_out ( DEBUG2 ) << "Output chanel "<< j << " signal is " << signal << endl;
				    }
				}
			    }
			    // push back this channel
			    outputvec.push_back(background);
			    if ( background > 0.0 )
			    {
				streamlog_out ( DEBUG3 ) << "Output chanel "<< j << " final signal is " << background << endl;
			    }
			}

			datavec = outputvec;

		    }
		    */

		    std::vector < int > clusterNumber;
		    
		    
		    
		    if ( _zsmode == 0 )
		    {
			clusterNumber.assign ( datavec.size ( ), 0 );

			for ( size_t ichan = 0; ichan < datavec.size ( ); ichan++ )
			{
			    int cluleft = 0;
			    int cluright = 0;
			    int clusize = 0;

			    if ( datavec[ichan] * 1.0 > 0 )
			    {
				streamlog_out ( DEBUG0 ) << "Chan " << ichan << " over threshold" << endl;
				cluleft = ichan;
				cluright = ichan;
				bool search = true;
				nClusters++;
				clusize++;
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["Hitmap_" + to_string ( _outputSensorID + i ) ] ) -> fill ( ichan );
				while ( search == true )
				{
				    if ( ( ichan + 1 ) < ( datavec.size ( ) - 1 ) )
				    {
					ichan++;
					if ( datavec[ichan] > 0 )
					{
					    if ( clusize < _maxclustersize )
					    {
						cluright = ichan;
						clusize++;
						dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["Hitmap_" + to_string ( _outputSensorID + i ) ] ) -> fill ( ichan );
						streamlog_out ( DEBUG1 ) << "Chan " << ichan << " added to cluster, size now " << clusize << endl;
					    }
					    else
					    {
						streamlog_out ( DEBUG4 ) << "Cluster size larger than allowed max cluster size of " << _maxclustersize << "! Discarding cluster in event " << anEvent -> getEventNumber ( ) << "!" << endl;
						search = false;
					    }
					}
					else
					{
					    search = false;
					}
				    }
				    else
				    {
					search = false;
				    }
				}
				for ( int k = cluleft; k <= cluright; k++ )
				{
				    clusterNumber.at ( k ) = nClusters;
				}
			    }
			}

		    }

		    if ( _zsmode > 0 )
		    {
			streamlog_out ( DEBUG4 ) << "Using ZS mode, sensor " << i << endl;
			clusterNumber.assign ( _chancount, 0 );
			int clupos = -2;
			int temppos = -2;
			int clusize = 0;
			if ( datavec.size ( ) > 0 )
			{
			    if ( _nonsensitiveaxis == "y" )
			    {
				for ( size_t ix = 0; ix <= ( datavec.size ( ) - 4 ); ix = ix + 4 )
				{
				    clupos = datavec[ix];
				    if ( ( ( clupos - temppos ) == 1 ) && clusize < _maxclustersize )
				    {
					clusterNumber.at ( clupos ) = nClusters;
					temppos = clupos;
					clusize++;
					streamlog_out ( DEBUG4 ) << "Enlarged cluster at " << clupos << endl;
					dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["Hitmap_" + to_string ( _outputSensorID + i ) ] ) -> fill ( clupos );
				    }
				    else
				    {
					nClusters++;
					clusterNumber.at ( clupos ) = nClusters;
					streamlog_out ( DEBUG4 ) << "Found new cluster at " << clupos << endl;
					temppos = clupos;
					clusize = 1;
					dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["Hitmap_" + to_string ( _outputSensorID + i ) ] ) -> fill ( clupos );
				    }
				}
			    }
			    if ( _nonsensitiveaxis == "x" )
			    {
				for ( size_t ix = 0; ix <= ( datavec.size ( ) - 4 ); ix = ix + 4 )
				{
				    clupos = datavec[ix + 1];
				    if ( ( ( clupos - temppos ) == 1 ) && clusize < _maxclustersize )
				    {
					clusterNumber.at ( clupos ) = nClusters;
					temppos = clupos;
					clusize++;
				    }
				    else
				    {
					nClusters++;
					clusterNumber.at ( clupos ) = nClusters;
					temppos = clupos;
					clusize = 1;
				    }
				}
			    }
			}
			datavec.clear ( );
			for ( size_t ii = 0; ii < clusterNumber.size ( ); ii++ )
			{
			    datavec.push_back ( clusterNumber.at ( ii ) );
			}
			
		    }

		    if ( nClusters > _maxclusters )
		    {
			clusterNumber.assign ( _chancount, 0 );
			streamlog_out ( DEBUG4 ) << "Found " << nClusters << " clusters in event " << anEvent -> getEventNumber ( ) << "! Discarding all of them!" << endl;
		    }

		    // now output the clusters we have found
		    CellIDEncoder < TrackerPulseImpl > zsDataEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, clusterCollection );
		    CellIDEncoder < TrackerDataImpl > idClusterEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );

		    // decode datavec into pixelvec, each chan goes into a pixel
		    //std::vector < eutelescope::EUTelGenericSparsePixel* > PixelVec;
		    std::vector<EUTelGenericSparsePixel> PixelVec;
		    eutelescope::EUTelGenericSparsePixel Pixel;
		    std::unique_ptr < EUTelTrackerDataInterfacerImpl < EUTelGenericSparsePixel > > pixelData = std::make_unique < EUTelTrackerDataInterfacerImpl < EUTelGenericSparsePixel > > ( trkdata );
		    for ( unsigned int iPixel = 0; iPixel < datavec.size ( ); iPixel++ )
		    {
			if ( _nonsensitiveaxis == "x" )
			{
			    Pixel.setXCoord ( 0.0 );
			    Pixel.setYCoord ( iPixel );
			}
			else
			{
			    Pixel.setXCoord ( iPixel );
			    Pixel.setYCoord ( 0.0 );
			}
			Pixel.setSignal ( datavec[iPixel] * 1.0 );
			PixelVec.push_back ( Pixel );
			streamlog_out ( DEBUG0 ) << "Decoded pixel with coordinates: ( " << Pixel.getXCoord ( ) << " , " << Pixel.getYCoord ( ) << " ) and charge: " << datavec[iPixel] << " !" << endl;
			streamlog_out ( DEBUG0 ) << "Input was: " << PixelVec.at ( iPixel ) << endl;
		    }

		    std::vector < int > ::iterator it;
		    for ( it = clusterNumber.begin ( ); it != clusterNumber.end ( ); ++it )
		    {
			lcio::TrackerPulseImpl * pulseFrame = new lcio::TrackerPulseImpl ( );
			lcio::TrackerDataImpl * clusterFrame = new lcio::TrackerDataImpl ( );
			eutelescope::EUTelSparseClusterImpl < eutelescope::EUTelGenericSparsePixel > *pixelCluster = new eutelescope::EUTelSparseClusterImpl < eutelescope::EUTelGenericSparsePixel > ( clusterFrame );
			
			for ( size_t j = 0; j < datavec.size ( ); j++ )
			{
			    if ( clusterNumber.at ( j ) == *it ) 
			    {
				if ( clusterNumber.at ( j ) > 0 )
				{
				    // put only these pixels in that ClusterCollection that belong to that cluster
				    auto temppix = PixelVec.at(i);
				    pixelCluster -> push_back ( temppix );
				    streamlog_out ( DEBUG1 ) << "Evt " << anEvent -> getEventNumber ( )<<" Adding channel " << j << " to cluster " << clusterNumber.at ( j ) << endl;
				    if ( it != clusterNumber.end ( ) - 1 )
				    {
					++it;
				    }
				}
			    }
			}

			// now we have pixelClusters with the corresponding pixels (which have the xy and q info) in them
			// this if stops making the "zero" cluster with all strips not in a cluster...
			unsigned int clumin = 1;
			unsigned int clumax = 99;
			if ( ( pixelCluster -> size ( ) >= clumin ) && ( pixelCluster -> size ( ) <= clumax ) && ( pixelCluster -> getTotalCharge ( ) >= 1 ) ) 
			{

			    // make a frame of each cluster and give it position, charge, etc. Then push back into clusterCollection and sparseClusterCollectionVec.
			    float x, y = 0;
			    int xsize, ysize = 0;
			    float charge = 0.0;
			    charge = pixelCluster -> getTotalCharge ( );
			    pulseFrame -> setCharge ( charge );
			    pixelCluster -> getCenterOfGravity ( x, y );
			    pixelCluster -> getClusterSize ( xsize, ysize );

			    streamlog_out( DEBUG1 ) << "Cluster: " << *it << ", Q: " << charge << " , x: " << x << " , y: " << y << " , dx: " << xsize << " , dy: " << ysize << " in event: " << anEvent -> getEventNumber ( ) << endl;

			    dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ClusterCharge_" + to_string ( _outputSensorID + i ) ] ) -> fill ( charge );
			    if ( _nonsensitiveaxis == "x" )
			    {
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ClusterSize_" + to_string ( _outputSensorID + i ) ] ) -> fill ( ysize );
			    }
			    if ( _nonsensitiveaxis == "y" )
			    {
				dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ClusterSize_" + to_string ( _outputSensorID + i ) ] ) -> fill ( xsize );
			    }

			    zsDataEncoder["sensorID"] = _outputSensorID + i;
			    zsDataEncoder["xSeed"] = static_cast < long > ( x );
			    zsDataEncoder["ySeed"] = static_cast < long > ( y );
			    zsDataEncoder["xCluSize"] = xsize;
			    zsDataEncoder["yCluSize"] = ysize;
			    zsDataEncoder["type"] = static_cast < int > ( kEUTelSparseClusterImpl );
			    zsDataEncoder["quality"] =  0;
			    zsDataEncoder.setCellID ( pulseFrame );
			    pulseFrame -> setTrackerData ( clusterFrame );
			    clusterCollection -> push_back ( pulseFrame );

			    idClusterEncoder["sensorID"] = _outputSensorID + i;
			    idClusterEncoder["sparsePixelType"] = 2;
			    idClusterEncoder["quality"] = 0;
			    idClusterEncoder.setCellID ( clusterFrame );
			    sparseClusterCollectionVec -> push_back ( clusterFrame );
			}
			else // clean up
			{
			    delete pulseFrame;
			    delete clusterFrame;
			}

			delete pixelCluster;

		    } // done clusterNumber iteration

		}
	    }
	    catch ( lcio::DataNotAvailableException& )
	    {

	    }

	anEvent->addCollection( sparseClusterCollectionVec, _cbcDataOutputCollectionName );
	anEvent->addCollection( clusterCollection, _cbcPulseOutputCollectionName );

}


void CBCClustering::check ( LCEvent * /* evt */ )
{

}


void CBCClustering::end ( )
{
    streamlog_out ( MESSAGE4 ) << "Successfully finished!" << endl;

}


void CBCClustering::fillHistos ( )
{

}


void CBCClustering::bookHistos ( )
{
    for ( int i = 0; i < 2; i++ )
    {
	string basePath = "Clustering_" + to_string ( _outputSensorID + i);
	AIDAProcessor::tree ( this ) -> mkdir ( basePath.c_str ( ) );
	basePath.append ( "/" );

	AIDA::IHistogram1D * csizeHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "ClusterSize_" + to_string ( _outputSensorID + i) ).c_str ( ), 10, -0.5, 9.5 );
	_aidaHistoMap.insert ( make_pair ( "ClusterSize_" + to_string ( _outputSensorID + i ), csizeHist ) );
	csizeHist -> setTitle ( "Cluster Size;Cluster Size;Entries" );

	AIDA::IHistogram1D * cchargeHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "ClusterCharge_" + to_string ( _outputSensorID + i ) ).c_str ( ), 10, -0.5, 9.5 );
	_aidaHistoMap.insert ( make_pair ( "ClusterCharge_" + to_string ( _outputSensorID + i ), cchargeHist ) );
	cchargeHist -> setTitle ( "Cluster Charge;Cluster Charge;Entries" );

	AIDA::IHistogram1D * hitmapHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "Hitmap_" + to_string ( _outputSensorID + i ) ).c_str ( ), _chancount + 1, -0.5, _chancount - 0.5 );
	_aidaHistoMap.insert ( make_pair ( "Hitmap_" + to_string ( _outputSensorID + i ), hitmapHist ) );
	hitmapHist -> setTitle ( "Hitmap;Channel;Entries" );
    }
}
