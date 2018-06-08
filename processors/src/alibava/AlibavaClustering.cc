/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"

// alibava includes ".h"
#include "AlibavaClustering.h"
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"
#include "AlibavaPedNoiCalIOManager.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// aida includes <.h>
#if defined ( USE_AIDA ) || defined ( MARLIN_USE_AIDA )
    #include <marlin/AIDAProcessor.h>
    #include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCEventImpl.h>

// ROOT includes ".h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>

// system includes <>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace alibava;
using namespace eutelescope;

AlibavaClustering::AlibavaClustering ( ) : AlibavaBaseProcessor ( "AlibavaClustering" ),
_clusterCollectionName ( ALIBAVA::NOTSET ),
_clustercharge ( ),
_clustercount ( 0 )
{

    // modify processor description
    _description = "AlibavaClustering does lots of things. It clusters non-zero supressed strip sensor data and can perform FIR filtering beforehand.";

    // first register the input collection
    registerInputCollection ( LCIO::TRACKERDATA, "InputCollectionName", "Input data collection name", _inputCollectionName, string ( "data" ) );

    // noise input
    registerProcessorParameter ( "NoiseInputFile", "The filename where the final noise is stored", _pedestalFile, string ( "finalnoise.slcio" ) );

    // seed cut
    registerProcessorParameter ( "SeedCut", "The seed SNR cut", _seedcut , float ( 3.0 ) );

    // cluster cut
    registerProcessorParameter ( "ClusterCut", "The cluster SNR cut", _clustercut , float ( 1.5 ) );

    // the polarity 
    registerProcessorParameter ( "Polarity", "The sensor polarity: -1 for negative cluster signals (p-type sensor), 1 for positive cluster signals (n-type sensor)", _polarity,  -1 );

    // now the optional parameters
    registerProcessorParameter ( "NoiseCollectionName", "The noise collection name", _noiseCollectionName, string ( "finalnoise" ) );

    // the name of the cluster collection
    registerProcessorParameter ( "ClusterCollectionName", "The cluster collection name", _clusterCollectionName, string ( "clustercollection" ) );

    // the unsensitve axis of the strip sensor
    registerProcessorParameter ( "UnsensitiveAxis", "The unsensitive axis of our strip sensor", _nonsensitiveaxis, string ( "x" ) );

    // the name of the sparse cluster collection
    registerProcessorParameter ( "SparseClusterCollectionName", "The sparse cluster collection name, this needs to be original_zsdata for hitmaker", _sparseclusterCollectionName, string ( "original_zsdata" ) );

    registerOptionalParameter ( "MinClustersize", "The minimum accepted clustersize in the sensitive Alibava direction. This must be larger or equal to 1!", _clusterminsize,  1 );

    registerOptionalParameter ( "MaxClustersize", "The maximum accepted clustersize in the sensitive Alibava direction. This should be larger or equal to MinClustersize!", _clustermaxsize, 99 );

    registerOptionalParameter ( "UseFIRFilter", "A FIR (finite impulse response) filter can be applied to the input data to minimise crosstalk. This switches the filter on.", _usefir, false );

    registerOptionalParameter ( "WriteFIRCoefficients", "From the eta distribution coefficients for filtering are calculated. This writes them to disk.", _writecoefficients, false );

    registerOptionalParameter ( "WriteZeroCoefficients", "For compatibility, two zero coefficients can be written to file.", _writezero, false );

    registerOptionalParameter ( "ReadFIRCoefficients", "FIR filter coefficients from a previous iteration can be read if this is switched on.", _readcoefficients, false );

    registerOptionalParameter ( "FIRCoefficientFile", "The filename to read/write coefficients to", _filterFileName, string ( "filtercoefficients.txt" ) );

    registerOptionalParameter ( "FIRCollectionName", "If filtering is used, it is saved into this collection", _filteredCollectionName, string ( "filteredcollection" ) );

}

void AlibavaClustering::init ( )
{

    streamlog_out ( MESSAGE4 ) << "Running init" << endl;

    printParameters ( );

    if ( Global::parameters -> isParameterSet ( ALIBAVA::CHANNELSTOBEUSED ) )
    {
	Global::parameters -> getStringVals ( ALIBAVA::CHANNELSTOBEUSED, _channelsToBeUsed );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The Global Parameter " << ALIBAVA::CHANNELSTOBEUSED << " is not set! All channels will be used!" << endl;
    }

    if ( Global::parameters -> isParameterSet ( ALIBAVA::SKIPMASKEDEVENTS ) )
    {
	_skipMaskedEvents = bool ( Global::parameters -> getIntVal ( ALIBAVA::SKIPMASKEDEVENTS ) );
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "The Global Parameter " << ALIBAVA::SKIPMASKEDEVENTS << " is not set! Masked events will be used!" << endl;
    }

    _readcoefficient1 = 0.0;
    _readcoefficient2 = 0.0;

    // I assume this is good for about 3m of flatband ribbon between alibava daughterboard and alibava motherboard in the DESY testbeam
    // this might differ for other setups, number was found by hand :-)
    _initcoefficient1 = 0.0373;
    _initcoefficient2 = 0.0162;

    if ( _usefir == true && _readcoefficients == true )
    {

	ifstream fileRead;
	fileRead.open ( _filterFileName.c_str ( ) );

	if ( fileRead.is_open ( ) )
	{
	    fileRead >> _readcoefficient1 >> _readcoefficient2;
	    streamlog_out ( MESSAGE4 ) << "Filter coefficients successfully loaded and set to " << _readcoefficient1 << " and " << _readcoefficient2 << " !" << endl;
	}
	else
	{
	    streamlog_out ( ERROR1 ) << "Unable to open file " << _filterFileName << " !" << endl;
	}
    }
}

void AlibavaClustering::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    auto arunHeader = std::make_unique < AlibavaRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    setChipSelection( arunHeader -> getChipSelection ( ) );

    setChannelsToBeUsed ( );

    setPedestals ( );

    bookHistos ( );

}

void AlibavaClustering::processEvent ( LCEvent * anEvent )
{
    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( anEvent );

    // the collection we read
    LCCollectionVec * inputCollectionVec;

    // the collection we output
    LCCollectionVec * clusterCollection;

    // the sparse collection 
    LCCollectionVec * sparseClusterCollectionVec = nullptr;
    sparseClusterCollectionVec =  new LCCollectionVec ( LCIO::TRACKERDATA );

    // the filtered collection
    filteredcollectionVec = new LCCollectionVec ( LCIO::TRACKERDATA );

    try
    {
	clusterCollection = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( _clusterCollectionName ) );
	streamlog_out ( DEBUG5 ) << "clusterCollection exists..." <<  endl;
    }
    catch ( lcio::DataNotAvailableException& e )
    {
	clusterCollection = new LCCollectionVec ( LCIO::TRACKERPULSE );
	streamlog_out ( DEBUG5 ) << "clusterCollection doesn't exist..." <<  endl;
    }

    // find the seed clusters on our data
    try
    {
	// give the collection vec its data
	inputCollectionVec = dynamic_cast < LCCollectionVec * > ( alibavaEvent -> getCollection ( getInputCollectionName ( ) ) ) ;

	// loop over detectors
	int noOfDetector = inputCollectionVec -> getNumberOfElements ( );
	for ( int i = 0; i < noOfDetector; ++i ) 
	{
	    TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( inputCollectionVec -> getElementAt ( i ) ) ;

	    // Only call the actual clustering if the event is not masked (if skipping masked events is activated).
	    // If this is skipped, the outputcollection is still written ( then empty), so the following processors stay in sync!
	    if ( ( _skipMaskedEvents && ( alibavaEvent -> isEventMasked ( ) == false ) ) || _skipMaskedEvents == false )
	    {
		streamlog_out ( DEBUG4 ) << "Calling clustering function in event: " << anEvent -> getEventNumber ( ) << " !" << endl;

		// let the clustering begin!
		findSeedClusters ( trkdata, clusterCollection, sparseClusterCollectionVec, alibavaEvent );
	    }
	}

	// Write an output event, even if there is no cluster in this event, this makes hitmaking and merging easier
	// If we skipped the call to findSeedClusters above because an event is masked, we still need to define the cellidencoding for the empty collection we are writing...
	if ( _skipMaskedEvents && ( alibavaEvent -> isEventMasked ( ) ) )
	{
	    CellIDEncoder < TrackerPulseImpl > pulseEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, clusterCollection );
	    CellIDEncoder < TrackerDataImpl > dataEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );
	}

	alibavaEvent -> addCollection ( sparseClusterCollectionVec, _sparseclusterCollectionName );
	alibavaEvent -> addCollection ( clusterCollection, getClusterCollectionName ( ) );

	if ( _usefir == true )
	{
	    alibavaEvent -> addCollection ( filteredcollectionVec, _filteredCollectionName );
	}

    }
    catch ( lcio::DataNotAvailableException& )
    {
	// do nothing
	streamlog_out ( ERROR5 ) << "Collection (" << getInputCollectionName ( ) << ") not found in event " << anEvent -> getEventNumber ( ) << " ! " << endl;
    }
}

void AlibavaClustering::findSeedClusters(TrackerDataImpl * trkdata, LCCollectionVec * clusterCollection, LCCollectionVec * sparseClusterCollectionVec, AlibavaEventImpl * alibavaEvent )
{
    // magic is done here
    streamlog_out ( DEBUG4 ) << "Find Seed loop " << alibavaEvent -> getEventNumber ( ) << endl;

    // the TDC of this event
    float tdc = alibavaEvent -> getEventTime ( );

    // the chip number
    int chipnum = getChipNum ( trkdata );
    streamlog_out ( DEBUG4 ) << "Chip " << chipnum << endl;

    // channel spacer
    int dchip = 0;
    if ( chipnum == 1 )
    {
	dchip = ALIBAVA::NOOFCHANNELS;
    }

    // the data
    FloatVec datavec;
    datavec = trkdata -> getChargeValues ( );

    // output the filtered data too
    FloatVec newdatavec;

    // do FIR filtering beforehand
    if ( _usefir == true )
    {

	// the output dataimpl
	TrackerDataImpl * newdataImpl = new TrackerDataImpl ( );

	// also decoders
	CellIDEncoder < TrackerDataImpl > chipIDEncoder ( ALIBAVA::ALIBAVADATA_ENCODE, filteredcollectionVec );
	chipIDEncoder[ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM] = chipnum;
	chipIDEncoder.setCellID ( newdataImpl );

	// read buffer
	double input_buffer[ALIBAVA::NOOFCHANNELS];
	double chargein = 0.0;
	double chargeout = 0.0;

	// set the input buffer by reading in the channels
	for ( int i = 0; i <= ( ALIBAVA::NOOFCHANNELS - 1 ); i++ )
	{
	    input_buffer[i] = datavec.at ( i );
	    chargein += datavec.at ( i );
	}

	// as the beetle chip reads channels 0, 1, 2, ... we have to subtract crosstalk the other way around
	// read reverse and subtract the filtercoefficients
	// the _initcoefficient1/2 always have to be subtracted in later iterations, as the original data in the datavec collection is not altered
	// if this is not done, later iterations will subtract too little crosstalk, as the coefficients will be smaller
	for ( int i = ( ALIBAVA::NOOFCHANNELS - 1 ); i >= 0; i-- )
	{
	    // magic lines:

	    if ( i == 0 || i == 1 || i == ( ALIBAVA::NOOFCHANNELS - 2 ) || i == ( ALIBAVA::NOOFCHANNELS - 1 ) )
	    {
		input_buffer[i] = input_buffer[i];
	    }
	    if ( i > 1 && i < 126 )
	    {
		double c1 = ( _readcoefficient1 + _initcoefficient1 ) * input_buffer[i-1];
		double c2 = ( _readcoefficient2 + _initcoefficient2 ) * input_buffer[i-2];
		input_buffer[i] = input_buffer[i] - c1 - c2;
		input_buffer[i - 1] = input_buffer[i - 1] + c1;
		input_buffer[i - 2] = input_buffer[i - 2] + c2;
	    }

	}

	// write out again
	for ( int i = 0; i < ( ALIBAVA::NOOFCHANNELS - 1 ); i++ )
	{
	    datavec.at ( i ) = input_buffer[i];
	    newdatavec.push_back ( input_buffer[i] );
	    chargeout += input_buffer[i];
	}

	if ( ( chargein - chargeout ) > 1 )
	{
	    streamlog_out ( ERROR1 ) << "Charge loss in filtering is : " << chargein - chargeout << " ADCs!" << endl;
	}

	// save output
	newdataImpl -> setChargeValues ( newdatavec );
	filteredcollectionVec -> push_back ( newdataImpl );
    }

    // cluster count in this event
    int nClusters = 0;

    // a vector with the cluster info, 0 for no cluster and nClusters for the individual cluster,
    // e.g. {0,0,0,0,0,0}, then  {0,1,0,0,0,0} after finding the first cluster, {0,1,0,2,2,0} after finding the second, etc.
    std::vector < int > clusterNumber;
    clusterNumber.assign ( datavec.size ( ), 0 );

    // go through all channels
    for ( size_t ichan = 0; ichan < datavec.size ( );ichan++ )
    {
	streamlog_out ( DEBUG3 ) << "Channel loop " << ichan << endl;

	float noise = 0.0;
	float seedlevel = 0.0;
	float clusterlevel = 0.0;

	// this is the noise level we get from the previous steps
	if ( isMasked ( chipnum, ichan ) == false )
	{
	    // get the noise levels
	    noise = getNoiseAtChannel ( chipnum, ichan );
	    streamlog_out ( DEBUG4 ) << "Channel " << ichan << " has a noise level of | " << noise << " | ADCs" << endl;
	    streamlog_out ( DEBUG4 ) << "Channel " << ichan << " has a signal level of " << datavec[ichan] << " ADCs" << endl;

	    // the level for a seed
	    seedlevel = _seedcut * noise;
	    streamlog_out ( DEBUG4 ) << "Seed level is " << seedlevel * _polarity << endl;

	    // the level for a cluster
	    clusterlevel = _clustercut * noise;
	    streamlog_out ( DEBUG4 ) << "Cluster level is " << clusterlevel * _polarity << endl;

	}
	else
	{
	    // all levels should be 0
	    streamlog_out ( DEBUG3 ) << "Channel " << ichan << " is masked. Setting noise level to " << noise << endl;
	}

	// if a channel passes the seed cluster clut, use only positive ADCs to check! The noise levels from previous steps will always be positive!
	// we also disalow a seed next to a bad channel
	if ( ( datavec[ichan] * _polarity ) >= seedlevel && ( isMasked ( chipnum, ichan ) == false ) && ( isMasked ( chipnum, ichan + 1 ) == false ) && ( isMasked ( chipnum, ichan - 1 ) == false ) )
	{

	    streamlog_out ( DEBUG4 ) << "Found seed on channel " << ichan << " with seed of " << datavec[ichan] * _polarity << endl;

	    // now check if any neighbours pass the cluster cut
	    // the cluster starts out with size 1, so nothing left or right
	    int posclustersize = 0;
	    int negclustersize = 0;

	    // look right
	    bool posdir = true;
	    // does this channel exist?
	    while ( ( ichan + posclustersize + 1 ) < datavec.size ( ) && posdir == true )
	    {
		if ( ( datavec[ichan + posclustersize + 1] * _polarity ) >= ( _clustercut * getNoiseAtChannel ( chipnum, ichan + posclustersize + 1 ) ) && ( isMasked ( chipnum, ichan + posclustersize + 1 ) == false ) )
		{

		    posclustersize += 1;
		    streamlog_out ( DEBUG3 ) << "Found large cluster in event: " << alibavaEvent -> getEventNumber ( ) << endl;
		    // if this happens, we have found a second seed in this cluster with a higher ADC count
		    // since ichan goes to the right, this feature can only be found here.

		    // don't compare ACDs, look for ADC/seedlevel
		    if ( ( ( datavec[ichan + posclustersize + 1] * _polarity ) / ( getNoiseAtChannel ( chipnum, ichan + posclustersize + 1 ) * _seedcut )  ) > ( ( datavec[ichan] * _polarity ) / seedlevel ) )
		    {

			streamlog_out ( DEBUG4 ) << "Found multiple seeds in a cluster in event " << alibavaEvent -> getEventNumber ( ) << " ... splitting!" << endl;
			streamlog_out ( DEBUG4 ) << ( ( datavec[ichan + posclustersize + 1] * _polarity ) / ( getNoiseAtChannel ( chipnum, ichan + posclustersize + 1 ) * _seedcut )  ) << " > " << ( ( datavec[ichan] * _polarity ) / seedlevel ) << endl;
			// two options: if they touch, we simply move the seed position to the right one and inc negativesize, if they don't touch, we split:
			if ( posclustersize == 1 )
			{

			    streamlog_out ( DEBUG5 ) << "Moving the seed in event: " << alibavaEvent -> getEventNumber ( ) << endl;

			    // move ichan and negclustersize
			    ichan++;
			    negclustersize++;

			    // now posclustersize has to go down, since we shifted the cluster seed position 
			    posclustersize--;
			}

			// if there are clusters below the seed between two seeds, we give them to the highest and split
			if ( posclustersize >= 2 )
			{
			    streamlog_out ( DEBUG5 ) << "Splitting the cluster between seeds in event: " << alibavaEvent -> getEventNumber ( ) << endl;

			    // setting this to zero, let the higher seed find the inbetween strip when looping
			    posclustersize = 0;
			    posdir = false;
			}

			if ( posclustersize >= 3 )
			{
			    // just report 
			    streamlog_out ( WARNING5 ) << "Multiple clusters between seeds in event: " << alibavaEvent -> getEventNumber ( ) << " ! The noise levels could be false!" << endl;
			}

		    }
		}
		else
		{
		    posdir = false;
		}
	    }

	    // look left
	    bool negdir = true;
	    // does this channel exist?
	    while ( negdir == true )
	    {
		// the channel has to fullfill 3 conditions: a) be above cluster level, b) not be masked and c) not be above the seed cut. c) is so that we don't refind clusters after splitting.
		if ( ( datavec[ichan - negclustersize - 1] * _polarity ) >= ( _clustercut * getNoiseAtChannel ( chipnum, ichan - negclustersize - 1 ) ) && ( isMasked ( chipnum, ichan - negclustersize - 1 ) == false ) && ( datavec[ichan - negclustersize - 1] * _polarity ) < ( _seedcut * getNoiseAtChannel ( chipnum, ichan - negclustersize - 1 ) ) )
		{
		    streamlog_out ( DEBUG3 ) << "Found large cluster in event: " << alibavaEvent -> getEventNumber ( ) << endl;
		    negclustersize += 1;
		}
		else
		{
		    negdir = false;
		}
	    }

	    // save the cluster information
	    int lowerlimit = ichan - negclustersize;
	    int upperlimit = ichan + posclustersize;
	    int clusize = upperlimit - lowerlimit + 1;

	    // allow a cut on the clustersize
	    if ( ( clusize >= _clusterminsize ) && ( clusize <= _clustermaxsize ) )
	    {

		// anything past here is an accepted cluster
		streamlog_out ( DEBUG5 ) << "Found a cluster in event: " << alibavaEvent -> getEventNumber ( ) << endl;

		// increment cluster count
		nClusters++;

		// debug output
		if ( clusize > 1 )
		{
		    streamlog_out ( DEBUG4 ) << "Clustersize: " << clusize << endl;
		}

		// fill the clusterNumber vector with the position of the individual clusters
		for ( int icluster = lowerlimit; icluster <= upperlimit; icluster++ )
		{
		    clusterNumber.at ( icluster ) = nClusters;
		    streamlog_out ( DEBUG3 ) << "Wrote cluster to clusterNumber at position(s): " <<  icluster << endl;
		}

		// record the charge to the left and right of the seed to spot asymetries:
		if ( isMasked ( chipnum, ichan - 2 ) == false && isMasked ( chipnum, ichan - 1 ) == false && isMasked ( chipnum, ichan + 1 ) == false && isMasked ( chipnum, ichan + 2 ) == false )
		{
		    _clustercharge[0] += fabs ( datavec[ichan - 2] );
		    _clustercharge[1] += fabs ( datavec[ichan - 1] );
		    _clustercharge[2] += fabs ( datavec[ichan] );
		    _clustercharge[3] += fabs ( datavec[ichan + 1] );
		    _clustercharge[4] += fabs ( datavec[ichan + 2] );
		    fillChargeDistHisto ( datavec[ichan - 3] * _polarity, datavec[ichan - 2] * _polarity, datavec[ichan - 1] * _polarity, datavec[ichan] * _polarity, datavec[ichan + 1] * _polarity, datavec[ichan + 2] * _polarity, datavec[ichan + 3] * _polarity );
		}

		// let's calculate the eta distribution, only works for clustersize >1
		// actually this is not the official definition of eta, but this helps nonetheless
		if ( clusize > 1 )
		{
		    float etaleft = 0;
		    float etaright = 0;
		    float cogpos = 0;
		    float totalsignal = 0;
		    float weight = 0;
		    for ( int i = lowerlimit; i <= upperlimit; i++ )
		    {
			totalsignal += fabs ( datavec[i] );
			weight += fabs ( datavec[i] * i );
		    }
		    cogpos = weight / totalsignal;
		    for ( int j = lowerlimit; j <= upperlimit; j++ )
		    {
			if ( float ( j ) < cogpos )
			{
			    etaleft += fabs ( datavec[j] );
			}
			else if ( float ( j ) > cogpos )
			{
			    etaright += fabs ( datavec[j] );
			}
		    }
		    float etaratio = etaleft / ( etaleft + etaright );
		    fillEtaHisto ( etaratio );

		    // CoG control plot:
		    if ( clusize >= 2 )
		    {
			fillCogHisto ( cogpos - ichan );
		    }
		}

		// do the real calculation of eta:
		// require good channels left and right so we don't bias
		if ( isMasked ( chipnum, ichan - 1 ) == false && isMasked ( chipnum, ichan + 1 ) == false )
		{
		    float etaleft_a = 0;
		    float etaright_a = 0;
		    float etaleft_b = 0;
		    float etaright_b = 0;
		    float etaratio = -1;

		    // divison by noise can be added by removeing comments
		    etaleft_a = datavec[ichan - 1] * _polarity;// / getNoiseAtChannel(chipnum,ichan-1);
		    etaright_a = datavec[ichan] * _polarity;// / getNoiseAtChannel(chipnum,ichan);

		    etaright_b = datavec[ichan + 1] * _polarity;// / getNoiseAtChannel(chipnum,ichan+1);
		    etaleft_b = datavec[ichan] * _polarity;// / getNoiseAtChannel(chipnum,ichan);

		    if ( etaleft_a > etaright_b )
		    {
			etaratio = etaleft_a / ( etaleft_a + etaright_a );
			fillEtaHisto2 ( etaratio );
			fillEtaHisto2TDC ( etaratio, tdc );
			fillEtaHistoPos ( etaratio, ichan + dchip );
			streamlog_out ( DEBUG2 ) << "Eta2: left side ratio written: " << etaratio << endl; 
		    }
		    else if ( etaright_b > etaleft_a )
		    {
			etaratio = etaleft_b / ( etaleft_b + etaright_b );
			fillEtaHisto2 ( etaratio );
			fillEtaHisto2TDC ( etaratio, tdc );
			fillEtaHistoPos ( etaratio, ichan + dchip );
			streamlog_out ( DEBUG2 ) << "Eta2: right side ratio written:" << etaratio << endl;
		    }
		    else
		    {
			streamlog_out ( DEBUG2 ) << "Eta2: no eta here (" << alibavaEvent -> getEventNumber ( ) << ") because fail! ela:" << etaleft_a << " era: " << etaright_a << " elb: "<< etaleft_b << " erb: " << etaright_b << endl;
		    }
		}

		// fill the cluster signal and snr histos
		float clustersignal = 0.0;
		float clusternoise = 0.0;
		for ( int i = ( lowerlimit ); i <= ( upperlimit ); i++ )
		{
		    clustersignal += _polarity * datavec[i];
		    clusternoise += getNoiseAtChannel ( chipnum, i );
		}
		fillSignalHisto ( clustersignal );
		fillSNRHisto ( clustersignal / clusternoise );

		// fill the seed charge histogram
		fillSeedChargeHisto ( _polarity * datavec[ichan] );

		// the overall cluster count in this run
		_clustercount++;

		// fill the clustersize histo
		fillClusterHisto ( clusize );

		// fill the hitmap histogram
		fillHitmapHisto ( ichan + dchip, negclustersize, posclustersize );

		// fill the seed histogram
		fillSeedHisto ( ichan + dchip );

		// final: move the ichan var so that we dont include a seed channel in multiple clusters
		if ( ( ichan + posclustersize + 1 ) < datavec.size ( ) )
		{
		    // this channel is the first one next to a previous cluster and below
		    // the clustercut so it will fail the  seed check in the next check
		    ichan = ichan + posclustersize + 1;
		}

	    }
	    else
	    {
		streamlog_out ( DEBUG5 ) << "Failed clustersize cut! Limits are: " << _clusterminsize << " to " << _clustermaxsize << " , this cluster has size: " << clusize << " !" << endl;
	    } // done clustersize cut

	} // done going over all seeds
    } // done going over all channels in this event

    // more debug output
    if ( nClusters >= 1 )
    {
	streamlog_out ( DEBUG4 ) << "Clusters in this event: " << nClusters << endl;
    }

    fillclusterspereventhisto ( nClusters, alibavaEvent -> getEventNumber ( ) );

    // now output the clusters we have found
    // encoders:
    CellIDEncoder < TrackerPulseImpl > pulseEncoder ( eutelescope::EUTELESCOPE::PULSEDEFAULTENCODING, clusterCollection );
    CellIDEncoder < TrackerDataImpl > dataEncoder ( eutelescope::EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );
    //CellIDEncoder<TrackerDataImpl> dataEncoder (EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );

    // What have we found?
    for ( unsigned int icluster = 0; icluster < clusterNumber.size ( ); icluster++ )
    {
	streamlog_out ( DEBUG3 ) << "Pos: " << icluster << " , Value: " << clusterNumber.at ( icluster ) << endl;
    }

    // Decode datavec into pixelvec, each strip goes into a Pixel
    // Move to positive ADCs only, since identification has already happend. This is for EUTelSimpleSparsePixel.

    //std::vector<eutelescope::EUTelGenericSparsePixel* > PixelVec;
    std::vector < EUTelGenericSparsePixel > PixelVec;
    eutelescope::EUTelGenericSparsePixel Pixel;

    for ( unsigned int iPixel = 0; iPixel < datavec.size ( ); iPixel++ )
    {

	// The missing coordinate, as passed through by the struct
	float missingvalue = 0.0;

	// chip 1 gets channels inc'ed by 128
	// select the sensor orientation, we give each "pixel" the missing coordinate on the unsensitive axis
	if ( _nonsensitiveaxis == "x" )
	{
	    if ( chipnum == 0 )
	    {
		Pixel.setXCoord ( missingvalue );
		Pixel.setYCoord ( iPixel );
	    }
	    if ( chipnum == 1 )
	    {
		Pixel.setXCoord ( missingvalue );
		Pixel.setYCoord ( iPixel + ALIBAVA::NOOFCHANNELS );
	    }
	}
	if ( _nonsensitiveaxis == "y" )
	{
	    if ( chipnum == 0 )
	    {
		Pixel.setXCoord ( iPixel );
		Pixel.setYCoord ( missingvalue );
	    }
	    if ( chipnum == 1 )
	    {
		Pixel.setXCoord ( iPixel + ALIBAVA::NOOFCHANNELS );
		Pixel.setYCoord ( missingvalue );
	    }
	}

	// Charge times 100 to get precision, since telescope cluster charge is of type short!
	// This may be changed in future eutelescope releases
	Pixel.setSignal ( datavec[iPixel] * _polarity * 100.0 );
	PixelVec.push_back ( Pixel );

	streamlog_out ( DEBUG3 ) << "Decoded pixel with coordinates: ( " << Pixel.getXCoord ( ) << " , " << Pixel.getYCoord ( ) << " ) and charge: " << datavec[iPixel] << " !" << endl;
	streamlog_out ( DEBUG3 ) << "Input was: " << PixelVec.at ( iPixel ) << endl;
    }

    // now let's move this out into trackerpulses and trackerdata
    // iterate over the clusters
    std::vector < int > ::iterator it;
    for ( it = clusterNumber.begin ( ); it != clusterNumber.end ( ); ++it )
    {
	std::unique_ptr < TrackerDataImpl > zsCluster = std::make_unique < TrackerDataImpl > ( );
        // prepare a reimplementation of sparsified cluster
	auto sparseCluster = std::unique_ptr < EUTelTrackerDataInterfacer > ( new EUTelTrackerDataInterfacerImpl < EUTelGenericSparsePixel > ( zsCluster.get ( ) ) );

	for ( unsigned int i = 0; i < datavec.size ( ); i++ )
	{
	    if ( ( clusterNumber.at ( i ) == *it ) && ( clusterNumber.at ( i ) > 0 ) )
	    {
		// put only these pixels in that ClusterCollection that belong to that cluster
		auto temppix = PixelVec.at ( i );
		sparseCluster -> push_back ( temppix );
		streamlog_out ( DEBUG6 ) << "Adding strip " << i << " to cluster " << clusterNumber.at ( i ) << endl;
		++it;
	    }
	}

	// now we have pixelClusters with the corresponding pixels (which have the xy and q info) in them
	// this if stops making the "zero" cluster with all strips not in a cluster...
	unsigned int clumin = 0;
	unsigned int clumax = 0;
	clumin = static_cast<unsigned int>( _clusterminsize );
	clumax = static_cast<unsigned int>( _clustermaxsize );
	if ( ( sparseCluster -> size ( ) >= clumin ) && ( sparseCluster -> size ( ) <= clumax ) )
	{

	    // this assumes six telescope planes, so the first chip will be plane 6, etc.
	    dataEncoder["sensorID"] = 6;
	    dataEncoder["sparsePixelType"] = static_cast < int > ( kEUTelSparseClusterImpl );
	    dataEncoder["quality"] =  0;
	    dataEncoder.setCellID ( zsCluster.get ( ) );
	    sparseClusterCollectionVec -> push_back ( zsCluster.get ( ) );

	    std::unique_ptr < TrackerPulseImpl > zsPulse = std::make_unique < TrackerPulseImpl > ( );
	    pulseEncoder["sensorID"] = 6;
	    pulseEncoder["type"] = static_cast < int > ( kEUTelSparseClusterImpl );
	    pulseEncoder.setCellID ( zsPulse.get ( ) );

	    zsPulse -> setTrackerData ( zsCluster.release ( ) );
	    clusterCollection -> push_back ( zsPulse.release ( ) );

	}
    }
}

void AlibavaClustering::check ( LCEvent * /* evt */ )
{
    // nothing to check here - could be used to fill check plots in reconstruction processor
}

void AlibavaClustering::fillclusterspereventhisto ( int clusters, int event )
{
    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap["ClustersVsEvents"] ) )
    {
	histo -> Fill ( event, clusters );
    }
}

void AlibavaClustering::fillHitmapHisto ( int ichan, int negclustersize, int posclustersize )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["ClusterHitmap"] ) )
    {
	for ( int i = ( ichan - negclustersize ); i <= ( ichan + posclustersize ); i++ )
	{
	    histo -> Fill ( i );
	}
    }
}

void AlibavaClustering::fillClusterHisto ( int clusize )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["ClusterSize"] ) )
    {
	histo -> Fill ( clusize );
    }
}

void AlibavaClustering::fillEtaHisto ( float etaratio )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["ClusterCOGEta"] ) )
    {
	histo -> Fill ( etaratio );
    }
}

void AlibavaClustering::fillEtaHisto2(float etaratio)
{
	if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["EtaDistribution"] ) )
	{
		histo->Fill(etaratio);
	}
}

void AlibavaClustering::fillEtaHisto2TDC ( float etaratio, float tdc )
{
    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap["EtaDistributionTDC"] ) )
    {
	histo -> Fill ( etaratio, tdc );
    }
}

void AlibavaClustering::fillEtaHistoPos ( float etaratio, int ichan )
{
    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap["EtaDistributionPos"] ) )
    {
	histo -> Fill ( etaratio, ichan );
    }
}

void AlibavaClustering::fillSeedHisto ( int ichan )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SeedHitmap"] ) )
    {
	histo -> Fill ( ichan );
    }
}

void AlibavaClustering::fillChargeDistHisto (float a, float b, float c, float d, float e, float f, float g )
{
    if ( TH2D * histo = dynamic_cast < TH2D* > ( _rootObjectMap["NeighbourChargeDistribution"] ) )
    {
	histo -> Fill ( -3.0, a / d );
	histo -> Fill ( -2.0, b / d );
	histo -> Fill ( -1.0, c / d );
	histo -> Fill ( 0.0, d / d );
	histo -> Fill ( 1.0, e / d );
	histo -> Fill ( 2.0, f / d );
	histo -> Fill ( 3.0, g / d );
    }

    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SignalLeft2"] ) )
    {
	histo -> Fill ( b / d );
    }

    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SignalLeft1"] ) )
    {
	histo -> Fill ( c / d );
    }

    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SignalRight1"] ) )
    {
	histo -> Fill ( e / d );
    }

    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SignalRight2"] ) )
    {
	histo -> Fill ( f / d );
    }

    // also fill some alternative histos: see thesis of Erik Butz
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["ZetaHisto"] ) )
    {
	histo -> Fill ( ( c + e ) / ( c + d + e ) );
    }
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SigmaHisto"] ) )
    {
	histo -> Fill ( ( c + e ) / ( 2 * d ) );
    }
}

void AlibavaClustering::fillSignalHisto ( float signal )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SignalfromClusters"] ) )
    {
	histo -> Fill ( signal );
    }
}

void AlibavaClustering::fillSNRHisto ( float signal )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SNRfromClusters"] ) )
    {
	histo -> Fill ( signal );
    }
}

void AlibavaClustering::fillSeedChargeHisto ( float signal )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["SignalfromSeeds"] ) )
    {
	histo -> Fill ( signal );
    }
}

void AlibavaClustering::fillCogHisto ( float cog )
{
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["RelativeCoGPosition"] ) )
    {
	histo -> Fill ( cog );
    }
}

void AlibavaClustering::end ( )
{
    // now that we have all the data, do a fit on the signal histos
    dolandaugausfit ( "SignalfromSeeds" );
    dolandaugausfit ( "SignalfromClusters" );

    // also do a fit on the neighbour charge histos
    TH1D * histol2 = dynamic_cast < TH1D* > ( _rootObjectMap["SignalLeft2"] );
    TF1 * fitl2 = dynamic_cast < TF1* > ( _rootObjectMap["Signal Left 2 Fit"] );
    histol2 -> Fit ( fitl2, "Q" );

    TH1D * histol1 = dynamic_cast < TH1D* > ( _rootObjectMap["SignalLeft1"] );
    TF1 * fitl1 = dynamic_cast < TF1* > ( _rootObjectMap["Signal Left 1 Fit"] );
    histol1 -> Fit ( fitl1, "Q" );

    TH1D * histor1 = dynamic_cast < TH1D* > ( _rootObjectMap["SignalRight1"] );
    TF1 * fitr1 = dynamic_cast < TF1* > ( _rootObjectMap["Signal Right 1 Fit"] );
    histor1 -> Fit ( fitr1, "Q" );

    TH1D * histor2 = dynamic_cast < TH1D* > ( _rootObjectMap["SignalRight2"] );
    TF1 * fitr2 = dynamic_cast < TF1* > ( _rootObjectMap["Signal Right 2 Fit"] );
    histor2 -> Fit ( fitr2, "Q" );

    // we can also do a fit on the neighbour eta to see if we still have asymetries...
    TH1D * etaHisto2 = dynamic_cast < TH1D* > ( _rootObjectMap["EtaDistribution"] );
    TF1 * etalfit = dynamic_cast < TF1* > ( _rootObjectMap["Eta Distribution left fit"] );
    etalfit -> SetRange ( -0.5, 0.5 );
    etaHisto2 -> Fit ( etalfit, "QR+" );
    TF1 * etarfit = dynamic_cast < TF1* > ( _rootObjectMap["Eta Distribution right fit"] );
    etarfit -> SetRange ( 0.5, 1.5 );
    etaHisto2 -> Fit ( etarfit, "QR+" );

    // make an integral of the track based eta distribution
    double counts = 0.0;
    double integral = 0.0;
    TH1D * etahisto = dynamic_cast < TH1D* > ( _rootObjectMap["EtaDistribution"] );
    TH1D * etaintegral = dynamic_cast < TH1D* > ( _rootObjectMap["EtaDistributionIntegral"] );
    for ( int i = 1; i < etahisto -> GetNbinsX ( ); i++ )
    {
	counts = etahisto -> GetBinContent ( i );
	integral += counts;
	etaintegral -> SetBinContent ( i, integral );
    }

    // some final output...
    streamlog_out ( MESSAGE4 ) << "===============================" << endl;
    double filtercoefficient1 = ( fitr1 -> GetParameter ( 1 ) - fitl1 -> GetParameter ( 1 ) );
    double filtercoefficient2 = ( fitr2 -> GetParameter ( 1 ) - fitl2 -> GetParameter ( 1 ) );

    streamlog_out ( MESSAGE4 ) << " Left 2:      " << fitl2 -> GetParameter ( 1 ) << endl;
    streamlog_out ( MESSAGE4 ) << " Left 1:      " << fitl1 -> GetParameter ( 1 ) << endl;
    streamlog_out ( MESSAGE4 ) << " Left shift:  " << etalfit -> GetParameter ( 1 ) << endl;
    streamlog_out ( MESSAGE4 ) << " Right shift: " << 1 - etarfit -> GetParameter ( 1 ) << endl;
    streamlog_out ( MESSAGE4 ) << " Right 1:     " << fitr1 -> GetParameter ( 1 ) << endl;
    streamlog_out ( MESSAGE4 ) << " Right 2:     " << fitr2 -> GetParameter ( 1 ) << endl;
    streamlog_out ( MESSAGE4 ) << "===============================" << endl;
    streamlog_out ( MESSAGE4 ) << "Clustercount:" << _clustercount << endl;
    streamlog_out ( MESSAGE4 ) << "===============================" << endl;
    streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;

    // from the eta distribution we obtain the coefficients for a possible re-run of this processor
    // they are written to a file here
    if ( _writecoefficients == true )
    {
	streamlog_out ( MESSAGE4 ) << "Writing filter coefficients to disk..." << endl;
	if ( _writezero == true )
	{
	    filtercoefficient1 = 0.0;
	    filtercoefficient2 = 0.0;
	}
	ofstream filterFile;
	filterFile.open ( _filterFileName.c_str ( ) );
	filterFile << filtercoefficient1 << endl;
	filterFile << filtercoefficient2 << endl;
	filterFile.close ( );
	streamlog_out ( MESSAGE4 ) << "File " << _filterFileName << " written." << endl;
    }
}

void AlibavaClustering::fillHistos ( )
{

}

void AlibavaClustering::bookHistos ( )
{
    AIDAProcessor::tree ( this ) -> cd ( this -> name ( ) );

    TH2D * cluevent = new TH2D ( "ClustersVsEvents", "", 5000, 0, 500000, 10, 0, 9 );
    _rootObjectMap.insert ( make_pair ( "ClustersVsEvents", cluevent ) );
    cluevent -> SetTitle ( "Clusters per Event vs. Event Nr;Event Nr.;Clusters in this Event");

    // a histogram showing the clustersize
    string tempHistoName = "ClusterSize";
    stringstream tempHistoTitle;
    tempHistoTitle << tempHistoName << ";Clustersize;NumberofEntries";

    TH1D * clusterHisto = new TH1D ( tempHistoName.c_str ( ), "", 4, 1, 5 );
    _rootObjectMap.insert ( make_pair ( tempHistoName, clusterHisto ) );
    string tmp_string = tempHistoTitle.str ( );
    clusterHisto -> SetTitle ( tmp_string.c_str ( ) );

    // a histogram showing the eta distribution
    string tempHistoName2 = "ClusterCOGEta";
    stringstream tempHistoTitle2;
    tempHistoTitle2 << tempHistoName2 << ";Eta;NumberofEntries";

    TH1D * etaHisto = new TH1D ( tempHistoName2.c_str ( ), "", 100, 0, 1 );
    _rootObjectMap.insert ( make_pair ( tempHistoName2, etaHisto ) );
    string tmp_string2 = tempHistoTitle2.str ( );
    etaHisto -> SetTitle ( tmp_string2.c_str ( ) );

    // a hitmap histo
    string tempHistoName3 = "ClusterHitmap";
    stringstream tempHistoTitle3;
    tempHistoTitle3 << tempHistoName3 << ";Channel;NumberofEntries";

    TH1D * hitmapHisto = new TH1D ( tempHistoName3.c_str ( ), "", 256, 0, 255 );
    _rootObjectMap.insert ( make_pair ( tempHistoName3, hitmapHisto ) );
    string tmp_string3 = tempHistoTitle3.str ( );
    hitmapHisto -> SetTitle ( tmp_string3.c_str ( ) );

    // a seed histo
    string tempHistoName4 = "SeedHitmap";
    stringstream tempHistoTitle4;
    tempHistoTitle4 << tempHistoName4 << ";Channel;NumberofEntries";

    TH1D * seedHisto = new TH1D ( tempHistoName4.c_str ( ), "", 256, 0, 255 );
    _rootObjectMap.insert ( make_pair ( tempHistoName4, seedHisto ) );
    string tmp_string4 = tempHistoTitle4.str ( );
    seedHisto -> SetTitle ( tmp_string4.c_str ( ) );

    // a histogram showing the eta distribution - alternative calculation
    TH1D * etaintegral = new TH1D ( "EtaDistributionIntegral", "", 120, -1, 2 );
    _rootObjectMap.insert ( make_pair ( "EtaDistributionIntegral", etaintegral ) );
    etaintegral -> SetTitle ( "Eta Distribution Integral;Eta;NumberofEntries" );

    // an integral over this plot
    TH1D * etaHisto2 = new TH1D ( "EtaDistribution", "", 120, -1, 2 );
    _rootObjectMap.insert ( make_pair ( "EtaDistribution", etaHisto2 ) );
    etaHisto2 -> SetTitle ( "Eta Distribution;Eta;NumberofEntries" );

    // two fits for this eta plot
    TF1 * etalfit = new TF1 ( "Eta Distribution left fit", "gaus" );
    _rootObjectMap.insert ( make_pair ( "Eta Distribution left fit", etalfit ) );

    TF1 * etarfit = new TF1 ( "Eta Distribution right fit", "gaus" );
    _rootObjectMap.insert ( make_pair ( "Eta Distribution right fit", etarfit ) );

    // a histogram showing the charge distribution
    string tempHistoName6 = "NeighbourChargeDistribution";
    stringstream tempHistoTitle6;
    tempHistoTitle6 << tempHistoName6 << ";Distance to seed;Strip charge / Seed charge";

    TH2D * chargehisto = new TH2D ( tempHistoName6.c_str ( ), "", 7, -3.5, 3.5, 1000, -1, 1 );
    _rootObjectMap.insert ( make_pair ( tempHistoName6, chargehisto ) );
    string tmp_string6 = tempHistoTitle6.str ( );
    chargehisto -> SetTitle ( tmp_string6.c_str ( ) );

    // a histogram showing the cluster signals
    string tempHistoName7 = "SignalfromClusters";
    stringstream tempHistoTitle7;
    tempHistoTitle7 << tempHistoName7 << ";Cluster signal (ADCs) * (-1);Number of Entries";

    TH1D * clustersignalhisto = new TH1D ( tempHistoName7.c_str ( ), "", 1000, 0, 100 );
    _rootObjectMap.insert ( make_pair ( tempHistoName7, clustersignalhisto ) );
    string tmp_string7 = tempHistoTitle7.str ( );
    clustersignalhisto -> SetTitle ( tmp_string7.c_str ( ) );

    // a histogram showing the cluster signal to noise ratio
    string tempHistoName8 = "SNRfromClusters";
    stringstream tempHistoTitle8;
    tempHistoTitle8 << tempHistoName8 << ";Cluster SNR;Number of Entries";

    TH1D * clustersnrhisto = new TH1D ( tempHistoName8.c_str ( ), "", 500, 0, 50 );
    _rootObjectMap.insert ( make_pair ( tempHistoName8, clustersnrhisto ) );
    string tmp_string8 = tempHistoTitle8.str ( );
    clustersnrhisto -> SetTitle ( tmp_string8.c_str ( ) );

    // a histogram showing the seed charge
    string tempHistoName9 = "SignalfromSeeds";
    stringstream tempHistoTitle9;
    tempHistoTitle9 << tempHistoName9 << ";Seed charge (ADCs) * (-1);Number of Entries";

    TH1D * seedchargehisto = new TH1D ( tempHistoName9.c_str ( ), "", 1000, 0, 100 );
    _rootObjectMap.insert ( make_pair ( tempHistoName9, seedchargehisto ) );
    string tmp_string9 = tempHistoTitle9.str ( );
    seedchargehisto -> SetTitle ( tmp_string9.c_str ( ) );

    // a cog control plot
    TH1D * cogplot = new TH1D ( "RelativeCoGPosition", "", 1000, -2, 2 );
    _rootObjectMap.insert ( make_pair ( "RelativeCoGPosition", cogplot ) );
    cogplot -> SetTitle ( "Relative CoG Position;Relative CoG;Number of Entries" );

    // eta vs tdc
    TH2D * etatdc = new TH2D ( "EtaDistributionTDC", "", 100, 0, 1, 20, 0, 100 );
    _rootObjectMap.insert ( make_pair ( "EtaDistributionTDC", etatdc ) );
    etatdc -> SetTitle ( "Eta distribution vs. Event TDC;Eta;TDC time [ns]" );

    // eta vs pos
    TH2D * etapos = new TH2D ( "EtaDistributionPos", "", 100, 0, 1, 256, 0, 255 );
    _rootObjectMap.insert ( make_pair ( "EtaDistributionPos", etapos ) );
    etapos -> SetTitle ( "Eta distribution vs. Seed Position;Eta;Seed Strip" );

    // charge distribution plots for 2 neighbours left and right of a seed
    TH1D * neighbour2left = new TH1D ( "SignalLeft2", "", 100, -2, 2 );
    _rootObjectMap.insert ( make_pair ( "SignalLeft2", neighbour2left ) );
    neighbour2left -> SetTitle ( "Signal Left 2;Relative Charge to Seed;Number of Entries" );

    TH1D * neighbour1left = new TH1D ( "SignalLeft1", "", 100, -2, 2 );
    _rootObjectMap.insert ( make_pair ( "SignalLeft1", neighbour1left ) );
    neighbour1left -> SetTitle ( "Signal Left 1;Relative Charge to Seed;Number of Entries" );

    TH1D * neighbour1right = new TH1D ("SignalRight1", "", 100, -2, 2 );
    _rootObjectMap.insert ( make_pair ( "SignalRight1", neighbour1right ) );
    neighbour1right -> SetTitle ( "Signal Right 1;Relative Charge to Seed;Number of Entries" );

    TH1D * neighbour2right = new TH1D ( "SignalRight2", "", 100, -2, 2 );
    _rootObjectMap.insert ( make_pair ( "SignalRight2", neighbour2right ) );
    neighbour2right-> SetTitle ( "Signal Right 2;Relative Charge to Seed;Number of Entries" );

    TF1 * nl2 = new TF1 ( "Signal Left 2 Fit", "gaus" );
    _rootObjectMap.insert ( make_pair ( "Signal Left 2 Fit", nl2 ) );

    TF1 * nl1 = new TF1 ( "Signal Left 1 Fit", "gaus" );
    _rootObjectMap.insert ( make_pair ( "Signal Left 1 Fit", nl1 ) );

    TF1 * nr1 = new TF1 ( "Signal Right 1 Fit", "gaus" );
    _rootObjectMap.insert ( make_pair ( "Signal Right 1 Fit", nr1 ) );

    TF1 * nr2 = new TF1 ( "Signal Right 2 Fit", "gaus" );
    _rootObjectMap.insert ( make_pair ( "Signal Right 2 Fit", nr2 ) );

    TH1D * zetahisto = new TH1D ( "ZetaHisto", "", 600, -3, 3 );
    _rootObjectMap.insert ( make_pair ( "ZetaHisto", zetahisto ) );
    zetahisto -> SetTitle ( "Zeta Distribution;Zeta;Number of Entries" );

    TH1D * sigmahisto = new TH1D ( "SigmaHisto", "", 600, -3, 3 );
    _rootObjectMap.insert ( make_pair ( "SigmaHisto", sigmahisto ) );
    sigmahisto -> SetTitle ( "Sigma Distribution;Sigma;Number of Entries" );

    streamlog_out ( MESSAGE1 )  << "End of Booking histograms. " << endl;
}

// getter and setter for _clusterCollectionName
void AlibavaClustering::setClusterCollectionName ( std::string clusterCollectionName )
{
    _clusterCollectionName = clusterCollectionName;
}

std::string AlibavaClustering::getClusterCollectionName ( )
{
    return _clusterCollectionName;
}

// Landau-Gauss fitting is done here:
void AlibavaClustering::dolandaugausfit ( string tempHistoName )
{
    streamlog_out ( DEBUG5 ) << "Fitting landau gaus on histo " << tempHistoName << endl;
    if ( TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap[tempHistoName] ) )
    {
	// Setting fit range and start values
	Double_t fr[2];
	Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
	fr[0] = 0.7 * histo -> GetMean ( );
	fr[1] = 5.0 * histo -> GetMean ( );
	pllo[0] = 0.05;
	pllo[1] = 0.50;
	pllo[2] = 0.1;
	pllo[3] = 0.04;
	plhi[0] = 50.0;
	plhi[1] = 500.0;
	plhi[2] = 1000000.0;
	plhi[3] = 50.0;
	sv[0] = 1.8;
	sv[1] = 20.0;
	sv[2] = 50000.0;
	sv[3] = 3.0;

	Double_t chisqr;
	Int_t ndf;
	TF1 *fitsnr = langaufit ( histo, fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf );
	histo -> Fit ( fitsnr, "QR" );

	Double_t SNRPeak, SNRFWHM;
	langaupro ( fp, SNRPeak, SNRFWHM );

	streamlog_out ( MESSAGE4 ) << "===============================" << endl;
	streamlog_out ( DEBUG5 ) << "Done fitting landau gaus on histo" << tempHistoName << endl;
	streamlog_out ( MESSAGE4 ) << "Landau-Gaus Peak for : " << tempHistoName << " is: " << SNRPeak << endl;
	streamlog_out ( MESSAGE4 ) << "Landau-Gaus Sigma for: " << tempHistoName << " is: " << SNRFWHM << endl;
	streamlog_out ( MESSAGE4 ) << "===============================" << endl;
    }

}

Double_t langaufun ( Double_t *x, Double_t *par )
{
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation), 
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // Numeric constants
    Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
    Double_t mpshift = -0.22278298; // Landau maximum location

    // Control constants
    Double_t np = 100.0; // number of convolution steps
    Double_t sc = 5.0; // convolution extends to +-sc Gaussian sigmas

    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0]; 

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = ( xupp - xlow ) / np;

    // Convolution integral of Landau and Gaussian by sum
    for ( i = 1.0; i <= np / 2; i++ )
    {
	xx = xlow + ( i - 0.5 ) * step;
	fland = TMath::Landau ( xx, mpc,par[0] ) / par[0];
	sum += fland * TMath::Gaus ( x[0], xx, par[3] );
	xx = xupp - ( i - 0.5 ) * step;
	fland = TMath::Landau ( xx, mpc,par[0] ) / par[0];
	sum += fland * TMath::Gaus ( x[0], xx, par[3] );
    }

    return ( par[2] * step * sum * invsq2pi / par[3] );
}

TF1 *langaufit ( TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF )
{
    // Variables for langaufit call:
    //   his             histogram to fit
    //   fitrange[2]     lo and hi boundaries of fit range
    //   startvalues[4]  reasonable start values for the fit
    //   parlimitslo[4]  lower parameter limits
    //   parlimitshi[4]  upper parameter limits
    //   fitparams[4]    returns the final fit parameters
    //   fiterrors[4]    returns the final fit errors
    //   ChiSqr          returns the chi square
    //   NDF             returns ndf

    Int_t i;
    Char_t FunName[100];

    sprintf ( FunName, "Fitfcn_%s", his -> GetName ( ) );

    TF1 *ffitold = static_cast<TF1*>( gROOT -> GetListOfFunctions ( ) -> FindObject ( FunName ) );
    if ( ffitold ) delete ffitold;

    TF1 *ffit = new TF1 ( FunName, langaufun, fitrange[0], fitrange[1], 4 );
    ffit -> SetParameters ( startvalues );
    ffit -> SetParNames ( "Width", "MP", "Area", "GSigma" );

    for ( i = 0; i < 4; i++ )
    {
	ffit -> SetParLimits ( i, parlimitslo[i], parlimitshi[i] );
    }

    // fit within specified range, use ParLimits, do not plot
    his -> Fit ( FunName, "RB0Q" );

    // obtain fit parameters
    ffit -> GetParameters ( fitparams );
    for ( i = 0; i < 4; i++ )
    {
	// obtain fit parameter errors
	fiterrors[i] = ffit -> GetParError ( i );
    }
    // obtain chi^2
    ChiSqr[0] = ffit -> GetChisquare ( );
    // obtain ndf
    NDF[0] = ffit -> GetNDF ( );
    // return fit function
    return ( ffit );
}

Int_t langaupro ( Double_t *params, Double_t &maxx, Double_t &FWHM )
{
    // Seaches for the location (x value) at the maximum of the 
    // Landau-Gaussian convolute and its full width at half-maximum.
    // The search is probably not very efficient, but it's a first try.

    Double_t p, x, fy, fxr, fxl;
    Double_t step;
    Double_t l, lold;
    Int_t i = 0;
    Int_t MAXCALLS = 10000;

    // Search for maximum
    p = params[1] - 0.1 * params[0];
    step = 0.05 * params[0];
    lold = -2.0;
    l = -1.0;

    while ( ( l != lold ) && ( i < MAXCALLS ) )
    {
	i++;
	lold = l;
	x = p + step;
	l = langaufun ( &x, params );

	if ( l < lold )
	{
	    step = -step / 10;
	}

	p += step;
    }

    if ( i == MAXCALLS )
    {
	return ( -1 );
    }

    maxx = x;
    fy = l / 2;

    // Search for right x location of fy
    p = maxx + params[0];
    step = params[0];
    lold = -2.0;
    l = -1e300;
    i = 0;

    while ( ( l != lold ) && ( i < MAXCALLS ) )
    {
	i++;
	lold = l;
	x = p + step;
	l = TMath::Abs ( langaufun ( &x, params ) - fy );

	if ( l > lold )
	{
	    step = -step / 10;
	}

	p += step;
    }

    if ( i == MAXCALLS )
    {
	return ( -2 );
    }

    fxr = x;

    // Search for left x location of fy
    p = maxx - 0.5 * params[0];
    step = -params[0];
    lold = -2.0;
    l = -1e300;
    i = 0;

    while ( ( l != lold ) && ( i < MAXCALLS ) )
    {
	i++;

	lold = l;
	x = p + step;
	l = TMath::Abs ( langaufun ( &x, params ) - fy );

	if ( l > lold )
	{
	    step = -step / 10;
	}

	p += step;
    }

    if ( i == MAXCALLS )
    {
	return ( -3 );
    }

    fxl = x;
    FWHM = fxr - fxl;
    return ( 0 );
}
