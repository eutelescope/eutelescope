/*
 * Rewritten and alibava features added by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 *
 *  Author  A.F.Zarnecki
 * (2007 University of Warsaw)
 *
 *  email:zarnecki@fuw.edu.pl
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined(USE_GEAR) && ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

// alibava includes
#include "AlibavaRunHeaderImpl.h"
#include "AlibavaEventImpl.h"
#include "ALIBAVA.h"

// eutelescope inlcudes
#include "EUTelFitTupleAlibava.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelSparseClusterImpl.h"


#include "EUTelSimpleVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

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
#include "gear/BField.h"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include "TGraph.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TProfile.h"

#define PI 3.14159265

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;
using namespace alibava;

// definition of static members mainly used to name histograms
std::string EUTelFitTupleAlibava::_FitTupleName  = "EUFit";

EUTelFitTupleAlibava::EUTelFitTupleAlibava ( ) : Processor ( "EUTelFitTupleAlibava" )
{
    // modify processor description
    _description = "Prepare n-tuple with track fit results. Do some DUT analysis." ;

    // register steering parameters: name, description, class-variable, default value, input collection first:
    registerInputCollection ( LCIO::TRACK, "InputCollectionName", "Name of the input track collection.", _inputColName, std::string ( "testfittracks" ) );

    registerInputCollection ( LCIO::TRACKERHIT, "InputDUTCollectionName", "Name of the input DUT hit collection.", _inputDUTColName, std::string ( "hit" ) );

    // processor parameters:
    registerProcessorParameter ( "AlibavaOutput", "Is there an Alibava DUT? If so, the alibava header data and reco data will also be written into the ntuple.", _alibava,  true );

    registerProcessorParameter ( "DistMax_X", "Maximum allowed distance in mm between fit and matched DUT hit in X.", _distMax_X, 0.1 );

    registerProcessorParameter ( "DistMax_Y", "Maximum allowed distance in mm between fit and matched DUT hit in Y.", _distMax_Y, 0.1 );

    std::vector < float > initAlign;
    initAlign.push_back ( 0.0 );
    initAlign.push_back ( 0.0 );
    initAlign.push_back ( 0.0 );
    initAlign.push_back ( 0.0 );
    initAlign.push_back ( 0.0 );
    initAlign.push_back ( 0.0 );

    registerProcessorParameter ( "DUTalignment", "These coordinates will be read from the gear file and used to position the DUT in the telescope and determine the track impact position. 0 to continue with a coordinate as it is in the input file collection, 1 to load the corresponding shift or rotation from the gear file. If the DUT is rotated, the coordinate should be switched on! Coordinates are: (X, Y, Z, ZY, ZX, XY).", _DUTalign, initAlign );

    registerProcessorParameter ( "ManualDUTid", "The telescope layer ID of the DUT.", _manualDUTid, 0 );

    registerProcessorParameter ( "MissingValue", "Value written to the ntuple for missing measurements.", _missingValue, -100.0 );

    registerProcessorParameter ( "UseManualDUT", "Set to true to use any DUT.", _useManualDUT, false );

    // optional parameters:
    registerOptionalParameter ( "AlibavaClusterCollectionName", "To compare matched hit positions against reconstructed clusters we need the name of the Alibava cluster collection.", _alibavaClusterCollectionName, std::string ( "alibava_clusters" ) );

    registerOptionalParameter ( "AlibavaRecoCollectionName", "If we are outputing the alibava data, this is the collection it is stored in.", _alibavaCollectionName, std::string ( "recodata_cmmd" ) );

    registerOptionalParameter ( "AlibavaUnfilteredCollectionName", "Alibava data from before a filtering step", _alibavaUnfilteredCollectionName, std::string ( "recodata_cmmd" ) );

    EVENT::StringVec _alignmentCollectionNameExamples;

    registerOptionalParameter ( "AlignmentCollectionName", "To reconstruct the local DUT coordinates, all alignment collection names which have been applied must be entered here. The order should correspond to their application, i.e.: alignment1 alignment2 etc.", _alignmentCollectionName, _alignmentCollectionNameExamples );

    registerOptionalParameter ( "AlignedHitPrefix", "The hit collection name of previously applied alignments, e.g. 'AlignedHit'", _alignedHitPrefix, std::string ( "AlignedHit" ) );

    registerOptionalParameter ( "DoAlignment","This processor can provide additional alignment for the DUT from its residuals. This is the flag to switch this on.", _doAlignment, false );

    registerOptionalParameter ( "DoGamma","If the DoAlignment flag is set, set this to true to provide alignment in XY-direction (gamma). False will provide alignment in ZY-direction (alpha), if x is the unsensitive axis or in ZX-direction (beta), if y is the unsensitive axis.", _doGamma, true );

    std::vector < float > initFiducial;
    initFiducial.push_back ( -99.0 );
    initFiducial.push_back ( 99.0 );
    initFiducial.push_back ( -99.0 );
    initFiducial.push_back ( 99.0 );

    registerOptionalParameter ( "Fiducial", "Additional residual histograms are created with these fiducial cuts: Xmin, Xmax, Ymin, Ymax. All are in global coordinates and in mm.", _fiducialcut, initFiducial );

    registerOptionalParameter ( "MaxTracksPerEvent", "The maximum number of tracks considered per event.", _tracklimit , 99 );

    registerOptionalParameter ( "OriginalHitCollection", "The original hit collection name before any alignment is added.", _originalCollectionName, std::string ( "PreAlignedHit2" ) );

    registerOptionalParameter ( "OriginalPreAlignment", "Any pre-Prealignment applied to the original hit collection.", _originalPreAlignment, std::string ( "prealignment" ) );

    registerOptionalParameter ( "OutputAlignment", "If we want to do our own DUT alignment here, this is the filename it will be stored in.", _outputalignment, std::string ( "hackalignment.slcio" ) );

    registerOptionalParameter ( "Polarity", "The polarity of our DUT. -1 for p-type, 1 for n-type", _polarity , -1 );

    EVENT::StringVec _pre_alignmentCollectionNameExamples;

    registerOptionalParameter ( "PreAlignmentCollectionName", "To reconstruct the local DUT coordinates, all prealignment collection names which have been applied must be entered here. The order should correspond to their application, i.e.: prealignment1 prealignment2 etc.", _pre_alignmentCollectionName, _pre_alignmentCollectionNameExamples );

    registerOptionalParameter ( "PreAlignedHitPrefix", "The hit collection name of previously applied prealignments, e.g. 'PreAlignedHit'", _prealignedHitPrefix, std::string ( "PreAlignedHit" ) );

    registerOptionalParameter ( "RecheckDealignment", "If this flag is set to true, an aligned hit from the DUT will be dealigned with the loaded alignment and prealignment collections and compared to the original DUT hit. This is only to check the correct loading and application of DUT alignment constants. DEBUG2 or lower should be set", _checkdealignment, false );

    registerOptionalParameter ( "ReferenceCollectionName", "The name of the reference collection. This is needed to get the correct Z DUT position after alignments.", _referencecollectionname, std::string ( "referenceHit" ) );

    registerOptionalParameter ( "UnAlignedHitPrefix", "The hit collection name without any alignment, e.g. 'hit'", _unalignedcollectionname, std::string ( "hit" ) );

    registerOptionalParameter ( "UnsensitiveAxis", "If we have a strip sensor DUT, set the unsensitive axis here (x or y). This will be used for the additional alignment (if you switch it on).", _unsensitiveaxis, std::string ( "x" ) );

    registerOptionalParameter ( "UseOriginalPreAlignment", "Subtract any pre-Prealignment applied to the original hit collection?", _useOriginalPreAlignment, true );

    registerOptionalParameter ( "UseTrackFit", "Use an existing track fit to estimate impact position? If false, a guess will be used.", _useTrackFit, false );

    registerOptionalParameter ( "UseTriplet", "Use triplet from upstream telescope to estimate impact position? If false, result of UseTrackFit will be used", _useTriplet, false );

}

void EUTelFitTupleAlibava::init ( )
{

    const gear::BField& Bfield = Global::GEAR -> getBField ( );
    gear::Vector3D vectorGlobal ( 0.1, 0.1, 0.1 );
    const double Bx = ( Bfield.at ( vectorGlobal ) .x ( ) );
    const double By = ( Bfield.at ( vectorGlobal ) .y ( ) );
    const double Bz = ( Bfield.at ( vectorGlobal ) .z ( ) );
    TVector3 B ( Bx, By, Bz );
    if ( fabs ( Bx ) > 1e-6 || fabs ( By ) > 1e-6 || fabs ( By ) > 1e-6 )
    {
	streamlog_out ( MESSAGE4 ) << "Running WITH magnetic field! " << Bx << " Bx, " << By << " By, " << Bz << " Bz!" <<  endl;
	streamlog_out ( MESSAGE4 ) << "The residual cuts will be multiplied by 10!" << endl;
	_distMax_X = _distMax_X * 10.0;
	_distMax_Y = _distMax_Y * 10.0;
    }
    else
    {
	streamlog_out ( MESSAGE4 ) << "Running WITHOUT magnetic field!" << endl;
    }

    // usually a good idea to
    printParameters ( );

    _nRun = 0;
    _nEvt = 0;
    _matchedhits = 0;

    // check if the GEAR manager pointer is not null!
    if ( Global::GEAR == nullptr )
    {
	streamlog_out ( ERROR5 ) << "The GearMgr is not available, for an unknown reason." << endl;
	exit ( -1 );
    }

    // Read geometry information from GEAR
    streamlog_out ( MESSAGE5 ) <<  "Reading telescope geometry description from GEAR " << endl;

    _siPlanesParameters  = const_cast < gear::SiPlanesParameters* > ( & ( Global::GEAR -> getSiPlanesParameters ( ) ) );
    _siPlanesLayerLayout = const_cast < gear::SiPlanesLayerLayout* > ( & ( _siPlanesParameters -> getSiPlanesLayerLayout ( ) ) );

    // Take all layers defined in GEAR geometry
    _nTelPlanes = _siPlanesLayerLayout -> getNLayers ( );

    // Check for DUT
    if ( _siPlanesParameters -> getSiPlanesType ( ) == _siPlanesParameters -> TelescopeWithDUT )
    {
	_iDUT = _nTelPlanes;
	_nTelPlanes++;
    }
    else
    {
	_iDUT = -1;
    }

    // Read position in Z (for sorting)
    _planeSort = new int[_nTelPlanes];
    _planePosition = new double[_nTelPlanes];

    for ( int ipl = 0; ipl < _siPlanesLayerLayout -> getNLayers ( ); ipl++ )
    {
	_planePosition[ipl] = _siPlanesLayerLayout -> getLayerPositionZ ( ipl );
	_planeSort[ipl] = ipl;
    }

    if ( _iDUT > 0 )
    {
	_planePosition[_iDUT] = _siPlanesLayerLayout -> getDUTPositionZ ( );
	_planeSort[_iDUT] = _iDUT;
    }

    // Binary sorting
    bool sorted;
    do
    {
	sorted = false;
	for ( int iz = 0; iz < _nTelPlanes - 1; iz++ )
	{
	    if ( _planePosition[iz] > _planePosition[iz + 1] )
	    {
		double _posZ = _planePosition[iz];
		_planePosition[iz] = _planePosition[iz + 1];
		_planePosition[iz + 1] = _posZ;
		int _idZ = _planeSort[iz];
		_planeSort[iz] = _planeSort[iz + 1];
		_planeSort[iz + 1] = _idZ;
		sorted = true;
	    }
	}
    }
    while ( sorted );

    // Book local geometry arrays
    _planeID = new int[_nTelPlanes];
    _isActive = new bool[_nTelPlanes];

    // Fill remaining layer parameters
    for ( int iz = 0; iz < _nTelPlanes; iz++ )
    {
	int ipl=_planeSort[iz];
	double resolution;

	if ( ipl != _iDUT )
	{
	    _planeID[iz] = _siPlanesLayerLayout -> getID ( ipl );
	    resolution = _siPlanesLayerLayout -> getSensitiveResolution ( ipl );
	}
	else
	{
	    _planeID[iz] = _siPlanesLayerLayout -> getDUTID ( );
	    resolution = _siPlanesLayerLayout -> getDUTSensitiveResolution ( );
	}
	_isActive[iz] = ( resolution > 0 );
    }

    // Get new DUT position (after sorting)
    for ( int iz = 0; iz < _nTelPlanes; iz++ )
    {
	if ( _planeSort[iz] == _iDUT )
	{
	    _iDUT = iz;
	    break;
	}
    }

    // DUT position can be changed by processor parameter
    if ( _useManualDUT )
    {
	bool _manualOK = false;
	for ( int iz = 0; iz < _nTelPlanes; iz++ )
	{
	    if ( _planeID[iz] == _manualDUTid )
	    {
		_iDUT = iz;
		_manualOK = true;
	    }
	}

	if ( !_manualOK )
	{
	    streamlog_out ( ERROR5 ) <<  "Manual DUT flag set, layer not found ID = " << _manualDUTid << " . Program will terminate! Correct geometry description!" << endl;
	    exit ( -1 );
	}
    }

    _zDUT=_planePosition[_iDUT];

    // Print out geometry information
    streamlog_out ( MESSAGE5 ) << "Telescope configuration with " << _nTelPlanes << " planes" << endl;

    for ( int ipl = 0; ipl < _nTelPlanes; ipl++ )
    {
	stringstream st;
	if ( ipl == _iDUT )
	{
	    st << "D.U.T.  plane";
	}
	else
	{
	    if ( _isActive[ipl] )
	    {
		st << "Active  plane";
	    }
	    else
	    {
		st << "Passive plane";
	    }
	}
	st << "  ID = " << _planeID[ipl] << "  at Z [mm] = " << _planePosition[ipl];
	streamlog_out ( MESSAGE5 ) <<  st.str ( ) << endl;
    }

    // Allocate arrays for track fitting
    _isMeasured = new bool[_nTelPlanes];
    _isFitted = new bool[_nTelPlanes];
    _measuredX = new double[_nTelPlanes];
    _measuredY = new double[_nTelPlanes];
    _measuredZ = new double[_nTelPlanes];
    _measuredQ = new double[_nTelPlanes];
    _fittedX = new double[_nTelPlanes];
    _fittedY = new double[_nTelPlanes];
    _fittedZ = new double[_nTelPlanes];
    _telescopeX = new double[_nTelPlanes];
    _telescopeY = new double[_nTelPlanes];

    // Book histograms
    bookHistos ( );

    // Take shifts and roations from gear, if so wanted
    // caveat: not all implemented, see bottom.
    // assumption: gear has mm and deg
    if ( _DUTalign.at ( 0 ) != 0 )
    {
	_DUTalign.at ( 0 ) = _siPlanesLayerLayout -> getLayerPositionX ( _iDUT );
	streamlog_out ( MESSAGE5 ) << "Set DUT X shift according to gear file: " << _DUTalign.at ( 0 ) << endl;
    }
    if ( _DUTalign.at ( 1 ) != 0 )
    {
	_DUTalign.at ( 1 ) = _siPlanesLayerLayout -> getLayerPositionY ( _iDUT );
	streamlog_out ( MESSAGE5 ) << "Set DUT Y shift according to gear file: " << _DUTalign.at ( 1 ) << endl;
    }
    if ( _DUTalign.at ( 2 ) != 0 )
    {
	_DUTalign.at ( 2 ) = _siPlanesLayerLayout -> getLayerPositionZ ( _iDUT );
	streamlog_out ( MESSAGE5 ) << "Set DUT Z shift according to gear file: " << _DUTalign.at ( 2 ) << endl;
    }
    if ( _DUTalign.at ( 3 ) != 0 )
    {
	_DUTalign.at ( 3 ) = _siPlanesLayerLayout -> getLayerRotationZY ( _iDUT );
	streamlog_out ( MESSAGE5 ) << "Set DUT ZY rotation according to gear file: " << _DUTalign.at ( 3 ) << endl;
    }
    if ( _DUTalign.at ( 4 ) != 0 )
    {
	_DUTalign.at ( 4 ) = _siPlanesLayerLayout -> getLayerRotationZX ( _iDUT );
	streamlog_out ( MESSAGE5 ) << "Set DUT ZX rotation according to gear file: " << _DUTalign.at ( 4 ) << endl;
    }
    if ( _DUTalign.at ( 5 ) != 0 )
    {
	_DUTalign.at ( 5 ) = _siPlanesLayerLayout -> getLayerRotationXY ( _iDUT );
	streamlog_out ( MESSAGE5 ) << "Set DUT XY rotation according to gear file: " << _DUTalign.at ( 5 ) << endl;
    }

    // load alignment from file
    _alignmentloaded = false;
    _pre_alignmentloaded = false;
    _referenceloaded = false;

    // allocate arrays for alignment constants
    _dut_align_x = new double[_alignmentCollectionName.size ( )];
    _dut_align_y = new double[_alignmentCollectionName.size ( )];
    _dut_align_z = new double[_alignmentCollectionName.size ( )];
    _dut_align_a = new double[_alignmentCollectionName.size ( )];
    _dut_align_b = new double[_alignmentCollectionName.size ( )];
    _dut_align_c = new double[_alignmentCollectionName.size ( )];
    _dut_align_x_error = new double[_alignmentCollectionName.size ( )];
    _dut_align_y_error = new double[_alignmentCollectionName.size ( )];
    _dut_align_z_error = new double[_alignmentCollectionName.size ( )];
    _dut_align_a_error = new double[_alignmentCollectionName.size ( )];
    _dut_align_b_error = new double[_alignmentCollectionName.size ( )];
    _dut_align_c_error = new double[_alignmentCollectionName.size ( )];
    _dut_pre_align_x = new double[_pre_alignmentCollectionName.size ( )];
    _dut_pre_align_y = new double[_pre_alignmentCollectionName.size ( )];
    _dut_pre_align_z = new double[_pre_alignmentCollectionName.size ( )];
    _dut_pre_align_a = new double[_pre_alignmentCollectionName.size ( )];
    _dut_pre_align_b = new double[_pre_alignmentCollectionName.size ( )];
    _dut_pre_align_c = new double[_pre_alignmentCollectionName.size ( )];

    // get the DUT parameters from gear
    _pitchx = _siPlanesLayerLayout -> getSensitivePitchX ( _iDUT );
    _pitchy = _siPlanesLayerLayout -> getSensitivePitchY ( _iDUT );
    _pixelx = _siPlanesLayerLayout -> getSensitiveNpixelX ( _iDUT );
    _pixely = _siPlanesLayerLayout -> getSensitiveNpixelY ( _iDUT );

}

void EUTelFitTupleAlibava::getReference ( )
{

    // get the DUT z position from the reference file
    // there might be a difference between this z and the one loaded from the gear file

    _x_refhit = 0.0;
    _y_refhit = 0.0;
    _z_refhit = 0.0;
    _a_refhit = 0.0;
    _b_refhit = 0.0;
    _c_refhit = 0.0;

    _x_refhit = geo::gGeometry ( ) .getPlaneXPosition ( _manualDUTid );
    _y_refhit = geo::gGeometry ( ) .getPlaneYPosition ( _manualDUTid );
    _z_refhit = geo::gGeometry ( ) .getPlaneZPosition ( _manualDUTid ) + 0.5 * geo::gGeometry ( ) .getPlaneZSize ( _manualDUTid );

    double refVec[3];
    refVec[0] = 0.0;
    refVec[1] = 0.0;
    refVec[2] = 1.0;

    double gRotation[3] = { 0.0, 0.0, 0.0 };

    gRotation[0] = geo::gGeometry ( ) .getPlaneZRotationDegrees ( _manualDUTid );
    gRotation[1] = geo::gGeometry ( ) .getPlaneYRotationDegrees ( _manualDUTid );
    gRotation[2] = geo::gGeometry ( ) .getPlaneXRotationDegrees ( _manualDUTid );
    gRotation[0] = gRotation[0] * 3.1415926 / 180.0;
    gRotation[1] = gRotation[1] * 3.1415926 / 180.0;
    gRotation[2] = gRotation[2] * 3.1415926 / 180.0;

    TVector3 _RotatedVector ( refVec[0], refVec[1], refVec[2] );
    TVector3 _Xaxis ( 1.0, 0.0, 0.0 );
    TVector3 _Yaxis ( 0.0, 1.0, 0.0 );
    TVector3 _Zaxis ( 0.0, 0.0, 1.0 );

    if ( TMath::Abs ( gRotation[2] ) > 1e-6 )
    {
	_RotatedVector.Rotate ( gRotation[2], _Xaxis );
    }
    if ( TMath::Abs ( gRotation[1] ) > 1e-6 )
    {
	_RotatedVector.Rotate ( gRotation[1], _Yaxis );
    }
    if ( TMath::Abs ( gRotation[0] ) > 1e-6 )
    {
	_RotatedVector.Rotate ( gRotation[0], _Zaxis );
    }

    _a_refhit = _RotatedVector[0];
    _b_refhit = _RotatedVector[1];
    _c_refhit = _RotatedVector[2];

    streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
    streamlog_out ( DEBUG7 ) << "Reference Collection loaded:" << endl;
    streamlog_out ( DEBUG7 ) << "DUT reference X is:     " << _x_refhit << endl;
    streamlog_out ( DEBUG7 ) << "DUT reference Y is:     " << _y_refhit << endl;
    streamlog_out ( DEBUG7 ) << "DUT reference Z is:     " << _z_refhit << endl;
    streamlog_out ( DEBUG7 ) << "DUT reference Alpha is: " << _a_refhit << endl;
    streamlog_out ( DEBUG7 ) << "DUT reference Beta is:  " << _b_refhit << endl;
    streamlog_out ( DEBUG7 ) << "DUT reference Gamma is: " << _c_refhit << endl;
}

void EUTelFitTupleAlibava::getAlignment ( LCEvent * event )
{

    // Get the alignment from file
    // assumes rotations in rad!

    LCCollectionVec * alignmentCollection;

    for ( size_t i = 0; i < _alignmentCollectionName.size ( ); i++ )
    {

	try
	{
	    alignmentCollection = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( _alignmentCollectionName[i] ) );

	    for ( int j = 0; j < _nTelPlanes; j++ )
	    {
		EUTelAlignmentConstant * alignment = static_cast < EUTelAlignmentConstant * > ( alignmentCollection -> getElementAt ( j ) );

		if ( alignment -> getSensorID ( ) == _manualDUTid )
		{
		    _dut_align_x[i] = alignment -> getXOffset ( );
		    _dut_align_y[i] = alignment -> getYOffset ( );
		    _dut_align_z[i] = alignment -> getZOffset ( );
		    _dut_align_a[i] = alignment -> getAlpha ( );
		    _dut_align_b[i] = alignment -> getBeta ( );
		    _dut_align_c[i] = alignment -> getGamma ( );
		    _dut_align_x_error[i] = alignment -> getXOffsetError ( );
		    _dut_align_y_error[i] = alignment -> getYOffsetError ( );
		    _dut_align_z_error[i] = alignment -> getZOffsetError ( );
		    _dut_align_a_error[i] = alignment -> getAlphaError ( );
		    _dut_align_b_error[i] = alignment -> getBetaError ( );
		    _dut_align_c_error[i] = alignment -> getGammaError ( );
		}
	    }
	}
	catch ( lcio::DataNotAvailableException& e )
	{
	    streamlog_out ( ERROR5 ) << "No alignment collection with name " << _alignmentCollectionName[i] << " !" << endl;
	}

	streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
	streamlog_out ( DEBUG7 ) << "Alignment loading, iteration " << i+1 << ":" << endl;
	streamlog_out ( DEBUG7 ) << "DUT alignment: X is:     " << _dut_align_x[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT alignment: Y is:     " << _dut_align_y[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT alignment: Z is:     " << _dut_align_z[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT alignment: Alpha is: " << _dut_align_a[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT alignment: Beta is:  " << _dut_align_b[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT alignment: Gamma is: " << _dut_align_c[i] << endl;
    }
}

void EUTelFitTupleAlibava::getPreAlignment ( LCEvent * event )
{

    // Get the alignment from file
    // assumes rotations in rad!

    LCCollectionVec * pre_alignmentCollection;

    for ( size_t i = 0; i < _pre_alignmentCollectionName.size ( ); i++ )
    {
	try
	{
	    pre_alignmentCollection = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( _pre_alignmentCollectionName[i] ) );

	    for ( int j = 0; j < _nTelPlanes; j++ )
	    {
		EUTelAlignmentConstant * pre_alignment = static_cast < EUTelAlignmentConstant * > ( pre_alignmentCollection -> getElementAt ( j ) );

		if ( pre_alignment -> getSensorID ( ) == _manualDUTid )
		{
		    _dut_pre_align_x[i] = pre_alignment -> getXOffset ( );
		    _dut_pre_align_y[i] = pre_alignment -> getYOffset ( );
		    _dut_pre_align_z[i] = pre_alignment -> getZOffset ( );
		    _dut_pre_align_a[i] = pre_alignment -> getAlpha ( );
		    _dut_pre_align_b[i] = pre_alignment -> getBeta ( );
		    _dut_pre_align_c[i] = pre_alignment -> getGamma ( );
		}
	    }
	}
	catch ( lcio::DataNotAvailableException& e )
	{
	    streamlog_out ( ERROR5 ) << "No pre_alignment collection with name " << _pre_alignmentCollectionName[i] << " !" << endl;
	}

	streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
	streamlog_out ( DEBUG7 ) << "Prealignment loading, iteration " << i+1 << ":" << endl;
	streamlog_out ( DEBUG7 ) << "DUT pre_alignment: X is:     " << _dut_pre_align_x[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT pre_alignment: Y is:     " << _dut_pre_align_y[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT pre_alignment: Z is:     " << _dut_pre_align_z[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT pre_alignment: Alpha is: " << _dut_pre_align_a[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT pre_alignment: Beta is:  " << _dut_pre_align_b[i] << endl;
	streamlog_out ( DEBUG7 ) << "DUT pre_alignment: Gamma is: " << _dut_pre_align_c[i] << endl;
    }

    // get original prealignment
    if ( _useOriginalPreAlignment == true )
    {

	_dut_original_pre_align_x = 0.0;
	_dut_original_pre_align_y = 0.0;
	_dut_original_pre_align_z = 0.0;
	_dut_original_pre_align_a = 0.0;
	_dut_original_pre_align_b = 0.0;
	_dut_original_pre_align_c = 0.0;

	try
	{
	    pre_alignmentCollection = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( _originalPreAlignment ) );

	    for ( int j = 0; j < _nTelPlanes; j++ )
	    {
		EUTelAlignmentConstant * pre_alignment = static_cast < EUTelAlignmentConstant * > ( pre_alignmentCollection -> getElementAt ( j ) );

		if ( pre_alignment -> getSensorID ( ) == _manualDUTid )
		{
		    _dut_original_pre_align_x = pre_alignment -> getXOffset ( );
		    _dut_original_pre_align_y = pre_alignment -> getYOffset ( );
		    _dut_original_pre_align_z = pre_alignment -> getZOffset ( );
		    _dut_original_pre_align_a = pre_alignment -> getAlpha ( );
		    _dut_original_pre_align_b = pre_alignment -> getBeta ( );
		    _dut_original_pre_align_c = pre_alignment -> getGamma ( );
		}
	    }
	}
	catch ( lcio::DataNotAvailableException& e )
	{
	    streamlog_out ( ERROR5 ) << "No original pre_alignment collection with name " << _originalPreAlignment << " !" << endl;
	}

	streamlog_out ( DEBUG7 ) << "**************************************************" << endl;
	streamlog_out ( DEBUG7 ) << "Original Prealignment loading:" << endl;
	streamlog_out ( DEBUG7 ) << "DUT original pre_alignment: X is:     " << _dut_original_pre_align_x << endl;
	streamlog_out ( DEBUG7 ) << "DUT original pre_alignment: Y is:     " << _dut_original_pre_align_y << endl;
	streamlog_out ( DEBUG7 ) << "DUT original pre_alignment: Z is:     " << _dut_original_pre_align_z << endl;
	streamlog_out ( DEBUG7 ) << "DUT original pre_alignment: Alpha is: " << _dut_original_pre_align_a << endl;
	streamlog_out ( DEBUG7 ) << "DUT original pre_alignment: Beta is:  " << _dut_original_pre_align_b << endl;
	streamlog_out ( DEBUG7 ) << "DUT original pre_alignment: Gamma is: " << _dut_original_pre_align_c << endl;
    }
}

void EUTelFitTupleAlibava::processRunHeader ( LCRunHeader* runHeader )
{
    auto eutelHeader = std::make_unique < EUTelRunHeaderImpl > ( runHeader );
    eutelHeader -> addProcessor ( type ( ) );
    _nRun++;

    // Decode and print out Run Header information - just a check
    _runNr = runHeader -> getRunNumber ( );
    streamlog_out ( MESSAGE5 ) << "Processing run header " << _nRun << ", run nr " << _runNr << endl;

    const std::string detectorName = runHeader -> getDetectorName ( );
    const std::string detectorDescription = runHeader -> getDescription ( );
    const std::vector < std::string > * subDets = runHeader -> getActiveSubdetectors ( );

    streamlog_out ( MESSAGE5 ) << detectorName << " : " << detectorDescription << endl;

    int nDet = subDets -> size ( );
    if ( nDet )
    {
	streamlog_out ( MESSAGE5 ) << nDet << " subdetectors defined :" << endl;
    }
    for ( int idet = 0; idet < nDet; idet++ )
    {
	streamlog_out ( MESSAGE5 ) << idet + 1 << " : " << subDets -> at ( idet ) << endl;
    }

}

void EUTelFitTupleAlibava::processEvent ( LCEvent * event )
{

    // load the alignment only once
    if ( _alignmentloaded == false )
    {
	getAlignment ( event );
	_alignmentloaded = true;
	streamlog_out ( MESSAGE5 ) <<  "Alignment loaded..." << endl;
	streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
	streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
    }

    // also load prealignment only once
    if ( _pre_alignmentloaded == false )
    {
	getPreAlignment ( event );
	_pre_alignmentloaded = true;
	streamlog_out ( MESSAGE5 ) <<  "Prealignment loaded..." << endl;
	streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
	streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
    }

    // also load reference collection only once
    if ( _referenceloaded == false )
    {
	getReference ( );
	_referenceloaded = true;
	streamlog_out ( MESSAGE5 ) <<  "Reference loaded..." << endl;
	streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
	streamlog_out ( MESSAGE5 ) <<  "**************************************************" << endl;
    }

    // stop if this is the last event
    EUTelEventImpl * euEvent = static_cast < EUTelEventImpl* > ( event );
    if ( euEvent -> getEventType ( ) == kEORE )
    {
	streamlog_out ( DEBUG5 ) <<  "EORE found: nothing else to do." << endl;
	return;
    }

    // the number of this event
    _nEvt ++;
    _evtNr = event -> getEventNumber ( );
    streamlog_out ( DEBUG2 ) << "Processing record " << _nEvt << " == event " << _evtNr << endl;

    // the track input collection of the event
    LCCollection* col;
    try
    {
	col = event -> getCollection ( _inputColName );
    }
    catch ( lcio::DataNotAvailableException& e )
    {
	streamlog_out ( DEBUG1 ) << "Not able to get track input collection " << _inputColName << " from event " << event -> getEventNumber ( ) << " in run " << event -> getRunNumber ( ) << endl;
	return;
    }

    // the hit input collection of the event
    LCCollection* hitcol = nullptr;
    bool _DUTok = true;
    try
    {
	hitcol = event -> getCollection ( _inputDUTColName );
    }
    catch ( lcio::DataNotAvailableException& e )
    {
	// this can happen more often, as there could be different DAQs or event counts, etc.
	// _DUTok is the flag for this
	streamlog_out ( DEBUG3 ) << "No able to get DUT hit input collection " << _inputDUTColName << " from event " << event -> getEventNumber ( ) << " in run " << event -> getRunNumber ( ) << endl;
	_DUTok = false;
    }

    // Loop over tracks in input collections
    int nTrack = col -> getNumberOfElements ( );
    streamlog_out ( DEBUG1 ) << "Total of " << nTrack << " tracks in input collection " << _inputColName << " !" << endl;

    int nInputHits = 0;

    if ( _DUTok )
    {
	nInputHits = hitcol -> getNumberOfElements ( );
    }
    streamlog_out ( DEBUG1 ) << "Total of " << nInputHits << " track hits in input collection " << _inputDUTColName << " !" << endl;

    // a matching DUT hit has not been found yet...
    _foundDUTHit = false;

    // Track limit per event
    int maxTracks = 0;
    // only if we have a track at all
    if ( nTrack > 0 )
    {
	if ( nTrack > _tracklimit )
	{
	    maxTracks = _tracklimit;
	}
	else if ( nTrack <= _tracklimit )
	{
	    maxTracks = nTrack;
	}
    }


    // first, loop over the tracks in the event
    for ( int itrack = 0; itrack < maxTracks; itrack++ )
    {
	Track * fittrack = dynamic_cast < Track* > ( col -> getElementAt ( itrack ) );

	// hit list assigned to track
	std::vector < EVENT::TrackerHit* > trackhits = fittrack -> getTrackerHits ( );

	// copy hits assigned to the track to a local table
	// assign hits to sensor planes
	int nHit = trackhits.size ( );
	double trackchi2 = 0.0;
	trackchi2 = fittrack -> getChi2 ( );
	double trackndf = 0.0;
	trackndf = fittrack -> getNdf ( );

	streamlog_out ( DEBUG7 ) << "Track " << itrack << " of " << nTrack << " (limit: " << maxTracks << " ) with " << nHit << " hits, Chi2 = " << trackchi2 << "/" << trackndf << " !" << endl;

	// clear plane tables
	for ( int ipl = 0; ipl < _nTelPlanes; ipl++ )
	{
	    _isMeasured[ipl] = false;
	    _isFitted[ipl] = false;

	    _measuredX[ipl] = _missingValue;
	    _measuredY[ipl] = _missingValue;
	    _measuredZ[ipl] = _missingValue;
	    _measuredQ[ipl] = _missingValue;

	    _fittedX[ipl] = _missingValue;
	    _fittedY[ipl] = _missingValue;
	    _fittedZ[ipl] = _missingValue;
	}

	// loop over hits and fill hit tables

	CellIDDecoder < TrackerHit >  hitCellDecoder ( EUTELESCOPE::HITENCODING );

	for ( int ihit = 0; ihit < nHit; ihit++ )
	{
	    TrackerHit * meshit = trackhits.at ( ihit );

	    // hit position
	    const double * pos = meshit -> getPosition ( );

	    // we find the plane number of the hit by looking at the Z position...
	    // FIXME this might be a problem in rotated sensors
	    double distMin = 5.0;
	    int hitPlane = -1;

	    for ( int ipl = 0; ipl < _nTelPlanes; ipl++ )
	    {
		double dist =  pos[2] - _planePosition[ipl];

		if ( dist * dist < distMin * distMin )
		{
		    hitPlane = ipl;
		    distMin = dist;
		}
	    }

	    // ignore hits not matched to any plane
	    if ( hitPlane < 0 )
	    {
		streamlog_out ( DEBUG3 ) << "Hit outside telescope plane at z [mm] = "  << pos[2] << " ! Maybe the gear settings for the DUT are incorrect!" << endl;
		hitPlane = _iDUT;
		continue;
	    }

	    if ( ( hitCellDecoder ( meshit ) ["properties"] & kFittedHit ) == 0 )
	    {
		// these are measured hits and have rawdata assigned
		_isMeasured[hitPlane] = true;
		_measuredX[hitPlane] = pos[0];
		_measuredY[hitPlane] = pos[1];
		_measuredZ[hitPlane] = pos[2];

		// get the cluster charge
		_measuredQ[hitPlane] = 0.0;

		EVENT::LCObjectVec rawdata =  meshit -> getRawHits ( );
		if ( rawdata.size ( ) > 0 && rawdata.at ( 0 ) != nullptr )
		{
		    EUTelVirtualCluster * cluster = new EUTelFFClusterImpl ( static_cast < TrackerDataImpl* > ( rawdata.at ( 0 ) ) );
		    _measuredQ[hitPlane] = cluster -> getTotalCharge ( );
		}

		streamlog_out ( DEBUG3 ) << "Measured hit in plane " << hitPlane << " at  X = " << pos[0] << ", Y = " << pos[1] << ", Q = " << _measuredQ[hitPlane] << endl;

	    }
	    else
	    {

		// fitted hits don't have rawdata, as they come from a track fit and not from a cluster
		_isFitted[hitPlane] = true;
		_fittedX[hitPlane] = pos[0];
		_fittedY[hitPlane] = pos[1];
		_fittedZ[hitPlane] = pos[2];

		streamlog_out ( DEBUG3 ) << "Fitted  hit  in plane " << hitPlane << " at  X = " << pos[0] << ", Y = " << pos[1] << ", Z = " << pos[2] << endl;

	    }

	    _telescopeX[hitPlane] = _fittedX[hitPlane];
	    _telescopeY[hitPlane] = _fittedY[hitPlane];

	} // end hit loop

	// fill the n-tuple
	// first 0-4 rows with event and track information
	_FitTuple -> fill ( 0, _nEvt );
	_FitTuple -> fill ( 1, _runNr );
	_FitTuple -> fill ( 2, _evtNr );
	_FitTuple -> fill ( 3, fittrack -> getNdf ( ) );
	_FitTuple -> fill ( 4, fittrack -> getChi2 ( ) );

	// now with the information of the hits
	int icol = 5;
	for ( int ipl = 0; ipl < _nTelPlanes; ipl++ )
	{
	    _FitTuple -> fill ( icol++, _measuredX[ipl] );
	    _FitTuple -> fill ( icol++, _measuredY[ipl] );
	    _FitTuple -> fill ( icol++, _measuredZ[ipl] );
	    _FitTuple -> fill ( icol++, _measuredQ[ipl] );
	    _FitTuple -> fill ( icol++, _fittedX[ipl] );
	    _FitTuple -> fill ( icol++, _fittedY[ipl] );
	    _FitTuple -> fill ( icol++, _fittedZ[ipl] );
	}

	// clear dut variables
	double dutHitR = _missingValue;
	double dutHitQ = _missingValue;

	// if there is no DUT in our setup, then we're done here and just write out _missingValue to preserve data structure in the ntuple
	if ( _iDUT == -1 )
	{
	    // track xyz in global
	    _FitTuple -> fill ( icol++, _missingValue );
	    _FitTuple -> fill ( icol++, _missingValue );
	    _FitTuple -> fill ( icol++, _missingValue );

	    // track xyz in local
	    _FitTuple -> fill ( icol++, _missingValue );
	    _FitTuple -> fill ( icol++, _missingValue );
	    _FitTuple -> fill ( icol++, _missingValue );

	    // hit xyz in global
	    _FitTuple -> fill ( icol++, _missingValue );
	    _FitTuple -> fill ( icol++, _missingValue );
	    _FitTuple -> fill ( icol++, _missingValue );

	    // hit information
	    _FitTuple -> fill ( icol++, dutHitR );
	    _FitTuple -> fill ( icol++, dutHitQ );
	}

	// if there is a DUT, we want to write it to the ntuple and do some analysis
	// require that there is no DUT hit in the track fit...
	// && _measuredX[_iDUT] == _missingValue
	if ( _iDUT != -1 )
	{

	    // define axis for rotations:
	    TVector3 _Xaxis ( 1.0, 0.0, 0.0 );
	    TVector3 _Yaxis ( 0.0, 1.0, 0.0 );
	    TVector3 _Zaxis ( 0.0, 0.0, 1.0 );

	    // the global coordinates of the track hit on the DUT which we want to find
	    double dutTrackX_global = 0.0;
	    double dutTrackY_global = 0.0;
	    double dutTrackZ_global = 0.0;

	    // some output
	    streamlog_out ( DEBUG6 ) << " " << endl;
	    streamlog_out ( DEBUG6 ) << "===============================================================================================================================================" << endl;
	    streamlog_out ( DEBUG6 ) << "Starting DUT magic in event " << _evtNr << " - dumpevent " << _evtNr + 1 << " !" << endl;
	    streamlog_out ( DEBUG6 ) << "===============================================================================================================================================" << endl;

	    // Several possibilities for guessing the track point on the DUT:

	    // ****************** //
	    // 0: straight line interpolation between track fits from DUT + 1 and DUT - 1
	    // very rough guess at the expected hit position on the DUT... do not use for real analysis
	    // ****************** //

	    //double dutTrackX =  _fittedX[_iDUT-1] + (_fittedX[_iDUT+1] - _fittedX[_iDUT-1]) / (_planePosition[_iDUT+1] - _planePosition[_iDUT-1]) * (_zDUT - _planePosition[_iDUT-1]);
	    //double dutTrackY =  _fittedY[_iDUT-1] + (_fittedY[_iDUT+1] - _fittedY[_iDUT-1]) / (_planePosition[_iDUT+1] - _planePosition[_iDUT-1]) * (_zDUT - _planePosition[_iDUT-1]);


	    // ****************** //
	    // 1: Plane and line intersection, with alignments
	    // this is an improvement of the following triplet method:
	    // - a vector is calculated from the hits in planes 0 and 3
	    // - the dut is (pre)aligned and rotated according to gear
	    // - the intersection of the vector and the (rotated and shifted) plane is calculated and defined as track impact point
	    // ****************** //

	    double planeintersectX = 0.0;
	    double planeintersectY = 0.0;
	    double planeintersectZ = 0.0;

	    // only calculate if mode is selected
	    if ( _useTrackFit == false && _useTriplet == false )
	    {
		streamlog_out ( DEBUG6 ) << "Using plane and vector intersection! " << endl;

		// the vector from preceeding planes
		double tr_x = _telescopeX[_iDUT - 1] - _telescopeX[_iDUT - 2];
		double tr_y = _telescopeY[_iDUT - 1] - _telescopeY[_iDUT - 2];
		double tr_z = _planePosition[_iDUT - 1] - _planePosition[_iDUT - 2];

		// DUT plane start vector should get all gear shifts, prealignment steps and alignment steps
		// the initial dut plane vector is it's z position:
		TVector3 dutplane ( 0.0, 0.0, 0.0 );

		// gear shifts get applied to start vector...
		TVector3 gearshift ( _DUTalign.at ( 0 ), _DUTalign.at ( 1 ), _DUTalign.at ( 2 ) );
		dutplane = dutplane + gearshift;

		// now the normal vector on the DUT plane
		// it gets rotations from gear and also the alignment rotations, as the normal vector now no longer only points along z:
		// as it points along z, it should be invariant to all rotations around gamma -> FIXME check this!
		// as a normal vector it is invariant against shifts, so the gear shifts and alignment shifts are not added

		// the gear rotations are passed from here
		double alpha = _DUTalign.at ( 3 );
		double beta = _DUTalign.at ( 4 );
		double gamma = _DUTalign.at ( 5 );

		// gamma should not be set
		if ( gamma != 0 )
		{
		    streamlog_out ( ERROR5 ) << "DUT GEAR alignment in gamma is switched on! This is not implemented!" << endl;
		}

		// the components of the plane normal vector rotated by gear angles
		// if all angles are zero, this should be ( 0 | 0 | 1 )
		double n_x = sin ( beta * PI / 180.0 );
		double n_y = -sin ( alpha * PI / 180.0 );
		double n_z = cos ( alpha * PI / 180.0 ) * cos ( beta * PI / 180.0 );
		TVector3 normalvector ( n_x, n_y, n_z );

		// now we have gear shifts and rotations for both plane start vector and plane normal vector
		// both need to be rotated with the alignment rotations
		// the start vector also needs the shifts, the normal vector should be invariant against this

		// the prealignment gets added...
		for ( unsigned int i = 0; i < _pre_alignmentCollectionName.size ( ); i++ )
		{
		    streamlog_out ( DEBUG2 ) << "Performing prealignment on DUT plane vector, step " << i + 1 << " !" << endl;

		    // add rotation -> negative rotation
		    dutplane.Rotate ( - _dut_pre_align_a[i], _Xaxis );
		    dutplane.Rotate ( - _dut_pre_align_b[i], _Yaxis );
		    dutplane.Rotate ( - _dut_pre_align_c[i], _Zaxis );

		    normalvector.Rotate ( - _dut_pre_align_a[i], _Xaxis );
		    normalvector.Rotate ( - _dut_pre_align_b[i], _Yaxis );
		    normalvector.Rotate ( - _dut_pre_align_c[i], _Zaxis );

		    // add shift -> subtract it
		    TVector3 alignshift ( _dut_pre_align_x[i], _dut_pre_align_y[i], _dut_pre_align_z[i] );
		    dutplane = dutplane - alignshift;
		}

		// the alignment gets added
		for ( unsigned int i = 0; i < _alignmentCollectionName.size ( ); i++ )
		{
		    streamlog_out ( DEBUG2 ) << "Performing alignment on DUT plane vector, step " << i + 1 << " !" << endl;

		    // add rotation -> negative rotation
		    dutplane.Rotate ( - _dut_align_a[i], _Xaxis );
		    dutplane.Rotate ( - _dut_align_b[i], _Yaxis );
		    dutplane.Rotate ( - _dut_align_c[i], _Zaxis );

		    normalvector.Rotate ( - _dut_align_a[i] ,_Xaxis );
		    normalvector.Rotate ( - _dut_align_b[i] ,_Yaxis );
		    normalvector.Rotate ( - _dut_align_c[i] ,_Zaxis );

		    // add shift -> subtract it
		    TVector3 alignshift ( _dut_align_x[i], _dut_align_y[i], _dut_align_z[i] );
		    dutplane = dutplane - alignshift;

		}

		TVector3 initialpos ( 0.0, 0.0, _z_refhit );
		dutplane = dutplane + initialpos;

		// finally the track length scalar
		double scalar = ( ( dutplane.X ( ) - _telescopeX[_iDUT - 2] ) * normalvector.X ( ) + ( dutplane.Y ( ) - _telescopeY[_iDUT - 2] ) * normalvector.Y ( ) + ( dutplane.Z ( ) - _planePosition[_iDUT - 2] ) * normalvector.Z ( ) ) / ( tr_x * normalvector.X ( ) + tr_y * normalvector.Y ( ) + tr_z * normalvector.Z ( ) );

		// And the final track position:
		planeintersectX = _telescopeX[_iDUT - 2] + scalar * tr_x;
		planeintersectY = _telescopeY[_iDUT - 2] + scalar * tr_y;
		planeintersectZ = _planePosition[_iDUT - 2] + scalar * tr_z;

		dutTrackX_global = planeintersectX;
		dutTrackY_global = planeintersectY;
		dutTrackZ_global = planeintersectZ;
	    }

	    // ****************** //
	    // 2: existing track fit!
	    // these variables will be set to _missingValue if the DUT was NOT in a previous track fit...
	    // preceeding daffitter processors should be run with settings to include the DUT in the trackfit
	    // ****************** //
	    if ( _useTrackFit == true && _useTriplet == false )
	    {
		dutTrackX_global = _fittedX[_iDUT];
		dutTrackY_global = _fittedY[_iDUT];
		dutTrackZ_global = _fittedZ[_iDUT];
		streamlog_out ( DEBUG6 ) << "Using DAF/GBL track fit!" << endl;
	    }

	    // ****************** //
	    // 3: triplet pointing: Assumes DUT is between 3 telescope planes!
	    // old: fit from planes 0 and 2
	    // new: fit from planes 1 and 2
	    // no real difference...
	    // ****************** //
	    //double old_triplet_X = ( _fittedX[_iDUT - 1] - _fittedX[_iDUT - 3] ) * _zDUT / _planePosition[_iDUT - 1]  + _fittedX[_iDUT - 3];
	    //double old_triplet_Y = ( _fittedY[_iDUT - 1] - _fittedY[_iDUT - 3] ) * _zDUT / _planePosition[_iDUT - 1]  + _fittedY[_iDUT - 3];
	    double triplet_X = ( _zDUT - _planePosition[_iDUT - 2] ) * ( _telescopeX[_iDUT - 1] - _telescopeX[_iDUT - 2] ) / ( _planePosition[_iDUT - 1] - _planePosition[_iDUT - 2] ) + _telescopeX[_iDUT - 2];
	    double triplet_Y = ( _zDUT - _planePosition[_iDUT - 2] ) * ( _telescopeY[_iDUT - 1] - _telescopeY[_iDUT - 2] ) / ( _planePosition[_iDUT - 1] - _planePosition[_iDUT - 2] ) + _telescopeY[_iDUT - 2];

	    // fill the difference between methods into a histogram
	    fillestimationhisto ( planeintersectX, planeintersectY, triplet_X, triplet_Y, _fittedX[_iDUT], _fittedY[_iDUT] );

	    if ( _useTrackFit == false && _useTriplet == true )
	    {
		dutTrackX_global = triplet_X;
		dutTrackY_global = triplet_Y;
		dutTrackZ_global = _planePosition[_iDUT];
		streamlog_out ( DEBUG6 ) << "Using triplet! " << endl;
	    }

	    streamlog_out ( DEBUG6 ) << "Intersect:         ( " << planeintersectX << " | " << planeintersectY << " | " << planeintersectZ << " )" << endl;
	    streamlog_out ( DEBUG6 ) << "DAF/GBL Track fit: ( " << _fittedX[_iDUT] << " | " << _fittedY[_iDUT] << " | " <<  _fittedZ[_iDUT] << " )" << endl;
	    streamlog_out ( DEBUG6 ) << "Triplet:           ( " << triplet_X << " | " << triplet_Y << " | " << _planePosition[_iDUT] << " )" << endl;
	    streamlog_out ( DEBUG6 ) << " " << endl;
	    streamlog_out ( DEBUG6 ) << "Track point:       ( " << dutTrackX_global << " | " << dutTrackY_global << " | " << dutTrackZ_global << " )" << endl;

	    // fill it into the ntuple
	    // this is the expected track impact position of the aligned and rotated DUT in the global coordinate system!
	    _FitTuple -> fill ( icol++, dutTrackX_global );
	    _FitTuple -> fill ( icol++, dutTrackY_global );
	    _FitTuple -> fill ( icol++, dutTrackZ_global );

	    // also dump this into a hitmap histogram
	    fillHitmapHisto ( dutTrackX_global, dutTrackY_global, 0 );

	    // we have found the global position of the track on the DUT
	    // we want the local coordinates (pixels) of the track for further analysis

	    // store the global position in a root 3vector
	    TVector3 iterationvector ( dutTrackX_global, dutTrackY_global, dutTrackZ_global );
	    TVector3 testvector ( 0.0, 0.0, 0.0 );

	    // as a check of the de-aligning process, we load an aligned hit, dialign it with the loaded collections, and compare it to the original hit it came from
	    int hittocheck = 0;
	    bool checkthisevent = false;
	    if ( _checkdealignment == true )
	    {

		LCCollectionVec * testCollection2;
		try
		{
		    testCollection2 = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( _inputDUTColName ) );
		}
		catch ( lcio::DataNotAvailableException& e )
		{
		    streamlog_out ( DEBUG5 ) << "No aligned hit with name " << testCollection2 << " ! Check _inputDUTColName! " << endl;
		    return;
		}
		int asize = testCollection2 -> getNumberOfElements ( );
		const double * pos2;
		for ( int j = 0; j < asize; j++ )
		{
		    TrackerHit * testhit2 = dynamic_cast < TrackerHit* > ( testCollection2 -> getElementAt ( j ) );
		    pos2 = testhit2 -> getPosition ( );
		    double distMin = 5.0;
		    int hitPlane = -1;
		    for ( int ipl = 0; ipl < _nTelPlanes; ipl++ )
		    {
			double dist =  pos2[2] - _zDUT;

			if ( dist * dist < distMin * distMin )
			{
			    hitPlane = 99;
			    distMin = dist;
			}
		    }
		    if ( hitPlane == 99 )
		    {
			streamlog_out ( DEBUG2 ) << "Test hit point:   ( " << pos2[0] << " | " << pos2[1] << " | " << pos2[2] << " )" << endl;
			// this is also the real point before dealignment :)
			streamlog_out ( DEBUG2 ) << "Real hit point:   ( " << pos2[0] << " | " << pos2[1] << " | " << pos2[2] << " )" << endl;
			TVector3 tempvector ( pos2[0], pos2[1], pos2[2] );
			testvector = testvector + tempvector;
			hittocheck = j;
			checkthisevent = true;
			break;
		    }
		}
	    }

	    // the z position is set to 0 due to the way the reference collection works...
	    if ( checkthisevent == true )
	    {
		iterationvector.SetZ ( testvector.Z ( ) );
	    }

	    // dealign the reference hit too
	    double refx = 0.0;
	    double refy = 0.0;
	    double refz = 0.0;

	    // the central coordinates of our DUT:
	    double dutcenter_x = 0.0;
	    double dutcenter_y = 0.0;
	    double dutcenter_z = 0.0;

	    // undo the alignment
	    for ( int i = _alignmentCollectionName.size ( ); i > 0; i-- )
	    {
		streamlog_out ( DEBUG2 ) << "Undoing alignment, step " << i << " !" << endl;

		// i-1 since the array starts at 0!

		// align the reference
		refx = _x_refhit + _dut_align_x[i - 1];
		refy = _y_refhit + _dut_align_y[i - 1];
		refz = _z_refhit + _dut_align_z[i - 1];

		TVector3 _RotatedVector( _a_refhit, _b_refhit, _c_refhit );
		_RotatedVector.RotateZ ( _dut_align_a[i - 1] );
		_RotatedVector.RotateY ( _dut_align_b[i - 1] );
		_RotatedVector.RotateX ( _dut_align_c[i - 1] );

		// align our global track (iterationvector)
		// relative to the reference
		double inputPosition[3] = { iterationvector.X ( ), iterationvector.Y ( ), iterationvector.Z ( ) };
		
		dutcenter_x = refx;
		dutcenter_y = refy;
		dutcenter_z = refz;

		dutcenter_x -= _dut_align_x[i - 1];
		dutcenter_y -= _dut_align_y[i - 1];
		dutcenter_z -= _dut_align_z[i - 1];

		inputPosition[0] = inputPosition[0] - dutcenter_x;
		inputPosition[1] = inputPosition[1] - dutcenter_y;
		inputPosition[2] = inputPosition[2] - dutcenter_z;

		double outputPosition[3] = { 0.0, 0.0, 0.0 };

		outputPosition[0] = dutcenter_x;
		outputPosition[1] = dutcenter_y;
		outputPosition[2] = dutcenter_z;

		TVector3 iCenterOfSensorFrame ( inputPosition[0], inputPosition[1], inputPosition[2] );

		iCenterOfSensorFrame.RotateZ ( _dut_align_c[i - 1] );
		iCenterOfSensorFrame.RotateY ( _dut_align_b[i - 1] );
		iCenterOfSensorFrame.RotateX ( _dut_align_a[i - 1] );

		outputPosition[0] += iCenterOfSensorFrame ( 0 );
		outputPosition[1] += iCenterOfSensorFrame ( 1 );
		outputPosition[2] += iCenterOfSensorFrame ( 2 );

		// second the shift
		outputPosition[0] += _dut_align_x[i - 1];
		outputPosition[1] += _dut_align_y[i - 1];
		outputPosition[2] += _dut_align_z[i - 1];

		iterationvector.SetX ( outputPosition[0] );
		iterationvector.SetY ( outputPosition[1] );
		iterationvector.SetZ ( outputPosition[2] );

		if ( _checkdealignment == true && checkthisevent == true )
		{

		    // align our global hit (testvector)
		    // relative to the reference
		    double inputPositionTest[3] = { testvector.X ( ), testvector.Y ( ), testvector.Z ( ) };

		    inputPositionTest[0] = inputPositionTest[0] - dutcenter_x;
		    inputPositionTest[1] = inputPositionTest[1] - dutcenter_y;
		    inputPositionTest[2] = inputPositionTest[2] - dutcenter_z;

		    double outputPositionTest[3] = { 0.0, 0.0, 0.0 };

		    outputPositionTest[0] = dutcenter_x;
		    outputPositionTest[1] = dutcenter_y;
		    outputPositionTest[2] = dutcenter_z;

		    TVector3 iCenterOfSensorFrameTest ( inputPositionTest[0], inputPositionTest[1], inputPositionTest[2] );

		    iCenterOfSensorFrameTest.RotateZ ( _dut_align_c[i - 1] );
		    iCenterOfSensorFrameTest.RotateY ( _dut_align_b[i - 1] );
		    iCenterOfSensorFrameTest.RotateX ( _dut_align_a[i - 1] );

		    outputPositionTest[0] += iCenterOfSensorFrameTest ( 0 );
		    outputPositionTest[1] += iCenterOfSensorFrameTest ( 1 );
		    outputPositionTest[2] += iCenterOfSensorFrameTest ( 2 );

		    // second the shift
		    outputPositionTest[0] += _dut_align_x[i - 1];
		    outputPositionTest[1] += _dut_align_y[i - 1];
		    outputPositionTest[2] += _dut_align_z[i - 1];

		    testvector.SetX ( outputPositionTest[0] );
		    testvector.SetY ( outputPositionTest[1] );
		    testvector.SetZ ( outputPositionTest[2] );

		    streamlog_out ( DEBUG2 ) << "De-aligned test hit is: ( " << testvector.X ( ) << " | " << testvector.Y ( ) << " | " << testvector.Z ( ) << " ) !" << endl;

		    LCCollectionVec * testCollection3;

		    std::string tempname = _alignedHitPrefix;
		    std::stringstream tempstream;
		    if ( i != 1 )
		    {
			tempstream << tempname << i - 1;
		    }
		    else if ( i == 1 )
		    {
			tempname = _prealignedHitPrefix;
			tempstream << tempname;
		    }

		    try
		    {
			testCollection3 = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( tempstream.str ( ) ) );
		    }
		    catch ( lcio::DataNotAvailableException& e )
		    {
			streamlog_out ( DEBUG5 ) << "Error in de-alignment iteration " << i << endl;
			return;
		    }
		    TrackerHit * testhit3 = dynamic_cast < TrackerHit* > ( testCollection3 -> getElementAt ( hittocheck ) );
		    const double * pos3 = testhit3 -> getPosition ( );
		    streamlog_out ( DEBUG2 ) << "De-aligned real hit is: ( " << pos3[0] << " | " << pos3[1] << " | " << pos3[2] << " ) !" << endl;
		}

		streamlog_out ( DEBUG2 ) << "De-aligned vector is:   ( " << iterationvector.X ( ) << " | " << iterationvector.Y ( ) << " | " << iterationvector.Z ( ) << " ) !" << endl;
	    }

	    // undo the prealignment
	    for ( int i = _pre_alignmentCollectionName.size ( ); i > 0; i-- )
	    {
		streamlog_out ( DEBUG2 ) << "Undoing prealignment, step " << i << " !" << endl;

		// i-1 since the array starts at 0!

		// align the reference
		refx = _x_refhit + _dut_pre_align_x[i - 1];
		refy = _y_refhit + _dut_pre_align_y[i - 1];
		refz = _z_refhit + _dut_pre_align_z[i - 1];

		TVector3 _RotatedVector ( _a_refhit, _b_refhit, _c_refhit );
		_RotatedVector.RotateZ ( _dut_pre_align_a[i - 1] );
		_RotatedVector.RotateY ( _dut_pre_align_b[i - 1] );
		_RotatedVector.RotateX ( _dut_pre_align_c[i - 1] );

		// align our global track (iterationvector)
		// relative to the reference
		double inputPosition[3] = { iterationvector.X ( ), iterationvector.Y ( ), iterationvector.Z ( ) };

		dutcenter_x = refx;
		dutcenter_y = refy;
		dutcenter_z = refz;

		dutcenter_x -= _dut_pre_align_x[i - 1];
		dutcenter_y -= _dut_pre_align_y[i - 1];
		dutcenter_z -= _dut_pre_align_z[i - 1];

		inputPosition[0] = inputPosition[0] - dutcenter_x;
		inputPosition[1] = inputPosition[1] - dutcenter_y;
		inputPosition[2] = inputPosition[2] - dutcenter_z;

		double outputPosition[3] = { 0.0, 0.0, 0.0 };

		outputPosition[0] = dutcenter_x;
		outputPosition[1] = dutcenter_y;
		outputPosition[2] = dutcenter_z;

		TVector3 iCenterOfSensorFrame ( inputPosition[0], inputPosition[1], inputPosition[2] );

		iCenterOfSensorFrame.RotateZ ( _dut_pre_align_c[i - 1] );
		iCenterOfSensorFrame.RotateY ( _dut_pre_align_b[i - 1] );
		iCenterOfSensorFrame.RotateX ( _dut_pre_align_a[i - 1] );

		outputPosition[0] += iCenterOfSensorFrame ( 0 );
		outputPosition[1] += iCenterOfSensorFrame ( 1 );
		outputPosition[2] += iCenterOfSensorFrame ( 2 );

		// second the shift
		outputPosition[0] += _dut_pre_align_x[i - 1];
		outputPosition[1] += _dut_pre_align_y[i - 1];
		outputPosition[2] += _dut_pre_align_z[i - 1];

		iterationvector.SetX ( outputPosition[0] );
		iterationvector.SetY ( outputPosition[1] );
		iterationvector.SetZ ( outputPosition[2] );

		if (_checkdealignment == true && checkthisevent == true )
		{

		    // align our global hit (testvector)
		    // relative to the reference
		    double inputPositionTest[3] = { testvector.X ( ), testvector.Y ( ), testvector.Z ( ) };

		    inputPositionTest[0] = inputPositionTest[0] - dutcenter_x;
		    inputPositionTest[1] = inputPositionTest[1] - dutcenter_y;
		    inputPositionTest[2] = inputPositionTest[2] - dutcenter_z;

		    double outputPositionTest[3] = { 0.0, 0.0, 0.0 };

		    outputPositionTest[0] = dutcenter_x;
		    outputPositionTest[1] = dutcenter_y;
		    outputPositionTest[2] = dutcenter_z;

		    TVector3 iCenterOfSensorFrameTest ( inputPositionTest[0], inputPositionTest[1], inputPositionTest[2] );

		    iCenterOfSensorFrameTest.RotateZ ( _dut_pre_align_c[i - 1] );
		    iCenterOfSensorFrameTest.RotateY ( _dut_pre_align_b[i - 1] );
		    iCenterOfSensorFrameTest.RotateX ( _dut_pre_align_a[i - 1] );

		    outputPositionTest[0] += iCenterOfSensorFrameTest ( 0 );
		    outputPositionTest[1] += iCenterOfSensorFrameTest ( 1 );
		    outputPositionTest[2] += iCenterOfSensorFrameTest ( 2 );

		    // second the shift
		    outputPositionTest[0] += _dut_pre_align_x[i - 1];
		    outputPositionTest[1] += _dut_pre_align_y[i - 1];
		    outputPositionTest[2] += _dut_pre_align_z[i - 1];

		    testvector.SetX ( outputPositionTest[0] );
		    testvector.SetY ( outputPositionTest[1] );
		    testvector.SetZ ( outputPositionTest[2] );

		    streamlog_out ( DEBUG2 ) << "De-prealigned test hit is: ( " << testvector.X ( ) << " | " << testvector.Y ( ) << " | " << testvector.Z ( ) << " ) !" << endl;

		    LCCollectionVec * testCollection3;

		    std::string tempname = _prealignedHitPrefix;
		    std::stringstream tempstream;
		    if ( i != 1 )
		    {
			tempstream << tempname << i - 1;
		    }
		    else if ( i == 1 )
		    {
			tempname = _unalignedcollectionname;
			tempstream << tempname;
		    }

		    try
		    {
			testCollection3 = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( tempstream.str ( ) ) );
		    }
		    catch ( lcio::DataNotAvailableException& e )
		    {
			streamlog_out ( DEBUG5 ) << "Error in de-prealignment iteration " << i << endl;
			return;
		    }
		    TrackerHit * testhit3 = dynamic_cast < TrackerHit* > ( testCollection3 -> getElementAt ( hittocheck ) );
		    const double * pos3 = testhit3 -> getPosition ( );
		    streamlog_out ( DEBUG2 ) << "De-prealigned real hit is: ( " << pos3[0] << " | " << pos3[1] << " | " << pos3[2] << " ) !" << endl;

		    // if we are at i = 1, then we have done the last de-alignment and de-prealignment step -> fill histos
		    // do a check on the z position to identify the sensor
		    // if there is no dut hit associated to the track, the first telescope hit will have been dealigned with the wrong alignment
		    if ( i == 1 && pos3[2] > 0.9 * _planePosition[_iDUT] && pos3[2] < 1.1 * _planePosition[_iDUT] )
		    {
			filldebughistos ( testvector.X ( ), pos3[0], testvector.Y ( ), pos3[1] );
		    }

		}

		streamlog_out ( DEBUG2 ) << "De-prealigned vector is:   ( " << iterationvector.X ( ) << " | " << iterationvector.Y ( ) << " | " << iterationvector.Z ( ) << " ) !" << endl;
	    }

	    // FIXME only x and y...
	    if (_useOriginalPreAlignment)
	    {
		iterationvector.SetX ( iterationvector.X ( ) + _dut_original_pre_align_x );
		iterationvector.SetY ( iterationvector.Y ( ) + _dut_original_pre_align_y );

		streamlog_out ( DEBUG2 ) << "Vector after undoing original prealignment is:   ( " << iterationvector.X ( ) << " | " << iterationvector.Y ( ) << " | " << iterationvector.Z ( ) << " ) !" << endl;
	    }

	    // iterationvector is now in local coordinates, still in mm
	    // revert gear rotations and shifts too
	    // FIXME for now only alpha and beta

	    float tempy = iterationvector.Y ( );
	    float roty = tempy / cos ( _DUTalign.at ( 3 ) * PI / 180.0 );
	    iterationvector.SetY ( roty );

	    float tempx = iterationvector.X ( );
	    float rotx = tempx / cos ( _DUTalign.at ( 4 ) * PI / 180.0 );
	    iterationvector.SetX ( rotx );

	    streamlog_out ( DEBUG6 ) << " " << endl;
	    streamlog_out ( DEBUG6 ) << "Local track hit:  ( " << iterationvector.X ( ) << " | " << iterationvector.Y ( ) << " | " << iterationvector.Z ( ) << " )" << endl;

	    // write the local information to the ntuple
	    _FitTuple -> fill ( icol++, iterationvector.X ( ) );
	    _FitTuple -> fill ( icol++, iterationvector.Y ( ) );
	    _FitTuple -> fill ( icol++, iterationvector.Z ( ) );

	    double dutTrackX_local_pix = ( _pixelx / 2.0 ) - 0.5 - iterationvector.X ( ) / _pitchx; 
	    double dutTrackY_local_pix = ( _pixely / 2.0 ) - 0.5 - iterationvector.Y ( ) / _pitchy;

	    streamlog_out ( DEBUG6 ) << "Points to pixels: ( " << dutTrackX_local_pix << " | " << dutTrackY_local_pix << " )" << endl;
	    streamlog_out ( DEBUG6 ) << " " << endl;

	    // write the pixel information to the ntuple
	    _FitTuple -> fill ( icol++, dutTrackX_local_pix );
	    _FitTuple -> fill ( icol++, dutTrackY_local_pix );
	    // local z not defined
	    _FitTuple -> fill ( icol++, _missingValue );

	    // for writing to the tuple: the local position of the matched hit in the collection
	    double foundhitlocalX = _missingValue;
	    double foundhitlocalY = _missingValue;
	    double foundhitlocalZ = _missingValue;

	    // for writing to the tuple: the pixels of the matched hit
	    double foundhitlocalPixX = _missingValue;
	    double foundhitlocalPixY = _missingValue;
	    double foundhitlocalPixZ = _missingValue;

	    // now we have the expected track impact position on the DUT in global and local coordinates.
	    // here we search for an actual DUT hit in the aligned track hits, this is done by looking at the z coordinate...
	    // _DUTok is a flag to see if the DUT hit collection has been loaded
	    if ( _DUTok )
	    {
		streamlog_out ( DEBUG2 ) << "Checking all hits if they are from the DUT..." << endl;
		double distPrevX = 99999.0;
		double distPrevY = 99999.0;
		const double * pos;
		double distX = 0.0;
		double distY = 0.0;
		int matchhit = -1;;
		double matchpos[3] = { 0.0 };

		// loop over the input hits 
		for ( int ihit = 0; ihit < nInputHits; ihit++ )
		{
		    TrackerHitImpl * meshit = dynamic_cast < TrackerHitImpl* > ( hitcol -> getElementAt ( ihit ) );

		    // hit position, is it from the dut?
		    pos = meshit -> getPosition ( );

		    UTIL::CellIDDecoder < TrackerHitImpl > hitDecoder ( EUTELESCOPE::HITENCODING );
		    const int tempid = hitDecoder ( meshit ) ["sensorID"];

		    if ( tempid == _manualDUTid )
		    {

			streamlog_out ( DEBUG3 ) << "Found possible DUT hit..." << endl;

			int clustersizeX = -1;
			int clustersizeY = -1;
			float clustercharge = 0.0;

			TrackerDataImpl* clusterVector = static_cast < TrackerDataImpl* > ( meshit -> getRawHits ( ) [0] );
			EUTelSimpleVirtualCluster * cluster = nullptr;
			cluster = new EUTelSparseClusterImpl < EUTelGenericSparsePixel > ( clusterVector );
			if ( cluster != nullptr )
			{
			    cluster -> getClusterSize ( clustersizeX, clustersizeY );
			    clustercharge = cluster -> getTotalCharge ( );
			}

			// for now we output via these vars to the ntuple...
			dutHitQ = clustercharge;
			dutHitR = clustersizeY;

			streamlog_out ( DEBUG3 ) << "Cluster size in x: " << clustersizeX << " size in y: " << clustersizeY << " charge: " << clustercharge << endl;

			// we compare against the fit...
			streamlog_out ( DEBUG4 ) << " " << endl;

			// we also fill a histogram of all hits on the dut for reference
			fillLocalHitmapHisto ( pos[0], pos[1], pos[2], 0 );

			streamlog_out ( DEBUG7 ) << "Comparing DUT hit  ( " << pos[0] << " | " << pos[1] << " | " << pos[2] << " )" << endl;
			streamlog_out ( DEBUG7 ) << "with DUT track hit ( " << dutTrackX_global << " | " << dutTrackY_global << " | " << dutTrackZ_global << ") " << endl;

			// the distance is between the hit and the fit!
			distX = ( pos[0] - dutTrackX_global );
			distY = ( pos[1] - dutTrackY_global );
			//dutHitR = sqrt(fabs(distX)*fabs(distY) );

			// are we within the maximum range to include the hit ?
			// also, if there are multiple hits in an event, only consider the hit nearest to this track
			// multiple tracks in an event can have different (or the same) hits
			if ( sqrt ( distX * distX ) < fabs ( _distMax_X ) && sqrt ( distY * distY ) < fabs ( _distMax_Y ) )
			{
			    if ( sqrt ( distX * distX ) < fabs ( distPrevX ) && sqrt ( distY * distY ) < fabs ( distPrevY ) )
			    {

				distPrevX = sqrt ( distX * distX );
				distPrevY = sqrt ( distY * distY );

				streamlog_out ( DEBUG5 ) << "###############################################################################################################################################" << endl;
				streamlog_out ( DEBUG6 ) << "Matched DUT hit : ( " << pos[0] << " | " << pos[1] << " )" << endl;
				streamlog_out ( DEBUG5 ) << "Distance in global mm is: X: " << distX << " , Y: " << distY << " -> " << dutHitR << " mm!" << endl;
				streamlog_out ( DEBUG5 ) << "###############################################################################################################################################" << endl;

				matchhit = ihit;
				matchpos[0] = pos[0];
				matchpos[1] = pos[1];
				matchpos[2] = pos[2];
			    }
			    else
			    {
				streamlog_out ( DEBUG5 ) << "Found closer DUT hit in event " << _evtNr << endl;
				streamlog_out ( DEBUG3 ) << " dX current: " << sqrt ( distX * distX ) << ", previous dX: " << fabs ( distPrevX ) <<", dY current: " << sqrt ( distY * distY ) << " previous dY: " << fabs ( distPrevY ) << endl;
			    }

			}
			else
			{
			    streamlog_out ( DEBUG4 ) << "DUT hit not matched!" << endl;
			    _foundDUTHit = false;
			}

		    } // z pos loop
	
		} // input hit loop

		if ( matchhit > -1 )
		{
		    // fill residual
		    fillresihistos ( matchpos[0], matchpos[1], dutTrackX_global, dutTrackY_global, _evtNr, trackchi2, trackndf );

		    // if we pass fiducial cuts fill an extra histogram:
		    if ( dutTrackX_global > _fiducialcut[0] && dutTrackX_global < _fiducialcut[1] && dutTrackY_global > _fiducialcut[2] && dutTrackY_global < _fiducialcut[3] )
		    {
			streamlog_out ( DEBUG2 ) << "Passed fiducial cut!" << endl;
			Filteredfillresihistos ( matchpos[0], matchpos[1], dutTrackX_global, dutTrackY_global );
		    }

		    _foundDUTHit = true;
		    _matchedhits++;

		    // if the DUT was not in the track fit, then _measuredX and _measuredY don't get filled in the first track loop, so we do this here
		    if ( _measuredX[_iDUT] == _missingValue || _measuredY[_iDUT] == _missingValue )
		    {
			_measuredX[_iDUT] = matchpos[0];
			_measuredY[_iDUT] = matchpos[1];
			_measuredZ[_iDUT] = matchpos[2];
		    }

		    LCCollection* foundhitcol = event -> getCollection ( _originalCollectionName );
		    TrackerHit * foundhit = dynamic_cast < TrackerHit* > ( foundhitcol -> getElementAt ( matchhit ) );

		    const double * foundpos = foundhit -> getPosition ( );
		    foundhitlocalX = foundpos[0];
		    foundhitlocalY = foundpos[1];
		    foundhitlocalZ = foundpos[2];

		    // subtract the original prealignment
		    // FIXME for now only x and y shifts
		    if ( _useOriginalPreAlignment )
		    {
			foundhitlocalX = foundhitlocalX + _dut_original_pre_align_x;
			foundhitlocalY = foundhitlocalY + _dut_original_pre_align_y;
			foundhitlocalZ = foundhitlocalZ + _dut_original_pre_align_z;
		    }

		    foundhitlocalPixX = ( _pixelx / 2.0 ) - 0.5 - foundhitlocalX / _pitchx; 
		    foundhitlocalPixY = ( _pixely / 2.0 ) - 0.5 - foundhitlocalY / _pitchy;
		    foundhitlocalPixZ = _missingValue;

		    streamlog_out ( DEBUG4 ) << "Found hit points to pixels: ( " << foundhitlocalPixX << " | " << foundhitlocalPixY << " ) !" << endl;

		    // fill a matched hitmap
		    fillLocalHitmapHisto ( matchpos[0], matchpos[1], matchpos[2], 1 );
		    fillHitmapHisto ( dutTrackX_global, dutTrackY_global, 1 );
		    fillInterStripHitmap ( foundhitlocalPixX, foundhitlocalPixY );

		}

		// some alibava specific output:
		// we want to see which channel (aka cluster seed) created the matched hit... this can be used as a check on the reconstructed track impact points...
		if ( _alibava == true )
		{

		    float locdistx = 0.0;
		    float locdisty = 0.0;

		    LCCollectionVec * tempCollectionVec = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( _alibavaClusterCollectionName ) );
		    int tempsize = tempCollectionVec -> getNumberOfElements ( );
		    CellIDDecoder < TrackerPulseImpl > inputDecoder ( tempCollectionVec );
		    for ( int i = 0; i < tempsize; i++ )
		    {
			lcio::TrackerPulseImpl * input  = dynamic_cast < lcio::TrackerPulseImpl * > ( tempCollectionVec -> getElementAt ( i ) );
			int xSeed = inputDecoder ( input ) ["xSeed"];
			int ySeed = inputDecoder ( input ) ["ySeed"];
			streamlog_out ( DEBUG5 ) << "Check strip / pixel hit of ( " << dutTrackX_local_pix << " | " << dutTrackY_local_pix << " ) against ( " << xSeed << " | " << ySeed << " ) !" << endl;
			locdistx = ( xSeed - dutTrackX_local_pix ) * _pitchx;
			locdisty = ( ySeed - dutTrackY_local_pix ) * _pitchy;

			streamlog_out ( DEBUG5 ) << "Distance in local mm is: X: " << locdistx << " , Y: " << locdisty << " -> " << sqrt ( fabs ( locdistx ) * fabs ( locdisty ) ) << " mm!" << endl;
			streamlog_out ( DEBUG5 ) << "###############################################################################################################################################" << endl;
			streamlog_out ( DEBUG5 ) << " " << endl;

			// there might be more clusters in an event, but we are too lazy to identify the one which made this hit...
			if ( tempsize < 2 )
			{
			    fillPrecisionHisto ( dutTrackX_local_pix - xSeed, dutTrackY_local_pix - ySeed );
			}
		    }

		    LCCollectionVec * alibavaCollection;
		    try
		    {
			alibavaCollection = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( _alibavaCollectionName ) );

			int noOfDetector = alibavaCollection -> getNumberOfElements ( );
			for ( int i = 0; i < noOfDetector; ++i )
			{

			    TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( alibavaCollection -> getElementAt ( i ) );
			    FloatVec datavec;
			    datavec = trkdata -> getChargeValues ( );
			    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( event );
			    float tdctime = alibavaEvent -> getEventTime ( );

			    // track eta plot: fill a histogram with the charge of two strips left and right of the track pointer
			    // this charge is also used for filling duthitq...
			    // FIXME: only one dut...
			    if ( _foundDUTHit == true )
			    {
				if ( _unsensitiveaxis == "x" && i == 0 )
				{
				    double tempvar = dutTrackY_local_pix;
				    int tempvar2 = static_cast < float > ( tempvar );

				    if ( tempvar2 > 1 && tempvar2 < 125 )
				    {
					fillEtaHisto ( datavec[tempvar2], datavec[tempvar2 + 1], datavec[tempvar2 - 1], datavec[tempvar2 + 2], datavec[tempvar2 - 2], datavec[tempvar2 + 3], dutTrackY_local_pix, _evtNr, tdctime );

				    }

				}
				else if ( _unsensitiveaxis == "y"  && i == 0 )
				{

				    double tempvar = dutTrackX_local_pix;
				    int tempvar2 = static_cast < float > ( tempvar );

				    if ( tempvar2 > 1 && tempvar2 < 125 )
				    {
					fillEtaHisto ( datavec[tempvar2], datavec[tempvar2 + 1], datavec[tempvar2 - 1], datavec[tempvar2 + 2], datavec[tempvar2 - 2], datavec[tempvar2 + 3], dutTrackX_local_pix, _evtNr, tdctime );

				    }

				}

				fillTDCResHisto ( _measuredX[_iDUT], _measuredY[_iDUT], dutTrackX_global, dutTrackY_global, tdctime );
			    }

			} // end noOfDetector loop

		    }
		    catch ( lcio::DataNotAvailableException& e )
		    {
			streamlog_out ( DEBUG5 ) << "No alibava collection with name " << _alibavaCollectionName << " ! Event " << event -> getEventNumber ( ) << " might be masked!" << endl;
			continue;
		    }

		} // end of if alibava

	    } // end of if(_DUTok)

	    // these are filled regardless of if there is a DUT or not to preserve structure
	    _FitTuple -> fill ( icol++, _measuredX[_iDUT] );
	    _FitTuple -> fill ( icol++, _measuredY[_iDUT] );
	    _FitTuple -> fill ( icol++, _measuredZ[_iDUT] );
	    _FitTuple -> fill ( icol++, foundhitlocalX );
	    _FitTuple -> fill ( icol++, foundhitlocalY );
	    _FitTuple -> fill ( icol++, foundhitlocalZ );
	    _FitTuple -> fill ( icol++, foundhitlocalPixX );
	    _FitTuple -> fill ( icol++, foundhitlocalPixY );
	    _FitTuple -> fill ( icol++, foundhitlocalPixZ );
	    _FitTuple -> fill ( icol++, dutHitR );
	    _FitTuple -> fill ( icol++, dutHitQ );

	    // if we have an Alibava DUT, we want the header info and the reco data in the ntuple
	    // the structure of the ntuple will have been expanded, so writing should be no problem
	    if ( _alibava == true )
	    {
		LCCollectionVec * alibavaCollection;
		try
		{
		    alibavaCollection = dynamic_cast < LCCollectionVec * > ( event -> getCollection ( _alibavaCollectionName ) );

		    // header: tdc and temperature
		    AlibavaEventImpl * alibavaEvent = static_cast < AlibavaEventImpl* > ( event );
		    float tdctime = alibavaEvent -> getEventTime ( );
		    float temperature = alibavaEvent -> getEventTemp ( );
		    _FitTuple -> fill ( icol++, tdctime );
		    _FitTuple -> fill ( icol++, temperature );

		    // loop over detectors
		    int noOfDetector = alibavaCollection -> getNumberOfElements ( );
		    for ( int i = 0; i < noOfDetector; ++i )
		    {
			TrackerDataImpl * trkdata = dynamic_cast < TrackerDataImpl * > ( alibavaCollection -> getElementAt ( i ) );
			FloatVec datavec;
			datavec = trkdata -> getChargeValues ( );

			// channel loop
			for ( size_t ichan = 0; ichan < datavec.size ( ); ichan++ )
			{
			    double alidata = 0.0;
			    alidata = datavec[ichan];
			    _FitTuple -> fill ( icol++, alidata );
			    streamlog_out ( DEBUG0 ) << "Filling alibava channel " << ichan << " with reconstructed data: " << alidata << " ADCs" << endl;
			}

		    }
		}
		catch ( lcio::DataNotAvailableException& e )
		{
		    streamlog_out ( DEBUG5 ) << "No alibava collection with name " << _alibavaCollectionName << " ! Event " << event -> getEventNumber ( ) << " might be masked!" << endl;
		    continue;
		}

	    } // end of if Alibava

	} // end of if DUT is present (_iDUT != -1)

	// now we should have all relevant information in the tuple
	_FitTuple->addRow();

    } // end of loop over tracks

    if ( event -> getEventNumber ( ) % 1000 == 0 )
    {
	streamlog_out ( MESSAGE5 ) << "Matched " << _matchedhits << " DUT hits at event " << event -> getEventNumber ( ) << endl;
    }

    return;
} // end of event loop

void EUTelFitTupleAlibava::check ( LCEvent * /* evt */ )
{

}

// finally, additional aligment is calculated if so wanted
// the resoltion histograms are fitted and some output is displayed
void EUTelFitTupleAlibava::end ( )
{

    // do lots of stuff for DUTs
    if ( _iDUT != -1 )
    {

	// now that we've filled the residual histograms, we can do a fit to obtain the resolutions...
	TH1D * histo = dynamic_cast < TH1D* > ( _rootObjectMap["ResidualX"] );
	TF1 * tempfit = dynamic_cast < TF1* > ( _rootObjectMap["ResidualXFit"] );
	histo -> Fit ( tempfit, "Q" );

	TH1D * histo2 = dynamic_cast < TH1D* > ( _rootObjectMap["ResidualY"] );
	TF1 * tempfit2 = dynamic_cast < TF1* > ( _rootObjectMap["ResidualYFit"] );
	histo2 -> Fit ( tempfit2, "Q" );

	// make an integral of the track based eta distribution
	if ( _alibava == true )
	{
	    double counts = 0.0;
	    double integral = 0.0;
	    TH1D * etahisto = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaCentral"] );
	    TH1D * etaintegral = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaIntegral"] );
	    for ( int i = 1; i < etahisto -> GetNbinsX ( ); i++ )
	    {
		counts = etahisto -> GetBinContent ( i );
		integral += counts;
		etaintegral -> SetBinContent ( i, integral );
	    }
	}

	// we also fit the profile plots to get some last-minute alignment...

	// FIXME until now only gamma (XY) rotations, as this is the easiest to align with a strip sensor DUT
	// for gamma, the XY rotation there can only be one significant residual, dX or dY
	double additionalZYrot = 0.0; // a
	double additionalZXrot = 0.0; // b
	double additionalXYrot = 0.0; // g

	if ( _unsensitiveaxis == "x" )
	{
	    TProfile* ResidualProfileDY_X = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDY_X"] );
	    TF1 * profilefitX = dynamic_cast < TF1* > ( _rootObjectMap["ProfileXFit"] );

	    // find the limits of the plot, so that we don't include the edges of the sensor...
	    float x_lobin = 0;
	    float x_hibin = 150;
	    for ( int i = 0; i < 150; i++ )
	    {
		if ( fabs ( ResidualProfileDY_X -> GetBinContent ( i ) ) > 0 )
		{
		    x_hibin = i;
		}
	    }
	    for ( int i = 149; i >= 0; i-- )
	    {
		if ( fabs ( ResidualProfileDY_X -> GetBinContent ( i ) ) > 0 && fabs ( ResidualProfileDY_X -> GetBinContent ( i ) ) < 0.1 )
		{
		    x_lobin = i;
		}
	    }
	    // -15.0 for the histo range, / 5.0 for the bin per um and 4/6 to disregard the edges
	    float x_min = -15.0 + x_lobin / 5.0 + 6.0;
	    float x_max = -15.0 + x_hibin / 5.0 - 4.0;
	    profilefitX -> SetRange ( x_min, x_max );
	    streamlog_out ( DEBUG1 ) << "Lower bin in x: " << x_lobin << " , upper bin in x: " << x_hibin << " !" << endl;
	    ResidualProfileDY_X -> Fit ( profilefitX, "QR" );

	    // get an additional XY rotation for alignment:
	    additionalXYrot = atan ( profilefitX -> GetParameter ( 1 ) );

	    TProfile* ResidualProfileDY_Y = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDY_Y"] );
	    TF1 * profilefitY = dynamic_cast < TF1* > ( _rootObjectMap["ProfileYFit"] );

	    // find the limits of the plot, so that we don't include the edges of the sensor...
	    x_lobin = 0;
	    x_hibin = 150;
	    for ( int i = 0; i < 150; i++ )
	    {
		if ( fabs ( ResidualProfileDY_Y -> GetBinContent ( i ) ) > 0 )
		{
		    x_hibin = i;
		}
	    }
	    for ( int i = 149; i >= 0; i-- )
	    {
		if ( fabs ( ResidualProfileDY_Y -> GetBinContent ( i ) ) > 0 && fabs ( ResidualProfileDY_Y -> GetBinContent ( i ) ) < 0.1 )
		{
		    x_lobin = i;
		}
	    }
	    // -15.0 for the histo range, / 5.0 for the bin per um and 0.5 to disregard the edges
	    x_min = -15.0 + x_lobin / 5.0 + 0.5;
	    x_max = -15.0 + x_hibin / 5.0 - 0.5;
	    profilefitY -> SetRange ( x_min, x_max );
	    streamlog_out ( DEBUG1 ) << "Lower bin in x: " << x_lobin << " , upper bin in x: " << x_hibin << " !" << endl;
	    streamlog_out ( DEBUG9 ) << "Lower pos in x: " << x_min << " , upper pos in x: " << x_max << " !" << endl;
	    ResidualProfileDY_Y -> Fit ( profilefitY, "QR" );

	    // get an additional ZY rotation for alignment:
	    additionalZYrot = atan ( profilefitY -> GetParameter ( 1 ) );

	}
	else if ( _unsensitiveaxis == "y" )
	{

	    // find the limits of the plot, so that we don't include the edges of the sensor...
	    TProfile* ResidualProfileDX_Y = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDX_Y"] );
	    TF1 * profilefitY = dynamic_cast < TF1* > ( _rootObjectMap["ProfileYFit"] );
	    float y_lobin = 0;
	    float y_hibin = 150;
	    for ( int i = 0; i < 150; i++ )
	    {
		if ( fabs ( ResidualProfileDX_Y -> GetBinContent ( i ) ) > 0 )
		{
		    y_hibin = i;
		}
	    }
	    for ( int i = 149; i >= 0; i-- )
	    {
		if ( fabs ( ResidualProfileDX_Y -> GetBinContent ( i ) ) > 0 )
		{
		    y_lobin = i;
		}
	    }
	    float y_min = -15.0 + y_lobin / 5.0 + 4.0;
	    float y_max = -15.0 + y_hibin / 5.0 - 4.0;
	    profilefitY -> SetRange ( y_min, y_max );
	    streamlog_out ( DEBUG1 ) << "Lower bin in y: " << y_lobin << " , upper bin in y: " << y_hibin << " !" << endl;
	    ResidualProfileDX_Y -> Fit ( profilefitY, "QR" );

	    // get an additional XY rotation for alignment:
	    additionalXYrot = atan ( profilefitY -> GetParameter ( 1 ) );

	    // find the limits of the plot, so that we don't include the edges of the sensor...
	    TProfile* ResidualProfileDX_X = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDX_X"] );
	    TF1 * profilefitX = dynamic_cast < TF1* > ( _rootObjectMap["ProfileXFit"] );
	    y_lobin = 0;
	    y_hibin = 150;
	    for ( int i = 0; i < 150; i++ )
	    {
		if ( fabs ( ResidualProfileDX_X -> GetBinContent ( i ) ) > 0 )
		{
		    y_hibin = i;
		}
	    }
	    for ( int i = 149; i >= 0; i-- )
	    {
		if ( fabs ( ResidualProfileDX_X -> GetBinContent ( i ) ) > 0 )
		{
		    y_lobin = i;
		}
	    }
	    y_min = -15.0 + y_lobin / 5.0 + 4.0;
	    y_max = -15.0 + y_hibin / 5.0 - 4.0;
	    profilefitX -> SetRange ( y_min, y_max );
	    streamlog_out ( DEBUG1 ) << "Lower bin in y: " << y_lobin << " , upper bin in y: " << y_hibin << " !" << endl;
	    ResidualProfileDX_X -> Fit ( profilefitX, "QR" );

	    // get an additional ZX rotation for alignment:
	    additionalZXrot = atan ( profilefitX -> GetParameter ( 1 ) );

	}

	// only output an additional rotation if we have enough statistics for this, otherwise just output zero
	if ( _matchedhits < 1000 )
	{
	    additionalXYrot = 0.0;
	    additionalZYrot = 0.0;
	    streamlog_out ( WARNING5 ) << "Not enough matched DUT hits ( " << _matchedhits << " < 1000 ) for additional alignment! Defaulting to zero!" << endl;
	}

	// some output of the plot fits:
	streamlog_out ( MESSAGE5 ) << " " << endl;
	streamlog_out ( MESSAGE5 ) << "===============================================================================================================================================" << endl;
	streamlog_out ( MESSAGE5 ) << "The profile plot suggests an additional ZY rotation of: " << additionalZYrot << " rad!" << endl;
	streamlog_out ( MESSAGE5 ) << "The profile plot suggests an additional ZX rotation of: " << additionalZXrot << " rad!" << endl;
	streamlog_out ( MESSAGE5 ) << "The profile plot suggests an additional XY rotation of: " << additionalXYrot << " rad!" << endl;
	streamlog_out ( MESSAGE5 ) << " " << endl;

	if ( _doAlignment == true )
	{
	    streamlog_out ( MESSAGE5 ) << "Will now write an additional rotation constant to file!" << endl;
	    streamlog_out ( MESSAGE5 ) << " " << endl;
	}

	streamlog_out ( MESSAGE5 ) << "===============================================================================================================================================" << endl;

	// with the residual plots we can now perform one last rotational alignment of the DUT
	// the alignment collection is written here
	if ( _doAlignment == true )
	{

	    // writer
	    LCWriter * lcWriter = LCFactory::getInstance ( ) -> createLCWriter ( );
	    try 
	    {
		    lcWriter -> open ( _outputalignment, LCIO::WRITE_NEW );
	    }
	    catch ( IOException& e )
	    {
		streamlog_out ( ERROR5 ) << e.what ( ) << endl;
		exit ( -1 );
	    }

	    // write an almost empty run header
	    LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
	    lcHeader -> setRunNumber ( 0 );
	    lcWriter -> writeRunHeader ( lcHeader );
	    delete lcHeader;

	    // an event:
	    LCEventImpl * event = new LCEventImpl;
	    event -> setRunNumber ( 0 );
	    event -> setEventNumber ( 0 );

	    // the alignment constant collection we want to write
	    // FIXME this assumes the DUT is in the middle of the telescope, with telescope planes 012 before it and 345 after it

	    LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );

	    // upstream:
	    for ( int i = 0; i < 3; i++ )
	    {
		EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
		constant -> setSensorID ( i );
		constant -> setXOffset ( 0.0 );
		constant -> setYOffset ( 0.0 );
		constant -> setZOffset ( 0.0 );
		constant -> setAlpha ( 0.0 );
		constant -> setBeta ( 0.0 );
		constant -> setGamma ( 0.0 );
		constantsCollection -> push_back ( constant );
	    }

	    // only g or (a or b), depending on unsensitive axis
	    if ( _doGamma )
	    {
		additionalZYrot = 0.0;
		additionalZXrot = 0.0;
	    }
	    else
	    {
		additionalXYrot = 0.0;
	    }

	    // DUT:
	    EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
	    constant -> setSensorID ( _manualDUTid );
	    constant -> setXOffset ( 0.0 );
	    constant -> setYOffset ( 0.0 );
	    constant -> setZOffset ( 0.0 );
	    constant -> setAlpha ( additionalZYrot );
	    constant -> setBeta ( additionalZXrot );
	    constant -> setGamma ( additionalXYrot );
	    constantsCollection -> push_back ( constant );

	    // downstream:
	    for ( int i = 3; i < _nTelPlanes - 1; i++ )
	    {
		EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;
		constant -> setSensorID ( i );
		constant -> setXOffset ( 0.0 );
		constant -> setYOffset ( 0.0 );
		constant -> setZOffset ( 0.0 );
		constant -> setAlpha ( 0.0 );
		constant -> setBeta ( 0.0 );
		constant -> setGamma ( 0.0 );
		constantsCollection -> push_back ( constant );
	    }

	    // collection name alignment is hardcoded, as it is in EuTelMille
	    // in future versions of EUTelescope, this might change
	    event -> addCollection ( constantsCollection, "alignment" );

	    // output all this
	    lcWriter -> writeEvent ( event );
	    delete event;

	    lcWriter -> close ( );

	} // end of _doAlignment

	dolandaugausfit ( "TrackSignal" );

	// some final output
	streamlog_out ( MESSAGE5 ) << " " << endl;
	streamlog_out ( MESSAGE5 ) << "N-tuple with " << _FitTuple -> rows ( ) << " rows created!" << endl;
	streamlog_out ( MESSAGE5 ) << " " << endl;
	streamlog_out ( MESSAGE5 ) << "===============================================================================================================================================" << endl;
	streamlog_out ( MESSAGE5 ) << " " << endl;
	streamlog_out ( MESSAGE5 ) << "DUT Resolution:               X:     " << tempfit -> GetParameter ( 2 ) << " +/- " << tempfit -> GetParError ( 2 ) << " mm" << endl;
	streamlog_out ( MESSAGE5 ) << "                              Y:     " << tempfit2 -> GetParameter ( 2 ) << " +/- " << tempfit2 -> GetParError ( 2 ) << " mm" << endl;
	streamlog_out ( MESSAGE5 ) << " " << endl;
	streamlog_out ( MESSAGE5 ) << "From " << _matchedhits << " matched hits!" << endl;
	streamlog_out ( MESSAGE5 ) << "===============================================================================================================================================" << endl;

	} // end of if DUT

	// clean memory
	delete [] _planeSort;
	delete [] _planePosition;
	delete [] _planeID;
	delete [] _isActive;
	delete [] _isMeasured;
	delete [] _isFitted;
	delete [] _measuredX;
	delete [] _measuredY;
	delete [] _measuredZ;
	delete [] _measuredQ;
	delete [] _fittedX;
	delete [] _fittedY;
	delete [] _telescopeX;
	delete [] _telescopeY;

}

// begin histogram filling functions
// a dut hitmap plot, 0 for all hits on the DUT, 1 for those that are matched
void EUTelFitTupleAlibava::fillLocalHitmapHisto ( float x, float y, float z, int map )
{
    if ( map == 0 )
    {
	if ( TH3D * hitmapallhits = dynamic_cast < TH3D* > ( _rootObjectMap["DUTHitmap"] ) )
	{
	    hitmapallhits -> Fill ( x, z - _planePosition[_iDUT], y );
	}
	if ( TH2D * hitmapallhitsmodpitch = dynamic_cast < TH2D* > ( _rootObjectMap["DUTHitmapModPitch"] ) )
	{
	    double fractpartx, fractparty, intpartx, intparty;
	    fractpartx = modf ( x / _pitchx , &intpartx );
	    fractparty = modf ( y / _pitchy , &intparty );
	    hitmapallhitsmodpitch -> Fill ( fabs ( fractpartx ), fabs ( fractparty ) );
	}
    }
    else if ( map == 1)
    {
	if ( TH3D * hitmapmatchedhits = dynamic_cast < TH3D* > ( _rootObjectMap["DUTHitmapMatched"] ) )
	{
	    hitmapmatchedhits -> Fill ( x, z - _planePosition[_iDUT], y );
	}
	if ( TH2D * hitmapallhitsmatchedmodpitch = dynamic_cast < TH2D* > ( _rootObjectMap["DUTHitmapMatchedModPitch"] ) )
	{
	    double fractpartx, fractparty, intpartx, intparty;
	    fractpartx = modf ( x / _pitchx , &intpartx );
	    fractparty = modf ( y / _pitchy , &intparty );
	    hitmapallhitsmatchedmodpitch -> Fill ( fabs ( fractpartx ), fabs ( fractparty ) );
	}
    }
}

// a track hitmap plot
// it is filled for dut hits which are matched to a track (1) or not (0)
// position mod pitch is also plotted to see if the fitting/matching somehow biases the track
void EUTelFitTupleAlibava::fillHitmapHisto ( float hitx, float hity, int match )
{
    if ( match == 0 )
    {
	if ( TH2D * hitmapHisto = dynamic_cast < TH2D* > ( _rootObjectMap["TrackHitmap"] ) )
	{
	    hitmapHisto -> Fill ( hitx, hity );
	}
	if ( TH2D * hitmapHistoModPitch = dynamic_cast < TH2D* > ( _rootObjectMap["TrackHitmapModPitch"] ) )
	{
	    double fractpartx, fractparty, intpartx, intparty;
	    fractpartx = modf ( hitx / _pitchx , &intpartx );
	    fractparty = modf ( hity / _pitchy , &intparty );
	    hitmapHistoModPitch -> Fill ( fabs ( fractpartx ), fabs ( fractparty ) );
	}
    }
    else if ( match == 1 )
    {
	if ( TH2D * hitmapHistoMatch = dynamic_cast < TH2D* > ( _rootObjectMap["TrackHitmapMatch"] ) )
	{
	    hitmapHistoMatch -> Fill ( hitx, hity );
	}
	if ( TH2D * hitmapHistoMatchModPitch = dynamic_cast < TH2D* > ( _rootObjectMap["TrackHitmapMatchModPitch"] ) )
	{
	    double fractpartx, fractparty, intpartx, intparty;
	    fractpartx = modf ( hitx / _pitchx , &intpartx );
	    fractparty = modf ( hity / _pitchy , &intparty );
	    hitmapHistoMatchModPitch -> Fill ( fabs ( fractpartx ), fabs ( fractparty ) );
	}
    }
}

// a plot to check reconstruction of local coordinates
// this gets called with the distance in local pitch units between track impact prediction and (Alibava) seed cluster position
void EUTelFitTupleAlibava::fillPrecisionHisto ( float x, float y )
{
    if ( TH2D * precHisto = dynamic_cast < TH2D* > ( _rootObjectMap["Precision"] ) )
    {
	precHisto -> Fill ( x, y );
    }
}

// some plots to compare the different track impact point estimation methods...
void EUTelFitTupleAlibava::fillestimationhisto ( float x1, float y1, float x2, float y2, float x3, float y3 )
{
    if ( TH2D * estimateHisto1 = dynamic_cast < TH2D* > ( _rootObjectMap["EstimationMethodDifferenceIntersectTriplet"] ) )
    {
	estimateHisto1 -> Fill ( x1 - x2, y1 - y2 );
    }
    if ( TH2D * estimateHisto2 = dynamic_cast < TH2D* > ( _rootObjectMap["EstimationMethodDifferenceIntersectDAF"] ) )
    {
	estimateHisto2 -> Fill ( x1 - x3, y1 - y3 );
    }
    if ( TH2D * estimateHisto3 = dynamic_cast < TH2D* > ( _rootObjectMap["EstimationMethodDifferenceDAFTriplet"] ) )
    {
	estimateHisto3 -> Fill ( x3 - x2, y3 - y2 );
    }
}

// resi vs tdc
void EUTelFitTupleAlibava::fillTDCResHisto ( float hitx, float hity, float trackx, float tracky, float tdctime )
{
    if ( TH2D * resiXHisto = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualXTDC"] ) )
    {
	resiXHisto -> Fill ( hitx - trackx, tdctime );
    }
    if ( TH2D * resiYHisto = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualYTDC"] ) )
    {
	resiYHisto -> Fill ( hity - tracky, tdctime );
    }
}

// the residual plots
// this is called with the matched DUT hit in X and Y, as well as with the track hit X and Y...
// all should be in mm and in the global coordinate system
void EUTelFitTupleAlibava::fillresihistos ( float hitx, float hity, float trackx, float tracky, int eventnr, double chi2, double ndf )
{
    if ( TH1D * resiXHisto = dynamic_cast < TH1D* > ( _rootObjectMap["ResidualX"] ) )
    {
	resiXHisto -> Fill ( hitx - trackx );
    }
    if ( TH1D * resiYHisto = dynamic_cast < TH1D* > ( _rootObjectMap["ResidualY"] ) )
    {
	resiYHisto -> Fill ( hity - tracky );
    }
    if ( TH2D * DYresiXHisto = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualDXX"] ) )
    {
	DYresiXHisto -> Fill ( hitx - trackx, trackx );
    }
    if ( TH2D * DYresiYHisto = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualDXY"] ) )
    {
	DYresiYHisto -> Fill ( hitx - trackx, tracky );
    }
    if ( TProfile* ResidualProfileDX_X = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDX_X"] ) )
    {
	ResidualProfileDX_X -> Fill ( trackx, hitx - trackx, 1 );
    }
    if ( TProfile* ResidualProfileDX_Y = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDX_Y"] ) )
    {
	ResidualProfileDX_Y -> Fill ( tracky, hitx - trackx, 1 );
    }
    if ( TH2D * DYresiXHisto = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualDYX"] ) )
    {
	DYresiXHisto -> Fill ( hity - tracky, trackx );
    }
    if ( TH2D * DYresiYHisto = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualDYY"] ) )
    {
	DYresiYHisto -> Fill ( hity - tracky, tracky );
    }
    if ( TProfile* ResidualProfileDY_X = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDY_X"] ) )
    {
	ResidualProfileDY_X -> Fill ( trackx, hity - tracky, 1 );
    }
    if ( TProfile* ResidualProfileDY_Y = dynamic_cast < TProfile* > ( _rootObjectMap["ResidualProfileDY_Y"] ) )
    {
	ResidualProfileDY_Y -> Fill ( tracky, hity - tracky, 1 );
    }
    if ( TH2D * resiXHistoTime = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualXTime"] ) )
    {
	resiXHistoTime -> Fill ( eventnr, hitx - trackx );
    }
    if ( TH2D * resiYHistoTime = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualYTime"] ) )
    {
	resiYHistoTime -> Fill ( eventnr, hity - tracky );
    }
    if ( TH2D * resiXHistoChi2 = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualXChi2"] ) )
    {
	resiXHistoChi2 -> Fill ( chi2 / ndf, hitx - trackx );
    }
    if ( TH2D * resiYHistoChi2 = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualYChi2"] ) )
    {
	resiYHistoChi2 -> Fill ( chi2 / ndf, hity - tracky );
    }
    double prob = TMath::Prob ( chi2, ndf );
    if ( TH2D * resiXHistoProb = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualXProb"] ) )
    {
	resiXHistoProb -> Fill ( prob, hitx - trackx );
    }
    if ( TH2D * resiYHistoProb = dynamic_cast < TH2D* > ( _rootObjectMap["ResidualYProb"] ) )
    {
	resiYHistoProb -> Fill ( prob, hity - tracky );
    }
}

// residual histograms after the fiducial cuts
// all should be in mm and in the global coordinate system
void EUTelFitTupleAlibava::Filteredfillresihistos ( float hitx, float hity, float trackx, float tracky )
{
    if ( TH1D * FilteredresiXHisto = dynamic_cast < TH1D* > ( _rootObjectMap["FilteredResidualX"] ) )
    {
	FilteredresiXHisto -> Fill ( hitx - trackx );
    }
    if ( TH1D * FilteredresiYHisto = dynamic_cast < TH1D* > ( _rootObjectMap["FilteredResidualY"] ) )
    {
	FilteredresiYHisto -> Fill ( hity - tracky );
    }
    if ( TH2D * FilteredDYresiXHisto = dynamic_cast < TH2D* > ( _rootObjectMap["FilteredResidualDYX"] ) )
    {
	FilteredDYresiXHisto -> Fill ( hity - tracky, trackx );
    }
    if ( TH2D * FilteredDYresiYHisto = dynamic_cast < TH2D* > ( _rootObjectMap["FilteredResidualDYY"] ) )
    {
	FilteredDYresiYHisto -> Fill ( hity - tracky, tracky );
    }
    if ( TProfile* FilteredResidualProfileDY_X = dynamic_cast < TProfile* > ( _rootObjectMap["FilteredResidualProfileDY_X"] ) )
    {
	FilteredResidualProfileDY_X -> Fill ( trackx, hity - tracky, 1 );
    }
    if ( TProfile* FilteredResidualProfileDY_Y = dynamic_cast < TProfile* > ( _rootObjectMap["FilteredResidualProfileDY_Y"] ) )
    {
	FilteredResidualProfileDY_Y -> Fill ( tracky, hity - tracky, 1 );
    }
}

// some control plots...
// if _checkdealignment is true, these histograms will be filled with the difference between applied alignment in reverse and our calculation...
void EUTelFitTupleAlibava::filldebughistos ( float x1, float x2, float y1, float y2 )
{
    if ( TH1D * calcX = dynamic_cast < TH1D* > ( _rootObjectMap["CalcX"] ) )
    {
	calcX -> Fill ( x1 - x2 );
    }
    if ( TH1D * calcY = dynamic_cast < TH1D* > ( _rootObjectMap["CalcY"] ) )
    {
	calcY -> Fill ( y1 - y2 );
    }

    if ( ( fabs ( x1 - x2 ) > 5 ) || ( fabs ( y1 - y2 ) > 5 ) )
    {
	streamlog_out ( WARNING5 ) << "Attention! Dealignment not working! X: " << x1 - x2 << " mm, Y: " << y1 - y2 << " mm!" << endl;
    }
}

// fill track based eta plots
// this is called if there is an Alibava DUT present to check the charge distribution relative to the track point
// left and right are used as the local track point reconstruction can be biased because of alignment errors
void EUTelFitTupleAlibava::fillEtaHisto ( float left, float right, float left1, float right1, float left2, float right2, float local, int event, float tdctime )
{
    if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaCentral"] ) )
    {
	eta -> Fill ( left / ( left + right ) );
    }
    if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaShiftedLeft"] ) )
    {
	eta -> Fill ( left1 / ( left1 + left ) );
    }
    if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaShiftedRight"] ) )
    {
	eta -> Fill ( right / ( right + right1 ) );
    }

    // if a track is actually "on" a strip
    double fractpart, intpart;
    fractpart = modf ( local , &intpart );

    // small: strip is left
    if ( fractpart < 0.25 )
    {
	if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaTrackOnStrip"] ) )
	{
	    eta -> Fill ( left1 / ( left1 + right ) );
	}
	if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaTrackOnStripLeftCharge"] ) )
	{
	    eta -> Fill ( left1 / left );
	}
	if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaTrackOnStripRightCharge"] ) )
	{
	    eta -> Fill ( right / left );
	}
    }

    // large: strip is right
    if ( fractpart > 0.75 )
    {
	if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaTrackOnStrip"] ) )
	{
	    eta -> Fill ( left / ( left + right1 ) );
	}
	if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaTrackOnStripLeftCharge"] ) )
	{
	    eta -> Fill ( left / right );
	}
	if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackEtaTrackOnStripRightCharge"] ) )
	{
	    eta -> Fill ( right1 / right );
	}
    }

    if ( TH2D * eta = dynamic_cast < TH2D* > ( _rootObjectMap["2DEta"] ) )
    {
	    eta -> Fill ( left / ( left + right ), fractpart );
    }

    if ( TH2D * etat = dynamic_cast < TH2D* > ( _rootObjectMap["TrackEtaOverEvents"] ) )
    {
	    etat -> Fill ( left / ( left + right ), event );
    }

    // signal of 5 strips under a track
    // depending on fractpart select left2 or right2
    float signal = 0.0;
    float leftneighbour = 0.0;
    float rightneighbour = 0.0;
    float seedsignal = 0.0;
    if ( fractpart < 0.5 )
    {
	signal = ( left2 + left1 + left + right + right1 ) * _polarity;
	seedsignal = left * _polarity;
	rightneighbour = right * _polarity;
	leftneighbour = left1 * _polarity;
    }
    if ( fractpart == 0.5 )
    {
	signal = ( ( left2 / 2.0 ) + left1 + left + right + right1 + ( right2 / 2.0 ) ) * _polarity;
	seedsignal = ( ( left + right ) / 2.0 ) * _polarity;
	rightneighbour = ( ( right + right1 ) / 2.0 ) * _polarity;
	leftneighbour = ( ( left1 + left ) / 2.0 ) * _polarity;
    }
    if ( fractpart > 0.5 )
    {
	signal = ( left1 + left + right + right1 + right2 ) * _polarity;
	seedsignal = right * _polarity;
	rightneighbour = right1 * _polarity;
	leftneighbour = left * _polarity;
    }

    if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackSignal"] ) )
    {
	eta -> Fill (signal);
    }

    // signal vs tdc
    if ( TH2D * eta = dynamic_cast < TH2D* > ( _rootObjectMap["TrackSignalTDC"] ) )
    {
	eta -> Fill ( tdctime, signal );
    }

    // seed and neighbour signal vs tdc
    if ( TH2D * eta = dynamic_cast < TH2D* > ( _rootObjectMap["SeedSignalTDC"] ) )
    {
	eta -> Fill ( tdctime, seedsignal );
    }
    if ( TH2D * eta = dynamic_cast < TH2D* > ( _rootObjectMap["LeftSignalTDC"] ) )
    {
	eta -> Fill ( tdctime, leftneighbour );
    }
    if ( TH2D * eta = dynamic_cast < TH2D* > ( _rootObjectMap["RightSignalTDC"] ) )
    {
	eta -> Fill ( tdctime, rightneighbour );
    }
}

// check how matched hits are distributed between strips
// this xy is the local pixel of the hit
void EUTelFitTupleAlibava::fillInterStripHitmap ( double x, double y )
{
    double fractpart = 0.0;
    double intpart = 0.0;
    if ( _unsensitiveaxis == "x" )
    {
	fractpart = modf ( y , &intpart );
    }
    if ( _unsensitiveaxis == "y" )
    {
	fractpart = modf ( x , &intpart );
    }

    if ( TH1D * eta = dynamic_cast < TH1D* > ( _rootObjectMap["TrackDistribution"] ) )
    {
	eta -> Fill ( fractpart );
    }
}

// begin histogram booking
// here we define the names of the ntuple colums we want to fill
// the order is important and should correspond to the one used by the processor above
// also all other histograms get booked here
void EUTelFitTupleAlibava::bookHistos ( )
{
    streamlog_out ( MESSAGE5 ) << "Booking fit n-tuple" << endl;

    std::vector < std::string > _columnNames;
    std::vector < std::string > _columnType;

    _columnNames.push_back ( "Event" );
    _columnType.push_back ( "int" );

    _columnNames.push_back ( "RunNr" );
    _columnType.push_back ( "int" );

    _columnNames.push_back ( "EvtNr" );
    _columnType.push_back ( "int" );

    _columnNames.push_back ( "Ndf" );
    _columnType.push_back ( "int" );

    _columnNames.push_back ( "Chi2" );
    _columnType.push_back ( "float" );

    const char * _varName[] = { "measX", "measY" , "measZ", "measQ", "fitX", "fitY", "fitZ" };

    for ( int ipl = 0; ipl < _nTelPlanes; ipl++ )
    {
	for ( int ivar = 0; ivar < 7; ivar++ )
	{
	    stringstream st;
	    st << _varName[ivar] << "_" << ipl;
	    _columnNames.push_back ( st.str ( ) );
	    _columnType.push_back ( "double" );
	}
    }

    // the DUT fields get filled by _missingValue if there is no DUT present

    _columnNames.push_back ( "dutTrackX_global" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackY_global" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackZ_global" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackX_local" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackY_local" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackZ_local" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackX_pixel" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackY_pixel" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutTrackZ_pixel" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitX_global" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitY_global" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitZ_global" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitX_local" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitY_local" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitZ_local" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitX_pixel" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitY_pixel" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitZ_pixel" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitR" );
    _columnType.push_back ( "double" );

    _columnNames.push_back ( "dutHitQ" );
    _columnType.push_back ( "double" );

    // these fields have the alibava information
    if ( _alibava == true && _iDUT != -1 )
    {

	// the header information
	_columnNames.push_back ( "alibava_tdc" );
	_columnType.push_back ( "float" );

	_columnNames.push_back ( "alibava_temp" );
	_columnType.push_back ( "float" );

	// FIXME assumes 2 chips with 128 channels each
	for (int i = 0; i < 256; i++ )
	{
	    stringstream st;
	    st << "alibava_reco_ch_" << i;
	    _columnNames.push_back ( st.str ( ) );
	    _columnType.push_back ( "double" );
	}
    }

    _FitTuple = AIDAProcessor::tupleFactory ( this ) -> create ( _FitTupleName, _FitTupleName, _columnNames, _columnType, "" );

    // done with the actual ntuple
    // now the histograms, order is no longer important

    // plot some DUT information/control plots
    if ( _iDUT != -1 )
    {
	TH2D * precHisto = new TH2D ( "Precision", "", 500, -50, 50, 500, -50, 50 );
	_rootObjectMap.insert ( make_pair ( "Precision", precHisto ) );
	precHisto -> SetTitle ( "Impact Prediction Precision;X_{reconstructed - clustered} [Pitch];Y_{reconstructed - clustered} [Pitch]" );

	TH2D * estimateHisto1 = new TH2D ( "EstimationMethodDifferenceIntersectTriplet", "", 500, -10, 10, 500, -10, 10 );
	_rootObjectMap.insert ( make_pair ( "EstimationMethodDifferenceIntersectTriplet", estimateHisto1 ) );
	estimateHisto1 -> SetTitle ( "Estimation Method Difference Between Intersect and Triplet;X_{Plane intersection - triplet} [mm];Y_{Plane intersection - triplet} [mm]" );

	TH2D * estimateHisto2 = new TH2D ( "EstimationMethodDifferenceIntersectDAF", "", 500, -10, 10, 500, -10, 10 );
	_rootObjectMap.insert ( make_pair ( "EstimationMethodDifferenceIntersectDAF", estimateHisto2 ) );
	estimateHisto2 -> SetTitle ( "Estimation Method Difference Between Intersect and DAF;X_{Plane intersection - DAF} [mm];Y_{Plane intersection - DAF} [mm]" );

	TH2D * estimateHisto3 = new TH2D ( "EstimationMethodDifferenceDAFTriplet", "", 500, -10, 10, 500, -10, 10 );
	_rootObjectMap.insert ( make_pair ( "EstimationMethodDifferenceDAFTriplet", estimateHisto3 ) );
	estimateHisto3 -> SetTitle ( "Estimation Method Difference Between DAF and Triplet;X_{Plane DAF - triplet} [mm];Y_{Plane DAF - triplet} [mm]" );

	TH1D * resiXHisto = new TH1D ( "ResidualX", "", 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualX", resiXHisto ) );
	resiXHisto -> SetTitle ( "DUT Residual in X;X_{hit - track} [mm];Events" );

	TH1D * resiYHisto = new TH1D ( "ResidualY", "", 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualY", resiYHisto ) );
	resiYHisto -> SetTitle ( "DUT Residual in Y;Y_{hit - track} [mm];Events" );

	TF1 * resiXFit = new TF1 ( "ResidualXFit", "gaus" );
	_rootObjectMap.insert ( make_pair ( "ResidualXFit", resiXFit) );

	TF1 * resiYFit = new TF1 ( "ResidualYFit", "gaus" );
	_rootObjectMap.insert ( make_pair ( "ResidualYFit", resiYFit) );

	TH2D * DXresiXHisto = new TH2D ( "ResidualDXX", "", 500, -0.5, 0.5, 500, -15.0, 15.0 );
	_rootObjectMap.insert ( make_pair ( "ResidualDXX", DXresiXHisto ) );
	DXresiXHisto -> SetTitle ( "DUT Residual in X vs. X;X_{hit - track} [mm];X_{track} [mm]" );

	TH2D * DXresiYHisto = new TH2D ( "ResidualDXY", "", 500, -0.5, 0.5, 500, -10.0, 10.0 );
	_rootObjectMap.insert ( make_pair ( "ResidualDXY", DXresiYHisto ) );
	DXresiYHisto -> SetTitle ( "DUT Residual in X vs. Y;X_{hit - track} [mm];Y_{track} [mm]" );

	TProfile * ResidualProfileDX_X = new TProfile ( "ResidualProfileDX_X", "", 150, -15, 15, -0.2, 0.2, "s" );
	_rootObjectMap.insert ( make_pair ( "ResidualProfileDX_X", ResidualProfileDX_X ) );
	ResidualProfileDX_X -> SetTitle ( "DUT Residual in X vs. X - Profile;X_{track} [mm];X_{hit - track} [mm]" );

	TProfile * ResidualProfileDX_Y = new TProfile ( "ResidualProfileDX_Y", "", 150, -15, 15, -0.2, 0.2, "s" );
	_rootObjectMap.insert ( make_pair ( "ResidualProfileDX_Y", ResidualProfileDX_Y ) );
	ResidualProfileDX_Y -> SetTitle ( "DUT Residual in X vs. Y - Profile;Y_{track} [mm];X_{hit - track} [mm]" );

	TH2D * DYresiXHisto = new TH2D ( "ResidualDYX", "", 500, -0.5, 0.5, 500, -15.0, 15.0 );
	_rootObjectMap.insert ( make_pair ( "ResidualDYX", DYresiXHisto ) );
	DYresiXHisto -> SetTitle ( "DUT Residual in Y vs. X;Y_{hit - track} [mm];X_{track} [mm]" );

	TH2D * DYresiYHisto = new TH2D ( "ResidualDYY", "", 500, -0.5, 0.5, 500, -10.0, 10.0 );
	_rootObjectMap.insert ( make_pair ( "ResidualDYY", DYresiYHisto ) );
	DYresiYHisto -> SetTitle ( "DUT Residual in Y vs. Y;Y_{hit - track} [mm];Y_{track} [mm]" );

	TProfile * ResidualProfileDY_X = new TProfile ( "ResidualProfileDY_X", "", 150, -15, 15, -0.2, 0.2, "s" );
	_rootObjectMap.insert ( make_pair ( "ResidualProfileDY_X", ResidualProfileDY_X ) );
	ResidualProfileDY_X -> SetTitle ( "DUT Residual in Y vs. X - Profile;X_{track} [mm];Y_{hit - track} [mm]" );

	TProfile * ResidualProfileDY_Y = new TProfile ("ResidualProfileDY_Y", "", 150, -15, 15, -0.2, 0.2, "s" );
	_rootObjectMap.insert ( make_pair ( "ResidualProfileDY_Y", ResidualProfileDY_Y ) );
	ResidualProfileDY_Y -> SetTitle ( "DUT Residual in Y vs. Y - Profile;Y_{track} [mm];Y_{hit - track} [mm]" );

	TH2D * resiXHistoTime = new TH2D ( "ResidualXTime", "", 1000, 0, 500000, 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualXTime", resiXHistoTime ) );
	resiXHistoTime -> SetTitle ( "DUT Residual in X over event number;Event Number;X_{hit - fit} [mm]" );

	TH2D * resiYHistoTime = new TH2D ( "ResidualYTime", "", 1000, 0, 500000, 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualYTime", resiYHistoTime) );
	resiYHistoTime -> SetTitle ( "DUT Residual in Y over event number;Event Number;Y_{hit - fit} [mm]" );

	TH2D * resiXHistoChi2 = new TH2D ( "ResidualXChi2", "", 100, 0, 20, 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualXChi2", resiXHistoChi2 ) );
	resiXHistoChi2 -> SetTitle ( "DUT Residual in X over Track Chi2/Ndf;Track Chi2/Ndf;X_{hit - fit} [mm]" );

	TH2D * resiYHistoChi2 = new TH2D ( "ResidualYChi2", "", 100, 0, 20, 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualYChi2", resiYHistoChi2 ) );
	resiYHistoChi2 -> SetTitle ( "DUT Residual in Y over Track Chi2/Ndf;Track Chi2/Ndf;Y_{hit - fit} [mm]" );

	TH2D * resiXHistoProb = new TH2D ( "ResidualXProb", "", 100, 0, 1, 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualXProb", resiXHistoProb ) );
	resiXHistoProb -> SetTitle ( "DUT Residual in X over Track Probability;Track Probability;X_{hit - fit} [mm]" );

	TH2D * resiYHistoProb = new TH2D ( "ResidualYProb", "", 100, 0, 1, 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "ResidualYProb", resiYHistoProb ) );
	resiYHistoProb -> SetTitle ( "DUT Residual in Y over Track Probability;Track Probability;Y_{hit - fit} [mm]" );

	TH1D * calcX = new TH1D ( "CalcX", "", 500, -5, 5 );
	_rootObjectMap.insert ( make_pair ( "CalcX", calcX ) );
	calcX -> SetTitle ( "De-alignment Calculation Difference in X;X_{Reconstruction - Alignment} [pitch];Events" );

	TH1D * calcY = new TH1D ( "CalcY", "", 500, -5, 5 );
	_rootObjectMap.insert ( make_pair ( "CalcY", calcY ) );
	calcY -> SetTitle ( "De-alignment Calculation Difference in Y;Y_{Reconstruction - Alignment} [pitch];Events" );

	TF1 * profilefitX = new TF1 ( "ProfileXFit", "pol1" );
	_rootObjectMap.insert ( make_pair ( "ProfileXFit", profilefitX ) );

	TF1 * profilefitY = new TF1 ( "ProfileYFit", "pol1" );
	_rootObjectMap.insert ( make_pair ( "ProfileYFit", profilefitY ) );

	TH1D * FilteredresiXHisto = new TH1D ( "FilteredResidualX", "", 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "FilteredResidualX", FilteredresiXHisto ) );
	FilteredresiXHisto -> SetTitle ( "Filtered DUT Residual in X;X_{hit - track} [mm];Events" );

	TH1D * FilteredresiYHisto = new TH1D ( "FilteredResidualY", "", 500, -0.5, 0.5 );
	_rootObjectMap.insert ( make_pair ( "FilteredResidualY", FilteredresiYHisto ) );
	FilteredresiYHisto -> SetTitle ( "Filtered DUT Residual in Y;Y_{hit - track} [mm];Events" );
	
	TH2D * FilteredDYresiXHisto = new TH2D ( "FilteredResidualDYX", "", 500, -0.5, 0.5, 500, -15.0, 15.0 );
	_rootObjectMap.insert ( make_pair ( "FilteredResidualDYX", FilteredDYresiXHisto ) );
	FilteredDYresiXHisto -> SetTitle ( "Filtered DUT Residual in Y vs. X;Y_{hit - track} [mm];X_{track} [mm]" );

	TH2D * FilteredDYresiYHisto = new TH2D ( "FilteredResidualDYY", "", 500, -0.5, 0.5, 500, -10.0, 10.0 );
	_rootObjectMap.insert ( make_pair ( "FilteredResidualDYY", FilteredDYresiYHisto ) );
	FilteredDYresiYHisto -> SetTitle ( "FilteredDUT Residual in Y vs. Y;Y_{hit - track} [mm];Y_{track} [mm]" );

	TProfile * FilteredResidualProfileDY_X = new TProfile ( "FilteredResidualProfileDY_X", "", 150, -15, 15, -0.2, 0.2, "s" );
	_rootObjectMap.insert ( make_pair ( "FilteredResidualProfileDY_X", FilteredResidualProfileDY_X ) );
	FilteredResidualProfileDY_X -> SetTitle ( "Filtered DUT Residual in Y vs. X - Profile;X_{track} [mm];Y_{hit - track} [mm]" );

	TProfile * FilteredResidualProfileDY_Y = new TProfile ( "FilteredResidualProfileDY_Y", "", 150, -15, 15, -0.2, 0.2, "s" );
	_rootObjectMap.insert ( make_pair ( "FilteredResidualProfileDY_Y", FilteredResidualProfileDY_Y ) );
	FilteredResidualProfileDY_Y -> SetTitle ( "Filtered DUT Residual in Y vs. Y - Profile;Y_{track} [mm];Y_{hit - track} [mm]" );

	TH2D * hitmapHisto = new TH2D ( "TrackHitmap", "", 4000, -20, 20, 4000, -20, 20 );
	_rootObjectMap.insert ( make_pair ( "TrackHitmap", hitmapHisto ) );
	hitmapHisto -> SetTitle ( "Track Point Hitmap on DUT;Track X [mm];Track Y [mm]" );

	TH2D * hitmapHistoMatch = new TH2D ( "TrackHitmapMatch", "", 4000, -20, 20, 4000, -20, 20 );
	_rootObjectMap.insert ( make_pair ( "TrackHitmapMatch", hitmapHistoMatch ) );
	hitmapHistoMatch -> SetTitle ( "Track Point Hitmap on DUT if Matched;Track X [mm];Track Y [mm]" );

	TH2D * hitmapHistoModPitch = new TH2D ( "TrackHitmapModPitch", "", 50, -1, 2, 50, -1, 2 );
	_rootObjectMap.insert ( make_pair ( "TrackHitmapModPitch", hitmapHistoModPitch ) );
	hitmapHistoModPitch -> SetTitle ( "Track Point Hitmap on DUT;Track X mod Pitch;Track Y mod Pitch" );

	TH2D * hitmapHistoMatchModPitch = new TH2D ( "TrackHitmapMatchModPitch", "", 50, -1, 2, 50, -1, 2 );
	_rootObjectMap.insert ( make_pair ( "TrackHitmapMatchModPitch", hitmapHistoMatchModPitch ) );
	hitmapHistoMatchModPitch -> SetTitle ( "Track Point Hitmap on DUT if Matched;Track X mod Pitch;Track Y mod Pitch" );

	TH3D * hitmapallhits = new TH3D ( "DUTHitmap", "", 200, -10, 10, 200, -10, 10, 200, -10, 10 );
	_rootObjectMap.insert ( make_pair ( "DUTHitmap", hitmapallhits ) );
	hitmapallhits -> SetTitle ( "DUT Hitmap of All Hits;X [mm];relative Z [mm];Y [mm]" );

	TH2D * hitmapallhitsmodpitch = new TH2D ( "DUTHitmapModPitch", "", 50, -1, 2, 50, -1, 2 );
	_rootObjectMap.insert ( make_pair ( "DUTHitmapModPitch", hitmapallhitsmodpitch ) );
	hitmapallhitsmodpitch -> SetTitle ( "DUT Hitmap of All Hits mod pitch;Hit X [Pitch];Hit Y [Pitch]" );

	TH3D * hitmapmatchedhits = new TH3D ( "DUTHitmapMatched", "", 200, -10, 10, 200, -10, 10, 200, -10, 10 );
	_rootObjectMap.insert ( make_pair ( "DUTHitmapMatched", hitmapmatchedhits ) );
	hitmapmatchedhits -> SetTitle ( "DUT Hitmap of All Matched Hits;X [mm];relative Z [mm];Y [mm]" );

	TH2D * hitmapmatchedhitsmodpitch = new TH2D ( "DUTHitmapMatchedModPitch", "", 50, -1, 2, 50, -1, 2 );
	_rootObjectMap.insert ( make_pair ( "DUTHitmapMatchedModPitch", hitmapmatchedhitsmodpitch ) );
	hitmapmatchedhitsmodpitch -> SetTitle ( "DUT Hitmap of All Matched Hits mod pitch;Hit X [Pitch];Hit Y [Pitch]" );

	// eta is only filled if there is an alibava present
	if ( _alibava == true )
	{

	    TH1D * etac = new TH1D ( "TrackEtaCentral", "", 120, -1, 2 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaCentral", etac ) );
	    etac -> SetTitle ( "Track Eta Central;Eta unit;Events" );

	    TH1D * etal = new TH1D ( "TrackEtaShiftedLeft", "", 120, -1, 2 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaShiftedLeft", etal ) );
	    etal -> SetTitle ( "Track Eta if track point is shifted left;Eta unit;Events" );

	    TH1D * etar = new TH1D ( "TrackEtaShiftedRight", "", 120, -1, 2 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaShiftedRight", etar ) );
	    etar -> SetTitle ( "Track Eta if track point is shifted right;Eta unit;Events" );

	    TH2D * eta2d = new TH2D ( "2DEta", "", 120, -1, 2, 100, -0.5, 1.5 );
	    _rootObjectMap.insert ( make_pair ( "2DEta", eta2d ) );
	    eta2d -> SetTitle ( "Eta 2D;Eta;Y mod pitch" );

	    TH1D * etac2 = new TH1D ( "TrackEtaTrackOnStrip", "", 120, -2, 3 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaTrackOnStrip", etac2 ) );
	    etac2 -> SetTitle ( "Track Eta Central;Eta unit;Events" );

	    TH1D * etal2 = new TH1D ( "TrackEtaTrackOnStripLeftCharge", "", 120, -2, 2 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaTrackOnStripLeftCharge", etal2 ) );
	    etal2 -> SetTitle ( "Track Eta Left;Relative Charge;Events" );

	    TH1D * etar2 = new TH1D ( "TrackEtaTrackOnStripRightCharge", "", 120, -2, 2 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaTrackOnStripRightCharge", etar2 ) );
	    etar2 -> SetTitle ( "Track Eta Right;Relative Charge;Events" );

	    TH2D * etat = new TH2D ( "TrackEtaOverEvents", "", 120, -1, 2, 1000, 0, 500000 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaOverEvents", etat ) );
	    etat -> SetTitle ( "Track Eta vs. Time;Eta unit;Eventnr." );

	    TH1D * etaintegral = new TH1D ( "TrackEtaIntegral", "", 120, -1, 2 );
	    _rootObjectMap.insert ( make_pair ( "TrackEtaIntegral", etaintegral ) );
	    etaintegral -> SetTitle ( "Track Eta Integral;Eta unit;Events" );

	    TH1D * trackdistri = new TH1D ( "TrackDistribution", "", 200, -0.5, 1.5 );
	    _rootObjectMap.insert ( make_pair ( "TrackDistribution", trackdistri ) );
	    trackdistri -> SetTitle ( "Matched Hit Distribution Between Strips;Pitch;Events" );

	    TH1D * tracksignal = new TH1D ( "TrackSignal", "", 2000, -500, 500 );
	    _rootObjectMap.insert ( make_pair ( "TrackSignal", tracksignal ) );
	    tracksignal -> SetTitle ( "Signal in 5 Strips Under Track;ADCs;Events" );

	    TH2D * tracksignaltdc = new TH2D ( "TrackSignalTDC", "", 101, 0, 100, 2000, -500, 500 );
	    _rootObjectMap.insert ( make_pair ( "TrackSignalTDC", tracksignaltdc ) );
	    tracksignaltdc -> SetTitle ( "Signal in 5 Strips Under Track vs. TDC Time;TDC [ns];Signal [ADCs]" );

	    TH2D * seedsignaltdc = new TH2D ( "SeedSignalTDC", "", 101, 0, 100, 2000, -500, 500 );
	    _rootObjectMap.insert ( make_pair ( "SeedSignalTDC", seedsignaltdc ) );
	    seedsignaltdc -> SetTitle ( "Signal in Seed Strip Under Track vs. TDC Time;TDC [ns];Signal [ADCs]" );

	    TH2D * leftsignaltdc = new TH2D ( "LeftSignalTDC", "", 101, 0, 100, 2000, -500, 500 );
	    _rootObjectMap.insert ( make_pair ( "LeftSignalTDC", leftsignaltdc ) );
	    leftsignaltdc -> SetTitle ( "Signal in Left Neighbour Strip Under Track vs. TDC Time;TDC [ns];Signal [ADCs]" );

	    TH2D * rightsignaltdc = new TH2D ( "RightSignalTDC", "", 101, 0, 100, 2000, -500, 500 );
	    _rootObjectMap.insert ( make_pair ( "RightSignalTDC", rightsignaltdc ) );
	    rightsignaltdc -> SetTitle ( "Signal in Right Neighbour Strip Under Track vs. TDC Time;TDC [ns];Signal [ADCs]" );

	    TH2D * resxtdc = new TH2D ( "ResidualXTDC", "", 500, -0.5, 0.5, 100, 0, 100 );
	    _rootObjectMap.insert ( make_pair ( "ResidualXTDC", resxtdc ) );
	    resxtdc -> SetTitle ( "DUT Residual in X over Cluster TDC;X_{hit - fit} [mm];TDC [ns]" );

	    TH2D * resytdc = new TH2D ( "ResidualYTDC", "", 500, -0.5, 0.5, 100, 0, 100 );
	    _rootObjectMap.insert ( make_pair ( "ResidualYTDC", resytdc ) );
	    resytdc -> SetTitle ( "DUT Residual in Y over Cluster TDC;Y_{hit - fit} [mm];TDC [ns]" );

	}

    }

    streamlog_out ( DEBUG5 ) << "Booking completed" << endl;

    return;
}

#endif // GEAR && AIDA

void EUTelFitTupleAlibava::dolandaugausfit ( string tempHistoName )
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
	TF1 *fitsnr = langaufit2 ( histo, fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf );
	histo -> Fit ( fitsnr, "QR" );

	Double_t SNRPeak, SNRFWHM;
	langaupro2 ( fp, SNRPeak, SNRFWHM );

	streamlog_out ( MESSAGE4 ) << "===============================" << endl;
	streamlog_out ( DEBUG5 ) << "Done fitting landau gaus on histo" << tempHistoName << endl;
	streamlog_out ( MESSAGE4 ) << "Landau-Gaus Peak for : " << tempHistoName << " is: " << SNRPeak << endl;
	streamlog_out ( MESSAGE4 ) << "Landau-Gaus Sigma for: " << tempHistoName << " is: " << SNRFWHM << endl;
	streamlog_out ( MESSAGE4 ) << "===============================" << endl;
    }

}

Double_t langaufun2 ( Double_t *x, Double_t *par )
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

TF1 *langaufit2 ( TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF )
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

    TF1 *ffitold = static_cast<TF1*>(gROOT -> GetListOfFunctions ( ) -> FindObject ( FunName ));
    if ( ffitold ) delete ffitold;

    TF1 *ffit = new TF1 ( FunName, langaufun2, fitrange[0], fitrange[1], 4 );
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

Int_t langaupro2 ( Double_t *params, Double_t &maxx, Double_t &FWHM )
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
	l = langaufun2 ( &x, params );

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
	l = TMath::Abs ( langaufun2 ( &x, params ) - fy );

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
	l = TMath::Abs ( langaufun2 ( &x, params ) - fy );

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


/*
 * 
 * Rescue Pig to the rescue!

                                 _
    _._ _..._ .-',     _.._(`) )
   '-. `     '  /-._.-'    ',/
      )         \            '.
     / _    _    |             \
    |  a    a    /              |
    \   .-.                     ;  
     '-('' ).-'       ,'       ;
        '-;           |      .'
           \           \    /
           | 7  .__  _.-\   \
           | |  |  ``/  /`  /
          /,_|  |   /,_/   /
             /,_/      '`-'
  
  
  */
