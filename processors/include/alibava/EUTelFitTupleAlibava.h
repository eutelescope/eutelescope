/*
 * Rewritten and alibava features added by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 *
 *  Author  A.F.Zarnecki
 * (2007 University of Warsaw)
 *
 *  email:zarnecki@fuw.edu.pl
 */

#ifndef EUTelFitTupleAlibava_h
#define EUTelFitTupleAlibava_h 1

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	#include <AIDA/IBaseHistogram.h>
	#include <AIDA/ITuple.h>
#endif

// ROOT includes <>
#include "TObject.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>

// system includes <>
#include <string>
#include <vector>
#include <map>

#include "TObject.h"

using namespace std;

    // This is used for the landau gaus function and fits. Details see: http://root.cern.ch/root/html/tutorials/fit/langaus.C.html
    Double_t langaufun2 ( Double_t *x, Double_t *par );
    TF1 *langaufit2 ( TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF );
    Int_t langaupro2 ( Double_t *params, Double_t &maxx, Double_t &FWHM );

namespace eutelescope
{

    // Processor for building n-tuple with track fit results

    class EUTelFitTupleAlibava : public marlin::Processor
    {

	public:

	    virtual Processor*  newProcessor ( )
	    {
		return new EUTelFitTupleAlibava;
	    }

	    EUTelFitTupleAlibava ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader* run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    std::map < std::string , TObject * > _rootObjectMap;

	    virtual void end ( );

	protected:

	    // Silicon plane parameters as described in GEAR. This object is provided by GEAR during the init() phase and stored here for local use.
	    gear::SiPlanesParameters * _siPlanesParameters;

	    // This is the real geometry description for each layer. This object is taken from _siPlanesParameters during the init() phase and stored for local use
	    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

	    // Input Track collection name
	    std::string _inputColName;

	    // Input TrackerHit collection name
	    std::string _inputDUTColName;

	    // Write our own alignment into this file
	    std::string _outputalignment;

	    bool _doAlignment;

	    bool _doGamma;

	    // Flag for manual DUT selection
	    bool _useManualDUT;

	    // Id of telescope layer which should be used as a DUT
	    int _manualDUTid;

	    //  Value to be used for missing measurements when filling the tuple
	    double _missingValue;

	    // Flag to recheck alignment of the DUT
	    bool _checkdealignment;

	    // Alibava output and usage on/off
	    bool _alibava;

	    // Alibava collection name
	    std::string _alibavaCollectionName;

	    // Alibava cluster collection name
	    std::string _alibavaClusterCollectionName;

	    // Alibava unfiltered collection name
	    std::string _alibavaUnfilteredCollectionName;

	    std::string _alignedHitPrefix;
	    std::string _prealignedHitPrefix;
	    std::string _unalignedcollectionname;

	    // Alignment collection name
	    std::vector<std::string> _alignmentCollectionName;

	    // Prealignment collection name
	    std::vector<std::string> _pre_alignmentCollectionName;

	    // Original hit collection name
	    std::string _originalCollectionName;

	    // Original prealignment applied to the original hit collection
	    std::string _originalPreAlignment;

	    bool _useOriginalPreAlignment;

	    // Unsensitive strip sensor axis
	    std::string _unsensitiveaxis;

	    // Telescope plane count
	    int _nTelPlanes;

	    // Telescope plane sorting
	    int * _planeSort;

	    // Telescope plane id
	    int * _planeID;

	    // Telescope plane position in z
	    double * _planePosition;

	    // Flag for active plane
	    bool   * _isActive;

	    // Flag for measurements...
	    bool   * _isMeasured;

	    // ... in X
	    double * _measuredX;

	    // ... in Y
	    double * _measuredY;

	    // ... in Z
	    double * _measuredZ;

	    // ... with charge
	    double * _measuredQ;

	    // Flag for fitted measurements...
	    bool   * _isFitted;

	    // ... in X
	    double * _fittedX;

	    // ... in Y
	    double * _fittedY;

	    // ... in Z
	    double * _fittedZ;

	    // The telescope coordinates for extrapolation
	    double * _telescopeX;

	    double * _telescopeY;

	    // The DUT plane
	    int _iDUT;

	    // DUT z position before rotation
	    double _zDUT;

	    // The maximum allowed distance for fitting a hit in X
	    double _distMax_X;

	    // The maximum allowed distance for fitting a hit in X
	    double _distMax_Y;

	    // Additional alignment and rotations from the GEAR file
	    std::vector < float > _DUTalign;

	    // Fiducial cuts for the residual
	    std::vector < float > _fiducialcut;

	    // The function to get the alignment from file
	    void getAlignment ( LCEvent * event );

	    // The function to get the prealignment from file
	    void getPreAlignment ( LCEvent * event );

	    // The function to fill the histogram comparing fitted track position and clusters
	    void fillPrecisionHisto ( float x, float y );

	    // The function to fill the histogram comparing different impact point methods
	    void fillestimationhisto ( float x1, float y1, float x2, float y2, float x3, float y3 );

	    // We fill histograms for the DUT residual
	    void fillTDCResHisto ( float hitx, float hity, float trackx, float tracky, float tdctime );

	    // We fill histograms for the DUT residual
	    void fillresihistos ( float hitx, float hity, float trackx, float tracky, int eventnr, double chi2, double ndf );

	    // A hitmap of the tracks
	    void fillHitmapHisto ( float hitx, float hity, int match );

	    // Residuals after cuts
	    void Filteredfillresihistos ( float hitx, float hity, float trackx, float tracky );

	    // Fill debug histo on calculations
	    void filldebughistos ( float x1, float x2, float y1, float y2 );

	    // Fill eta histo
	    void fillEtaHisto ( float left, float right, float left1, float right1, float left2, float right2, float local, int event, float tdctime );

	    // Flag if alignment is loaded
	    bool _alignmentloaded;

	    // Flag if prealignment is loaded
	    bool _pre_alignmentloaded;

	    // Flag if reference is loaded
	    bool _referenceloaded;

	    // The variables to store the alignment from the files
	    double * _dut_align_x;
	    double * _dut_align_y;
	    double * _dut_align_z;
	    double * _dut_align_a;
	    double * _dut_align_b;
	    double * _dut_align_c;
	    double * _dut_align_x_error;
	    double * _dut_align_y_error;
	    double * _dut_align_z_error;
	    double * _dut_align_a_error;
	    double * _dut_align_b_error;
	    double * _dut_align_c_error;
	    double * _dut_pre_align_x;
	    double * _dut_pre_align_y;
	    double * _dut_pre_align_z;
	    double * _dut_pre_align_a;
	    double * _dut_pre_align_b;
	    double * _dut_pre_align_c;
	    double _dut_original_pre_align_x;
	    double _dut_original_pre_align_y;
	    double _dut_original_pre_align_z;
	    double _dut_original_pre_align_a;
	    double _dut_original_pre_align_b;
	    double _dut_original_pre_align_c;

	    // The pitch and pixel count of the DUT
	    double _pitchx;
	    double _pitchy;
	    int _pixelx;
	    int _pixely;

	    // Counting run headers in a run
	    int _nRun;

	    // Counting events in a run
	    int _nEvt;

	    // Run number
	    int _runNr;

	    // Event number
	    int _evtNr;

	    // Use previous track fit
	    bool _useTrackFit;

	    // Use triplet guess
	    bool _useTriplet;

	    // Matched a DUT hit
	    bool _foundDUTHit;

	    // Track limit
	    int _tracklimit;

	    // The total count of matched hits 
	    int _matchedhits;

	    #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

	    // The Tuple name
	    static std::string _FitTupleName;
	    AIDA::ITuple * _FitTuple;

	    // DUT hitmap
	    void fillLocalHitmapHisto ( float x, float y, float z, int map );

	    void fillInterStripHitmap ( double x, double y );

	    // The reference hit collection
	    std::string _referencecollectionname;

	    // The function to get the reference hit collection
	    void getReference ( );

	    // DUT positions from the reference hit collection
	    double _x_refhit;
	    double _y_refhit;
	    double _z_refhit;
	    double _a_refhit;
	    double _b_refhit;
	    double _c_refhit;

	    int _polarity;

	    #endif

	    void dolandaugausfit ( string tempHistoName );

	} ;

	// A global instance of the processor.
	EUTelFitTupleAlibava aEUTelFitTupleAlibava ;

}

#endif
