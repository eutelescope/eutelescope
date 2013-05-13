#ifdef USE_GBL

// C++
#include <map>
#include <memory>
#include <string>
#include <vector>

// LCIO
#include <EVENT/LCCollection.h>

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

// AIDA
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMath.h"
#endif

// EUTELESCOPE
#include "EUTelProcessorTrackingGBLTrackFit.h"

#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelUtility.h"
#include "EUTelExhaustiveTrackFinder.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
// Cluster types
#include "EUTelSparseCluster2Impl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelFFClusterImpl.h"


using namespace lcio;
using namespace marlin;
using namespace eutelescope;

/**  EUTelProcessorTrackingGBLTrackFit
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates histograms.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Requires collection of track candidates hits.
 *
 *  <h4>Output</h4> 
 *  <li> Histograms.
 *  <li> Collection of fitted tracks.
 */

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelProcessorTrackingGBLTrackFit::_chi2GblFitHistName = "chi2GblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_probGblFitHistName = "probGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_residGblFitHistName = "ResidualsGblFit";
std::string EUTelProcessorTrackingGBLTrackFit::_residGblFitHistNameX = "ResidualsGblFit_x";
std::string EUTelProcessorTrackingGBLTrackFit::_residGblFitHistNameY = "ResidualsGblFit_y";
std::string EUTelProcessorTrackingGBLTrackFit::_resid2DGblFitHistNameXvsX = "Residuals2DGblFit_xVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_resid2DGblFitHistNameYvsX = "Residuals2DGblFit_yVSx";
std::string EUTelProcessorTrackingGBLTrackFit::_resid2DGblFitHistNameXvsY = "Residuals2DGblFit_xVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_resid2DGblFitHistNameYvsY = "Residuals2DGblFit_yVSy";
std::string EUTelProcessorTrackingGBLTrackFit::_kinkGblFitHistNameX = "KinksGblFit_x";
std::string EUTelProcessorTrackingGBLTrackFit::_kinkGblFitHistNameY = "KinksGblFit_y";
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

/** Default constructor */
EUTelProcessorTrackingGBLTrackFit::EUTelProcessorTrackingGBLTrackFit( ) :
 Processor( "EUTelProcessorTrackingGBLTrackFit" ),
 _trackCandidateHitsInputCollectionName( "TrackCandidateHitCollection" ),
 _tracksOutputCollectionName( "TrackCollection" ),
 _nProcessedRuns( 0 ),
 _nProcessedEvents( 0 ) {

    // Processor description
    _description = "EUTelProcessorTrackingGBLTrackFit performs track fits using GBL optionally writing data files for MILLEPEDE II.";

    // TrackerHit input collection
    registerInputCollection( LCIO::LCGENERICOBJECT,
            "TrackCandHitOutputCollectionName",
            "Input track candidates hits collection name",
            _trackCandidateHitsInputCollectionName,
            std::string( "TrackCandidateHitCollection" ) );

    // Track output collection
    registerOutputCollection( LCIO::TRACK,
            "TracksOutputCollectionName",
            "Output tracks collection name",
            _tracksOutputCollectionName,
            std::string( "TrackCollection" ) );



    // Optional processor parameters that define finder settings

    registerOptionalParameter( "BinaryFilename", "Name of the Millepede binary file.", _binaryFilename, std::string( "mille.bin" ) );
    
    registerOptionalParameter( "GeometryFilename", "Name of the TGeo geometry definition file.", _tgeoFileName , std::string( "TELESCOPE.root" ) );

}

void EUTelProcessorTrackingGBLTrackFit::init( ) {

    streamlog_out( DEBUG ) << "EUTelProcessorTrackingGBLTrackFit::init( )" << std::endl;

    // usually a good idea to
    printParameters( );

    // Reset counters
    _nProcessedRuns = 0;
    _nProcessedEvents = 0;


    // Getting access to geometry description
    _geometry = new EUTelGeometryTelescopeGeoDescription;
    _geometry -> initializeTGeoDescription( _tgeoFileName );

    // Instantiate track fitter. This is a working horse of the processor.
    {
        streamlog_out( DEBUG ) << "Initialisation of track fitter" << std::endl;

        EUTelGBLFitter* Fitter = new EUTelGBLFitter( "myGBLFitter" );
        Fitter->SetMilleBinary( _milleGBL );
        _trackFitter = Fitter;

        if ( _trackFitter == 0 ) {
            streamlog_out( ERROR ) << "Can't allocate an instance of EUTelGBLFitter. Stopping ..." << std::endl;
            throw UnknownDataTypeException( "Track finder was not created" );
        }
    }

    // Instantiate millipede output. 
    {
        streamlog_out( DEBUG ) << "Initialising Mille..." << std::endl;

        const unsigned int reserveSize = 80000;
        _milleGBL = new gbl::MilleBinary( _binaryFilename, reserveSize );

        if ( _milleGBL == 0 ) {
            streamlog_out( ERROR ) << "Can't allocate an instance of mMilleBinary. Stopping ..." << std::endl;
            throw UnknownDataTypeException( "MilleBinary was not created" );
        }
    }
    // Book histograms
    bookHistograms( );
}

void EUTelProcessorTrackingGBLTrackFit::processRunHeader( LCRunHeader* run ) {

    auto_ptr<EUTelRunHeaderImpl> header( new EUTelRunHeaderImpl( run ) );
    header->addProcessor( type( ) );


    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, warn the user.

    if ( header->getGeoID( ) == 0 )
        streamlog_out( WARNING0 ) << "The geometry ID in the run header is set to zero." << endl
            << "This may mean that the GeoID parameter was not set" << endl;


    if ( header->getGeoID( ) != _geometry->_siPlanesParameters->getSiPlanesID( ) ) {
        streamlog_out( WARNING5 ) << "Error during the geometry consistency check: " << endl
                << "The run header says the GeoID is " << header->getGeoID( ) << endl
                << "The GEAR description says is     " << _geometry->_siPlanesParameters->getSiPlanesID( ) << endl;
    }

    _nProcessedRuns++;
}

void EUTelProcessorTrackingGBLTrackFit::processEvent( LCEvent * evt ) {

    if ( isFirstEvent( ) ) {
        ;
    }

    EUTelEventImpl * event = static_cast < EUTelEventImpl* > ( evt );

    // Do not process last events
    if ( event->getEventType( ) == kEORE ) {
        streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
        return;
    } else if ( event->getEventType( ) == kUNKNOWN ) {
        streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber( ) << " in run " << event->getRunNumber( )
                << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }

    // Try to access collection

    LCCollection* col = NULL;
    try {
        col = evt->getCollection( _trackCandidateHitsInputCollectionName );
    } catch ( DataNotAvailableException e ) {
        streamlog_out( WARNING ) << _trackCandidateHitsInputCollectionName << " collection not available" << std::endl;
        throw marlin::SkipEventException( this );
    }

    // this will only be entered if the collection is available
    if ( col != NULL ) {
        if ( _nProcessedEvents % 10000 == 1 ) streamlog_out( MESSAGE5 ) << "EUTelProcessorTrackingGBLTrackFit" << endl;

        //        TrackerHitImpl * hit = static_cast < TrackerHitImpl* > ( collection->getElementAt( iHit ) );

        // Perform fit for all found track candidates
        // ------------------------------------------
        unsigned int numData;
        TVectorD residual( 200 );
        TVectorD measErr( 200 );
        TVectorD residualErr( 200 );
        TVectorD downWeight( 200 );
        //        if (_nTracks != 0 && _nTracks == 1) {
        //            _trackFitter->SetTrackCandidates(trackCandidates);
        //            _trackFitter->FitTracks();
        ////
        //            double chi2Trk = 0.;
        //            int ndfTrk = 0;
        //
        //            IMPL::LCCollectionVec* fittrackvec;
        //            fittrackvec = static_cast<EUTelGBLFitter*> (_trackFitter)->GetFitTrackVec();
        //            IMPL::LCCollectionVec::const_iterator itFitTrack;
        //
        //            int iCounter = 0;
        //            for (itFitTrack = fittrackvec->begin(); itFitTrack != fittrackvec->end(); ++itFitTrack) {
        //                chi2Trk = static_cast<TrackImpl*> (*itFitTrack)->getChi2();
        //                ndfTrk = static_cast<TrackImpl*> (*itFitTrack)->getNdf();
        //
        //                static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _chi2GblFitHistName ]) -> fill(chi2Trk);
        //                if( chi2Trk < 17 ) static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _probGblFitHistName ]) -> fill(TMath::Prob(chi2Trk, ndfTrk));
        //
        //
        //                std::map< int, gbl::GblTrajectory* > gblTracks = static_cast<EUTelGBLFitter*> (_trackFitter)->GetGblTrackCandidates();
        //
        //                std::stringstream sstr;
        //                gbl::GblTrajectory* gblTraj = gblTracks[ iCounter ];
        //		//gblTraj->printTrajectory();
        //		//gblTraj->printPoints();
        //		//gblTraj->printData();
        //                std::vector< gbl::GblPoint > gblPointVec = static_cast<EUTelGBLFitter*> (_trackFitter)->GetGblTracksPoints()[iCounter];
        //                std::vector< gbl::GblPoint >::const_iterator itGblPoint = gblPointVec.begin();
        //                int iPlane = 0; // wrong in case of missing planes
        //                for (; itGblPoint != gblPointVec.end(); ++itGblPoint) {
        //                    if (iPlane > 5) continue;
        //		    //if ( itGblPoint->getLabel() < 1000 )
        //		    if ( itGblPoint->getLabel() % 3 == 1 ) 
        //		    {
        //			    streamlog_out(DEBUG0) << std::setw(15) << itGblPoint->getLabel() << std::endl;
        //                            // spatial residuals
        //			    gblTraj->getMeasResults(itGblPoint->getLabel(), numData, residual, measErr, residualErr, downWeight);
        //			    sstr << _residGblFitHistNameX << iPlane;
        //			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[0] << std::setw(15) << std::setprecision(5) << residualErr[0] << std::endl;
        //			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[0] / residualErr[0], downWeight[0]);
        //			    sstr.str(std::string());
        //			    sstr << _residGblFitHistNameY << iPlane;
        //			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[1] << std::setw(15) << std::setprecision(5) << residualErr[1] << std::endl;
        //			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[1] / residualErr[1], downWeight[1]);
        //			    sstr.str(std::string());
        //                            // kinks
        //                            gblTraj->getScatResults (itGblPoint->getLabel(), numData, residual, measErr, residualErr, downWeight);
        //                            sstr << _kinkGblFitHistNameX << iPlane;
        //			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[0] << std::setw(15) << std::setprecision(5) << residualErr[0] << std::endl;
        //			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[0], downWeight[0]);
        //			    sstr.str(std::string());
        //			    sstr << _kinkGblFitHistNameY << iPlane;
        //			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[1] << std::setw(15) << std::setprecision(5) << residualErr[1] << std::endl;
        //			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[1], downWeight[1]);
        //			    sstr.str(std::string());
        //                            
        //                            // 2D histograms
        //                            const double* hitpos = trackCandidates.begin()->second[iPlane]->getPosition();              // wrong in case of empty planes
        //                            sstr << _resid2DGblFitHistNameXvsX << iPlane;
        //                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[0], hitpos[0], downWeight[0]);
        //                            sstr.str(std::string());
        //                            sstr << _resid2DGblFitHistNameXvsY << iPlane;
        //                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[0], hitpos[1], downWeight[0]);
        //                            sstr.str(std::string());
        //                            sstr << _resid2DGblFitHistNameYvsX << iPlane;
        //                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[1], hitpos[0], downWeight[1]);
        //                            sstr.str(std::string());
        //                            sstr << _resid2DGblFitHistNameYvsY << iPlane;
        //                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[1], hitpos[1], downWeight[1]);
        //                            sstr.str(std::string());
        //                            
        //			    if ( itGblPoint->getLabel() < 1000 )++iPlane;
        //		    }
        //                }
        //
        //                IMPL::LCCollectionVec::const_iterator itFitTrack;
        //                iCounter++;
        //            }
        //        }
    }

    _nProcessedEvents++;
    
    if ( isFirstEvent( ) ) _isFirstEvent = false;
}

void EUTelProcessorTrackingGBLTrackFit::check( LCEvent * evt ) {
    // nothing to check here
}

void EUTelProcessorTrackingGBLTrackFit::end( ) {

    delete _geometry;
    delete _milleGBL;

    streamlog_out( MESSAGE ) << "EUTelProcessorTrackingGBLTrackFit::end()  " << name( )
            << " processed " << _nProcessedEvents << " events in " << _nProcessedRuns << " runs "
            << std::endl;

}

void EUTelProcessorTrackingGBLTrackFit::bookHistograms( ) {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    try {
        streamlog_out( DEBUG ) << "Booking histograms..." << std::endl;

        const int chi2NBin = 1000;
        const double chi2Min = 0.;
        const double chi2Max = 1000.;

        const int probNBin = 1000;
        const double probMin = 0.;
        const double probMax = 1.;

        int NBinX = 4000;
        double MinX = -20.;
        double MaxX = 20.;
        int NBinY = 4000;
        double MinY = -20.;
        double MaxY = 20.;


        // GBL fits
        AIDA::IHistogram1D * chi2GblFit =
                marlin::AIDAProcessor::histogramFactory( this )->createHistogram1D( _chi2GblFitHistName, chi2NBin, chi2Min, chi2Max );
        if ( chi2GblFit ) {
            chi2GblFit->setTitle( "#chi^{2} of track candidates; #chi^{2};N Tracks" );
            _aidaHistoMap1D.insert( std::make_pair( _chi2GblFitHistName, chi2GblFit ) );
        } else {
            streamlog_out( ERROR2 ) << "Problem booking the " << ( _chi2GblFitHistName ) << std::endl;
            streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }

        AIDA::IHistogram1D * probGblFit =
                marlin::AIDAProcessor::histogramFactory( this )->createHistogram1D( _probGblFitHistName, probNBin, probMin, probMax );
        if ( probGblFit ) {
            probGblFit->setTitle( "Probability of track fit; Prob;N Tracks" );
            _aidaHistoMap1D.insert( std::make_pair( _probGblFitHistName, probGblFit ) );
        } else {
            streamlog_out( ERROR2 ) << "Problem booking the " << ( _probGblFitHistName ) << std::endl;
            streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
        }

        // Residuals after fit
        std::stringstream sstm;
        std::string residGblFitHistName;
        std::string histTitle;
        for ( int iPlane = 0; iPlane < 6; iPlane++ ) {
            sstm << _residGblFitHistNameX << iPlane;
            residGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Normalised residuals. Plane " << iPlane << "X direction; r; N hits";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram1D( residGblFitHistName, NBinX, MinX, MaxX );
            if ( residGblFit ) {
                residGblFit->setTitle( histTitle );
                _aidaHistoMap1D.insert( std::make_pair( residGblFitHistName, residGblFit ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( residGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );
        }

        for ( int iPlane = 0; iPlane < 6; iPlane++ ) {
            sstm << _residGblFitHistNameY << iPlane;
            residGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Normalised residuals. Plane " << iPlane << "Y direction; r; N hits";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram1D( residGblFitHistName, NBinX, MinX, MaxX );
            if ( residGblFit ) {
                residGblFit->setTitle( histTitle );
                _aidaHistoMap1D.insert( std::make_pair( residGblFitHistName, residGblFit ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( residGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );
        }

        // 2D histograms
        NBinX = 100;
        NBinY = 100;
        MinX = -0.001;
        MaxX = 0.001;
        std::string resid2DGblFitHistName;
        for ( int iPlane = 0; iPlane < 6; iPlane++ ) {
            sstm << _resid2DGblFitHistNameXvsX << iPlane;
            resid2DGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Normalised residuals. Plane " << iPlane << "; x (mm); rx";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram2D( resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY );
            if ( residGblFit1 ) {
                residGblFit1->setTitle( histTitle );
                _aidaHistoMap2D.insert( std::make_pair( resid2DGblFitHistName, residGblFit1 ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( resid2DGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );

            sstm << _resid2DGblFitHistNameXvsY << iPlane;
            resid2DGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Normalised residuals. Plane " << iPlane << "; y (mm); rx";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram2D( resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY );
            if ( residGblFit2 ) {
                residGblFit2->setTitle( histTitle );
                _aidaHistoMap2D.insert( std::make_pair( resid2DGblFitHistName, residGblFit2 ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( resid2DGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );
        }

        for ( int iPlane = 0; iPlane < 6; iPlane++ ) {
            sstm << _resid2DGblFitHistNameYvsX << iPlane;
            resid2DGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Normalised residuals. Plane " << iPlane << "; x (mm); ry";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram2D( resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY );
            if ( residGblFit1 ) {
                residGblFit1->setTitle( histTitle );
                _aidaHistoMap2D.insert( std::make_pair( resid2DGblFitHistName, residGblFit1 ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( resid2DGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );

            sstm << _resid2DGblFitHistNameYvsY << iPlane;
            resid2DGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Normalised residuals. Plane " << iPlane << "; y (mm); ry";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram2D( resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY );
            if ( residGblFit2 ) {
                residGblFit2->setTitle( histTitle );
                _aidaHistoMap2D.insert( std::make_pair( resid2DGblFitHistName, residGblFit2 ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( resid2DGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );
        }

        // Kink angles after fit
        MinX = -0.001;
        MaxX = 0.001;
        std::string kinkGblFitHistName;
        for ( int iPlane = 0; iPlane < 6; iPlane++ ) {
            sstm << _kinkGblFitHistNameX << iPlane;
            kinkGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Kink angles. Plane " << iPlane << "X direction; kink (rad); N hits";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram1D( kinkGblFitHistName, NBinX, MinX, MaxX );
            if ( residGblFit ) {
                residGblFit->setTitle( histTitle );
                _aidaHistoMap1D.insert( std::make_pair( kinkGblFitHistName, residGblFit ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( kinkGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );
        }

        for ( int iPlane = 0; iPlane < _geometry->_nPlanes; iPlane++ ) {
            sstm << _kinkGblFitHistNameY << iPlane;
            kinkGblFitHistName = sstm.str( );
            sstm.str( std::string( ) );
            sstm << "Kink angles. Plane " << iPlane << "Y direction; kink (rad); N hits";
            histTitle = sstm.str( );
            sstm.str( std::string( "" ) );
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory( this )->createHistogram1D( kinkGblFitHistName, NBinX, MinX, MaxX );
            if ( residGblFit ) {
                residGblFit->setTitle( histTitle );
                _aidaHistoMap1D.insert( std::make_pair( kinkGblFitHistName, residGblFit ) );
            } else {
                streamlog_out( ERROR2 ) << "Problem booking the " << ( kinkGblFitHistName ) << std::endl;
                streamlog_out( ERROR2 ) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            }
            sstm.str( std::string( "" ) );
        }
    } catch ( lcio::Exception& e ) {
        streamlog_out( WARNING2 ) << "Can't allocate histgrams. Continue without histogramming" << endl;
    }
#endif // defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
}

#endif // USE_GBL