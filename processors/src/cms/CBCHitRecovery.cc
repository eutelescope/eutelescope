/*
 * Created by Mykyta Haranko
 *  (2018 DESY)
 *
 *  email:mykyta.haranko@desy.de
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
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
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
#include <cmath>
#include <algorithm>

// eutelescope includes ""
#include "EUTELESCOPE.h"
#include "anyoption.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"

// processor include
#include "CBCHitRecovery.h"

#define PI 3.14159265

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace IMPL;
using namespace eutelescope;


CBCHitRecovery::CBCHitRecovery ( ) : Processor ( "CBCHitRecovery" ),
_aidaHistoMap ( )
{

    _description = "CBCHitRecovery recovers the real hit positions from the virtual dut and track.";

    registerProcessorParameter ( "InputTrackCollectionName", "The name of the Tracks collection we want to read", _InputTrackCollectionName, string ( "tracks_input" ) );

    registerProcessorParameter ( "InputFitHitsCollectionName", "The name of the fithits collection we want to read", _InputFitHitsCollectionName, string ( "fithits_input" ) );
    
    registerProcessorParameter ( "CBCInputCollectionName", "The name of the CBC collection we want to read", _cbcInputCollectionName, string ( "cbc_input" ) );

    registerProcessorParameter ( "CBCDataOutputCollectionName", "The name of the CBC data collection we want to write", _cbcDataOutputCollectionName, string ( "cbc_data_output" ) );

    registerProcessorParameter ( "CBCVirtualDUTId", "The id of the virtual CBC DUT", _cbcVirtualDUTId, 0 );
    
    registerProcessorParameter ( "CBCRealDUTsVec", "The ids of the real DUTs (has to be two)", _cbcRealDUTsVec, std::vector < int > () );
}


void CBCHitRecovery::init ( )
{
        streamlog_out ( MESSAGE4 ) << "Running init" << endl;
        printParameters ( );

        _siPlanesParameters  = const_cast < gear::SiPlanesParameters* > ( & ( Global::GEAR -> getSiPlanesParameters ( ) ) );
        _siPlanesLayerLayout = const_cast < gear::SiPlanesLayerLayout* > ( & ( _siPlanesParameters -> getSiPlanesLayerLayout ( ) ) );

        for ( int i = 0; i < _siPlanesParameters -> getSiPlanesNumber ( ); i++ )
        {
                int cPlaneID = _siPlanesLayerLayout -> getID (i);
                if ( std::find(_cbcRealDUTsVec.begin(), _cbcRealDUTsVec.end(), cPlaneID) != _cbcRealDUTsVec.end())
                {
                        // get the angles
                        double alpha = _siPlanesLayerLayout -> getLayerRotationZY ( i );
                        double beta  = _siPlanesLayerLayout -> getLayerRotationZX ( i );

                        // we need to calculate the normal vectors 
                        double *normalvector = new double[3];
                        normalvector[0] = sin ( beta * PI / 180.0 );
                        normalvector[1] = -sin ( alpha * PI / 180.0 );
                        normalvector[2] = cos ( alpha * PI / 180.0 ) * cos ( beta * PI / 180.0 );

                        // put it to the map now
                        _dutNormalMap.insert ( make_pair ( cPlaneID, normalvector ) );

                        // we need to calculate the normal vectors 
                        double *posvector = new double[3];
                        posvector[0] = _siPlanesLayerLayout -> getLayerPositionX( i );
                        posvector[1] = _siPlanesLayerLayout -> getLayerPositionY( i );
                        posvector[2] = _siPlanesLayerLayout -> getLayerPositionZ( i ) + 0.5 *_siPlanesLayerLayout -> getSensitiveThickness( i );

                        // put it to the map now
                        _dutPosMap.insert ( make_pair ( cPlaneID, posvector ) );

                        // test if it was saved properly
                        double *pos_test = _dutPosMap.at(cPlaneID);
                        streamlog_out ( MESSAGE2 ) << "The position vector of the DUT" << cPlaneID << " is: X = " << pos_test[0] << ", Y = " << pos_test[1] << ", Z = " << pos_test[2] << endl;
                        double *normal_test = _dutNormalMap.at(cPlaneID);
                        streamlog_out ( MESSAGE2 ) << "The normal vector of the DUT" << cPlaneID << " is: X = " << normal_test[0] << ", Y = " << normal_test[1] << ", Z = " << normal_test[2] << endl;
                }               
        }
}


void CBCHitRecovery::processRunHeader ( LCRunHeader * rdr )
{
    streamlog_out ( MESSAGE4 ) << "Running processRunHeader" << endl;

    auto arunHeader = std::make_unique < EUTelRunHeaderImpl > ( rdr );
    arunHeader -> addProcessor ( type ( ) );

    bookHistos ( );

}


void CBCHitRecovery::processEvent ( LCEvent * anEvent )
{

    if ( anEvent -> getEventNumber ( ) % 1000 == 0 )
    {
        streamlog_out ( MESSAGE4 ) << "Looping events " << anEvent -> getEventNumber ( ) << endl;
    }

        // the collection we read
        //LCCollectionVec * inputTrackVec;
        LCCollectionVec * inputFitPointVec;
        LCCollectionVec * inputHitsVec;
        
        // the collection we output
        LCCollectionVec * outputCollectionVec;

        try
        {
            outputCollectionVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcDataOutputCollectionName ) );
        }
        catch ( lcio::DataNotAvailableException& e )
        {
            outputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);
        }

        // find process the tracks
        try
        {
                // give the collection vec its data
                //inputTrackVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _InputTrackCollectionName ) );
                inputFitPointVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _InputFitHitsCollectionName ) );
                inputHitsVec = dynamic_cast < LCCollectionVec * > ( anEvent -> getCollection ( _cbcInputCollectionName ) );
                
                // getting the DUT normal vector in the first event
                if( anEvent -> getEventNumber ( ) == 0 ) {
                        // get dut normal
                        LCCollectionVec* dutnormalvec_collection = dynamic_cast < LCCollectionVec * > (anEvent->getCollection("dutnormal"));
                        // get object
                        LCGenericObjectImpl* dutnormalvec = dynamic_cast < LCGenericObjectImpl * > (dutnormalvec_collection->at(0));

                        // allocare
                        _VirtualDutNormal_zplus = new double[3];
                        _VirtualDutNormal_zminus = new double[3];
                        _VirtualDutPos = new double[3];

                        // dut normal vec                       
                        for(int i = 0; i < 3; i++) {
                                _VirtualDutNormal_zplus[i] = dutnormalvec->getDoubleVal(i+3);   
                                _VirtualDutNormal_zminus[i] = dutnormalvec->getDoubleVal(i+3);  
                                _VirtualDutPos[i] = dutnormalvec->getDoubleVal(i)/1000.0;       
                        }
                        // we need to get the opposite normal vec for pointing upstream
                        _VirtualDutNormal_zminus[2] = -1*_VirtualDutNormal_zminus[2];

                        streamlog_out ( MESSAGE4 ) << "Loaded the DUT Pos Vector: X = " << _VirtualDutPos[0] << ", Y = " << _VirtualDutPos[1] << ", Z = " << _VirtualDutPos[2] << endl;
                        streamlog_out ( MESSAGE4 ) << "Loaded the DUT Normal Vector: X = " << _VirtualDutNormal_zplus[0] << ", Y = " << _VirtualDutNormal_zplus[1] << ", Z = " << _VirtualDutNormal_zplus[2] << endl;
                
                        //delete dutnormalvec;
                        //delete dutnormalvec_collection;
                }

                
                /*// loop over the tracks in the event
                int nTracks = inputTrackVec -> getNumberOfElements ( );
                for ( int i = 0; i < nTracks; ++i ) 
                {
                    TrackImpl * track = dynamic_cast < TrackImpl * > ( inputTrackVec -> getElementAt ( i ) );

                    //delete track;
                }*/

                // variables
                const double *cTelescope2Pos = nullptr;
                const double *cVirtualHitPos = nullptr;
                const double *cTelescope3Pos = nullptr;

                // find the fit hit points of the track (needed to calculate the fir hits in the dut
                int nEntries = inputFitPointVec -> getNumberOfElements ( );
                for ( int i = 0; i < nEntries; ++i ) 
                {
                    TrackerHitImpl * hit = dynamic_cast < TrackerHitImpl * > ( inputFitPointVec -> getElementAt ( i ) );

                    // sensor id decoder
                    CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( inputFitPointVec );
                    int sensorID = inputCellIDDecoder ( hit ) ["sensorID"];
                    // check that we are on the virtual cbc hit
                    if (sensorID == _cbcVirtualDUTId) {
                        cVirtualHitPos = hit->getPosition();
                        
                        dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualHitPosX"] ) -> fill (cVirtualHitPos[0]);
                        dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualHitPosY"] ) -> fill (cVirtualHitPos[1]);
                        dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualHitPosZ"] ) -> fill (cVirtualHitPos[2]);
                    } else if (sensorID == 2) {
                        cTelescope2Pos = hit->getPosition();
                    } else if (sensorID == 3) {
                        cTelescope3Pos = hit->getPosition();
                    }

                    //delete hit;
                }

                // check that there are tracks in the event
                if (nEntries != 0) {
                        // and that we have succeeded to get the track positions
                        if (cTelescope2Pos == nullptr || cVirtualHitPos == nullptr || cTelescope3Pos == nullptr) {
                                streamlog_out ( WARNING ) << "No fit hit for some of the planes in the event " <<  anEvent -> getEventNumber ( ) << endl;
                        } else {
                                // track from upstream (it pointed from dut to the telescope plane 2, then we will reconstruct sensor 60 hit)
                                double cTrackVecUpstream[3];
                                for(int i = 0; i < 3; i++) cTrackVecUpstream[i] = cTelescope2Pos[i] - cVirtualHitPos[i];

                                // track to downstream
                                double cTrackVecDownstream[3];
                                for(int i = 0; i < 3; i++) cTrackVecDownstream[i] = cTelescope3Pos[i] - cVirtualHitPos[i];

                                // print more for debugging
                                /*if (anEvent->getEventNumber() < 50) {
                                        streamlog_out (DEBUG1) << "vec up: " << cTrackVecUpstream[0] << ", " << cTrackVecUpstream[1] << ", " << cTrackVecUpstream[2] << endl;
                                        streamlog_out (DEBUG1) << "vec down: " << cTrackVecDownstream[0] << ", " << cTrackVecDownstream[1] << ", " << cTrackVecDownstream[2] << endl;
                                }*/

                                // cos phi between the vectors
                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["cosPhi"] ) -> fill (cos_alpha(cTrackVecUpstream,cTrackVecDownstream));

                                // fit hit pos in dut 0
                                double *fit_hit_pos_dut0 = this->hit_pos(cTrackVecUpstream, _VirtualDutNormal_zminus, cVirtualHitPos, -2.0);
                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosX_DUT0"] ) -> fill (fit_hit_pos_dut0[0]);
                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosY_DUT0"] ) -> fill (fit_hit_pos_dut0[1]);
                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosZ_DUT0"] ) -> fill (fit_hit_pos_dut0[2]);


                                // fit hit pos in dut 1
                                double *fit_hit_pos_dut1 = this->hit_pos(cTrackVecDownstream, _VirtualDutNormal_zplus, cVirtualHitPos, 2.0);
                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosX_DUT1"] ) -> fill (fit_hit_pos_dut1[0]);
                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosY_DUT1"] ) -> fill (fit_hit_pos_dut1[1]);
                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["FitHitPosZ_DUT1"] ) -> fill (fit_hit_pos_dut1[2]);

                                // now calculate the residuals
                                int nHits = inputHitsVec -> getNumberOfElements ( );
                                for ( int iHit = 0; iHit < nHits; ++iHit ) {
                                        TrackerHitImpl * hit = dynamic_cast < TrackerHitImpl * > ( inputHitsVec -> getElementAt ( iHit ) );

                                        // sensor id decoder
                                        CellIDDecoder < TrackerHitImpl > inputCellIDDecoder ( inputFitPointVec );
                                        int sensorID = inputCellIDDecoder ( hit ) ["sensorID"];

                                        // virtual
                                        if (sensorID == _cbcVirtualDUTId) {
                                                const double *pos = hit->getPosition();
                                                double resX = (cVirtualHitPos[0] - pos[0])*1000;
                                                //double resY = (cVirtualHitPos[1] - pos[1])*1000;
                                                //double resZ = (cVirtualHitPos[2] - pos[2])*1000;

                                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["VirtualResidualX"] ) -> fill (resX);
                                        }

                                        // first sensor
                                        if (sensorID == _cbcRealDUTsVec.at(0)) {
                                                const double *pos = hit->getPosition();
                                                double resX = (fit_hit_pos_dut0[0] - pos[0])*1000;
                                                //double resY = (fit_hit_pos_dut0[1] - pos[1])*1000;
                                                //double resZ = (fit_hit_pos_dut0[2] - pos[2])*1000;

                                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ResidualX_DUT0"] ) -> fill (resX);
                                        }

                                        // second sensor
                                        if (sensorID == _cbcRealDUTsVec.at(1)) {
                                                const double *pos = hit->getPosition();
                                                double resX = (fit_hit_pos_dut1[0] - pos[0])*1000;
                                                //double resY = (fit_hit_pos_dut1[1] - pos[1])*1000;
                                                //double resZ = (fit_hit_pos_dut1[2] - pos[2])*1000;

                                                dynamic_cast < AIDA::IHistogram1D* > ( _aidaHistoMap["ResidualX_DUT1"] ) -> fill (resX);
                                        }
                                }

                                // delete the hist positions
                                delete fit_hit_pos_dut0;
                                delete fit_hit_pos_dut1;

                                // clear pos vectors
                                //delete cTelescope2Pos;
                                //delete cVirtualHitPos;
                                //delete cTelescope3Pos;
                        }
                }

        }
        catch ( lcio::DataNotAvailableException& )
        {

        }

        anEvent->addCollection( outputCollectionVec, _cbcDataOutputCollectionName );

}


void CBCHitRecovery::check ( LCEvent * /* evt */ )
{

}


void CBCHitRecovery::end ( )
{
    streamlog_out ( MESSAGE4 ) << "Successfully finished!" << endl;

}


void CBCHitRecovery::fillHistos ( )
{

}


void CBCHitRecovery::bookHistos ( )
{
        string basePath = "Input";
        AIDAProcessor::tree ( this ) -> mkdir ( basePath.c_str ( ) );
        basePath.append ( "/" );

        AIDA::IHistogram1D * cVirtualHitXHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "VirtualHitPosX" ).c_str ( ), 100, -100, 100 );
        _aidaHistoMap.insert ( make_pair ( "VirtualHitPosX", cVirtualHitXHist ) );
        cVirtualHitXHist -> setTitle ( "Virtual Hit Position X;X [mm];Entries" );

        AIDA::IHistogram1D * cVirtualHitYHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "VirtualHitPosY" ).c_str ( ), 100, -100, 100 );
        _aidaHistoMap.insert ( make_pair ( "VirtualHitPosY", cVirtualHitYHist ) );
        cVirtualHitYHist -> setTitle ( "Virtual Hit Position Y;Y [mm];Entries" );

        AIDA::IHistogram1D * cVirtualHitZHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "VirtualHitPosZ" ).c_str ( ), 100, 365, 375 );
        _aidaHistoMap.insert ( make_pair ( "VirtualHitPosZ", cVirtualHitZHist ) );
        cVirtualHitZHist -> setTitle ( "Virtual Hit Position Z;Z [mm];Entries" );

        AIDA::IHistogram1D * cCosPhiHist = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePath + "cosPhi" ).c_str ( ), 500, -1.000001, -0.99999 );
        _aidaHistoMap.insert ( make_pair ( "cosPhi", cCosPhiHist ) );
        cCosPhiHist -> setTitle ( "cos between inbound and outbound track vectors;cos(phi);Entries" );

        string basePathOut = "Output";
        AIDAProcessor::tree ( this ) -> mkdir ( basePathOut.c_str ( ) );
        basePathOut.append ( "/" );

        AIDA::IHistogram1D * cDUT0FitHitPosX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosX_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 100, -100, 100 );
        _aidaHistoMap.insert ( make_pair ( "FitHitPosX_DUT0", cDUT0FitHitPosX ) );
        cDUT0FitHitPosX -> setTitle ( "Fit Hit Pos X - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";X [mm];Entries" );

        AIDA::IHistogram1D * cDUT0FitHitPosY = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosY_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 100, -100, 100 );
        _aidaHistoMap.insert ( make_pair ( "FitHitPosY_DUT0", cDUT0FitHitPosY ) );
        cDUT0FitHitPosY -> setTitle ( "Fit Hit Pos Y - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";Y [mm];Entries" );

        AIDA::IHistogram1D * cDUT0FitHitPosZ = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosZ_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 100, 365, 375 );
        _aidaHistoMap.insert ( make_pair ( "FitHitPosZ_DUT0", cDUT0FitHitPosZ ) );
        cDUT0FitHitPosX -> setTitle ( "Fit Hit Pos Z - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";Z [mm];Entries" );

        AIDA::IHistogram1D * cDUT1FitHitPosX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosX_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 100, -100, 100 );
        _aidaHistoMap.insert ( make_pair ( "FitHitPosX_DUT1", cDUT1FitHitPosX ) );
        cDUT1FitHitPosX -> setTitle ( "Fit Hit Pos X - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";X [mm];Entries" );

        AIDA::IHistogram1D * cDUT1FitHitPosY = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosY_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 100, -100, 100 );
        _aidaHistoMap.insert ( make_pair ( "FitHitPosY_DUT1", cDUT1FitHitPosY ) );
        cDUT1FitHitPosY -> setTitle ( "Fit Hit Pos Y - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";Y [mm];Entries" );

        AIDA::IHistogram1D * cDUT1FitHitPosZ = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "FitHitPosZ_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 100, 365, 375 );
        _aidaHistoMap.insert ( make_pair ( "FitHitPosZ_DUT1", cDUT1FitHitPosZ ) );
        cDUT1FitHitPosX -> setTitle ( "Fit Hit Pos Z - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";Z [mm];Entries" );

        AIDA::IHistogram1D * cVirtualResidualX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "VirtualResidualX" ).c_str ( ), 100, -500, 500 );
        _aidaHistoMap.insert ( make_pair ( "VirtualResidualX", cVirtualResidualX ) );
        cVirtualResidualX -> setTitle ( "Residual X - Virtual DUT;X [um];Entries" );

        AIDA::IHistogram1D * cDUT0ResidualX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "ResidualX_DUT" + to_string(_cbcRealDUTsVec.at(0))).c_str ( ), 600, -1500, 1500 );
        _aidaHistoMap.insert ( make_pair ( "ResidualX_DUT0", cDUT0ResidualX ) );
        cDUT0ResidualX -> setTitle ( "Residual X - DUT " + to_string(_cbcRealDUTsVec.at(0)) + ";X [um];Entries" );

        AIDA::IHistogram1D * cDUT1ResidualX = AIDAProcessor::histogramFactory ( this ) -> createHistogram1D ( ( basePathOut + "ResidualX_DUT" + to_string(_cbcRealDUTsVec.at(1))).c_str ( ), 600, -1500, 1500 );
        _aidaHistoMap.insert ( make_pair ( "ResidualX_DUT1", cDUT1ResidualX ) );
        cDUT1ResidualX -> setTitle ( "Residual X - DUT " + to_string(_cbcRealDUTsVec.at(1)) + ";X [um];Entries" );

}

// get the module of vector
double CBCHitRecovery::vector_get_length(const double *vec) {
        double length = 0;
        for(int i = 0; i < 3; i++) length += vec[i]*vec[i];
        return sqrt(length);
}

// set certain length of the vector
void CBCHitRecovery::vector_set_length(double *& vec, double length) {
        // we also need to normalize it
        double current_length = vector_get_length(vec);
        // now set the new length
        for(int i = 0; i < 3; i++) vec[i] = vec[i]*length/current_length;
}

// angle between two vector
double CBCHitRecovery::cos_alpha(const double *vec1, const double *vec2) {
        // variable 
        double angle = 0;
        // calculate the scalar
        for(int i = 0; i < 3; i++) angle += vec1[i]*vec2[i];
        // now divide by modules
        angle = angle/(vector_get_length(vec1)*vector_get_length(vec2));
        // return
        return angle;   
}

// finds the position of fithit in real plane
double* CBCHitRecovery::hit_pos(const double track[],const double normal[],const double virtual_hit[], double distance) {
        // init vector
        double *pos = new double[3];
        // module of distance
        double abs_distance = distance;
        if(abs_distance < 0) abs_distance = -abs_distance;
        // create copies of the vectors to modify them
        double *track1 = new double[3];
        for(int i = 0; i < 3; i++) {
                track1[i] = track[i];
        }
        // set the length of the track vector to distance/cos(alpha), alpha - angle between normal and track
        vector_set_length(track1, abs_distance/cos_alpha(track, normal));
        // add two to get the position
        for(int i = 0; i < 3; i++) pos[i] = track1[i] + virtual_hit[i];
        // clear
        delete track1;
        // return
        return pos;
}
