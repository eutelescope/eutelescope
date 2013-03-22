//! Contact: Denys Lontkovskyi, DESY <mailto:denys.lontkovskyi@desy.de>
//
// Version: $Id: EUTelMilleGBL.cc 2299 2013-01-19 20:01:17Z hamnett $
/*!
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if GEAR and MARLINUTIL are used; require USE_GBL flag as well
#if defined(USE_GEAR) && defined(USE_MARLINUTIL) && defined(USE_GBL)

// eutelescope includes ".h"
#include "EUTelMilleGBL.h"
#include "EUTELESCOPE.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelPStream.h"
#include "EUTelReferenceHit.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelSparseCluster2Impl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelExhaustiveTrackFinder.h"
#include "EUTelGBLFitter.h"
#include "EUTelVirtualCluster.h"

// marlin includes ".h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"


// marlin util includes
#include "mille/Mille.h"

// GBL
#include "include/MilleBinary.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>


// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITree.h>
#endif


// lcio includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IO/LCWriter.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>


// ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TVector3.h>
#endif


// system includes <>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelMilleGBL::_numberTracksCandidatesHistName = "NumberTracksCandidates";
std::string EUTelMilleGBL::_chi2GblFitHistName = "chi2GblFit";
std::string EUTelMilleGBL::_probGblFitHistName = "probGblFit";
std::string EUTelMilleGBL::_residGblFitHistName = "ResidualsGblFit";
std::string EUTelMilleGBL::_residGblFitHistNameX = "ResidualsGblFit_x";
std::string EUTelMilleGBL::_residGblFitHistNameY = "ResidualsGblFit_y";
std::string EUTelMilleGBL::_resid2DGblFitHistNameXvsX = "Residuals2DGblFit_xVSx";
std::string EUTelMilleGBL::_resid2DGblFitHistNameYvsX = "Residuals2DGblFit_yVSx";
std::string EUTelMilleGBL::_resid2DGblFitHistNameXvsY = "Residuals2DGblFit_xVSy";
std::string EUTelMilleGBL::_resid2DGblFitHistNameYvsY = "Residuals2DGblFit_yVSy";
std::string EUTelMilleGBL::_kinkGblFitHistNameX = "KinksGblFit_x";
std::string EUTelMilleGBL::_kinkGblFitHistNameY = "KinksGblFit_y";
#endif

EUTelMilleGBL::EUTelMilleGBL() : Processor("EUTelMilleGBL") {

    //some default values
    FloatVec MinimalResidualsX;
    FloatVec MinimalResidualsY;
    FloatVec MaximalResidualsX;
    FloatVec MaximalResidualsY;

    FloatVec PedeUserStartValuesX;
    FloatVec PedeUserStartValuesY;

    FloatVec PedeUserStartValuesGamma;

    FloatVec SensorZPositions;

    FloatVec SensorXShifts;
    FloatVec SensorYShifts;

    FloatVec SensorGamma;

    FloatVec SensorAlpha;

    FloatVec SensorBeta;

    //maybe one has to chose a larger value than 6?
    for (int i = 0; i < 6; i++) {
        MinimalResidualsX.push_back(0.0);
        MinimalResidualsY.push_back(0.0);
        MaximalResidualsX.push_back(0.0);
        MaximalResidualsY.push_back(0.0);

        PedeUserStartValuesX.push_back(0.0);
        PedeUserStartValuesY.push_back(0.0);

        PedeUserStartValuesGamma.push_back(0.0);

        float zpos = 20000.0 + 20000.0 * (float) i;
        SensorZPositions.push_back(zpos);

        SensorXShifts.push_back(0.0);
        SensorYShifts.push_back(0.0);

        SensorGamma.push_back(0.0);
        SensorAlpha.push_back(0.0);
        SensorBeta.push_back(0.0);
    }



    // modify processor description
    _description =
            "EUTelMilleGBL uses the MILLE program to write data files for MILLEPEDE II.";

    // choose input mode
    registerOptionalParameter("InputMode", "Selects the source of input hits."
            "\n0 - hits read from hitfile with simple trackfinding. "
            "\n1 - hits read from output of tracking processor. "
            "\n2 - Test mode. Simple internal simulation and simple trackfinding. "
            "\n3 - Mixture of a track collection from the telescope and hit collections for the DUT (only one DUT layer can be used unfortunately)",
            _inputMode, static_cast<int> (0));

    registerOptionalParameter("AllowedMissingHits", "Set how many hits (=planes) can be missing on a track candidate.",
            _allowedMissingHits, static_cast<int> (0));

    registerOptionalParameter("MimosaClusterChargeMin", "Remove Mimosa26 clusters with a charge (i.e. number of fired pixels in cluster) below or equal to this value",
            _mimosa26ClusterChargeMin, static_cast<int> (1));


    // input collections
    std::vector<std::string > HitCollectionNameVecExample;
    HitCollectionNameVecExample.push_back("corrhits");

    registerInputCollections(LCIO::TRACKERHIT, "HitCollectionName",
            "Hit collections name",
            _hitCollectionName, HitCollectionNameVecExample);

    registerInputCollection(LCIO::TRACK, "TrackCollectionName",
            "Track collection name",
            _trackCollectionName, std::string("fittracks"));

    // parameters

    registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection to be saved into the output slcio file",
            _hotPixelCollectionName, static_cast<std::string> ("hotpixel_apix"));


    registerOptionalParameter("DistanceMax", "Maximal allowed distance between hits entering the fit per 10 cm space between the planes.",
            _distanceMax, static_cast<float> (2000.0));

    registerOptionalParameter("DistanceMaxVec", "Maximal allowed distance between hits entering the fit per 10 cm space between the planes. One value for each neighbor planes. DistanceMax will be used for each pair if this vector is empty.",
            _distanceMaxVec, std::vector<float> ());


    registerOptionalParameter("ExcludePlanes", "Exclude planes from fit according to their sensor ids.", _excludePlanes_sensorIDs, std::vector<int>());

    registerOptionalParameter("FixedPlanes", "Fix sensor planes in the fit according to their sensor ids.", _FixedPlanes_sensorIDs, std::vector<int>());



    registerOptionalParameter("MaxTrackCandidatesTotal", "Maximal number of track candidates (Total).", _maxTrackCandidatesTotal, static_cast<int> (10000000));

    registerOptionalParameter("MaxTrackCandidates", "Maximal number of track candidates.", _maxTrackCandidates, static_cast<int> (2000));



    registerOptionalParameter("TelescopeResolution", "Resolution of the telescope for Millepede (sigma_x=sigma_y.", _telescopeResolution, static_cast<float> (3.0));

    registerOptionalParameter("OnlySingleHitEvents", "Use only events with one hit in every plane.", _onlySingleHitEvents, static_cast<int> (0));

    registerOptionalParameter("OnlySingleTrackEvents", "Use only events with one track candidate.", _onlySingleTrackEvents, static_cast<int> (0));

    registerOptionalParameter("AlignMode", "Number of alignment constants used. Available mode are: "
            "\n1 - shifts in the X and Y directions and a rotation around the Z axis,"
            "\n2 - only shifts in the X and Y directions"
            "\n3 - (EXPERIMENTAL) shifts in the X,Y and Z directions and rotations around all three axis",
            _alignMode, static_cast<int> (1));

    registerOptionalParameter("UseResidualCuts", "Use cuts on the residuals to reduce the combinatorial background. 0 for off, 1 for on", _useResidualCuts,
            static_cast<int> (0));

    registerOptionalParameter("AlignmentConstantLCIOFile", "This is the name of the LCIO file name with the output alignment"
            "constants (add .slcio)", _alignmentConstantLCIOFile, static_cast<std::string> ("alignment.slcio"));

    registerOptionalParameter("AlignmentConstantCollectionName", "This is the name of the alignment collection to be saved into the slcio file",
            _alignmentConstantCollectionName, static_cast<std::string> ("alignment"));



    registerOptionalParameter("ResidualsXMin", "Minimal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _residualsXMin, MinimalResidualsX);

    registerOptionalParameter("ResidualsYMin", "Minimal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _residualsYMin, MinimalResidualsY);

    registerOptionalParameter("ResidualsXMax", "Maximal values of the hit residuals in the X direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _residualsXMax, MaximalResidualsX);

    registerOptionalParameter("ResidualsYMax", "Maximal values of the hit residuals in the Y direction for a track. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _residualsYMax, MaximalResidualsY);



    registerOptionalParameter("ResolutionX", "X resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _resolutionX, std::vector<float> (static_cast<int> (6), 10.));

    registerOptionalParameter("ResolutionY", "Y resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _resolutionY, std::vector<float> (static_cast<int> (6), 10.));

    registerOptionalParameter("ResolutionZ", "Z resolution parameter for each plane. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _resolutionZ, std::vector<float> (static_cast<int> (6), 10.));



    registerOptionalParameter("ReferenceCollection", "reference hit collection name ", _referenceHitCollectionName, static_cast<std::string> ("reference"));

    registerOptionalParameter("ApplyToReferenceCollection", "Do you want the reference hit collection to be corrected by the shifts and tilts from the alignment collection? (default - false )", _applyToReferenceHitCollection, static_cast<bool> (false));


    registerOptionalParameter("BinaryFilename", "Name of the Millepede binary file.", _binaryFilename, std::string("mille.bin"));
    
    registerOptionalParameter("PedeConstraintsFilename", "Name of the Millepede constraints file.", _pedeConstraintsFilename, std::string("constraints.txt"));

    registerOptionalParameter("FixParameter", "Fixes the given alignment parameters in the fit if alignMode==3 is used. For each sensor an integer must be specified (If no value is given, then all parameters will be free). bit 0 = x shift, bit 1 = y shift, bit 2 = z shift, bit 3 = alpha, bit 4 = beta, bit 5 = gamma. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.", _FixParameter, std::vector<int> (static_cast<int> (6), 24));

    registerOptionalParameter("GeneratePedeSteerfile", "Generate a steering file for the pede program.", _generatePedeSteerfile, static_cast<int> (0));

    registerOptionalParameter("PedeSteerfileName", "Name of the steering file for the pede program.", _pedeSteerfileName, std::string("steer_mille.txt"));

    registerOptionalParameter("RunPede", "Execute the pede program using the generated steering file.", _runPede, static_cast<int> (0));

    registerOptionalParameter("UsePedeUserStartValues", "Give start values for pede by hand (0 - automatic calculation of start values, 1 - start values defined by user).", _usePedeUserStartValues, static_cast<int> (0));

    registerOptionalParameter("PedeUserStartValuesX", "Start values for the alignment for shifts in the X direction.", _pedeUserStartValuesX, PedeUserStartValuesX);

    registerOptionalParameter("PedeUserStartValuesY", "Start values for the alignment for shifts in the Y direction.", _pedeUserStartValuesY, PedeUserStartValuesY);

    registerOptionalParameter("PedeUserStartValuesZ", "Start values for the alignment for shifts in the Z direction.", _pedeUserStartValuesZ, std::vector<float> (static_cast<int> (6), 0.0));

    registerOptionalParameter("PedeUserStartValuesAlpha", "Start values for the alignment for the angle alpha.", _pedeUserStartValuesAlpha, std::vector<float> (static_cast<int> (6), 0.0));

    registerOptionalParameter("PedeUserStartValuesBeta", "Start values for the alignment for the angle beta.", _pedeUserStartValuesBeta, std::vector<float> (static_cast<int> (6), 0.0));

    registerOptionalParameter("PedeUserStartValuesGamma", "Start values for the alignment for the angle gamma.", _pedeUserStartValuesGamma, PedeUserStartValuesGamma);

    registerOptionalParameter("TestModeSensorResolution", "Resolution assumed for the sensors in test mode.", _testModeSensorResolution, static_cast<float> (3.0));

    registerOptionalParameter("TestModeXTrackSlope", "Width of the track slope distribution in the x direction", _testModeXTrackSlope, static_cast<float> (0.0005));

    registerOptionalParameter("TestModeYTrackSlope", "Width of the track slope distribution in the y direction", _testModeYTrackSlope, static_cast<float> (0.0005));

    registerOptionalParameter("TestModeSensorZPositions", "Z positions of the sensors in test mode.", _testModeSensorZPositions, SensorZPositions);

}

void EUTelMilleGBL::init() {

    // check if the GEAR manager pointer is not null!
    if (marlin::Global::GEAR == 0x0) {
        streamlog_out(ERROR2) << "The GearMgr is not available, for an unknown reason." << std::endl;
        throw InvalidGeometryException("GEAR manager is not initialised");
    }


    // sensor-planes in geometry navigation:
    _siPlanesParameters = const_cast<gear::SiPlanesParameters*> (&(marlin::Global::GEAR->getSiPlanesParameters()));
    _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> (&(_siPlanesParameters->getSiPlanesLayerLayout()));

    // clear the sensor ID vector
    _sensorIDVec.clear();

    // clear the sensor ID map
    _sensorIDVecMap.clear();
    _sensorIDtoZOrderMap.clear();

    // clear the sensor ID vector (z-axis order)
    _sensorIDVecZOrder.clear();

    // copy-paste from another class (should be ideally part of GEAR!)
    double* keepZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];

    for (int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++) {
        int sensorID = _siPlanesLayerLayout->getID(iPlane);
        _sensorIDVec.push_back(sensorID);
        _sensorIDVecMap.insert(std::make_pair(sensorID, iPlane));

        // count number of the sensors to the left of the current one:
        int sensorsToTheLeft = 0;
        keepZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
        for (int jPlane = 0; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++)
            if (_siPlanesLayerLayout->getLayerPositionZ(jPlane) + 1e-06 < keepZPosition[ iPlane ]) sensorsToTheLeft++;

        _sensorIDVecZOrder.push_back(sensorsToTheLeft);
        _sensorIDtoZOrderMap.insert(std::make_pair(sensorID, sensorsToTheLeft));
    }

    delete [] keepZPosition;

    _histogramSwitch = true;

    _referenceHitVec = 0;


    //lets guess the number of planes
    if (_inputMode == 0) {
        // the number of planes is got from the GEAR description and is
        // the sum of the telescope reference planes and the DUT (if any)
        _nPlanes = _siPlanesParameters->getSiPlanesNumber();
        if (_siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT) ++_nPlanes;

        if (_useSensorRectangular.empty()) {
            streamlog_out(MESSAGE4) << "No rectangular limits on pixels of sensorplanes applied" << std::endl;
        } else {
            if (_useSensorRectangular.size() % 5 != 0) {
                streamlog_out(WARNING2) << "Wrong number of arguments in RectangularLimits! Ignoring this cut!" << std::endl;
            } else {
                streamlog_out(MESSAGE4) << "Reading in SensorRectangularCuts: " << std::endl;
                int sensorcuts = _useSensorRectangular.size() / 5;
                for (int i = 0; i < sensorcuts; ++i) {
                    int sensor = _useSensorRectangular.at(5 * i + 0);
                    int A = _useSensorRectangular.at(5 * i + 1);
                    int B = _useSensorRectangular.at(5 * i + 2);
                    int C = _useSensorRectangular.at(5 * i + 3);
                    int D = _useSensorRectangular.at(5 * i + 4);
                    Utility::SensorRectangular r(sensor, A, B, C, D);
                    r.print();
                    _rect.addRectangular(r);
                }
            }
        }
    } else if (_inputMode == 1) _nPlanes = _siPlanesParameters->getSiPlanesNumber();

    else {
        streamlog_out(ERROR2) << "unknown input mode " << _inputMode << std::endl;
        throw InvalidParameterException("unknown input mode");
    }


    // an associative map for getting also the sensorID ordered
    std::map< double, int > sensorIDMap;

    // create an array with the z positions of each layer
    for (int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++) {
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getLayerPositionZ(iPlane));
        sensorIDMap.insert(std::make_pair(_siPlanesLayerLayout->getLayerPositionZ(iPlane), _siPlanesLayerLayout->getID(iPlane)));
    }

    if (_siPlanesParameters->getSiPlanesType() == _siPlanesParameters->TelescopeWithDUT) {
        _siPlaneZPosition.push_back(_siPlanesLayerLayout->getDUTPositionZ());
        sensorIDMap.insert(std::make_pair(_siPlanesLayerLayout->getDUTPositionZ(), _siPlanesLayerLayout->getDUTID()));
    }

    // sort the array with increasing z
    sort(_siPlaneZPosition.begin(), _siPlaneZPosition.end());


    //the user is giving sensor ids for the planes to be excluded. this
    //sensor ids have to be converted to a local index according to the
    //planes positions along the z axis.
    for (size_t i = 0; i < _FixedPlanes_sensorIDs.size(); i++) {
        std::map< double, int >::iterator iter = sensorIDMap.begin();
        int counter = 0;
        while (iter != sensorIDMap.end()) {
            if (iter->second == _FixedPlanes_sensorIDs[i]) {
                _FixedPlanes.push_back(counter);
                break;
            }
            ++iter;
            ++counter;
        }
    }
    for (size_t i = 0; i < _excludePlanes_sensorIDs.size(); i++) {
        std::map< double, int >::iterator iter = sensorIDMap.begin();
        int counter = 0;
        while (iter != sensorIDMap.end()) {
            if (iter->second == _excludePlanes_sensorIDs[i]) {
                _excludePlanes.push_back(counter);
                break;
            }
            ++iter;
            ++counter;
        }
    }

    // strip from the map the sensor id already sorted.
    std::map< double, int >::iterator iter = sensorIDMap.begin();
    unsigned int counter = 0;
    while (iter != sensorIDMap.end()) {
        bool excluded = false;
        for (size_t i = 0; i < _excludePlanes.size(); i++) {
            if (_excludePlanes[i] == counter) {
                excluded = true;
                break;
            }
        }
        if (!excluded)
            _orderedSensorID_wo_excluded.push_back(iter->second);
        _orderedSensorID.push_back(iter->second);

        ++iter;
        ++counter;
    }
    //


    //consistency
    if (_siPlaneZPosition.size() != _nPlanes) {
        streamlog_out(ERROR2) << "the number of detected planes is " << _nPlanes << " but only " << _siPlaneZPosition.size() << " layer z positions were found!" << std::endl;
        throw InvalidParameterException("number of layers and layer z positions mismatch");
    }

    // this method is called only once even when the rewind is active
    // usually a good idea to
    printParameters();

    // set to zero the run and event counters
    _iRun = 0;
    _iEvt = 0;

    // Initialize number of excluded planes
    _nExcludePlanes = _excludePlanes.size();

    streamlog_out(MESSAGE2) << "Number of planes excluded from the alignment fit: " << _nExcludePlanes << std::endl;

    // Initialise Mille statistics
    _nMilleDataPoints = 0;
    _nMilleTracks = 0;

    _waferResidX = new double[_nPlanes];
    _waferResidY = new double[_nPlanes];
    _waferResidZ = new double[_nPlanes];


    _xFitPos = new double[_nPlanes];
    _yFitPos = new double[_nPlanes];

    _telescopeResolX = new double[_nPlanes];
    _telescopeResolY = new double[_nPlanes];
    _telescopeResolZ = new double[_nPlanes];

    // fill resolution arrays
    for (size_t help = 0; help < _nPlanes; help++) {
        _telescopeResolX[help] = _telescopeResolution;
        _telescopeResolY[help] = _telescopeResolution;
    }

    // booking histograms
    bookHistos();

    streamlog_out(MESSAGE2) << "Initialising Mille..." << std::endl;

    unsigned int reserveSize = 80000;
    _milleGBL = new gbl::MilleBinary(_binaryFilename, reserveSize);

    for (int i = 0; i < _maxTrackCandidates; i++) {
        _xPos.push_back(std::vector<double>(_nPlanes, 0.0));
        _yPos.push_back(std::vector<double>(_nPlanes, 0.0));
        _zPos.push_back(std::vector<double>(_nPlanes, 0.0));
    }

    if (!_distanceMaxVec.empty()) {
        if (_distanceMaxVec.size() != static_cast<unsigned int> (_nPlanes)) {
            streamlog_out(WARNING2) << "Consistency check of the DistanceMaxVec array failed. Its size is different compared to the number of planes! Will now use _distanceMax for each pair of planes." << std::endl;
            _distanceMaxVec.clear();
            for (size_t i = 0; i < _nPlanes - 1; i++) {
                _distanceMaxVec.push_back(_distanceMax);
            }
        }
    } else {
        _distanceMaxVec.clear();
        for (size_t i = 0; i < _nPlanes - 1; i++) {
            _distanceMaxVec.push_back(_distanceMax);
        }
    }

    streamlog_out(MESSAGE) << "Initialisation of track finder" << std::endl;

    EUTelExhaustiveTrackFinder* myFinder = new EUTelExhaustiveTrackFinder("MyFinder", getAllowedMissingHits(), getMaxTrackCandidates());
    myFinder->SetMode(2);
    myFinder->SetDistanceMaxVec(_distanceMaxVec);
    myFinder->SetResidualsYMax(_residualsYMax);
    myFinder->SetResidualsYMin(_residualsYMin);
    myFinder->SetResidualsXMax(_residualsXMax);
    myFinder->SetResidualsXMin(_residualsXMin);
    _theFinder = myFinder;

    if (_theFinder == 0) {
        streamlog_out(ERROR) << "Can't allocate an instance of myTrackFinder. Stopping ..." << std::endl;
        //      return ;
    }

    streamlog_out(MESSAGE) << "Initialisation of track fitter" << std::endl;

    EUTelGBLFitter* myFitter = new EUTelGBLFitter("myGBLFitter");
    myFitter->SetMilleBinary(_milleGBL);
    _theFitter = myFitter;

    if (_theFitter == 0) {
        streamlog_out(ERROR) << "Can't allocate an instance of myTrackFitter. Stopping ..." << std::endl;
        //      return ;
    }


    streamlog_out(MESSAGE2) << "end of initialisation" << std::endl;
}

void EUTelMilleGBL::processRunHeader(LCRunHeader * rdr) {

    std::auto_ptr<EUTelRunHeaderImpl> header(new EUTelRunHeaderImpl(rdr));
    header->addProcessor(type());

    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, instead of barely
    // quitting ask the user what to do.

    if (header->getGeoID() != _siPlanesParameters->getSiPlanesID()) {
        streamlog_out(ERROR2) << "Error during the geometry consistency check: " << std::endl;
        streamlog_out(ERROR2) << "The run header says the GeoID is " << header->getGeoID() << std::endl;
        streamlog_out(ERROR2) << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesNumber() << std::endl;

#ifdef EUTEL_INTERACTIVE
        string answer;
        while (true) {
            streamlog_out(ERROR2) << "Type Q to quit now or C to continue using the actual GEAR description anyway [Q/C]" << std::endl;
            cin >> answer;
            // put the answer in lower case before making the comparison.
            transform(answer.begin(), answer.end(), answer.begin(), ::tolower);
            if (answer == "q") {
                exit(-1);
            } else if (answer == "c") {
                break;
            }
        }
#endif

    }

    // increment the run counter
    ++_iRun;
}

void EUTelMilleGBL::processEvent(LCEvent * event) {

    if (isFirstEvent()) {

        _hotPixelMap = Utility::FillHotPixelMap(event, _hotPixelCollectionName);

        if (_applyToReferenceHitCollection) _referenceHitVec = dynamic_cast<LCCollectionVec*> (event->getCollection(_referenceHitCollectionName));
    }

    if (_iEvt % 100 == 0 && _iEvt % 1000 != 0) {
        streamlog_out(MESSAGE4) << "Processing event "
                << std::setw(6) << std::setiosflags(std::ios::right) << event->getEventNumber() << " in run "
                << std::setw(6) << std::setiosflags(std::ios::right) << std::setfill('0') << event->getRunNumber() << std::setfill(' ')
                << " (Total = " << std::setw(10) << _iEvt << ")" << std::resetiosflags(std::ios::left) << std::endl;
    }

    if (_iEvt % 1000 == 0) {
        streamlog_out(MESSAGE6) << "Processing event "
                << std::setw(6) << std::setiosflags(std::ios::right) << event->getEventNumber() << " in run "
                << std::setw(6) << std::setiosflags(std::ios::right) << std::setfill('0') << event->getRunNumber() << std::setfill(' ')
                << " (Total = " << std::setw(10) << _iEvt << ")" << std::resetiosflags(std::ios::left) << std::endl;
    }

    // stop when necessary amount of tracks reached
    if (_nMilleTracks > _maxTrackCandidatesTotal) throw marlin::StopProcessingException(this);

    EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

    if (evt->getEventType() == kEORE) {
        streamlog_out(DEBUG2) << "EORE found: nothing else to do." << std::endl;
        return;
    }

    std::vector< EVENT::TrackerHitVec > _allHitsArray(_nPlanes, EVENT::TrackerHitVec());
    std::vector<int> indexconverter(_nPlanes, -1);

    Utility::FillNotExcludedPlanesIndices(indexconverter, _excludePlanes, _nPlanes);

    // check if running in input mode 0 or 2
    if (_inputMode == 0) {

        for (size_t i = 0; i < _hitCollectionName.size(); i++) {

            LCCollection* collection;
            try {
                collection = event->getCollection(_hitCollectionName[i]);
            } catch (DataNotAvailableException& e) {
                streamlog_out(WARNING2) << "No input collection " << _hitCollectionName[i] << " found for event " << event->getEventNumber()
                        << " in run " << event->getRunNumber() << std::endl;
                throw marlin::SkipEventException(this);
            }
            int layerIndex = -1;

            // loop over all hits in collection
            for (int iHit = 0; iHit < collection->getNumberOfElements(); iHit++) {

                TrackerHitImpl * hit = static_cast<TrackerHitImpl*> (collection->getElementAt(iHit));
                ///code hoes here;
                if (Utility::HitContainsHotPixels(hit, _hotPixelMap)) {
                    streamlog_out(DEBUG3) << "Hit " << i << " contains hot pixels; skip this one. " << std::endl;
                    continue;
                }

                LCObjectVec clusterVector = hit->getRawHits();

                EUTelVirtualCluster * cluster;
                cluster = Utility::GetClusterFromHit(hit);

                //Warning! It is very bad idea to do like this. This code must ideally be implemented in GetClusterFromHit
                if (hit->getType() == kEUTelSparseClusterImpl) {
                    // ok the cluster is of sparse type, but we also need to know
                    // the kind of pixel description used. This information is
                    // stored in the corresponding original data collection.
                    LCCollectionVec * sparseClusterCollectionVec = dynamic_cast<LCCollectionVec *> (evt->getCollection("original_zsdata"));

                    TrackerDataImpl * oneCluster = dynamic_cast<TrackerDataImpl*> (sparseClusterCollectionVec->getElementAt(0));
                    CellIDDecoder<TrackerDataImpl > anotherDecoder(sparseClusterCollectionVec);
                    SparsePixelType pixelType = static_cast<SparsePixelType> (static_cast<int> (anotherDecoder(oneCluster)["sparsePixelType"]));

                    // now we know the pixel type. So we can properly create a new
                    // instance of the sparse cluster
                    if (pixelType == kEUTelSimpleSparsePixel) {
                        cluster = new EUTelSparseClusterImpl< EUTelSimpleSparsePixel > (static_cast<TrackerDataImpl *> (clusterVector[ 0 ]));
                    } else if (pixelType == kEUTelAPIXSparsePixel) {
                        cluster = new EUTelSparseClusterImpl<EUTelAPIXSparsePixel > (static_cast<TrackerDataImpl *> (clusterVector[ 0 ]));
                    } else {
                        streamlog_out(ERROR4) << "Unknown pixel type.  Sorry for quitting." << std::endl;
                        throw UnknownDataTypeException("Pixel type unknown");
                    }
                }

                if (hit->getType() == kEUTelAPIXClusterImpl) {
                    TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> (clusterVector[0]);
                    EUTelSparseClusterImpl< EUTelAPIXSparsePixel > *apixCluster = new EUTelSparseClusterImpl< EUTelAPIXSparsePixel > (clusterFrame);
                    int sensorID = apixCluster->getDetectorID();
                    bool skipHit = false;
                    for (size_t iPixel = 0; iPixel < apixCluster->size(); ++iPixel) {
                        EUTelAPIXSparsePixel apixPixel;
                        apixCluster->getSparsePixelAt(iPixel, &apixPixel);
                        int pixelX = apixPixel.getXCoord();
                        int pixelY = apixPixel.getYCoord();
                        skipHit = !_rect.isInside(sensorID, pixelX, pixelY);
                    }

                    if (skipHit) {
                        streamlog_out(MESSAGE4) << "Skipping cluster due to rectangular cut!" << std::endl;
                        continue;
                    }
                    if (skipHit) {
                        streamlog_out(MESSAGE4) << "TEST: This message should never appear!" << std::endl;
                    }
                }

                if (hit->getType() == kEUTelDFFClusterImpl ||
                        hit->getType() == kEUTelFFClusterImpl ||
                        hit->getType() == kEUTelSparseClusterImpl) {
                    if (cluster->getTotalCharge() <= getMimosa26ClusterChargeMin()) {
                        streamlog_out(DEBUG) << " Thin cluster (charge <=" << getMimosa26ClusterChargeMin()
                                << ") found and removed (hit type " << hit->getType()
                                << " on detector w/ id " << cluster->getDetectorID() << ")" << std::endl;
                        delete cluster;
                        continue;
                    }
                }

                unsigned int localSensorID = Utility::GuessSensorID(hit);
                layerIndex = _sensorIDVecMap[ localSensorID ];
                _allHitsArray[layerIndex].push_back(hit);

                delete cluster; // <--- destroying the cluster
            } // end loop over all hits in collection

        }
    }

    int _nTracks = 0;
    int _nGoodTracks = 0;

    std::map< int, EVENT::TrackerHitVec > trackCandidates;

    if (_inputMode == 0) {

        streamlog_out(DEBUG1) << "Event #" << _iEvt << std::endl;
        streamlog_out(DEBUG1) << "Initialising hits for _theFinder..." << std::endl;
        TrackerHitVec::iterator itHit;
        for (int iLayer = 0; iLayer < 6; iLayer++) {
            streamlog_out(DEBUG0) << iLayer << std::endl;
            for (itHit = _allHitsArray[iLayer].begin(); itHit != _allHitsArray[iLayer].end(); ++itHit) {
                streamlog_out(DEBUG0) << *itHit << std::endl;
                streamlog_out(DEBUG0) << std::setw(10) << (*itHit)->id() << std::endl;
                streamlog_out(DEBUG0) << std::setw(10) << (*itHit)->getPosition()[0]
                        << std::setw(10) << (*itHit)->getPosition()[1]
                        << std::setw(10) << (*itHit)->getPosition()[2] << std::endl;
            }
        }

        _theFinder->SetAllHits(_allHitsArray);
        streamlog_out(DEBUG1) << "Searching tracks..." << std::endl;
        EUTelTrackFinder::SearchResult search_result = _theFinder->SearchTracks();
        streamlog_out(DEBUG1) << "Search results retrieved..." << std::endl;
        streamlog_out(DEBUG1) << "Search results = " << (int) search_result << std::endl;
        streamlog_out(DEBUG1) << "Retrieving track candidates..." << std::endl;
        trackCandidates = _theFinder->GetTrackCandidates();
        std::vector<std::vector<int> > indexarray;

        _nTracks = (int) trackCandidates.size();
        static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _numberTracksCandidatesHistName ]) -> fill(_nTracks);
        //        
        streamlog_out(DEBUG1) << "Event #" << _iEvt << std::endl;
        streamlog_out(DEBUG1) << "Track finder " << _theFinder->GetName() << " found  " << _nTracks << std::endl;

        // Perform fit for all found track candidates
        // ------------------------------------------
        unsigned int numData;
        TVectorD residual(200);
        TVectorD measErr(200);
        TVectorD residualErr(200);
        TVectorD downWeight(200);
        if (_nTracks != 0 && _nTracks == 1) {
            _theFitter->SetTrackCandidates(trackCandidates);
            _theFitter->FitTracks();
//
            double chi2Trk = 0.;
            int ndfTrk = 0;

            IMPL::LCCollectionVec* fittrackvec;
            fittrackvec = static_cast<EUTelGBLFitter*> (_theFitter)->GetFitTrackVec();
            IMPL::LCCollectionVec::const_iterator itFitTrack;

            int iCounter = 0;
            for (itFitTrack = fittrackvec->begin(); itFitTrack != fittrackvec->end(); ++itFitTrack) {
                chi2Trk = static_cast<TrackImpl*> (*itFitTrack)->getChi2();
                ndfTrk = static_cast<TrackImpl*> (*itFitTrack)->getNdf();

                static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _chi2GblFitHistName ]) -> fill(chi2Trk);
                if( chi2Trk < 17 ) static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ _probGblFitHistName ]) -> fill(TMath::Prob(chi2Trk, ndfTrk));


                std::map< int, gbl::GblTrajectory* > gblTracks = static_cast<EUTelGBLFitter*> (_theFitter)->GetGblTrackCandidates();

                std::stringstream sstr;
                gbl::GblTrajectory* gblTraj = gblTracks[ iCounter ];
		//gblTraj->printTrajectory();
		//gblTraj->printPoints();
		//gblTraj->printData();
                std::vector< gbl::GblPoint > gblPointVec = static_cast<EUTelGBLFitter*> (_theFitter)->GetGblTracksPoints()[iCounter];
                std::vector< gbl::GblPoint >::const_iterator itGblPoint = gblPointVec.begin();
                int iPlane = 0; // wrong in case of missing planes
                for (; itGblPoint != gblPointVec.end(); ++itGblPoint) {
                    if (iPlane > 5) continue;
		    //if ( itGblPoint->getLabel() < 1000 )
		    if ( itGblPoint->getLabel() % 3 == 1 ) 
		    {
			    streamlog_out(DEBUG0) << std::setw(15) << itGblPoint->getLabel() << std::endl;
                            // spatial residuals
			    gblTraj->getMeasResults(itGblPoint->getLabel(), numData, residual, measErr, residualErr, downWeight);
			    sstr << _residGblFitHistNameX << iPlane;
			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[0] << std::setw(15) << std::setprecision(5) << residualErr[0] << std::endl;
			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[0] / residualErr[0], downWeight[0]);
			    sstr.str(std::string());
			    sstr << _residGblFitHistNameY << iPlane;
			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[1] << std::setw(15) << std::setprecision(5) << residualErr[1] << std::endl;
			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[1] / residualErr[1], downWeight[1]);
			    sstr.str(std::string());
                            // kinks
                            gblTraj->getScatResults (itGblPoint->getLabel(), numData, residual, measErr, residualErr, downWeight);
                            sstr << _kinkGblFitHistNameX << iPlane;
			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[0] << std::setw(15) << std::setprecision(5) << residualErr[0] << std::endl;
			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[0], downWeight[0]);
			    sstr.str(std::string());
			    sstr << _kinkGblFitHistNameY << iPlane;
			    streamlog_out(DEBUG0) << std::setw(15) << std::setprecision(5) << residual[1] << std::setw(15) << std::setprecision(5) << residualErr[1] << std::endl;
			    static_cast<AIDA::IHistogram1D*> (_aidaHistoMap1D[ sstr.str() ]) -> fill(residual[1], downWeight[1]);
			    sstr.str(std::string());
                            
                            // 2D histograms
                            const double* hitpos = trackCandidates.begin()->second[iPlane]->getPosition();              // wrong in case of empty planes
                            sstr << _resid2DGblFitHistNameXvsX << iPlane;
                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[0], hitpos[0], downWeight[0]);
                            sstr.str(std::string());
                            sstr << _resid2DGblFitHistNameXvsY << iPlane;
                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[0], hitpos[1], downWeight[0]);
                            sstr.str(std::string());
                            sstr << _resid2DGblFitHistNameYvsX << iPlane;
                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[1], hitpos[0], downWeight[1]);
                            sstr.str(std::string());
                            sstr << _resid2DGblFitHistNameYvsY << iPlane;
                            static_cast<AIDA::IHistogram2D*> (_aidaHistoMap2D[ sstr.str() ]) -> fill(residual[1], hitpos[1], downWeight[1]);
                            sstr.str(std::string());
                            
			    if ( itGblPoint->getLabel() < 1000 )++iPlane;
		    }
                }

                IMPL::LCCollectionVec::const_iterator itFitTrack;
                iCounter++;
            }
        }
        _theFinder->Reset(); // Breaks Fitter track candidates reference???!!!
    }

    // count events
    ++_iEvt;
    if (isFirstEvent()) _isFirstEvent = false;

}

int EUTelMilleGBL::guessSensorID(double * hit) {

    int sensorID = -1;
    double minDistance = std::numeric_limits< double >::max();

    if (_referenceHitVec == 0 || _applyToReferenceHitCollection == false) {
        // use z information of planes instead of reference vector
        for (int iPlane = 0; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane) {
            double distance = std::abs(hit[2] - _siPlaneZPosition[ iPlane ]);
            if (distance < minDistance) {
                minDistance = distance;
                sensorID = _siPlanesLayerLayout->getID(iPlane);
            }
        }
        if (minDistance > 30) {
            // advice the user that the guessing wasn't successful 
            streamlog_out(WARNING3) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
                    "Please check the consistency of the data with the GEAR file: hitPosition[2]=" << hit[2] << std::endl;
        }

        return sensorID;
    }

    for (size_t ii = 0; ii < (unsigned int) _referenceHitVec->getNumberOfElements(); ii++) {
        EUTelReferenceHit* refhit = static_cast<EUTelReferenceHit*> (_referenceHitVec->getElementAt(ii));

        TVector3 hit3d(hit[0], hit[1], hit[2]);
        TVector3 hitInPlane(refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane(refhit->getAlpha(), refhit->getBeta(), refhit->getGamma());

        double distance = fabs(norm2Plane.Dot(hit3d - hitInPlane));
        if (distance < minDistance) {
            minDistance = distance;
            sensorID = refhit->getSensorID();
        }

    }

    return sensorID;
}

void EUTelMilleGBL::end() {

    delete [] _telescopeResolY;
    delete [] _telescopeResolX;
    delete [] _telescopeResolZ;
    delete [] _yFitPos;
    delete [] _xFitPos;
    delete [] _waferResidY;
    delete [] _waferResidX;
    delete [] _waferResidZ;

    delete _theFinder;
    delete _theFitter;

    // if write the pede steering file
    if (_generatePedeSteerfile) {
        streamlog_out ( MESSAGE4 ) << "\n Generating the steering file for the pede program..." << std::endl;

        std::string tempHistoName;
        double *meanX = new double[_nPlanes];
        double *meanY = new double[_nPlanes];
        double *meanZ = new double[_nPlanes];

//        // loop over all detector planes
//        for(unsigned int iDetector = 0; iDetector < _nPlanes; iDetector++ ) {
//
//          int sensorID = _orderedSensorID.at( iDetector );
//
//          if ( _histogramSwitch ) {
//            tempHistoName =  _residualXLocalname + "_d" + to_string( sensorID );
//            if ( AIDA::IHistogram1D* residx_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
//              meanX[iDetector] = residx_histo->mean();
//            else {
//              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualXLocalname << endl;
//              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
//              _histogramSwitch = false;
//            }
//          }
//
//          if ( _histogramSwitch ) {
//            tempHistoName =  _residualYLocalname + "_d" + to_string( sensorID );
//            if ( AIDA::IHistogram1D* residy_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
//              meanY[iDetector] = residy_histo->mean();
//            else {
//              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualYLocalname << endl;
//              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
//              _histogramSwitch = false;
//            }
//          }
//
//          if ( _histogramSwitch ) {
//            tempHistoName =  _residualZLocalname + "_d" + to_string( sensorID );
//            if ( AIDA::IHistogram1D* residz_histo = dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[tempHistoName.c_str()]) )
//              meanZ[iDetector] = residz_histo->mean();
//            else {
//              streamlog_out ( ERROR2 ) << "Not able to retrieve histogram pointer for " << _residualZLocalname << endl;
//              streamlog_out ( ERROR2 ) << "Disabling histogramming from now on" << endl;
//              _histogramSwitch = false;
//            } 
//          }
//        } // end loop over all detector planes

        std::ofstream steerFile;
        steerFile.open(_pedeSteerfileName.c_str());

        if (steerFile.is_open()) {

          // find first and last excluded plane
          unsigned int firstnotexcl = _nPlanes;
          unsigned int lastnotexcl = 0;

          // loop over all planes
          for (unsigned int help = 0; help < _nPlanes; help++) 
          {

            int excluded = 0;

            // loop over all excluded planes
            for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) 
            {
              if (help == _excludePlanes[helphelp]) 
              {
    //            excluded = 1;
              }
            } // end loop over all excluded planes

            if (excluded == 0 && firstnotexcl > help) 
            {
              firstnotexcl = help;
            }

            if (excluded == 0 && lastnotexcl < help) 
            {
              lastnotexcl = help;
            }
          } // end loop over all planes

          // calculate average
          double averageX = (meanX[firstnotexcl] + meanX[lastnotexcl]) / 2;
          double averageY = (meanY[firstnotexcl] + meanY[lastnotexcl]) / 2;
          double averageZ = (meanZ[firstnotexcl] + meanZ[lastnotexcl]) / 2;

          steerFile << "fortranfiles" << std::endl;
          steerFile << _pedeConstraintsFilename << std::endl;
          steerFile << std::endl;
          
          steerFile << "Cfiles" << std::endl;
          steerFile << _binaryFilename << std::endl;
          steerFile << std::endl;

          steerFile << "Parameter" << std::endl;

          int counter = 0;

          // loop over all planes
          for (unsigned int help = 0; help < _nPlanes; help++) {

            int excluded = 0; // flag for excluded planes

            // loop over all excluded planes
            for (int helphelp = 0; helphelp < _nExcludePlanes; helphelp++) {

              if (help == _excludePlanes[helphelp]) {
    //            excluded = 1;
              }

            } // end loop over all excluded planes

            // if plane not excluded
            if (excluded == 0) {

              bool fixed = false;
              for(size_t i = 0;i< _FixedPlanes.size(); i++)
                {
                  if(_FixedPlanes[i] == static_cast< int >(help))
                    fixed = true;
                }

              if( fixed || (_FixedPlanes.empty() && (help == firstnotexcl || help == lastnotexcl) ) )
                {
                  if (_alignMode == 1) {
                    steerFile << (counter * 3 + 1) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 3 + 2) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 3 + 3) << " 0.0 -1.0" << std::endl;
                  } else if (_alignMode == 2) {
                    steerFile << (counter * 2 + 1) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 2 + 2) << " 0.0 -1.0" << std::endl;
                  } else if (_alignMode == 3) {
                    steerFile << (counter * 6 + 1) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 6 + 2) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 6 + 3) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 6 + 4) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 6 + 5) << " 0.0 -1.0" << std::endl;
                    steerFile << (counter * 6 + 6) << " 0.0 -1.0" << std::endl;
                  } else if (_alignMode == 4) {
                        if(help==0)
                        {       
                            steerFile << 1 << " " << _pedeUserStartValuesX[help] << " -1.0" << std::endl;
                            steerFile << 2 << " " << _pedeUserStartValuesY[help] << " -1.0" << std::endl;
                        }
                  }
                } else {

                if (_alignMode == 1) {

                  if (_usePedeUserStartValues == 0) {
                    steerFile << (counter * 3 + 1) << " " << (averageX - meanX[help]) << " 0.0" << std::endl;
                    steerFile << (counter * 3 + 2) << " " << (averageY - meanY[help]) << " 0.0" << std::endl;
                    steerFile << (counter * 3 + 3) << " " << " 0.0 0.0" << std::endl;
                  } else {
                    steerFile << (counter * 3 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << std::endl;
                    steerFile << (counter * 3 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << std::endl;
                    steerFile << (counter * 3 + 3) << " " << _pedeUserStartValuesGamma[help] << " 0.0" << std::endl;
                  }

                } else if (_alignMode == 2) {

                  if (_usePedeUserStartValues == 0) {
                    steerFile << (counter * 2 + 1) << " " << (averageX - meanX[help]) << " 0.0" << std::endl;
                    steerFile << (counter * 2 + 2) << " " << (averageY - meanY[help]) << " 0.0" << std::endl;
                  } else {
                    steerFile << (counter * 2 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << std::endl;
                    steerFile << (counter * 2 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << std::endl;
                  }

                } else if (_alignMode == 3) {
                  if (_usePedeUserStartValues == 0)
                    {
                      if(_FixParameter[help] & (1 << 0))
                        steerFile << (counter * 6 + 1) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 1) << " " << (averageX - meanX[help]) << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 1))
                        steerFile << (counter * 6 + 2) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 2) << " " << (averageY - meanY[help]) << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 2))
                        steerFile << (counter * 6 + 3) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 3) << " " << (averageZ - meanZ[help]) << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 3))
                        steerFile << (counter * 6 + 4) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 4) << " 0.0 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 4))
                        steerFile << (counter * 6 + 5) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 5) << " 0.0 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 5))
                        steerFile << (counter * 6 + 6) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 6) << " 0.0 0.0" << std::endl;
                    } else {
                      if(_FixParameter[help] & (1 << 0))
                        steerFile << (counter * 6 + 1) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 1) << " " << _pedeUserStartValuesX[help] << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 1))
                        steerFile << (counter * 6 + 2) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 2) << " " << _pedeUserStartValuesY[help] << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 2))
                        steerFile << (counter * 6 + 3) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 3) << " " << _pedeUserStartValuesZ[help] << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 3))
                        steerFile << (counter * 6 + 4) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 4) << " " << _pedeUserStartValuesAlpha[help] << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 4))
                        steerFile << (counter * 6 + 5) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 5) << " " << _pedeUserStartValuesBeta[help] << " 0.0" << std::endl;

                      if(_FixParameter[help] & (1 << 5))
                        steerFile << (counter * 6 + 6) << " 0.0 -1.0" << std::endl;
                      else
                        steerFile << (counter * 6 + 6) << " " << _pedeUserStartValuesGamma[help] << " 0.0" << std::endl;
                    }
                } else if (_alignMode == 4) {
                      if (_usePedeUserStartValues == 0) {
                        if(help==0)
                        {       
                            steerFile << 1 << " " << _pedeUserStartValuesX[help] << " -1.0" << std::endl;
                            steerFile << 2 << " " << _pedeUserStartValuesY[help] << " -1.0" << std::endl;
                        }
                      } 
                }
              }

              counter++;

            } // end if plane not excluded

          } // end loop over all planes

          steerFile << std::endl;
          steerFile << "chiscut 10.0 5.0" << std::endl;
          steerFile << "! outlierdownweighting 4" << std::endl;
          steerFile << std::endl;
          steerFile << "method inversion 10 0.001" << std::endl;
          steerFile << std::endl;
          steerFile << "! histprint" << std::endl;
          steerFile << std::endl;
          steerFile << "end" << std::endl;

          steerFile.close();

          streamlog_out ( MESSAGE5 ) << "File " << _pedeSteerfileName << " written." << std::endl;

        } else {

          streamlog_out ( ERROR2 ) << "Could not open steering file." << std::endl;

        }
        //cleaning up
        delete [] meanX;
        delete [] meanY;
        delete [] meanZ;
    } // end if write the pede steering file

    streamlog_out(MESSAGE2) << std::endl;
    streamlog_out(MESSAGE2) << "Number of data points used: " << _nMilleDataPoints << std::endl;
    streamlog_out(MESSAGE2) << "Number of tracks used: " << _nMilleTracks << std::endl;

    // if running pede using the generated steering file
    if (_runPede == 1) {
        // check if steering file exists
        if (_generatePedeSteerfile == 1) {

          std::string command = "pede " + _pedeSteerfileName;

          streamlog_out ( MESSAGE5 ) << "Starting pede...: " << command.c_str() << std::endl;

          bool encounteredError = false;

          // run pede and create a streambuf that reads its stdout and stderr
          redi::ipstream pede( command.c_str(), redi::pstreams::pstdout|redi::pstreams::pstderr ); 

          if (!pede.is_open()) {
              streamlog_out( ERROR5 ) << "Pede cannot be executed: command not found in the path" << std::endl;
              encounteredError = true;	  
          } else {

            // output multiplexing: parse pede output in both stdout and stderr and echo messages accordingly
            char buf[1024];
            std::streamsize n;
            std::stringstream pedeoutput; // store stdout to parse later
            std::stringstream pedeerrors;
            bool finished[2] = { false, false };
            while (!finished[0] || !finished[1])
              {
                if (!finished[0])
                  {
                    while ((n = pede.err().readsome(buf, sizeof(buf))) > 0){
                      streamlog_out( ERROR5 ).write(buf, n).flush();
                      std::string error (buf, n);
                      pedeerrors << error;
                      encounteredError = true;
                    }
                    if (pede.eof())
                      {
                        finished[0] = true;
                        if (!finished[1])
                          pede.clear();
                      }
                  }

                if (!finished[1])
                  {
                    while ((n = pede.out().readsome(buf, sizeof(buf))) > 0){
                      streamlog_out( MESSAGE4 ).write(buf, n).flush();
                      std::string output (buf, n);
                      pedeoutput << output;
                    }
                    if (pede.eof())
                      {
                        finished[1] = true;
                        if (!finished[0])
                          pede.clear();
                      }
                  }
              }

            // pede does not return exit codes on some errors (in V03-04-00)
            // check for some of those here by parsing the output
            {
              const char * pch = strstr(pedeoutput.str().data(),"Too many rejects");
              if (pch){
                streamlog_out ( ERROR5 ) << "Pede stopped due to the large number of rejects. " << std::endl;
                encounteredError = true;
              }
            }

            {
              const char* pch0 = strstr(pedeoutput.str().data(),"Sum(Chi^2)/Sum(Ndf) = ");
              if (pch0 != 0){
                streamlog_out ( DEBUG5 ) << " Parsing pede output for final chi2/ndf result.. " << std::endl;
                // search for the equal sign after which the result for chi2/ndf is stated within the next 80 chars 
                // (with offset of 22 chars since pch points to beginning of "Sum(..." string just found)
                char* pch = (char*)(memchr (pch0+22, '=', 180));
                if (pch!=NULL){
                  char str[16];
                  // now copy the numbers after the equal sign
                  strncpy ( str, pch+1, 15 );
                  str[15] = '\0';   /* null character manually added */
                  // monitor the chi2/ndf in CDash when running tests
//                  CDashMeasurement meas_chi2ndf("chi2_ndf",atof(str));  cout << meas_chi2ndf; // output only if DO_TESTING is set
                  streamlog_out ( MESSAGE6 ) << "Final Sum(Chi^2)/Sum(Ndf) = " << str << std::endl;
                }	    
              }
            }

            // wait for the pede execution to finish
            pede.close();

            // check the exit value of pede / react to previous errors
            if ( pede.rdbuf()->status() == 0 && !encounteredError) 
            {
              streamlog_out ( MESSAGE7 ) << "Pede successfully finished" << std::endl;
            } else {
              streamlog_out ( ERROR5 ) << "Problem during Pede execution, exit status: " << pede.rdbuf()->status() << ", error messages (repeated here): " << std::endl;
              streamlog_out ( ERROR5 ) << pedeerrors.str() << std::endl;
              // TODO: decide what to do now; exit? and if, how?
              streamlog_out ( ERROR5 ) << "Will exit now" << std::endl;
              //exit(EXIT_FAILURE); // FIXME: can lead to (ROOT?) seg faults - points to corrupt memory? run valgrind...
              return; // does fine for now
            }

            // reading back the millepede.res file and getting the
            // results.
            std::string millepedeResFileName = "millepede.res";

            streamlog_out ( MESSAGE6 ) << "Reading back the " << millepedeResFileName << std::endl
                                       << "Saving the alignment constant into " << _alignmentConstantLCIOFile << std::endl;

            // open the millepede ASCII output file
            ifstream millepede( millepedeResFileName.c_str() );


            // reopen the LCIO file this time in append mode
            LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

            try 
            {
              lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
            }
            catch ( IOException& e ) 
            {
              streamlog_out ( ERROR4 ) << e.what() << std::endl
                                       << "Sorry for quitting. " << std::endl;
              exit(-1);
            }


            // write an almost empty run header
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


            if ( millepede.bad() || !millepede.is_open() ) 
            {
              streamlog_out ( ERROR4 ) << "Error opening the " << millepedeResFileName << std::endl
                                       << "The alignment slcio file cannot be saved" << std::endl;
            }
            else 
            {
              std::vector<double > tokens;
              std::stringstream tokenizer;
              std::string line;
              double buffer;

              // get the first line and throw it away since it is a
              // comment!
              getline( millepede, line );

              int counter = 0;

              while ( ! millepede.eof() ) {

                EUTelAlignmentConstant * constant = new EUTelAlignmentConstant;

                bool goodLine = true;
                unsigned int numpars = 0;
                if(_alignMode != 3)
                  numpars = 3;
                else
                  numpars = 6;

                bool _nonzero_tokens = false;

                for ( unsigned int iParam = 0 ; iParam < numpars ; ++iParam ) 
                {
                  getline( millepede, line );

                  if ( line.empty() ) {
                    goodLine = false;
                    continue;
                  }

                  tokens.clear();
                  tokenizer.clear();
                  tokenizer.str( line );

                  // check that all parts of the line are non zero
                  while ( tokenizer >> buffer ) {
                    tokens.push_back( buffer ) ;
                    if(buffer> 1e-12) _nonzero_tokens = true;
                  }

                  if ( ( tokens.size() == 3 ) || ( tokens.size() == 6 ) || (tokens.size() == 5) ) {
                    goodLine = true;
                  } else goodLine = false;

                  bool isFixed = ( tokens.size() == 3 );
                  if(_alignMode != 3 && _alignMode == 4) //GBL case
                    {
                     if ( iParam == 0 ) {
                        constant->setXOffset( tokens[1] );
			streamlog_out(DEBUG5)<<"X offset:"<<std::endl;
			streamlog_out(DEBUG5)<<tokens[1]<<std::endl;
                        if ( ! isFixed ) constant->setXOffsetError( tokens[4] ) ;
                      }
                      if ( iParam == 1 ) {
			streamlog_out(DEBUG5)<<"Y offset:"<<std::endl;
			streamlog_out(DEBUG5)<<tokens[1]<<std::endl;
                        constant->setYOffset( tokens[1] ) ;
                        if ( ! isFixed ) constant->setYOffsetError( tokens[4] ) ;
                      }
                      if ( iParam == 2 ) {
			streamlog_out(DEBUG5)<<"Gamma:"<<std::endl;
			streamlog_out(DEBUG5)<<tokens[1]<<std::endl;
                        constant->setGamma( tokens[1]  ) ;
                        if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
                      }
                    }
                  else
                    {
                     if ( iParam == 0 ) {
                        constant->setXOffset( tokens[1] / 1000 );
                        if ( ! isFixed ) constant->setXOffsetError( tokens[4] / 1000 ) ;                    
                      }
                      if ( iParam == 1 ) {
                        constant->setYOffset( tokens[1] / 1000 ) ;
                        if ( ! isFixed ) constant->setYOffsetError( tokens[4] / 1000 ) ;
                      }
                      if ( iParam == 2 ) {
                        constant->setZOffset( tokens[1] / 1000 ) ;
                        if ( ! isFixed ) constant->setZOffsetError( tokens[4] / 1000 ) ;
                      }
                      if ( iParam == 3 ) {
                        constant->setAlpha( tokens[1]  ) ;
                        if ( ! isFixed ) constant->setAlphaError( tokens[4] ) ;
                      } 
                      if ( iParam == 4 ) {
                        constant->setBeta( tokens[1]  ) ;
                        if ( ! isFixed ) constant->setBetaError( tokens[4] ) ;
                      } 
                      if ( iParam == 5 ) {
                        constant->setGamma( tokens[1]  ) ;
                        if ( ! isFixed ) constant->setGammaError( tokens[4] ) ;
                      } 

                    }

                }


                // right place to add the constant to the collection
                if ( goodLine  ) {
    //               constant->setSensorID( _orderedSensorID_wo_excluded.at( counter ) );
                  constant->setSensorID( _orderedSensorID.at( counter ) );
                  ++ counter;
                  constantsCollection->push_back( constant );
                  streamlog_out ( MESSAGE0 ) << (*constant) << std::endl;
                }
                else delete constant;
              }

            }



            event->addCollection( constantsCollection, _alignmentConstantCollectionName );
            lcWriter->writeEvent( event );
            delete event;

            lcWriter->close();

            millepede.close();

          }
        } else {

          streamlog_out ( ERROR2 ) << "Unable to run pede. No steering file has been generated." << std::endl;

        }
    } // end if running pede using the generated steering file

    streamlog_out(MESSAGE2) << std::endl;
    streamlog_out(MESSAGE2) << "Successfully finished" << std::endl;
}

void EUTelMilleGBL::bookHistos() {


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    try {
        streamlog_out(MESSAGE2) << "Booking histograms..." << std::endl;

        const int tracksNBin = 20;
        const double tracksMin = -0.5;
        const double tracksMax = 19.5;

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
        

        // Track candidate search
        AIDA::IHistogram1D * numberTracksCandidates =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_numberTracksCandidatesHistName, tracksNBin, tracksMin, tracksMax);
        if (numberTracksCandidates) {
            numberTracksCandidates->setTitle("Number of track candidates;N tracks;N events");
            _aidaHistoMap1D.insert(std::make_pair(_numberTracksCandidatesHistName, numberTracksCandidates));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_numberTracksCandidatesHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            _histogramSwitch = false;
        }

        // GBL fits
        AIDA::IHistogram1D * chi2GblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_chi2GblFitHistName, chi2NBin, chi2Min, chi2Max);
        if (chi2GblFit) {
            chi2GblFit->setTitle("#chi^{2} of track candidates; #chi^{2};N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_chi2GblFitHistName, chi2GblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_chi2GblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            _histogramSwitch = false;
        }

        AIDA::IHistogram1D * probGblFit =
                marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(_probGblFitHistName, probNBin, probMin, probMax);
        if (probGblFit) {
            probGblFit->setTitle("Probability of track fit; Prob;N Tracks");
            _aidaHistoMap1D.insert(std::make_pair(_probGblFitHistName, probGblFit));
        } else {
            streamlog_out(ERROR2) << "Problem booking the " << (_probGblFitHistName) << std::endl;
            streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
            _histogramSwitch = false;
        }

        // Residuals after fit
        std::stringstream sstm;
        std::string residGblFitHistName;
        std::string histTitle;
        for (int iPlane = 0; iPlane < 6; iPlane++) {
            sstm << _residGblFitHistNameX << iPlane;
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << iPlane << "X direction; r; N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
        }

        for (int iPlane = 0; iPlane < 6; iPlane++) {
            sstm << _residGblFitHistNameY << iPlane;
            residGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << iPlane << "Y direction; r; N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(residGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
        }
        
        // 2D histograms
        NBinX = 100;
        NBinY = 100;
        MinX = -0.001;
        MaxX = 0.001;
        std::string resid2DGblFitHistName;
        for (int iPlane = 0; iPlane < 6; iPlane++) {
            sstm << _resid2DGblFitHistNameXvsX << iPlane;
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << iPlane << "; x (mm); rx";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit1) {
                residGblFit1->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit1));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
            
            sstm << _resid2DGblFitHistNameXvsY << iPlane;
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << iPlane << "; y (mm); rx";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit2) {
                residGblFit2->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit2));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
        }

        for (int iPlane = 0; iPlane < 6; iPlane++) {
            sstm << _resid2DGblFitHistNameYvsX << iPlane;
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << iPlane << "; x (mm); ry";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram2D * residGblFit1 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit1) {
                residGblFit1->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit1));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
            
            sstm << _resid2DGblFitHistNameYvsY << iPlane;
            resid2DGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Normalised residuals. Plane " << iPlane << "; y (mm); ry";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram2D * residGblFit2 =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D(resid2DGblFitHistName, NBinX, MinX, MaxX, NBinY, MinY, MaxY);
            if (residGblFit2) {
                residGblFit2->setTitle(histTitle);
                _aidaHistoMap2D.insert(std::make_pair(resid2DGblFitHistName, residGblFit2));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (resid2DGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
        }
        
        // Kink angles after fit
        MinX = -0.001;
        MaxX = 0.001;
        std::string kinkGblFitHistName;
        for (int iPlane = 0; iPlane < 6; iPlane++) {
            sstm << _kinkGblFitHistNameX << iPlane;
            kinkGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Kink angles. Plane " << iPlane << "X direction; kink (rad); N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(kinkGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(kinkGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (kinkGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
        }

        for (int iPlane = 0; iPlane < 6; iPlane++) {
            sstm << _kinkGblFitHistNameY << iPlane;
            kinkGblFitHistName = sstm.str();
            sstm.str(std::string());
            sstm << "Kink angles. Plane " << iPlane << "Y direction; kink (rad); N hits";
            histTitle = sstm.str();
            sstm.str(std::string(""));
            AIDA::IHistogram1D * residGblFit =
                    marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(kinkGblFitHistName, NBinX, MinX, MaxX);
            if (residGblFit) {
                residGblFit->setTitle(histTitle);
                _aidaHistoMap1D.insert(std::make_pair(kinkGblFitHistName, residGblFit));
            } else {
                streamlog_out(ERROR2) << "Problem booking the " << (kinkGblFitHistName) << std::endl;
                streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
                _histogramSwitch = false;
            }
            sstm.str(std::string(""));
        }

    } catch (lcio::Exception& e) {
#ifdef EUTEL_INTERACTIVE
        streamlog_out(ERROR2) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << std::endl;
        string answer;
        while (true) {
            streamlog_out(ERROR2) << "[q]/[c]" << std::endl;
            cin >> answer;
            transform(answer.begin(), answer.end(), answer.begin(), ::tolower);
            if (answer == "q") {
                exit(-1);
            } else if (answer == "c")
                _histogramSwitch = false;
            break;
        }
#else
        streamlog_out(WARNING2) << "No AIDAProcessor initialized. Continue without histogramming" << std::endl;
#endif
    }
#endif
}
#endif // USE_GEAR

