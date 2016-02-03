// Version: $Id$
// Author Havard Gjersdal, UiO(haavagj@fys.uio.no)
/*!
 * This is a track fitting processor for the Eutelescope package. 
 *
 * It preforms track finding and fitting on a supplied hit collection.
 *
 * The track finder works by propagating all hits to plane 0, currently assuming straight
 * line fits, then running a cluster finder. Hit clusters above some set value are considered
 * track candidates.
 *
 * This track candidate is then fitted using a implementation of a Deterministic Annealing
 * Filter (DAF), that in short is a Kalman Filter running iteratively over a set of weighted
 * measurements, reweighing the measurements after each fit based on the residuals and a
 * supplied chi2 cut off.
 *
 * This package uses the Eigen library for linear algebra. This package is very quick when
 * compiled properly, but very slow when compiled for debugging. Make sure to compile
 * properly before running productions.
 *
 * Running 'cmake -i' inside the build folder, and then when it asks
 * Variable Name: CMAKE_CXX_FLAGS_RELEASE
 * Description: Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).                          
 *
 * enter:
 * New Value (Enter to keep current value): -O3 -msse2 -ftree-vectorize -DNDEBUG
 *
 * When it asks
 * Variable Name: CMAKE_BUILD_TYPE
 * enter:
 * New Value (Enter to keep current value): Release
 *
 * If youc cpu supports it, you could try -msse4 or -msse3 aswell.
 */

// built only if GEAR and MARLINUTIL are used
#if defined(USE_GEAR)
// eutelescope includes ".h"
#include "EUTelDafAlign.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelVirtualCluster.h"
#include "EUTelExceptions.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelPStream.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "marlin/AIDAProcessor.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#endif

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


EUTelDafAlign::EUTelDafAlign ()
: EUTelDafBase("EUTelDafAlign"),
  _runPede(false), 
  _pedeSteerfileName(""),
  _binaryFilename(""),
  _alignmentConstantLCIOFile(""),
  _alignmentConstantCollectionName(""),
  _translate(),
  _translateX(),
  _translateY(),
  _zRot(),
  _scale(),
  _scaleX(),
  _scaleY(),
  _resXMin(),
  _resXMax(),
  _resYMin(),
  _resYMax(),
  _mille(NULL),
  _resX(),
  _resY(),
  _dutMatches()
{
    //Child spesific params and description
  dafParams();
}

void EUTelDafAlign::dafParams(){
  _description = "This processor preforms track reconstruction. The tracks are used for alignment of a partially aligned system.";

  //Millepede options
  registerOptionalParameter("RunPede","Build steering file, binary input file, and execute the pede program.",_runPede, static_cast <bool> (true));
  registerOptionalParameter("PedeSteerfileName","Name of the steering file for the pede program.",_pedeSteerfileName, string("steer_mille.txt"));
  registerOptionalParameter("BinaryFilename","Name of binary input file for Millepede.",_binaryFilename, string ("mille.bin"));
  registerOptionalParameter("AlignmentConstantLCIOFile","Name of LCIO db file where alignment constantds will be stored", 
			    _alignmentConstantLCIOFile, std::string( "alignment.slcio" ) );
  registerOptionalParameter("AlignmentConstantCollectionName", "This is the name of the alignment collection to be saved into the slcio file",
                            _alignmentConstantCollectionName, std::string( "alignment" ));
  //Alignment parameter options
  registerOptionalParameter("Translate", "List of sensor IDs to where translations should be free", _translate, std::vector<int>());
  registerOptionalParameter("TranslateX", "List of sensor IDs to where translations should be free in the x direction", _translateX, std::vector<int>());
  registerOptionalParameter("TranslateY", "List of sensor IDs to where translations should be free in the y direction", _translateY, std::vector<int>());
  registerOptionalParameter("ZRotate", "List of sensor IDs to where rotation over z-axis should be free", _zRot, std::vector<int>());
  registerOptionalParameter("Scale", "List of sensor IDs to where scales of X and Y axis should be free", _scale, std::vector<int>());
  registerOptionalParameter("ScaleX", "List of sensor IDs to where scales of X axis should be free", _scaleX, std::vector<int>());
  registerOptionalParameter("ScaleY", "List of sensor IDs to where scales of Y axis should be free", _scaleY, std::vector<int>());
  //Track cuts
  registerOptionalParameter("ResidualXMin", "Min X for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resXMin, std::vector<float>());
  registerOptionalParameter("ResidualXMax", "Max X for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resXMax, std::vector<float>());
  registerOptionalParameter("ResidualYMin", "Min Y for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resYMin, std::vector<float>());
  registerOptionalParameter("ResidualYMax", "Max Y for residual cuts. All sensors shoule be included, same order as parameter DutPlanes",
			    _resYMax, std::vector<float>());
  //Region of cell
  registerOptionalParameter("MinCol", "Minimum allowed local X of hit. All clusters containing hits outside will be discarded from alignment.", _colMin, std::vector<int>());
  registerOptionalParameter("MaxCol", "Maximum allowed local X of hit. All clusters containing hits outside will be discarded from alignment.", _colMax, std::vector<int>());
  registerOptionalParameter("MinRow", "Minimum allowed local Y of hit. All clusters containing hits outside will be discarded from alignment.", _rowMin, std::vector<int>());
  registerOptionalParameter("MaxRow", "Maximum allowed local Y of hit. All clusters containing hits outside will be discarded from alignment.", _rowMax, std::vector<int>());
}

void EUTelDafAlign::dafInit() {
  for(size_t ii = 0; ii < _translate.size(); ii++){
    _translateX.push_back( _translate.at(ii) );
    _translateY.push_back( _translate.at(ii) );
  }
  for(size_t ii = 0; ii < _scale.size(); ii++){
    _scaleX.push_back( _scale.at(ii) );
    _scaleY.push_back( _scale.at(ii) );
  }
  for(size_t ii = 0; ii < _translateX.size(); ii++){
    
  }
  if(_runPede){
    _mille = new Mille(_binaryFilename.c_str());
    streamlog_out ( MESSAGE5 ) << "The filename for the mille binary file is: " << _binaryFilename.c_str() << endl;
  }
  
  for(size_t ii = 0; ii < _dutPlanes.size(); ii++){
    int iden = _dutPlanes.at(ii);
    float xMin;
    if(_resXMin.size() > ii){
      xMin = _resXMin.at(ii);
    }  else{
      xMin = -9999999;
    }
    float xMax;
    if(_resXMax.size() > ii){
      xMax = _resXMax.at(ii);
    }  else{
      xMax = 9999999;
    }
    float yMin;
    if(_resYMin.size() > ii){
      yMin = _resYMin.at(ii);
    }  else{
      yMin = -9999999;
    }
    float yMax;
    if(_resYMax.size() > ii){
      yMax = _resYMax.at(ii);
    }  else{
      yMax = 9999999;
    }
    streamlog_out ( MESSAGE5 ) << "xMin " << xMin << " " << xMax << endl;
    streamlog_out ( MESSAGE5 ) << "yMin " << yMin << " " << yMax << endl;
    _resX[iden] = make_pair(xMin, xMax);
    _resY[iden] = make_pair(yMin, yMax);
  }
}

int EUTelDafAlign::checkDutResids(daffitter::TrackCandidate<float,4>& track){
  int nHits(0);
  _dutMatches.clear();
  for( size_t ii = 0; ii < _system.planes.size() ; ii++){
    daffitter::FitPlane<float>& plane = _system.planes.at(ii);
    int iden = plane.getSensorID();
    if( find(_dutPlanes.begin(), _dutPlanes.end(), iden) == _dutPlanes.end()){ continue; }
    std::pair<float, float> resX = _resX[iden];
    std::pair<float, float> resY = _resY[iden];

    for(size_t w = 0; w < plane.meas.size(); w++){
      float measWeight = 1.0f;
      daffitter::Measurement<float>& meas = plane.meas.at(w);
      daffitter::TrackEstimate<float,4>& estim = track.estimates.at(ii);
      //Resids 
      if( (estim.getX() - meas.getX()) < resX.first) { measWeight = 0.0f;}
      if( (estim.getX() - meas.getX()) > resX.second){ measWeight = 0.0f;}
      if( (estim.getY() - meas.getY()) < resY.first) { measWeight = 0.0f;}
      if( (estim.getY() - meas.getY()) > resY.second){ measWeight = 0.0f;}
      if(measWeight > 0.5){ 
	nHits++; 
	if( not meas.goodRegion() ){  measWeight = 0.0f;  }
      }
      track.weights.at(ii)(w) = measWeight;
      if( measWeight > 0.5 ){
	_dutMatches.push_back(ii);
      }
    }
  }
  return(nHits);
}

void EUTelDafAlign::dafEvent (LCEvent* /*event*/) {
  //Check found tracks
  for(size_t ii = 0; ii < _system.getNtracks(); ii++ ){
    //run track fitter
    _nCandidates++;
    _system.fitPlanesInfoDaf(_system.tracks.at(ii));
    //Check resids, intime, angles
    if(not checkTrack( _system.tracks.at(ii))) { continue;};
    //This guy includes DUT planes and adds weights to measurements based on resid cuts
    if( checkDutResids ( _system.tracks.at(ii)) >= _nDutHits ) {
      //Fill plots
      if(_histogramSwitch){ 
	fillPlots( _system.tracks.at(ii) ); 
	fillDetailPlots( _system.tracks.at(ii) ); 
      }
      _nTracks++;
      if(_runPede){
	for(size_t dut = 0; dut < _dutMatches.size(); dut++){
	  _system.planes.at( _dutMatches.at(dut) ).include();
	  _system.weightToIndex(_system.tracks.at(ii));
	  _system.fitPlanesInfoBiased(_system.tracks.at(ii));
	  //Add to mille bin file
	  addToMille( _system.tracks.at(ii));
	  _system.planes.at( _dutMatches.at(dut) ).exclude();
	}
      }
    }
  }
}

void EUTelDafAlign::addToMille(daffitter::TrackCandidate<float,4>& track){
  const int nLC = 4; //number of local parameters
  const int nGL = _system.planes.size() * 5; // number of global parameters
  
  float *derLC = new float[nLC]; // array of derivatives for local parameters
  float *derGL = new float[nGL]; // array of derivatives for global parameters
  int *label = new int[nGL]; // array of parameter labels
  
  for(int ii = 0; ii < nGL; ii++){ 
    label[ii] = ii + 1;
    derGL[ii] = 0;
  }
  for(int ii = 0; ii < nLC; ii++){ derLC[ii] = 0; }

  for(size_t ii = 0; ii < _system.planes.size(); ii++){
    daffitter::FitPlane<float>& pl = _system.planes.at(ii);
    if(pl.isExcluded() ) { continue; }
    int index = track.indexes.at(ii);
    //index < 0 means plane is excluded
    if( index < 0) { continue; }
    daffitter::Measurement<float>& meas = pl.meas.at(index);
    daffitter::TrackEstimate<float,4>& estim = track.estimates.at(ii);
    derGL[(ii * 5)    ] = -1; //Derivatives of residuals w.r.t. shift in x
    derGL[(ii * 5) + 2] = meas.getY(); //Derivatives of residuals w.r.t. z rotations
    derGL[(ii * 5) + 3] = meas.getX(); //Derivatives of residuals w.r.t. scale of x axis
    derLC[0] = 1; //Derivatives of fit pos w.r.t. x
    derLC[2] = pl.getMeasZ(); //Derivatives of fit pos w.r.t. dx/dz
    _mille->mille(nLC, derLC, nGL, derGL, label, estim.getX()- meas.getX(), pl.getSigmaX());

    derGL[(ii * 5)]     = 0;
    derGL[(ii * 5) + 2] = 0;
    derGL[(ii * 5) + 3] = 0;
    derLC[0] = 0;
    derLC[2] = 0;

    derGL[(ii * 5) + 1] = -1; //Derivatives of residuals w.r.t. shift in y
    derGL[(ii * 5) + 2] = -1 * meas.getX(); //Derivatives of residuals w.r.t. z rotations
    derGL[(ii * 5) + 4] = meas.getY();//Derivatives of residuals w.r.t. scales of y axis

    derLC[1] = 1; //Derivatives of fit pos w.r.t. y
    derLC[3] = pl.getMeasZ(); //Derivatives of fit pos w.r.t. dy/dz
    
    _mille->mille(nLC, derLC, nGL, derGL, label, estim.getY() - meas.getY(), pl.getSigmaY());
    
    derGL[(ii * 5) + 1] = 0;
    derGL[(ii * 5) + 2] = 0;
    derGL[(ii * 5) + 4] = 0;
    derLC[1] = 0;
    derLC[3] = 0;
  }

  delete [] derLC;
  delete [] derGL;
  delete [] label;

  _mille->end();
}


void EUTelDafAlign::generatePedeSteeringFile(){
  ofstream steerFile;
  steerFile.open(_pedeSteerfileName.c_str());
  if (not steerFile.is_open()) {
    streamlog_out ( ERROR2 ) << "Unable to open pede steering file, quitting" << endl;
    throw runtime_error("Unable to open file " + _pedeSteerfileName);
  }
  steerFile << "Cfiles" << endl;
  steerFile << _binaryFilename << endl;
  steerFile << endl;
  steerFile << "Parameter" << endl;
  for(size_t ii = 0; ii < _system.planes.size(); ii++){
    daffitter::FitPlane<float>& pl = _system.planes.at(ii);
    int iden = pl.getSensorID();
    steerLine(steerFile, (ii * 5) + 1, iden, _translateX);
    steerLine(steerFile, (ii * 5) + 2, iden, _translateY);
    steerLine(steerFile, (ii * 5) + 3, iden, _zRot);
    steerLine(steerFile, (ii * 5) + 4, iden, _scaleX);
    steerLine(steerFile, (ii * 5) + 5, iden, _scaleY);
  }
  steerFile << endl;
  steerFile << "! chiscut 5.0 2.5" << endl;
  steerFile << "! outlierdownweighting 4" << endl;
  steerFile << endl;
  steerFile << "method inversion 10 0.001" << endl;
  steerFile << endl;
  steerFile << "histprint" << endl;
  steerFile << endl;
  steerFile << "end" << endl;
  steerFile.close();
  streamlog_out ( MESSAGE5 ) << "File " << _pedeSteerfileName << " written." << endl;
}

void EUTelDafAlign::steerLine(ofstream &steerFile, int label, int iden, std::vector<int> idens){
  if( find(idens.begin(), idens.end(), iden) != idens.end()){
    steerFile << label << " 0.0 0.0" << endl; 
    streamlog_out( MESSAGE5 ) << "Steering line: " << label << " 0.0 0.0" << endl; 
  } else {
    steerFile << label << " 0.0 -1.0" << endl; 
  }
}

void EUTelDafAlign::runPede(){
  std::string command = "pede " + _pedeSteerfileName;
  // create a new process
  redi::ipstream which("which pede");
  // wait for the process to finish
  which.close();
 
  if (  which.rdbuf()->status() == 255 ) {
    streamlog_out( ERROR5 ) << "Cannot find pede program in path. Nothing to do." << endl;
    throw runtime_error("Cannot find pede in path.");
  }
  streamlog_out ( MESSAGE5 ) << "Starting pede..." << endl;
  redi::ipstream pede( command.c_str() );
  string output;
  while ( getline( pede, output ) ) { streamlog_out( MESSAGE5 ) << output << endl; }
  // wait for the pede execution to finish
  pede.close();
  // check the exit value of pede
  if ( pede.rdbuf()->status() == 0 ) {
    streamlog_out ( MESSAGE5 ) << "Pede successfully finished" << endl;
  } else {
    throw runtime_error("Pede exitted abmormally.");
  }
  // reading back the millepede.res file and getting the results. 
  string millepedeResFileName = "millepede.res";

  streamlog_out ( MESSAGE5 ) << "Reading back " << millepedeResFileName << endl
			     << "Saving the alignment constant into " << _alignmentConstantLCIOFile << endl;

  // open the millepede ASCII output file
  ifstream millepede( millepedeResFileName.c_str() );

  if(not millepede.is_open() ){
    throw runtime_error("Unable to open " + millepedeResFileName + ".");
  }
  //Open the alignment db file
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
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
  
  vector<double > values;
  stringstream linestream;
  string line;
  double value;

  // get the first line and throw it away since it is a comment!
  getline( millepede, line );
  while ( not millepede.eof() ) {
    EUTelAlignmentConstant* constant = NULL;
    getline( millepede, line );
   
    values.clear(); linestream.clear(); linestream.str( line );
    
    while ( linestream >> value ) { values.push_back( value ); }
    
    int nValues = values.size();
    
    if(( nValues != 3 ) and ( nValues != 5) and (nValues != 6)){ continue; }
    bool isFixed = ( nValues == 3 );
    int label = int(values.at(0) - 1);//Possible loss of data...
    int plane = label / 5;
    switch( label % 5){
    case 0:
      //First line, we need a new alignment constant thingy.
      constant = new EUTelAlignmentConstant();
      constant->setXOffset ( values.at(1) / 1000.0);
      if( not isFixed){  constant->setXOffsetError( values.at(4)/ 1000.0);}
      break;
    case 1:
      constant->setYOffset ( values.at(1) / 1000.0);
      if( not isFixed){  constant->setYOffsetError( values.at(4)/ 1000.0);}
      break;
    case 2:
      constant->setGamma( values.at(1));
      if( not isFixed){  constant->setGammaError( values.at(4));}
      break;
    case 3:
      constant->setAlpha( values.at(1));
      if( not isFixed){  constant->setAlphaError( values.at(4));}
      break;
    case 4:
      constant->setBeta( values.at(1));
      if( not isFixed){  constant->setBetaError( values.at(4));}
      //Last line, add constant to collection
      constant->setSensorID( _system.planes.at(plane).getSensorID() );
      constantsCollection->push_back( constant );
      streamlog_out ( MESSAGE5 ) << (*constant) << endl;
      break;
    default:
      throw runtime_error("Error parsing millepede.res");
    }
  }
  event->addCollection( constantsCollection, _alignmentConstantCollectionName );
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
  millepede.close();
}

void EUTelDafAlign::dafEnd() {
  if(not _runPede){ return; }
  delete _mille;
  generatePedeSteeringFile();
  runPede();
}
#endif // USE_GEAR
