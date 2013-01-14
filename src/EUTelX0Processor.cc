#include "EUTelX0Processor.h"
#include <cmath>
using namespace marlin;
using namespace eutelescope;
using namespace std;

std::string EUTelX0Processor::_histoResidualX = "ResidualX";
std::string EUTelX0Processor::_histoResidualXPlane1 = "ResidualXPlane1";
std::string EUTelX0Processor::_histoResidualXPlane2 = "ResidualXPlane2";
std::string EUTelX0Processor::_histoResidualXPlane3 = "ResidualXPlane3";
std::string EUTelX0Processor::_histoResidualXPlane4 = "ResidualXPlane4";
std::string EUTelX0Processor::_histoResidualYPlane1 = "ResidualYPlane1";
std::string EUTelX0Processor::_histoResidualYPlane2 = "ResidualYPlane2";
std::string EUTelX0Processor::_histoResidualYPlane3 = "ResidualYPlane3";
std::string EUTelX0Processor::_histoResidualYPlane4 = "ResidualYPlane4";
std::string EUTelX0Processor::_histoResidualY = "ResidualY";
std::string EUTelX0Processor::_histoResidualXZ = "ResidualXZ";
std::string EUTelX0Processor::_histoResidualYZ = "ResidualYZ";
std::string EUTelX0Processor::_histoResidualXY = "ResidualYZ";

EUTelX0Processor::EUTelX0Processor():Processor("EUTelX0Processor")
{
}

void EUTelX0Processor::init()
{
  _debug = false; 
  int nobins = 100, nobinsangle = 100;//Number of bins in the histograms
  double minbin = -0.4, maxbin = 0.4;//Maximum and minimum bin values
  double minbinangle = -25, maxbinangle = 25, minbinalpha = -5, maxbinalpha = 5;
  std::vector<double> empty;  
  
  AIDA::IHistogram2D * ResidualXZ = AIDAProcessor::histogramFactory(this)->createHistogram2D(_histoResidualXZ, nobins, minbin, maxbin, nobins, minbin, maxbin );//Creates a profile histogram which is going to store the residuals as a function of X and Y
  ResidualXZ->setTitle(_histoResidualXZ);//Set the title for this histogram
  _histoThing.insert(make_pair(_histoResidualXZ, ResidualXZ));//Insert the histogram into the map
  
  AIDA::IHistogram2D * ResidualYZ = AIDAProcessor::histogramFactory(this)->createHistogram2D(_histoResidualYZ , nobins, minbin, maxbin, nobins, minbin, maxbin );//Creates a profile histogram which is going to store the residuals as a function of X and Y
  ResidualYZ->setTitle(_histoResidualYZ );//Set the title for this histogram
  _histoThing.insert(make_pair(_histoResidualYZ , ResidualYZ));//Insert the histogram into the map
  
  AIDA::IHistogram2D * ResidualXY = AIDAProcessor::histogramFactory(this)->createHistogram2D(_histoResidualXY , nobins, minbin, maxbin, nobins, minbin, maxbin );//Creates a profile histogram which is going to store the residuals as a function of X and Y
  ResidualXY->setTitle(_histoResidualXY);//Set the title for this histogram
  _histoThing.insert(make_pair(_histoResidualXY , ResidualXY));//Insert the histogram into the map
  
  AIDA::IHistogram2D * ResidualXwrtX = AIDAProcessor::histogramFactory(this)->createHistogram2D("ResidualXwrtX",nobins,minbin,maxbin,nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualXwrtX->setTitle("ResidualXwrtX");
  _histoThing.insert(make_pair("ResidualXwrtX",ResidualXwrtX));
  
  AIDA::IHistogram2D * ResidualXwrtY = AIDAProcessor::histogramFactory(this)->createHistogram2D("ResidualXwrtY",nobins,minbin,maxbin,nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualXwrtY->setTitle("ResidualXwrtY");
  _histoThing.insert(make_pair("ResidualXwrtY",ResidualXwrtY));

  AIDA::IHistogram2D * ResidualYwrtX = AIDAProcessor::histogramFactory(this)->createHistogram2D("ResidualYwrtX",nobins,minbin,maxbin,nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualYwrtX->setTitle("ResidualYwrtX");
  _histoThing.insert(make_pair("ResidualYwrtX",ResidualYwrtX));
  
  AIDA::IHistogram2D * ResidualYwrtY = AIDAProcessor::histogramFactory(this)->createHistogram2D("ResidualYwrtY",nobins,minbin,maxbin,nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualYwrtY->setTitle("ResidualYwrtY");
  _histoThing.insert(make_pair("ResidualYwrtY",ResidualYwrtY));

  AIDA::IHistogram1D * ResidualX = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualX",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualX->setTitle("ResidualX");
  _histoThing.insert(make_pair("ResidualX",ResidualX));
  _histoData["ResidualX"] = empty;
 
  AIDA::IHistogram1D * ResidualXPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualXPlane1",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualXPlane1->setTitle("ResidualXPlane1");
  _histoThing.insert(make_pair("ResidualXPlane1",ResidualXPlane1));
  _histoData["ResidualXPlane1"] = empty;
 
  AIDA::IHistogram1D * ResidualXPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualXPlane2",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualXPlane2->setTitle("ResidualXPlane2");
  _histoThing.insert(make_pair("ResidualXPlane2",ResidualXPlane2));
  _histoData["ResidualXPlane2"] = empty;
 
  AIDA::IHistogram1D * ResidualXPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualXPlane3",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualXPlane3->setTitle("ResidualXPlane3");
  _histoThing.insert(make_pair("ResidualXPlane3",ResidualXPlane3));
  _histoData["ResidualXPlane3"] = empty;
 
  AIDA::IHistogram1D * ResidualXPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualXPlane4",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualXPlane4->setTitle("ResidualXPlane4");
  _histoThing.insert(make_pair("ResidualXPlane4",ResidualXPlane4));
  _histoData["ResidualXPlane4"] = empty;
 
  AIDA::IHistogram1D * ResidualYPlane1 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualYPlane1",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualYPlane1->setTitle("ResidualYPlane1");
  _histoThing.insert(make_pair("ResidualYPlane1",ResidualYPlane1));
  _histoData["ResidualYPlane1"] = empty;
 
  AIDA::IHistogram1D * ResidualYPlane2 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualYPlane2",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualYPlane2->setTitle("ResidualYPlane2");
  _histoThing.insert(make_pair("ResidualYPlane2",ResidualYPlane2));
  _histoData["ResidualYPlane2"] = empty;
 
  AIDA::IHistogram1D * ResidualYPlane3 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualYPlane3",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualYPlane3->setTitle("ResidualYPlane3");
  _histoThing.insert(make_pair("ResidualYPlane3",ResidualYPlane3));
  _histoData["ResidualYPlane3"] = empty;
 
  AIDA::IHistogram1D * ResidualYPlane4 = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualYPlane4",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualYPlane4->setTitle("ResidualYPlane4");
  _histoThing.insert(make_pair("ResidualYPlane4",ResidualYPlane4));
  _histoData["ResidualYPlane4"] = empty;
  
  AIDA::IHistogram1D * ResidualY = AIDAProcessor::histogramFactory(this)->createHistogram1D("ResidualY",nobins,minbin,maxbin);//Create a histogram for the residual
  ResidualY->setTitle("ResidualY");
  _histoThing.insert(make_pair("ResidualY",ResidualY));
  _histoData["ResidualY"] = empty;
 
  AIDA::IHistogram1D * Theta012 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Theta012",nobinsangle,minbinangle,maxbinangle);//Create a histogram for the residual
  Theta012->setTitle("Theta012");
  _histoThing.insert(make_pair("Theta012",Theta012));
  _histoData["Theta012"] = empty;
 
  
  AIDA::IHistogram1D * Phi012 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Phi012",nobinsangle,minbinangle,maxbinangle);//Create a histogram for the residual
  Phi012->setTitle("Phi012");
  _histoThing.insert(make_pair("Phi012",Phi012));
  _histoData["Phi012"] = empty;
 
  
  AIDA::IHistogram1D * Alpha012 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Alpha012",nobinsangle,minbinalpha,maxbinalpha);//Create a histogram for the residual
  Alpha012->setTitle("Alpha012");
  _histoThing.insert(make_pair("Alpha012",Alpha012));
  _histoData["Alpha012"] = empty;
 
  
  AIDA::IHistogram1D * Theta432 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Theta432",nobinsangle,minbinangle,maxbinangle);//Create a histogram for the residual
  Theta432->setTitle("Theta432");
  _histoThing.insert(make_pair("Theta432",Theta432));
  _histoData["Theta432"] = empty;
 
  
  AIDA::IHistogram1D * Phi432 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Phi432",nobinsangle,minbinangle,maxbinangle);//Create a histogram for the residual
  Phi432->setTitle("Phi432");
  _histoThing.insert(make_pair("Phi432",Phi432));
  _histoData["Phi432"] = empty;
 
  
  AIDA::IHistogram1D * Alpha432 = AIDAProcessor::histogramFactory(this)->createHistogram1D("Alpha432",nobinsangle,minbinalpha,maxbinalpha);//Create a histogram for the residual
  Alpha432->setTitle("Alpha432");
  _histoThing.insert(make_pair("Alpha432",Alpha432));
  _histoData["Alpha432"] = empty;
 

  _inputCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);//Used to store the values of the hit events
  // readin parameter for input hit collection name (default alignedHit)
  //registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName", "The name of the input hit collection", _inputHitCollectionName, string("correctedHit"));//Used to store the values of the hit events
  
  registerInputCollection(LCIO::TRACKERHIT,"AlignedHitCollectionName",
                           "Collection name for corrected particle positions",
                           _correctedHitColName, string ("alignedHit"));
  /*
  registerInputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _inputTrackColName, string ("testfittracks"));

  registerInputCollection(LCIO::TRACKERHIT,"CorrectedHitCollectionName",
                           "Collection name for corrected particle positions",
                           _correctedHitColName, string ("corrfithits"));

  registerInputCollection(LCIO::TRACKERHIT,"OutputHitCollectionName",
                           "Collection name for fitted particle hits (positions)",
                           _inputHitColName, string ("testfithits"));
*/
  registerOptionalParameter("ReferenceCollection","This is the name of the reference it collection (init at 0,0,0)", _referenceHitCollectionName, static_cast< string > ( "refhit" ) );//Necessary for working out which layer the particle is detected in
  registerProcessorParameter("CutValue","Used to determine cuts in the system", _cutValue1, static_cast< double > (10.0));
  registerProcessorParameter ("DebugEventCount", "Print out every DebugEnevtCount event", _debugCount, static_cast < int > (100));//Not sure if I need this or not...


}

void EUTelX0Processor::processRunHeader(LCRunHeader *run)
{
}

void EUTelX0Processor::processEvent(LCEvent *evt)
{
  map< int, vector< TVector3> > hitInfoTemp;//Temporary map for storing vectors of hit information pre-cut
  LCCollection* col = evt->getCollection("alignedHit");//Create the collection of alignedHits for this event
  _referenceHitVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_referenceHitCollectionName));//Create the reference hit vector (used for figuring out what layer the hit is in)
  threePointResolution(col);
/*  try{
    LCCollection* col = evt->getCollection("track");//Create the collection of alignedHits for this event
    _referenceHitVec = dynamic_cast < LCCollectionVec * > (evt->getCollection(_referenceHitCollectionName));//Create the reference hit vector (used for figuring out what layer the hit is in)
    testtrack(col);
    createResiduals(col);
  }
  catch(...){
    cerr << "track collection did not exist in this event, moving to next event" << endl;
  }
*/

  
}

void EUTelX0Processor::createResiduals(LCCollection *trackCollection){
  int sizeElement = trackCollection->getNumberOfElements();
  for(int i = 0; i < sizeElement; ++i){
    Track* trackElement = dynamic_cast< Track* >(trackCollection->getElementAt(i));
    const vector< TrackerHit* > hits = trackElement->getTrackerHits();
    int sizeHits = hits.size();
    map< int, TVector3 > vecHits;
    for(int j = 0; j < sizeHits; ++j){
      if(hits[j]->getType() >= 32){  //Insist on a fitted hit
        TVector3 temp(hits[j]->getPosition()[0],hits[j]->getPosition()[1],hits[j]->getPosition()[2]);
        int layer = guessSensorID(hits[j]->getPosition());
        vecHits[layer] = temp;
      }
    }
    //Now we have all the hits of a track, we have to work out what the residual of the track is
    //Start with the forward residuals (TODO: (Phillip Hamnett) add the backward residuals too at some point)
    for(int j = 0; j < _noLayers; ++j){
      double dxf,dyf,dzf;
      if(j == 0){
        dxf = vecHits[j+1].x() - vecHits[j].x();
        dyf = vecHits[j+1].y() - vecHits[j].y();
        dzf = vecHits[j+1].z() - vecHits[j].z();
      }
      else if(j == 5){
        dxf = dyf = dzf = 0;
      }
      else{
        dxf = vecHits[j+1].x() - vecHits[j].x();
        dyf = vecHits[j+1].y() - vecHits[j].y();
        dzf = vecHits[j+1].z() - vecHits[j].z();
      }
      if(dxf != 0 || dyf != 0 || dzf != 0){
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualX"]))->fill(dxf/dzf);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualY"]))->fill(dyf/dzf);
//        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing[_histoResidualXY]))->fill(dxf/dyf);
      }
    }
  }
}

void EUTelX0Processor::testtrack(LCCollection *trackCollection){
  if(_debug){
    streamlog_out(DEBUG) << "trackCollection contains " << trackCollection->getNumberOfElements() << " elements" << endl;
    streamlog_out(DEBUG) << "The elements are as follows: " << endl;
    for(int i = 0; i < trackCollection->getNumberOfElements(); ++i){
      Track* trackElement = dynamic_cast< Track* >(trackCollection->getElementAt(i));
      streamlog_out(DEBUG) << "Element " << i << " contains the following values:" << endl;
      streamlog_out(DEBUG) << "D0 = " << trackElement->getD0() << endl
           << "Omega = " << trackElement->getOmega() << endl
           << "Phi = " << trackElement->getPhi() << endl
           << "Z0 = " << trackElement->getZ0() << endl
           << "Tan Lambda = " << trackElement->getTanLambda() << endl
           << "Reference Point = " << trackElement->getReferencePoint()[0] << ", " <<  trackElement->getReferencePoint()[1] << ", " << trackElement->getReferencePoint()[2] << endl
           << "Type = " << trackElement->getType() << endl
           << "Chi2 = " << trackElement->getChi2() << endl 
           << "Degrees of Freedom = " << trackElement->getNdf() << endl 
           << "The covariance matrix is as follows:" << endl;
      const vector<float> covMatrix = trackElement->getCovMatrix();
      for(size_t j = 0; j < covMatrix.size(); ++j){
        streamlog_out(DEBUG) << covMatrix[j] << endl;
      }
      streamlog_out(DEBUG) <<  "The track hits are finally listed here:" << endl;
      const vector<TrackerHit*> hits = trackElement->getTrackerHits();
      for(size_t j = 0; j < hits.size(); ++j){
        streamlog_out(DEBUG) << "Hit " << j << " = " << hits[j]->getPosition()[0] << ", " <<  hits[j]->getPosition()[1] << ", " << hits[j]->getPosition()[2]  << endl;
      }
      streamlog_out(DEBUG) << endl;
    }
  }
}

int EUTelX0Processor::guessSensorID(const double * hit )
{
  int sensorID = -1;
  double minDistance = numeric_limits< double >::max() ;
  streamlog_out(DEBUG) << "referencehit collection: " << _referenceHitCollectionName << " at "<< _referenceHitVec << endl;
  if( _referenceHitVec == 0)
  {
    streamlog_out(MESSAGE) << "_referenceHitVec is empty" << endl;
    return 0;
  }

  if( isFirstEvent() )
  {
    // message<DEBUG > ( log() << "number of elements : " << _referenceHitVec->getNumberOfElements() << endl );
  }

  for(int ii = 0 ; ii < _referenceHitVec->getNumberOfElements(); ii++)
  {
    EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
    if(refhit == 0 ) continue;
     TVector3 hit3d( hit[0], hit[1], hit[2] );
    TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
    TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
    double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
    streamlog_out(DEBUG) << " hit " << hit[0] << " "<< hit[1] << " " << hit[2] << endl;
    streamlog_out(DEBUG) << " " << refhit->getXOffset() << " " << refhit->getYOffset() << " " << refhit->getZOffset() << endl;
    streamlog_out(DEBUG) << " " << refhit->getAlpha() << " " << refhit->getBeta() << " " << refhit->getGamma() << endl;
    streamlog_out(DEBUG) << " distance " << distance << endl;
    if ( distance < minDistance )
    {
      minDistance = distance;
      sensorID = refhit->getSensorID();
    }
  }
  //Some useful debug printouts:
  bool debug = ( _debugCount>0 );

  if(debug)
  {
    for(int ii = 0 ; ii < _referenceHitVec->getNumberOfElements(); ii++)
    {
      EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
      if(refhit == 0 ) continue;
      // if( sensorID != refhit->getSensorID() ) continue;
      // streamlog_out(DEBUG) << " _referenceHitVec " << _referenceHitVec << " " << _referenceHitCollectionName.c_str() << " " << refhit << " at "
      // << refhit->getXOffset() << " " << refhit->getYOffset() << " " << refhit->getZOffset() << " "
      //<< refhit->getAlpha() << " " << refhit->getBeta() << " " << refhit->getGamma() << endl ;
      //message<DEBUG> ( log() << "iPlane " << refhit->getSensorID() << " hitPos: [" << hit[0] << " " << hit[1] << " " << hit[2] << "] distance: " << minDistance << endl );
    }
  }
  //cout << "Layer Number of this hit is " << sensorID << endl;
  return sensorID;
}
 
void EUTelX0Processor::end()
{
  //calculateX0();
  _hitInfo.clear();  //This stores the hit position in a TVector3. If there are multiple hits then they are all stored in the vector of TVector3's. The int key refers to the layer number
  _projectedHits.clear();  //int refers to 1 or 2, with 1 being the projection from the 01 plane and 2 being the projection from the 43 plane
  _histoThing.clear();  //This map contains all the histograms that are going to be plotted in this processor
  _residual.clear();  //The pair of doubles refers to the x and y coordinates in layer 2. The vector of TVector3's refers to the positions of the projected hits
  _residualAngle.clear();  //As above but the TVector3 doesn't contain xyz but instead theta, phi and alpha
  _residualProfile.clear(); //TODO(Phillip Hamnett): Can this be joined with _residual? //Used as above but for created a profile histogram
  _inputCollectionVec = NULL;  //Stores the information being brought in from the Marlin process which contains information about post-aligned hits
  _referenceHitVec = NULL;  
}

double EUTelX0Processor::calculateX0()
{
  double X0(0);
  //Two methods for figuring out the material budget, X0, are used. The loss of energy and the scattering angle due to multiple scattering.
  //Potentially might use both a Forward-backward Kalmin filter and a Global Linear Estimator to work out the two above values, might give a more accurate meaurement of X0.
  double sigma_msi(0);  //Highland formula
  double sigma_measi(0);  //Sum of the variance of histogram alpha012 and alpha432
  double E = 5;  //The total particle energy (GeV) TODO: (Phillip Hamnett) is there some way to make this information (and the mass) into the steering file?
  double m = 0.000511;  //Electron Mass (GeV)
  double p = sqrt(E*E - m*m);  //The particle momentum (GeV)
  double theta012(0);  //The polar angle between the momentum vector and the z-axis
  double beta012(0);  //The difference between the azimuthal angle of the momentum vector and the azimuthal angle of a point on the track
  double theta432(0);
  double beta432(0);
  X0 = 100;  //An estimate of X0
  double thickness = 0.00005;  //Thickness of mimosa layer is 50 micrometers
  double X = X0/thickness;
  double BetterX = DBL_MAX;
  double currentlikelihood = DBL_MIN;
  double lastX0 = DBL_MIN;
  while(X0 != BetterX){
    sigma_measi = 5;  //HACK - This is actually the sum of the variances of theta in the forward and backward direction. 
    size_t size = _histoData["Theta012"].size();
    double loglikelihoodX0(0);
    cout << "Size = " << size << endl;
    for(size_t i = 0; i < size; ++i){
      theta012 = _histoData["Theta012"][i]*3.1415/180.0;
      theta432 = _histoData["Theta432"][i]*3.1415/180.0;
      beta012 = _histoData["Phi012"][i]*3.1415/180.0;
      beta432 = _histoData["Phi432"][i]*3.1415/180.0;
      double deltatheta = sqrt((theta432 - theta012)*(theta432 - theta012));  //Take the positive difference between the two angles in theta and beta
      double deltabeta = sqrt((beta432 - beta012)*(beta432 - beta012));
      sigma_msi = (0.015*E/(p*p))*sqrt(X/(sin(deltatheta)*cos(deltabeta)))*(1 + 0.038*std::log(X/(sin(deltatheta)*cos(deltabeta))));
      if(sigma_msi != sigma_msi){
        cerr << "Code failed with following values: " << endl;
        cout << "E = " << E << endl << "p = " << p << endl << "X = " << X << endl << "theta012 = " << theta012 << endl << "beta012 = " << beta012 << endl << "theta432 = " << theta432 << endl << "beta432 = " << beta432 << endl << "deltatheta = " << deltatheta << endl << "deltabeta = " << deltabeta << endl << "sigma_msi = " << sigma_msi << endl << "sigma_measi = " << sigma_measi << endl << "loglikelihoodX0 = " << loglikelihoodX0 << endl << endl;
      }
      //Maximum log-likelihood:
      loglikelihoodX0 += -std::log(sigma_msi*sigma_msi + sigma_measi*sigma_measi) - (theta012*theta012)/(sigma_msi*sigma_msi + sigma_measi*sigma_measi); 

      cout << "E = " << E << endl << "p = " << p << endl << "X = " << X << endl << "theta012 = " << theta012 << endl << "beta012 = " << beta012 << endl << "sigma_msi = " << sigma_msi << endl << "sigma_measi = " << sigma_measi << endl << "loglikelihoodX0 = " << loglikelihoodX0 << endl << endl;
    }
    if(loglikelihoodX0 > currentlikelihood){
      BetterX = X0;
      currentlikelihood = loglikelihoodX0;
      if(lastX0 < X0){
        X++;
      }
      else{
        X--;
      }
    }
    else{
      if(lastX0 < X0){
        X--;
      }
      else{
        X++;
      }
    }
    X = X0/thickness;
    loglikelihoodX0 = 0;
  }

  cout << "The MAX Log Likelihood is " << currentlikelihood << endl;
  return X0;
}

void EUTelX0Processor::threePointResolution(LCCollection *alignedHitCollection){
  //cout << "Running threePointResolution function" << endl;
  map< int, vector< TVector3> > hitvectortemp;//Temporary map for storing vectors of hit information pre-cut
  int collectionsize = alignedHitCollection->getNumberOfElements();
  for(int i = 0; i < collectionsize; ++i){
    //cout << "Make a collection of the hits: " << i << endl;
    TrackerHit* tHit = dynamic_cast<TrackerHit*>(alignedHitCollection->getElementAt(i));//Get the hit from each element of the collection
    const double* pos = tHit->getPosition();//Get the position of the hit in x, y and z
    int layernumber = guessSensorID(pos);//Get the layer number of the hit
    Double_t X = pos[0], Y = pos[1], Z = pos[2];//Store the hit positions as Double_t's to keep ROOT happy
    TVector3 tempvec(X,Y,Z);//Add the positions to a TVector3
    hitvectortemp[layernumber].push_back(tempvec);//Then put it in the temporary map
    if(_debug == true){
      cout << "In for loop 'for(int i = 0; i < collectionsize; ++i)':" << endl
           << "layernumber = " << layernumber << ", X = " << X << ", Y = " << Y << ", Z = " << Z << endl
           << "hitvectortemp[" << layernumber << "].size() = " << hitvectortemp[layernumber].size() << endl;
    }
  }
  for(int i = 0; i < _noLayers - 1; ++i){
    size_t hitvecisize = hitvectortemp[i].size();
    size_t hitveci2size = hitvectortemp[i+2].size();
    size_t hitveci1size = hitvectortemp[i+1].size();
    //cout << "Size of each hitvec is: " << endl <<
    //  "hitvecisize = " << hitvecisize << endl << "hitveci1size = " << hitveci1size << endl << "hitveci2size = " << hitveci2size << endl;
    for(size_t j = 0; j < hitvecisize; ++j){
      TVector3 hiti = hitvectortemp[i][j];
      for(size_t k = 0; k < hitveci2size; ++k){
        TVector3 hiti2 = hitvectortemp[i+2][k];
        for(size_t l = 0; l < hitveci1size; ++l){
          //printf("%i,%i,%i,%i\n",i,j,k,l);
          TVector3 hiti1 = hitvectortemp[i+1][l];
          double averagex = (hiti.x() + hiti2.x())/2.0;
          double averagey = (hiti.y() + hiti2.y())/2.0;
          if(sqrt(pow((averagex-hiti1.x()),2) + pow((averagey-hiti1.y()),2)) < _cutValue1){
            (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualX"]))->fill(averagex - hiti1.x());
            (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualY"]))->fill(averagex - hiti1.x());
            if(i == 0){
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualXPlane1"]))->fill(averagex - hiti1.x());
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualYPlane1"]))->fill(averagex - hiti1.x());
            }
            else if(i == 1){
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualXPlane2"]))->fill(averagex - hiti1.x());
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualYPlane2"]))->fill(averagex - hiti1.x());
            }
            else if(i == 2){
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualXPlane3"]))->fill(averagex - hiti1.x());
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualYPlane3"]))->fill(averagex - hiti1.x());
            }
            else if(i == 3){
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualXPlane4"]))->fill(averagex - hiti1.x());
              (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualYPlane4"]))->fill(averagex - hiti1.x());
            }
          }
        }
      }
    }
  }
}

void EUTelX0Processor::basicFitter(LCCollection *alignedHitCollection){
  
  map< int, vector< TVector3> > hitvectortemp;//Temporary map for storing vectors of hit information pre-cut
  int collectionsize = alignedHitCollection->getNumberOfElements();
  for(int i = 0; i < collectionsize; ++i){
    TrackerHit* tHit = dynamic_cast<TrackerHit*>(alignedHitCollection->getElementAt(i));//Get the hit from each element of the collection
    const double* pos = tHit->getPosition();//Get the position of the hit in x, y and z
    int layernumber = guessSensorID(pos);//Get the layer number of the hit
    Double_t X = pos[0], Y = pos[1], Z = pos[2];//Store the hit positions as Double_t's to keep ROOT happy
    TVector3 tempvec(X,Y,Z);//Add the positions to a TVector3
    hitvectortemp[layernumber].push_back(tempvec);//Then put it in the temporary map
    if(_debug == true){
      cout << "In for loop 'for(int i = 0; i < collectionsize; ++i)':" << endl
           << "layernumber = " << layernumber << ", X = " << X << ", Y = " << Y << ", Z = " << Z << endl
           << "hitvectortemp[" << layernumber << "].size() = " << hitvectortemp[layernumber].size() << endl;
    }
  }//End of: for(int i = 0; i < colSize; ++i)
  if(_debug == true){
    cout << "So now, hitvectortemp contains: " << endl;
    for(map< int, vector< TVector3 > >::iterator it = hitvectortemp.begin(); it != hitvectortemp.end(); ++it){
      int layer = it->first;
      cout << "Layer " << layer << ":" << endl;
      for(vector< TVector3 >::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2){
        cout  << it2->x() << ", " << it2->y() << ", " << it2->z() << endl;
      }
    }
  }
  double layerdifference = 0.1;  //This is the radius of a circle in which the two particles must lie
  vector< pair< TVector3, TVector3 > > relatedhits01;
  if(_debug == true){
    cout << "The following hits are within each others acceptable position difference (radius) of " << layerdifference << ":" << endl;
  }
  for(vector< TVector3 >::iterator hit0 = hitvectortemp[0].begin(); hit0 != hitvectortemp[0].end(); ++hit0){
    for(vector< TVector3 >::iterator hit1 = hitvectortemp[1].begin(); hit1 != hitvectortemp[1].end(); ++hit1){
      if(sqrt(pow(hit0->x() - hit1->x(),2) + pow(hit0->y() - hit1->y(),2)) < layerdifference){
        if(_debug == true){
          cout << "Layer 0: " << hit0->x() << ", " << hit0->y() << ", " << hit0->z() << endl;
          cout << "Layer 1: " << hit1->x() << ", " << hit1->y() << ", " << hit1->z() << endl;
        }
        pair< TVector3, TVector3 > relatedhits(*hit0,*hit1);
        relatedhits01.push_back(relatedhits);
      }
    }
  }
  if(_debug == true){
    cout << "The number of related hits in Layers 0 and 1 is:  " << relatedhits01.size() << endl;
  }
  //Make a projected point on layer 2 based on all related hits
  vector< TVector3 > projectedhits01;
  for(vector< pair< TVector3, TVector3 > >::iterator it = relatedhits01.begin(); it != relatedhits01.end(); ++it){
    TVector3 projectedhit(2*it->second.x() - it->first.x(), 2*it->second.y() - it->first.y(), 2*it->second.z());
    projectedhits01.push_back(projectedhit);
  }
  if(_debug == true){
    cout << "These are the predicted positions in layer 2 based on layers 0 and 1:" << endl;
    for(vector< TVector3 >::iterator it = projectedhits01.begin(); it != projectedhits01.end(); ++it){
      cout << it->x() << ", " << it->y() << ", " << it->z() << endl;
    }
  }
  //Clean up
  relatedhits01.clear();
  //Now to repeat on the other side, and then compare results
  vector< pair< TVector3, TVector3 > > relatedhits43;
  if(_debug == true){
    cout << "The following hits are within each others acceptable position difference (radius) of " << layerdifference << ":" << endl;
  }
  for(vector< TVector3 >::iterator hit4 = hitvectortemp[4].begin(); hit4 != hitvectortemp[4].end(); ++hit4){
    for(vector< TVector3 >::iterator hit3 = hitvectortemp[3].begin(); hit3 != hitvectortemp[3].end(); ++hit3){
      if(sqrt(pow(hit4->x() - hit3->x(),2) + pow(hit4->y() - hit3->y(),2)) < layerdifference){
        if(_debug == true){
          cout << "Layer 4: " << hit4->x() << ", " << hit4->y() << ", " << hit4->z() << endl;
          cout << "Layer 3: " << hit3->x() << ", " << hit3->y() << ", " << hit3->z() << endl;
        }
        pair< TVector3, TVector3 > relatedhits(*hit4,*hit3);
        relatedhits43.push_back(relatedhits);
      }
    }
  }
  if(_debug == true){
    cout << "The number of related hits in Layers 4 and 3 is:  " << relatedhits43.size() << endl;
  }
  //Make a projected point on layer 2 based on all related hits
  vector< TVector3 > projectedhits43;
  for(vector< pair< TVector3, TVector3 > >::iterator it = relatedhits43.begin(); it != relatedhits43.end(); ++it){
    TVector3 projectedhit(2*it->second.x() - it->first.x(), 2*it->second.y() - it->first.y(), 2*it->second.z());
    projectedhits43.push_back(projectedhit);
  }
  if(_debug == true){
    cout << "These are the predicted positions in layer 2 based on layers 4 and 3:" << endl;
    for(vector< TVector3 >::iterator it = projectedhits43.begin(); it != projectedhits43.end(); ++it){
      cout << it->x() << ", " << it->y() << ", " << it->z() << endl;
    }
  }
  //Clean up
  relatedhits43.clear();
  cout << "projectedhits01.size() = " << projectedhits01.size() << endl;
  cout << "projectedhits43.size() = " << projectedhits43.size() << endl;
  //Now we have projected hits in the 2nd layer from both the forward and backward direction. Now we make a residual between the two assuming that anything within a radius of 'layerdifference' of each other should be taken account of only
  for(vector< TVector3 >::iterator hits012 = projectedhits01.begin(); hits012 != projectedhits01.end(); ++hits012){
    for(vector< TVector3 >::iterator hits432 = projectedhits43.begin(); hits432 != projectedhits43.end(); ++hits432){
      if(sqrt(pow(hits012->x() - hits432->x(),2) + pow(hits012->y() - hits432->y(),2)) < layerdifference){
        cout << "Hit recorded in histogram" << endl;
        double newx = hits012->x() - hits432->x();
        double newy = hits012->y() - hits432->y();
        double newz = 20;  //HACK - find a better way of deducing the z position of layer 2
        Double_t theta012 = (atan(hits012->x()/newz)/3.14)*180;
        Double_t phi012 = (atan(hits012->y()/newz)/3.14)*180;
        Double_t alpha012 = tan(theta012)*tan(phi012)/sqrt(tan(theta012)*tan(theta012) + tan(phi012)*tan(phi012));
        Double_t theta432 = (atan(hits432->x()/newz)/3.14)*180;
        Double_t phi432 = (atan(hits432->y()/newz)/3.14)*180;
        Double_t alpha432 = tan(theta432)*tan(phi432)/sqrt(tan(theta432)*tan(theta432) + tan(phi432)*tan(phi432));
        if(_debug == true){
          cout << newx << ", " << newy << ", " << newz << endl;
          cout << theta012 << ", " << phi012 << ", " << alpha012 << endl;
        }
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualX"]))->fill(newx);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["ResidualY"]))->fill(newy);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["Theta012"]))->fill(theta012);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["Phi012"]))->fill(phi012);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["Alpha012"]))->fill(alpha012);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["Theta432"]))->fill(theta432);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["Phi432"]))->fill(phi432);
        (dynamic_cast< AIDA::IHistogram1D* > (_histoThing["Alpha432"]))->fill(alpha432);
        _histoData["ResidualX"].push_back(newx);
        _histoData["ResidualY"].push_back(newy);
        _histoData["Theta012"].push_back(theta012);
        _histoData["Phi012"].push_back(phi012);
        _histoData["Alpha012"].push_back(alpha012);
        _histoData["Theta432"].push_back(theta432);
        _histoData["Phi432"].push_back(phi432);
        _histoData["Alpha432"].push_back(alpha432);
      }
    }
  }
  projectedhits01.clear();
  projectedhits43.clear();

}

