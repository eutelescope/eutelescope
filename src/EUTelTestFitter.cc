// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version: $Id: EUTelTestFitter.cc,v 1.1 2007-06-04 06:45:03 bulgheroni Exp $
// Date 2007.06.04

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#include "EUTelTestFitter.h"
#include <iostream>
#include <fstream>
#include <cmath>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
// #include <EVENT/TrackerHit.h>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/LCFlagImpl.h"


using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;


EUTelTestFitter::EUTelTestFitter() : Processor("EUTelTestFitter") {
  
  // modify processor description
  _description = "Analytical track fitting processor for EUDET telescope" ;
  

  // register steering parameters: 
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACKERHIT,
			   "InputCollectionName" , 
			   "Name of the input TrackerHit collection"  ,
			   _inputColName ,
			   std::string("meshit") ) ;

  // compulsory parameters:

  registerProcessorParameter ("DebugEventCount",
			      "Print out every DebugEnevtCount event",
			      _debugCount,  static_cast < int > (100));

  registerProcessorParameter ("GeometryFileName",
                              "Name of the geometry description file",
                              _geometryFileName,
                              string ("EUTelTestFitter.geom"));

  registerProcessorParameter ("AllowMissingHits",
			      "Allowed number of missing hits in the track",
			      _allowMissingHits,  static_cast < int > (0));

  registerProcessorParameter ("AllowSkipHits",
			      "Allowed number of hits removed from the track",
			      _allowSkipHits,  static_cast < int > (0));

  registerProcessorParameter ("MaxPlaneHits",
			      "Maximum number of considered hits per plane",
			      _maxPlaneHits,  static_cast < int > (5));

  registerProcessorParameter ("MissingHitPenalty",
			      "Chi2 penalty for missing hit in the track",
			      _missingHitPenalty,  static_cast < double > (0.));

  registerProcessorParameter ("SkipHitPenalty",
			      "Chi2 penalty for removing hit from the track",
			      _skipHitPenalty,  static_cast < double > (100.));

  registerProcessorParameter ("Chi2Max",
			      "Maximum Chi2 for accepted track fit",
			      _chi2Max,  static_cast < double > (1000.));

  registerProcessorParameter ("UseNominalResolution",
			      "Flag for using nominal resolution instead of position errors",
			      _useNominalResolution,  static_cast < bool > (false));

  registerProcessorParameter ("UseDUT",
                              "Flag for including DUT measurement in the fit",
			      _useDUT,  static_cast < bool > (false));

  registerProcessorParameter ("Ebeam",
			      "Beam energy [GeV]",
			      _eBeam,  static_cast < double > (6.0));

  // optional parameters

  registerOptionalParameter ("OutputTrackCollectionName",
                             "Collection name for fitted tracks",
                             _outputTrackColName, string ("testfittracks"));

  registerOptionalParameter ("OutputHitCollectionName",
                             "Collection name for fitted particle hits (positions)",
                             _outputHitColName, string ("testfithits"));

  registerOptionalParameter ("UseBeamConstraint",
			     "Flag for using beam direction constraint in the fit",
                             _useBeamConstraint,  static_cast < bool > (false));

  registerOptionalParameter ("BeamSpread",
                             "Assumed angular spread of the beam [rad]",
                             _beamSpread,  static_cast < double > (0.1));

  registerOptionalParameter ("SearchMultipleTracks",
			     "Flag for searching multiple tracks in events with multiple hits",
                             _searchMultipleTracks,  static_cast < bool > (false));

}


void EUTelTestFitter::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // Read geometry information from file

  ifstream geometryFile;
  geometryFile.exceptions(ifstream::failbit | ifstream::badbit);

  try {
    geometryFile.open(_geometryFileName.c_str(),ios::in);
  }
  catch (exception& e) {
    cerr << "IO exception " << e.what() << " with " 
	 << _geometryFileName << ".\nExiting." << endl;
    return;
  }

  cout << "Reading telescope geometry description from " 
       << _geometryFileName << endl ;

  geometryFile >> _nTelPlanes >> _iDUT ;
  _iDUT-- ;

  _planePosition   = new double[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];
  _nActivePlanes = 0 ;

  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      int iActive;
      double resolution;

      // All dimensions should be given in mm !!!

      geometryFile >> _planePosition[ipl]
                   >> _planeThickness[ipl]
                   >> _planeX0[ipl]
                   >> iActive 
                   >> resolution ;

      if(iActive && (ipl != _iDUT || _useDUT ))
	{
	  _isActive[ipl] = true ;
	  _planeResolution[ipl]=resolution;
	  _nActivePlanes++ ;
	}
      else
	{
          _isActive[ipl] = false ;
          _planeResolution[ipl]=0.;
	}
    }

  // Print out geometry information

  cout << "Telescope configuration with " << _nTelPlanes << " planes" <<  endl;

  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      if(ipl == _iDUT)
	cout << "D.U.T.  plane at" ;
      else
	if(_isActive[ipl])
	  cout << "Active  plane at" ;
	else
	  cout << "Passive plane at" ; 

      cout << "  Z [mm] = " << _planePosition[ipl] 
           << " dZ [um] = " << _planeThickness[ipl]*1000. ;

      if(_isActive[ipl])
	cout << "  Res [um] = " << _planeResolution[ipl]*1000. ;

      cout << endl;
    }

  cout << "Total of " << _nActivePlanes << " active sensor planes " << endl;

  // Allocate arrays for track fitting

  _planeX  = new double[_nTelPlanes];
  _planeEx = new double[_nTelPlanes];
  _planeY  = new double[_nTelPlanes];
  _planeEy = new double[_nTelPlanes];

  _planeDist = new double[_nTelPlanes];
  _planeScat = new double[_nTelPlanes];

  _fitX  = new double[_nTelPlanes];
  _fitEx = new double[_nTelPlanes];
  _fitY  = new double[_nTelPlanes];
  _fitEy = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];
  _nominalFitArray = new double[arrayDim];

  _nominalError = new double[_nTelPlanes];

  // Fill nominal fit matrices and 
  // calculate expected precision of track fitting

  // Planes have to be ordered in position along the beam line !
  // This is not checked !!!

  for(int ipl=0; ipl<_nTelPlanes ; ipl++)
    {
      if(ipl>0)
	_planeDist[ipl-1]=1./(_planePosition[ipl] - _planePosition[ipl-1]) ;

      _planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl])
	* (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl])) ;

      if(ipl==0 && _useBeamConstraint)
	_planeScat[ipl]= 1./(_planeScat[ipl]*_planeScat[ipl]+ _beamSpread*_beamSpread) ; 
      else
	_planeScat[ipl]= 1./(_planeScat[ipl] * _planeScat[ipl]) ; 

      _fitX[ipl] =_fitY[ipl] = 0. ;
      _nominalError[ipl]= _planeResolution[ipl];
    }

  // Fit with nominal parameters

  int status = DoAnalFit(_fitX,_nominalError);

  if(status)
    cerr << "\n Fit with nominal geometry failed !?!" << endl ;

  // Store fit matrix

  for(int imx=0; imx<arrayDim; imx++)
    _nominalFitArray[imx] = _fitArray[imx];

  cout << "Expected position resolutions [um]: " << endl ;
  for(int ipl=0; ipl<_nTelPlanes ; ipl++)
    cout << _nominalError[ipl]*1000. << "  " ;
  cout << endl;
}

void EUTelTestFitter::processRunHeader( LCRunHeader* runHeader) { 

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  int runNr = runHeader->getRunNumber();

  cout << "\n Processing run header " << _nRun 
       << ", run nr " << runNr << endl ;

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  cout << detectorName << " : " << detectorDescription << endl;

  int nDet = subDets->size();

  cout << nDet << " subdetectors defined :" << endl;
  for(int idet=0;idet<nDet;idet++)
    cout << idet+1 << " : " << subDets->at(idet) << endl;

} 

void EUTelTestFitter::processEvent( LCEvent * event ) { 

  bool debug = ( _debugCount>0 && _nEvt%_debugCount == 0) ;

  _nEvt ++ ;
  int evtNr = event->getEventNumber();

  if(debug)
    cout << "\n EUTelTestFitter: Processing " << _nEvt 
	 << " event record, event nr " << evtNr << endl ;

  LCCollection* col = event->getCollection( _inputColName ) ;

  if( col == 0 ) 
    {
      cerr << "Missing collection " << _inputColName 
	   << " in event " <<   _nEvt << " / " << evtNr << endl ;

      return ;
    }

  // Copy hits to local table
  // Assign hits to sensor planes
  // =============================


  int nHit = col->getNumberOfElements()  ;

  if(debug)
    cout << "Total of " << nHit << " tracker hits in input collection " << endl ;

  if(nHit + _allowMissingHits < _nActivePlanes)
    {
      if(debug) 
	cout << "Not enough hits to perform the fit, exiting... " << endl ;
      return;
    }


  double * hitX  = new double[nHit];
  double * hitEx = new double[nHit];
  double * hitY  = new double[nHit];
  double * hitEy = new double[nHit];
  double * hitZ  = new double[nHit];
  int * hitPlane = new int[nHit];

  int * nPlaneHits = new int[_nTelPlanes];
  int * nPlaneChoice = new int[_nTelPlanes];
  IntVec * planeHitID = new IntVec[_nTelPlanes];

  double * bestX  = new double[_nTelPlanes];
  double * bestEx = new double[_nTelPlanes];
  double * bestY  = new double[_nTelPlanes];
  double * bestEy = new double[_nTelPlanes];

  // Loop over hits

  int nGoodHit = 0;

  for(int ihit=0; ihit< nHit ; ihit++)
    {
      TrackerHit * meshit = dynamic_cast<TrackerHit*>( col->getElementAt(ihit) ) ;

      // Hit position

      const double * pos = meshit->getPosition();

      hitX[ihit] = pos[0];
      hitY[ihit] = pos[1];
      hitZ[ihit] = pos[2];

      // Plane ID should be stored as hit type 
      //      hitPlane[ihit] = meshit->getType() - 1 ;


      // If not, we have to find Plane ID of the hit
      // by looking at the Z position

      double distMin = 1.;
      hitPlane[ihit] = -1 ;

      for(int ipl=0;ipl<_nTelPlanes;ipl++)
	{
	  double dist =  hitZ[ihit] - _planePosition[ipl] ;

	  if(dist*dist < distMin*distMin)
	    {
	      hitPlane[ihit]=ipl;
	      distMin=dist;
	    }
	}

      // Ignore hits not matched to any plane

      if(hitPlane[ihit]<0)
	{
	  cerr << "Reconstructed hit outside sensor planes z [mm] = " 
	       << hitZ[ihit] << endl ;
	  continue;
	}

      // Ignore hit, if plane not declared as active (i.e. not used in the fit)

      if(! _isActive[hitPlane[ihit]])
	continue ;

      // Ignore hit also, if maximum number of hits already matched to this plane

      if(_maxPlaneHits>0 && 
         _maxPlaneHits-planeHitID[hitPlane[ihit]].size()<=0) 
	continue ;

      // Add hit to hit list for given plane - to be used in track selection

      planeHitID[hitPlane[ihit]].push_back(ihit);
      nGoodHit++;

      // Position uncertainty. Use nominal resolution if not properly defined

      const EVENT::FloatVec cov = meshit->getCovMatrix();

      if(cov.at(0)>0.)
	hitEx[ihit]=sqrt(cov.at(0));
      else
	hitEx[ihit]=_planeResolution[hitPlane[ihit]];

      if(cov.at(2)>0.)
	hitEy[ihit]=sqrt(cov.at(2));
      else
	hitEy[ihit]=_planeResolution[hitPlane[ihit]];


      if(debug)
        cout << "Hit " << ihit
             << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]  
             << "   Y = " << hitY[ihit] << " +/- " << hitEy[ihit]  
             << "   Z = " << hitZ[ihit] << " (plane" << hitPlane[ihit] << ")" << endl ;
    }


  // Main analysis loop: finding multiple tracks (if allowed)

  // Define output track and hit collections
  LCCollectionVec     * fittrackvec = new LCCollectionVec(LCIO::TRACK);
  LCCollectionVec     * fitpointvec = new LCCollectionVec(LCIO::TRACKERHIT);

  // Set flag for storing track hits in track collection

  LCFlagImpl flag(fittrackvec->getFlag()); 
  flag.setBit( LCIO::TRBIT_HITS );
  fittrackvec->setFlag(flag.getFlag());


  int nFittedTracks = 0 ;
  bool firstTrack = true ;
  int ibest;

  do{

    // Count planes active in this event and number of fit possibilities

    int nFiredPlanes = 0;
    int nChoice = 1;
    ibest=-1;

    for(int ipl=0;ipl<_nTelPlanes;ipl++)
      {
	nPlaneHits[ipl] = planeHitID[ipl].size() ;

	if(nPlaneHits[ipl]>0)
	  {
	    nFiredPlanes++;

	    if(_allowSkipHits>0)
	      nPlaneChoice[ipl]=nPlaneHits[ipl]+1;
	    else
	      nPlaneChoice[ipl]=nPlaneHits[ipl];
	  }
	else
	  nPlaneChoice[ipl]=1;

	nChoice*=nPlaneChoice[ipl];

	if(debug && _isActive[ipl] && firstTrack)
	  {
	    cout << "Plane " << ipl << "  " << nPlaneHits[ipl] << " hit(s), hit IDs : " ;
	    for( int ihit=0; ihit < (int) planeHitID[ipl].size() ; ihit ++)
	      cout << planeHitID[ipl].at(ihit) << " " ;
	    cout << endl;
	  }
      }
 
    // Check if fit can be done

    if(nFiredPlanes + _allowMissingHits < _nActivePlanes)
      {
        if(debug && firstTrack)
	  cout << "Not enough planes hit to perform the fit, exiting... " << endl ;

        break;
      }


    if(debug && firstTrack)
      cout << nFiredPlanes << " active sensor planes hit, checking "
	   << nChoice << " fit possibilities " << endl ;

    // Check all track possibilities

    double chi2best = _chi2Max;
    double bestPenalty = 0;
    int nBestFired = 0;

    // Loop over fit possibilities

    for(int ichoice=0; ichoice<nChoice; ichoice++)
      {
	int modchoice=ichoice;

	int nChoiceFired=0;
	double choiceChi2=-1.;

	// Fill position and error arrays for this hit configuration

	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  {
	    _planeX[ipl]=_planeY[ipl]=_planeEx[ipl]=_planeEy[ipl]=0.;

	    if(_isActive[ipl])
	      {
		int ihit=modchoice%nPlaneChoice[ipl];
		modchoice/=nPlaneChoice[ipl];

		if(ihit<nPlaneHits[ipl])
		  {
		    int jhit = planeHitID[ipl].at(ihit);
		    _planeX[ipl]=hitX[jhit];
		    _planeY[ipl]=hitY[jhit];
		    _planeEx[ipl]=(_useNominalResolution)?_planeResolution[ipl]:hitEx[jhit];
		    _planeEy[ipl]=(_useNominalResolution)?_planeResolution[ipl]:hitEy[jhit];
		    nChoiceFired++;
		  }
	      }
	  }

	// Check number of selected hits

	if(nChoiceFired + _allowMissingHits < _nActivePlanes ||
	   nChoiceFired + _allowSkipHits    < nFiredPlanes )
	  continue;

	// Select fit method
	// "Nominal" fit only if all active planes used

	if(_useNominalResolution && (nChoiceFired == _nActivePlanes))
          choiceChi2 = NominalFit();
	else
	  if(_useNominalResolution)
	    choiceChi2 = SingleFit();
	  else
	    choiceChi2 = MatrixFit();

	// Fit failed ?

	if(choiceChi2 < 0.)continue ;

	// Penalty for missing or skiped hits

	double penalty = (_nActivePlanes-nFiredPlanes)*_missingHitPenalty
	  +   (nFiredPlanes-nChoiceFired)*_skipHitPenalty ;

	// Best fit ?

	if(choiceChi2+penalty<chi2best)
          {
	    chi2best=choiceChi2+penalty;
	    bestPenalty=penalty;
	    ibest=ichoice;
	    nBestFired=nChoiceFired;
	    for(int ipl=0;ipl<_nTelPlanes;ipl++)
	      {
		bestX[ipl]=_fitX[ipl];
		bestY[ipl]=_fitY[ipl];
		bestEx[ipl]=_fitEx[ipl];
		bestEy[ipl]=_fitEy[ipl];
	      }
	  }
      }
    // End of loop over track possibilities

    if(ibest<0 && debug && firstTrack)
      cout << "No track fulfilling search criteria found ! " << endl ;

    if(ibest>=0)
      {
	if(debug)
	  {
	    cout << "Best fit from choice " << ibest << " : " << nBestFired << " hits " 
		 << " chi2 = " << chi2best << " including penalties of " << bestPenalty << endl ;

	    for(int ipl=0;ipl<_nTelPlanes;ipl++)
	      cout << "  X = " << bestX[ipl] << " +/- " << bestEx[ipl] 
		   << "  Y = " << bestY[ipl] << " +/- " << bestEy[ipl] 
		   << "  at Z = " << _planePosition[ipl] << endl ;
	  }

	// Write fit result out

	TrackImpl * fittrack = new TrackImpl();

	// Following parameters are not used for Telescope
	// and are set to zero (just in case)
	fittrack->setOmega(0.);     // curvature of the track
	fittrack->setD0(0.);        // impact paramter of the track in (r-phi)
	fittrack->setZ0(0.);        // impact paramter of the track in (r-z)
	fittrack->setPhi(0.);       // phi of the track at reference point
	fittrack->setTanLambda(0.); // dip angle of the track at reference point

	// Used class members

	fittrack->setChi2(chi2best);  // Chi2 of the fit (including penalties)
	fittrack->setNdf(nBestFired); // Number of planes fired (!)

	// Fitted position at DUT is stored as Reference Point
	// (see below in loop over telescope planes)
	// If no DUT present: use position in the first plane !

	fittrack->setIsReferencePointPCA(false);  
	float refpoint[3];

	// Track points fitted in each plane are stored as track hits !!!

	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  {
	    TrackerHitImpl * fitpoint = new TrackerHitImpl;

	    // Plane number stored as hit type

	    fitpoint->setType(ipl+1);

	    // fitted position in a plane

	    double pos[3];

	    pos[0]=bestX[ipl];
	    pos[1]=bestY[ipl];
	    pos[2]=_planePosition[ipl];

	    fitpoint->setPosition(pos);

	    // Covariance matrix of the position 
	    // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).

	    float cov[TRKHITNCOVMATRIX];

	    cov[0]=bestEx[ipl]*bestEx[ipl];
	    cov[1]=0.;
	    cov[2]=bestEy[ipl]*bestEy[ipl];
	    cov[3]=cov[4]=0.;
	    cov[5]=_planeThickness[ipl]*_planeThickness[ipl]/12.;

	    fitpoint->setCovMatrix(cov);

	    // store fit point 

	    fitpointvec->push_back(fitpoint);

	    //   add point to track

	    fittrack->addHit(fitpoint);

	    // Use position at DUT as a track reference point.
	    // If no DUT present: use position in the first plane !

	    if(ipl==_iDUT || (_iDUT<0 && ipl==0))
	      for(int iref=0;iref<3;iref++)
		refpoint[iref]=pos[iref];
	  }

	// Store track reference point.
	fittrack->setReferencePoint(refpoint);

	fittrackvec->addElement(fittrack);
	nFittedTracks++;
      }

    // If multiple tracks allowed: remove hits from fitted track from the list

    if(ibest>=0 && _searchMultipleTracks)
      {
	int modchoice=ibest;

	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	  if(_isActive[ipl])
	    {
	      int ihit=modchoice%nPlaneChoice[ipl];
	      modchoice/=nPlaneChoice[ipl];

	      if(ihit<nPlaneHits[ipl])
		{
                  planeHitID[ipl].erase(planeHitID[ipl].begin()+ihit);
                  nGoodHit--;
		}
	    }
      }

    // Suppress output for secondary tracks

    // firstTrack=false;
  }
  while(_searchMultipleTracks && ibest>=0 && nGoodHit + _allowMissingHits >= _nActivePlanes);

  // End of track loop

  if(nFittedTracks > 0 )
    {
      event->addCollection(fittrackvec,_outputTrackColName);
      event->addCollection(fitpointvec,_outputHitColName);
    }
  else
    {
      delete fittrackvec;
      delete fitpointvec;
    }

  // Clear all working arrays

  delete hitX;
  delete hitEx;
  delete hitY;
  delete hitEy;
  delete hitZ;
  delete hitPlane;

  delete nPlaneHits;
  delete nPlaneChoice;
  delete planeHitID;

  delete bestX;
  delete bestEx;
  delete bestY;
  delete bestEy;

  return;
}



void EUTelTestFitter::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelTestFitter::end(){ 
  
  //   std::cout << "EUTelTestFitter::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;


  // Clean memory 

  delete _planePosition ;
  delete _planeThickness  ;
  delete _planeX0  ;
  delete _planeResolution ;
  delete _planeDist ;
  delete _planeScat ;
  delete _isActive ;

  delete _planeX ;
  delete _planeEx ;
  delete _planeY ;
  delete _planeEy ;

  delete _fitX  ;
  delete _fitEx ;
  delete _fitY ;
  delete _fitEy ;
  delete _fitArray ;

  delete _nominalFitArray ;
  delete _nominalError ;
}

//
// ===============================================================================
//
//  Private function members
//


double EUTelTestFitter::MatrixFit()
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitX[ipl]=_planeX[ipl];
      _fitEx[ipl]=_planeEx[ipl];
      _fitY[ipl]=_planeY[ipl];
      _fitEy[ipl]=_planeEy[ipl];
    }

  int status = DoAnalFit(_fitX,_fitEx);

  if(status)return -1. ;

  status = DoAnalFit(_fitY,_fitEy);

  if(status)return -1. ;

  double chi2=GetFitChi2();

  return chi2 ;
}

double EUTelTestFitter::SingleFit()
{

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitX[ipl]=_planeX[ipl];
      _fitEx[ipl]=_planeEx[ipl];
    }

  int status = DoAnalFit(_fitX,_fitEx);

  if(status)return -1. ;

  // Use same matrix to solve equation in Y

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitEy[ipl]=_fitEx[ipl];

      _fitY[ipl]=0. ;
      for(int jpl=0; jpl<_nTelPlanes;jpl++)
	if(_planeEy[jpl]>0.)
	  _fitY[ipl]+=_fitArray[ipl+jpl*_nTelPlanes]*_planeY[jpl]/_planeEy[jpl]/_planeEy[jpl];
    }

  double chi2=GetFitChi2();

  return chi2 ;
}

double EUTelTestFitter::NominalFit()
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      _fitEy[ipl]=_fitEx[ipl]=_nominalError[ipl];

      _fitX[ipl]=0. ;
      _fitY[ipl]=0. ;
      for(int jpl=0; jpl<_nTelPlanes;jpl++)
	{
	  if(_planeEx[jpl]>0.)
	    _fitX[ipl]+=_nominalFitArray[ipl+jpl*_nTelPlanes]*_planeX[jpl]/_planeEx[jpl]/_planeEx[jpl];
	  if(_planeEy[jpl]>0.)
	    _fitY[ipl]+=_nominalFitArray[ipl+jpl*_nTelPlanes]*_planeY[jpl]/_planeEy[jpl]/_planeEy[jpl];
	}
    }

  double chi2=GetFitChi2();

  return chi2 ;
}


int EUTelTestFitter::DoAnalFit(double * pos, double *err)
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      if(_isActive[ipl] && err[ipl]>0)
        err[ipl]=1./err[ipl]/err[ipl] ;
      else
        err[ipl] = 0. ;

      pos[ipl]*=err[ipl];
    }


  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    for(int jpl=0; jpl<_nTelPlanes;jpl++)
      {
	int imx=ipl+jpl*_nTelPlanes;

	_fitArray[imx] = 0.;

	if(jpl==ipl-2)
	  _fitArray[imx] += _planeDist[ipl-2]*_planeDist[ipl-1]*_planeScat[ipl-1] ;

	if(jpl==ipl+2)
	  _fitArray[imx] += _planeDist[ipl]*_planeDist[ipl+1]*_planeScat[ipl+1] ;

	if(jpl==ipl-1)
	  {
	    if(ipl>0 &&  ipl < _nTelPlanes-1)
              _fitArray[imx] -= _planeDist[ipl-1]*(_planeDist[ipl]+_planeDist[ipl-1])*_planeScat[ipl] ;
	    if(ipl>1)
              _fitArray[imx] -= _planeDist[ipl-1]*(_planeDist[ipl-1]+_planeDist[ipl-2])*_planeScat[ipl-1] ;
	  }

	if(jpl==ipl+1)
	  {
	    if(ipl>0 && ipl < _nTelPlanes-1)
	      _fitArray[imx] -= _planeDist[ipl]*(_planeDist[ipl]+_planeDist[ipl-1])*_planeScat[ipl] ;
	    if(ipl < _nTelPlanes-2)
	      _fitArray[imx] -= _planeDist[ipl]*(_planeDist[ipl+1]+_planeDist[ipl])*_planeScat[ipl+1] ;
	  }

	if(jpl==ipl)
	  {
            _fitArray[imx] += err[ipl] ;

            if(ipl>0 && ipl<_nTelPlanes-1)
	      _fitArray[imx] += _planeScat[ipl]*(_planeDist[ipl]+_planeDist[ipl-1])*(_planeDist[ipl]+_planeDist[ipl-1]) ;

            if(ipl > 1 )
	      _fitArray[imx] += _planeScat[ipl-1]*_planeDist[ipl-1]*_planeDist[ipl-1] ;

            if(ipl < _nTelPlanes-2)
	      _fitArray[imx] += _planeScat[ipl+1]*_planeDist[ipl]*_planeDist[ipl] ;
	  }

	// For beam constraint

	if(ipl==jpl && ipl<2 && _useBeamConstraint)
	  _fitArray[imx] += _planeScat[1]*_planeDist[1]*_planeDist[1] ;

	if(ipl+jpl==1 && _useBeamConstraint) 
	  _fitArray[imx] -= _planeScat[1]*_planeDist[1]*_planeDist[1] ;


      }

  int status=GaussjSolve(_fitArray,pos,_nTelPlanes) ;

  if(status)
    {
      cerr << "Singular matrix in track fitting algorithm ! " << endl;
      for(int ipl=0;ipl<_nTelPlanes;ipl++)
	err[ipl]=0. ;
    }
  else
    for(int ipl=0;ipl<_nTelPlanes;ipl++)
      err[ipl]=sqrt(_fitArray[ipl+ipl*_nTelPlanes]);

  return status ;
}



double EUTelTestFitter::GetFitChi2()
{
  double chi2=0. ;

  // Measurements

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    if(_isActive[ipl])
      {
	if(_planeEx[ipl]>0.)
	  chi2+=(_fitX[ipl]-_planeX[ipl])*(_fitX[ipl]-_planeX[ipl])/_planeEx[ipl]/_planeEx[ipl] ;

	if(_planeEy[ipl]>0.)
	  chi2+=(_fitY[ipl]-_planeY[ipl])*(_fitY[ipl]-_planeY[ipl])/_planeEy[ipl]/_planeEy[ipl] ;

      }

  // Scattering angles
  // Use approximate formulas, corresponding to the approximation 
  // used in fitting algorithm

  for(int ipl=1; ipl<_nTelPlanes-1;ipl++)
    {
      double th1,th2,dth;
    
      th2=(_fitX[ipl+1]-_fitX[ipl])*_planeDist[ipl] ;
      th1=(_fitX[ipl]-_fitX[ipl-1])*_planeDist[ipl-1] ;
      //    dth=atan(th2)-atan(th1) ;
      dth= th2 - th1 ;
      chi2 += _planeScat[ipl] * dth * dth;

      th2=(_fitY[ipl+1]-_fitY[ipl])*_planeDist[ipl] ;
      th1=(_fitY[ipl]-_fitY[ipl-1])*_planeDist[ipl-1] ;
      //    dth=atan(th2)-atan(th1) ;
      dth= th2 - th1 ;
      chi2 += _planeScat[ipl] * dth * dth;
    }

  // Beam constraint

  if(_useBeamConstraint)
    {
      double dth;
    
      //    dth=atan((_fitX[1]-_fitX[0])*_planeDist[0]) ;
      dth=(_fitX[1]-_fitX[0])*_planeDist[0] ;
      chi2 += _planeScat[0] * dth * dth;

      //    dth=atan((_fitY[1]-_fitY[0])*_planeDist[0]) ;
      dth=(_fitY[1]-_fitY[0])*_planeDist[0] ;
      chi2 += _planeScat[0] * dth * dth;
    }


  return chi2 ;
}



int EUTelTestFitter::GaussjSolve(double *alfa,double *beta,int n)
{
  int *ipiv;
  int *indxr;
  int *indxc;
  int i,j,k;
  int irow=0;
  int icol=0;
  double abs,big,help,pivinv;

  ipiv = new int[n];
  indxr = new int[n];
  indxc = new int[n];

  for(i=0;i<n;ipiv[i++]=0);

  for(i=0;i<n;i++)
    {
      big=0.;
      for(j=0;j<n;j++)
	{
	  if(ipiv[j]==1)continue;
	  for(k=0;k<n;k++)
	    {
	      if(ipiv[k]!=0)continue;
	      abs=fabs(alfa[n*j+k]);
	      if(abs>big)
		{
		  big=abs;
		  irow=j;
		  icol=k;
		}
	    }
	}
      ipiv[icol]++;

      if(ipiv[icol]>1)
	return 1;

      if(irow!=icol)
	{
	  help=beta[irow];
	  beta[irow]=beta[icol];
	  beta[icol]=help;
	  for(j=0;j<n;j++)
	    {
	      help=alfa[n*irow+j];
	      alfa[n*irow+j]=alfa[n*icol+j];
	      alfa[n*icol+j]=help;
	    }
	}
      indxr[i]=irow;
      indxc[i]=icol;

      if(alfa[n*icol+icol]==0.)
	return 1;

      help=alfa[n*icol+icol];
      pivinv=1./help;
      alfa[n*icol+icol]=1.;
      for(j=0;j<n;alfa[n*icol+(j++)]*=pivinv);
      beta[icol]*=pivinv;

      for(j=0;j<n;j++)
	{
	  if(j==icol)continue;
	  help=alfa[n*j+icol];
	  alfa[n*j+icol]=0.;
	  for(k=0;k<n;k++)
	    alfa[n*j+k]-=alfa[n*icol+k]*help;
	  beta[j]-=beta[icol]*help;
	}
    }

  for(i=n-1;i>=0;i--)
    {
      if(indxr[i]==indxc[i])continue;
      for(j=0;j<n;j++)
	{
	  help=alfa[n*j+indxr[i]];
	  alfa[n*j+indxr[i]]=alfa[n*j+indxc[i]];
	  alfa[n*j+indxc[i]]=help;
	}
    }

  delete [] ipiv;
  delete [] indxr;
  delete [] indxc;

  return 0;
}
