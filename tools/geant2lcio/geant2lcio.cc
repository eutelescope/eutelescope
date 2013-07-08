// This test program reads GEANT simuation output from text file
// and write it out in LCIO format
//
// Compile with: 
//  g++ -o geant2lcio geant2lcio.cc -lgsl -lgslcblas -lm -llcio -lsio -lz
//
// after setting proper include path (for LCIO) and lib paths (for LCIO and GSL)
//
// A.F.Zarnecki   March 2007
// updated January 2008 for use with new simulation results
// (backside hits only in the ascii input file)

// system includes

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>

// lcio includes

#include "lcio.h"
#include "IMPL/LCRunHeaderImpl.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/SimTrackerHitImpl.h"
#include "IMPL/LCCollectionVec.h"

// gsl includes

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace lcio;

int main(int argc, char ** argv) {

  // Input parameters are input and putput file names

  if(argc<3){
    cerr << "Parameters missing: [input file] [output file] {geometry file}"  << endl;
    return -1;
  }

  // input, output and geometry file names

  string inputFileName  = argv[1]; 
  string outputFileName = argv[2];
  string geometryFileName = (argc>3)?argv[3]:"geant2lcio.geom";

  cerr << "Converting " << inputFileName.c_str()
       << " to "        <<  outputFileName.c_str() <<  endl;

  cerr << "Using geometry description from " << geometryFileName.c_str() <<  endl;

  // Initialize input stream; exeption handling taken from AB

  ifstream inputFile;
  inputFile.exceptions(ifstream::failbit | ifstream::badbit);
  
  // open the input file
  try {
    inputFile.open(inputFileName.c_str(),ios::in);
  }
  catch (exception& e) {
    cerr << "IO exception " << e.what() << " with " 
	 << inputFileName << ".\nExiting." << endl;
    return -1;
  }
  
  // Open geometry description file

  ifstream geometryFile;
  geometryFile.open(geometryFileName.c_str(),ios::in);


  // Now prepare the output slcio file.
  
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

  // open the file
  try {
    lcWriter->open(outputFileName.c_str(),LCIO::WRITE_NEW);
  }
  catch (IOException& e) {
    cerr << e.what() << endl;
    return -1;
  }

/* create GSL generator chosen by the 
   environment variable GSL_RNG_TYPE */

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

 // Prepare a run header

  int runNumber       = 1;
  int eventNumber     = 0;
  string detectorName = "Eutelescope";
  string detectorDescription = "EUDET telescope, WN-WW configuration";

  LCRunHeaderImpl * runHeader = new LCRunHeaderImpl();
  runHeader->setRunNumber(runNumber);
  runHeader->setDetectorName(detectorName);
  runHeader->setDescription(detectorDescription);

  // Read detector data from geometry description file

  int Ndet,beamaxis,ix,iy,iz;
  double pileup,noise,beamspot;

  geometryFile >> Ndet >> pileup >> noise >> beamaxis >> beamspot ;

// set axis directions for decoding input file
// beamaxis = 1..3

  iz=beamaxis-1;
  ix=(iz+1)%3;
  iy=(iz+2)%3;

  cerr << "Telescope setup with " << Ndet << " layers" <<  endl;

  double * zdet = new double[Ndet];
  double * xmin = new double[Ndet];
  double * xmax = new double[Ndet];
  double * ymin = new double[Ndet];
  double * ymax = new double[Ndet];
  double * resol = new double[Ndet];
  double * effi = new double[Ndet];

  for(int idet=0; idet < Ndet; idet++)
     {
      string detName;
      geometryFile >> zdet[idet] >> xmin[idet] >> xmax[idet] 
                   >> ymin[idet] >> ymax[idet] >> resol[idet] >> effi[idet];
      
      // change resolution units to mm

      resol[idet]/=1000.;

      std::getline(geometryFile,detName,'\n');

     // Add plane names to run header

      runHeader->addActiveSubdetector(detName);
      }

  // Check subdetector list

  const std::vector<std::string> * subDets = 
                               runHeader->getActiveSubdetectors();

  Ndet = subDets->size();

  cerr << Ndet << " subdetectors defined :" << endl;
  for(int idet=0;idet<Ndet;idet++)
    cerr << idet+1 << " : " << subDets->at(idet) << endl;

  // write the header to the output file

  lcWriter->writeRunHeader(runHeader);

  // delete the run header since not used anymore

  delete runHeader;


// Input form Geant

//  int Nsub = 2*Ndet+2;  // front and back side for each sensor plane
  
  int Nsides = 1;            // back side only
  int Nsub = Nsides*Ndet+2;   
  double rawpos[3];

// Beamspot offset

  double offset[3]; 

// MC point tables:

  double * xgen = new double[Ndet];
  double * ygen = new double[Ndet];
  double * zgen = new double[Ndet];

// Measured points tables:

  double * xmes = new double[Ndet];
  double * ymes = new double[Ndet];
  double * zmes = new double[Ndet];
  bool  * fired = new bool[Ndet];

// Event loop

  while ( true ) {
    


    // Prepare event header and collections 

     eventNumber++;

        if (eventNumber%1000 == 0) 
              cout << "Converting event number " << eventNumber << endl;

     LCEventImpl * event = new LCEventImpl();

        event->setDetectorName(detectorName);

        event->setRunNumber(runNumber);

        event->setEventNumber(eventNumber);


    // prepare a collection to store generated points 

    LCCollectionVec     * simhitvec = new LCCollectionVec(LCIO::SIMTRACKERHIT);

    // prepare a collection to store measured points

    LCCollectionVec     * meshitvec = new LCCollectionVec(LCIO::TRACKERHIT);


   // Loop over geant events to simulate pileup

   // Number of tracks to include

    unsigned int npile=1;

    npile += gsl_ran_poisson(r,pileup);

    for(int ipile=0; ipile<npile;ipile++)
      {
      int nread=0;

// Beam spot

   offset[0]=gsl_ran_gaussian(r,beamspot);
   offset[1]=gsl_ran_gaussian(r,beamspot);
   offset[2]=0.;

// Read one Geant event from file

    for(int idet=0;idet<Ndet;idet++)
             {
             xgen[idet]=offset[0];
             ygen[idet]=offset[1];
             zgen[idet]=offset[2];
             }

    try{ 
      
      for(int isub=0 ; isub <  Nsub ; isub++ )
        {
// Geant4 input: beam along Z axis
	  inputFile >> rawpos[2] >> rawpos[0] >> rawpos[1] ;

// Input file is in um -> convert to mm

          int idet=(isub-1)/Nsides;

          if(isub>0 && idet<Ndet)
             {
             xgen[idet]+=rawpos[0]/1000./Nsides;
             ygen[idet]+=rawpos[1]/1000./Nsides;
             zgen[idet]+=rawpos[2]/1000./Nsides;
             }

	  ++nread ;
	}

    }      
    catch(exception& e){

      if( !inputFile.eof() ) 
	cerr << " a read exception occured : " << e.what() << endl  ;

      if( nread != 0 && nread != Nsub ){
	cerr << " less than " << Nsub << " points read - event is incomplete ! " 
	     << endl ;
      }
      break ;
    }    

// calculate measured points

    for(int idet=0;idet<Ndet;idet++)
        {

       // Apply Gaussian smearing

        xmes[idet]=xgen[idet]+gsl_ran_gaussian(r,resol[idet]);
        ymes[idet]=ygen[idet]+gsl_ran_gaussian(r,resol[idet]);
        zmes[idet]=zgen[idet];


        fired[idet]=true;

        // check detector range

        if(xmes[idet]<xmin[idet] || xmes[idet]>xmax[idet])fired[idet]=false;
        if(ymes[idet]<ymin[idet] || ymes[idet]>ymax[idet])fired[idet]=false;

        // Detector efficiency

        if(gsl_rng_uniform(r)>effi[idet])fired[idet]=false;

        }


// Add generated points to collection

    for(int idet=0;idet<Ndet;idet++)
        {

    // Fill and store single MC points 
    // Each point has to remain in memory until event is written out !

        SimTrackerHitImpl * simhit = new SimTrackerHitImpl;

    // Cell ID is just the plane number  ID=1...Ndet
        
        simhit->setCellID(idet+1);

    // Get position, change axis according to the beam direction 

        double pos[3];

        pos[ix]=xgen[idet];
        pos[iy]=ygen[idet];
        pos[iz]=zgen[idet];

        simhit->setPosition(pos);

     // store this hit

       simhitvec->push_back(simhit);

        }

// Add measured points to collection

    for(int idet=0;idet<Ndet;idet++)
        {

       // take only valid hits:

        if(fired[idet])
            {
       // Fill and store measured points 
       // Each point has to remain in memory until event is written out !

            TrackerHitImpl * meshit = new TrackerHitImpl;

            if(meshit == NULL) 
                  cerr << "TrackerHitImpl allocation failed" << endl;

       // Store plane number  as hit type:

             meshit->setType(idet+1);

      // store measured position

             double pos[3];

             pos[ix]=xmes[idet];
             pos[iy]=ymes[idet];
             pos[iz]=zmes[idet];

             meshit->setPosition(pos);

      // Covariance matrix of the position 
      // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).

            float cov[TRKHITNCOVMATRIX];

            cov[ix+ix*ix]=resol[idet]*resol[idet];
            cov[iy+iy*iy]=resol[idet]*resol[idet];
            cov[iz+iz*iz]=0.;

            cov[1]=cov[3]=cov[4]=0.;

            meshit->setCovMatrix(cov);


      // store measured hit 

             meshitvec->push_back(meshit);

             }       // if(fired[idet])
       }    // idet loop


// Just in case of EOF

  if (inputFile.eof()) break;

      }    // pileup event loop


// Add noise:

    for(int idet=0;idet<Ndet;idet++)
        {

        unsigned int nnoise=gsl_ran_poisson(r,noise);

        for(int inoise=0;inoise<nnoise;inoise++)
            {
       // Fill and store noise points 

             TrackerHitImpl * meshit = new TrackerHitImpl;

       // Store plane number  as hit type:

             meshit->setType(idet+1);

      // noise position : uniform close to true hit

             double pos[3];

             pos[ix]=gsl_rng_uniform(r)*(xmax[idet]-xmin[idet])+xmin[idet];
             pos[iy]=gsl_rng_uniform(r)*(ymax[idet]-ymin[idet])+ymin[idet];
             pos[iz]=zdet[idet];

             meshit->setPosition(pos);

      // Covariance matrix of the position 
      // (stored as lower triangle matrix, i.e.  cov(xx),cov(y,x),cov(y,y) ).

            float cov[TRKHITNCOVMATRIX];

            cov[ix+ix*ix]=resol[idet]*resol[idet];
            cov[iy+iy*iy]=resol[idet]*resol[idet];
            cov[iz+iz*iz]=0.;

            cov[1]=cov[3]=cov[4]=0.;

            meshit->setCovMatrix(cov);


      // store measured hit 

             meshitvec->push_back(meshit);

             }       // noise loop

        }    // idet loop




    //  if(npile>2)cerr << "Pileup of " << npile << " MC particles in event " << eventNumber << endl;


// add hit collections to an event

    event->addCollection(simhitvec,"simhit");


// add measured hit collection to an event

    event->addCollection(meshitvec,"meshit");


    // write event out

    lcWriter->writeEvent(event);
    
    // clean
    // deleting an event also delets everything what was put into this event...

    delete event;

   // if the file is at the end ... break the loop!

  if (inputFile.eof()) break;

// End of event loop
}


// That is all! Close all streams...

  gsl_rng_free (r);
  
  lcWriter->close();
  inputFile.close();
  geometryFile.close();

 
  return 0;
}

