#ifndef EUTELSUCIMAIMAGERREADER_H
#define EUTELSUCIMAIMAGERREADER_H

// personal includes ".h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"
//#include "TH2.h"
//#include "TF1.h"
//#include "TFile.h"

// lcio includes <.h>

// system includes <>
namespace eutelescope
{

  class DEPFETReader:public marlin::DataSourceProcessor 
   {
     public:

     //! Default constructor
      DEPFETReader();
     virtual DEPFETReader * newProcessor ();
     virtual void readDataSource (int Ntrig);
     virtual void init ();
     virtual void end ();
 
   protected:
 
       std::string _fileName;
       int _noOfXPixel;
       int _noOfYPixel;
       int _runNumber;
       int _yearflag;
       int _SkipStartGate;
       short *_buffer;
       std::vector<int >_excludePlane;
       std::vector<int >_setupPlane;

   };DEPFETReader  gDEPFETReader;
 
 
}                               // end namespace eutelescope
#endif
