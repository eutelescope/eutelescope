/*=====================================================================*/
/*          DEPFET file converter (RAW->LCIO)                          */
/*          Author: Julia Furletova                                    */
/*                (julia@mail.desy.de or furletova@physik.uni-bonn.de) */
/*          Created   02 nov 2007                                      */
/*          Modified                                                   */
/*=====================================================================*/

// user includes
#include "DEPFET_FR.h"
#include "Mod_DEPF.h"
#include "DEPFETReader.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"


// marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// lcio includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
// #include <UTIL/LCTOOLS.h>

// system includes
#include <fstream>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

using namespace std;
using namespace marlin;
using namespace DEPFET;
using namespace eutelescope;


/*=====================================================================*/
DEPFETReader::DEPFETReader ():DataSourceProcessor  ("DEPFETReader") {
/*=====================================================================*/
 #define SETUP_MAX 7

  _description =
    "Reads data files and creates LCEvent with TrackerRawData collection.\n"
    "Make sure to not specify any LCIOInputFiles in the steering in order to read DEPFET files.";


  IntVec excludePlanes;

  IntVec setupPlanes;
  for(int j=0;j<SETUP_MAX;j++){
      excludePlanes.push_back(0);
      setupPlanes.push_back(0);
  }


  registerProcessorParameter ("FileName", "Input file",
                              _fileName, std::string ("input.dat"));
  registerProcessorParameter ("runNumber", "RunNumber",
                              _runNumber, static_cast < int >(01));

  registerProcessorParameter ("NoOfXPixel", "Number of pixels along X",
                              _noOfXPixel, static_cast < int >(64));
  registerProcessorParameter ("NoOfYPixel", "Number of pixels along Y",
                              _noOfYPixel, static_cast < int >(128));
  registerProcessorParameter ("Excludeplane", "Exclude one or more planes: 9 14 1 4 (module ID)",
                              _excludePlane,excludePlanes);
  registerProcessorParameter ("SETUP", "Setup starting from beam 9 14 1 4 (module ID)",
                              _setupPlane,setupPlanes);
  registerProcessorParameter ("YEAR", " year",
                              _yearflag,static_cast < int >(07));

  registerProcessorParameter ("SkipStartGate",
                              "Skip 4 rows after StartGate",
                              _SkipStartGate, static_cast < int >(0));


}

/*=====================================================================*/
 DEPFETReader * DEPFETReader::newProcessor () {
/*=====================================================================*/

   return new DEPFETReader;
}

/*=====================================================================*/
 void DEPFETReader::init () {
/*=====================================================================*/

    printParameters ();
}

/*=====================================================================*/
 void DEPFETReader::readDataSource (int Ntrig) {
/*=====================================================================*/
    EUTelEventImpl *event = NULL;
    LCCollectionVec *rawData = NULL;

    TrackerRawDataImpl *rawMatrix[20];
    int rc,eventNumber;
    static int First=1;
    char filename[80];
    int Nmodtot,iprint=0,i,idet;
    int Nex,Ndet,kk;
    int DATA[64][128];
    Mod_DEPF *DEPFET[10];
    DEPFET_FR  *myfile;
    static int MOD2ID[32], ID2MOD[32];

    sprintf(filename,"input.dat");
    myfile=new DEPFET_FR(filename,Ntrig,&rc);
    printf("RC=%d\n",rc);
    if (rc<0) {printf("problem with opening file\n"); exit (-1); };
   /*=====================================================================*/
   /*                     READ Group HEADER                               */
   /*=====================================================================*/
    if(_yearflag==9 ) {Nmodtot=6; goto year09a;}

     Nmodtot=myfile->READ_GROUP_HEADER(iprint);
 //    int _runNumber=1802;
     if(_yearflag==06 ) Nmodtot=5;
     // if(_yearflag==07) Nmodtot=2;
    year09a:;
     Nex=0; Ndet=0;
     for(i=0;i<SETUP_MAX;i++) {
         if (_setupPlane[i]<=0) break;
         ID2MOD[_setupPlane[i]]=-1;
         printf("-> ALL::  i=%d  Mod=%d  iMOD=%d \n"
                ,i,_setupPlane[i],ID2MOD[i]);
         for(int ii=0;ii<SETUP_MAX;ii++) {
             if(_excludePlane[ii]==_setupPlane[i]) { Nex+=1; goto excl; }
         }
         ID2MOD[Ndet]=_setupPlane[i];
         MOD2ID[_setupPlane[i]]=Ndet;
         printf("   SET::  i=%d  Mod=%d MOD2ID=%d iMOD=%d \n"
                ,i,_setupPlane[i],MOD2ID[_setupPlane[i]], ID2MOD[Ndet]);

         DEPFET[Ndet] = new Mod_DEPF(ID2MOD[Ndet]);
         Ndet++;
     excl:;
     }
     printf (" DEPFET  Init:: Ndet=%d Nex=%d\n",Ndet, Nex);
     if(_yearflag==9) goto year09b;
     kk=0;
     for( i=0;i<Nmodtot;i++) {
         myfile->READ_HEADER(iprint);
         printf("DEPFET_FR :: %d \n",  myfile->header.ModuleNo);
         if (MOD2ID[myfile->header.ModuleNo]<0) continue;
         kk++;
         printf("MMMMM=%d \n",myfile->header.ModuleNo);
         printf(".i=%d, kk=%d, id2mod=%d, Modid=%d, mod2id=%d........\n"
                ,i,kk,ID2MOD[kk], myfile->header.ModuleNo,MOD2ID[myfile->header.ModuleNo]);
     }
 year09b:;
     printf ("================================== \n");
     printf ("Number of Modules %d , Run Number=%d\n",Nmodtot,_runNumber);
     printf ("================================== EOF HEADER\n");
     eventNumber=0;

   /*=====================================================================*/
   /*                 READ DATA                                           */
   /*=====================================================================*/

     while (true)  {  idet=0; kk=-1;

     if(iprint) printf ("================================== START READ DATA\n");

     rc=myfile->READ_EVENT_HEADER(iprint); //---read 2 words from event(Header+Trigger)

     /*-------------Nevents=Ntrig...go to end---------*/
     if(rc==-1) goto aaa;
     /*-------------Event Group HEADER----------------*/
     if(rc==100)  { // --goto skip_1;

     if(!First) {
         event->addCollection (rawData, "rawdata");
         ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (event));
         delete event;
     }
     /*---------BOR event---------*/
     else { First=0;
         printf("FIRST event:: start fill HEADER\n");
         // in the case it is the first run so we need to process and
         // write out the run header.

         auto_ptr<IMPL::LCRunHeaderImpl> lcHeader  ( new IMPL::LCRunHeaderImpl );
         auto_ptr<EUTelRunHeaderImpl>    runHeader ( new EUTelRunHeaderImpl (lcHeader.get()) );
         runHeader->addProcessor( type() );
         runHeader->lcRunHeader()->setDescription(" Events read from DEPFET input file: " + _fileName);
         runHeader->lcRunHeader()->setRunNumber (_runNumber);
         runHeader->setHeaderVersion (0.0001);
         runHeader->setDataType (EUTELESCOPE::CONVDATA);
         runHeader->setDateTime ();
         runHeader->addIntermediateFile (_fileName);
         runHeader->addProcessor (_processorName);
         // this is a mistake here only for testing....
         runHeader->setNoOfEvent(Ntrig);
         //////////////////////////////////////////////
         //////////////////////////////////////////////
         runHeader->setNoOfDetector(Ndet);
         runHeader->setMinX(IntVec(Ndet, 0));
         runHeader->setMaxX(IntVec(Ndet, _noOfXPixel - 1));
         runHeader->setMinY(IntVec(Ndet, 0));
         runHeader->setMaxY(IntVec(Ndet, _noOfYPixel - 1));

         runHeader->lcRunHeader()->setDetectorName("DEPFET");
         //UTIL::LCTOOLS::dumpRunHeader(runHeader);

         // process the run header
         ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> ( lcHeader.release()) );

         // end of first event
         _isFirstEvent = false;
         printf("END BOR event:: \n");

     };    /*---------END BOR event---------*/

     event = new EUTelEventImpl;
     event->setDetectorName("DEPFET");
     event->setEventType(kDE);
     LCTime * now = new LCTime;
     event->setTimeStamp(now->timeStamp());
     delete now;
     rawData = new LCCollectionVec (LCIO::TRACKERRAWDATA); //<---- Only for RAW data

     eventNumber++;
     if(eventNumber%10==0) printf("Event number=%d  Trigger=%d \n",eventNumber,myfile->header.Triggernumber);
     if(eventNumber%16==0) { printf("SKIP Event number=%d  Trigger=%d \n",eventNumber,myfile->header.Triggernumber); goto myskip;}
//       eventNumber=myfile->header.Triggernumber;
     if(iprint) printf("EventNumber=%d\n",eventNumber);
     event->setRunNumber (_runNumber);
     event->setEventNumber (eventNumber);


     } // skip_1:  //--  if rc!=100  --> end group header

     kk=0;
     do { idet++;
         if(iprint) printf ("===============> idet= %d  rc=%d\n",idet,rc);
         rc=myfile->READ_EVENT_HEADER(iprint);
         if(iprint)printf("READ_DATA rc=%d \n",rc);
         /*-------------Module DATA HEADER----------------*/
         if(rc==1) goto aaa;
         if(rc==2) {   // goto skip_2;

         if(myfile->READ_DEPFET_EVENT(iprint,_SkipStartGate,DATA)==1) goto aaa;
         int iMOD=MOD2ID[myfile->header.ModuleNo];
         if (iMOD<0) continue;
         if (iprint)
             printf("Start %d %d,%d Trig=%d\n",iMOD, DEPFET[iMOD]->ModID(),myfile->header.ModuleNo,myfile->header.Triggernumber);

         rawMatrix[iMOD] = new TrackerRawDataImpl;
         CellIDEncoder < TrackerRawDataImpl > idEncoder (EUTELESCOPE::MATRIXDEFAULTENCODING, rawData);
//       idEncoder["sensorID"] = DEPFET[iMOD]->ModID();
         idEncoder["sensorID"] = iMOD;
         idEncoder["xMin"] = 0;
         idEncoder["xMax"] = _noOfXPixel - 1;
         idEncoder["yMin"] = 0;
         idEncoder["yMax"] = _noOfYPixel - 1;
         idEncoder.setCellID (rawMatrix[iMOD]);


         for (int yPixel = 0; yPixel < _noOfYPixel; yPixel++) {
             for (int xPixel = 0; xPixel < _noOfXPixel; xPixel++) {
                 rawMatrix[iMOD]->adcValues ().push_back (DATA[xPixel][yPixel]);
                 if(iprint==6)printf("Det=%d ADCvalue=%d,  ipix=%d %d,\n",i,DATA[xPixel][yPixel],xPixel,yPixel);
             }
         }
//       rawData->push_back (rawMatrix[iMOD]);

         }  // skip_2:  //end eventheader rc==2

     myskip:;
     } while (idet!=Nmodtot);
     printf("-------------------------------\n");
     for (int ii=0;ii<Ndet;ii++)  rawData->push_back (rawMatrix[ii]);

     }; // end while(true)
 aaa:
     event = new EUTelEventImpl;
     event->setDetectorName("DEPFET");
     LCTime * now = new LCTime;
     event->setTimeStamp(now->timeStamp());
     delete now;
     event->setRunNumber (_runNumber);
     event->setEventNumber (eventNumber++);
     event->setEventType(kEORE);
     ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (event));
     delete event;
 }


/*=====================================================================*/

 void DEPFETReader::end () {
   message<MESSAGE> ("Successfully finished") ;
 }

/*=====================================================================*/


