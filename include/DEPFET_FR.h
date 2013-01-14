#ifndef DEPFET_FR_H_
#define DEPFET_FR_H_
#include <stdio.h>
#include <string.h>
#include <unistd.h>
//#include "shmem_depfet.h"


namespace DEPFET {

#define MAXDATA 10000   //-- in words; 8195 * 4 = 32780 bytes -> size of one depfet
#define MAXMOD  6

    struct Header {
	unsigned int    EventSize: 20;
	unsigned short   flag0: 1;
	unsigned short   flag1: 1;
	unsigned short  EventType: 2;
	unsigned short  ModuleNo:  4;
	unsigned short  DeviceType: 4;
	unsigned int    Triggernumber;  
    } ;
//.........................................................  
//    unsigned short succframes:       2;
//    unsigned short startgate:       6; //jf new for 128x128 
//    unsigned short succframes:       4;
    struct InfoWord  {
	unsigned int framecnt: 10; // number of Bits
	unsigned int startgate: 10; //jf new for 128x128 
	unsigned int zerosupp: 1;
	unsigned int startgate_ver: 1; 
	unsigned int temperature: 10;                     
    };
//.........................................................  
    struct Hit  {
	unsigned short data: 16; 
	unsigned short col : 6;
	unsigned short row : 7;
	unsigned short nil : 3;    
    };
    
    typedef struct
    { 
	unsigned int HEADER;
	unsigned int Trigger; 
	unsigned int Startgate;
	unsigned int DATA[MAXDATA];
    }
    event_DEPFET;
    
    //! DEPFET File reader
  /*! A class to read and decode a DEPFET file
   *
   *  @author Yulia Furletova, Uni-Bonn <mailto:yulia@mail.cern.ch>
   *  @version $Id$
   */
   
    
    class DEPFET_FR 
    {
	
    private:
	
	FILE *DF;
	int Ntrig;
	int NtrigFlag;
	
    public:
	
	static const int NUM_COL = 64;
	static const int NUM_ROW = 128;
        struct Hit hit;
        struct InfoWord infoword;
        struct Header header;
	struct InfoWord *infoword1;
	
	DEPFET_FR(const char FileName[80],int Ntriggers,int *rc);
	~DEPFET_FR() { };
        void skip_n_words(int n);
	int READ_DEPFET_EVENT(int,int);
	int READ_EVENT_HEADER(int);
	int READ_GROUP_HEADER(int);
	int READ_HEADER(int);
//	int READ_DEPFET_EVENT(int iprint,int DATA1[NUM_COL][NUM_ROW]);
	int READ_DEPFET_EVENT(int iprint,int SkipStartGate, int *DATA1);
	int READ_DEPFET_EVENT1(int iprint, int SkipStartGate, event_DEPFET *depfet_mod,int *eventsize);
	int READ_DEPFET_EVENT_DCD(int iprint, int SkipStartGate, event_DEPFET *depfet_mod,int *eventsize);
	
	
	
	
    };
}
#endif /*DEPFET_FR_H_*/
