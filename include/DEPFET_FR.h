#ifndef DEPFET_FR_H_
#define DEPFET_FR_H_
#include <stdio.h>
#include <string.h>
#include <unistd.h>

struct Header {
    unsigned int    EventSize: 20;
     unsigned short   flag0: 1;
     unsigned short   flag1: 1;
    unsigned short  EventType: 2;
    unsigned short  ModuleNo: 4;
    unsigned short  DeviceType: 4;
    unsigned int    Triggernumber;  
} ;
//.........................................................  
struct InfoWord  {
    unsigned short framecnt:        10; // number of Bits
    unsigned short startgate:        6;
    unsigned short succframes:       4;
    unsigned short zerosupp:         1;
    unsigned short startgate_ver:    1; 
    unsigned short temperature:     10;                     
};
//.........................................................  
struct Hit  {
    unsigned short data: 16; 
    unsigned short col : 6;
    unsigned short row : 7;
    unsigned short nil : 3;    
};



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


	DEPFET_FR(const char FileName[80],int Ntriggers,int *rc);
	~DEPFET_FR() { };

       void skip_n_words(int n);
 	int READ_DEPFET_EVENT(int,int);
	int READ_EVENT_HEADER(int);
	int READ_GROUP_HEADER(int);
	int READ_HEADER(int);
	int READ_DEPFET_EVENT(int iprint,int SkipStartGate,int DATA1[NUM_COL][NUM_ROW]);

	//    int READ_DEPFET_EVENT1(int iprint,event_DEPFET *depfet_mod,int *eventsize);

};
#endif /*DEPFET_FR_H_*/
