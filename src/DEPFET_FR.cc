#include "DEPFET_FR.h"

using namespace DEPFET;

DEPFET_FR::DEPFET_FR(const char FileName[80],int Ntriggers,int *rc)
{
  *rc=0;
  printf("DEPFET_FR::FileName=%s\n",FileName);
  if((DF = fopen(FileName,"r")) == NULL)
    { printf("Can not open file %s \n",FileName);  *rc=-1; }
  Ntrig=Ntriggers;
  if(Ntrig>=0) { NtrigFlag=1;} else {NtrigFlag=0;}
  if(Ntrig>0)  printf("Trigers to run =%d \n",Ntrig);
  *rc=1;
}

/*=====================================================================*/
int DEPFET_FR::READ_HEADER(int iprint){
  aaa:    fread(&header,sizeof(header),1,DF);
  if(iprint) { 
    printf( "header.DeviceType=0x%x\n",header.DeviceType);
    printf( "header.EventSize=%d\n",header.EventSize);
    printf( "header.Triggernumber=%d\n",header.Triggernumber);
    printf( "header.ModuleNo=%d\n",header.ModuleNo);
    printf( "header.EventType=0x%x\n",header.EventType);
  }
  if( header.DeviceType ==10 || header.DeviceType==11 ||header.DeviceType==5 || header.DeviceType==14 || header.DeviceType==13) {
    fseek(DF,(header.EventSize-2)*4,SEEK_CUR);
    if(iprint)printf("skip devices: %d  ModId=%d\n",header.DeviceType,header.ModuleNo);
    goto aaa;
  }

  return 0;
}
/*=====================================================================*/
int DEPFET_FR::READ_GROUP_HEADER(int iprint){
  int Nmodtot = 0;
  printf("==================================\n");

  READ_HEADER(iprint);
  if(header.EventType==0x0) {
    Nmodtot=int ((header.EventSize-2)/2) ;
    printf("GROUP HEADER:: BORE :: Nmodules=%d \n",Nmodtot);
  } else if (header.EventType==0x2) {

    return -1;
  };


  printf("GROUP HEADER:: BORE done\n");
  printf("==================================\n");

  return Nmodtot;
}
/*=====================================================================*/
int DEPFET_FR::READ_EVENT_HEADER(int iprint) {
  static int first=1;

  READ_HEADER(iprint);

  switch (header.DeviceType) {
//...................................................
  case (0x0):
    //  printf("Group event \n");
    if(header.Triggernumber%10==0&& !first) {
      // printf( "Triggernumber=%d\n",header.Triggernumber);
    }
    // printf("ModTot=%d\n",Nmodtot);
    if(first) { first=0;
//              BOOK_HIST();
    }
    if(NtrigFlag==1 && header.Triggernumber>(unsigned)Ntrig) return -1;
    return 100;
    break;
//....................................................
  case (0x2):
    //printf("DEPFET event \n");
    if(header.EventType==0x0) {
      printf("Mod %d BORE\n",header.ModuleNo);}

    else if(header.EventType==0x1){
      printf("Mod %d EORE\n",header.ModuleNo); return 1;}
    else if(header.EventType==0x2){

      if(iprint) printf("Mod %d data\n",header.ModuleNo); return 2;}

    break;
//....................................................
  case (0x3): // DEPFET Event 128x128
    //printf("DEPFET event \n");
    if(header.EventType==0x0) { 
       printf("Mod %d BORE\n",header.ModuleNo);}
      
    else if(header.EventType==0x1){
       printf("Mod %d EORE\n",header.ModuleNo); return 1;}
    else if(header.EventType==0x2){
	
         if(iprint) printf("Mod %d data\n",header.ModuleNo); return 2;}
    break;

//....................................................
  case(0xd):// TLU event
    if(header.EventType==0x2 ) { printf("TLU EVENT\n");
    fseek(DF,(header.EventSize-2)*4,SEEK_CUR);

    return 3;
    };
    break;
//....................................................
  case(0xe):// TLU event
    fseek(DF,(header.EventSize-2)*4,SEEK_CUR);

    return 4;

    break;

//....................................................

  }// end switch



  return 0;
}



/*=====================================================================*/
/*=====================================================================*/
int  DEPFET_FR::READ_DEPFET_EVENT(int iprint,int SkipStartGate,int DATA1[NUM_COL][NUM_ROW] ){
  int rc,i;
  int icol,irow,startgate = 0;

   unsigned int Startgate;
    struct InfoWord *infoword1;
    struct InfoWord infoword2;

    rc=fread(&Startgate,sizeof(Startgate),1,DF);
    //  printf("unsign =%d struct =%d \n",sizeof(Startgate),sizeof(infoword2));
    //printf("Infoword2 =%d \n",sizeof(struct InfoWord));
    //  rc=fread(&infoword1,sizeof(infoword1),1,DF);
   infoword1=(struct InfoWord *)&Startgate;

//  rc=fread(&infoword,sizeof(infoword),1,DF);
  if(rc==0) {printf("EOF \n"); return 1;};
  if(iprint) { 
    printf("InfoWord.framecnt=0x%x\n",infoword1->framecnt);
    printf("InfoWord.startgate=0x%x\n",infoword1->startgate);
    printf("InfoWord.succframes=0x%x\n",infoword1->succframes);
    printf("InfoWord.zerosupp=0x%x \n",infoword1->zerosupp);
    printf("InfoWord.startgate_ver=0x%x\n",infoword1->startgate_ver);
    printf("InfoWord.temperature=0x%x\n",infoword1->temperature);
  };

  int eventsize =  header.EventSize-3 ;
  if(iprint) printf("DEPFET_FR:: eventsize=%d \n",eventsize);

    Hit * depfet_hits = new Hit[eventsize];

  rc=fread(depfet_hits,sizeof(Hit)*eventsize,1,DF);
  // printf("DEPFET_FR:: rc=%d \n",rc);
  if(infoword.zerosupp) {
    for (icol=0;icol<NUM_COL;icol++) {
      for (irow=0;irow<NUM_ROW;irow++) {
        DATA1[icol][irow]=-10000;
      }
    }
  }

  if(iprint) {printf("DEPFET_FR:: Read depfet_hits done \n");} 
  for (i=0;i<eventsize;i++) {
      Hit *tmp;
      tmp=&depfet_hits[i];
      //     memcpy(&tmp,&(depfet_hits[i]),sizeof(Hit));
             if (i==0) startgate=depfet_hits[i].row;
      //if (i==0) startgate=(*tmp).row;
        int z= depfet_hits[i].data & 0x0000ffff;
	//   printf("z=%d\n",z);

//     DATA[depfet_hits[i].col][depfet_hits[i].row]=z;
    DATA1[depfet_hits[i].col][depfet_hits[i].row]=z;

  };
  if(iprint) {printf("DEPFET_FR:: decoding done \n");} 
  if( SkipStartGate) {
    if (startgate/2<62) {
      for(icol=0;icol<NUM_COL;icol++) {
        for (irow=0;irow<4;irow++) {
          DATA1[icol][startgate+irow]=15000;
        };
      };
    } else {startgate=0;};
  };

  delete [] depfet_hits;

  return 0;
}

/*=====================================================================*/
