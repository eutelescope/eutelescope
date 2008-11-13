#ifndef Mod_DEPF_A
#define Mod_DEPF_A 1
#include <stdio.h>
#include <string.h>
#include <unistd.h>


namespace DEPFET {

  class Mod_DEPF{

    int ID;

  public:

    static const int NUM_COL = 64;
    static const int NUM_ROW = 128;

    static const int PITCHX = 33;
    static const int PITCHY = 22;

    struct InfoWord DEPF_InfoWord;
    struct Header DEPF_Header;
    struct Hit DEPF_Hit;
    int DATA1[NUM_COL][NUM_ROW];
    /*............................................*/
    void fillDATA(int icol,int irow, int idata) {
      DATA1[icol][irow]=idata;
    }
    /*............................................*/
    int ModID(){ return ID; }
    /*............................................*/

    /*............................................*/
    Mod_DEPF(int ii) {
      ID=ii;
      printf("Mod_DEPF:: ModID=%d \n",ID);
    }

    /*............................................*/
    ~Mod_DEPF() {}
    /*............................................*/
  };

}
#endif
