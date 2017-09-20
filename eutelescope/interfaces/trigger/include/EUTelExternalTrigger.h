#ifndef EUTELEXTERNALTRIGGER_H
#define EUTELEXTERNALTRIGGER_H

#include "EUTELESCOPE.h"

namespace eutelescope {

  /* Basic class for external triggers provided 
     in addition to TLU triggers
     @auther Dorothea vom Bruch, <mailto: vombruch@uni-mainz.de>
  */

  class EUTelExternalTrigger {

  public:
  
    // Default constructor with all arguments
    EUTelExternalTrigger(long long unsigned timestamp, short label);
    
    // Default constructor without arguments
    EUTelExternalTrigger() = default;

    // Default destructor
    ~EUTelExternalTrigger() { ; }

    /*
      Get time stamp of trigger
     */
    
    long long unsigned getTimestamp() const; 
    
    /*
      Get trigger label
     */
    short getLabel() const;

    /*
      Set time stamp of trigger
     */
    void setTimestamp(long long unsigned timestamp);
    
    /*
      Set trigger label
     */
    void setLabel(short label);

    /*
      Get number of elements in lcio data structure
     */
    unsigned int GetNoOfElements() const;

    
  private:
  
  //! time stamp of trigger
  long long unsigned _timestamp;
  
  //! trigger label
  short _label;
  
  //! number of 32-bit wide elements in data structure
  unsigned int _nElement;

  };

}

#endif
