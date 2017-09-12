#ifndef EUTELEXTERNALTRIGGER_H
#define EUTELEXTERNALTRIGGER_H

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
    EUTelExternalTrigger();

    // Default destructor
    ~EUTelExternalTrigger() { ; }

    /*
      Get time stamp of trigger
     */
    
    long long unsigned getTimestamp() const { return _timestamp; }
    
    /*
      Get trigger label
     */
    short getLabel() const { return _label; }

    /*
      Set time stamp of trigger
     */
    void setTimestamp(long long unsigned timestamp) { _timestamp = timestamp; }
    
    /*
      Set trigger label
     */
    void setLabel(short label) { _label = label; }

    /*
      Get number of elements in lcio data structure
     */
    unsigned int GetNoOfElements() const { return _nElement; }

    
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
