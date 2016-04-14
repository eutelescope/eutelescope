
/*
  Basic class for external triggers provided 
  in addition to TLU triggers
  
  @auther Dorothea vom Bruch, <mailto: vombruch@uni-mainz.de>
 */


// personal includes
#include "EUTelExternalTrigger.h"

using namespace eutelescope;

EUTelExternalTrigger::EUTelExternalTrigger() : 
  _timestamp(0),
  _label(-1),
  _nElement(3)
{
  
}

EUTelExternalTrigger::EUTelExternalTrigger(long timestamp, short label) :
  _timestamp(timestamp),
  _label(label),
  _nElement(3)
{

}
