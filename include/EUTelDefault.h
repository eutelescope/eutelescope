// LCIO
#include <EVENT/LCCollection.h>
#include "lcio.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;
using namespace std;


namespace eutelescope {

	class  EUTelProcessorDefault : public Processor {

  	private:
  	DISALLOW_COPY_AND_ASSIGN( EUTelProcessorDefault);      // prevent users from making (default) copies of processors
        
    public:

    virtual Processor* newProcessor() {
    	return new  EUTelProcessorDefault;
    }

    EUTelProcessorDefault();

	}

    EUTelProcessorDefault gEUTelProcessorDefault;

}
