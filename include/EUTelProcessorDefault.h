// LCIO
#include <EVENT/LCCollection.h>
#include "lcio.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

//EUTelescope
#include "EUTelUtility.h"

namespace eutelescope {

	class  EUTelProcessorDefault : public marlin::Processor {

  	private:
  	DISALLOW_COPY_AND_ASSIGN(EUTelProcessorDefault)
        
    public:

    virtual Processor* newProcessor() {
    	return new  EUTelProcessorDefault;
    }

    EUTelProcessorDefault();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every run.*/
    virtual void processRunHeader(LCRunHeader* run);

    /** Called for every event - the working horse.*/
    virtual void processEvent(LCEvent * evt);

    virtual void check(LCEvent * evt);

   	/** Called after data processing for clean up. **/
		virtual void end();

	};

    EUTelProcessorDefault gEUTelProcessorDefault;

}
