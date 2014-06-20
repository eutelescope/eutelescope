#ifndef EUTELPROCESSORMILLEALIGN_H
#define	EUTELPROCESSORMILLEALIGN_H


// C++
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>

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
#include "EUTelRunHeaderImpl.h"
#include "EUTelExceptions.h"


using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

 class  EUTelProcessorMilleAlign : public Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorMilleAlign);     // prevent users from making (default) copies of processors
        
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorMilleAlign;
        }

        EUTelProcessorMilleAlign();
        
    public:
        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init();

        /** Called for every run.
         */
        virtual void processRunHeader(LCRunHeader* run);

        /** Called for every event - the working horse.
         */
        virtual void processEvent(LCEvent * evt);

        virtual void check(LCEvent * evt);

        /** Called after data processing for clean up. **/
	
			  virtual void end();

		protected: 

				std::string _milleBinaryFilename;
				std::string _milleSteeringFilename;
				std::string _milleResultFileName;
				std::string _gear_aligned_file;

};

    /** A global instance of the processor */
    EUTelProcessorMilleAlign gEUTelProcessorMilleAlign;

}

#endif	/* EUTELPROCESSORMILLEALIGN_H */
