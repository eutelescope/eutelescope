#ifndef EUTELESCOPEPROCESSORGBLFITCANDIDATES_H
#define	EUTELESCOPEPROCESSORGBLFITCANDIDATES_H

#ifdef USE_GBL



// LCIO
#include "lcio.h"

#include "marlin/Processor.h"

using namespace lcio;
using namespace marlin;
using namespace std;

namespace eutelescope {

 class EUTelProcessorGBLFitCandidates : public Processor {

    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelProcessorGBLFitCandidates)     // prevent users from making (default) copies of processors
        
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorGBLFitCandidates;
        }

        EUTelProcessorGBLFitCandidates();
        
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

        /** Called after data processing for clean up.


#endif // USE_GBL

#endif	/* EUTELESCOPEPROCESSORTRACKINGGBLTRACKFIT_H */
