#ifdef USE_GBL    // Not sure where this is defined. However it is used in all the other GBL processors will use it.

// C++
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>

// LCIO
#include <EVENT/LCCollection.h>

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

// AIDA
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IProfile2D.h>
#endif // MARLIN_USE_AIDA

#include "EUTelProcessorGBLFitCandidates.h"

EUTelProcessorGBLFitCandidates::EUTelProcessorGBLFitCandidates() :
Processor("EUTelProcessorGBLFitCandidates"){}

void EUTelProcessorGBLFitCandidates::init() {}

void EUTelProcessorGBLFitCandidates::processRunHeader(LCRunHeader * run) {}

void check(LCEvent * evt){}

void EUTelProcessorGBLFitCandidates::processEvent(LCEvent * evt){}

void EUTelProcessorTrackingGBLTrajectory::end() {]

#endif // USE_GBL
