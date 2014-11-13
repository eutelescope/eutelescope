#include "EUTelProcessorDefault.h"

using namespace eutelescope;

EUTelProcessorDefault::EUTelProcessorDefault() :
Processor("EUTelProcessorDefault"){}


void EUTelProcessorDefault::init(){}

void EUTelProcessorDefault::processRunHeader(LCRunHeader * run) {}

void EUTelProcessorDefault::check(LCEvent * evt){}

void EUTelProcessorDefault::processEvent(LCEvent * evt){}

void EUTelProcessorDefault::end(){}
