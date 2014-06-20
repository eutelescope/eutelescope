
#include "EUTelProcessorMilleAlign.h"

using namespace eutelescope;

EUTelProcessorMilleAlign::EUTelProcessorMilleAlign() :
Processor("EUTelProcessorGBLFitCandidates"),
_milleBinaryFilename("mille.bin"),
_milleSteeringFilename("pede-steer.txt"),
_milleResultFileName("millepede.res"),
_gear_aligned_file("gear-00001-aligned.xml"){

	registerOptionalParameter("MilleBinaryFilename", "Name of the Millepede binary file", _milleBinaryFilename, std::string("mille.bin"));

	registerOptionalParameter("MilleSteeringFilename", "Name of the Millepede steering file to be created", _milleSteeringFilename, std::string("pede-steer.txt"));
    
	registerOptionalParameter("MilleResultFilename", "Name of the Millepede result file", _milleResultFileName, std::string("millepede.res"));
    
	registerOptionalParameter("GearAlignedFile", "Suffix to add to the new Gear with alignment corrections", _gear_aligned_file, std::string("gear-00001-aligned.xml"));




}

void EUTelProcessorMilleAlign::init() {}

void EUTelProcessorMilleAlign::processRunHeader(LCRunHeader * run) {}

void EUTelProcessorMilleAlign::processEvent(LCEvent * evt){}

void EUTelProcessorMilleAlign::check(LCEvent * evt){}

void EUTelProcessorMilleAlign::end(){}
