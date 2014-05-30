//Writen by Alexander Morton <alexander.morton@desy.de>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

//Can only use this processor if AIDA is defined.
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

//Histogram header files
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram2D.h>

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

//Deals with all the geometry of the telescope
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "EUTelProcessorHitPlotter.h"

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;


//Define constructor. This will take collection name as input 
EUTelProcessorHitPlotter::EUTelProcessorHitPlotter() : Processor("EUTelProcessorHitPlotter"),
_hitCollectionNameInput(){ //Variables here initialized with the constructor of this class

	// modify processor description
	_description ="EUTelLocaltoGlobalHitMaker is responsible to change local coordinates to global. This is done using the GEAR geometry description";

  registerInputCollection(LCIO::TRACKERHIT,"hitCollectionNameInput",
                           "Hit collection name",
                           _hitCollectionNameInput, string ( "" ));
}

void EUTelProcessorHitPlotter::init() {}

void EUTelProcessorHitPlotter::processRunHeader (LCRunHeader * rdr) {}

void EUTelProcessorHitPlotter::check(LCEvent *event){}

void EUTelProcessorHitPlotter::processEvent (LCEvent * event) {}

void EUTelProcessorHitPlotter::end(){}






#endif
