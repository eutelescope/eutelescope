/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPROCESSORTRANSFORMFROMGEAR_H
#define EUTELPROCESSORTRANSFORMFROMGEAR_H

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// system includes <>
#include <string>
#include <map>
#include <vector>

//EIGEN includes
#include <Eigen/Core>

namespace eutelescope {


class EUTelProcessorTransformFromGEAR: public marlin::Processor
{
  public:
	virtual Processor* newProcessor() 
	{
		return new EUTelProcessorTransformFromGEAR;
	}

	//! Default constructor
	EUTelProcessorTransformFromGEAR();
	virtual void init();
	virtual void processRunHeader(LCRunHeader * run);
	virtual void processEvent(LCEvent * evt);
	virtual void end();

  private:
	std::map<int, Eigen::Matrix3d> _rotMat;
	std::map<int, Eigen::Vector3d> _offVec;

  protected:
	std::string _inputHitCollectionName;
	std::string _outputHitCollectionName;
	
	int _iRun;
	int _iEvt;

	int _initialOutputCollectionSize;
};

//! A global instance of the processor
 EUTelProcessorTransformFromGEAR gEUTelProcessorTransformFromGEAR;
}
#endif
