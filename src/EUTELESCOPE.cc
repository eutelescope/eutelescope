// Author: Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id: EUTELESCOPE.cc,v 1.5 2007-02-22 08:09:36 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#include "EUTELESCOPE.h"

using namespace eutelescope;

const char *   EUTELESCOPE::HEADERVERSION       = "HeaderVersion";
const char *   EUTELESCOPE::NOOFEVENT           = "NoOfEvent";
const char *   EUTELESCOPE::DATATYPE            = "DataType";
const char *   EUTELESCOPE::DATETIME            = "DateTime";
const char *   EUTELESCOPE::DAQHWNAME           = "DAQHWName";
const char *   EUTELESCOPE::DAQHWVERSION        = "DAQHWVersion";
const char *   EUTELESCOPE::DAQSWNAME           = "DAQSWName";
const char *   EUTELESCOPE::DAQSWVERSION        = "DAQSWVersion";
const char *   EUTELESCOPE::SIMULSWNAME         = "SimulSWName";
const char *   EUTELESCOPE::SIMULSWVERSION      = "SimulSWVersion";
const char *   EUTELESCOPE::GEOID               = "GeoID";
const char *   EUTELESCOPE::NOOFDETECTOR        = "NoOfDetector";
const char *   EUTELESCOPE::MINX                = "MinX";
const char *   EUTELESCOPE::MAXX                = "MaxX";
const char *   EUTELESCOPE::MINY                = "MinY";
const char *   EUTELESCOPE::MAXY                = "MaxY";
const char *   EUTELESCOPE::INTERMEDIATEFILE    = "IntermediateFile";
const char *   EUTELESCOPE::APPLIEDPROCESSOR    = "AppliedProcessor";
const char *   EUTELESCOPE::DAQDATA             = "DAQData";
const char *   EUTELESCOPE::SIMULDATA           = "SimulData";
const char *   EUTELESCOPE::CONVDATA            = "ConvData";
const char *   EUTELESCOPE::EUDRB               = "EUDRB";
const char *   EUTELESCOPE::IPHCIMAGER          = "IPHCIMAGER";
const char *   EUTELESCOPE::SUCIMAIMAGER        = "SUCIMAIMAGER";
const char *   EUTELESCOPE::EUDAQ               = "EUDAQ";
const int      EUTELESCOPE::GOODPIXEL           =  0;
const int      EUTELESCOPE::BADPIXEL            =  1;
const int      EUTELESCOPE::HITPIXEL            = -1;
const char *   EUTELESCOPE::ABSOLUTENOISEVALUE  = "AbsoluteNoiseValue";
const char *   EUTELESCOPE::NOISEDISTRIBUTION   = "NoiseDistribution";
const char *   EUTELESCOPE::MEANRMS             = "MeanRMS";
const char *   EUTELESCOPE::AIDAPROFILE         = "AIDAProfile";
const char *   EUTELESCOPE::FIXEDFRAME          = "FixedFrame";
const char *   EUTELESCOPE::MATRIXDEFAULTENCODING  = "sensorID:5,xMin:12,xMax:12,yMin:12,yMax:12";
const char *   EUTELESCOPE::CLUSTERDEFAULTENCODING = "sensorID:5,clusterID:8,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5";
const char *   EUTELESCOPE::FIXEDWEIGHT         = "FixedWeight";
