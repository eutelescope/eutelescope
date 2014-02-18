// Version: $Id$
/*========================================================================*/
/*          CMSPixel CalibrateEvent (Calibration of raw detector data)    */
/*          Author: Simon Spannagel                                       */
/*                (simon.spannagel@student.kit.edu or s.spannagel@cern.ch)*/
/*          Created       23 feb 2012                                     */
/*          Last modified 04 apr 2012                                     */
/*========================================================================*/

/*
*   This source code is part of the Eutelescope package of Marlin.
*   You are free to use this source files for your own development as
*   long as it stays in a public research context. You are not
*   allowed to use it for commercial purpose. You must put this
*   header with author names in all development based on this file.
*
*/
#ifndef CMSPIXELCALIBRATEEVENT_H
#define CMSPIXELCALIBRATEEVENT_H 1

// Marlin includes
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// AIDA includes
#include <AIDA/IBaseHistogram.h>
#endif

// GEAR includes
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// system includes <>
#include <string>
#include <map>
#include <vector>

#define MAXROC 16
#define MAXEVENTS 4160

namespace eutelescope {

    typedef struct {
        double par0;
        double par1;
        double par2;
        double par3;
    } cal_param;


    class CMSPixelCalibrateEventProcessor:public marlin::Processor {

        public:

            virtual Processor * newProcessor() {
                return new CMSPixelCalibrateEventProcessor;
            }

            //! Default constructor
            CMSPixelCalibrateEventProcessor ();

            virtual void init ();
            virtual void processRunHeader (LCRunHeader * run);
            virtual void processEvent (LCEvent * evt);
            virtual void check (LCEvent * evt);
            virtual void end();

            void initializeGeometry();            
            void initializeCalibration() throw ( marlin::StopProcessingException );
            void initializeGaintanhCalibration() throw ( marlin::StopProcessingException );
            
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            void fillHistos (int raw, int calibrated, int sensorID);
            void bookHistos();
        	static std::string _pulseVcalHistoName;
        	static std::string _pulseHeightHistoName;
        	static std::string _NANMapHistoName;
#endif            

        protected:

            std::string _sparseDataCollectionName;
            std::string _calibratedDataCollectionName;
            std::string _calibrationFile;
            int _pixelType;
            int _iRun;
            bool _fillHistos;
            bool _phCalibration;
            
            unsigned int _noOfXPixel;
            unsigned int _noOfYPixel;
            unsigned int _noOfROC;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        	std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;
#endif        	
        	bool _isGeometryReady;
        	gear::SiPlanesParameters* _siPlanesParameters;
        	gear::SiPlanesLayerLayout* _siPlanesLayerLayout;  
        	std::map< int , int > _layerIndexMap;        	      	
             
        private:

            bool calTanH(double &corr, double y, double p0, double p1, double p2, double p3);
	    bool calWeibull(double &corr, double y);
	    bool calLinear(double &corr, double y);
            std::vector< std::vector< cal_param > > calibration;

    };

    //! A global instance of the processor
    CMSPixelCalibrateEventProcessor gCMSPixelCalibrateEventProcessor;

}
#endif
