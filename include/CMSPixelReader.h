/*========================================================================*/
/*          CMSPixel file converter (RAW->LCIO)                           */
/*          Author: Simon Spannagel (s.spannagel@cern.ch)                 */
/*          Created       23 feb 2012                                     */
/*          Last modified 19 jul 2012                                     */
/*========================================================================*/

#ifndef EUTELCMSPIXELREADER_H
#define EUTELCMSPIXELREADER_H

// EUTelescope includes
#include "EUTELESCOPE.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// AIDA includes
#include <AIDA/IBaseHistogram.h>
#endif

// Marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"

// GEAR includes
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// System includes
#include <vector>
#include <map>
#include <string.h>

namespace eutelescope
{
    class CMSPixelReader: public marlin::DataSourceProcessor 
    {
        public:
            /*! This method returns an new instance of the this processor. It
            *  is called by Marlin execution framework and it shouldn't be
            *  called/used by the final user.
            *
            *  @return a new CMSPixelReader.
            */
            virtual Processor * newProcessor() {
                return new CMSPixelReader;
            }

            //! Default constructor
            CMSPixelReader();
            virtual void readDataSource (int Ntrig);
            virtual void init ();
            virtual void end ();

            void initializeGeometry();            

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            void fillHistos (int xCoord, int yCoord, int value, int sensorID);
            void bookHistos();
        	static std::string _hitMapHistoName;
        	static std::string _pulseHeightHistoName;
#endif            
    
        protected:
     
            std::string _fileName;
            std::string _levelsFile;
            std::string _sparseDataCollectionName;

            unsigned int _noOfXPixel;
            unsigned int _noOfYPixel;
            std::string _srunNumber;
            int _runNumber;
            unsigned int _noOfROC;
            bool _writeEmptyEvents;
            bool _digitalROC;
            IntVec _shufflePlanes;
            bool _lazyDecoding;
            int _event_selection;
            bool _haveTBM;
            bool _useIPBus;
            short *_buffer;
            std::vector<int >_excludePlane;
            std::vector<int >_setupPlane;
            //! DEBUG
            int _debugSwitch;

            bool _fillHistos;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        	std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;
#endif        	
        	bool _isGeometryReady;
        	gear::SiPlanesParameters* _siPlanesParameters;
        	gear::SiPlanesLayerLayout* _siPlanesLayerLayout;  
        	std::map< int , int > _layerIndexMap;        	      	
           
        private:
            bool _isFirstEvent;
            int eventNumber;
            unsigned int iROC;
            int status;
            int flags;
            char filename[80];
            char levelsFile[80];
            char exception[80];

    };
    CMSPixelReader gCMSPixelReader;
    
} // end namespace eutelescope

#endif
