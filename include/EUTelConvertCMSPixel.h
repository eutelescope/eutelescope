// Version: $Id: EUTelConvertCMSPixel.h 2348 2013-02-06 13:54:11Z spanns $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCMSPIXELREADER_H
#define EUTELCMSPIXELREADER_H 1

//Include the Decoder:
#include "CMSPixelDecoder.h"

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

  //! EUTelescope Native Reader processor for the CMS Pixel Detector (PSI46) read-out chips
  /*! This processor can be used to translate native detector data from the PSI46
   * pixel read-out chips (ROCs) into the common LCIO format to process the data with the
   * EUTelescope framework. The reader uses an instance of the CMSPixelDecoder to
   * decode the raw data stream into single events, hits and the respective the pixel
   * address and pulse height from the data stream. Furthermore different checks on the
   * decoded pixel hit or complete event can be performed using a set of decoding flags.
   * @see CMSPixelDecoder
   * 
   * The processor is able to decode data streams from any number of both analog 
   * PSI46v2, PSI46xdb or digital PSI46dig ROCs. The data stream can contain TBM
   * signatures or can be recorded without. In addition a set of parameters and flags
   * can be used to influence the verbosity and decoding behavior:
   * 
   * @param sparseDataCollection specifies the LCIO output collection name
   * @param FileName is the native raw data input file.
   * @param runNumber is the DUT run number to be embedded into the LCIO run header.
   * @param HistogramFilling can be used to turn the creation of histograms on and off.
   * @param addressLevelsFile is the file which contains the address level binning for
   * TBM and ROCs. This parameter is only needed for analog ROCs.
   * @param debugDecoder sets the verbosity level of the raw data decoder class.
   * @param 
   * @param haveTBMheaders switch between data streams recorded with TBM or TBM emulation and
   * streams without any TBM data.
   * @param writeEmptyEvents Enable or disable the writing of events which contain no
   * pixel hits. This is mandatory if using the ROC as DUT in a telescope to maintain
   * event correlation.
   * @param digitalROC lets you choose between the different readout flavors of PSI46 chips,
   * analog PSI46v2 and PSI46xdb or digital PSI46dig.
   * @param useIPBus swtiches the data input format from the native raw file between the
   * readout with the PSI testboard over USB interface and the RAL testboatd with IPBus over
   * optical Ethernet.
   * @param shufflePlanes gives the possibility to re-order the ROCs coming from the native
   * raw data. This is useful if the readout order of the ROCs is not identical with the
   * order of the ROC telescope planes.
   *
   */

    class EUTelConvertCMSPixel: public marlin::DataSourceProcessor 
    {
        public:
            /*! This method returns an new instance of the this processor. It
            *  is called by Marlin execution framework and it shouldn't be
            *  called/used by the final user.
            *
            *  @return a new CMSPixelReader.
            */
            virtual Processor * newProcessor() {
                return new EUTelConvertCMSPixel;
            }

            //! Default constructor
            EUTelConvertCMSPixel();

	    //! Creates LCIO events from the native raw data
	    /*! This is the real method. This is looping over all the events
	     *  containied into the native file and calling the appropriate
	     *  decoder for each of the subevent found into the current event.
	     *
	     *  @param Ntrig The number of events to be read out.
	     */
            virtual void readDataSource (int Ntrig);

	    //! Called at the job beginning.
	    /*! This is executed only once in the whole execution. It prints
	     *  out the processor parameters and reset all needed data
	     *  members.
	     */
            virtual void init ();

	    //! Called after data processing.
	    /*! This method is called when the loop on events is finished. It
	     *  prints only a goodbye message
	     */
            virtual void end ();

	    //! Initialize the geometry information
	    /*! This method is called to initialize the geometry information,
	     *  namely the total number of sensors and the boundaries for each sensors.
	     *
	     *  @throw In case the geometry file does not contain all the needed
	     *  information, a InvalidGeometryException is thrown.
	     */
            void initializeGeometry();            

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	    //! Book histograms
	    /*! This method is used to prepare the needed directory structure
	     *  within the current ITree folder and books all required
	     *  histograms. Histogram pointers are stored into
	     *  EUTelHistogramMaker::_aidaHistoMap so that they can be
	     *  recalled and filled from anywhere in the code.
	     */
            void bookHistos();

	    //! Fill histograms
	    /*! This method is called for each event and the cluster
	     *  information are inserted into specific AIDA histograms.
	     *
	     *  @param xCoord The current pixel x coordinate
	     *  @param xCoord The current pixel y coordinate
	     *  @param xCoord The current pixel pulse height
	     *  @param xCoord The current sensor ID
	     */
            void fillHistos (int xCoord, int yCoord, int value, int sensorID,
			     int64_t timestamp, int triggerphase, int evt);

	    // FIXME probably move to protected?
	    static std::string _triggerPhaseHistoName;
	    static std::string _triggerPhaseHitHistoName;
	    static std::string _triggerPhaseHitCutHistoName;
	    static std::string _dcolMonitorHistoName;
	    static std::string _dcolMonitorEvtHistoName;
	    //! Histogram name of the hit map
	    static std::string _hitMapHistoName;
	    static std::string _hitMapTrigHistoName;
	    static std::string _hitMapCutHistoName;
	    //! Histogram name of the pulse height distribution
            static std::string _pulseHeightHistoName;
#endif            
    
        protected:
	    //! Native rae data file name
            std::string _fileName;
	    
	    //! Address levels file name for analog PSI46 chips
            std::string _levelsFile;

	    //! Output collection name
            std::string _sparseDataCollectionName;

	    //! Number of pixels in x (columns) as read from GEAR file
            unsigned int _noOfXPixel;
	    //! Number of pixels in y (rows) as read from GEAR file
            unsigned int _noOfYPixel;
            
	    //! Run number string
	    std::string _srunNumber;

	    //! Timestamp of the first event recorded
	    int64_t timestamp_event1;

	    //! Run number
            int _runNumber;

	    //! Number of ROCs as read from GEAR file
            unsigned int _noOfROC;

	    //! Writing empty events for event correlation
            bool _writeEmptyEvents;
	    
	    //! The ROC type. Can be: psi46v2, psi46xdb, psi46dig_trig, psi46dig, psi46digv2_b, psi46digv2.
            int _ROC_type;

	    //! Testboard type. Can be: PSI_ATB, PSI_DTB, RAL.
	    std::string _TB_type;

	    //! Reordering of telescope planes w.r.t. readout order
            IntVec _shufflePlanes;

            IntVec _cutHitmap;

	    //! TBM data stream enabled
            bool _haveTBM;

	    //! IPBus data stream enabled
            bool _useIPBus;

            short *_buffer;
            std::vector<int >_excludePlane;
            std::vector<int >_setupPlane;

            //! Decoder verbosity level
	    std::string _debugSwitch;

	    //! Enable histogram writing
            bool _fillHistos;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	    //! Map of histograms
	    std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;
#endif        	

	    //! Geometry correctly read from GEAR file
	    bool _isGeometryReady;

	    //! Planes parameters read from GEAR file
	    gear::SiPlanesParameters* _siPlanesParameters;

	    //! Layer layout read from GEAR file
	    gear::SiPlanesLayerLayout* _siPlanesLayerLayout;  

	    //! Layer index map of the ROC planes
	    std::map< int , int > _layerIndexMap;        	      	
           
        private:
	    //! First procesed event
            bool _isFirstEvent;

	    //! Event number in current run
            int eventNumber;

	    
            unsigned int iROC;
            int status;
            
	    //! Flags to be sent to decoder instance
	    int flags;

            char filename[80];
            char levelsFile[80];
            char exception[80];

    };
    //! A global instance of the processor
    EUTelConvertCMSPixel gEUTelConvertCMSPixel;
    
} // end namespace eutelescope

#endif
