/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELSTRASMIMOTELREADER_H
#define EUTELSTRASMIMOTELREADER_H

// personal includes ".h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>
#include <string>
#include <iostream>

#define RUN_HEADER_RESERVED_W32    1040
#define EVENT_HEADER_RESERVED_W32    19


namespace eutelescope {
  
  //!  Strasbourg DAQ Mimotel converter
  /*!  This data reader is used to read the output file produced by
   *   the Strasbourg DAQ system and put them into the LCIO format
   *   according to the EUTelescope data model. 
   *   
   *   In principle this reader can be used for each detector read by
   *   this DAQ, but for simplicity for the time being some values are
   *   hardcoded to properly decode MimoTel information.
   *
   *   <h4>Input collection</h4>
   *   None
   *
   *   <h4>Output</h4>
   *   LCEvent with TrackerRawData collection
   *
   *   @param LEPSIRunNumber Integer number corresponding to the run
   *   number
   *
   *   @param XNoOfPixel Number of pixels along the x axis
   *
   *   @param YNoOfPixel Number of pixels along the x axis
   *  
   *   @param Frame0CollectionName Name of the Frame0 collection
   * 
   *   @param Frame1CollectionName Name of the Frame1 collection
   *
   *   @param CDSCollectionName Name of the CDS collection
   *
   *   @author  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *   @version $Id$
   *
   */

  class EUTelStrasMimoTelReader : public marlin::DataSourceProcessor   {

  public:

    //! Default constructor
    EUTelStrasMimoTelReader ();

    //! New processor
    /*! Return a new instance of a EUTelStrasMimoTelReader. It is
     *  called by the Marlin execution framework and shouldn't be used
     *  by the final user.
     */
    virtual EUTelStrasMimoTelReader *newProcessor ();

    //! Creates events from a SUCIMA Imager ASCII data file
    /*! This methods reads the input ASCII file, produced by the
     *  SUCIMA Imager DAQ and converts it into a LCIO file with a run
     *  header and the event structure.
     */
    virtual void readDataSource (int numEvents);

    //! Init method
    /*! It is called at the beginning of the cycle and it prints out
     *  the parameters.
     */
    virtual void init ();

    //! End method
    /*! It deletes the EUTelStrasMimoTelReader::_buffer array
     */
    virtual void end ();

    //! Add an End Of Run Event
    /*! When the end of the input files is reached,
     *  then add a end of file event to the output LCIO file
     */
    void addEORE();

    //! Class to store run header information
    /*! This is an attempt to write a class mimicking the TRunHeader
     *  struct defined by Gilles in the DAQ code. For the time being
     *  it is just having a reduced set of information.
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @author Gilles CLAUS, LEPSI <mailto:claus@lepsi.in2p3.fr>
     *  @version $Id$
     */ 
    class StrasRunHeader {

    public:
      //! LittleEndian
      /*! Set to 1 if the DAQ system save data in little endian
       */
      int LittleEndian;
      
      //! ByteSex
      /*! The DAQ system always write 0x3615ABCD or 0 in this field
       */ 
      int ByteSex;
      
      
      //! FileNo
      /*! Don't care about this field
       */
      int FileNo;
      
      //! FileEvNb
      /*! Event number per file
       */
      int FileEvNb;
      
      //! TotEvNb
      /*! Total number of event in the run. It could be less if the
       *  run is stopped before the end.
       */ 
      int TotEvNb;

      //! TotalSz
      /*! Event total size. The sum of the header, the data and the
       *  trailer sizes. This number is in <b>bytes</b>.
       */ 
      int TotalSz;

      //! DataSz
      /*! Data size = DBD_VFAS_RAM_SZ_W32 times ADC number (@see
       *  VFasRunNb). This number is in <b>bytes</b> and not in words.
       */
      int DataSz;

      //! RunNo
      /*! Copy of the run number
       */
      int RunNo;

      //! RunType
      /*! Run type = Xrc Run Type parameter (???)
       */
      int RunType;
      
      //! VFasPresentNb
      /*! Number of ADC boards installed in the system
       */
      int VFasPresentNb;
      
      //! VFasRunNb
      /*! Number of ADC used in the current run
       */
      int VFasRunNb;

      //! VFasChanRunNb
      /*! Number of ADC input channels used in the current run
       */
      int VFasChanRunNb;
      
      //! VFasCalcChanNb
      /*! Used for on-line monitoring: should be 0
       */
      int VFasCalcChanNb;

      //! Other fields
      /*! There are here below other field I cannot decode because
       *  their size is defined by preprocessor macro I don't have. Up
       *  to here 52 bytes have been used and since the information
       *  file is in total 4212 bytes, I'm here below defining a
       *  reserved region of the right size
       * 
       */
      int Reserved[RUN_HEADER_RESERVED_W32];

      // there are actually many other fields here below but with a
      // size defined into a preprocessor macro that I don't have.
      // sorry...
      
    }  _runHeader;

    //! Streamer for the StrasRunHeader
    inline friend std::ostream& operator<<(std::ostream& os, const StrasRunHeader& runHeader) { 
      os << "===============================================\n" 
	 << "| LittleEndian   " << runHeader.LittleEndian << "\n"
	 << "| ByteSex        " << std::hex << "0x"<< runHeader.ByteSex  << std::dec  << "\n"
	 << "| FileNo         " << runHeader.FileNo       << "\n"
	 << "| FileEvNb       " << runHeader.FileEvNb     << "\n"
	 << "| TotEvNb        " << runHeader.TotEvNb      << "\n"
	 << "| TotalSz        " << runHeader.TotalSz        << "\n"
	 << "| DataSz         " << runHeader.DataSz       << "\n"
	 << "| RunNo          " << runHeader.RunNo        << "\n"
	 << "| RunType        " << runHeader.RunType      << "\n"
	 << "| VFasPresentNb  " << runHeader.VFasPresentNb << "\n"
	 << "| VFasRunNb      " << runHeader.VFasRunNb << "\n"
	 << "| VFasChanRunNb  " << runHeader.VFasChanRunNb << "\n"
	 << "| VFasCalcChanNb " << runHeader.VFasCalcChanNb << "\n"
	 << "===============================================";
      return os;
    }

    //! Class to store the event header information
    /*! This is again an attempt to interpret the documentation
     *  provided by Gilles for the event header in the LEPSI DAQ
     *  system output format.
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @author Gilles CLAUS, LEPSI <mailto:claus@lepsi.in2p3.fr>
     *  @version $Id$
     */ 
    class StrasEventHeader {
      
    public:
      //! Event trigger
      /*! Increment from 0 at each trigger
       */
      unsigned int EvTrig;
      
      //! Event number
      /*! Increment from 0 at each event taken
       */
      unsigned int EvNo;

      //! Event position
      /*! Event position in file: increment from 0 after each event
       *  write
       */
      unsigned int EvPos;

      //! Event Tag
      /*! Event tag: for future use, should be 0
       */
      int EvTag;

      //! Event date
      /*! Event date = 0
       */
      int EvDate;
      
      //! Event time
      /*! Event time = SSSS ssss (S = second s = ms - 10 ms
       *  resolution) 
       */
      int EvTime;

      //! TrigCnt
      /*! Total triggers received by acquisition system: accepted or
       *  not. 
       */ 
      unsigned int TrigCnt;
      
      //! EvVmeTime
      /*! Time to read VME boards in ms [10 ms res]
       */ 
      unsigned int EvVmeTime;
      
      //! VFasCnt
      /*! Hardware ADC strip counter
       * 
       *  AB: this is supposed to be an array of a fixed size
       *  NB_MAX_VFAS, that is preprocessor macro I don't know. Since,
       *  according to my knowledge the only usable information is the
       *  one contained into VFasCnt[0] so I can forget about the rest
       *  of the size. But this is not allowing me to continue a
       *  proper decoding.
       *
       *  According to something I have found looking around the
       *  NB_MAX_VFAS should be defined as 4. Since for the time being
       *  the only information really needed of this register is the
       *  first, I don't want to touch 
       */
      int VFasCnt;
      
      //! Reserved space
      /*! The total file header size is 112 bytes. Up to here only 36
       *  bytes have been used so I'm defining a reserved space 112 -
       *  36 bytes long
       */ 
      int Reserved[EVENT_HEADER_RESERVED_W32];

    }  _eventHeader;


    //! Streamer operator for the event reader
    friend std::ostream& operator<<(std::ostream& os, const StrasEventHeader& eventHeader) {
      os << "===============================================\n" 
	 << "| EvTrig         " << eventHeader.EvTrig << "\n"
	 << "| EvNo           " << eventHeader.EvNo << "\n"
	 << "| EvPos          " << eventHeader.EvPos << "\n"
	 << "| EvTag          " << eventHeader.EvTag << "\n"
	 << "| EvDate         " << eventHeader.EvDate << "\n"
	 << "| EvTime         " << eventHeader.EvTime << "\n"
	 << "| TrigCnt        " << eventHeader.TrigCnt << "\n"
	 << "| EvVmeTime      " << eventHeader.EvVmeTime << "\n"
	 << "| VFasCnt        " << eventHeader.VFasCnt << "\n"
	 << "===============================================";    
      return os;
    }


    //! Class to store the trailer information
    /*! This is again an attempt to interpret the documentation
     *  provided by Gilles for the trailer in the LEPSI DAQ
     *  system output format.
     * 
     *  Gilles in the documentation file "das_data_format_d7.txt"
     *  claims it should be 0xDEADFACE, but it looks like to be 0x89ABCDEF.
     *
     *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
     *  @author Gilles CLAUS, LEPSI <mailto:claus@lepsi.in2p3.fr>
     *  @version $Id$
     */ 
    class StrasEventTrailer {
      
      //! End of record. 
      /*! Gilles in the documentation file "das_data_format_d7.txt"
       *  claims it should be 0xDEADFACE, but it looks like to be 0x89ABCDEF. 
       */
    public:
      
      int Eor;
    } _eventTrailer;

    //! Streamer for the Strasbourg even trailer
    friend std::ostream& operator<<(std::ostream& os, const StrasEventTrailer& trailer) {
      os << std::hex << "0x" << trailer.Eor << std::dec;
      return os;
    }

  protected:

    //! Data block
    /*! This is used to read the data block in between the event
     *  header and the trailer. It has a size that is variable and
     *  defined once forever in the StrasRunHeader::DataSz.
     *
     */ 
    int * _dataBuffer;

    //! Input run number
    /*! This is the run number and it is read directly from the
     *  steering file. This number is used as a replacement of the
     *  name, since the files have are all contained into a directory
     *  named after this number. Marlin has to be executed from the
     *  data folder and the directory named after the run should be
     *  there.
     *
     */ 
    int _runNumber;

    //! Run header file name
    /*! Also known as info file, this is in the run folder and is
     *  named according to the following rule:
     *  @c RUN_XXX_i.rz
     *  where XXX is the run number.
     */ 
    std::string _runHeaderFileName;
    
    //! Data file base name
    /*! The DAQ output data are split in several files having all in
     *  principle the same size but the last one that can have a
     *  smaller size if not properly closed the run. Files are named
     *  according to the following rule:
     *  @c RUN_XXX_Y.rz
     *  where XXX is the run number and Y is the file number.
     */ 
    static std::string _dataFileBaseName;

    //! File name extension
    static const std::string _fileNameExt;
    
    //! Number of pixels along X
    /*! This parameter can be selected by the user via the steering
     *  file
     */ 
    int _noOfXPixel;

    //! Number of pixels along Y
    /*! This parameter can be selected by the user via the steering
     *  file
     */ 
    int _noOfYPixel;

    //! Name of the frame 0 collection
    /*! This is the name of the first frame collection. 
     */
    std::string _frame0CollectionName;
    
    //! Name of the frame 1 collection
    /*! This is the name of the second frame collection. 
     */
    std::string _frame1CollectionName;   

    //! Name of the frame 1 collection
    /*! This is the name of the second frame collection. 
     */
    std::string _cdsCollectionName; 

    //! Static variable == 1 for the time being
    static int _noOfSubMatrix;
    
  };

  //! A global instance of the processor
  EUTelStrasMimoTelReader gEUTelStrasMimoTelReader;

  

}                               // end namespace eutelescope
#endif
