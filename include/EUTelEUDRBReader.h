/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifdef EXPERIMENTAL
#ifndef EUTELEUDRBREADER_H 
#define EUTELEUDRBREADER_H 1

// personal includes ".h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>


namespace eutelescope {

  //! This is the file header. 
  /*! There is a structure like this at the beginning of the file and
   *  contains many information about the current setup.
   *  
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  struct EUDRBFileHeader {
    
    //! The total number of events stored in the file
    /*! This number represents how many events are saved into the
     *  file. This number is not terribly important from the point of
     *  view of this file structure, since the input file is read
     *  within a while loop until the EOF is reached, but it is
     *  important for the following analysis steps where it is
     *  important to know the number of events in the file. For the
     *  time being neither the BORE not the EORE are implemented (see
     *  <a
     *  href=http://forum.linearcollider.org/index.php?t=tree&th=295&rid=239&S=a6eba1d4660c56539adfcabc48017ce9#page_top>the
     *  linear collider forum</a>)
     *
     *  <code>
     *  sizeof(EUDRBFileHeader) = 32
     *  </code>
     */ 
    int  numberOfEvent;  //  4 bytes
    
    //! The number of separate detector 
    /*! For the time being this number is going to be equal to one; in
     *  fact only one detector is read with one EUDRB board. As soon
     *  as multiple detectors will be read, this number can be
     *  different by 1, but for that time I hope we will have already
     *  the final data format (LCIO based) not using anymore this
     *  basic debug format
     */ 
    int  numberOfDetector; // 4 bytes

    //! The number of pixel along x 
    /*! This is the number of pixel along the x direction for each
     *  single channel. So for example this is 66 for a MimoTel
     *  detector
     */ 
    int nXPixel; // 4 bytes

    //! The number of pixel along y
    /*! This is the number of pixel along the y direction for each
     *  single channel. So for example this is 256 for a MimoTel
     *  detector
     */ 
    int nYPixel; // 4 bytes
    
    //! The total event size
    /*! This is the number of bytes contained in one event
     *  structure. This is actually equivalent to: <code>
     *  sizeof(EUDRBEventHeader) + sizeof(EUDRBDataBlock) +
     *  sizeof(EUDRBTrailer)</code>. For the time being this is
     *  equivalent to:
     *
     *  \li @c sizeof(EUDRBEventHeader) 8 bytes
     *  \li @c dataSize
     *  \li @c sizeof(EUDRBTrailer) 4 bytes
     */ 
    int eventSize;
    
    //! The total data size
    /*! This is the size in bytes of the data block. It is equal to
     *  the following: <code> 4 channel * nXPixel * nYPixel * 3 frame
     *  / 2 samples per record </code>
     */
    int dataSize;

    //! Data bit-mask for channel A and C
    /*! This integer number is used to mask the data part for channels
     *  A and C in the transferred bus. This mask has to applied to
     *  each record and then a right shift must be applied to obtained
     *  the ADC value.
     *
     *  For the time being this bit mask is 0x0FFF0000;
     */ 
    int chACBitMask; // 4 bytes

    //! Right shift for the channel A/C data
    /*! To obtain the ADC value for channels A and C from the
     *  transferred buffer record, first the chACBitMask has to be
     *  applied, and then a right shift of @a chACRightShift must be
     *  applied
     *
     *  For the time being this number is 16
     */
    int chACRightShift; // 4 bytes

    //! Data bit-mask for channel B and C
    /*! This integer number is used to mask the data part for channels
     *  B and D in the transferred bus. This mask has to applied to
     *  each record and then a right shift must be applied to obtained
     *  the ADC value.
     *
     *  For the time being this bit mask is 0x00000FFF;
     */ 
    int chBDBitMask; // 4 bytes

    //! Right shift for the channel B/D data
    /*! To obtain the ADC value for channels B and D from the
     *  transferred buffer record, first the chBDBitMask has to be
     *  applied, and then a right shift of @a chBDRightShift must be
     *  applied
     *
     *  For the time being this number is 0
     */
    int chBDRightShift; // 4 bytes

  };

  //! This is the event file header.
  /*! There is a structure like this at the beginning of each
   *  event. The total size is 8 bytes.
   * 
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */ 
  struct EUDRBEventHeader {
    
    //! The current event number (starting from 0)
    int eventNumber;     //  4 bytes
    
    //! The current trigger number if available 
    int triggerNumber;   //  4 bytes
    
  };


  //! This is the EUDRB trailer
  /*! This is the trailer appended at the end of each event.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$   
   */
  struct EUDRBTrailer {
    //! The trailer
    unsigned int trailer;  // 4 bytes
  };
  
  
  
  //!  Reads test data set written with the EUDRB board
  /*!  During the debug phase of the EUDRB in non zero suppressed
   *   mode, data are saved on disk as they are coming from the VME
   *   bus without any modifications. This is required to avoid any
   *   corruption due to data manipulation by the very preliminary DAQ
   *   software.
   *
   *   The result of the MBLT data transfer from the VME bus is a
   *   buffer of integers in which the value of two pixel samples are
   *   saved into a 32 bit integer number with the following structure:
   *
   *   <code>
   *   31 . . 28 27 . . 24 23 . . 20 19 . . 16 15 . . 12 11 . . 8 7 . . 4 3 . . 0<br>
   *   <-  c1 -> <-----    pixel A / C   ----> <-  c2 -> <----   pixel B / D --->
   *   </code>
   * 
   *   Odd records contain pixels from channels A and B, while even
   *   records contain pixels from channels C and D. @c c1 and @c c2
   *   are two control fields not yet considered in the data
   *   processing.
   *
   *   At the beginning of the file there is EUDRBFileHeader that is a
   *   struct contains important detector and setup information.  At
   *   the beginning of each event there is instead a EUDRBEventHeader
   *   that is used only to store the event number; this is followed
   *   by the data block which size corresponds to the following
   *   formula:
   *
   *   <code>
   *   xPixel * yPixel * nChannel * nFrame / pixelPerRecord
   *   </code>
   *
   *   At the end of each event a trailer containing 0x89ABCDEF is
   *   written.
   *
   *   The user can select which algorithm has to be applied to the
   *   incoming data. Those are the available choices: \li
   *   <b>CDS32</b>. The TrackerRawData contains the difference from
   *   frame 3 and frame 2. \li <b>CDS21</b>. As above but the
   *   TrackerRawData contains the difference from frame 2 and frame
   *   1. \li <b>LFn</b> with n=1, 2, 3. The TrackerRawData contains
   *   the frame specified by the integer number.
   *   
   *   <h4>Input - Prerequisites</h4> None
   *
   *   <h4>Output</h4>
   *   The event with the TrackerRawData collection
   *   
   *   @param FileName The input file name in the EUDRB format
   *
   *   @param CalculationAlgorithm The algorithm to be used to fill
   *   the TrackerRawData
   *
   *   @author  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *   @version $Id$
   *
   */
  
  class EUTelEUDRBReader: public marlin::DataSourceProcessor {
    
  public:
    
    //! Default constructor
    EUTelEUDRBReader ();
    
    //! New processor
    /*! Return a new instance of a EUTelEUDRBReader. It is
     *  called by the Marlin execution framework and shouldn't be used
     *  by the final user.
     */
    virtual EUTelEUDRBReader *newProcessor ();
    
    //! Creates events from the debug data stream
    /*! At the beginning of the development of a new hardware device,
     *  it is crucial to test the quality of the transmitted data
     *  stream taking care of reducing as much as possible the number
     *  of manipulations the DAQ software has to perform on the data
     *  stream. This method is taking as input the data stream as it
     *  comes from the VME bus with no manipulation at all and
     *  reconstruct the LCIO data structure needed for basic data
     *  analysis like the pedestal and noise distribution.
     *
     *  @param numEvents This is the total number of events the should
     *  be processed. This value it is actually passed to the
     *  DataSourceProcessor by the ProcessorMgr
     */
    virtual void readDataSource (int numEvents);
    
    //! Init method
    /*! It is called at the beginning of the cycle and it prints out
     *  the parameters.
     */
    virtual void init ();
    
    //! End method
    /*! It deletes the EUTelEUDRBReader::_buffer array
     */
    virtual void end ();
    
  protected:
    
    //! Input file name 
    std::string _fileName;
    
    //! Calculation algorithm
    std::string _algo;
    
  private:
    
    //! A EUDRBFileHeader instance
    /*! This object is used to read the file header from the input
     *  file and the content is used to keep all the useful
     *  information for the data processing
     */ 
    EUDRBFileHeader * _fileHeader;
    
    //! The buffer container
    int * _buffer;
        
  };

  
}                               // end namespace eutelescope
#endif
#endif
