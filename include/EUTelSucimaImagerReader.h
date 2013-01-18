/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELSUCIMAIMAGERREADER_H
#define EUTELSUCIMAIMAGERREADER_H

// personal includes ".h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>


namespace eutelescope
{

   //!  Reads SUCIMA Imager ASCII file
   /*!  This Marlin processor is used to convert existing data file
    *   acquired using the SUCIMA Imager DAQ system. The output data
    *   format is a very inconvenient one in terms of data storage and
    *   I/O performance, but, at the same time, it is very practical
    *   from the debug point of view. Pixel signals (short integers)
    *   are saved in the ASCII file in a single row and are tab
    *   separated.
    *   
    *   <h4>Input - Prerequisites</h4>
    *   SUCIMA Imager ASCII data file
    *   Number of pixels in the horizontal direction
    *   Number of pixels in the vertical direction
    *
    *   <h4>Output</h4>
    *   LCEvent with TrackerRawData collection
    *
    *   @param   SUCIMAImagerFileName  name of the input file
    *   @param   NoOfXPixel number of pixels along X
    *   @param   NoOfYPixel number of pixels along Y
    *   @author  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
    *   @version $Id$
    *
    */

   class EUTelSucimaImagerReader:public marlin::DataSourceProcessor
   {

    public:

      //! Default constructor
      EUTelSucimaImagerReader ();

      //! New processor
      /*! Return a new instance of a EUTelSucimaImagerReader. It is
       *  called by the Marlin execution framework and shouldn't be used
       *  by the final user.
       */
      virtual EUTelSucimaImagerReader *newProcessor ();

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
      /*! It deletes the EUTelSucimaImagerReader::_buffer array
       */
      virtual void end ();

    protected:

      //! Input file name 
        std::string _fileName;

      //! Number of pixels along X
      int _noOfXPixel;

      //! Number of pixels along Y
      int _noOfYPixel;

      //! The buffer to store data from file
      /*! This array of short is used to temporary store data from the
       *  disk before moving them to the TrackerRawData. This is a good
       *  attitude because, if something goes wrong with the data
       *  reading (most likely a I/O error, you still have the chance of
       *  save the other data. This array is dynamically allocated in
       *  the EUTelSucimaImagerReader::readDataSource(int) only when
       *  processing the first event. It is deleted in the
       *  EUTelSucimaImagerReader::end()
       */
      short *_buffer;

   };

  //! A global instance of the processor
  EUTelSucimaImagerReader gEUTelSucimaImagerReader;

}                               // end namespace eutelescope
#endif
