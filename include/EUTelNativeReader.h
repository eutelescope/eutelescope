// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifdef USE_EUDAQ 
#ifndef EUTELNATIVEREADER_H
#define EUTELNATIVEREADER_H 1

// personal includes ".h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// eudaq includes <.h>

// lcio includes <.h>

// system includes <>
#include <string>

namespace eutelescope {

  //!  Reads the data produced by the EUDRB boards
  /*!  This Marlin reader is taking as an input the output file of
   *   the eudaq software and converting it to a LCIO event. This is
   *   linking against libeudaq and using directly the data structure
   *   used in the DAQ software. This has to be thought as a sort of
   *   link between the native DAQ raw format and the LCIO data model
   *   used for the telescope data description.
   *
   *   This processor is automatically guessing both the kind of data
   *   contained into the raw file, i.e. RAW2, RAW3, ZS and also the
   *   detector type (MimoTel and Mimosa18).
   *
   *   The main goal of this data reader is to test the quality of
   *   the data output from the hardware but it is going to disappear
   *   soon being the LCIO output already produced by the online
   *   system.
   *
   *   @see @ref cds3frame
   *
   *   The marker removal procedure is done a very general way
   *   allowing to copy paste the same code for a detector with a
   *   different configuration of markers.
   *
   *   <h4>Input collection</h4>
   *   None
   *
   *   <h4>Output collections</h4>
   * 
   *   <b>First frame</b> A collection of TrackerRawData containing
   *   the first decoded frame
   *
   *   <b>Second frame</b> A collection of TrackerRawData containing
   *   the second decoded frame
   *
   *   <b>Third frame</b> A collection of TrackerRawData containing
   *   the third decoded frame
   *
   *   <b>CDS</b> A collection of TrackerRawData containing
   *   the CDS make taking into account the pivot pixel.
   *
   *   @param FirstFrameCollectionName The name of the first frame
   *   collection.
   *
   *   @param SecondFrameCollectionName The name of the second frame collection.
   *
   *   @param ThirdFrameCollectionName The name of the third frame
   *   collection.
   *
   *   @param CDSCollection The name of the CDS collection.
   *
   *   @param CDS Switch to enable / disable the CDS calculation.
   *   
   *   @param InputDataFileName The input file to convert.
   *
   *   @param SignalPolarity The expected signal polarity.
   *
   *   @param GeoID This is the identification number of the used
   *   geometry. It has to be the same here and in the corresponding
   *   XML geometry description.
   *
   *   @param RemoveMarker This boolean is used to remove (true) the
   *   markers from the output data stream. 
   *
   *   @param MarkerPosition This is a vector of integer containing
   *   the marker position in pixel number start counting from 0.
   *
   *   @author  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *   @version $Id: EUTelNativeReader.h,v 1.1 2008-08-06 20:37:00 bulgheroni Exp $
   *
   */
  
  class EUTelNativeReader : public marlin::DataSourceProcessor    {
    
  public:
     
    //! Default constructor
    EUTelNativeReader ();
     
    //! New processor
    /*! Return a new instance of a EUTelNativeReader. It is
     *  called by the Marlin execution framework and shouldn't be used
     *  by the final user.
     */
    virtual EUTelNativeReader * newProcessor ();
     
    //! Creates events from the eudaq software
    /*! This method reads a certain number of events from the input
     *  raw data file and generates LCIO events with at least three
     *  collections. This processor is very specific for the MimoTel
     *  setup, so it cannot be used in general to read the output of
     *  the EUDRB board. This limitation is due to the fact that the
     *  main goal of this data reader is to check the quality of the
     *  data producer and soon the LCIO output will be provided
     *  directly from the DAQ software.
     *
     *  @param numEvents The number of events to be read out.
     */
    virtual void readDataSource (int numEvents);

    //! Init method
    /*! It is called at the beginning of the cycle and it prints out
     *  the parameters.
     */
    virtual void init ();

    //! End method
    /*! It prints out a good bye message 
     */
    virtual void end ();

  protected:

    //! The input file name
    /*! It is set as a Marlin parameter in the constructor
     */ 
    std::string _fileName;

    //! The first frame collection name
    std::string _firstFrameCollectionName;

    //! The second frame collection name
    std::string _secondFrameCollectionName;

    //! The third frame collection name
    std::string _thirdFrameCollectionName;

    //! The CDS collection name
    std::string _cdsCollectionName;

    //! The zero suppressed frame collection name
    std::string _zsFrameCollectionName;

    //! The CDS enable flag
    /*! This flag is true if the converter should perform also online
     *  CDS calculation. It is false otherwise.
     */
    bool _cdsCalculation;

    //! Signal polarity 
    /*! This is used to change the signal polarity in the CDS
     *  calculation. 
     */
    int _polarity;

    //! Geometry ID
    /*! This is the unique identification number of the telescope
     *  geometry. This identification number is saved in the run
     *  header and then crosscheck against the XML geometry
     *  description during the reconstruction phase. 
     *  
     *  In the future, this ID can be used to browse a geometry database.
     */ 
    int _geoID;

    //! Marker removal switch
    /*! This boolean is used to set the removal of markers from the
     *  TrackerRawData output collections.
     *  When this switch is set to true, the output collection will
     *  not contain the marker information.
     */
    bool _removeMarkerSwitch;

    //! Marker position vector
    /*! This vector of integer contains the position of all the
     *  markers that have to stripped away from the output
     *  collection. 
     *  Markers are assumed to occur always at the same position in
     *  the row.
     *
     */
    std::vector<int > _markerPositionVec;

    //! Type of sparsified pixel
    /*! Which information of the pixel passing the zero suppression
     *  can be stored into different data structure. The user can
     *  select which one in the steering file using the
     *  SparsePixelType enumerator
     */
    int _pixelType;

    //! Maximum number of pixel along the x direction
    /*! This variable can assume two values depending if the user
     *  wants to remove the markers of not
     */ 
    int _xMax;

    //! Maximum number of pixel along the y direction
    /*! Since on the MimoTel there are only columns of markers (and no
     *  row), the maximum number of pixel along y is constant.
     */ 
    int _yMax;

    //! Activate / Deactivate out of synch skipping
    /*! In normal conditions all the boards should work perfectly
     *  synchronized, but due to several reasons the synchronization
     *  can be lost and the corresponding data might be
     *  corrupted. Activating this switch, it is possible to remove
     *  from the output file all events not properly synchronized.
     * 
     */
    bool _skipOutOfSynch;

    //! Out of synch threshold
    /*! The definition of an out of synch event is based upon the
     *  value of pivot pixel address recorded by all the boards. An
     *  event is declared de-synchronized if the difference between the
     *  maximum and the minimum pivot pixel address is exceed this threshold.
     *
     */ 
    int _outOfSynchThr;
  };

  //! A global instance of the processor
  EUTelNativeReader gEUTelNativeReader;

}                               // end namespace eutelescope
#endif
#endif
