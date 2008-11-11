// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMAPSDIGI_H
#define EUTELMAPSDIGI_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif


// system includes <>
#include <string>
#include <vector>
#include <map>


typedef std::map<unsigned int, double> type_pixelChargeMap;

namespace eutelescope {

  //! MAPS digitization processor
  /*! Realistic model based on the beam test results is used to
   *  calculate the expected response of MAPS sensors for the Geant4
   *  simulated hits, as returned by Mokka.
   *
   *  First step: conversion from global reference frame to local
   *  reference frame of individual sensor (based on the GEAR geometry
   *  description) was adopted from EUTelHitMaker.
   *
   *  Full documentation will be added with the first public version
   *  of the code.
   *
   *  @author Aleksander Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
   *  @version $Id: EUTelMAPSdigi.h,v 1.1 2008-11-11 09:18:11 bulgheroni Exp $
   *
   */

  class EUTelMAPSdigi : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelMAPSdigi
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelMAPSdigi.
     */
    virtual Processor * newProcessor() {
      return new EUTelMAPSdigi;
    }

    //! Default constructor
    EUTelMAPSdigi ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file.
     *
     *  Each simulated hit is translated from the global to the local
     *  frame of reference  thanks to the GEAR geometry description.
     *
     *  The track segment is then divided into smaller fragments and
     *  the expected charge sharing between pixels is calculated based
     *  on the realistic parametrization.
     *
     *  Finally, all pixels fired are put into the output collection
     *  in the format corresponding to sparcified raw data.
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();


    //! Histogram booking
    /*! Some control histograms are filled during this procedure in
     *  order to be able to perform easy check on the quality of the
     *  output hits and also to understand if the frame of reference
     *  conversion has been properly done. Of course this method is
     *  effectively doing something only in the case MARLIN_USE_AIDA.
     */
    void bookHistos();


  protected:

    //
    // Processor parameter definitions start here
    //

    //! SimTrackerHit collection name
    /*! This is the name of the collection holding the simulated hits
     *  from Mokka.
     */
    std::string _simhitCollectionName;

    //! Maxmimu energy per step [MeV]
    /*! Maximum energy for single charge sharing calculation.
     *  If more energy is deposited in single Mokka step, the step
     *  is divided into smaller ones.
     */
    double _stepdEmax;

    //! Maximum step length
    /*! Longer Mokka hits have to be divided into smaller steps
     */
    double  _stepLmax;

    //! Maximum diffusion range in pixels
    /*! Charge sharing is calculated for all pixels closer
     * in X and Y to the seed pixel than the maximum given
     */

    int _spreadRange;

    //! Attenuation length
    /*! Charge attenuation in diffusion
     */

    double _attenLength;

    //! Single Mokka hit output mode
    /*! Flag controling if pixel list should be written to output
     *  file after each Mokka hit is processed or only after
     *  the whole event
     */

    bool _singleHitOutput;

    //! Debug Event Count
    /*! Print out debug and information
     *  messages only for one out of given number of events. If
     *  zero, no debug information is printed.
     */
    int _debugCount ;


  private:

    // Local functions used in charge sharing calculations

    //! Checks track segment start and end point
    /*! Track segment start and end points are compared with sensor
     *  dimensions. Variables _pathStart and _pathEnd are set to
     *  values providing that the trakck is inside the sensitive
     *  volume.
     */

    double CheckPathLimits();


    //! Get energy deposited in one step
    /*! Divides energy deposited in one Mokka hit into nStep deposits
     *   in the sensor. Fluctuations can be taken into account
     *   (not implemented yet).
     */
    double GetStepEnergy(int nStep);


    //! Main routine calculating charge distribution

    double CalculateChargeDistribtion(double stepCenter, double stepDeposit);


    //! Storing calculated charge in pixel map

    void AddPixelCharge(int ix, int iy, double dQ);


    //! Number of stored (fired) pixels
    /* Returns number of pixels stored in pixel charge map
     */

    inline unsigned int GetPixelNumber()
      {
        return _pixelChargeMap->size();
      }


    //! Pixel charge
    /*! Returns charge attributed to given pixel.
     *   Pixels are numbered with position in the pixel map, default
     *   iterator is used.
     */

    double GetPixelCharge(int iPixel);


    //! Pixel index
    /*! Returns pixel index, coding its position in X-Y
     *   Pixels are numbered with position in the pixel map, default
     *   iterator is used.
     */

    int GetPixelIndex(int iPixel);


    //! Pixel position in X

    int GetPixelX(int iPixel);


    //! Pixel position in Y

    int GetPixelY(int iPixel);


    //! Code pixel index

    inline int PixelIndex(int ix, int iy)
      {
        return (ix<<16) + iy;
      }


    //
    // Local variables
    //


    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

    //! Conversion ID map.
    /*! In the data file, each cluster is tagged with a detector ID
     *  identify the sensor it belongs to. In the geometry
     *  description, there are along with the sensors also "passive"
     *  layers and other stuff. Those are identify by a layerindex. So
     *  we need a conversion table to go from the detectorID to the
     *  layerindex.
     */
    std::map< int, int > _conversionIdMap;

    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    //! Name of the local hit map histo
    /*! The histogram pointed by this name is a 2D histo. The x and y
     *  axes correspond to the pixel detector axes in its own local
     *  frame of reference already converted in millimeter. There is
     *  a histo like this for each detector in the geometry.
     */
    static std::string _hitHistoLocalName;

    //! Name of the hit map histo
    /*! The histogram pointed by this name is a 2D histo. The x and y
     *  axes correspond to the pixel detector axes already in the
     *  telescope frame of reference. There is a histo like this for
     *  each detector in the geometry.
     */
    static std::string _hitHistoTelescopeName;

    //! Name of the pixel map histo
    /*! The histogram pointed by this name is a 2D histo. The x and y
     *  axes correspond to the pixel number in its own local frame. There is
     *  a histo like this for each detector in the geometry.
     */
    static std::string _pixelHistoName;

    //! Name of the density plot
    /*! This is a very nice plot showing in a 3D frame where all hits
     *  have been found. If the run is sufficiently long and the hits
     *  are uniformly distributed on the sensor surface, this plot
     *  should recall the shape of the telescope itself.
     */

#endif

    //! Fill histogram switch
    /*! This boolean switch was initially introduced for debug reason
     *  but then we realized that it could stay there and protect
     *  against missing AIDA::Processor.
     *
     */
    bool _histogramSwitch;

    /*!
     *  Following variables contain information about Mokka hit
     *  transformed to the local (sensor) coordinate frame
     */

    //! Mokka hit position in local frame

    double _localPosition[3];

    //! Particle momentum (as stored by Mokka) in local frame

    double _localMomentum[3];

    //! Particle direction in local frame (normalized to 1)

    double _localDirection[3];

    //! Path length corresponding to Mokka hit

    double _mokkaPath;

    //! Start position on the path
    /*! Due to multiple scattering Mokka path can be longer than
     *  the detector thickness. Start position tells us at which track
     *  point the track segment enters the sensor. Should be 0 for
     *  most hits. This is calculated in CheckPathLimits routine
     */

    double _pathStart;

    //! Start position on the path

    double _pathEnd;

    //! Energy deposit corresponding to Mokka hit

    double _mokkaDeposit;

    /*!
     *  Following variables contain information about sensor
     *  dimensions in the local coordinate frame
     */

    //! Sensor size in local reference frame
    /*! Sensor dimensions can be different in local reference frame
     *  as X can swap with Y due to rotation.
     */

    double _localSize[3];

    //! Sensor pitch in XY in local reference frame
    /*! Sensor pitch can be different in local reference frame
     *  X can swap with Y due to rotation
     */

    double _localPitch[2];

    //! Pixel Charge map
    /*! Pixel charges calculated in subsequent steps are added to this
     *  map. Each (fired) pixel is identified by an index coding X-Y
     *  pixel position.
     */

    type_pixelChargeMap *_pixelChargeMap;

    type_pixelChargeMap::iterator _pixelIterator;

    std::map< int,  type_pixelChargeMap *>  _pixelChargeMapCollection;

  };

  //! A global instance of the processor
  EUTelMAPSdigi gEUTelMAPSdigi;

}
#endif
#endif
