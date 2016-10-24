/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if USE_GEAR
// #if defined(USE_GEAR)
#ifndef EUTELAPPLYALIGNMENTPROCESSOR_H
#define EUTELAPPLYALIGNMENTPROCESSOR_H 1

// eutelescope includes ".h"
#include "EUTelAlignmentConstant.h"
#include "EUTelEventImpl.h"
#include "EUTelReferenceHit.h"
#include "EUTelExceptions.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
//#include "marlin/EventModifier.h"


// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// ROOT includes
#include "TString.h"

// lcio includes
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
//#include <TrackerHitImpl2.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <IMPL/TrackImpl.h>

// system includes <>
#include <iostream>
#include <string>
#include <map>
#include <vector>


namespace eutelescope {

  //! Apply the alignment correction to the input hit collection
  /*! The main goal of this processor is to apply the alignment
   *  correction to the input hit collection.
   *
   *  The alignment constants are, as well, contained into collections
   *  of EUTelAlignmentConstant objects and can be loaded from an
   *  external LCIO file or from a DB query.
   *
   *  <h4>Input collections</h4>
   *  <br><b>InputHit</b>.
   *  This is a collection of TrackerHit containing all the hits that
   *  need to be corrected for misalignment.
   *
   *  <br><b>AlignmentConstant</b>.
   *  This is a collection of EUTelAlignmentConstant. This is a
   *  LCGenericObject  containing all the needed alignment constants
   *  calculated by previous processors.
   *
   *  <h4>Output</h4>
   *
   *  <br><b>OutputHit</b>.
   *  This is a collection of TrackerHit with the correct hits.
   *
   *  @param CorrectionMethod There are actually several different
   *  methods to apply the alignment constants. Here below a list of
   *  available methods:
   *     0. <b>Shifts only</b> No rotations angle will be taken into
   *  account.
   *     1. <b>Rotation first</b> The rotational matrix is applied
   *  before the translation.
   *     2. <b>Translation first</b> The hits are first shifted and
   *  then rotated.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   *
   *
   */

  class EUTelApplyAlignmentProcessor:public marlin::Processor {

  public:

//    virtual void modifyEvent( LCEvent * evt ) ;
//    virtual const std::string & name() const { return Processor::name() ; }

    //! Returns a new instance of EUTelApplyAlignmentProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelApplyAlignmentProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelApplyAlignmentProcessor;
    }

    //! Default constructor
    EUTelApplyAlignmentProcessor ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and reset all needed data
     *  members.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run the LCRunHeader of the this current run
     */
    void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. A few consistency
     *  checks are done on the first event to see if the alignment
     *  constants and the input collection are compatible.
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Apply Alignment
    /*!
     *
     */
    virtual void Direct(LCEvent *event);
    
    //! Apply Alignment in reverse direction
    /*!
     *
     */
    virtual void Reverse(LCEvent *event);
 
    //! Apply GEAR shits and rotations
    /*!
     *
     */
    virtual void ApplyGear6D(LCEvent *event);
 
    //! Revert GEAR shits and rotations
    /*!
     *
     */
    virtual void RevertGear6D(LCEvent *event);
 
    //! Apply alignment to a reference hit collection
    /*!
     *
     */
    virtual void AlignReferenceHit(EUTelEventImpl *evt,  EUTelAlignmentConstant * alignment );
   
    //!
    /*
     */
    virtual inline int GetApplyAlignmentDirection(){return _applyAlignmentDirection;}
     
    void TransformToLocalFrame(TrackerHitImpl* outputHit);
    void revertAlignment(double & x, double & y, double & z) ;
 
    //! Perform Euler rotations
    void _EulerRotation(double* _telPos, double* _gRotation);

    //! Perform Euler rotations backwards
    void _EulerRotationInverse(double* _telPos, double* _gRotation);

    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. For the time being there is
     *  nothing to check and do in this slot.
     *
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check (LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is finished. It
     *  prints only a goodbye message
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

    virtual    LCCollectionVec* CreateDummyReferenceHitCollection();

    virtual void CheckIOCollections(LCEvent* event);

  private:
    //! Conversion ID map.
    /*! In the data file, each cluster is tagged with a detector ID
     *  identify the sensor it belongs to. In the geometry
     *  description, there are along with the sensors also "passive"
     *  layers and other stuff. Those are identify by a layerindex. So
     *  we need a conversion table to go from the sensorID to the
     *  layerindex.
     */
    std::map< int, int > _conversionIdMap;

  	DISALLOW_COPY_AND_ASSIGN(EUTelApplyAlignmentProcessor)

  protected:

    // Internal processor variables
    int _nRun ;
    int _nEvt ;

    int _iDUT;
	int	_indexDUT;
	double	_xPitch, _yPitch, _rot00, _rot01, _rot10, _rot11;
  
    // class internal strings and pointers
    std::string      internal_inputHitCollectionName;
    LCCollectionVec* internal_inputCollectionVec      ;

    std::string      internal_referenceHitCollectionName;
    LCCollectionVec* internal_referenceHitVec;    

    //! Input collection name.
    /*! This is the name of the input hit collection.
     */
    std::string      _inputHitCollectionName;
    LCCollectionVec* _inputCollectionVec      ;

    //! Alignment constant collection name
    /*! This is the name of the collection containing the alignment
     *  constants. This should be the results of the execution of a
     *  EUTelMille and pede.
     */
    std::string      _alignmentCollectionName;
    LCCollectionVec* _alignmentCollectionVec  ;

    //! Output collection name.
    /*! This is the name of the output hit collection.
     */
    std::string      _outputHitCollectionName;
    LCCollectionVec* _outputCollectionVec     ;

    //! reference HitCollection name 
    /*!
     */
    bool        _applyToReferenceHitCollection;
 
    std::string _referenceHitCollectionName;
    LCCollectionVec* _referenceHitVec;    

    std::string _outputReferenceHitCollectionName;
    LCCollectionVec* _outputReferenceHitVec;    

    //! Correction method
    /*! There are actually several different
     *  methods to apply the alignment constants. Here below a list of
     *  available methods:
     *     0. <b>Shifts only</b> No rotations angle will be taken into
     *  account.
     *     1. <b>Rotation first</b> The rotational matrix is applied
     *  before the translation.
     *     2. <b>Translation first</b> The hits are first shifted and
     *  then rotated.
     */
    int _correctionMethod;

    //! apply alignment direction modes
    /*! There are basically two methods.
     *  Here below a list of available methods:
     *     0. <b>Direct </b> Normal alignment.
     *     1. <b>Reverse </b> Unalign everything back to how it was. 
     */
    int _applyAlignmentDirection;

    // 18 January 2011
    EVENT::StringVec		_alignmentCollectionNames;
    EVENT::StringVec		_alignmentCollectionSuffixes;

    // 10.09.2012
    EVENT::StringVec		_hitCollectionNames;
    EVENT::StringVec		_hitCollectionSuffixes;

    EVENT::StringVec		_refhitCollectionNames;
    EVENT::StringVec		_refhitCollectionSuffixes;

    //! DoGear
    bool            _doGear; 

    //! DoAlignCollection
    bool            _doAlignCollection; 

    //! Ignore _doGear and _doAlignCollection flags
    /*! set _doAlignmentInOneGo
     *  If you want to do all (anti)alignment steps in one go
     */ 
    bool            _doAlignmentInOneGo; 


    //! DEBUG
    bool            _debugSwitch; 
    //! rotation in ZY
    //
    double			_alpha;
    //! rotation in ZX
    //
    double			_beta;
    //! rotation in XY
    //
    double			_gamma;

    //! common DEBUG1 information
    /* set the number of events to have additional DEBUG1 information
     */
    int 			_printEvents;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Look Up Table for the sensor ID
    //    std::map< int, int > _lookUpTable;
    std::map< std::string, std::map< int, int > > _lookUpTable;

    //! boolean to mark the first processed event
    bool _fevent;

#if (defined(USE_AIDA) || defined(MARLIN_USE_AIDA))
    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;


    static std::string _densityPlotBeforeAlignName;
    static std::string _densityPlotAfterAlignName ;
    static std::string _hitHistoBeforeAlignName;
    static std::string _hitHistoAfterAlignName;

#endif

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

    //! An array with the Z position of planes
    double * _siPlaneZPosition;

    //! The ordered sensor ID vector
    /*! This vector contains the sensorID of all the detectors in the
     *  input collection in the same order as they appear. This vector
     *  has to be used to number the histogram booking and filling.
     */
    std::vector< int > _orderedSensorIDVec;

    //! Fill histogram switch
    /*! This boolean switch was initially introduced for debug reason
     *  but then we realized that it could stay there and protect
     *  against missing AIDA::Processor.
     *
     */
    bool _histogramSwitch;

  };

  //! A global instance of the processor
  EUTelApplyAlignmentProcessor gEUTelApplyAlignmentProcessor;

}
#endif
//#endif // GEAR
