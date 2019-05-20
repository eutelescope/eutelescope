/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

//GEAR is required for this processor
#if defined(USE_GEAR)

// eutelescope includes ".h"
#include "EUTelCorrelator.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTELESCOPE.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelVirtualCluster.h"

// marlin includes ".h"
#include "marlin/Global.h"
#include "marlin/Processor.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <marlin/AIDAProcessor.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <Exceptions.h>

// system includes <>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace eutelescope;

EUTelCorrelator::EUTelCorrelator ():Processor ("EUTelCorrelator"),
_sensorIDVec ()
{

  _description = "EUTelCorrelator fills histograms with correlation plots";

  EVENT::StringVec _clusterCollectionVecExample;

  registerInputCollections (LCIO::TRACKERPULSE,
			    "InputClusterCollections",
			    "List of cluster collections",
			    _clusterCollectionVec,
			    _clusterCollectionVecExample);

  registerInputCollection (LCIO::TRACKERHIT,
			   "InputHitCollectionName",
			   "Hit collection name",
			   _inputHitCollectionName, std::string ("hit"));

  registerProcessorParameter ("ClusterChargeMinimum",
			      "Minimum allowed cluster charge to be taken into "
			      "account for the correlation plots (default = 2)",
			      _clusterChargeMin, 2);

  registerProcessorParameter ("RequiredEvents",
			      "How many events are needed to get "
			      "reasonable correlation plots (and "
			      "Offset DB)? (default=1000)",
			      _requiredEvents, 1000);

  registerOptionalParameter ("FixedPlane",
			     "SensorID of fixed plane", _fixedPlaneID, 0);

  registerOptionalParameter ("ExcludedPlanes",
			     "The list of sensor IDs that shall be excluded (e.g. passive plane).",
			     _excludedPlaneIDVec, std::vector < int >());

  registerOptionalParameter ("ResidualsXMin",
			     "Minimal values of the hit residuals in the X direction "
			     "for a correlation band (ordered along z of sensors).",
			     _residualsXMin, std::vector < float >(6, -10.));

  registerOptionalParameter ("ResidualsYMin",
			     "Minimal values of the hit residuals in the Y direction "
			     "for a correlation band (ordered along z of sensors).",
			     _residualsYMin, std::vector < float >(6, -10.));

  registerOptionalParameter ("ResidualsXMax",
			     "Maximal values of the hit residuals in the X direction "
			     "for a correlation band (ordered along z of sensors).",
			     _residualsXMax, std::vector < float >(6, 10.));

  registerOptionalParameter ("ResidualsYMax",
			     "Maximal values of the hit residuals in the Y direction "
			     "for a correlation band (ordered along z of sensors).",
			     _residualsYMax, std::vector < float >(6, 10.));

  registerOptionalParameter ("MinNumberOfCorrelatedHits",
			     "If there are more then this number of correlated "
			     "hits (planes->track candidate) (default=5)",
			     _minNumberOfCorrelatedHits, 5);
}

void
EUTelCorrelator::init ()
{

  //usually a good idea to do
  printParameters ();

  //initialize geometry
  geo::gGeometry ().initializeTGeoDescription (EUTELESCOPE::GEOFILENAME,
					       EUTELESCOPE::DUMPGEOROOT);

  //get sensor ID vector and subtract excluded IDs
  _sensorIDVec = geo::gGeometry ().sensorIDsVec ();
for (auto excludeID:_excludedPlaneIDVec)
    {
      _sensorIDVec.
	erase (std::
	       remove (_sensorIDVec.begin (), _sensorIDVec.end (), excludeID),
	       _sensorIDVec.end ());
    }

  //fill sensorIDtoZ map
  for (std::vector < int >::iterator it = _sensorIDVec.begin ();
       it != _sensorIDVec.end (); it++)
    {
      _sensorIDtoZ.insert (std::make_pair (*it,
					   static_cast <
					   int >(it -
						 _sensorIDVec.begin ())));
    }

  //reset run and event counters
  _iRun = 0;
  _iEvt = 0;

  //set initalization flag
  _isInitialize = false;
}

void
EUTelCorrelator::processRunHeader (LCRunHeader * rdr)
{

  EUTelRunHeaderImpl *runHeader = new EUTelRunHeaderImpl (rdr);

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file. If the numbers are different, instead of barely
  // quitting ask the user what to do.
  if (runHeader->getGeoID () == 0)
    streamlog_out (WARNING0)
      << "The geometry ID in the run header is set to zero." << std::endl
      << "This may mean that the GeoID parameter was not set" << std::endl;

  if (static_cast < unsigned int >(runHeader->getGeoID ()) !=
      geo::gGeometry ().getLayoutID ())
    {
      streamlog_out (WARNING5)
	<< "Error during the geometry consistency check: " << std::endl
	<< "The run header says the GeoID is " << runHeader->
	getGeoID () << std::
	endl << "The GEAR description says is     " << geo::gGeometry ().
	getLayoutID () << std::endl;
    }

  delete runHeader;

  //increment run counter
  ++_iRun;
}

void
EUTelCorrelator::processEvent (LCEvent * event)
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  //required events reached, then stop
  if (_iEvt > _requiredEvents)
    return;

  //increment event counter
  ++_iEvt;

  //check event type
  EUTelEventImpl *evt = static_cast < EUTelEventImpl * >(event);
  if (evt->getEventType () == kEORE)
    {
      streamlog_out (DEBUG4) << "EORE found: nothing else to do." << std::
	endl;
      return;
    }

  //intialize booleans
  _hasClusterCollection = false;
  _hasHitCollection = false;

  //[START] loop over cluster collections
  for (size_t icoll = 0; icoll < _clusterCollectionVec.size (); icoll++)
    {

      std::string _inputClusterCollectionName = _clusterCollectionVec[icoll];
      try
      {
	//check if cluster collections exist
	event->getCollection (_inputClusterCollectionName);
	_hasClusterCollection = true;
	streamlog_out (DEBUG5) << "found cluster collection " << icoll
	  << " with name "
	  << _inputClusterCollectionName.c_str () << std::endl;

      }
      catch (lcio::Exception & e)
      {
	_hasClusterCollection = false;
	streamlog_out (WARNING) << "NOT found cluster collection " << icoll
	  << " with name "
	  << _inputClusterCollectionName.c_str () << std::endl;
	break;
      }
    }				//[END] loop over cluster collections

  try
  {
    //check if hit collections exist
    event->getCollection (_inputHitCollectionName);
    _hasHitCollection = true;
    streamlog_out (DEBUG5) << "found hit collection with name "
      << _inputHitCollectionName.c_str () << std::endl;

  }
  catch (lcio::Exception & e)
  {

    _hasHitCollection = false;
    streamlog_out (DEBUG5) << "NOT found hit collection with name "
      << _inputHitCollectionName.c_str () << std::endl;
  }

  //if there is no collection at all, stop for this event
  if (!_hasClusterCollection && !_hasHitCollection)
    {
      return;
    }

  //if not initialized, book histograms
  if (!_isInitialize)
    {
      bookHistos ();
      _isInitialize = true;
    }

  //[IF] hasCluster
  if (_hasClusterCollection && !_hasHitCollection)
    {

      //[START] loop over collection (external)
      for (size_t eCol = 0; eCol < _clusterCollectionVec.size (); eCol++)
	{

	  std::string externalInputClusterCollectionName =
	    _clusterCollectionVec[eCol];

	  LCCollectionVec *externalInputClusterCollection =
	    static_cast <
	    LCCollectionVec *
	    >(event->getCollection (externalInputClusterCollectionName));
	  CellIDDecoder < TrackerPulseImpl >
	    pulseCellDecoder (externalInputClusterCollection);

	  //[START] loop over cluster (external)
	  for (size_t iExt = 0;
	       iExt < externalInputClusterCollection->size (); ++iExt)
	    {

	      TrackerPulseImpl *externalPulse =
		static_cast <
		TrackerPulseImpl *
		>(externalInputClusterCollection->getElementAt (iExt));

	      EUTelVirtualCluster *externalCluster;

	      ClusterType type =
		static_cast < ClusterType > (static_cast <
					     int
					     >((pulseCellDecoder
						(externalPulse)["type"])));

	      //check that the type of cluster is ok
	      if (type == kEUTelDFFClusterImpl)
		{
		  externalCluster =
		    new EUTelDFFClusterImpl (static_cast <
					     TrackerDataImpl *
					     >(externalPulse->
					       getTrackerData ()));
		}
	      else if (type == kEUTelBrickedClusterImpl)
		{
		  externalCluster =
		    new EUTelBrickedClusterImpl (static_cast <
						 TrackerDataImpl *
						 >(externalPulse->
						   getTrackerData ()));
		}
	      else if (type == kEUTelFFClusterImpl)
		{
		  externalCluster =
		    new EUTelFFClusterImpl (static_cast <
					    TrackerDataImpl *
					    >(externalPulse->
					      getTrackerData ()));
		}
	      else if (type == kEUTelSparseClusterImpl)
		{
		  externalCluster =
		    new EUTelSparseClusterImpl < EUTelGenericSparsePixel >
		    (static_cast <
		     TrackerDataImpl * >(externalPulse->getTrackerData ()));
		  if (externalCluster != nullptr
		      && externalCluster->getTotalCharge () <
		      _clusterChargeMin)
		    {
		      delete externalCluster;
		      continue;
		    }
		}
	      else
		{
		  continue;
		}

	      int externalSensorID =
		pulseCellDecoder (externalPulse)["sensorID"];

	      streamlog_out (DEBUG1) << "externalSensorID : " <<
		externalSensorID << " externalCluster=" << externalCluster <<
		std::endl;

	      //get coordinates of external seed
	      float externalXCenter = 0.;
	      float externalYCenter = 0.;
	      externalCluster->getCenterOfGravity (externalXCenter,
						   externalYCenter);

	      //check minimal charge requirement
	      if (externalCluster->getTotalCharge () <= _clusterChargeMin)
		{
		  delete externalCluster;
		  continue;
		}

	      //[START] loop over collection (internal)
	      for (size_t iCol = 0; iCol < _clusterCollectionVec.size ();
		   iCol++)
		{

		  std::string internalInputClusterCollectionName =
		    _clusterCollectionVec[iCol];

		  LCCollectionVec *internalInputClusterCollection =
		    static_cast <
		    LCCollectionVec *
		    >(event->
		      getCollection (internalInputClusterCollectionName));
		  CellIDDecoder < TrackerPulseImpl >
		    pulseCellDecoder (internalInputClusterCollection);

		  //[START] loop over cluster (internal)
		  for (size_t iInt = 0;
		       iInt < internalInputClusterCollection->size (); ++iInt)
		    {

		      TrackerPulseImpl *internalPulse =
			static_cast <
			TrackerPulseImpl *
			>(internalInputClusterCollection->
			  getElementAt (iInt));

		      EUTelVirtualCluster *internalCluster;

		      ClusterType type =
			static_cast < ClusterType > (static_cast <
						     int
						     >((pulseCellDecoder
							(internalPulse)
							["type"])));

		      //check that the type of cluster is ok
		      if (type == kEUTelDFFClusterImpl)
			{
			  internalCluster =
			    new EUTelDFFClusterImpl (static_cast <
						     TrackerDataImpl *
						     >(internalPulse->
						       getTrackerData ()));
			}
		      else if (type == kEUTelBrickedClusterImpl)
			{
			  internalCluster =
			    new EUTelBrickedClusterImpl (static_cast <
							 TrackerDataImpl *
							 >(internalPulse->
							   getTrackerData
							   ()));
			}
		      else if (type == kEUTelFFClusterImpl)
			{
			  internalCluster =
			    new EUTelFFClusterImpl (static_cast <
						    TrackerDataImpl *
						    >(internalPulse->
						      getTrackerData ()));
			}
		      else if (type == kEUTelSparseClusterImpl)
			{
			  internalCluster =
			    new EUTelSparseClusterImpl <
			    EUTelGenericSparsePixel > (static_cast <
						       TrackerDataImpl *
						       >(internalPulse->
							 getTrackerData ()));
			  if (internalCluster != nullptr
			      && internalCluster->getTotalCharge () <
			      _clusterChargeMin)
			    {
			      delete internalCluster;
			      continue;
			    }
			}
		      else
			{
			  continue;
			}

		      //check charge requirement
		      if (internalCluster->getTotalCharge () <
			  _clusterChargeMin)
			{
			  delete internalCluster;
			  continue;
			}

		      int internalSensorID =
			pulseCellDecoder (internalPulse)["sensorID"];

		      if ((internalSensorID != getFixedPlaneID () &&
			   externalSensorID == getFixedPlaneID ()) ||
			  (_sensorIDtoZ.at (internalSensorID) >
			   _sensorIDtoZ.at (externalSensorID)))
			{

			  //get coordinates of internal seed
			  float internalXCenter = 0.;
			  float internalYCenter = 0.;
			  internalCluster->
			    getCenterOfGravity (internalXCenter,
						internalYCenter);

			  streamlog_out (DEBUG5) << "Filling histo for "
			    << "extID " << externalSensorID << " and intID "
			    << internalSensorID << std::endl;

			  //input coordinates in correlation matrix (for X and Y)
			  _clusterXCorrelationMatrix[externalSensorID]
			    [internalSensorID]->fill (externalXCenter,
						      internalXCenter);
			  _clusterYCorrelationMatrix[externalSensorID]
			    [internalSensorID]->fill (externalYCenter,
						      internalYCenter);
			  streamlog_out (MESSAGE1) << " ex " <<
			    externalSensorID << " = [" << externalXCenter <<
			    ":" << externalYCenter << "]" << " in " <<
			    internalSensorID << " = [" << internalXCenter <<
			    ":" << internalYCenter << "]" << std::endl;
			}

		      delete internalCluster;
		    }		//[END] loop over cluster (internal)
		}		//[END] loop over collection (internal)

	      delete externalCluster;
	    }			//[END] loop over cluster (external)
	}			//[END] loop over collection (external)
    }				//[ENDIF] hasCluster

  //[IF] hasCollection
  if (_hasHitCollection)
    {

      LCCollectionVec *inputHitCollection =
	static_cast <
	LCCollectionVec * >(event->getCollection (_inputHitCollectionName));
      UTIL::CellIDDecoder < TrackerHitImpl >
	hitDecoder (EUTELESCOPE::HITENCODING);

      streamlog_out (MESSAGE2) << "inputHitCollection "
	<< _inputHitCollectionName.c_str () << std::endl;

      //[START] loop over collection (external) 
      for (size_t iExt = 0; iExt < inputHitCollection->size (); ++iExt)
	{

	  std::vector < double >trackXVec;
	  std::vector < double >trackYVec;
	  std::vector < int >planeIDVec;
	  trackXVec.clear ();
	  trackYVec.clear ();
	  planeIDVec.clear ();

	  //get external hit
	  TrackerHitImpl *externalHit =
	    static_cast <
	    TrackerHitImpl * >(inputHitCollection->getElementAt (iExt));
	  double *externalPosition =
	    const_cast < double *>(externalHit->getPosition ());
	  int externalSensorID = hitDecoder (externalHit)["sensorID"];
	  double etrackPointLocal[] =
	    { externalPosition[0], externalPosition[1],
	    externalPosition[2]
	  };
	  double etrackPointGlobal[] =
	    { externalPosition[0], externalPosition[1],
	    externalPosition[2]
	  };

	  //check for coordinate system
	  if (hitDecoder (externalHit)["properties"] != kHitInGlobalCoord)
	    {
	      //transfer to global frame
	      geo::gGeometry ().local2Master (externalSensorID,
					      etrackPointLocal,
					      etrackPointGlobal);
	    }
	  else
	    {
	      //do nothing, already in global telescope frame
	    }

	  trackXVec.push_back (etrackPointGlobal[0]);
	  trackYVec.push_back (etrackPointGlobal[1]);
	  planeIDVec.push_back (externalSensorID);

	  streamlog_out (MESSAGE2)
	    << "extPlane:" << externalSensorID << " at local position: " <<
	    etrackPointLocal[0] << " " << etrackPointLocal[1] <<
	    " and global position: " << etrackPointGlobal[0] << " " <<
	    etrackPointGlobal[1] << std::endl;

	  //[START] loop over collection (internal)
	  for (size_t iInt = 0; iInt < inputHitCollection->size (); ++iInt)
	    {

	      TrackerHitImpl *internalHit =
		static_cast <
		TrackerHitImpl * >(inputHitCollection->getElementAt (iInt));

	      double *internalPosition =
		const_cast < double *>(internalHit->getPosition ());
	      int internalSensorID = hitDecoder (internalHit)["sensorID"];
	      double itrackPointLocal[] =
		{ internalPosition[0], internalPosition[1],
		internalPosition[2]
	      };
	      double itrackPointGlobal[] =
		{ internalPosition[0], internalPosition[1],
		internalPosition[2]
	      };

	      //check for coordinate system
	      if (hitDecoder (internalHit)["properties"] != kHitInGlobalCoord)
		{
		  //transfer to global frame
		  geo::gGeometry ().local2Master (internalSensorID,
						  itrackPointLocal,
						  itrackPointGlobal);
		}
	      else
		{
		  //do nothing, already in global telescope frame
		}

	      //[IF] check planes
	      if ((internalSensorID != getFixedPlaneID () &&
		   externalSensorID == getFixedPlaneID ()) ||
		  (_sensorIDtoZ.at (internalSensorID) >
		   _sensorIDtoZ.at (externalSensorID)))
		{

		  int iz = _sensorIDtoZ.at (internalSensorID);

		  //[IF] check residual requirement
		  if (((etrackPointGlobal[0] - itrackPointGlobal[0]) <
		       _residualsXMax[iz]) &&
		      (_residualsXMin[iz] <
		       (etrackPointGlobal[0] - itrackPointGlobal[0])) &&
		      ((etrackPointGlobal[1] - itrackPointGlobal[1]) <
		       _residualsYMax[iz]) &&
		      (_residualsYMin[iz] <
		       (etrackPointGlobal[1] - itrackPointGlobal[1])))
		    {

		      trackXVec.push_back (itrackPointGlobal[0]);
		      trackYVec.push_back (itrackPointGlobal[1]);
		      planeIDVec.push_back (internalSensorID);

		      streamlog_out (MESSAGE2) << "intPlane:" <<
			internalSensorID << " at loc position: " <<
			itrackPointLocal[0] << " " << itrackPointLocal[1] <<
			" and global position: " << itrackPointGlobal[0] <<
			" " << itrackPointGlobal[1] << std::endl;
		    }		//[ENDIF] check residual requirement
		}		//[ENDIF] check planes
	    }			//[END] loop over collection (internal)

	  std::vector < int >planeID_unique = planeIDVec;

	  //remove duplicates
	  std::vector < int >::iterator p_end;

	  p_end = unique (planeID_unique.begin (), planeID_unique.end ());
	  //FIXME: should there be real erasing, change function... (jha92)
	  //planeID_unique.erase(p_end,planeID_unique.end());

	  //[IF] check for minimal number of correlated hits
	  if (static_cast < int >(planeID_unique.size ()) >
	      _minNumberOfCorrelatedHits
	      && trackXVec.size () == trackYVec.size ())
	    {

	      size_t indexPlane = 0;

	      for (size_t i = 0; i < trackXVec.size (); i++)
		{
		  if (i == indexPlane)
		    continue;	//skip as this one is not booked

		  _hitXCorrelationMatrix[planeIDVec[indexPlane]][planeIDVec
								 [i]]->
		    fill (trackXVec[indexPlane], trackXVec[i]);
		  _hitYCorrelationMatrix[planeIDVec[indexPlane]][planeIDVec
								 [i]]->
		    fill (trackYVec[indexPlane], trackYVec[i]);
		  //assumption: all rotations were done in hitmaker processor
		  _hitXCorrShiftMatrix[planeIDVec[indexPlane]][planeIDVec
							       [i]]->
		    fill (trackXVec[indexPlane],
			  trackXVec[indexPlane] - trackXVec[i]);
		  _hitYCorrShiftMatrix[planeIDVec[indexPlane]][planeIDVec
							       [i]]->
		    fill (trackYVec[indexPlane],
			  trackYVec[indexPlane] - trackYVec[i]);
		}
	    }			//[ENDIF]
	}			//[END] loop over collection (external)
    }				//[ENDIF] hasCollection
#endif
}

void
EUTelCorrelator::end ()
{

  streamlog_out (MESSAGE4) << "Successfully finished" << std::endl;
}

void
EUTelCorrelator::bookHistos ()
{

  //if no cluster and hits, no histograms needed
  if (!_hasClusterCollection && !_hasHitCollection)
    return;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  try
  {
    streamlog_out (DEBUG5) << "Booking histograms" << std::endl;

    //declare variables
    int xBin, yBin;
    double xMin, xMax;
    double yMin, yMax;

    //[IF] cluster correlation
    if (_hasClusterCollection && !_hasHitCollection)
      {

	//create directories
	marlin::AIDAProcessor::tree (this)->mkdir ("ClusterX");
	marlin::AIDAProcessor::tree (this)->mkdir ("ClusterY");

	//[START] loop over sensors (from)
      for (auto fromID:_sensorIDVec)
	  {

	    std::map < unsigned int, AIDA::IHistogram2D * >innerMapXCluster;
	    std::map < unsigned int, AIDA::IHistogram2D * >innerMapYCluster;

	    //[START] loop over sensors (to)
	  for (auto toID:_sensorIDVec)
	      {

		if ((toID != getFixedPlaneID ()
		     && fromID == getFixedPlaneID ())
		    || (_sensorIDtoZ.at (toID) > _sensorIDtoZ.at (fromID)))
		  {

		    //create 2D histogram: cluster correlation X
		    std::string histName_clusterXCorr =
		      "ClusterX/ClusterXCorrelation_d" + to_string (fromID) +
		      "_d" + to_string (toID);

		    xBin = geo::gGeometry ().getPlaneNumberOfPixelsX (fromID);
		    xMin = 0.;
		    xMax = geo::gGeometry ().getPlaneNumberOfPixelsX (fromID);
		    yBin = geo::gGeometry ().getPlaneNumberOfPixelsX (toID);
		    yMin = 0.;
		    yMax = geo::gGeometry ().getPlaneNumberOfPixelsX (toID);

		    AIDA::IHistogram2D * hist2D_clusterXCorr =
		      marlin::AIDAProcessor::histogramFactory (this)->
		      createHistogram2D (histName_clusterXCorr, xBin, xMin,
					 xMax, yBin, yMin, yMax);

		    hist2D_clusterXCorr->
		      setTitle ("Cluster correlation in X (d" +
				std::to_string (fromID) + "->d" +
				std::to_string (toID) + "); X_d" +
				std::to_string (fromID) + " [mm]; X_d" +
				std::to_string (toID) + " [mm]");

		    innerMapXCluster[toID] = hist2D_clusterXCorr;

		    //create 2D histogram: cluster correlation Y
		    std::string histName_clusterYCorr =
		      "ClusterY/ClusterYCorrelation_d" + to_string (fromID) +
		      "_d" + to_string (toID);

		    xBin = geo::gGeometry ().getPlaneNumberOfPixelsY (fromID);
		    xMin = 0.;
		    xMax = geo::gGeometry ().getPlaneNumberOfPixelsY (fromID);
		    yBin = geo::gGeometry ().getPlaneNumberOfPixelsY (toID);
		    yMin = 0.;
		    yMax = geo::gGeometry ().getPlaneNumberOfPixelsY (toID);

		    AIDA::IHistogram2D * hist2D_clusterYCorr =
		      marlin::AIDAProcessor::histogramFactory (this)->
		      createHistogram2D (histName_clusterYCorr, xBin, xMin,
					 xMax, yBin, yMin, yMax);

		    hist2D_clusterYCorr->
		      setTitle ("Cluster correlation in Y (d" +
				std::to_string (fromID) + "->d" +
				std::to_string (toID) + "); Y_d" +
				std::to_string (fromID) + " [mm]; Y_d" +
				std::to_string (toID) + " [mm]");

		    innerMapYCluster[toID] = hist2D_clusterYCorr;
		  }
		else
		  {
		    innerMapXCluster[toID] = nullptr;
		    innerMapYCluster[toID] = nullptr;
		  }
	      }			//[END] loop over sensors (to)
	    _clusterXCorrelationMatrix[fromID] = innerMapXCluster;
	    _clusterYCorrelationMatrix[fromID] = innerMapYCluster;
	  }			//[END] loop over sensors (from)
      }				//[ENDIF] cluster correlation

    //[IF] hit correlation
    if (_hasHitCollection)
      {

	//create directories
	marlin::AIDAProcessor::tree (this)->mkdir ("HitX");
	marlin::AIDAProcessor::tree (this)->mkdir ("HitY");
	marlin::AIDAProcessor::tree (this)->mkdir ("HitXShift");
	marlin::AIDAProcessor::tree (this)->mkdir ("HitYShift");

	//[START] loop over sensors (from)
      for (auto fromID:_sensorIDVec)
	  {

	    std::map < unsigned int, AIDA::IHistogram2D * >innerMapXHit;
	    std::map < unsigned int, AIDA::IHistogram2D * >innerMapYHit;
	    std::map < unsigned int, AIDA::IHistogram2D * >innerMapXHitShift;
	    std::map < unsigned int, AIDA::IHistogram2D * >innerMapYHitShift;

	    //[START] loop over sensors (to)
	  for (auto toID:_sensorIDVec)
	      {

		if ((toID != getFixedPlaneID ()
		     && fromID == getFixedPlaneID ())
		    || (_sensorIDtoZ.at (toID) > _sensorIDtoZ.at (fromID)))
		  {

		    //create 2D histogram: hit correlation X
		    std::string histName_hitXCorr = "HitX/HitXCorrelation_d" +
		      to_string (fromID) + "_d" + to_string (toID);

		    xBin = 100;
		    xMin = -0.5 * geo::gGeometry ().getPlaneXSize (fromID);
		    xMax = 0.5 * geo::gGeometry ().getPlaneXSize (fromID);
		    yBin = 100;
		    yMin = -0.5 * geo::gGeometry ().getPlaneXSize (toID);
		    yMax = 0.5 * geo::gGeometry ().getPlaneXSize (toID);

		    AIDA::IHistogram2D * hist2D_hitXCorr =
		      marlin::AIDAProcessor::histogramFactory (this)->
		      createHistogram2D (histName_hitXCorr, xBin, xMin, xMax,
					 yBin, yMin, yMax);

		    hist2D_hitXCorr->setTitle ("Hit correlation in X (d" +
					       std::to_string (fromID) +
					       "->d" + std::to_string (toID) +
					       "); X_d" +
					       std::to_string (fromID) +
					       " [mm]; X_d" +
					       std::to_string (toID) +
					       " [mm]");

		    innerMapXHit[toID] = hist2D_hitXCorr;

		    //create 2D histogram: hit correlation X shift
		    std::string histName_hitXCorrShift =
		      "HitXShift/HitXCorrShift_d" + to_string (fromID) +
		      "_d" + to_string (toID);

		    AIDA::IHistogram2D * hist2D_hitXCorrShift =
		      marlin::AIDAProcessor::histogramFactory (this)->
		      createHistogram2D (histName_hitXCorrShift, xBin, xMin,
					 xMax, yBin, yMin, yMax);

		    hist2D_hitXCorrShift->
		      setTitle ("Hit correlation shift in X (d" +
				std::to_string (fromID) + "->d" +
				std::to_string (toID) + "); X_d" +
				std::to_string (fromID) + " [mm]; X_d" +
				std::to_string (fromID) + "-X_d" +
				std::to_string (toID) + " [mm]");

		    innerMapXHitShift[toID] = hist2D_hitXCorrShift;

		    //create 2D histogram: hit correlation Y
		    std::string histName_hitYCorr = "HitY/HitYCorrelation_d" +
		      to_string (fromID) + "_d" + to_string (toID);
		    xBin = 100;
		    xMin = -0.5 * geo::gGeometry ().getPlaneYSize (fromID);
		    xMax = 0.5 * geo::gGeometry ().getPlaneYSize (fromID);
		    yBin = 100;
		    yMin = -0.5 * geo::gGeometry ().getPlaneYSize (toID);
		    yMax = 0.5 * geo::gGeometry ().getPlaneYSize (toID);

		    AIDA::IHistogram2D * hist2D_hitYCorr =
		      marlin::AIDAProcessor::histogramFactory (this)->
		      createHistogram2D (histName_hitYCorr, xBin, xMin, xMax,
					 yBin, yMin, yMax);

		    hist2D_hitYCorr->setTitle ("Hit correlation in Y (d" +
					       std::to_string (fromID) +
					       "->d" + std::to_string (toID) +
					       "); Y_d" +
					       std::to_string (fromID) +
					       " [mm]; Y_d" +
					       std::to_string (toID) +
					       " [mm]");

		    innerMapYHit[toID] = hist2D_hitYCorr;

		    //create 2D histogram: hit correlation Y shift
		    std::string histName_hitYCorrShift =
		      "HitYShift/HitYCorrShift_d" + to_string (fromID) +
		      "_d" + to_string (toID);

		    AIDA::IHistogram2D * hist2D_hitYCorrShift =
		      marlin::AIDAProcessor::histogramFactory (this)->
		      createHistogram2D (histName_hitYCorrShift, xBin, xMin,
					 xMax, yBin, yMin, yMax);

		    hist2D_hitYCorrShift->
		      setTitle ("Hit correlation shift in Y (d" +
				std::to_string (fromID) + "->d" +
				std::to_string (toID) + "); Y_d" +
				std::to_string (fromID) + " [mm]; Y_d" +
				std::to_string (fromID) + "-Y_d" +
				std::to_string (toID) + " [mm]");

		    innerMapYHitShift[toID] = hist2D_hitYCorrShift;
		  }
		else
		  {

		    innerMapXHit[toID] = nullptr;
		    innerMapYHit[toID] = nullptr;
		    innerMapXHitShift[toID] = nullptr;
		    innerMapYHitShift[toID] = nullptr;
		  }
	      }			//[END] loop over sensors (to)

	    _hitXCorrelationMatrix[fromID] = innerMapXHit;
	    _hitYCorrelationMatrix[fromID] = innerMapYHit;
	    _hitXCorrShiftMatrix[fromID] = innerMapXHitShift;
	    _hitYCorrShiftMatrix[fromID] = innerMapYHitShift;
	  }			//[END] loop over sensors (from)
      }				//[ENDIF] hit correlation

  }
  catch (lcio::Exception & e)
  {

    streamlog_out (ERROR1)
      << "No AIDAProcessor initialized. Sorry for quitting..." << endl;
    exit (-1);
  }
#endif
}

std::vector < double >
EUTelCorrelator::guessSensorOffset (int internalSensorID,
				    int externalSensorID,
				    std::vector < double >cluCenter)
{
  double internalXCenter = cluCenter.at (0);
  double internalYCenter = cluCenter.at (1);
  double externalXCenter = cluCenter.at (2);
  double externalYCenter = cluCenter.at (3);

  double xDet_in =
    internalXCenter * geo::gGeometry ().getPlaneXPitch (internalSensorID);
  double yDet_in =
    internalYCenter * geo::gGeometry ().getPlaneYPitch (internalSensorID);
  double xDet_ex =
    externalXCenter * geo::gGeometry ().getPlaneXPitch (externalSensorID);
  double yDet_ex =
    externalYCenter * geo::gGeometry ().getPlaneYPitch (externalSensorID);

  double xCoo_in = internalXCenter;
  double yCoo_in = internalYCenter;

  //get rotated sensors coordinates (in mm or um)
  double xPos_in = xDet_in;
  double yPos_in = yDet_in;
  double xPos_ex = xDet_ex;
  double yPos_ex = yDet_ex;

  double xCooPos_in = xCoo_in;
  double yCooPos_in = yCoo_in;

  //get rotated sensor coordinates (only in pixel number: col number)
  xPos_in += geo::gGeometry ().getPlaneXPosition (internalSensorID) +
    geo::gGeometry ().getPlaneXSize (internalSensorID) / 2.;

  xCooPos_in = xCooPos_in;

  yPos_in += geo::gGeometry ().getPlaneYPosition (internalSensorID) +
    geo::gGeometry ().getPlaneYSize (internalSensorID) / 2.;

  yCooPos_in = yCooPos_in;

  xPos_ex += geo::gGeometry ().getPlaneXPosition (externalSensorID) +
    geo::gGeometry ().getPlaneXSize (externalSensorID) / 2.;

  yPos_ex += geo::gGeometry ().getPlaneYPosition (externalSensorID) +
    geo::gGeometry ().getPlaneYSize (externalSensorID) / 2.;

  std::vector < double >cluster_offset;

  cluster_offset.push_back (-(xPos_in - xPos_ex));
  cluster_offset.push_back (-(yPos_in - yPos_ex));

  //add also internal sensor X and Y coord
  cluster_offset.push_back (xCooPos_in);
  cluster_offset.push_back (yCooPos_in);

  return cluster_offset;
}

#endif // USE_GEAR
