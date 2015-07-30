*Use gear_250x50_SLAC.xml gearfile for converter and mapping
*Use gear_500x25_wrongsetup_SLAC.xml gearfile for clustering onward


*Run Converter using, FiringFreqCutAPIX = 0.015
*Run Mapping
*Switch Gearfile to 500x25_wrongsetup_SLAC.xml and run clustering
 
*Run hitlocal with,

ResidualsXMax =  100.  100.  100.  100.  100.  100.  100.  100.
ResidualsXMin = -100. -100. -100. -100. -100. -100. -100. -100.
ResidualsYMax =  100.  100.  100.  100.  100.  100.  100.  100.
ResidualsYMin = -100. -100. -100. -100. -100. -100. -100. -100.


*switch to gear_500x25_wrongsetup_SLAC_pre.xml

*change gear_500x25_wrongsetup_SLAC_pre.xml so as plane Z positions
correspond to the SLAC run 442 set up.

	   Plane ID:   Z Position:
	   0                0
	   1		    20
	   2		    40
	   21		    250
	   3		    370
	   4		    390
	   5		    410
	   20		    470

*change DUT ladder ID and sensitive ID to match. (both 20 or 21)

*add 180' rotation to DUT 20 in the ZX plane.

*change sensitive thickness for all 6 mimosa's to 50 microns.

*saved as gear_500x25_wrongsetup_SLAC_pre_adjusted.xml (change runlist to this gear file)

*run hitlocal and check correlation.

*run patternRecognition with DUT 21 excluded:

     TripletConnectDistCut =  0.5 0.5 
     TripletSlopeCuts = 0.01 0.01
     DoubletCenDistCut = 0.5 0.5
     DoubletDistCut = 0.5 0.5
     DUTWindow= 3
     excludeplanes=  21


*run GBLAlign with:

     FixXshifts=0  5
     FixYshifts=0  5
     FixZshifts=0  5
     FixXrot=   0 1 2 3 4 5 20 21  
     FixYrot=   0 1 2 3 4 5 20 21
     FixZrot=   0  5 

     r = 0.5

     dutX=0.150
     dutY=0.01

     xResolutionPlane        =  %(r)s %(r)s %(r)s %(dutX)s %(r)s %(r)s %(r)s %(dutX)s
     yResolutionPlane        =  %(r)s %(r)s %(r)s %(dutY)s %(r)s %(r)s %(r)s %(dutY)s 

*GBLAlign output should be:

	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The alignment parameters to update with these corrections:
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1  Type:  positionX  Value:  -0.81842E-02  Error:  0.73867E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1  Type:  rotationXY  Value:  0.59302E-02  Error:  0.18755E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2  Type:  rotationXY  Value:  0.47027E-02  Error:  0.18335E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3  Type:  positionY  Value:  0.45508E-01  Error:  0.13290E-01
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3  Type:  rotationXY  Value:  0.56014E-02  Error:  0.18182E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4  Type:  positionX  Value:  0.17893E-01  Error:  0.74819E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4  Type:  positionY  Value:  -0.22856E-01  Error:  0.13789E-01
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4  Type:  rotationXY  Value:  -0.23442E-02  Error:  0.18581E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20  Type:  positionX  Value:  -0.51281  Error:  0.66959E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20  Type:  positionY  Value:  -0.13511  Error:  0.89264E-02
	  17:22:40 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20  Type:  rotationXY  Value:  -0.92887E-02  Error:  0.15770E-02

