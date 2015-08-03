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

#jobsub -c config/config_SLAC.cfg  -csv runlist/runlist.csv patternRecognition 442

#TripletConnectDistCut =  0.5 0.5 
#TripletSlopeCuts = 0.01 0.01
#DoubletCenDistCut = 0.5 0.5
#DoubletDistCut = 0.5 0.5
#DUTWindow= 4
#excludeplanes= 21  

#jobsub -c config/config_SLAC.cfg  -csv runlist/runlist.csv GBLAlign 442

FixXshifts=0  5 21   
FixYshifts=0  5 21  
FixZshifts=0  5 21  
FixXrot=   0 1 2 3 4 5 20 21  
FixYrot=   0 1 2 3 4 5 20 21
FixZrot=   0  5 21 
##Estimate resolution to pass to millepede.
#This will have to be increased if you have a large number of rejects.
rm26 = 0.003
dutXirr=100000
dutYirr=10000
dutXref=1
dutYref=1


18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  20  mode:  positionZ  result  -0.42950  error  0.99802
Removed!
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  20  mode:  rotationXY  result  -0.70009E-03  error
0.23742E-02     Removed!
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The alignment
parameters to update with these corrections:
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  positionX  Value:  -0.11755E-01  Error:  0.43936E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  positionY  Value:  -0.42912E-02  Error:  0.45064E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  positionZ  Value:  2.9373  Error:  0.70570E-01
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  rotationXY  Value:  0.57113E-02  Error:  0.10732E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionX  Value:  0.24282E-03  Error:  0.47205E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionY  Value:  -0.96906E-03  Error:  0.47164E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionZ  Value:  4.7031  Error:  0.75842E-01
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  rotationXY  Value:  0.43869E-02  Error:  0.11533E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionX  Value:  -0.10864E-01  Error:  0.47307E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionY  Value:  0.36426E-01  Error:  0.83289E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionZ  Value:  1.6135  Error:  0.75842E-01
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  rotationXY  Value:  0.53017E-02  Error:  0.11439E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionX  Value:  0.13338E-01  Error:  0.44433E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionY  Value:  -0.31936E-01  Error:  0.78601E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionZ  Value:  0.97644  Error:  0.70570E-01
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  rotationXY  Value:  -0.27258E-02  Error:  0.10635E-04
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  positionX  Value:  -0.55879  Error:  0.10184E-01
18:54:32 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  positionY  Value:  -0.13857  Error:  0.13959E-01

gear: alignedGear-iter1-run000442.xml
FixXshifts=0  5 21   
FixYshifts=0  5 21  
FixZshifts=0  5 21  
FixXrot=   0 1 2 3 4 5 20 21  
FixYrot=   0 1 2 3 4 5 20 21
FixZrot=   0  5 21 
##Estimate resolution to pass to millepede.
#This will have to be increased if you have a large number of rejects.
rm26 = 0.003
dutXirr=100000
dutYirr=10000
dutXref=0.5
dutYref=0.5




18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  1  mode:  positionX  result  -0.96437E-05  error  0.43193E-04
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  1  mode:  positionZ  result  -0.41158E-01  error  0.67841E-01
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  1  mode:  rotationXY  result  0.16643E-05  error  0.10538E-04
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  2  mode:  positionY  result  0.29232E-04  error  0.45848E-04
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  2  mode:  rotationXY  result  -0.24247E-06  error
0.11187E-04     Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  3  mode:  positionY  result  0.59751E-04  error  0.83882E-04
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  3  mode:  rotationXY  result  0.93424E-05  error  0.11512E-04
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  4  mode:  positionY  result  0.47489E-04  error  0.78633E-04
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  4  mode:  rotationXY  result  0.45537E-05  error  0.10661E-04
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  20  mode:  positionY  result  0.13626E-02  error  0.65130E-02
Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] INFO change smaller
than error. ID:  20  mode:  rotationXY  result  -0.92002E-03  error
0.11837E-02     Removed!
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The alignment
parameters to update with these corrections:
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  positionY  Value:  0.51932E-04  Error:  0.44316E-04
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionX  Value:  0.51190E-04  Error:  0.45838E-04
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionZ  Value:  -0.13309  Error:  0.72027E-01
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionX  Value:  0.14338E-03  Error:  0.47602E-04
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionZ  Value:  -0.14354  Error:  0.74709E-01
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionX  Value:  0.68917E-04  Error:  0.44649E-04
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionZ  Value:  -0.74608E-01  Error:  0.69243E-01
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  positionX  Value:  -0.74296E-02  Error:  0.51381E-02
18:58:34 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  positionZ  Value:  -1.0970  Error:  0.99174



gear: alignedGear-iter2-run000442.xml

#To continue this we need to think about the 180 degree rotation and the
corrections given by millepede.
Still have large offsets: alignedGear-iter2-fix20and21Shifts-run000442_pre.xml


Run pat rec again:

TripletConnectDistCut =  0.5 0.5 
TripletSlopeCuts = 0.01 0.01
DoubletCenDistCut = 0.5 0.5
DoubletDistCut = 0.5 0.5
DUTWindow= 0.5
excludeplanes=     


Now GBLAlign:

FixXshifts=0  5    
FixYshifts=0  5   
FixZshifts=0  5   
FixXrot=   0 1 2 3 4 5 20  21  
FixYrot=   0 1 2 3 4 5 20  21 
FixZrot=   0  5  
##Estimate resolution to pass to millepede.
#This will have to be increased if you have a large number of rejects.
rm26 = 0.003
dutXirr=0.08
dutYirr=0.03
dutXref=0.08
dutYref=0.03



19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  positionX  Value:  -0.16197E-02  Error:  0.28107E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  positionY  Value:  -0.15260E-02  Error:  0.31960E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  positionZ  Value:  0.37923  Error:  0.41125E-01
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  1
Type:  rotationXY  Value:  0.40000E-04  Error:  0.58509E-05
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionX  Value:  0.37781E-02  Error:  0.32125E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionY  Value:  0.29768E-03  Error:  0.35119E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  positionZ  Value:  0.65765  Error:  0.46954E-01
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  2
Type:  rotationXY  Value:  0.43979E-04  Error:  0.66865E-05
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionX  Value:  0.76799E-02  Error:  0.39137E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionY  Value:  0.14973E-01  Error:  0.53795E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  positionZ  Value:  0.51007  Error:  0.57018E-01
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  3
Type:  rotationXY  Value:  -0.82243E-05  Error:  0.80793E-05
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionX  Value:  0.40594E-02  Error:  0.30838E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionY  Value:  0.20415E-01  Error:  0.42514E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  positionZ  Value:  0.30072  Error:  0.44706E-01
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  4
Type:  rotationXY  Value:  -0.13004E-04  Error:  0.63165E-05
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  positionX  Value:  -0.72074E-02  Error:  0.54973E-03
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  positionY  Value:  0.16035  Error:  0.35425E-03
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  positionZ  Value:  12.705  Error:  0.37896
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  20
Type:  rotationXY  Value:  0.18256E-01  Error:  0.57677E-04
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  21
Type:  positionX  Value:  0.12685  Error:  0.11144E-02
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  21
Type:  positionY  Value:  0.14014E-01  Error:  0.54344E-03
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  21
Type:  positionZ  Value:  49.278  Error:  0.67998
19:12:55 jobsub.GBLAlign(INFO): [ MESSAGE9 "TrackAlign"] The sensor ID:  21
Type:  rotationXY  Value:  0.66527E-03  Error:  0.18042E-03

Gear output:alignedGear-iter3-run000442.xml
Do another iteration and final gear:alignedGear-iter4-run000442.xml
Irradiated senso does not look good clearly not random matching. Tighter cut
needed here most likely.


