export RUNRANGE=33-68
#export RUNRANGE=68

../../jobsub.py -c config.cfg -csv runlist.csv converter  $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv clustering $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv hitmaker   $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv align      $RUNRANGE
../../jobsub.py -c config.cfg -csv runlist.csv fitter     $RUNRANGE

