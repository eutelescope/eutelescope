#jobsub -c config.cfg -csv runlist.csv -g noisypixelmasker 501
#jobsub -c config.cfg -csv runlist.csv -g clustering 501
#jobsub -c config.cfg -csv runlist.csv -g hitmaker 501

#jobsub -c config.cfg -csv runlist.csv -g noisypixelmasker 502
jobsub -c config.cfg -csv runlist.csv -g clustering 502
jobsub -c config.cfg -csv runlist.csv -g hitmaker 502

#jobsub -c config.cfg -csv runlist.csv -g noisypixelmasker 503
jobsub -c config.cfg -csv runlist.csv -g clustering 503
jobsub -c config.cfg -csv runlist.csv -g hitmaker 503

#jobsub -c config.cfg -csv runlist.csv -g noisypixelmasker 504
jobsub -c config.cfg -csv runlist.csv -g clustering 504
jobsub -c config.cfg -csv runlist.csv -g hitmaker 504

exit 1

jobsub -c config.cfg -csv runlist.csv -g alignGBL1 501
jobsub -c config.cfg -csv runlist.csv -g alignGBL2 501
jobsub -c config.cfg -csv runlist.csv -g alignGBL3 501
jobsub -c config.cfg -csv runlist.csv -g fitGBL 501

jobsub -c config.cfg -csv runlist.csv -g noisypixelmasker 498
jobsub -c config.cfg -csv runlist.csv -g clustering 498
jobsub -c config.cfg -csv runlist.csv -g hitmaker 498
jobsub -c config.cfg -csv runlist.csv -g alignGBL1 498
jobsub -c config.cfg -csv runlist.csv -g alignGBL2 498
jobsub -c config.cfg -csv runlist.csv -g alignGBL3 498
jobsub -c config.cfg -csv runlist.csv -g fitGBL 498

jobsub -c config.cfg -csv runlist.csv -g noisypixelmasker 499
jobsub -c config.cfg -csv runlist.csv -g clustering 499
jobsub -c config.cfg -csv runlist.csv -g hitmaker 499
jobsub -c config.cfg -csv runlist.csv -g alignGBL1 499
jobsub -c config.cfg -csv runlist.csv -g alignGBL2 499
jobsub -c config.cfg -csv runlist.csv -g alignGBL3 499
jobsub -c config.cfg -csv runlist.csv -g fitgGBL 499
