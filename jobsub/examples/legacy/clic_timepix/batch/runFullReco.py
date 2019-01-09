import os 



runs = range(3209,3215)


for run in runs :
	os.system("./jobsub.py -c examples/clic_timepix/configFebruary2014.cfg converter  %i"%run)
	os.system("./jobsub.py -c examples/clic_timepix/configFebruary2014.cfg clusearch  %i"%run)
	os.system("./jobsub.py -c examples/clic_timepix/configFebruary2014.cfg hitmaker  %i"%run)
	os.system("./jobsub.py -c examples/clic_timepix/configFebruary2014.cfg align  %i"%run)	
	os.system("./jobsub.py -c examples/clic_timepix/configFebruary2014.cfg fitter  %i"%run)
	
	
	
	
	
	
