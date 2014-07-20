import os 



runs = [1110,1174,1231]


for run in runs :
	os.system("./jobsub.py -c examples/clic_timepix/configOctober2013_telalone.cfg converter  %i"%run)
	os.system("./jobsub.py -c examples/clic_timepix/configOctober2013_telalone.cfg clusearch  %i"%run)
	os.system("./jobsub.py -c examples/clic_timepix/configOctober2013_telalone.cfg hitmaker  %i"%run)
	os.system("./jobsub.py -c examples/clic_timepix/configOctober2013_telalone.cfg align  %i"%run)	
	os.system("./jobsub.py -c examples/clic_timepix/configOctober2013_telalone.cfg fitter  %i"%run)
	
	
	
	
	
	
