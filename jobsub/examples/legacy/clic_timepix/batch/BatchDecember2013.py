import os
from os import listdir
from os import walk
import time,os.path

def FetchIncompleteRuns(folder,BadRuns=[]):
    runs = []
    for (dirpath, dirnames, filenames) in walk(folder):     
        for file in filenames :
		if( os.path.getsize("%s/%s"%(folder,file)) < 15000) :
			tmp = int(''.join(x for x in file if x.isdigit()))
			if tmp not in BadRuns : 
				runs.append(tmp)

        break
    return runs


def GetRunNumber(file) : 
    return int(''.join(x for x in file if x.isdigit()))
	
def CreateConfigTelAlone(run) :
	f = open("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_telalone_batch.cfg")
	f2 =open("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_telalone_batch_run%06i.cfg"%run,"w")
	lines = f.readlines()
	for line in lines : 
		line=line.replace("@RunNumber@","%06i"%run)				
		f2.write(line)
	
	f.close()
	f2.close()
def CreateConfig(run) :
	f = open("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_batch.cfg")
	f2 =open("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_batch_run%06i.cfg"%run,"w")
	lines = f.readlines()
	for line in lines : 
		line=line.replace("@RunNumber@","%06i"%run)				
		f2.write(line)
	
	f.close()
	f2.close()


#runs = range(2032-2042) + range(2084,2105) + range(2159,2184)
runs = range(2159,2184) + range(2296,2312) + range(2356,2364) + range(2365,2380)
#runs = [3005]
print  runs

queue = "1nd"
batch_folder = "/afs/cern.ch/user/b/benoit/batch"
os.system("mkdir %s"%batch_folder)
os.system("rm -fr  %s/*"%batch_folder)
os.system("mkdir %s/launch"%batch_folder)
launch_folder = "%s/launch"%batch_folder
log_folder = "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_December2013_results/logs"

batch = []



for run in runs :
	CreateConfigTelAlone(run)
	CreateConfig(run)
	filename="%s/Batch_Run%i.sh"%(launch_folder,run)
	batch.append(filename)
	f=open(filename,'w')
	#f.write("cd /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation \n")
						
	#f.write("source /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/castor_test.sh \n")
	
	f.write("source /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/setup_eutelescope.sh \n")
	#f.write("cd /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub \n")	
		
	f.write("ls \n")	
	f.write("pwd \n")
		
	#File structure for castor 
	f.write("mkdir run%06i \n"%run)
	f.write("mkdir run%06i/db \n"%run)	
	f.write("mkdir run%06i/logs \n"%run)	
	f.write("mkdir run%06i/results \n"%run)
	f.write("mkdir run%06i/histo \n"%run)
	f.write("mkdir run%06i/lcio-raw \n"%run)
	f.write("mkdir run%06i/raw \n"%run)
	
	f.write("rfcp -v2 /castor/cern.ch/clicdet/CLIC_Vertex_TB_Data/DESY_TB_DATA_December2013/run%06i.raw run%06i/raw/run%06i.raw \n"%(run,run,run))	
	#f.write("cp /afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_December2013/run%06i.raw run%06i/raw \n"%(run,run))
	#f.write("cp /afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_December2013_results/db/run%06i-prealign-telalone-db.slcio run%06i/db \n"%(run,run))
	#f.write("cp /afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_December2013_results/db/run%06i-align-telalone-db.slcio run%06i/db \n"%(run,run))	
	f.write("ls \n")	
	
	#Align Telescope
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_telalone_batch_run%06i.cfg  converter %i \n"%(run,run))	
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_telalone_batch_run%06i.cfg  clusearch %i \n"%(run,run))
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_telalone_batch_run%06i.cfg  hitmaker %i \n"%(run,run))
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_telalone_batch_run%06i.cfg  align %i \n"%(run,run))
	
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_batch_run%06i.cfg  converter %i \n"%(run,run))	
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_batch_run%06i.cfg  clusearch %i \n"%(run,run))
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_batch_run%06i.cfg  hitmaker %i \n"%(run,run))
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configDecember2013_batch_run%06i.cfg  fitter %i \n"%(run,run))
	f.write("cp run%06i/histo/tbtrackrun%06i.root /afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_December2013_results/tbtrack \n"%(run,run))
	f.write("tar -pczvf run%06i.tar.gz run%06i \n"%(run,run))	
	
	f.write("rfcp run%06i.tar.gz /castor/cern.ch/clicdet/CLIC_Vertex_TB_results/CLIC_Vertex_TB_December2013_results/run%06i.tar.gz \n "%(run,run))
	#f.write("rm -fr run%06i.tar.gz run%06i/  \n"%(run,run))		

	f.write("ls \n")		
		
	os.system("chmod u+x %s"%filename)
	f.close()
	print filename

	

for job in batch :
    os.system("cd ~/batch/launch")
    run = GetRunNumber(job)
    log = "%s/Run%06i"%(log_folder,run)
    os.system("mkdir %s"%log)
    os.system("bsub -o %s/STDOUT -q %s %s"%(log,queue,job))
