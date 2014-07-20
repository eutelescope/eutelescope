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
	

def CreateConfig(run) :
	f = open("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configOctober2013_batch.cfg")
	f2 =open("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configOctober2013_batch_run%06i.cfg"%run,"w")
	lines = f.readlines()
	for line in lines : 
		line=line.replace("@RunNumber@","%06i"%run)
		#line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001006-prealign-db.slcio")
		#line=line.replace("@Alignement@"," /afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001006-align-db.slcio")			
		
		if run in range(1005,1110) : 
			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001006-prealign-db.slcio")
			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001006-align-db.slcio")		
		
		elif run in range(1110,1172) : 
			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001110-prealign-db.slcio")
			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001110-align-db.slcio")		
		
		elif run in range(1172,1229) : 		
			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001174-prealign-db.slcio")
			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001174-align-db.slcio")		
		
		elif run in range(1231,1244) : 		
			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001231-prealign-db.slcio")
			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run001231-align-db.slcio")	

#		elif run in range(220,274) : 		
#			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000220-prealign-db.slcio")
#			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000220-align-db.slcio")					
#
#		elif run in range(274,293) : 	
#			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000274-prealign-db.slcio")
#			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000274-align-db.slcio")		
#
#		elif run in range(298,334) : 	
#			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000298-prealign-db.slcio")
#			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000298-align-db.slcio")		
#
#		elif run in range(404,411) : 	
#			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000404-prealign-db.slcio")
#			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000404-align-db.slcio")
#
#		elif run in range(413,866) : 	
#			line=line.replace("@Prealignment@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000413-prealign-db.slcio")
#			line=line.replace("@Alignement@","/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/db/run000413-align-db.slcio")			
				
		f2.write(line)
	
	f.close()
	f2.close()


#runs = range(48,400)
#runs = range(30,65) + range(97,154) + range(155,184) + range(195,220) + range(220,274)+ range(274,293) + range(298,334) + range(335,400)
#runs = [52,53,58,61,97,98,127,130,132,135,150,152,154,155,159,160,161,164,165,171,174,175,176,177,178,183,184,195,196,197,199,200,201,202,204,205,206,207,208,220,223,225,226,227,228,229,230,231,232,233,234,236,237,239,240,241,242,243,244,245,246,247,249,250,252,253,256,257,258,259,260,261,262,263] + range(269,278) + [280] + range(282,289) + [291,292] + range(298,310) + range(311,323)
#runs = runs + [324,325,326] + range(329,337) + [339,341,343,344] + range(347,356) + [358,359,360,362,363,365,366,367,369,370,371,372,374] + range(377,383) + [384,385,387,388,389,390,392,393,394,395,396,397,399] 

#runs = range(404,866)
#BadRuns = [1,2,3,4,5] + range(13,30) + range(31,38) + [40,41,42,44,46,47,99,128,129,147,148,149,156,157,158,163,169,179,180,181,182,406,407,409,412]+ range(65,96) + range(101,127) +range(185,195)+ range(209,219) +range(293,298) +range(434,444) + range(445,900)

BadRuns = [1017] + range(1037,1044) + [1050,1052,1056,1057] + range(1082,1108) + [1150] + range(1194,1203) 
#runs = FetchIncompleteRuns("/afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/tbtrack",BadRuns)
runs = range(1006,1244)
runs = [run for run in runs if run not in BadRuns]

print  runs

queue = "1nd"
batch_folder = "/afs/cern.ch/user/m/mbenoit/batch"
os.system("mkdir %s"%batch_folder)
os.system("rm -fr  %s/*"%batch_folder)
os.system("mkdir %s/launch"%batch_folder)
launch_folder = "%s/launch"%batch_folder
log_folder = "/afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_October2013_results/logs"

batch = []



for run in runs :
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
	
	f.write("rfcp -v2 /castor/cern.ch/clicdet/CLIC_Vertex_TB_Data/DESY_TB_DATA_October2013/run%06i.raw run%06i/raw/run%06i.raw \n"%(run,run,run))	
	f.write("ls \n")	

	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configOctober2013_batch_run%06i.cfg  converter %i \n"%(run,run))	
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configOctober2013_batch_run%06i.cfg  clusearch %i \n"%(run,run))
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configOctober2013_batch_run%06i.cfg  hitmaker %i \n"%(run,run))
	f.write("/afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/jobsub.py -c /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub/examples/clic_timepix/configOctober2013_batch_run%06i.cfg  fitter %i \n"%(run,run))
	f.write("cp run%06i/histo/tbtrackrun%06i.root /afs/cern.ch/eng/clic/TBData/DESY_TB_DATA_October2013_results/tbtrack \n"%(run,run))
	f.write("tar -pczf run%06i.tar.gz run%06i \n"%(run,run))	
	
	f.write("rfcp run%06i.tar.gz /castor/cern.ch/clicdet/CLIC_Vertex_TB_results_October2013/run%06i.tar.gz \n "%(run,run))
	#f.write("rm -fr run%06i.tar.gz run%06i/  \n"%(run,run))		

	f.write("ls \n")		
		
	os.system("chmod u+x %s"%filename)
	f.close()
	print filename


#for run in runs :
#	CreateConfig(run)
#	filename="%s/Batch_Run%i.sh"%(launch_folder,run)
#	batch.append(filename)
#	f=open(filename,'w')
#	f.write("cd /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation \n")
#						
#	f.write("source /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/setup_eutelescope.sh \n")
#	f.write("cd /afs/cern.ch/user/m/mbenoit/TestBeam_Analysis/LCD/EUTelescopeInstallation/Eutelescope/v00-08-00/jobsub \n")	
#	f.write("rfcp -v2 /castor/cern.ch/clicdet/CLIC_Vertex_TB_results/run%06i.tar.gz  run%06i.tar.gz \n"%(run,run))
#	f.write("tar -xvzf run%06i.tar.gz\n"%run)
#	f.write("./jobsub.py -c examples/clic_timepix/configAugust2013_batch_run%06i.cfg  fitter %i \n"%(run,run))		
#	f.write("cp run%06i/histo/tbtrack%06i.root /afs/cern.ch/eng/clic/data/DESY_TB_DATA_October2013_results/tbtrack/"%(run,run))
#	f.write("rfcp -v2 run%06i/histo/tbtrack%06i.root /castor/cern.ch/clicdet/CLIC_Vertex_TB_results \n "%(run,run))
#	f.write("rm -fr run%06i.tar.gz run%06i/  \n"%(run,run))		
#	
#	os.system("chmod u+rwx %s"%filename)
#	f.close()
	
	

for job in batch :
    os.system("cd ~/batch/launch")
    run = GetRunNumber(job)
    log = "%s/Run%06i"%(log_folder,run)
    os.system("mkdir %s"%log)
    os.system("bsub -o %s/STDOUT -q %s %s"%(log,queue,job))
