#! /bin/zsh
#
###Make sure to enforce the right shell if not default
#$ -S /bin/zsh
#
###=======================================================
####Initial job configuration:
#
### Email address:
#$ -M eda.yildirim@desy.de
#
### When to send email (e.g on begin, end, suspension)
##$ -m beas
#$ -m as
#
### File for standard output and standard error
### stderr and stdout are merged together to stdout:
#$ -j y
### execute the job from the current directory and not relative to your home directory:
#$ -cwd
#### Your job name 
#$ -N RunFORMATTEDRUNNUM_itITERATION
#
# where to put the logfiles
#$ -o /nfs/dust/atlas/user/yeda/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/alibava/gblPython/output/runFORMATTEDRUNNUM/qsublogs/
### Job requirements:
#
### necessary cpu time:
#####$ -l h_cpu=80:00:00
#$ -l h_rt=05:30:00
#
### maximum necessary memory:
#$ -l h_vmem=8G
#
## filesize
#$ -l h_fsize=3G
#
### Choosing the operating system
#$ -l os=sld6
###======================================================================
basePath="/nfs/dust/atlas/user/yeda/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/alibava"
gblPythonPath="${basePath}/gblPython"
outputPath="${gblPythonPath}/output/runFORMATTEDRUNNUM"
qsubLogPath="${gblPythonPath}/output/runFORMATTEDRUNNUM/qsublogs"

mkdir -p ${qsubLogPath}
cd "${gblPythonPath}"
source PYROOT
/usr/bin/python "mytelmain.py" -g "GEARFILE" -r RUNNUM -e BEAMENERGY -c "GBLCUTS"
#/usr/bin/python "telmain.py" -g "${gblPythonPath}/GEARFILE" -r RUNNUM

cd "${outputPath}"
${basePath}/MillepedeII/pede steer.txt
/usr/bin/python "${basePath}/gearUpdate.py" -og "GEARFILE" -ng "${outputPath}/GEARFILE_NEW"

mv steer.txt steer_ITERATION.txt
mv cuts.root cuts_ITERATION.root
mv dutTree.root dutTree_ITERATION.root
mv milleBinary.dat milleBinary_ITERATION.dat
mv millepede.res millepede_ITERATION.res

