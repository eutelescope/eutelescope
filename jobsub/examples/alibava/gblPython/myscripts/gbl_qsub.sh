
runnum="$1"
gearfile="$2"
beamenergy="$3"
gblcuts="$4"
iteration="$5"

gbldir="/nfs/dust/atlas/user/yeda/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/alibava/gblPython"
mkdir -p "${gbldir}/output/run${(l:6::0:)runnum}/qsublogs"
qsubfile="${gbldir}/qsubtemplates/qsub_${runnum}_it${iteration}.sh"

cp "${gbldir}/qsubtemplates/qsubtemplate.sh" ${qsubfile}
sed -i "s@FORMATTEDRUNNUM@${(l:6::0:)runnum}@g" ${qsubfile}   
sed -i "s@GEARFILE_NEW@gear_${runnum}_it${iteration}.xml@g" ${qsubfile}   
sed -i "s@GEARFILE@${gearfile}@g" ${qsubfile}
sed -i "s@RUNNUM@${runnum}@g" ${qsubfile}
sed -i "s@ITERATION@${iteration}@g" ${qsubfile}
sed -i "s@BEAMENERGY@${beamenergy}@g" ${qsubfile}
sed -i "s@GBLCUTS@${gblcuts}@g" ${qsubfile}
qsub ${qsubfile}

