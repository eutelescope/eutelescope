 
export ANALYSIS=/nfs/dust/cms/user/hjansen/ilcsoft/v01-17-05/Eutelescope/trunk/jobsub/examples/datura-noDUT/analysis-20mm/
cd $ANALYSIS
mkdir -p $ANALYSIS/output/histograms  && mkdir -p $ANALYSIS/output/database \
           && mkdir -p $ANALYSIS/output/logs \
           && mkdir -p $ANALYSIS/output/lcio
