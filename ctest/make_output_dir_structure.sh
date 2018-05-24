# Run the second argument in the directory created from the first argument
mkdir -p $1_test/output/logs
mkdir -p $1_test/output/histograms
mkdir -p $1_test/output/database
mkdir -p $1_test/output/lcio
mkdir -p $1_test/gear
cp $EUTELESCOPE/jobsub/examples/$1/gear/* $1_test/gear 
cd $1_test
exec $2
