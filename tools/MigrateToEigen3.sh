extpath=$EUTELESCOPE/external
tempdir="$extpath/packages.tmp" # where to store the .tar.gz files and company

if [ ! -d $tempdir ]; then mkdir $tempdir; fi;

echo "Removing old Eigen version..."
#rm old Eigen
rm $extpath/Eigen/*

eigenversion="3.2.2"
wget --no-check-certificate --output-document="$tempdir/${eigenversion}.tar.gz" 'http://bitbucket.org/eigen/eigen/get/'${eigenversion}.tar.gz
echo "Extracting tar archive..."
tar --strip-components 1 -C "$extpath/Eigen" -xzf "$tempdir/${eigenversion}.tar.gz" 
echo "... done with Eigen library"

# clean up
rm $tempdir/*
rmdir $tempdir
