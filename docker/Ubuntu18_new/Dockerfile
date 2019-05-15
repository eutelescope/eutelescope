FROM ubuntu:18.04

ENV ILCSOFT /ilcsoft

#install packages
RUN apt update && apt-get -y install libeigen3-dev build-essential gfortran cmake subversion default-jdk libgsl0-dev libxpm-dev libxft-dev libx11-dev libxext-dev python-dev git wget nano 

#get ilcInstall & installation
RUN mkdir -p $ILCSOFT
RUN git clone -b master https://github.com/eutelescope/ilcinstall $ILCSOFT/ilcinstall
RUN cd $ILCSOFT/ilcinstall && git status
RUN $ILCSOFT/ilcinstall/ilcsoft-install -i -v $ILCSOFT/ilcinstall/release-local-standalone.cfg

#RUN mkdir -p /data/TestExampleAconiteQuadFEI4
#RUN ["/bin/bash", "-c", "cd /data/TestExampleAconiteQuadFEI4 && wget -O run001085.raw http://cernbox.cern.ch/index.php/s/mEgF1o96bBIpwei/download && ls -l && cd -"]
#RUN sed -i 's@NativePath.*= \/afs\/desy\.de\/group\/telescopes\/EutelTestData\/TestExampleAconiteQuadFEI4@NativePath = \/data\/TestExampleAconiteQuadFEI4@g' $EUTELESCOPE/jobsub/examples/aconite-4chipLoc#al/config.cfg

RUN mkdir -p /data/TestExampleGBLLocal
RUN ["/bin/bash", "-c", "cd /data/TestExampleGBLLocal && wget -O run000117.raw https://cernbox.cern.ch/index.php/s/U68owtYMk7pdjcg/download && ls -l && cd -"]
RUN sed -i 's@NativePath.*= \/afs\/desy\.de\/group\/telescopes\/EutelTestData\/TestExampleDaturaNoDUT_GBL@NativePath = \/data\/TestExampleGBLLocal@g' $EUTELESCOPE/jobsub/examples/gbl_local/config.cfg

WORKDIR ${EUTELESCOPE}/build/
RUN make clean
CMD /bin/bash -c "echo $PR_NO; git fetch origin refs/pull/$PR_NO/merge; git checkout FETCH_HEAD; cmake ..; make install -j$(nproc); ctest -V -R gbl_local_test/*"
