FROM ubuntu:16.04

ENV ILCSOFT /ilcsoft

RUN apt update && apt-get -y install libeigen3-dev build-essential gfortran cmake subversion default-jdk libgsl0-dev libxpm-dev libxft-dev libx11-dev libxext-dev python-dev git wget 

RUN mkdir -p $ILCSOFT
RUN git clone -b dev-base https://github.com/eutelescope/ilcinstall $ILCSOFT/ilcinstall
RUN cd $ILCSOFT/ilcinstall && git status
RUN $ILCSOFT/ilcinstall/ilcsoft-install -i -v $ILCSOFT/ilcinstall/examples/eutelescope/release-standalone.cfg

ENV EUTELESCOPE /ilcsoft/v01-19-02/Eutelescope/master/
ENV EUDAQ /ilcsoft/v01-19-02/Eutelescope/master/external/eudaq/v1.8.0/
ENV ILCSOFT_CMAKE_ENV /ilcsoft/v01-19-02/ILCSoft.cmake.env.sh
ENV MILLEPEDEII /ilcsoft/v01-19-02/Eutelescope/master/external/millepede2/tags/V04-03-09
ENV MILLEPEDEII_VERSION tags/V04-03-09
ENV GEAR /ilcsoft/v01-19-02/gear/v01-06-eutel-pre
ENV MARLIN /ilcsoft/v01-19-02/Marlin/v01-11
ENV MARLIN_DLL ${EUTELESCOPE}/lib/libEutelescope.so:${EUTELESCOPE}/lib/libEutelProcessors.so:${EUTELESCOPE}/lib/libEutelReaders.so:${EUDAQ}/lib/libNativeReader.so:${MARLIN_DLL}
ENV ROOTSYS /ilcsoft/v01-19-02/root/6.13.02
ENV GBL /ilcsoft/v01-19-02/GBL/V02-01-03
ENV CMAKE /ilcsoft/v01-19-02/CMake/3.11.2
ENV PATH ${CMAKE}/bin:${ROOTSYS}/bin:${MARLIN}/bin:${MILLEPEDEII}:${EUTELESCOPE}/bin:${GEAR}/tools:${GEAR}/bin:${PATH}
ENV LD_LIBRARY_PATH ${EUTELESCOPE}/lib:${ROOTSYS}/lib:${EUDAQ}/lib:${GEAR}/lib:${GBL}/lib:${LD_LIBRARY_PATH}

RUN mkdir -p /data/TestExampleAconiteQuadFEI4
RUN ["/bin/bash", "-c", "cd /data/TestExampleAconiteQuadFEI4 && wget -O run001085.raw http://cernbox.cern.ch/index.php/s/mEgF1o96bBIpwei/download && ls -l && cd -"]
RUN sed -i 's@NativePath.*= \/afs\/desy\.de\/group\/telescopes\/EutelTestData\/TestExampleAconiteQuadFEI4@NativePath = \/data\/TestExampleAconiteQuadFEI4@g' $EUTELESCOPE/jobsub/examples/aconite-4chipLocal/config.cfg

RUN mkdir -p /data/TestExampleGBLLocal
RUN ["/bin/bash", "-c", "cd /data/TestExampleGBLLocal && wget -O run000117.raw https://cernbox.cern.ch/index.php/s/U68owtYMk7pdjcg/download && ls -l && cd -"]
RUN sed -i 's@NativePath.*= \/afs\/desy\.de\/group\/telescopes\/EutelTestData\/TestExampleDaturaNoDUT_GBL@NativePath = \/data\/TestExampleGBLLocal@g' $EUTELESCOPE/jobsub/examples/gbl_local/config.cfg

WORKDIR ${EUTELESCOPE}/build/
RUN make clean
CMD /bin/bash -c "echo $PR_NO; git fetch origin refs/pull/$PR_NO/merge; git checkout FETCH_HEAD; cmake ..; make install -j 4; ctest -V -R aconite-4chipLocal_test/*; ctest -V -R gbl_local_test/*"
