FROM cern/slc6-base

ENV ILCSOFT /ilcsoft
RUN yum -y install sl-release-scl
RUN yum -y install git rh-git29 python27 devtoolset-6 eigen3-devel java-1.8.0-openjdk java-1.8.0-openjdk-devel gsl-devel libX11-devel libXpm-devel libXft-devel libXext-devel wget subversion tar 

RUN mkdir -p $ILCSOFT
RUN git clone -b dev-findpython https://github.com/eutelescope/ilcinstall $ILCSOFT/ilcinstall
RUN cd $ILCSOFT/ilcinstall && git status

ENV PATH $ILCSOFT/v01-19-02/CMake/3.6.3/bin/:${PATH}

RUN ["/bin/bash", "-c", "source /opt/rh/python27/enable && source /opt/rh/rh-git29/enable && source /opt/rh/devtoolset-6/enable && $ILCSOFT/ilcinstall/ilcsoft-install -i -v $ILCSOFT/ilcinstall/examples/eutelescope/release-standalone.cfg"]

ENV EUTELESCOPE /ilcsoft/v01-19-02/Eutelescope/master/
ENV EUDAQ /ilcsoft/v01-19-02/Eutelescope/master/external/eudaq/v1.7-dev/
ENV ILCSOFT_CMAKE_ENV /ilcsoft/v01-19-02/ILCSoft.cmake.env.sh
ENV MILLEPEDEII /ilcsoft/v01-19-02/Eutelescope/master/external/millepede2/tags/V04-03-03
ENV MILLEPEDEII_VERSION tags/V04-03-03
ENV GEAR /ilcsoft/v01-19-02/gear/v01-06-eutel-pre
ENV MARLIN /ilcsoft/v01-19-02/Marlin/v01-09
ENV MARLIN_DLL ${EUTELESCOPE}/lib/libEutelescope.so:${EUTELESCOPE}/lib/libEutelProcessors.so:${EUTELESCOPE}/lib/libEutelReaders.so:${EUDAQ}/lib/libNativeReader.so:${MARLIN_DLL}
ENV ROOTSYS /ilcsoft/v01-19-02/root/6.10.04
ENV GBL /ilcsoft/v01-19-02/GBL/V02-01-02
ENV PATH ${ROOTSYS}/bin:${MARLIN}/bin:${MILLEPEDEII}:${EUTELESCOPE}/bin:${GEAR}/tools:${GEAR}/bin:${PATH}
ENV LD_LIBRARY_PATH ${EUTELESCOPE}/lib:${ROOTSYS}/lib:${EUDAQ}/lib:${GEAR}/lib:${GBL}/lib:${LD_LIBRARY_PATH}

WORKDIR ${EUTELESCOPE}/jobsub/examples/aconite-4chipLocal/

RUN mkdir -p output/lcio output/logs output/histograms output/database
RUN mkdir -p /afs/
#CMD jobsub -c config.cfg -csv runlist.csv -g converter 1085
CMD /bin/bash -c "echo $PR_NO; git fetch origin refs/pull/$PR_NO/merge; git checkout FETCH_HEAD; cd $EUTELESCOPE/build; cmake ..; make install -j 4"

