FROM fedora:24

ENV ILCSOFT /ilcsoft

RUN dnf -y install wget automake
RUN wget ftp://rpmfind.net/linux/fedora/linux/releases/23/Everything/x86_64/os/Packages/g/gsl-devel-1.16-17.fc23.x86_64.rpm
RUN wget ftp://rpmfind.net/linux/fedora/linux/releases/23/Everything/x86_64/os/Packages/g/gsl-1.16-17.fc23.x86_64.rpm

RUN rpm -U gsl-1.16-17.fc23.x86_64.rpm
RUN rpm -U gsl-devel-1.16-17.fc23.x86_64.rpm

RUN dnf -y install gcc-c++ gcc-gfortran subversion java-1.8.0-openjdk-devel libX11-devel libXpm-devel libXft-devel libXext-devel python-devel tar make cmake which git 

RUN mkdir -p $ILCSOFT
RUN svn co https://svnsrv.desy.de/public/ilctools/ilcinstall/branches/v01-17-10-eutel-dev/ $ILCSOFT/ilcinstall_v01-17-10-eutel 
RUN $ILCSOFT/ilcinstall_v01-17-10-eutel/ilcsoft-install -i -v $ILCSOFT/ilcinstall_v01-17-10-eutel/examples/eutelescope/release-standalone.cfg

ENV EUTELESCOPE /ilcsoft/v01-17-10/Eutelescope/master/
ENV EUDAQ /ilcsoft/v01-17-10/Eutelescope/master/external/eudaq/v1.7-dev/
ENV ILCSOFT_CMAKE_ENV /ilcsoft/v01-17-10/ILCSoft.cmake.env.sh
ENV MILLEPEDEII /ilcsoft/v01-17-10/Eutelescope/master/external/millepede2/tags/V04-00-02
ENV GEAR /ilcsoft/v01-17-10/gear/v01-06-eutel-pre
ENV MARLIN /ilcsoft/v01-17-10/Marlin/v01-09
ENV MARLIN_DLL /ilcsoft/v01-17-10/Eutelescope/master/lib/libEutelescope.so:${EUTELESCOPE}/lib/libEutelProcessors.so:${EUTELESCOPE}/lib/libEutelReaders.so:${EUDAQ}/lib/libNativeReader.so:${MARLIN_DLL}

ENV ROOTSYS /ilcsoft/v01-17-10/root/6.08.00
ENV PATH /ilcsoft/v01-17-10/CMake/3.5.2/bin:${ROOTSYS}/bin:${MARLIN}/bin:${MILLEPEDEII}:${EUTELESCOPE}/bin:${GEAR}/tools:${GEAR}/bin:${PATH}
ENV LD_LIBRARY_PATH ${EUTELESCOPE}/lib:${ROOTSYS}/lib:${EUDAQ}/lib:${GEAR}/lib:${LD_LIBRARY_PATH}

WORKDIR ${EUTELESCOPE}/jobsub/examples/aconite-4chipLocal/

RUN mkdir -p output/lcio output/logs output/histograms output/database
RUN mkdir -p /afs/
#CMD jobsub -c config.cfg -csv runlist.csv -g converter 1085
CMD /bin/bash -c "echo $PR_NO; git fetch origin refs/pull/$PR_NO/merge; git checkout FETCH_HEAD; cd $EUTELESCOPE/build; cmake ..; make install -j 4"

