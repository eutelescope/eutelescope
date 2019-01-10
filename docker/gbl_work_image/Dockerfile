FROM ubuntu:16.04
SHELL ["bash", "-c"]

ENV ILCSOFT /ilcsoft

#install packages
RUN apt update && apt-get -y install libeigen3-dev build-essential gfortran cmake subversion default-jdk libgsl0-dev libxpm-dev libxft-dev libx11-dev libxext-dev python-dev git wget emacs vim

#get ilcInstall
RUN mkdir -p $ILCSOFT
RUN git clone -b dev-base https://github.com/eutelescope/ilcinstall $ILCSOFT/ilcinstall
RUN cd $ILCSOFT/ilcinstall && git status

#add adapted script
ADD eutelescope.py eutelescope.py

#install EUTelescope
RUN mv eutelescope.py $ILCSOFT/ilcinstall/ilcsoft/eutelescope.py
RUN $ILCSOFT/ilcinstall/ilcsoft-install -i -v $ILCSOFT/ilcinstall/examples/eutelescope/release-standalone.cfg

ENV EUTELESCOPE /ilcsoft/v01-19-02/Eutelescope/master/
WORKDIR ${EUTELESCOPE}

#checkout gbl_work branch
RUN git checkout gbl_work && git pull

#build it
RUN chmod +x build_env.sh
RUN source build_env.sh && cd build && cmake .. && cmake --build . -- -j$(nproc) && make install

#change input path for files in examples
##GBL_noDUT
WORKDIR /ilcsoft/v01-19-02/Eutelescope/master/jobsub/examples/GBL_noDUT
RUN sed -i 's/\/afs\/desy/\/input\/desy/g' config.cfg
##GBL_SUT
WORKDIR /ilcsoft/v01-19-02/Eutelescope/master/jobsub/examples/GBL_SUT
RUN sed -i 's/\/afs\/desy/\/input\/desy/g' config.cfg
RUN sed -i 's/\/afs\/desy/\/input\/desy/g' config_DESY2017.cfg
##GBL_DUT
WORKDIR /ilcsoft/v01-19-02/Eutelescope/master/jobsub/examples/GBL_DUT
RUN sed -i 's/\/afs\/desy/\/input\/desy/g' config.cfg
##mixed-mode
WORKDIR /ilcsoft/v01-19-02/Eutelescope/master/jobsub/examples/mixed-mode
RUN sed -i 's/\/afs\/desy/\/input\/desy/g' duplicated/config.cfg
RUN sed -i 's/\/afs\/desy/\/input\/desy/g' standard/config.cfg


#RUN mkdir /afs /output /input
RUN mkdir /output /input
RUN chmod +rwx /
#RUN chmod +rwx /afs
RUN chmod +rwx /output
RUN chmod +rwx /dummy

#start environment (only on local machine possible)
CMD bash -c " cd /ilcsoft/v01-19-02/Eutelescope/master/ ; source build_env.sh ; cd - ; bash"