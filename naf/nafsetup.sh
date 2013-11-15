mkdir ErrorLogs OutputLogs XmlFiles Jobs JobSubmission
cp $EUTELESCOPE/naf/.setupfiles/MakeAllJobs.sh Jobs/
cp $EUTELESCOPE/naf/.setupfiles/CreateJobScripts.cc Jobs/
cp $EUTELESCOPE/naf/.setupfiles/CLEANUP JobSubmission/
cp $EUTELESCOPE/naf/.setupfiles/CodeToCreateSubmitJobs.cc JobSubmission/
cp $EUTELESCOPE/naf/.setupfiles/MakeXmlFiles XmlFiles/
cd Jobs
g++ CreateJobScripts.cc -O3 -o WriteJobFiles
cd ../JobSubmission
g++ CodeToCreateSubmitJobs.cc -O3 -o CreateSubmitJobs
cd ..
