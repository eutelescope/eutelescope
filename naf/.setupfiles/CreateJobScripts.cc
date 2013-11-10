#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;


vector< int > FillRuns(string runrange);
int StringToInteger(string thestring);
void WriteFile(vector< int > runs, string analysistype);


int main(int argc, char *argv[]){
  if(argc != 3){
    //explain how to use program
    cerr << "Argument 1: Analysis Type (converter, hitmaker, fitter etc.)" << endl << "Argument 2: Run Range (1-10,12,13,25-36 etc.)" << endl;
    terminate();
  }
  cout << "There are " << argc << " arguments passed." << endl;
  cout << "The analysis step to be turned into a job is: " << argv[1] << endl;
  cout << "To be carried out on the runs in the range: " << argv[2] << endl;
  
  string analysistype = argv[1];
  string runrange = argv[2];
  vector< int > runs = FillRuns(runrange);
  WriteFile(runs,analysistype);
}

vector< int > FillRuns(string runrange){
  vector< int > runs;
  size_t currentpos(0);
  do{  
    size_t poshyphen = runrange.find('-',currentpos);
    size_t poscomma = runrange.find(',',currentpos);
    cout << "currentpos = " << currentpos << endl << "poscomma = " << poscomma << endl << "poshyphen = " << poshyphen << endl;
    if(poscomma < poshyphen && poscomma != static_cast< size_t >(-1)){
      cout << "Running through poscomma < poshyphen && poscomma != static_cast< size_t >(-1)" << endl;
      if(poscomma <= 0 || poscomma <= currentpos){
        cerr << "Run range specified is wrong, terminating because poscomma <= 0 || poscomma [" << poscomma << "] <= [" << currentpos << "] currentpos" << endl;
        terminate();
      }
      string stringrun = runrange.substr(currentpos,poscomma-currentpos);
      int intrun = StringToInteger(stringrun);
      cout << "Extracting run " << intrun << endl;
      runs.push_back(intrun);
      currentpos = poscomma + 1;
      cout << "Updated values:" << endl << "currentpos = " << currentpos << endl << endl;
    } else if(poshyphen < poscomma && poshyphen != static_cast< size_t >(-1)){
      cout << "Running through poshyphen < poscomma && poshyphen != static_cast< size_t >(-1)" << endl;
      if(poshyphen <= 0 || poshyphen <= currentpos){
        cerr << "Run range specified is wrong, terminating because poshyphen <= 0 || poshyphen <= currentpos" << endl;
        terminate();
      }
      string stringrun1 = runrange.substr(currentpos,poshyphen-currentpos);
      size_t nextcomma = runrange.find(',',poshyphen+1);
      size_t nexthyphen = runrange.find('-',poshyphen+1);
      cout << "stringrun1 = " << stringrun1 << endl;
      cout << "nextcomma = " << nextcomma << endl;
      cout << "nexthyphen = " << nexthyphen << endl;
      if(nexthyphen < nextcomma){
        cerr << "Can not have two hyphens next to each other like this 1-3-5" << endl;
        terminate();
      }
      string stringrun2("");
      if(nextcomma == static_cast< size_t >(-1)){
        stringrun2 = runrange.substr(poshyphen+1,string::npos);
      } else{
        stringrun2 = runrange.substr(poshyphen+1,poshyphen+1+nextcomma);
      }
      cout << "stringrun2 = " << stringrun2 << endl;
      int run1 = StringToInteger(stringrun1);
      int run2 = StringToInteger(stringrun2);
      for(int i = run1; i <= run2; ++i){
        cout << "Pushing back run " << i << " into the vector" << endl;
        runs.push_back(i);
      }
      currentpos = nextcomma+1;
    } else if(poscomma == static_cast< size_t >(-1) && poshyphen == static_cast< size_t >(-1)){
      string run1 = runrange.substr(currentpos,string::npos);
      cout << "Final run to be read in is " << run1 << endl;
      int finalrun = StringToInteger(run1);
      runs.push_back(finalrun);
      break;
    } else{
      cerr << "Somehow the position of the comma is in the same place as the hyphen, terminating" << endl;
      terminate();
    }   
  } while (currentpos < runrange.size() && currentpos != 0);
  return runs;
}     
      
int StringToInteger(string thestring){
  stringstream x;
  x << thestring; 
  int result(0);
  x >> result;
  return result;
}     
    
void WriteFile(vector< int > runs, string analysistype){
  ofstream file;
  for(vector< int >::iterator i = runs.begin(); i != runs.end(); ++i){
    stringstream ss;
    ss << "Run" << *i << analysistype;
    string numberextension("");
    if(*i < 10) numberextension = "00000";
    else if(*i < 100) numberextension = "0000";
    else if(*i < 1000) numberextension = "000";
    else if(*i < 10000) numberextension = "00";
    else if(*i < 100000) numberextension = "0";
    file.open(ss.str().c_str(), ofstream::trunc);
    file << "#!/bin/bash" << endl;
    file << "uname -a" << endl;
    file << "if [ -z \"$EUTELESCOPE\" ]; then source /nfs/dust/atlas/user/hamnett/ilcsoft/v01-17-02/Eutelescope/v00-08-03/build_env.sh; fi" << endl;
    if(analysistype == "FullChain"){
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv converter  " << *i << endl;
      file << "Marlin converter-" << numberextension.c_str() << *i << ".xml" << endl;
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv clustering  " << *i << endl;
      file << "Marlin clustering-" << numberextension.c_str() << *i << ".xml" << endl;
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv filter  " << *i << endl;
      file << "Marlin filter-" << numberextension.c_str() << *i << ".xml" << endl;
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv hitmaker  " << *i << endl;
      file << "Marlin hitmaker-" << numberextension.c_str() << *i << ".xml" << endl;
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv aligndaf  " << *i << endl;
      file << "Marlin aligndaf-" << numberextension.c_str() << *i << ".xml" << endl;
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv fitter  " << *i << endl;
      file << "Marlin fitter-" << numberextension.c_str() << *i << ".xml" << endl;
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv x0  " << *i << endl;
      file << "Marlin x0-" << numberextension.c_str() << *i << ".xml" << endl;
    } else{
//      file << "jobsub -c $EUTELESCOPE/jobsub/examples/X0/configNoDUTNEW.cfg -csv $EUTELESCOPE/jobsub/examples/X0/runs.csv  " << analysistype.c_str() << " " << *i << endl;
      file << "Marlin " << analysistype.c_str() << "-" << numberextension.c_str() << *i << ".xml" << endl;
    }
    file.close();
  }   
}

