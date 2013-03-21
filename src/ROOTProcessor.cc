// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)

// marlin includes ".h"
#include "ROOTProcessor.h"
#include "marlin/Processor.h"

// root includes <.h>
#include <TObject.h>
#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <TDirectory.h>
#include <TObjArray.h>
#include <TObjString.h>

// system includes <>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;
using namespace eutelescope;
using namespace marlin;

namespace eutelescope {

  //! A global instance
  ROOTProcessor aROOTProcessor ;

  //! The ROOTProcessor singleton
  /*! Like for the AIDAProcessor, also the ROOTProcessor is a
   *  singleton class instanciated only once for each execution. This
   *  is to avoid of having multiple ROOT histogram interface running
   *  together and offer the advantage of interoperability among the
   *  different active processors.
   */
  ROOTProcessor* ROOTProcessor::_me ;
}

ROOTProcessor::ROOTProcessor() : Processor("ROOTProcessor"),
                                 _rootOutputFile(NULL)  {

  _listOfList = new TList();

  _description = "Processor that handles ROOT files. Creates one directory per processor. "
    " Within process you can create all TObjects you wish and you only need to use ROOTProcessor::addTObject() to have it save"
    " Needs to be the first ActiveProcessor" ;

  registerProcessorParameter( "FileName" ,
                              "The output ROOT file name"  ,
                              _fileName ,
                              std::string("file.root") ) ;
}



Processor*  ROOTProcessor::newProcessor() {
  if( ! _me )
    _me = new ROOTProcessor ;
  return _me ;
}

void ROOTProcessor::init() {

  printParameters() ;

}


void ROOTProcessor::processRunHeader( LCRunHeader* /*run*/) { /*NO-OP*/ ; }

void ROOTProcessor::processEvent( LCEvent * /* evt */ )  { /*NO-OP*/ ; }

void ROOTProcessor::check( LCEvent * /* evt */  ) { /*NO-OP*/ ; }

void ROOTProcessor::end(){

  _rootOutputFile = TFile::Open(_fileName.c_str(), "RECREATE");

  TIter nextList(_listOfList);
  while ( TList * list = dynamic_cast<TList*> (nextList()) ) {

    // try to see if the full name has some "/" inside if so, it
    // means that this list has to be saved into a separate folder
    TString fullname = list->GetName();
    if ( fullname.First("/") != -1 ) {
      // ok it needs to be saved into a separate folder
      TDirectory * baseDir = gDirectory;
      TObjArray *  array   = fullname.Tokenize("/");
      for (int iDir = 0; iDir < array->GetEntriesFast(); iDir++) {
        TString      dir     = ( dynamic_cast< TObjString* >( array->At(iDir)))->GetString();
        // if the directory doens't exist, create it!
        if ( !baseDir->cd(dir) ) {
          baseDir->mkdir(dir);
          baseDir->cd(dir);
        }
        baseDir = gDirectory;
      }
      list->Write();
      _rootOutputFile->cd(_rootOutputFile->GetPath());
    } else {
      _rootOutputFile->mkdir(list->GetName());
      _rootOutputFile->cd(list->GetName());
      list->Write();
      _rootOutputFile->cd("..");
    }

    _rootOutputFile->Close();

  }
}

void ROOTProcessor::addTObject(const marlin::Processor * proc, TObject * obj) {

  _me->addTObject(proc,"",obj);

}

void ROOTProcessor::addTObject(const marlin::Processor * proc, const char * subfolder, TObject * obj) {

  // first of all I need to prepare the the proper list name that is proc->name() + "/" + subfolder
  stringstream ss;
  if ( string(subfolder).size() == 0 ) {
    ss << proc->name() ;
  } else {
    ss << proc->name() << "/" << subfolder;
  }
  string listName = ss.str();

  TList * list = dynamic_cast<TList*> (_me->_listOfList->FindObject(listName.c_str()));

  if ( list == 0x0 ) {
    list = new TList();
    list->SetName(listName.c_str());
    _me->_listOfList->Add(list);
  }

  list->Add(obj);

}



TObject * ROOTProcessor::getTObject(const marlin::Processor * proc, const char * objName) {

  return _me->getTObject(proc, "", objName);
}


TObject * ROOTProcessor::getTObject(const marlin::Processor * proc, const char * subfolder, const char * objName) {

  // first of all I need to prepare the the proper list name that is proc->name() + "/" + subfolder
  stringstream ss;
  if ( string(subfolder).size() == 0 ) {
    ss << proc->name() ;
  } else {
    ss << proc->name() << "/" << subfolder;
  }
  string listName = ss.str();

  TList * list = dynamic_cast<TList*> (_me->_listOfList->FindObject(listName.c_str()));

  if ( list == 0x0 ) {
    // something is really going badly wrong, you are looking for an
    // object but it is probably from another processor, or not
    // created at all. Returning 0x0, hoping you are check the
    // return value before going on. Consider the possibility to
    // throw an exception
    return 0x0;
  }

  return list->FindObject(objName);
}


TList * ROOTProcessor::getTList(const marlin::Processor * proc)  {

  TList * list = dynamic_cast<TList*> (_me->_listOfList->FindObject(proc->name().c_str())) ;

  if ( list == 0x0 ) {
    // it means that this is the first time, ok no problem, create a
    // new one, make its name the same as the processor name and add
    // it to _listOfList, so next time will be already there...
    list = new TList();
    list->SetName(proc->name().c_str());
    _me->_listOfList->Add(list);
  }

  return list;
}





#endif




