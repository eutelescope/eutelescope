// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: ROOTProcessor.cc,v 1.1 2007-02-26 09:28:53 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifdef MARLIN_USE_ROOT

// marlin includes ".h"
#include "ROOTProcessor.h"
#include "marlin/Processor.h"

// root includes <.h>
#include <TObject.h>
#include <TList.h>
#include <TFile.h>

// system includes <>
#include <iostream>

namespace marlin { 
  
  ROOTProcessor aROOTProcessor ;
  
  
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
   
  ROOTProcessor* ROOTProcessor::_me ;

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
      _rootOutputFile->mkdir(list->GetName());
      _rootOutputFile->cd(list->GetName());
      list->Write();
      _rootOutputFile->cd("..");
    }
 
    _rootOutputFile->Close();
    
  }
  
  void ROOTProcessor::addTObject(const Processor * proc, TObject * obj) {
    
    // first of all check if this is the first time this processor was calling this method
    TList * list = dynamic_cast<TList*> (_me->_listOfList->FindObject(proc->name().c_str())) ;
    
    if ( list == 0x0 ) {
      // it means that this is the first time, ok no problem, create a
      // new one, make its name the same as the processor name and add
      // it to _listOfList, so next time will be already there...
      list = new TList();
      list->SetName(proc->name().c_str());
      _me->_listOfList->Add(list);
    }

    // add the object! That easy!!!
    list->Add(obj);

  }

  TObject * ROOTProcessor::getTObject(const Processor * proc, const char * objName) {

    // let me get the right list first
    TList * list = dynamic_cast<TList*> (_me->_listOfList->FindObject(proc->name().c_str())) ;
    
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

  TList * ROOTProcessor::getTList(const Processor * proc)  {
    
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

 
    

}; // namespace

#endif




