// -*- C++ -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifdef MARLIN_USE_ROOT

#ifndef ROOTPROCESSOR_H
#define ROOTPROCESSOR_H

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <lcio.h>

// forward declaration
class TList;
class TFile;
class TObject;


namespace marlin {


  //! Implementation of the ROOT histogramming interface
  /*! This Marlin processor is playing the role of a ROOT interface in
   *  the meaning that allows the user to use any ROOT objects and
   *  having them saved in the output ROOT file, grouped in folder
   *  named after the processor where the object have been generated.
   *
   *  It is offering the same functionalities of the AIDAProcessor
   *  even if ROOT is not making use of Factories owning all the
   *  produced objects.
   *
   *  This processor is actually built with Marlin only the user
   *  exported the enviromental variable MARLIN_USE_ROOT before
   *  starting. This is changing the building option in the top level
   *  GNUmakefile.
   *
   *  To use ROOT as histogramming package you need to add this
   *  processor to the list of active processors and it has to be
   *  first processor requiring access to ROOT TObject. Putting this
   *  processor as the first one has to be considered as a good
   *  attitude.
   *
   *  When the init() method is called, the ROOTProcessor singleton
   *  will open the ROOT output file for writing with the user
   *  specified name.  When the user wants to save a ROOT object in
   *  the file, he/she has to take care of creating the object and to
   *  use the addTObject(Processor*,TObject*) to add it to the list of
   *  objects being saved at the end.
   *
   *  <h4>Input - Prerequisites</h4>
   *  None
   * 
   *  <h4>Output</h4> 
   *  A ROOT file containing all the objects added by
   *  other processors.
   *
   *  @param ROOT file name
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: ROOTProcessor.h,v 1.1 2007-02-26 09:28:53 bulgheroni Exp $
   */

  class ROOTProcessor : public marlin::Processor {
    
  public:

    //! Default constructor
    /*! This is the default constructor. It is hidded among private
     *  members to avoid singleton breaking
     */ 
    ROOTProcessor();
    
    //! Return a new ROOTProcessor
    /*! ROOTProcessor needs to be a singleton, so in fact this is not
     *  returning a new processor, but returning the only static
     *  one. This is created right here if it doesn't exist.
     */ 
    virtual Processor*  newProcessor() ;
    

    //! Init method
    /*! This is called at the beginning, for the time being it is just
     *  printing out the processor parameters
     */
    virtual void init() ;

    //! Process run header
    /*! Nothing to do.
     */
    virtual void processRunHeader( LCRunHeader* /* run */ ) ;
    
    //! Process event 
    /*! Nothing to do.
     */
    virtual void processEvent( LCEvent * /* evt */) ; 

    //! Check
    /*! Nothing to do.
     */
    virtual void check( LCEvent * /* evt */ ) ; 


    //! End method
    /*! This is called at the end of the data processing. This is the
     *  place where all TObjects are saved
     */ 
    virtual void end() ;

    //! Add TObject to the output file
    /*! This static member is used to add any kind of TObject to the
     *  list of objects are going to be saved into the output ROOT
     *  files. 
     * 
     *  @param proc The processor is calling the addTObject method. So
     *  usually this.
     *  @param obj The TObject we want to add
     */ 
    static void addTObject(const Processor * proc, TObject * obj);

    //! Get a TObject
    /*! This static member is used to get back the pointer of a
     *  certain objcted. The object is recognized by its name.
     *
     *  @param proc A pointer to the processor calling the method
     *  @param name The name of the object we want to retrieve
     *  @return A pointer to the object we would like to have
     *
     *  @todo Consider the possibility to define and throw an
     *  exception if the object the user is looking for doesn't exist
     */ 
    static TObject * getTObject(const Processor * proc, const char* objName) ;

    //! Get the TList of TObjects
    /*! This static member is used to get back a pointer to the list
     *  of objects created within a certain processor
     *  
     *  @param proc A pointer to the processor calling the method
     *  @return A pointer to the TList containing all objects added by
     *  this processor
     */
    static TList * getTList(const Processor * proc) ;

  protected:
    
    //! The output ROOT file name 
    /*! This is the name of the output file. If the user wants to save
     *  the file into a different folder, the full path has to be
     *  specified in the file name. If the output file is saved on a
     *  remote host via a supported networking protocol (rootd,
     *  proofd, xrootd), then it has to specified in the file name as
     *  well.
     */ 
    std::string _fileName ;
  
private:
    
    //! Static instance of the singleton
    static ROOTProcessor* _me ;

    //! A pointer to the output file
    /*! This a pointer to the output ROOT file.
     */ 
    TFile * _rootOutputFile;

    //! A list of list
    /*! This is the global list of TObjects. It contains a list for
     *  each processor that is callig the addTObject static method.
     */
    TList * _listOfList;


} ;

} // end namespace marlin 
#endif

#endif



