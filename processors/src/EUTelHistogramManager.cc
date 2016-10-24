// Author:  Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version: $Id$

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// #ifdef EXPERIMENTAL
// personal includes ".h"
#include "EUTelHistogramManager.h"

// marlin includes ".h"
#include "marlin/Exceptions.h"
#include "marlin/tinyxml.h"

// lcio includes <.h>

// system includes
#include <string>
#include <map>
#include <exception>
#include <iostream>


using namespace std;
using namespace marlin;
using namespace eutelescope;




EUTelHistogramManager::~EUTelHistogramManager() {

  map< string, EUTelHistogramInfo * >::iterator iter = _histoInfoMap.begin();
  while ( iter != _histoInfoMap.end() ) {
    delete iter->second;
    ++iter;
  }
  _histoInfoMap.clear();

}

bool EUTelHistogramManager::init() throw( std::exception, marlin::ParseException ) {

  
  TiXmlDocument * doc = new TiXmlDocument;

  if ( doc->LoadFile( _histoInfoFileName ) ) {

    TiXmlHandle    hDoc(doc);
    TiXmlHandle    hRoot(0);
    TiXmlElement * pElem;
    
    pElem = hDoc.FirstChildElement().Element();
    if ( !pElem ) {
      delete doc;
      throw ParseException( string( "EUTelHistogramManager::init: no root tag <HistogramManager> ... </HistogramManager> "
				    "found in ") + _histoInfoFileName );
    } else {
      hRoot = TiXmlHandle(pElem);
    }
	

    _histoInfoMap.clear();

    TiXmlNode * nHistosBlock = hRoot.FirstChild( "histos" ).ToNode();
    if ( !nHistosBlock ) {
      delete doc;
      throw ParseException( string( "EUTelHistogramManager::init: no <histos> ... </histos> block found in ")  + _histoInfoFileName);
    }

    TiXmlElement * pHistoNode = hRoot.FirstChild( "histos" ).FirstChild( "histo" ).Element();
    while ( pHistoNode ) {

      EUTelHistogramInfo * histoInfo = new EUTelHistogramInfo;
      histoInfo->_name  = pHistoNode->Attribute("name");
      
      if ( pHistoNode->Attribute("type") == NULL ) {
	delete doc;
	delete histoInfo;
	throw ParseException( string( "EUTelHistogramManager::init: no type found for " +  string(pHistoNode->Attribute("name"))) );
      } else   histoInfo->_type  = pHistoNode->Attribute("type");
      
      if ( pHistoNode->Attribute("title") == NULL  ) histoInfo->_title = "";
      else histoInfo->_title = pHistoNode->Attribute("title");
      
      if ( ( histoInfo->_type != string("C1D") ) &&
	   ( histoInfo->_type != string("C2D") ) &&
	   ( histoInfo->_type != string("C3D") ) ) {
	pHistoNode->QueryIntAttribute("xBin", &(histoInfo->_xBin));
	pHistoNode->QueryFloatAttribute("xMin", &(histoInfo->_xMin));
	pHistoNode->QueryFloatAttribute("xMax", &(histoInfo->_xMax));
	
	if ( ( histoInfo->_type == string("H2D") ) || 
	     ( histoInfo->_type == string("H3D") ) ) {
	  pHistoNode->QueryIntAttribute("yBin", &(histoInfo->_yBin));
	  pHistoNode->QueryFloatAttribute("yMin", &(histoInfo->_yMin));
	  pHistoNode->QueryFloatAttribute("yMax", &(histoInfo->_yMax));
	}
	
	if ( histoInfo->_type == string("H3D") ) {
	  pHistoNode->QueryIntAttribute("zBin", &(histoInfo->_zBin));
	  pHistoNode->QueryFloatAttribute("zMin", &(histoInfo->_zMin));
	  pHistoNode->QueryFloatAttribute("zMax", &(histoInfo->_zMax));
	}
	
	if ( histoInfo->_type == string("P1D")) {
	  pHistoNode->QueryFloatAttribute("yMin", &(histoInfo->_yMin));
	  pHistoNode->QueryFloatAttribute("yMax", &(histoInfo->_yMax));
	}
	
	if ( histoInfo->_type == string("P2D")) {
	  pHistoNode->QueryIntAttribute("yBin", &(histoInfo->_yBin));
	  pHistoNode->QueryFloatAttribute("yMin", &(histoInfo->_yMin));
	  pHistoNode->QueryFloatAttribute("yMax", &(histoInfo->_yMax));
	  pHistoNode->QueryFloatAttribute("zMin", &(histoInfo->_zMin));
	  pHistoNode->QueryFloatAttribute("zMax", &(histoInfo->_zMax));
	}
      }
      _histoInfoMap.insert( make_pair( histoInfo->_name, histoInfo ) );
      pHistoNode = pHistoNode->NextSiblingElement();
    }

  } else {
    
    stringstream ss;
    ss << "EUTelHistogramManager::init error in file " << _histoInfoFileName 
       << ", row: " << doc->ErrorRow() << ", col: " << doc->ErrorCol() << " : " 
       << doc->ErrorDesc();
    delete doc;
    throw ParseException( ss.str() );
  }

  delete doc;

  return true;
}

EUTelHistogramInfo * EUTelHistogramManager::getHistogramInfo(std::string histoName) const {

  std::map< std::string , EUTelHistogramInfo *>::const_iterator iter = _histoInfoMap.find(histoName);
  if ( iter == _histoInfoMap.end() ) return 0x0;
  return iter->second;
  

}


// #endif 
