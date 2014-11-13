#include "EUTelHistogram.h"



namespace eutelescope {

	EUTelHistogram::EUTelHistogram(std::string name, std::string histoName){
    streamlog_out( DEBUG5) << " Running book constructor " << std::endl;
		setHistogramInfoName(histoName);
		book(name);
};

	//This will book any histograms given within the xml brackets (name)
	int EUTelHistogram::book(std::string name){
/*
		//Load the xml file and then check that the brackets histogramManagere exist////////////////////////////////BEGIN
		TiXmlDocument * doc = new TiXmlDocument; //open document

  	if ( doc->LoadFile( _histoInfoFileName ) ){ //load and if was ok begin

  	  TiXmlHandle    hDoc(doc);
  	  TiXmlHandle    hRoot(0);
  	  TiXmlElement * pElem;
    
  	  pElem = hDoc.FirstChild("HistogramManager").Element();
  	  if ( !pElem ) {
  	    delete doc;
  	    throw ParseException( string( "EUTelHistogramManager::init: no root tag <HistogramManager> ... </HistogramManager> found in ") + _histoInfoFileName );
    	} else {
       	streamlog_out( DEBUG5) << " Found HistogramManager Node! " << std::endl;
      	hRoot = TiXmlHandle(pElem);
    

			///////////////////////////////////////////////////Check if brackets histos is inside xml file ///////////////////BEGIN
  		TiXmlNode * nHistosBlock = hRoot.FirstChild( "histos" ).ToNode();
			if ( !nHistosBlock ){
  		delete doc;
  		throw ParseException( string( "EUTelHistogramManager::init: no <histos> ... </histos> block found in ")  + _histoInfoFileName);
  		}else{
       	streamlog_out( DEBUG5) << " Found Histos Node! " << std::endl;


				//////////////////////////Look for the tag HistoProcessorName////////////////////////////////////BEGIN
				TiXmlElement * pHistoNode = NULL;
    		pHistoNode = hRoot.FirstChild( "histos" ).FirstChild("ProcessName").Element(); //This is where we open the actual histo
    		if( pHistoNode != NULL ){
					if( pHistoNode->Attribute("Name") = name.c_str()){
					streamlog_out( DEBUG5) << " Name of processor found: "<<  name.c_str() <<std::endl;
					TiXmlHandle process = TiXmlHandle(pHistoNode);
    			TiXmlElement * detector = process.FirstChild( "Directory" ).Element(); 
					//for( detector; detector; detector=detector->NextSiblingElement()){ //Loop thorugh all directories of the same level
						TiXmlHandle detectorHandle = TiXmlHandle(detector);
						std::string detectorString = detector->Attribute("List");
						streamlog_out( DEBUG5) << " The detector String: "<< detectorString <<std::endl;
						std::vector<std::string> ListDetector;
						stringSplitting(detectorString, ListDetector);
						streamlog_out( DEBUG5) << " The List of detectors size: "<< ListDetector.size() <<" The first addition "<< ListDetector[0] <<" The last "<<ListDetector[ListDetector.size() -1]<< std::endl;			
						for(int i = 0; i<ListDetector.size() -1; i++){ //Create 
						TiXmlElement* histogram1D = detectorHandle.FirstChild( "H1D" ).ToElement();
						TiXmlHandle histogram1DHandle = TiXmlHandle(histogram1D);
						TiXmlElement* HistoName = histogram1DHandle.FirstChild( "HistoName" ).ToElement();
						for( HistoName; HistoName; HistoName=HistoName->NextSiblingElement() ){
							std::string HistoNameString = HistoName->Attribute("name");
							std::string HistoTitleString = HistoName->Attribute("title");
							int HistoBinX = atoi(HistoName->Attribute("xBin"));
							float HistoxMin = atof(HistoName->Attribute("xMin"));
							float HistoxMax = atof(HistoName->Attribute("xMax"));
							streamlog_out( DEBUG5) << " This is the name of this histogram: "<< HistoNameString << " title:  " <<HistoTitleString << " bin, min, max: " <<HistoBinX <<" , " <<HistoxMin << " , " <<HistoxMax <<std::endl;
							TiXmlHandle HistoNameHandle = TiXmlHandle(HistoName);
							TiXmlElement* HistoList = HistoNameHandle.FirstChild( "HistoNameList" ).ToElement();
							std::string HistoListNameString = HistoList->Attribute("List");
							std::string HistoListTitleString = HistoList->Attribute("titleList");
							std::vector<std::string> listName;
							std::vector<std::string> titleList;
							stringSplitting(HistoListNameString, listName);
							stringSplitting(HistoListTitleString, titleList);
							if(listName.size() == titleList.size()){
								streamlog_out( DEBUG5) << " List name size "<< listName.size() <<" The first addition "<< listName[0] <<" The last "<<listName[listName.size() -1]<< std::endl;
								streamlog_out( DEBUG5) << " titleList name size "<< titleList.size() <<" The first addition "<< titleList[0] <<" The last "<<titleList[titleList.size() -1]<< std::endl;

								for(int j = 0 ; j <listName.size(); j++){ 
								//	std::string total_string_ID = ListDetector[i] + HistoNameString + 

								}//end of list name 
							
				 
						
								}else{
								streamlog_out( DEBUG5) << " Title list number and identification do not match: "<< listName.size()<< " and " << titleList.size() <<std::endl; 
								} //end of are list title and list name same size
						}//END of histogram name loop
					} //END of detector size loop
					
					
				}else{
  				throw ParseException( string( "Did not find bracket with name  ")  + name);
				} 
			}//END OF HISTOS IF STATEMENT 
 		}//END OF HISTOMANAGER IF STATEMENT
	
	}else{		
		streamlog_out( DEBUG5) << " No histogram info found at: "<< _histoInfoFileName<< std::endl;	
	}	//END OF DOC IS OPEN IF STATEMENT  */
}//END OF FUNCTION BOOK HISTOGRAMS


void EUTelHistogram::stringSplitting(std::string input, std::vector<std::string> &output){
	std::stringstream stringStream(input);
	std::string line;
	while(std::getline(stringStream, line)) 
	{
	    std::size_t prev = 0, pos;
	    while ((pos = line.find_first_of(";", prev)) != std::string::npos)
	    {
	        if (pos > prev)
	            output.push_back(line.substr(prev, pos-prev));
	        prev = pos+1;
	    }
	    if (prev < line.length())
	        output.push_back(line.substr(prev, std::string::npos));	
	}
}

}

void EUTelHistogram::pushStringTogether(std::string input1, std::string input2, std::string & output){

}
