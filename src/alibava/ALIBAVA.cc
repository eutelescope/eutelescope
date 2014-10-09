/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


#include "ALIBAVA.h"

// system includes
#include <algorithm>
#include <string>

using namespace alibava;
using namespace std;

// Alibava Readout system parameters
// If you do change these, you might need to change things in AlibavaConverter processor
// const int      ALIBAVA::NOOFCHIPS            = 2;
// const int      ALIBAVA::NOOFCHANNELS         = 128;

// Global Alibava Processor parameters
const char *   ALIBAVA::CHANNELSTOBEUSED    = "ChannelsToBeUsed";
const char *   ALIBAVA::SKIPMASKEDEVENTS    = "SkipMaskedEvents";


// General Parameters
const char *   ALIBAVA::GEOID               = "GeoID";
const char *   ALIBAVA::NOTSET              = "notset";

// Optional Parameters
const char *   ALIBAVA::TILTANGLE           = "TiltAngle";
const char *   ALIBAVA::SENSORTEMPERATURE   = "SensorTemperature";

// Run Header Parameters
const char *   ALIBAVA::HEADER              = "Header";
const char *   ALIBAVA::HEADERVERSION       = "HeaderVersion";
const char *   ALIBAVA::NOOFEVENT           = "NoOfEvent";
const char *   ALIBAVA::DATATYPE            = "DataType";
const char *   ALIBAVA::DATETIME            = "DateTime";
const char *   ALIBAVA::HEADERPEDESTAL      = "HeaderPedestal";
const char *   ALIBAVA::HEADERNOISE         = "HeaderNoise";
const char *   ALIBAVA::APPLIEDPROCESSOR         = "AppliedProcessor";
const char *   ALIBAVA::ISTHERESELECTION         = "IsThereChipSelection";
const char *   ALIBAVA::SELECTEDCHIPNUM         = "SelectedChipNumbers";

// Event Parameters
const char *   ALIBAVA::EVENTTYPE           = "EventType";
const char *   ALIBAVA::EVENTVALUE          = "EventValue";
const char *   ALIBAVA::EVENTSIZE           = "EventSize";
const char *   ALIBAVA::EVENTCODE           = "EventCode";
const char *   ALIBAVA::EVENTTIME           = "EventTime";
const char *   ALIBAVA::EVENTTEMP           = "EventTemp";
const char *   ALIBAVA::CALCHARGE           = "CalCharge";
const char *   ALIBAVA::CALDELAY            = "CalDelay";
const char *   ALIBAVA::EVENTMASK           = "EventMask";

// CellIDEncoding for AlibavaData
const char *   ALIBAVA::ALIBAVADATA_ENCODE           = "ChipNum:1,";
const char *   ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM         = "ChipNum";

// CellIDEncoding for AlibavaClusters
const char *   ALIBAVA::ALIBAVACLUSTER_ENCODE  = "ClusterID:6,IsSensitiveAxisX:1,IsSignalNegative:1,ChipNum:1,ClusterSize:8,Seed:16";
const char *   ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERID = "ClusterID";
const char *   ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSENSITIVEAXISX = "IsSensitiveAxisX";
const char *   ALIBAVA::ALIBAVACLUSTER_ENCODE_ISSIGNALNEGATIVE = "IsSignalNegative";
const char *   ALIBAVA::ALIBAVACLUSTER_ENCODE_CHIPNUM = "ChipNum";
const char *   ALIBAVA::ALIBAVACLUSTER_ENCODE_SEED = "Seed";
const char *   ALIBAVA::ALIBAVACLUSTER_ENCODE_CLUSTERSIZE = "ClusterSize";


std::string trim_right(const std::string &s)
{
	std::string b=" \t\n";
	std::string str = s;
	return str.erase(str.find_last_not_of(b) +1);
}

std::string trim_left(const std::string &s)
{
	std::string b=" \t\n";
	std::string str = s;
	return str.erase( 0, str.find_first_not_of(b) );
}

std::string trim_str(const std::string &s)
{
	std::string str = s;
	return trim_left(trim_right(str) );
}

string getSubStringUpToChar(string search_string, const char* search_char, size_t start_position){
	size_t found_position = search_string.find(search_char,start_position);
	if (found_position==string::npos) found_position = search_string.size();
	
	size_t nchar = found_position - start_position;
	string s = search_string.substr(start_position,nchar);
	return s;
	
}

