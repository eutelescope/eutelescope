//////////////////////
//
//	Analysis of Testbeam data
//
//	This loops over some ntuples created with EUTelescope and then does magic
//
//	Thomas Eichhorn 2015, updated 2016
//	(thomas.eichhorn@desy.de)
//
//////////////////////

// start with:		root +x -q -l -b 'main.cc("runselection")'
// compile with:	g++ -I `root-config --incdir` -o asdf main.cc `root-config --libs` -std=c++0x -pedantic -Wextra


/*
 *
 * Rescue Pig to the rescue!

                                 _
    _._ _..._ .-',     _.._(`))
   '-. `     '  /-._.-'    ',/
      )         \            '.
     / _    _    |             \
    |  a    a    /              |
    \   .-.                     ;  
     '-('' ).-'       ,'       ;
        '-;           |      .'
           \           \    /
           | 7  .__  _.-\   \
           | |  |  ``/  /`  /
          /,_|  |   /,_/   /
             /,_/      '`-'


*/


//C++ headers
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <utility>
#include <map>
#include <sstream>
#include <memory>

//Root headers
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTree.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TStopwatch.h"
#include "TLatex.h"

#define PI 3.14159265

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	CUTS
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// the number of runs we read from the runlist -> 230 max for all runs

// 91 for all 0deg
// 11 for good n
// 33 for 0deg p
// 25 for 0deg y
// 33 for 0deg n (all)

// 54 for 0,25 p
// 41 for 0,25 y
# define _runcount 231

// the evaluation step for the x position in mm, must be > 0.4
double _xevalstep = 2.5;

bool _doxposcheck = false;

// the conversion from ADCs to electrons
double _electronconversion = 183.0;

// the error on this
double _electronconversionerror = 1.0;

// the maximum number of tracks we allow in an event
int _cut_maxtracksperevent = 100;

// only singletrack events?
bool _cut_onlysingletrackevents = false;

// do we require the track signal to be highest in the middle 3 strips?
bool _cut_highchannelmiddle = true;

// the +- distance for a track to be considered "on" a strip
// 0.25 * 0.08 pitch -> 20 um
float _cut_onstrip = 0.25;

// the debug level of this software
int _debug = -5;

// the maximum number of telescope tracks in a run
int _maxtotaltracks = 99999999;
//int _maxtotaltracks = 50;

// the telescope resolution in um - this is used to plot a line in the residual vs voltage plot, this is not really a cut
float _cut_telescoperesolution = 5.0;

// the noise range to plot in ADCs
float _cut_minnoise = 0.0;
float _cut_maxnoise = 50.0;

// the residual range to plot in mm
float _cut_minXresidual = 0.0;
float _cut_maxXresidual = 0.5;
float _cut_minYresidual = 0.01;
float _cut_maxYresidual = 0.04;

// the signal range to plot in ADCs
float _cut_minsignal = 0.0;
float _cut_maxsignal = 200.0;

// do we want an extra comparison between the eta integral of matched and unmatched etas?
bool _cut_drawetaintegralcomparison = false;

// the eta value above (and 1- below) we consider as chargesharing
float _cut_eta_chargeshare = 0.2;

// the maximum channel noise value to calculate rghs in ADCs
// if lots of channels in a run are over this noise then rgh calculation makes no sense, since the sensor is too noisy
float _cut_maxrghnoise = 75.0;

// the maximum rgh percent value to plot
float _cut_maxrghplot = 10.0;

// the maximum cluster count per track we want to plot
float _cut_maxclustercount = 1.0;

// the cut times noise to define charge sharing
float _cut_chargesharing = 1.5;

// apply a temperature correction?
bool _cut_applytempcorrection = true;

// do a check of the TDC to find highest area?
bool _dotdccheck = false;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Global settings
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// the sensor features
// pitch in mm
float _pitchx = 0.018402778;
float _pitchy = 0.08;

// absolute error on the voltage in V:
float _voltageerror = 5.0;

// absolute error on the clustercount per run:
float _clustercounterror = 0.0;

// the runlist if nothing is selected
string _runfile = "../runlist.csv";

// the suffix of our root files
//string _filesuffix = "-alibava-tracking-3.root";
string _filesuffix = "-alibava-tracking-gbl.root";

// the suffix of the root cluster files
string _clusterfilesuffix = "-alibava-clustering-2.root";

// the vectors to store runlist information in
// _total gives the number of distinct entries in the vectors for plotting etc.
std::vector<int> _counter;
std::vector<int> _runnumber;
std::vector<int> _pedestal;
std::vector<int> _thickness;
std::vector<int> _thickness_list;
int _thickness_total = 0;
std::vector<int> _dutrotation;
std::vector<int> _dutrotation_list;
int _dutrotation_total = 0;
std::vector<int> _biasvoltage;
std::vector<int> _biasvoltage_list;
int _biasvoltage_total = 0;
std::vector<int> _temperature;
std::vector<int> _temperature_list;
int _temperature_total = 0;
std::vector<double> _sensorcurrent;
std::vector<double> _sensorcurrent_list;
int _sensorcurrent_total = 0;
std::vector<int> _polarity;
std::vector<int> _polarity_list;
int _polarity_total = 0;
std::vector<string> _sensorname;
std::vector<string> _sensorname_list;
int _sensorname_total = 0;
std::vector<char> _type;
std::vector<char> _type_list;
int _type_total = 0;
std::vector<float> _irradfluence;
std::vector<float> _irradfluence_list;
int _irradfluence_total = 0;
int _channels[_runcount][128];
float _noise[_runcount][128];
std::vector<int> _goodtracks;
std::vector<float> _sensorminx;
std::vector<float> _sensormaxx;
std::vector<float> _sensormintdc;
std::vector<float> _sensormaxtdc;
std::vector<float> _sensoralignx;
std::vector<float> _sensoraligny;
std::vector<float> _sensoralignz;
std::vector<float> _sensoraligna;
std::vector<float> _sensoralignb;
std::vector<float> _sensoralignc;
std::vector<int> _uselevel;
std::vector<string> _runcomments;

// the value which is written to the tuple by EUTelescope if there is no entry
int _missingvalue = -999;

// map for the root objects
std::map < string , TObject * > _rootObjectMap;

// the global run time
TStopwatch _totaltime;

// how many we have read from file
int _actual_runcount = 0;

int _bookinglimit = 200;

string _thefilename = "failnoselection.root";

TFile* _outputFile;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations

void readrunlist();

void bookhistos();

void openfile();

TF1* lanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit);

TF1* gausLanGausFitFixGausNoise(TH1D* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// let's go!

int main(int argc, char** argv)
{

	string  typetorun = "fail";
	stringstream astream;
	if (argc>1)
	{
		astream << argv[1];
		typetorun = astream.str();
		cout << " " << endl;
		cout << " Your input: " << typetorun << endl;
		cout << " " << endl;
	} else {
		cout << " " << endl;
		cout << "Please select a type to run!" << endl;
		cout << " " << endl;
		return 0;
	}

	// root pre stuff
	//TApplication* test = new TApplication();

	gROOT->SetBatch();
	gStyle->SetLabelSize(0.035, "x");
	gStyle->SetLabelSize(0.035, "y");
	gStyle->SetTitleSize(0.05, "x");
	gStyle->SetTitleSize(0.05, "y");
	gStyle->SetTitleOffset(0.95, "x");
	gStyle->SetTitleOffset(0.95, "y");
	gStyle->SetOptFit(1111);

	_totaltime.Start();

	cout << "##############################################" << endl;
	cout << " " << endl;
	cout << "Hello!" << endl;
	cout << " " << endl;

	// get the bulkselection from user:
	if (typetorun == "n")
	{
		cout << "You selected all n-types!" << endl;
		_runfile = "../runlists/n.csv";
		_thefilename = "epi_n.root";

	} else if (typetorun == "p") {

		cout << "You selected all p-stop-types!" << endl;
		_runfile = "../runlists/p.csv";
		_thefilename = "epi_p.root";

	} else if (typetorun == "y") {

		cout << "You selected all p-spray-types!" << endl;
		_runfile = "../runlists/y.csv";
		_thefilename = "epi_y.root";

	} else if (typetorun == "all") {

		cout << "You selected all runs!" << endl;
		_runfile = "../runlist.csv";
		_thefilename = "epi_all.root";

	} else if (typetorun == "t") {

		cout << "You selected test!" << endl;
		_runfile = "../testlist.csv";
		_thefilename = "test.root";

	} else if (typetorun == "d15") {

		cout << "You selected all d15 runs!" << endl;
		_runfile = "../runlists/d15.csv";
		_thefilename = "epi_d15.root";

	} else if (typetorun == "d16") {

		cout << "You selected all d16 runs!" << endl;
		_runfile = "../runlists/d16.csv";
		_thefilename = "epi_d16.root";

	} else if (typetorun == "a") {

		cout << "You selected angular scan runs!" << endl;
		_runfile = "../runlists/a.csv";
		_thefilename = "epi_angular.root";

	} else {
		cout << "You selected nothing!" << endl;
	}

	// the output
	_outputFile = new TFile(_thefilename.c_str(), "RECREATE");

	// read the runlist into the vectors
	readrunlist();

	// book ALL histograms
	bookhistos();

	// open the root files and fill the histograms
	openfile();

	// done
	_totaltime.Stop();

	cout << " " << endl;
	cout << "##############################################" << endl;
	cout << " " << endl;
	cout << "Time elapsed: " << _totaltime.RealTime() << " s!" << endl;
	cout << " " << endl;
	cout << "Goodbye!" << endl;
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// this reads the runlist to get filenames of tuples
void readrunlist()
{

	cout << "##############################################" << endl;
	cout << " " << endl;
	cout << "Reading " << _runfile << " !" << endl;
	cout << " " << endl;

	// the strings we parse into
	string fail;
	string sensor;
	string fluence;
	string rotation;
	string voltage;
	string temp;
	string current;
	string minx;
	string maxx;
	string tdcmin;
	string tdcmax;
	string bonds;
	string pede;

	// the ints
	int ped;
	int run;
	int pol;

	// the floats
	float flu;
	float minxval;
	float maxxval;
	float mintdcval;
	float maxtdcval;

	// the line to read
	string line;

	// the stream
	ifstream fileRead;

	// open
	fileRead.open(_runfile.c_str());

	// count on the found sensors
	int i = 0;

	// flag for finding information
	bool foundinfo = false;
	bool foundpede = false;

	// loop over the lines in the file
	while (std::getline(fileRead, line))
	{
		string startpart = line.substr(0,2);
		if (startpart == "#!")
		{

			// we found something
			foundinfo = true;

			// the input string of this line
			string input;
			std::istringstream iss(line);
			input = iss.str();

			// the delimiter of the comment (space)
			string delimiter = " ";

			// the position we found of the delimiter
			size_t pos = 0;

			// the iterator for counting elements
			std::vector<int>::iterator it;

			// first the #!
			pos = input.find(delimiter);
			fail = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());

			// now the sensor
			// push back into _polarity, sensor name, _thickness and _type
			pos = input.find(delimiter);
			sensor = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			int thick = 0;
			thick = atoi(sensor.c_str());
			_thickness.push_back(thick);
			// for now assume simulated data ('S') is also of positive polarity
			if (sensor.at(3) == 'N' || sensor.at(3) == 'S')
			{
				pol = 1;
			} else {
				pol = -1;
			}

			it = find (_polarity.begin(), _polarity.end(), pol);
			if (it != _polarity.end())
			{
				if (_debug <=1)
				{
					cout << "Polarity " << pol << " already in list!" << endl;
				}
			} else {
				_polarity_total++;
				_polarity_list.push_back(pol);
				if (_debug <=2)
				{
					cout << "Polarity " << pol << " not in list!" << endl;
					cout << " Total polarities now " << _polarity_total << endl;
				}
			}

			// string iterator
			std::vector<string>::iterator its;

			its = find (_sensorname.begin(), _sensorname.end(), sensor);
			if (its != _sensorname.end())
			{
				if (_debug <=1)
				{
					cout << "Sensor " << sensor << " already in list!" << endl;
				}
			} else {
				_sensorname_total++;
				_sensorname_list.push_back(sensor);
				if (_debug <=2)
				{
					cout << "Sensor " << sensor << " not in list!" << endl;
					cout << " Total sensors now " << _sensorname_total << endl;
				}
			}

			// char iterator
			std::vector<char>::iterator itc;

			itc = find (_type.begin(), _type.end(), sensor.at(3));
			if (itc != _type.end())
			{
				if (_debug <=1)
				{
					cout << "Type " << sensor.at(3) << " already in list!" << endl;
				}
			} else {
				_type_total++;
				_type_list.push_back(sensor.at(3));
				if (_debug <=2)
				{
					cout << "Type " << sensor.at(3) << " not in list!" <<endl;
					cout << " Total types now " << _type_total << endl;
				}
			}

			_polarity.push_back(pol);
			_sensorname.push_back(sensor);
			_type.push_back(sensor.at(3));

			// now the fluence
			pos = input.find(delimiter);
			fluence = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			if (fluence == "unirr,")
			{
				flu = 0.0;
			} else {
			    flu = atof(fluence.c_str());
			}

			// float iterator
			std::vector<float>::iterator itf;

			itf = find (_irradfluence.begin(), _irradfluence.end(), flu);
			if (itf != _irradfluence.end())
			{
				if (_debug <=1)
				{
					cout << "Irradiation " << flu << " already in list!" << endl;
				}
			} else {
				_irradfluence_total++;
				_irradfluence_list.push_back(flu);
				if (_debug <=2)
				{
					cout << "Irradiation " << flu << " not in list!" << endl;
					cout << " Total fluences now " << _irradfluence_total << endl;
				}
			}

			_irradfluence.push_back(flu);

			// now the rotation
			pos = input.find(delimiter);
			rotation = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			int rot = 0;
			rot = atoi(rotation.c_str());

			it = find (_dutrotation.begin(), _dutrotation.end(), rot);
			if (it != _dutrotation.end())
			{
				if (_debug <=1)
				{
					cout << "Rotation " << rot << "deg already in list!" << endl;
				}
			} else {
				_dutrotation_total++;
				_dutrotation_list.push_back(rot);
				if (_debug <=2)
				{
					cout << "Rotation " << rot << "deg not in list!" << endl;
					cout << " Total rotations now " << _dutrotation_total << endl;
				}
			}

			_dutrotation.push_back(rot);

			// now the voltage
			pos = input.find(delimiter);
			voltage = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			int vol = 0;
			vol = atoi(voltage.c_str());

			it = find (_biasvoltage.begin(), _biasvoltage.end(), vol);
			if (it != _biasvoltage.end())
			{
				if (_debug <=1)
				{
					cout << "Voltage " << vol << "V already in list!" << endl;
				}
			} else {
				_biasvoltage_total++;
				_biasvoltage_list.push_back(vol);
				if (_debug <=2)
				{
					cout << "Voltage " << vol << "V not in list!" << endl;
					cout << " Total voltages now " << _biasvoltage_total << endl;
				}
			}

			_biasvoltage.push_back(vol);

			// now the temp
			pos = input.find(delimiter);
			temp = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			int tem = 0;
			tem = atoi(temp.c_str());

			it = find (_temperature.begin(), _temperature.end(), tem);
			if (it != _temperature.end())
			{
				if (_debug <=1)
				{
					cout << "Temperature " << tem << "°C already in list!" << endl;
				}
			} else {
				_temperature_total++;
				_temperature_list.push_back(tem);
				if (_debug <=2)
				{
					cout << "Temperature " << tem << "°C not in list!" << endl;
					cout << " Total temperatures now " << _temperature_total << endl;
				}
			}

			_temperature.push_back(tem);

			// now the current
			pos = input.find(delimiter);
			current = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			float cur = 0;
			cur = atof(current.c_str());

			// double iterator
			std::vector<double>::iterator itd;

			itd = find (_sensorcurrent.begin(), _sensorcurrent.end(), cur);
			if (itd != _sensorcurrent.end())
			{
				if (_debug <=1)
				{
					cout << "Sensor Current " << cur << "mA already in list!" << endl;
				}
			} else {
				_sensorcurrent_total++;
				_sensorcurrent_list.push_back(cur);
				if (_debug <=2)
				{
					cout << "Sensor Current " << cur << "mA not in list!" << endl;
					cout << " Total sensor currents now " << _sensorcurrent_total << endl;
					cout << " " << endl;
				}
			}

			_sensorcurrent.push_back(cur);

			// now the min/max x

			pos = input.find(delimiter);
			minx = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			minxval = 0;
			minxval = atof(minx.c_str());
			_sensorminx.push_back(minxval);

			pos = input.find(delimiter);
			maxx = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			maxxval = 0;
			maxxval = atof(maxx.c_str());
			_sensormaxx.push_back(maxxval);

			// the min/max tdc
			pos = input.find(delimiter);
			tdcmin = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			mintdcval = 0;
			mintdcval = atof(tdcmin.c_str());
			_sensormintdc.push_back(mintdcval);

			pos = input.find(delimiter);
			tdcmax = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			maxtdcval = 0;
			maxtdcval = atof(tdcmax.c_str());
			_sensormaxtdc.push_back(maxtdcval);

			// now the aligment positions
			string tempstring;
			float tempfloat;
			pos = input.find(delimiter);
			tempstring = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			tempfloat = 0;
			tempfloat = atof(tempstring.c_str());
			_sensoralignx.push_back(tempfloat);

			pos = input.find(delimiter);
			tempstring = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			tempfloat = 0;
			tempfloat = atof(tempstring.c_str());
			_sensoraligny.push_back(tempfloat);

			pos = input.find(delimiter);
			tempstring = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			tempfloat = 0;
			tempfloat = atof(tempstring.c_str());
			_sensoralignz.push_back(tempfloat);

			pos = input.find(delimiter);
			tempstring = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			tempfloat = 0;
			tempfloat = atof(tempstring.c_str());
			_sensoraligna.push_back(tempfloat);

			pos = input.find(delimiter);
			tempstring = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			tempfloat = 0;
			tempfloat = atof(tempstring.c_str());
			_sensoralignb.push_back(tempfloat);

			pos = input.find(delimiter);
			tempstring = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			tempfloat = 0;
			tempfloat = atof(tempstring.c_str());
			_sensoralignc.push_back(tempfloat);

			pos = input.find(delimiter);
			tempstring = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());
			int tempint = 0;
			tempint = atoi(tempstring.c_str());
			_uselevel.push_back(tempint);

			// the remaining input string is the comment
			_runcomments.push_back(input);

			// output
			if (_debug <= 2)
			{
				cout << "Sensor " << i << ": " << sensor << " " << fluence << " " << rotation << " " << voltage << " " << temp << " " << current << endl;
				cout << " Range: X: " << minxval << " " << maxxval << " T: " << mintdcval << " " << maxtdcval << endl;
				cout << " Uselevel: " << tempint << ", comments: " << input << endl;
				cout << " " << endl;
			}

		} else {

			// if the line is not a special comment, we look if the bools are true.
			// first run, then ped, so we don't find things twice!

			if (foundpede == true)
			{

				// read run numbers
				std::istringstream iss(line);
				while (iss >> run)
				{
					if (_debug <= 3)
					{
						cout << "Data runnumber " << i << ": " << run << endl;
						cout << " " << endl;
					}
				}

				// set back to false
				foundpede = false;

				_counter.push_back(i);
				_runnumber.push_back(run);

				i++;
			}

			if (foundinfo == true)
			{
				string input;
				std::istringstream iss(line);
				input = iss.str();
				string delimiter = ",";
				size_t pos = input.find(delimiter);
				pede = input.substr(0, pos);
				input.erase(0, pos + delimiter.length());
				ped = atoi(pede.c_str());
				_pedestal.push_back(ped);

				if (_debug <= 1)
				{
					cout << "Pedestal runnumber " << i << ": " << ped << endl;
					cout << " " << endl;
				}

				delimiter = "'";
				pos = input.find(delimiter);
				fail = input.substr(0, pos);
				input.erase(0, pos + delimiter.length());

				pos = input.find(delimiter);
				bonds = input.substr(0, pos);
				input.erase(0, pos + delimiter.length());

				// set all channels to off
				for (int ii=0;ii<128;ii++)
				{
					_channels[i][ii] = 0;
				}

				int start, end;
				pos=0;
				while (pos<bonds.size()-1)
				{

					delimiter = ":";
					pos = bonds.find(delimiter);
					fail = bonds.substr(0, pos);
					bonds.erase(0, pos + delimiter.length());
					start = atoi(bonds.c_str());

					delimiter = "-";
					pos = bonds.find(delimiter);
					fail = bonds.substr(0, pos);
					bonds.erase(0, pos + delimiter.length());
					end = atoi(bonds.c_str());

					for (int ii=start;ii<=end;ii++)
					{
						_channels[i][ii] = 1;
					}

				}

				if (_debug <= 2)
				{
					cout << "Good Channels:" << endl;
					cout << "   0 - 31: ";
					for (int ii=0;ii<32;ii++)
					{
						cout << _channels[i][ii];
					}
					cout << endl << "  32 - 63: ";
					for (int ii=32;ii<64;ii++)
					{
						cout << _channels[i][ii];
					}
					cout << endl << "  64 - 95: ";
					for (int ii=64;ii<96;ii++)
					{
						cout << _channels[i][ii];
					}
					cout << endl << " 96 - 127: ";
					for (int ii=96;ii<128;ii++)
					{
						cout << _channels[i][ii];
					}
					cout << endl;
					cout << " " << endl;
				}

				// set back to false
				foundinfo = false;
				foundpede = true;

				_pedestal.push_back(ped);
			}

		}

		// limit the reading to save time...
		if (i == _runcount)
		{
			break;
		}
	}

	cout << "Done reading runlist after " << i << " found sensors!" << endl;
	cout << " " << endl;
	_actual_runcount = i;
	// protect against small runs, as voltagegraphs is encoded...
	if (_actual_runcount < 15)
	{
		_bookinglimit = 50;
	} else {
		_bookinglimit = _actual_runcount*5;
	}
	fileRead.close();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// this reads the nth file in the main _runnumber vector
void openfile()
{

	cout << "##############################################" << endl;
	cout << " " << endl;
	cout << "Opening root files and filling histograms!" << endl;
	cout << " " << endl;

	// the names and titles for the histograms
	char name[100];

	// the string for the name
	string histoname;
	stringstream sstream;

	// the fractions used later for mod pitch
	double fractpartx, fractparty, intpartx, intparty;

	// this will count the entries in a graph
	int tempcount = 0;

	// the vector of all graphs if we do vs voltage
	std::vector<int> voltagegraphs_total;
	// the individual ones
	std::vector<int> voltagegraphs;
	// the point count in each graph
	int voltagepoint[500];
	for (int i = 0; i< 500; i++)
	{
		voltagepoint[i] = 0;
	}

	// these global histos will get the alignments of each run
	histoname = "xshift";
	sstream << histoname;
	TH1D* xshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
	sstream.str(string());

	histoname = "yshift";
	sstream << histoname;
	TH1D* yshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
	sstream.str(string());

	histoname = "zshift";
	sstream << histoname;
	TH1D* zshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
	sstream.str(string());

	histoname = "ashift";
	sstream << histoname;
	TH1D* ashift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
	sstream.str(string());

	histoname = "bshift";
	sstream << histoname;
	TH1D* bshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
	sstream.str(string());

	histoname = "cshift";
	sstream << histoname;
	TH1D* cshift = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
	sstream.str(string());

	// now open each file!
	for ( int irun = 0 ; irun < _actual_runcount ; irun++ )
	{

		// set the noise map in this run to zero
		for (int ii=0;ii<128;ii++)
		{
			_noise[irun][ii] = 0.0;
		}

		// get the actual filenumber from the vector
		int filenumbertoopen = _runnumber.at(irun);

		// is this an x-sensitive run?
		bool xsensitive = false;
		if (filenumbertoopen == 813)
		{
			xsensitive = true;
		}

		// output
		if (_debug <= 3)
		{
			cout << " " << endl;
			cout << "**********************************************" << endl;
			cout << "Reading run " << irun << " of " << _actual_runcount - 1 << " , nr. " << _runnumber.at(irun) << endl;
		}

		// base string
		string filename("../output/histograms/00");
		string cfilename("../output/histograms/00");

		// leading zero
		if (filenumbertoopen <= 99)
		{
			filename+="0";
			cfilename+="0";
		}

		// leading zero again
		if (filenumbertoopen <= 999)
		{
			filename+="0";
			cfilename+="0";
		}

		// move the number to a string
		stringstream ss;
		ss << filenumbertoopen;
		string numbertoopen = ss.str();

		// construct the filename
		filename+=numbertoopen;
		filename+=_filesuffix;
		cfilename+=numbertoopen;
		cfilename+=_clusterfilesuffix;

		if (_debug <= 2)
		{
			cout << "Opening file: " << filename << endl;
		}

		// enter the file...
		TFile *f1 = TFile::Open(filename.c_str());
		f1->cd();

		TFile *f2 = TFile::Open(cfilename.c_str());
		f2->cd();

		if (_debug <= 2)
		{
			cout << "File open!" << endl;
		}

		// clone 3d hitmap
		TH3D* DUThitmap3D = (TH3D*) (f1->Get("Ntuple/DUTHitmap")); // valgrind complains
		DUThitmap3D->SetDirectory(0);

		// fail FIXME
		//get the final alignment from this
		/*
		float min3dx = 0.0;
		float max3dx = 0.0;
		float min3dy = 0.0;
		float max3dy = 0.0;
		float min3dz = 0.0;
		float max3dz = 0.0;
		*/
		/*
		// the hitmap is filled xzy!
		min3dx = DUThitmap3D->FindFirstBinAbove(0,1);
		max3dx = DUThitmap3D->FindLastBinAbove(0,1);
		min3dy = DUThitmap3D->FindFirstBinAbove(0,3);
		max3dy = DUThitmap3D->FindLastBinAbove(0,3);
		min3dz = DUThitmap3D->FindFirstBinAbove(0,2);
		max3dz = DUThitmap3D->FindLastBinAbove(0,2);
		
		cout << " " << endl;
		cout << "x: " << max3dx << " to " << min3dx << endl;
		cout << "yshift: " << max3dy << " to " << min3dy << endl;
		cout << "zshift: " << max3dz << " to " << min3dz << endl;

		cout << " " << endl;
		*/

		// since alpha is turned, subtract the nominal position to get the "alignment"
		xshift->Fill(_sensoralignx.at(irun));
		yshift->Fill(_sensoraligny.at(irun));
		zshift->Fill(_sensoralignz.at(irun));
		ashift->Fill(_sensoraligna.at(irun) - _dutrotation.at(irun));
		bshift->Fill(_sensoralignb.at(irun));
		cshift->Fill(_sensoralignc.at(irun));

		// get the cluster info
		TH1D* clustersizehisto = (TH1D*) (f2->Get("MyAlibavaClustering/ClusterSize"))->Clone(); // valgrind complains
		TH1D* clusteretadistribution = (TH1D*) (f2->Get("MyAlibavaClustering/EtaDistribution"))->Clone();
		TH1D* clustersignalhisto = (TH1D*) (f2->Get("MyAlibavaClustering/SignalfromClusters"))->Clone();

		// set the run vars we write out later on...

		double clustersize = 0.0;
		double clustersizeerror = 0.0;
		double clustercount = 0.0;
		double clustercounterror = 0.0;
		double signaltonoise = 0.0;
		double signaltonoiseerror = 0.0;

		clustersize = clustersizehisto->GetMean();
		clustersizeerror = clustersizehisto->GetRMS();

		// the temperature correction
		float tempcorscale= 1.0;
		if (_cut_applytempcorrection == true)
		{
			if (_temperature.at(irun) == 20)
			{
				tempcorscale = 1.19110;
				cout << "Applying temperature correction!" << endl;
			}
		}

		// declare all the vars
		// general
		int Event, RunNr, EvtNr, Ndf;
		float Chi2;

		// measurements and fits
		double measX_0, measY_0, measZ_0, measQ_0, fitX_0, fitY_0, fitZ_0;
		double measX_1, measY_1, measZ_1, measQ_1, fitX_1, fitY_1, fitZ_1;
		double measX_2, measY_2, measZ_2, measQ_2, fitX_2, fitY_2, fitZ_2;
		double measX_3, measY_3, measZ_3, measQ_3, fitX_3, fitY_3, fitZ_3;
		double measX_4, measY_4, measZ_4, measQ_4, fitX_4, fitY_4, fitZ_4;
		double measX_5, measY_5, measZ_5, measQ_5, fitX_5, fitY_5, fitZ_5;
		double measX_6, measY_6, measZ_6, measQ_6, fitX_6, fitY_6, fitZ_6;

		// the track point on the dut
		double dutTrackX_global, dutTrackY_global, dutTrackZ_global;
		double dutTrackX_local, dutTrackY_local, dutTrackZ_local;
		double dutTrackX_pixel, dutTrackY_pixel, dutTrackZ_pixel;

		// the hit on the dut if it was matched
		double dutHitX_global, dutHitY_global, dutHitZ_global;
		double dutHitX_local, dutHitY_local, dutHitZ_local;
		double dutHitX_pixel, dutHitY_pixel, dutHitZ_pixel;

		// misc hit info
		double dutHitR, dutHitQ;

		// alibava header stuff
		float alibava_tdc, alibava_temp;

		// and the alibava reco data
		double alibava_reco_ch_[128];

		// go to position...
		TTree *ttel = (TTree*)f1->Get("Ntuple/EUFit");  // valgrind complains

		// ... and read!
		ttel->SetBranchAddress("Event", &Event);
		ttel->SetBranchAddress("RunNr", &RunNr);
		ttel->SetBranchAddress("EvtNr", &EvtNr);
		ttel->SetBranchAddress("Ndf", &Ndf);
		ttel->SetBranchAddress("Chi2", &Chi2);
		ttel->SetBranchAddress("measX_0", &measX_0);
		ttel->SetBranchAddress("measY_0", &measY_0);
		ttel->SetBranchAddress("measZ_0", &measZ_0);
		ttel->SetBranchAddress("measQ_0", &measQ_0);
		ttel->SetBranchAddress("fitX_0", &fitX_0);
		ttel->SetBranchAddress("fitY_0", &fitY_0);
		ttel->SetBranchAddress("fitZ_0", &fitZ_0);
		ttel->SetBranchAddress("measX_1", &measX_1);
		ttel->SetBranchAddress("measY_1", &measY_1);
		ttel->SetBranchAddress("measZ_1", &measZ_1);
		ttel->SetBranchAddress("measQ_1", &measQ_1);
		ttel->SetBranchAddress("fitX_1", &fitX_1);
		ttel->SetBranchAddress("fitY_1", &fitY_1);
		ttel->SetBranchAddress("fitZ_1", &fitZ_1);
		ttel->SetBranchAddress("measX_2", &measX_2);
		ttel->SetBranchAddress("measY_2", &measY_2);
		ttel->SetBranchAddress("measZ_2", &measZ_2);
		ttel->SetBranchAddress("measQ_2", &measQ_2);
		ttel->SetBranchAddress("fitX_2", &fitX_2);
		ttel->SetBranchAddress("fitY_2", &fitY_2);
		ttel->SetBranchAddress("fitZ_2", &fitZ_2);
		ttel->SetBranchAddress("measX_3", &measX_3);
		ttel->SetBranchAddress("measY_3", &measY_3);
		ttel->SetBranchAddress("measZ_3", &measZ_3);
		ttel->SetBranchAddress("measQ_3", &measQ_3);
		ttel->SetBranchAddress("fitX_3", &fitX_3);
		ttel->SetBranchAddress("fitY_3", &fitY_3);
		ttel->SetBranchAddress("fitZ_3", &fitZ_3);
		ttel->SetBranchAddress("measX_4", &measX_4);
		ttel->SetBranchAddress("measY_4", &measY_4);
		ttel->SetBranchAddress("measZ_4", &measZ_4);
		ttel->SetBranchAddress("measQ_4", &measQ_4);
		ttel->SetBranchAddress("fitX_4", &fitX_4);
		ttel->SetBranchAddress("fitY_4", &fitY_4);
		ttel->SetBranchAddress("fitZ_4", &fitZ_4);
		ttel->SetBranchAddress("measX_5", &measX_5);
		ttel->SetBranchAddress("measY_5", &measY_5);
		ttel->SetBranchAddress("measZ_5", &measZ_5);
		ttel->SetBranchAddress("measQ_5", &measQ_5);
		ttel->SetBranchAddress("fitX_5", &fitX_5);
		ttel->SetBranchAddress("fitY_5", &fitY_5);
		ttel->SetBranchAddress("fitZ_5", &fitZ_5);
		ttel->SetBranchAddress("measX_6", &measX_6);
		ttel->SetBranchAddress("measY_6", &measY_6);
		ttel->SetBranchAddress("measZ_6", &measZ_6);
		ttel->SetBranchAddress("measQ_6", &measQ_6);
		ttel->SetBranchAddress("fitX_6", &fitX_6);
		ttel->SetBranchAddress("fitY_6", &fitY_6);
		ttel->SetBranchAddress("fitZ_6", &fitZ_6);
		ttel->SetBranchAddress("dutTrackX_global", &dutTrackX_global);
		ttel->SetBranchAddress("dutTrackY_global", &dutTrackY_global);
		ttel->SetBranchAddress("dutTrackZ_global", &dutTrackZ_global);
		ttel->SetBranchAddress("dutTrackX_local",  &dutTrackX_local);
		ttel->SetBranchAddress("dutTrackY_local",  &dutTrackY_local);
		ttel->SetBranchAddress("dutTrackZ_local",  &dutTrackZ_local);
		ttel->SetBranchAddress("dutTrackX_pixel",  &dutTrackX_pixel);
		ttel->SetBranchAddress("dutTrackY_pixel",  &dutTrackY_pixel);
		ttel->SetBranchAddress("dutTrackZ_pixel",  &dutTrackZ_pixel);
		ttel->SetBranchAddress("dutHitX_global",   &dutHitX_global);
		ttel->SetBranchAddress("dutHitY_global",   &dutHitY_global);
		ttel->SetBranchAddress("dutHitZ_global",   &dutHitZ_global);
		ttel->SetBranchAddress("dutHitX_local",    &dutHitX_local);
		ttel->SetBranchAddress("dutHitY_local",    &dutHitY_local);
		ttel->SetBranchAddress("dutHitZ_local",    &dutHitZ_local);
		ttel->SetBranchAddress("dutHitX_pixel",    &dutHitX_pixel);
		ttel->SetBranchAddress("dutHitY_pixel",    &dutHitY_pixel);
		ttel->SetBranchAddress("dutHitZ_pixel",    &dutHitZ_pixel);
		ttel->SetBranchAddress("dutHitR",          &dutHitR);
		ttel->SetBranchAddress("dutHitQ",          &dutHitQ);
		ttel->SetBranchAddress("alibava_tdc",      &alibava_tdc);
		ttel->SetBranchAddress("alibava_temp",     &alibava_temp);

		for(int i = 0; i < 128; ++i)
		{
			sprintf(name, "alibava_reco_ch_%i", i);
			ttel->SetBranchAddress(name, &alibava_reco_ch_[i]);
		}

		// set all to inactive, only the ones we need to 1:
		ttel->SetBranchStatus("*", 0);
		ttel->SetBranchStatus("RunNr", 1);
		ttel->SetBranchStatus("EvtNr", 1);
		//ttel->SetBranchStatus("Event", 1);
		ttel->SetBranchStatus("alibava*", 1);
		ttel->SetBranchStatus("dutTrackX_*", 1);
		ttel->SetBranchStatus("dutTrackY_*", 1);
		ttel->SetBranchStatus("dutTrackZ_*", 1);
		ttel->SetBranchStatus("dutHitX_global", 1);
		ttel->SetBranchStatus("dutHitY_global", 1);
		ttel->SetBranchStatus("dutHitZ_global", 1);
		//ttel->SetBranchStatus("dutHitY_pixel", 1);
		ttel->SetBranchStatus("dutHitR", 1);
		ttel->SetBranchStatus("dutHitQ", 1);
		//ttel->SetBranchStatus("measX_*", 1);
		//ttel->SetBranchStatus("measY_*", 1);
		//ttel->SetBranchStatus("measZ_*", 1);
		ttel->SetBranchStatus("fitX_0", 1);
		ttel->SetBranchStatus("fitY_0", 1);
		ttel->SetBranchStatus("fitZ_0", 1);
		ttel->SetBranchStatus("fitX_2", 1);
		ttel->SetBranchStatus("fitY_2", 1);
		ttel->SetBranchStatus("fitZ_2", 1);
		ttel->SetBranchStatus("fitX_4", 1);
		ttel->SetBranchStatus("fitY_4", 1);
		ttel->SetBranchStatus("fitZ_4", 1);
		ttel->SetBranchStatus("fitX_6", 1);
		ttel->SetBranchStatus("fitY_6", 1);
		ttel->SetBranchStatus("fitZ_6", 1);

		// how many entries in this tuple? -> tracks!
		long int tupleentrycount = ttel->GetEntries();

		// define the loop limit
		int tuplelimit = 0;
		if (tupleentrycount >= _maxtotaltracks)
		{
			tuplelimit = _maxtotaltracks; 
		}
		if (_maxtotaltracks > tupleentrycount)
		{
			tuplelimit = tupleentrycount;
		}

		// the run we are in
		int j = _runnumber.at(irun);

		// the histos to be filled

		TH1D* noise[128];
		TH1D* noiseelec[128];
		for (int ii=0;ii<128;ii++)
		{
			histoname = "noise_";
			sstream << histoname << j << "_chan_" << ii;
			noise[ii] = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());

			histoname = "noiseelec_";
			sstream << histoname << j << "_chan_" << ii;
			noiseelec[ii] = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
		}

		TH1D* tdceval[91];
		TH1D* tdcevalelec[91];
		for (int ii=0;ii<91;ii++)
		{
			histoname = "tdceval_";
			sstream << histoname << j << "_tdc_" << ii;
			tdceval[ii] = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());

			histoname = "tdcevalelec_";
			sstream << histoname << j << "_tdc_" << ii;
			tdcevalelec[ii] = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
		}

		TH1D* posxeval[100];
		for (int ii=0;ii<100;ii++)
		{
			histoname = "posxeval_";
			sstream << histoname << j << "_pos_" << ii;
			posxeval[ii] = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
		}

		histoname = "overalltdc_";
		sstream << histoname << j ;
		TGraphErrors* overalltdc = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "overallxpos_";
		sstream << histoname << j ;
		TGraphErrors* overallxpos = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "tdcdistri_";
		sstream << histoname << j ;
		TH1D* tdcdistri = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "tempdistri_";
		sstream << histoname << j ;
		TH1D* tempdistri = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "allnoise_";
		sstream << histoname << j ;
		TH1D* allnoise = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "allnoiseelec_";
		sstream << histoname << j ;
		TH1D* allnoiseelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "adcnoise_";
		sstream << histoname << j ;
		TH1D* adcnoise = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "adcnoiseelec_";
		sstream << histoname << j ;
		TH1D* adcnoiseelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "fivechannoise_";
		sstream << histoname << j ;
		TH1D* fivechannoise = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "fivechannoiseelec_";
		sstream << histoname << j ;
		TH1D* fivechannoiseelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsX_";
		sstream << histoname << j ;
		TH1D* residualsX = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsY_";
		sstream << histoname << j ;
		TH1D* residualsY = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsZ_";
		sstream << histoname << j ;
		TH1D* residualsZ = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsXY_";
		sstream << histoname << j ;
		TH2D* residualsXY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsXZ_";
		sstream << histoname << j ;
		TH2D* residualsXZ = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsYZ_";
		sstream << histoname << j ;
		TH2D* residualsYZ = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsYQ_";
		sstream << histoname << j ;
		TH2D* residualsYQ = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsYR_";
		sstream << histoname << j ;
		TH2D* residualsYR = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsYT_";
		sstream << histoname << j ;
		TH2D* residualsYT = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsY_Q_";
		sstream << histoname << j ;
		TH2D* residualsY_Q = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsY_clu1_";
		sstream << histoname << j ;
		TH1D* residualsY_clu1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsY_clu2_";
		sstream << histoname << j ;
		TH1D* residualsY_clu2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsY_clu3_";
		sstream << histoname << j ;
		TH1D* residualsY_clu3 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsY_clu4_";
		sstream << histoname << j ;
		TH1D* residualsY_clu4 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDXvsX_";
		sstream << histoname << j ;
		TH2D* residualsDXvsX = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDXvsY_";
		sstream << histoname << j ;
		TH2D* residualsDXvsY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDXvsZ_";
		sstream << histoname << j ;
		TH2D* residualsDXvsZ = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDYvsX_";
		sstream << histoname << j ;
		TH2D* residualsDYvsX = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDYvsY_";
		sstream << histoname << j ;
		TH2D* residualsDYvsY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDYvsZ_";
		sstream << histoname << j ;
		TH2D* residualsDYvsZ = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDZvsX_";
		sstream << histoname << j ;
		TH2D* residualsDZvsX = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDZvsY_";
		sstream << histoname << j ;
		TH2D* residualsDZvsY = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDZvsZ_";
		sstream << histoname << j ;
		TH2D* residualsDZvsZ = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDXvsSXU_";
		sstream << histoname << j ;
		TH2D* residualsDXvsSXU = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDXvsSYU_";
		sstream << histoname << j ;
		TH2D* residualsDXvsSYU = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDYvsSXU_";
		sstream << histoname << j ;
		TH2D* residualsDYvsSXU = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDYvsSYU_";
		sstream << histoname << j ;
		TH2D* residualsDYvsSYU = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDZvsSXU_";
		sstream << histoname << j ;
		TH2D* residualsDZvsSXU = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDZvsSYU_";
		sstream << histoname << j ;
		TH2D* residualsDZvsSYU = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDXvsSXD_";
		sstream << histoname << j ;
		TH2D* residualsDXvsSXD = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDXvsSYD_";
		sstream << histoname << j ;
		TH2D* residualsDXvsSYD = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDYvsSXD_";
		sstream << histoname << j ;
		TH2D* residualsDYvsSXD = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDYvsSYD_";
		sstream << histoname << j ;
		TH2D* residualsDYvsSYD = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDZvsSXD_";
		sstream << histoname << j ;
		TH2D* residualsDZvsSXD = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsDZvsSYD_";
		sstream << histoname << j ;
		TH2D* residualsDZvsSYD = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualProfileDXvsX_";
		sstream << histoname << j ;
		TProfile* residualProfileDXvsX = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualProfileDXvsY_";
		sstream << histoname << j ;
		TProfile* residualProfileDXvsY = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());


		histoname = "residualProfileDYvsX_";
		sstream << histoname << j ;
		TProfile* residualProfileDYvsX = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualProfileDYvsY_";
		sstream << histoname << j ;
		TProfile* residualProfileDYvsY = dynamic_cast<TProfile*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsXvsEvt_";
		sstream << histoname << j ;
		TH2D* residualsXvsEvt = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsYvsEvt_";
		sstream << histoname << j ;
		TH2D* residualsYvsEvt = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsZvsEvt_";
		sstream << histoname << j ;
		TH2D* residualsZvsEvt = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualMapTemp_";
		sstream << histoname << j ;
		TH3D* residualMapTemp = dynamic_cast<TH3D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapDUT_";
		sstream << histoname << j ;
		TH2D* hitmapDUT = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapDUTmodpitch_";
		sstream << histoname << j ;
		TH2D* hitmapDUTmodpitch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapTracks_";
		sstream << histoname << j ;
		TH2D* hitmapTracks = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapTracksmodpitch_";
		sstream << histoname << j ;
		TH2D* hitmapTracksmodpitch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapMatch_";
		sstream << histoname << j ;
		TH2D* hitmapMatch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapMatchmodpitch_";
		sstream << histoname << j ;
		TH2D* hitmapMatchmodpitch = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapMatchTracks_";
		sstream << histoname << j ;
		TH2D* hitmapMatchTracks = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapChargeShared0_";
		sstream << histoname << j ;
		TH2D* hitmapChargeShared0 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapChargeShared1_";
		sstream << histoname << j ;
		TH2D* hitmapChargeShared1 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapChargeShared2_";
		sstream << histoname << j ;
		TH2D* hitmapChargeShared2 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapChargeShared3_";
		sstream << histoname << j ;
		TH2D* hitmapChargeShared3 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapChargeShared4_";
		sstream << histoname << j ;
		TH2D* hitmapChargeShared4 = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapClusterSizeA_";
		sstream << histoname << j ;
		TH1D* hitmapClusterSizeA = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapClusterSize1_";
		sstream << histoname << j ;
		TH1D* hitmapClusterSize1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapClusterSize2_";
		sstream << histoname << j ;
		TH1D* hitmapClusterSize2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapClusterSize3_";
		sstream << histoname << j ;
		TH1D* hitmapClusterSize3 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "hitmapClusterSize4_";
		sstream << histoname << j ;
		TH1D* hitmapClusterSize4 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize1_";
		sstream << histoname << j ;
		TH1D* signalClusterSize1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize1elec_";
		sstream << histoname << j ;
		TH1D* signalClusterSize1elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize2_";
		sstream << histoname << j ;
		TH1D* signalClusterSize2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize2elec_";
		sstream << histoname << j ;
		TH1D* signalClusterSize2elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize3_";
		sstream << histoname << j ;
		TH1D* signalClusterSize3 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize3elec_";
		sstream << histoname << j ;
		TH1D* signalClusterSize3elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize4_";
		sstream << histoname << j ;
		TH1D* signalClusterSize4 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalClusterSize4elec_";
		sstream << histoname << j ;
		TH1D* signalClusterSize4elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "scatterX_";
		sstream << histoname << j ;
		TProfile2D* scatterX = dynamic_cast<TProfile2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "scatterY_";
		sstream << histoname << j ;
		TProfile2D* scatterY = dynamic_cast<TProfile2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "scatterXY_";
		sstream << histoname << j ;
		TProfile2D* scatterXY = dynamic_cast<TProfile2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "kinkX_";
		sstream << histoname << j ;
		TH1D* kinkX = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "kinkY_";
		sstream << histoname << j ;
		TH1D* kinkY = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "matchedEta_";
		sstream << histoname << j ;
		TH1D* matchedEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "matchedEta2D_";
		sstream << histoname << j ;
		TH2D* matchedEta2D = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "matchedEtaIntegral_";
		sstream << histoname << j ;
		TH1D* matchedEtaIntegral = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "unmatchedEta_";
		sstream << histoname << j ;
		TH1D* unmatchedEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "unmatchedEta2D_";
		sstream << histoname << j ;
		TH2D* unmatchedEta2D = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "unmatchedEtaIntegral_";
		sstream << histoname << j ;
		TH1D* unmatchedEtaIntegral = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "striphitEta_";
		sstream << histoname << j ;
		TH1D* striphitEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "striphitEta2D_";
		sstream << histoname << j ;
		TH2D* striphitEta2D = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "striphitEtaIntegral_";
		sstream << histoname << j ;
		TH1D* striphitEtaIntegral = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "noiseEta_";
		sstream << histoname << j ;
		TH1D* noiseEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "noisesubtractedEta_";
		sstream << histoname << j ;
		TH1D* noisesubtractedEta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaR10_";
		sstream << histoname << j ;
		TH1D* EtaR10 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaR20_";
		sstream << histoname << j ;
		TH1D* EtaR20 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaR30_";
		sstream << histoname << j ;
		TH1D* EtaR30 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaR40_";
		sstream << histoname << j ;
		TH1D* EtaR40 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaR50_";
		sstream << histoname << j ;
		TH1D* EtaR50 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaL10_";
		sstream << histoname << j ;
		TH1D* EtaL10 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaL20_";
		sstream << histoname << j ;
		TH1D* EtaL20 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaL30_";
		sstream << histoname << j ;
		TH1D* EtaL30 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaL40_";
		sstream << histoname << j ;
		TH1D* EtaL40 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "EtaL50_";
		sstream << histoname << j ;
		TH1D* EtaL50 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackSignal_";
		sstream << histoname << j ;
		TH1D* trackSignal = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackSignalelec_";
		sstream << histoname << j ;
		TH1D* trackSignalelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackSignalMap_";
		sstream << histoname << j ;
		TH2D* trackSignalMap = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackSignalMapelec_";
		sstream << histoname << j ;
		TH2D* trackSignalMapelec = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackSignalTDC_";
		sstream << histoname << j ;
		TH2D* trackSignalTDC = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackSignalTDCelec_";
		sstream << histoname << j ;
		TH2D* trackSignalTDCelec = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "tracksperevent_";
		sstream << histoname << j ;
		TH1D* tracksperevent = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "tracksvsevents_";
		sstream << histoname << j ;
		TH2D* tracksvsevents = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "highchanneldistribution_";
		sstream << histoname << j ;
		TH1D* highchanneldistribution = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalLeft2_";
		sstream << histoname << j ;
		TH1D* signalLeft2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalLeft2elec_";
		sstream << histoname << j ;
		TH1D* signalLeft2elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalLeft1_";
		sstream << histoname << j ;
		TH1D* signalLeft1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalLeft1elec_";
		sstream << histoname << j ;
		TH1D* signalLeft1elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalCenter_";
		sstream << histoname << j ;
		TH1D* signalCenter = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalCenterelec_";
		sstream << histoname << j ;
		TH1D* signalCenterelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalRight1_";
		sstream << histoname << j ;
		TH1D* signalRight1 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalRight1elec_";
		sstream << histoname << j ;
		TH1D* signalRight1elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalRight2_";
		sstream << histoname << j ;
		TH1D* signalRight2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalRight2elec_";
		sstream << histoname << j ;
		TH1D* signalRight2elec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalGoodEvents_";
		sstream << histoname << j ;
		TH1D* signalGoodEvents = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalGoodEventselec_";
		sstream << histoname << j ;
		TH1D* signalGoodEventselec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapA_";
		sstream << histoname << j ;
		TH1D* signalmapA = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapAelec_";
		sstream << histoname << j ;
		TH1D* signalmapAelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapB_";
		sstream << histoname << j ;
		TH1D* signalmapB = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapBelec_";
		sstream << histoname << j ;
		TH1D* signalmapBelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapC_";
		sstream << histoname << j ;
		TH1D* signalmapC = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapCelec_";
		sstream << histoname << j ;
		TH1D* signalmapCelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapD_";
		sstream << histoname << j ;
		TH1D* signalmapD = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalmapDelec_";
		sstream << histoname << j ;
		TH1D* signalmapDelec = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalareaplot_";
		sstream << histoname << j;
		TGraphErrors* signalareaplot = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalareaplotelec_";
		sstream << histoname << j;
		TGraphErrors* signalareaplotelec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "fiducial_discard_";
		sstream << histoname << j;
		TH1D* fiducial_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "goodchannel_discard_";
		sstream << histoname << j;
		TH1D* goodchannel_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "goodevent_discard_";
		sstream << histoname << j;
		TH1D* goodevent_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackselection_discard_";
		sstream << histoname << j;
		TH1D* trackselection_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "timecut_discard_";
		sstream << histoname << j;
		TH1D* timecut_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "highchannel_discard_";
		sstream << histoname << j;
		TH1D* highchannel_discard = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "fiducial_allow_";
		sstream << histoname << j;
		TH1D* fiducial_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "goodchannel_allow_";
		sstream << histoname << j;
		TH1D* goodchannel_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "nohittrack_";
		sstream << histoname << j;
		TH1D* nohittrack = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "goodevent_allow_";
		sstream << histoname << j;
		TH1D* goodevent_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "trackselection_allow_";
		sstream << histoname << j;
		TH1D* trackselection_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "timecut_allow_";
		sstream << histoname << j;
		TH1D* timecut_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "highchannel_allow_";
		sstream << histoname << j;
		TH1D* highchannel_allow = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalMapTemp_";
		sstream << histoname << j ;
		TH3D* signalMapTemp = dynamic_cast<TH3D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		// group some runs together!
		// find the right histo!
		std::vector<string>::iterator its;
		its = find(_sensorname_list.begin(), _sensorname_list.end(), _sensorname.at(irun));
		std::string::size_type pos1 = its - _sensorname_list.begin();

		std::vector<float>::iterator itf;
		itf = find(_irradfluence_list.begin(), _irradfluence_list.end(), _irradfluence.at(irun));
		size_t pos2 = itf - _irradfluence_list.begin();

		std::vector<int>::iterator it;
		it = find(_dutrotation_list.begin(), _dutrotation_list.end(), _dutrotation.at(irun));
		size_t pos3 = it - _dutrotation_list.begin();

		// if plotting quantity vs voltage the encoding is: sensor - fluence - rotation
		int vsvoltage_histonr = pos1*_irradfluence_total*_dutrotation_total + pos2*_dutrotation_total + pos3;

		// posN must be inc'ed for sensible output, as it starts at 0
		pos1++;
		pos2++;
		pos3++;

		if (_debug<=0)
		{
			cout << " " << endl;
			cout << "Sensor:   " << pos1 << " of " << _sensorname_total << endl;
			cout << "Fluence:  " << pos2 << " of " << _irradfluence_total << endl;
			cout << "Rotation: " << pos3 << " of " << _dutrotation_total << endl;
		}

		voltagegraphs_total.push_back(vsvoltage_histonr);

		it = find (voltagegraphs.begin(), voltagegraphs.end(), vsvoltage_histonr);
		if (it != voltagegraphs.end())
		{
			if (_debug <=0)
			{
				cout << "Voltage graph " << vsvoltage_histonr << " already in list!" << endl;
			}
		} else {
			voltagegraphs.push_back(vsvoltage_histonr);
			if (_debug <=0)
			{
				cout << "Voltage graph " << vsvoltage_histonr << " not in list!" << endl;
			}
		}

		// protection against crashes, FIXME
		if (vsvoltage_histonr > _bookinglimit)
		{
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " WARNING: Increase booking limit to at least "<< vsvoltage_histonr << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			break;
		}

		histoname = "residualsXvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* residualsXvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsYvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* residualsYvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "residualsYvoltageAngle_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* residualsYvoltageAngle = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* signalvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalvoltageelec_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* signalvoltageelec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "snvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* snvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "noisevoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* noisevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "noisevoltageelec_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* noisevoltageelec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "rghvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* rghvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "clustersizevoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "clustercountvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* clustercountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "channelcountvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* channelcountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "currentvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* currentvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "volumecurrent_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* volumecurrent = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "surfacecurrent_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* surfacecurrent = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "chargesharingvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* chargesharingvoltage= dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "etachargesharingvoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* etachargesharingvoltage= dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signaldistancevoltage_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* signaldistancevoltage= dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signaldistancevoltageelec_";
		sstream << histoname << vsvoltage_histonr;
		TGraphErrors* signaldistancevoltageelec= dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		histoname = "signalareamap_";
		sstream << histoname << vsvoltage_histonr;
		TH2D* signalareamap = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		/*
		 * 
		 * Loop to get the track multiplicity in an event
		 * 
		// var for the event nr so we can cut


		cout << "start multi loop " << endl;
		// track multiplicity loop
		int tempcurrentevent = -1;
		int trackmulti = -1;
		std::vector<int> multiplicity;
		for(int i = 0; i< tuplelimit;i++)
		{
			// get the event
			ttel->GetEntry(i);

			// easy case: same event as before, inc count
			if (EvtNr == tempcurrentevent)
			{
				trackmulti++;
			}

			// a new event
			if (EvtNr > tempcurrentevent )
			{
				// how many track-less events have passed?
				if (EvtNr > 0)
				{
					multiplicity.push_back(trackmulti);
					for (int j = 1; j<(EvtNr-tempcurrentevent);j++)
					{
						multiplicity.push_back(0);
					}

				}
				// new
				tempcurrentevent = EvtNr;
				trackmulti = 1;
			}
		}
		*/

		TStopwatch runtime;
		TStopwatch noisetime;

		if (_debug <= 2)
		{
			cout << " " << endl;
			cout << "Starting hit tuple loop..." << endl;
		}
		runtime.Start(true);

		float goodtrackcount = 0.0;
		float matchedhitcount = 0.0;

		int trackpereventcut = 0;
		int fiducialcut = 0;
		int pixelcut = 0;
		int finalcut = 0;

		int currentevent = -1;
		int tracksinthisevent = 1;
		bool isthisfirstevent = false;
		bool goodevent = false;

		// loop over the tuple contents and fill
		for(int i = 0; i< tuplelimit;i++)
		{

			ttel->GetEntry(i);
			// start filling histos

			float eventtemp = 0.0;
			eventtemp = alibava_temp;

			if (_temperature.at(irun) < 0 && _cut_applytempcorrection == true && RunNr < 4000)
			{
				float targetgain = 0.004700854;
				tempcorscale = targetgain / (0.00429259 + 7.9121e-5 * eventtemp - 3.82946e-6 * eventtemp*eventtemp);
			}
			if (_cut_applytempcorrection == true && RunNr > 4000)
			{
			     tempcorscale = 1 / (-0.0074953853 * eventtemp + 1.0763859491);
			}

			// first find out where we are
			if (EvtNr == currentevent)
			{
				// this is not the first track in the event
				tracksinthisevent++;

				// this is now false
				isthisfirstevent = false;

			} else {
				// this IS the first track in the event
				tdcdistri->Fill(alibava_tdc);
				tempdistri->Fill(eventtemp);

				// catch events with 0 tracks
				if (EvtNr > 0)
				{
					for (int j = 1; j<(EvtNr - currentevent); j++)
					{
						tracksperevent->Fill(0);
						tracksvsevents->Fill((currentevent+j),0);
					}
				}

				// tracksinthisevent still has the track count from the previous event, so we fill the histo now!
				tracksperevent->Fill(tracksinthisevent);
				tracksvsevents->Fill(currentevent,tracksinthisevent);

				// did all tracks in the previous event fullfill the cuts? then we have a good event!
				if (goodevent == true)
				{
					for (int j=0;j<tracksinthisevent;j++)
					{
						_goodtracks.push_back(i-1-j);
					}

					goodevent_allow->Fill(1.0);

				} else {
					goodevent_discard->Fill(1.0);
				}

				// first event: reset
				goodevent = true;

				// this is now true
				isthisfirstevent = true;

				// reset
				tracksinthisevent = 1;
				currentevent = EvtNr;

			}

			// adding 0.5 makes the rounding correct!
			// this will hold the integer channel nearest to the track

			float tempfloat = dutTrackY_pixel + 0.5;
			if (xsensitive == true)
			{
				tempfloat = dutTrackX_pixel + 0.5;
			}
			int pixelpos = static_cast <int>(tempfloat);

			// pixelbin: hit is between .0 and .999 -> this is used for eta left right
			// this will hold the integer channel before the decimal point
			int pixelbin = static_cast <int>(dutTrackY_pixel);
			if (xsensitive == true)
			{
				pixelbin = static_cast <int>(dutTrackX_pixel);
			}

			// the number of channels we have summed for 5chan noise_
			int channelssummed = 0;
			float sumnoise = 0.0;

			// noise calculation: only take the first track per event, exclude area under track +- 3
			// noise eta: additionally require the channel and its neighbour to be good

			// noise does not get temp corrected!
			if (isthisfirstevent == true)
			{
				noisetime.Start(false);
				for (int ii=0;ii<128;ii++)
				{
					if ((ii != (pixelpos-5)) && (ii != (pixelpos-4)) && (ii != (pixelpos-3)) && (ii != (pixelpos-2)) && (ii != (pixelpos-1)) && (ii != (pixelpos)) && (ii != (pixelpos+1)) && (ii != (pixelpos+2)) && (ii != (pixelpos+3)) && (ii != (pixelpos+4)) && (ii != (pixelpos+5)))
					{
						noise[ii]->Fill(alibava_reco_ch_[ii]);
						noiseelec[ii]->Fill(alibava_reco_ch_[ii]*_electronconversion);
						if (ii > 1 && ii < 127 && _channels[irun][ii] == 1 && _channels[irun][ii+1] == 1)
						{
							adcnoise->Fill(alibava_reco_ch_[ii]);
							adcnoiseelec->Fill(alibava_reco_ch_[ii]*_electronconversion);
							noiseEta->Fill(alibava_reco_ch_[ii+1]/(alibava_reco_ch_[ii]+alibava_reco_ch_[ii+1]));
						}
						if (channelssummed <5)
						{
							if (_channels[irun][ii] == 1)
							{
								sumnoise+=alibava_reco_ch_[ii];
								channelssummed++;
							}
						}
						if (channelssummed == 5)
						{
							fivechannoise->Fill(sumnoise*_polarity.at(irun));
							fivechannoiseelec->Fill(sumnoise*_polarity.at(irun)*_electronconversion);
							sumnoise = 0.0;
							channelssummed = 0;
						}
					}
				}
				noisetime.Stop();
			}

			double scatX = _missingvalue;
			double scatY = _missingvalue;

			// calculate the scattering of the track...
			scatX = ((fitX_6-fitX_4) / (fitZ_6-fitZ_4) - (fitX_2-fitX_0) / (fitZ_2-fitZ_0))*1E3;
			scatY = ((fitY_6-fitY_4) / (fitZ_6-fitZ_4) - (fitY_2-fitY_0) / (fitZ_2-fitZ_0))*1E3;
			kinkX->Fill(scatX);
			kinkY->Fill(scatY);
			scatterX->Fill(dutTrackX_global,dutTrackY_global,scatX);
			scatterY->Fill(dutTrackX_global,dutTrackY_global,scatY);
			scatterXY->Fill(dutTrackX_global,dutTrackY_global,sqrt(scatX*scatX+scatY*scatY));

			// cut: tracks per event
			if (tracksinthisevent <= _cut_maxtracksperevent)
			{

				// fiducial x on the track
				if (dutTrackX_global >= _sensorminx.at(irun) && dutTrackX_global <= _sensormaxx.at(irun))
				{

					// cut: existing and good channel pointed to - not really needed here...
					if (pixelpos >= 2 && pixelpos <= 125 && _channels[irun][pixelpos] == 1 && _channels[irun][pixelpos+1] == 1 && _channels[irun][pixelpos-1] == 1 && _channels[irun][pixelpos+2] == 1 && _channels[irun][pixelpos-2] == 1)
					{

						// dutHits are already preselected in time, if this is wide, then it will not affect the hits...
						// cut: in time
						//	    if (alibava_tdc > _sensormintdc.at(irun) && alibava_tdc < _sensormaxtdc.at(irun))
						//	    {

						// cut: matched hit
						// hit is from EUTel definition!!!
						// cut: not missingvalue!
						//  && (dutHitX_global) > _cut_minx && (dutHitX_global) < _cut_maxx
						if (dutHitX_global != _missingvalue && dutHitY_global != _missingvalue )//&& dutTrackX_global != _missingvalue && dutTrackY_global != _missingvalue )
						{
							// z scale
							double initzpos = 270.0;
							dutTrackZ_global = dutTrackZ_global - initzpos;
							dutHitZ_global = dutHitZ_global - initzpos;

							// Q scale, was introduced to get around integer storage in the clustering...
							double qscale = 1000.0;
							dutHitQ = dutHitQ / qscale;

							residualsX->Fill(dutHitX_global-dutTrackX_global);
							residualsY->Fill(dutHitY_global-dutTrackY_global);
							residualsZ->Fill(dutHitZ_global-dutTrackZ_global);
							residualsXY->Fill(dutHitX_global-dutTrackX_global,dutHitY_global-dutTrackY_global);
							residualsXZ->Fill(dutHitX_global-dutTrackX_global,dutHitZ_global-dutTrackZ_global);
							residualsYZ->Fill(dutHitY_global-dutTrackY_global,dutHitZ_global-dutTrackZ_global);

							residualsYQ->Fill(dutHitY_global-dutTrackY_global,dutHitQ);
							residualsYR->Fill(dutHitY_global-dutTrackY_global,dutHitR);
							residualsYT->Fill(dutHitY_global-dutTrackY_global,alibava_tdc);

							residualsDXvsX->Fill(dutHitX_global-dutTrackX_global,dutTrackX_global);
							residualsDXvsY->Fill(dutHitX_global-dutTrackX_global,dutTrackY_global);
							residualsDXvsZ->Fill(dutHitX_global-dutTrackX_global,dutTrackZ_global);
							residualsDYvsX->Fill(dutHitY_global-dutTrackY_global,dutTrackX_global);
							residualsDYvsY->Fill(dutHitY_global-dutTrackY_global,dutTrackY_global);
							residualsDYvsZ->Fill(dutHitY_global-dutTrackY_global,dutTrackZ_global);
							residualsDZvsX->Fill(dutHitZ_global-dutTrackZ_global,dutTrackX_global);
							residualsDZvsY->Fill(dutHitZ_global-dutTrackZ_global,dutTrackY_global);
							residualsDZvsZ->Fill(dutHitZ_global-dutTrackZ_global,dutTrackZ_global);

							residualsDXvsSXU->Fill(dutHitX_global-dutTrackX_global,(fitX_0-fitX_2)/(fitZ_0 - fitZ_2)*1E3);
							residualsDXvsSYU->Fill(dutHitX_global-dutTrackX_global,(fitY_0-fitY_2)/(fitZ_0 - fitZ_2)*1E3);
							residualsDYvsSXU->Fill(dutHitY_global-dutTrackY_global,(fitX_0-fitX_2)/(fitZ_0 - fitZ_2)*1E3);
							residualsDYvsSYU->Fill(dutHitY_global-dutTrackY_global,(fitY_0-fitY_2)/(fitZ_0 - fitZ_2)*1E3);
							residualsDZvsSXU->Fill(dutHitZ_global-dutTrackZ_global,(fitX_0-fitX_2)/(fitZ_0 - fitZ_2)*1E3);
							residualsDZvsSYU->Fill(dutHitZ_global-dutTrackZ_global,(fitY_0-fitY_2)/(fitZ_0 - fitZ_2)*1E3);

							residualsDXvsSXD->Fill(dutHitX_global-dutTrackX_global,(fitX_4-fitX_6)/(fitZ_4 - fitZ_6)*1E3);
							residualsDXvsSYD->Fill(dutHitX_global-dutTrackX_global,(fitY_4-fitY_6)/(fitZ_4 - fitZ_6)*1E3);
							residualsDYvsSXD->Fill(dutHitY_global-dutTrackY_global,(fitX_4-fitX_6)/(fitZ_4 - fitZ_6)*1E3);
							residualsDYvsSYD->Fill(dutHitY_global-dutTrackY_global,(fitY_4-fitY_6)/(fitZ_4 - fitZ_6)*1E3);
							residualsDZvsSXD->Fill(dutHitZ_global-dutTrackZ_global,(fitX_4-fitX_6)/(fitZ_4 - fitZ_6)*1E3);
							residualsDZvsSYD->Fill(dutHitZ_global-dutTrackZ_global,(fitY_4-fitY_6)/(fitZ_4 - fitZ_6)*1E3);

							residualsXvsEvt->Fill(EvtNr,dutHitX_global-dutTrackX_global);
							residualsYvsEvt->Fill(EvtNr,dutHitY_global-dutTrackY_global);
							residualsZvsEvt->Fill(EvtNr,dutHitZ_global-dutTrackZ_global);

							residualProfileDXvsX->Fill(dutTrackX_global,dutHitX_global-dutTrackX_global,1);
							residualProfileDXvsY->Fill(dutTrackY_global,dutHitX_global-dutTrackX_global,1);
							residualProfileDYvsX->Fill(dutTrackX_global,dutHitY_global-dutTrackY_global,1);
							residualProfileDYvsY->Fill(dutTrackY_global,dutHitY_global-dutTrackY_global,1);

							residualMapTemp->Fill(dutTrackX_global,dutTrackY_global,dutHitY_global-dutTrackY_global);

							hitmapMatch->Fill(dutHitX_global,dutHitY_global);
							hitmapMatchTracks->Fill(dutTrackX_global,dutTrackY_global);

							fractpartx = modf (dutHitX_global/_pitchx , &intpartx);
							fractparty = modf (dutHitY_global/_pitchy , &intparty);

							hitmapMatchmodpitch->Fill(fabs(fractpartx),fabs(fractparty));

							matchedEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]));
							fractparty = modf (dutTrackY_pixel, &intparty);
							matchedEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]),fractparty);

							matchedhitcount++;

							double tempsignal = tempcorscale*alibava_reco_ch_[pixelbin-2] + tempcorscale*alibava_reco_ch_[pixelbin-1] + tempcorscale*alibava_reco_ch_[pixelbin] + tempcorscale*alibava_reco_ch_[pixelbin+1] + tempcorscale*alibava_reco_ch_[pixelbin+2];
							fractparty = modf (dutTrackY_pixel, &intparty);
							if (dutHitR > 0.0)
							{
								hitmapClusterSizeA->Fill(fractparty);
							}
							if (dutHitR >= 1.0 && dutHitR < 2.0)
							{
								hitmapClusterSize1->Fill(fractparty);
								signalClusterSize1->Fill(tempsignal);
								signalClusterSize1elec->Fill(tempsignal*_electronconversion);
								residualsY_clu1->Fill(dutHitY_global-dutTrackY_global);
							}
							if (dutHitR >= 2.0 && dutHitR < 3.0)
							{
								hitmapClusterSize2->Fill(fractparty);
								signalClusterSize2->Fill(tempsignal);
								signalClusterSize2elec->Fill(tempsignal*_electronconversion);
								residualsY_clu2->Fill(dutHitY_global-dutTrackY_global);
							}
							if (dutHitR >= 3.0 && dutHitR < 4.0)
							{
								hitmapClusterSize3->Fill(fractparty);
								signalClusterSize3->Fill(tempsignal);
								signalClusterSize3elec->Fill(tempsignal*_electronconversion);
								residualsY_clu3->Fill(dutHitY_global-dutTrackY_global);
							}
							if (dutHitR >= 4.0)
							{
								hitmapClusterSize4->Fill(fractparty);
								signalClusterSize4->Fill(tempsignal);
								signalClusterSize4elec->Fill(tempsignal*_electronconversion);
								residualsY_clu4->Fill(dutHitY_global-dutTrackY_global);
							}

						// done EUTel hit... from now on we don't care about them...
						} else {
							finalcut++;
							nohittrack->Fill(1.0);
						}

						// cut: track passes through a strip
						fractparty = modf (dutTrackY_pixel, &intparty);
						if (fractparty <= _cut_onstrip)
						{
							striphitEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin-1]+tempcorscale*alibava_reco_ch_[pixelbin+1]));
							striphitEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin-1]+tempcorscale*alibava_reco_ch_[pixelbin+1]),fractparty);
						}
						if (fractparty >= _cut_onstrip+0.5)
						{
							striphitEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+2]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+2]));
							striphitEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+2]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+2]),fractparty);
						}

						goodtrackcount++;

						//	    } // done time cut, this does not kill a event!

						goodchannel_allow->Fill(1.0);

					} else { // done existing good channel

						goodevent = false;
						goodchannel_discard->Fill(1.0);
						pixelcut++;
					}

					fiducial_allow->Fill(1.0);

				} else {// done fiducial x cut

					goodevent = false;
					fiducial_discard->Fill(1.0);
					fiducialcut++;

				}

				// all tracks: no cuts at all:
				hitmapTracks->Fill(dutTrackX_global,dutTrackY_global);
				fractpartx = modf (dutTrackX_global/_pitchx , &intpartx);
				fractparty = modf (dutTrackY_global/_pitchy , &intparty);
				hitmapTracksmodpitch->Fill(fabs(fractpartx),fabs(fractparty));


				// cut on non-missingvalue: all dut hits, if matched or not
				if (dutHitX_global != _missingvalue && dutHitY_global != _missingvalue)
				{
					hitmapDUT->Fill(dutHitX_global,dutHitY_global);
					fractpartx = modf (dutHitX_global/_pitchx , &intpartx);
					fractparty = modf (dutHitY_global/_pitchy , &intparty);
					hitmapDUTmodpitch->Fill(fabs(fractpartx),fabs(fractparty));
				} // done missingvalue cut

				trackselection_allow->Fill(1.0);

			} else { // tracks per event cut

				goodevent = false;
				trackpereventcut++;
				trackselection_discard->Fill(1.0);

			}

			//////////
			//
			//	All cuts should be completed!
			//
			//////////


		} // done tuple contents loop

		runtime.Stop();
		if (_debug <= 2)
		{
			cout << "Done event and hit tuple loop after " << runtime.RealTime() << " s, noise calculation " << noisetime.RealTime() << " s!" << endl;
			cout << " " << endl;
			cout << "Initial track count from EUTelescope:        " << tuplelimit << endl;
			cout << " " << endl;
			cout << "Tracks failing max tracks per event cut:     " << trackpereventcut << endl;
			cout << "Tracks failing fiducial x cut:               " << fiducialcut << endl;
			cout << "Tracks failing good channel cut:             " << pixelcut << endl;
			cout << " " << endl;
			cout << "Discarded bad tracks:                        " << trackpereventcut + fiducialcut + pixelcut << endl;
			cout << "Tracks in good events:                       " << tuplelimit - (trackpereventcut + fiducialcut + pixelcut) << endl;
			cout << " " << endl;
			cout << "Tracks in good events with no DUT hit:       " << finalcut << endl;
			cout << "Total hit track accepted count:              " << matchedhitcount << endl;
			cout << " " << endl;
		}

		// deactivate no longer needed branches
		ttel->SetBranchStatus("fit*", 0);
		ttel->SetBranchStatus("dutHit*", 0);
		ttel->SetBranchStatus("dutTrackZ_*", 0);

		double temp1 = hitmapDUT->GetMean(1);
		double temp2 = hitmapDUT->GetRMS(1);
		double temp3 = 1.5*temp2;
		double temp4 = temp1 - temp3;
		double temp5 = temp1 + temp3;
		cout << "positions: " << temp4 << " " << temp5 << endl;
		cout << " " << endl;

		if (_debug <= 2)
		{
			cout << "Starting residual mapping..." << endl;
		}

		// call the FitSlicesZ method
		// since we want to save the output, no direct function call
		//residualMapTemp->FitSlicesZ(0,0,100,0,100,1, "R") ;
		TF1 * fit1 = 0;
		Int_t binminx = 0;
		Int_t binmaxx = 50;
		Int_t binminy = 0;
		Int_t binmaxy = 50;
		Int_t cut = 1;
		Option_t *option = "QR";

		Int_t nbinsx  = residualMapTemp->GetXaxis()->GetNbins();
		Int_t nbinsy  = residualMapTemp->GetYaxis()->GetNbins();
		Int_t nbinsz  = residualMapTemp->GetZaxis()->GetNbins();

		if (binminx < 1)
		{
			binminx = 1;
		}
		if (binmaxx > nbinsx)
		{
			binmaxx = nbinsx;
		}
		if (binmaxx < binminx)
		{
			binminx = 1;
			binmaxx = nbinsx;
		}
		if (binminy < 1)
		{
			binminy = 1;
		}
		if (binmaxy > nbinsy)
		{
			binmaxy = nbinsy;
		}
		if (binmaxy < binminy)
		{
			binminy = 1;
			binmaxy = nbinsy;
		}

		//default is to fit with a gaussian
		if (fit1 == 0)
		{
			fit1 = (TF1*)gROOT->GetFunction("gaus");
			if (fit1 == 0)
			{
				fit1 = new TF1("gaus","gaus",residualMapTemp->GetZaxis()->GetXmin(),residualMapTemp->GetZaxis()->GetXmax());
			} else {
				fit1->SetRange(residualMapTemp->GetZaxis()->GetXmin(),residualMapTemp->GetZaxis()->GetXmax());
			}
		}
		const char *fname = fit1->GetName();
		//Int_t npar = fit1->GetNpar();
		Int_t npar = 3;
		Double_t *parsave = new Double_t[npar];
		fit1->GetParameters(parsave);

		//Create one 2-d histogram for each function parameter
		Int_t ipar;
		char name[80], title[80];
		TH2D** hlist = new TH2D*[3];
		const TArrayD *xbins = residualMapTemp->GetXaxis()->GetXbins();
		const TArrayD *ybins = residualMapTemp->GetYaxis()->GetXbins();
		for (ipar=0;ipar<npar;ipar++)
		{
			snprintf(name,80,"%s_%d","ResidualMap_Fits",ipar);
			snprintf(title,80,"Fitted value of par[%d]=%s",ipar,fit1->GetParName(ipar));
			if (xbins->fN == 0)
			{
				hlist[ipar] = new TH2D(name, title,nbinsx, residualMapTemp->GetXaxis()->GetXmin(), residualMapTemp->GetXaxis()->GetXmax(),nbinsy, residualMapTemp->GetYaxis()->GetXmin(), residualMapTemp->GetYaxis()->GetXmax());
			} else {
				hlist[ipar] = new TH2D(name, title,nbinsx, xbins->fArray,nbinsy, ybins->fArray);
			}
			hlist[ipar]->GetXaxis()->SetTitle(residualMapTemp->GetXaxis()->GetTitle());
			hlist[ipar]->GetYaxis()->SetTitle(residualMapTemp->GetYaxis()->GetTitle());
		}
		snprintf(name,80,"%s_chi2","test2");
		TH2D *hchi2 = new TH2D(name,"chisquare", nbinsx, residualMapTemp->GetXaxis()->GetXmin(), residualMapTemp->GetXaxis()->GetXmax(), nbinsy, residualMapTemp->GetYaxis()->GetXmin(), residualMapTemp->GetYaxis()->GetXmax());

		//Loop on all cells in X,Y generate a projection along Z
		TH1D *hpz = new TH1D("R_temp","_temp",nbinsz, residualMapTemp->GetZaxis()->GetXmin(), residualMapTemp->GetZaxis()->GetXmax());
		Int_t bin,binx,biny,binz;
		for (biny=binminy;biny<=binmaxy;biny++)
		{
			Float_t y = residualMapTemp->GetYaxis()->GetBinCenter(biny);
			for (binx=binminx;binx<=binmaxx;binx++)
			{
				Float_t x = residualMapTemp->GetXaxis()->GetBinCenter(binx);
				hpz->Reset();
				Int_t nfill = 0;
				for (binz=1;binz<=nbinsz;binz++)
				{
					bin = residualMapTemp->GetBin(binx,biny,binz);
					Float_t w = residualMapTemp->GetBinContent(bin);
					if (w == 0)
					{
						continue;
					}
					hpz->Fill(residualMapTemp->GetZaxis()->GetBinCenter(binz),w);
					hpz->SetBinError(binz,residualMapTemp->GetBinError(bin));
					nfill++;
				}
				if (nfill < cut)
				{
					continue;
				}
				fit1->SetParameters(parsave);
				hpz->Fit(fname,option); // valgrind complains
				Int_t npfits = fit1->GetNumberFitPoints();
				if (npfits > npar && npfits >= cut)
				{
					for (ipar=0;ipar<npar;ipar++)
					{
						hlist[ipar]->Fill(x,y,fit1->GetParameter(ipar));
						hlist[ipar]->SetBinError(binx,biny,fit1->GetParError(ipar));
					}
					hchi2->SetBinContent(binx,biny,fit1->GetChisquare()/(npfits-npar));
				}
			}
		}

		delete [] parsave;
		delete hpz;
		delete fit1;
		delete hchi2;

		if (_debug <= 2)
		{
			cout << "Done residual mapping." << endl;
			cout << " " << endl;
		}

		// define cluster count as matched hits per good track
		//clustercount = clustersizehisto->GetEntries()/goodtrackcount;
		clustercount = matchedhitcount / goodtrackcount;
		clustercounterror = _clustercounterror;
		if (_debug <=0)
		{
			cout << "Matched cluster count / good track is: " << clustercount << " !" << endl;
			cout << " " << endl;
		}

		// set the current point and scale it
		double scaledcurrent = 0.0;
		double scaledcurrenterror = 0.0;
		double kboltz = 0.00008617;
		double reftemp = 253.0;
		double kelvintemp = 273.15 + _temperature.at(irun);

		// go from milli amp to micro amp:
		scaledcurrent = 1000*_sensorcurrent.at(irun)*_polarity.at(irun) * reftemp * reftemp / kelvintemp / kelvintemp * exp (-1.12 / 2 / kboltz * ((1/reftemp)-(1/kelvintemp)));

		// 10% as error:
		scaledcurrenterror = scaledcurrent * 0.1;

		double volcur = 0.0;
		double volcurerror = 0.0;

		double surfcur = 0.0;
		double surfcurerror = 0.0;

		// 1e4 for um to cm conversion
		volcur = scaledcurrent / (2.5*0.512*_thickness.at(irun)/1e4);
		surfcur = scaledcurrent / (2.5*0.512);

		if (_debug <= 0)
		{
			cout << "Scaled volume current:  " << volcur << " uA/cm^3 !" << endl;
			cout << "Scaled surface current: " << surfcur << " uA/cm^2 !" << endl;
			cout << " " << endl;
		}
		volcurerror = volcur * 0.1;
		surfcurerror = surfcur * 0.1;

		if ((_uselevel.at(irun) & (1<<4)) == 0)
		{

			tempcount = 0;
			tempcount = currentvoltage->GetN();

			currentvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),scaledcurrent);
			currentvoltage->SetPointError(tempcount,_voltageerror,scaledcurrenterror);

			tempcount = 0;
			tempcount = volumecurrent->GetN();

			volumecurrent->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),volcur);
			volumecurrent->SetPointError(tempcount,_voltageerror,volcurerror);

			tempcount = 0;
			tempcount = surfacecurrent->GetN();

			surfacecurrent->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),surfcur);
			surfacecurrent->SetPointError(tempcount,_voltageerror,surfcurerror);
		}

		if ((_uselevel.at(irun) & (1<<2)) == 0)
		{

			tempcount = 0;
			tempcount = clustersizevoltage->GetN();

			// set the cluster points
			clustersizevoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),clustersize);
			clustersizevoltage->SetPointError(tempcount,_voltageerror,clustersizeerror);

			tempcount = 0;
			tempcount = clustercountvoltage->GetN();

			if (clustercount < _cut_maxclustercount)
			{
				clustercountvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),clustercount);
				clustercountvoltage->SetPointError(tempcount,_voltageerror,clustercounterror);
			}
		}

		// calculate noise and rgh of this run

		// -2 and +1 as in matteo's analysis...
		TF1 * fivechanfit = new TF1("gausFit", "gaus", fivechannoise->GetMean() - 2 * fivechannoise->GetRMS(), fivechannoise->GetMean() + 1 * fivechannoise->GetRMS());
		fivechannoise->Fit(fivechanfit, "RQ");

		TF1 * fivechanfitelec = new TF1("gausFit", "gaus", fivechannoiseelec->GetMean() - 2 * fivechannoiseelec->GetRMS(), fivechannoiseelec->GetMean() + 1 * fivechannoiseelec->GetRMS());
		fivechannoiseelec->Fit(fivechanfitelec, "RQ");

		int channelsfornoise = 0;
		double noisecount = 0.0;
		double noiseerrorcount = 0.0;
		double sensornoise = 0.0;
		double sensornoisecenter = 0.0;
		double sensornoiseerror = 0.0;
		double sensornoiseelec = 0.0;
		double sensornoisecenterelec = 0.0;
		double sensornoiseerrorelec = 0.0;
		double globalrgh = 0.0;
		double globalnonrgh = 0.0;
		double rghratio = 0.0;
		bool rghfail = false;
		TF1 * noisefit[128];
		TF1 * noisefitelec[128];
		for (int ii=0;ii<128;ii++)
		{
			if (_channels[irun][ii] == 1)
			{
				noisefit[ii] = new TF1("noisefit","gaus",-200,200);
				noisefitelec[ii] = new TF1("noisefitelec","gaus",-30000,30000);
				noise[ii]->Fit(noisefit[ii],"QR"); // valgrind complains
				noiseelec[ii]->Fit(noisefitelec[ii],"QR"); // valgrind complains
				double tempnoise = 0.0;
				double tempnoiseerror = 0.0;
				double tempnoiseelec = 0.0;
				tempnoise = noisefit[ii]->GetParameter(2);
				tempnoiseerror = noisefit[ii]->GetParError(2);
				tempnoiseelec = noisefitelec[ii]->GetParameter(2);
				_noise[irun][ii] = tempnoise;
				allnoise->SetBinContent(ii,tempnoise);
				allnoiseelec->SetBinContent(ii,tempnoiseelec);
				if (tempnoise > 1)
				{
					channelsfornoise++;
					noisecount+=tempnoise;
					noiseerrorcount+=tempnoiseerror;
				}

				// calculate RGH ratio for the good channels
				// loop over this channel's noise and count the entries outside of 5 * sigma
				// 1000 bins
				double rgh = 0.0;
				double nonrgh = 0.0;

				// if the noise is far too high, something already went wrong in EUTel reconstruction
				// assume a high noise to estimate the (probably high) rgh ratio
				if (tempnoise >= _cut_maxrghnoise)
				{
					cout << "Fail in RGH! Too much noise in channel " << ii << " !" << endl;
					rghfail = true;
				}

				// only good channels
				// the maxnoise cut is redundant, but stays incase the above fixed high noise definition changes
				if (tempnoise < _cut_maxrghnoise)
				{
					// loop the 1k bins
					for (int jj=1; jj<1000;jj++)
					{
						// count the bin contents below center bin - 5sigma and above
						if (jj < (500.0-5.0*tempnoise) || jj > (500.0+5.0*tempnoise))
						{
							rgh += noise[ii]->GetBinContent(jj);
						} else {
							nonrgh += noise[ii]->GetBinContent(jj);
						}
					}

					globalrgh += rgh;
					globalnonrgh += nonrgh;
				}
			}
		}

		for (int ii=0;ii<128;ii++)
		{
			if (_channels[irun][ii] == 1)
			{
				delete noisefit[ii];
				delete noisefitelec[ii];
			}
		}

		TF1 * adcfit = new TF1("adcnoisefit","gaus",-200,200);
		adcnoise->Fit(adcfit,"QR"); // valgrind complains
		TF1 * adcfitelec = new TF1("adcnoisefitelec","gaus",-30000,30000);
		adcnoiseelec->Fit(adcfitelec,"QR"); // valgrind complains

		if (_debug <= 0)
		{
			cout << "RGH events: " << globalrgh << " , non-RGH events: " << globalnonrgh << " !" << endl;
			cout << " " << endl;
		}

		// catch div by zero
		if ((globalnonrgh + globalrgh) >0 && rghfail == false)
		{
			rghratio = globalrgh / (globalnonrgh + globalrgh);
		}
		if (globalnonrgh + globalrgh == 0)
		{
			cout << "FAIL in RGH estimation! No hits at all!" << endl;
			cout << " " << endl;
			rghratio = 1.0;
		}
		if (rghfail == true)
		{
			cout << "FAIL in RGH estimation! A channel had too much noise!" << endl;
			cout << " " << endl;
			rghratio = 1.0;
		}

		float rghpercent = rghratio * 100;
		
		if ((_uselevel.at(irun) & (1<<3)) == 0)
		{

			tempcount = 0;
			tempcount = rghvoltage->GetN();

			// set the point for the graph
			if (rghpercent <= _cut_maxrghplot)
			{
				rghvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),rghpercent);
				rghvoltage->SetPointError(tempcount,_voltageerror,rghpercent*0.1);
			}

			// the noise of the sensor, divided by the used channels
			sensornoise = noisecount / channelsfornoise;
			sensornoiseerror = noiseerrorcount / channelsfornoise;
			sensornoisecenter = adcfit->GetParameter(1);
			sensornoiseelec = sensornoise * _electronconversion;
			sensornoiseerrorelec = sensornoiseerror * _electronconversion;
			sensornoisecenterelec = adcfitelec->GetParameter(1);

			tempcount = 0;
			tempcount = noisevoltage->GetN();

			// cut for the histogram: only plot sensible values, drop others
			if (sensornoise > _cut_minnoise && sensornoise < _cut_maxnoise)
			{
				noisevoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),sensornoise);
				noisevoltage->SetPointError(tempcount,_voltageerror,sensornoiseerror);
				noisevoltageelec->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),sensornoiseelec);
				noisevoltageelec->SetPointError(tempcount,_voltageerror,sensornoiseerrorelec);
			} else {
				cout << "Failed cut in max sensor noise! The noise is " << sensornoise << " ADCs, but the cut is " << _cut_maxnoise << " ADCs!" << endl;
				cout << " " << endl;
			}

		}

		delete adcfit;
		delete adcfitelec;

		// done noise and rgh of this run

		if (_debug <= 2)
		{
			cout << "Starting track selection loop on " << _goodtracks.size() << " good events..." << endl;
		}
		runtime.Start(true);

		// limit to one track per event, find the one which gives more charge
		// can help with the 3GeV runs and kills some 0adc counts

		std::vector<int> selectedtracks;
		int tempeventnr = -1;
		float highsignal = 0.0;
		int trackselection = 0;
		bool singletrackevent = true;

		// loop on the good tracks
		for (unsigned int j=0;j<_goodtracks.size();j++)
		{

			float tempsignal = 0.0;
			ttel->GetEntry(_goodtracks.at(j));
			float tempfloat = dutTrackY_pixel + 0.5;
			if (xsensitive == true)
			{
				tempfloat = dutTrackX_pixel + 0.5;
			}
			int pixelpos = static_cast <int>(tempfloat);
			float eventtemp = alibava_temp;

			if (_temperature.at(irun) < 0 && _cut_applytempcorrection == true && RunNr < 4000)
			{
				float targetgain = 0.004700854;
				tempcorscale = targetgain / (0.00429259 + 7.9121e-5 * eventtemp - 3.82946e-6 * eventtemp*eventtemp);
			}
			if (_cut_applytempcorrection == true && RunNr > 4000)
			{
			     tempcorscale = 1 / (-0.0074953853 * eventtemp + 1.0763859491);
			}

			tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-2];
			tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-1];
			tempsignal += tempcorscale*alibava_reco_ch_[pixelpos];
			tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+1];
			tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+2];

			// only works for j!= 0
			if (j>0)
			{
				// if this track is from the same (aka previous) event
				if (tempeventnr == EvtNr)
				{
					singletrackevent = false;

					// trackselection is a shift to select the track with highest signal
					if (tempsignal*_polarity.at(irun) > highsignal*_polarity.at(irun))
					{
						trackselection = 0;
						highsignal = tempsignal;
					} else {
						trackselection++;
					}

				// this is a new event
				} else {

					if (_cut_onlysingletrackevents == false)
					{
						// push the previous selected track into the vector
						selectedtracks.push_back(_goodtracks.at(j-1-trackselection));
						trackselection = 0;
					}

					if (_cut_onlysingletrackevents == true && singletrackevent == true)
					{
						// push the previous selected track into the vector
						selectedtracks.push_back(_goodtracks.at(j-1-trackselection));
						trackselection = 0;
					}
					singletrackevent = true;
					highsignal = tempsignal;
				}
			}

			// set this to compare the next iteration against
			tempeventnr = EvtNr;

		} // done good track loop
		runtime.Stop();

		if (_debug <= 2)
		{
			cout << "Done track selection loop after " << runtime.RealTime() << " s with " << selectedtracks.size() << " selected good tracks!"<< endl;
			cout << " " << endl;
			cout << "Starting good track loop..." << endl;
		}
		runtime.Start(true);

		int goodtracksintimecut = 0;

		float signalL2 = 0.0;
		float signalL1 = 0.0;
		float signalC = 0.0;
		float signalR1 = 0.0;
		float signalR2 = 0.0;

		// counts where we dropped a track
		int pixelcheck = 0;
		int tdccheck = 0;
		int highchancheck = 0;

		// loop over selected good tracks
		for (unsigned int j=0;j<selectedtracks.size();j++)
		{

			ttel->GetEntry(selectedtracks.at(j));

			float tempfloat = dutTrackY_pixel + 0.5;
			if (xsensitive == true)
			{
				tempfloat = dutTrackX_pixel + 0.5;
			}
			int pixelpos = static_cast <int>(tempfloat);
			float tempsignal = 0.0;
			int pixelbin = static_cast <int>(dutTrackY_pixel);
			if (xsensitive == true)
			{
				pixelbin = static_cast <int>(dutTrackX_pixel);
			}

			float eventtemp = 0.0;
			eventtemp = alibava_temp;

			if (_temperature.at(irun) < 0 && _cut_applytempcorrection == true && RunNr < 4000)
			{
				float targetgain = 0.004700854;
				tempcorscale = targetgain / (0.00429259 + 7.9121e-5 * eventtemp - 3.82946e-6 * eventtemp*eventtemp);
			}
			if (_cut_applytempcorrection == true && RunNr > 4000)
			{
			     tempcorscale = 1 / (-0.0074953853 * eventtemp + 1.0763859491);
			}

			// add the signal of the hit channels
			if (pixelpos <=125 && pixelpos >= 2)
			{

				if (alibava_tdc >= _sensormintdc.at(irun) && alibava_tdc <= _sensormaxtdc.at(irun))
				{

					// cut: the strip with the highest pulse height is one of the center 3:
					int highchannel = -999;
					float highpulseheight = 0.0;
					for (int jj = -9;jj<=9;jj++)
					{
						// only consider good channels for this:
						//cout << "no seed found in evt " << EvtNr << endl; 
						if ( _channels[irun][pixelpos+jj] == 1 )
						{
							if ( ( tempcorscale*alibava_reco_ch_[pixelpos+jj] * _polarity.at(irun) ) > highpulseheight )
							{
								highpulseheight = tempcorscale*alibava_reco_ch_[pixelpos+jj] * _polarity.at(irun);
								highchannel = jj;
							} else {
							    //cout << " failed w/ sig " << tempcorscale*alibava_reco_ch_[pixelpos+jj] * _polarity.at(irun) << endl;
							}
						}
					}
					highchanneldistribution->Fill(highchannel);

					if (_cut_highchannelmiddle == false)
					{
						highchannel = 0;
					}

					if (highchannel == -1 || highchannel == 0 || highchannel == 1)
					{

						tempsignal = 0.0;
						// add the signal of the hit channels
						tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-2];
						tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-1];
						tempsignal += tempcorscale*alibava_reco_ch_[pixelpos];
						tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+1];
						tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+2];

						// positivise
						tempsignal = tempsignal*_polarity.at(irun);

						// now we can fill the positivised signal
						trackSignal->Fill(tempsignal);
						trackSignalMap->Fill(dutTrackY_global,tempsignal);
						trackSignalelec->Fill(tempsignal*_electronconversion);
						trackSignalMapelec->Fill(dutTrackY_global,tempsignal*_electronconversion);

						goodtracksintimecut++;

						signalGoodEvents->Fill(tempsignal);
						signalGoodEventselec->Fill(tempsignal*_electronconversion);
						
						signalMapTemp->Fill(dutTrackX_global,dutTrackY_global,tempsignal);

						if (_doxposcheck)
						{
							int ki=0;
							for (double kk=-10;kk<(20-_xevalstep);kk=kk+0.25)
							{
								if (dutTrackX_global > kk && dutTrackX_global <= (kk+_xevalstep))
								{
									posxeval[ki]->Fill(tempsignal);
								}
								ki++;
							}
						}

						signalL2 += tempcorscale*alibava_reco_ch_[pixelpos-2]*_polarity.at(irun);
						signalLeft2->Fill(tempcorscale*alibava_reco_ch_[pixelpos-2]*_polarity.at(irun));
						signalLeft2elec->Fill(tempcorscale*alibava_reco_ch_[pixelpos-2]*_polarity.at(irun)*_electronconversion);

						signalL1 += tempcorscale*alibava_reco_ch_[pixelpos-1]*_polarity.at(irun);
						signalLeft1->Fill(tempcorscale*alibava_reco_ch_[pixelpos-1]*_polarity.at(irun));
						signalLeft1elec->Fill(tempcorscale*alibava_reco_ch_[pixelpos-1]*_polarity.at(irun)*_electronconversion);

						signalC += tempcorscale*alibava_reco_ch_[pixelpos]*_polarity.at(irun);
						signalCenter->Fill(tempcorscale*alibava_reco_ch_[pixelpos]*_polarity.at(irun));
						signalCenterelec->Fill(tempcorscale*alibava_reco_ch_[pixelpos]*_polarity.at(irun)*_electronconversion);

						signalR1 += tempcorscale*alibava_reco_ch_[pixelpos+1]*_polarity.at(irun);
						signalRight1->Fill(tempcorscale*alibava_reco_ch_[pixelpos+1]*_polarity.at(irun));
						signalRight1elec->Fill(tempcorscale*alibava_reco_ch_[pixelpos+1]*_polarity.at(irun)*_electronconversion);

						signalR2 += tempcorscale*alibava_reco_ch_[pixelpos+2]*_polarity.at(irun);
						signalRight2->Fill(tempcorscale*alibava_reco_ch_[pixelpos+2]*_polarity.at(irun));
						signalRight2elec->Fill(tempcorscale*alibava_reco_ch_[pixelpos+2]*_polarity.at(irun)*_electronconversion);

						fractparty = modf (dutTrackY_pixel, &intparty);
						unmatchedEta->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]));
						unmatchedEta2D->Fill(tempcorscale*alibava_reco_ch_[pixelbin+1]/(tempcorscale*alibava_reco_ch_[pixelbin]+tempcorscale*alibava_reco_ch_[pixelbin+1]),fractparty);

						int temppixelbin = static_cast <int>(dutTrackY_pixel + 0.5);
						EtaR50->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel + 0.4);
						EtaR40->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel + 0.3);
						EtaR30->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel + 0.2);
						EtaR20->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel + 0.1);
						EtaR10->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel - 0.1);
						EtaL10->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel - 0.2);
						EtaL20->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel - 0.3);
						EtaL30->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel - 0.4);
						EtaL40->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));
						temppixelbin = static_cast <int>(dutTrackY_pixel - 0.5);
						EtaL50->Fill(tempcorscale*alibava_reco_ch_[temppixelbin+1]/(tempcorscale*alibava_reco_ch_[temppixelbin]+tempcorscale*alibava_reco_ch_[temppixelbin+1]));

						// divide the interstrip region into 8 parts and fill the signal accordingly
						if((fractparty>=0.0 && fractparty<0.125) || (fractparty >= 0.875 && fractparty < 1.0))
						{
							signalmapA->Fill(tempsignal);
							signalmapAelec->Fill(tempsignal*_electronconversion);
						}
						if((fractparty>=0.125 && fractparty<0.25) || (fractparty >= 0.75 && fractparty < 0.875))
						{
							signalmapB->Fill(tempsignal);
							signalmapBelec->Fill(tempsignal*_electronconversion);
						}
						if((fractparty>=0.25 && fractparty<0.375) || (fractparty >= 0.625 && fractparty < 0.75))
						{
							signalmapC->Fill(tempsignal);
							signalmapCelec->Fill(tempsignal*_electronconversion);
						}
						if(fractparty>=0.375 && fractparty<0.625)
						{
							signalmapD->Fill(tempsignal);
							signalmapDelec->Fill(tempsignal*_electronconversion);
						}

						// charge sharing map:
						// identify the shared channels:
						bool c = false;
						bool l1 = false;
						bool l2 = false;
						bool l3 = false;
						bool r1 = false;
						bool r2 = false;
						bool r3 = false;

						// if we have charge in the pointed-to chans:
						if (tempcorscale*alibava_reco_ch_[pixelpos]*_polarity.at(irun) > _noise[irun][pixelpos]*_cut_chargesharing)
						{
							c = true;
						}
						if (tempcorscale*alibava_reco_ch_[pixelpos-1]*_polarity.at(irun) > _noise[irun][pixelpos-1]*_cut_chargesharing)
						{
							l1 = true;
						}
						if (tempcorscale*alibava_reco_ch_[pixelpos-2]*_polarity.at(irun) > _noise[irun][pixelpos-2]*_cut_chargesharing)
						{
							l2 = true;
						}
						if (tempcorscale*alibava_reco_ch_[pixelpos-3]*_polarity.at(irun) > _noise[irun][pixelpos-3]*_cut_chargesharing)
						{
							l3 = true;
						}
						if (tempcorscale*alibava_reco_ch_[pixelpos+1]*_polarity.at(irun) > _noise[irun][pixelpos+1]*_cut_chargesharing)
						{
							r1 = true;
						}
						if (tempcorscale*alibava_reco_ch_[pixelpos+2]*_polarity.at(irun) > _noise[irun][pixelpos+2]*_cut_chargesharing)
						{
							r2 = true;
						}
						if (tempcorscale*alibava_reco_ch_[pixelpos+3]*_polarity.at(irun) > _noise[irun][pixelpos+3]*_cut_chargesharing)
						{
							r3 = true;
						}

						// fill maps
						if (c==true)
						{
							// 1:
							if (l1==false && r1==false)
							{
								hitmapChargeShared1->Fill(dutTrackX_global,dutTrackY_global);
							}

							// 2:
							if (l1==false && r1==true && r2==false)
							{
								hitmapChargeShared2->Fill(dutTrackX_global,dutTrackY_global);
							}
							if (l2==false && l1==true && r1==false)
							{
								hitmapChargeShared2->Fill(dutTrackX_global,dutTrackY_global);
							}

							// 3:
							if (l2==false && l1==true && r1==true && r2 == false)
							{
								hitmapChargeShared3->Fill(dutTrackX_global,dutTrackY_global);
							}
							if (l3==false && l2==true && l1==true && r1 == false)
							{
								hitmapChargeShared3->Fill(dutTrackX_global,dutTrackY_global);
							}
							if (l1==false && r1==true && r2==true && r3 == false)
							{
								hitmapChargeShared3->Fill(dutTrackX_global,dutTrackY_global);
							}

							// 4:
							if (l3==true && l2==true && l1==true && r1==false)
							{
								hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
							}
							if (l3==false && l2==true && l1==true && r1 == true && r2==false)
							{
								hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
							}
							if (l2==false && l1==true && r1==true && r2 == true && r3==false)
							{
								hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
							}
							if (l1==false && r1==true && r2==true && r3 == true)
							{
								hitmapChargeShared4->Fill(dutTrackX_global,dutTrackY_global);
							}

						} else {
							hitmapChargeShared0->Fill(dutTrackX_global,dutTrackY_global);
						} // done charge in pointed to channel

						highchannel_allow->Fill(1.0);

					} else {

						highchannel_discard->Fill(1.0);
						highchancheck++;

					}

					timecut_allow->Fill(1.0);

				} else { // done timecut

					timecut_discard->Fill(1.0);
					tdccheck++;

				}

				// this kills runtime!

				if (alibava_tdc>0.1)
				{

					tempsignal = 0.0;
					tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-2];
					tempsignal += tempcorscale*alibava_reco_ch_[pixelpos-1];
					tempsignal += tempcorscale*alibava_reco_ch_[pixelpos];
					tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+1];
					tempsignal += tempcorscale*alibava_reco_ch_[pixelpos+2];

					// positivise
					tempsignal = tempsignal*_polarity.at(irun);
					trackSignalTDC->Fill(alibava_tdc,tempsignal);
					trackSignalTDCelec->Fill(alibava_tdc,tempsignal*_electronconversion);

					if (_dotdccheck == true)
					{
						for (int itemp = 0; itemp < 91; itemp++)
						{
							if (alibava_tdc >= itemp && alibava_tdc <= (itemp+10))
							{
								int highchannel = -999;
								double highpulseheight = 0.0;
								for (int jj = -9;jj<=9;jj++)
								{
									if ( _channels[irun][pixelpos+jj] == 1 )
									{
										if ( ( tempcorscale*alibava_reco_ch_[pixelpos+jj] * _polarity.at(irun) ) > highpulseheight )
										{
											highpulseheight = tempcorscale*alibava_reco_ch_[pixelpos+jj] * _polarity.at(irun);
											highchannel = jj;
										}
									}
								}
								if (_cut_highchannelmiddle == false)
								{
									highchannel = 0;
								}
								if (highchannel == -1 || highchannel == 0 || highchannel == 1)
								{
									tdceval[itemp]->Fill(tempsignal);
									//tdcevalelec[itemp]->Fill(tempsignal*_electronconversion); //FIXME 
								}
							}
						}
					}

				} // 0.1 drops some entries to save time...

			} else { // pixelpos safety check
				pixelcheck++;
			}

		}
		runtime.Stop();
		if (_debug <= 2)
		{
			cout << "Done good track loop after " << runtime.RealTime() << " s!" << endl;
			cout << endl;
		}

		if (_dotdccheck == true)
		{
			double tdcevalsignal[91] = {0.0};
			double tdcevalsignalerror[91] = {0.0};
			int highbin = -1;
			double hightdc = -1;
			if (_debug <= 2)
			{
				cout << "Doing TDC checking!" << endl;
				cout << endl;
			}
			for (int itemp = 0; itemp < 91; itemp++)
			{
				if (tdceval[itemp]->GetEntries() > 100)
				{
					TF1 *tempfit = gausLanGausFitFixGausNoise(tdceval[itemp],1.0,7.0,fivechanfit->GetParameter(1),fivechanfit->GetParameter(2));
					//TF1 *tempfit = gausLanGausFitFixGausNoise(tdceval[itemp],1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
					tdcevalsignal[itemp] = tempfit->GetParameter(4);
					tdcevalsignalerror[itemp] = tempfit->GetParError(4);
					tempcount = 0;
					tempcount = overalltdc->GetN();
					overalltdc->SetPoint(tempcount,itemp+5,tdcevalsignal[itemp]);
					overalltdc->SetPointError(tempcount,5.0,tdcevalsignalerror[itemp]);
					if (_debug <= 0)
					{
						cout << "TDC bin " << itemp << " - " << itemp+10 << ": Signal " << tdcevalsignal[itemp]<< " ADCs!" << endl;
					}
					if (tdcevalsignal[itemp] > hightdc)
					{
						hightdc = tdcevalsignal[itemp];
						highbin = itemp;
					}
					delete tempfit;
				}
			}
			if (_debug <= 2)
			{
				cout << "Highest signal is from TDC "<< highbin << " to " << highbin+10 << " with signal " << hightdc << " ADCs!" << endl;
				cout << endl;
			}
		}

		if (_doxposcheck == true)
		{
			double posxevalsignal[100] = {0.0};
			double posxevalsignalerror[100] = {0.0};
			int highbin = -1;
			double highposx = -1;
			if (_debug <= 2)
			{
				cout << "Doing x position checking!" << endl;
				cout << endl;
			}
			for (int itemp = 0; itemp < 100; itemp++)
			{
				if (posxeval[itemp]->GetEntries() > 50)
				{
					TF1 *tempfit = gausLanGausFitFixGausNoise(posxeval[itemp],1.0,7.0,fivechanfit->GetParameter(1),fivechanfit->GetParameter(2));
					posxevalsignal[itemp] = tempfit->GetParameter(4);
					posxevalsignalerror[itemp] = tempfit->GetParError(4);
					tempcount = 0;
					tempcount = overallxpos->GetN();
					overallxpos->SetPoint(tempcount,(-10.0+(itemp*0.25)+(_xevalstep/2.0)),posxevalsignal[itemp]);
					overallxpos->SetPointError(tempcount,(_xevalstep/2.0),posxevalsignalerror[itemp]);
					if (_debug <= 0)
					{
						cout << "posx bin " << -10.0+(itemp*0.25) << " - " << -10.0+(itemp*0.25)+_xevalstep << ": Signal " << posxevalsignal[itemp]<< " ADCs!" << endl;
					}
					if (posxevalsignal[itemp] > highposx)
					{
						highposx = posxevalsignal[itemp];
						highbin = itemp;
					}
					delete tempfit;
				}
			}
			if (_debug <= 2)
			{
				cout << "Highest signal is from posx "<< -10.0+(highbin*0.25) << " to " << -10.0+(highbin*0.25)+_xevalstep << " with signal " << highposx << " ADCs!" << endl;
				cout << endl;
			}
		}

		// the amount of charge sharing
		float chargeshared = 0.0;
		float tracksshared = 0.0;
		tracksshared += hitmapChargeShared4->GetEntries();
		tracksshared += hitmapChargeShared3->GetEntries();
		tracksshared += hitmapChargeShared2->GetEntries();

		if (_debug <1)
		{
			cout << "Total used tracks:              " << goodtracksintimecut << " !" << endl;
			cout << "Tracks failing pixel check:     " << pixelcheck << endl;
			cout << "Tracks failing tdc check:       " << tdccheck << endl;
			cout << "Tracks failing high chan check: " << highchancheck << endl;
			cout << " " << endl;
		}

		if ((_uselevel.at(irun) & (1<<1)) == 0)
		{

			chargeshared = tracksshared / goodtracksintimecut * 100;
			tempcount = 0;
			tempcount = chargesharingvoltage->GetN();
			chargesharingvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),chargeshared);
			chargesharingvoltage->SetPointError(tempcount,_voltageerror,0.1*chargeshared);

			if (_debug <1)
			{
				cout << "Charge sharing, track method: " << chargeshared << " %!" << endl;
				cout << " " << endl;
			}
		}

		_goodtracks.clear();
		selectedtracks.clear();

		// integrate the eta distributions
		double counts = 0.0;
		double integral = 0.0;
		for (int ii = 1; ii<matchedEta->GetNbinsX(); ii++)
		{
			counts = matchedEta->GetBinContent(ii);
			integral += counts ;
			matchedEtaIntegral->SetBinContent(ii,integral);
		}

		counts = 0.0;
		integral = 0.0;
		for (int ii = 1; ii<striphitEta->GetNbinsX(); ii++)
		{
			counts = striphitEta->GetBinContent(ii);
			integral += counts ;
			striphitEtaIntegral->SetBinContent(ii,integral);
		}

		counts = 0.0;
		integral = 0.0;
		for (int ii = 1; ii<unmatchedEta->GetNbinsX(); ii++)
		{
			counts = unmatchedEta->GetBinContent(ii);
			integral += counts ;
			unmatchedEtaIntegral->SetBinContent(ii,integral);
		}

		// since we have calculated the eta distribution of noise, we can subtract it from unmatchedeta...

		// scale to -0.5 and +1.5, average and subtract
		float leftnoise = noiseEta->GetBinContent(20);
		float rightnoise = noiseEta->GetBinContent(100);
		float avgetanoise = (leftnoise+rightnoise)/2.0;
		float etaleftsig = unmatchedEta->GetBinContent(20);
		float etarightsig = unmatchedEta->GetBinContent(100);
		float avgetasig = (etaleftsig+etarightsig)/2.0;
		float etascale = avgetasig/avgetanoise;

		for (int ii=1;ii<unmatchedEta->GetNbinsX();ii++)
		{
			float tempeta = unmatchedEta->GetBinContent(ii);
			float tempeta2 = noiseEta->GetBinContent(ii);
			if ((tempeta - etascale*tempeta2) > 0)
			{
				noisesubtractedEta->SetBinContent(ii,(tempeta - etascale*tempeta2));
			} else {
				noisesubtractedEta->SetBinContent(ii,0);
			}
		}

		// define charge sharing out of eta:
		float totaleta = unmatchedEta->GetEntries();
		// noise are the counts <0 and >1
		float etanoise = 0.0;
		float etachargeshared = 0.0;

		// which eta distribution to use? unmatched or noisesubtractedEta ?

		for (int ii = 1; ii<unmatchedEta->GetNbinsX(); ii++)
		{
			if (unmatchedEta->GetBinCenter(ii) < 0 || unmatchedEta->GetBinCenter(ii) > 1.0)
			{
				etanoise+=unmatchedEta->GetBinContent(ii);
			}
			if (unmatchedEta->GetBinCenter(ii) > _cut_eta_chargeshare && unmatchedEta->GetBinCenter(ii) < (1.0 - _cut_eta_chargeshare) )
			{
				etachargeshared+=unmatchedEta->GetBinContent(ii);
			}
		}

		/*
		for (int ii = 1; ii<noisesubtractedEta->GetNbinsX(); ii++)
		{
			if (noisesubtractedEta->GetBinCenter(ii) < 0 || noisesubtractedEta->GetBinCenter(ii) > 1.0)
			{
				etanoise+=noisesubtractedEta->GetBinContent(ii);
			}
			if (noisesubtractedEta->GetBinCenter(ii) > _cut_eta_chargeshare && noisesubtractedEta->GetBinCenter(ii) < (1.0 - _cut_eta_chargeshare) )
			{
				etachargeshared+=noisesubtractedEta->GetBinContent(ii);
			}
		}
		*/

		//cout << "total eta counts: " << totaleta << endl;
		//cout << "shared counts: " << etachargeshared << endl;
		if ((totaleta) > 0)
		{
			//etanoise = etanoise / totaleta * 100.0;
			etachargeshared = etachargeshared / (totaleta) * 100.0;
		} else {
			etanoise = 0.0;
			etachargeshared = 0.0;
		}

		if ((_uselevel.at(irun) & (1<<1)) == 0)
		{
			tempcount = 0;
			tempcount = etachargesharingvoltage->GetN();
			etachargesharingvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),etachargeshared);
			etachargesharingvoltage->SetPointError(tempcount,_voltageerror,0.1*etachargeshared);

			if (_debug < 1)
			{
				cout << "Charge sharing, eta method: " << etachargeshared << " %!" << endl;
				cout << " " << endl;
			}
		}

		// residuals of this run

		if (_debug < 1)
		{
			cout << "Calculating residuals!" << endl;
			cout << " " << endl;
		}

		TF1 * residualsXfit = new TF1("residualsXfit","gaus + [3]", -0.5, 0.5);
		TF1 * residualsYfit = new TF1("residualsYfit","gaus + [3]", -0.5, 0.5);
		residualsXfit->SetParNames("Constant", "Mean", "Sigma", "Offset");
		residualsXfit->SetLineColor(kRed);
		residualsXfit->SetParameter(0, residualsX->GetMaximum());
		residualsXfit->SetParameter(1, residualsX->GetMean());
		residualsXfit->SetParameter(2, residualsX->GetRMS() * 0.5);
		residualsXfit->SetParameter(3, 0);
		residualsXfit->SetParLimits(0, 0, residualsX->GetEntries());
		residualsXfit->SetParLimits(1, residualsX->GetMean() - residualsX->GetRMS(), residualsX->GetMean() + residualsX->GetRMS());
		residualsXfit->SetParLimits(2, 0, 2 * residualsX->GetRMS());
		residualsXfit->SetParLimits(3, 0, residualsX->GetEntries());

		residualsYfit->SetParNames("Constant", "Mean", "Sigma", "Offset");
		residualsYfit->SetLineColor(kRed);
		residualsYfit->SetParameter(0, residualsY->GetMaximum());
		residualsYfit->SetParameter(1, residualsY->GetMean());
		residualsYfit->SetParameter(2, residualsY->GetRMS() * 0.5);
		residualsYfit->SetParameter(3, 0);
		residualsYfit->SetParLimits(0, 0, residualsY->GetEntries());
		residualsYfit->SetParLimits(1, residualsY->GetMean() - residualsY->GetRMS(), residualsY->GetMean() + residualsY->GetRMS());
		residualsYfit->SetParLimits(2, 0, 2 * residualsY->GetRMS());
		residualsYfit->SetParLimits(3, 0, residualsY->GetEntries());

		residualsX->Fit(residualsXfit,"QR"); // valgrind complains
		residualsY->Fit(residualsYfit,"QR"); // valgrind complains

		// set the graph points for resi vs voltage
		double fitparx1 = 0.0;
		double fitpary1 = 0.0;
		double fitparx2 = 0.0;
		double fitpary2 = 0.0;
		fitparx1 = residualsXfit->GetParameter(2);
		fitparx2 = residualsXfit->GetParError(2);
		fitpary1 = residualsYfit->GetParameter(2);
		fitpary2 = residualsYfit->GetParError(2);

		delete residualsXfit;
		delete residualsYfit;

		tempcount = 0;
		tempcount = residualsXvoltage->GetN();
		//cout << " " << endl;
		//cout << tempcount << " points in Residuals Y graph, adding one..." << endl;
		if (fitparx1 < _cut_maxXresidual && ((_uselevel.at(irun) & (1<<0)) == 0))
		{
			residualsXvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1);
			residualsXvoltage->SetPointError(tempcount,_voltageerror,fitparx2);
			if (_debug <= 2)
			{
				cout << "Residual in X is "<< fitparx1 << " mm!" << endl;
				cout << " " << endl;
			}
		} else {
			if (_debug <= 2)
			{
				cout << "Dropping X residual because it is too high or wrong uselevel!" << endl;
				cout << " " << endl;
			}
		}

		tempcount = 0;
		tempcount = residualsYvoltage->GetN();

		if (fitpary1 < _cut_maxYresidual && ((_uselevel.at(irun) & (1<<0)) == 0))
		{
			residualsYvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitpary1);
			residualsYvoltage->SetPointError(tempcount,_voltageerror,fitpary2);
			if (_debug <= 2)
			{
				cout << "Residual in Y is "<< fitpary1 << " mm!" << endl;
				cout << " " << endl;
			}
		} else {
			if (_debug <= 2)
			{
				cout << "Dropping Y residual because it is too high or wrong uselevel!" << endl;
				cout << " " << endl;
			}
		}

		tempcount = 0;
		tempcount = residualsYvoltageAngle->GetN();
		residualsYvoltageAngle->SetPoint(tempcount,_dutrotation.at(irun),fitpary1);
		residualsYvoltageAngle->SetPointError(tempcount,2.5,fitpary2);

		// done this run's residuals

		// signal  and fit of this run:

		// the clone-histos can only be written later on if they exist -> the (if tracksignal ) has to be true
		bool skipsignals = false;

		// limit this to good track signals -> at least 100 entires
		if (trackSignal->GetEntries() > 100)
		{

			if (_debug < 1)
			{
				cout << "Calculating signal!" << endl;
				cout << " " << endl;
			}

			// a gaus to get the first noise peak
			TF1 * gausnoise = new TF1("gausnoise","gaus",-150,150.0);

			// constrain the parameters:
			// first set all to "sensible" values: the constant to the bin entry at 0 ADCs (500th bin)
			// the mean gets set to the sensor noise mean and fixed
			// the sigma is assumed to be 1.5 times the channel sigma
			double tempval = 0.0;
			tempval = trackSignal->GetBinContent(500);
			// protect against empty bin
			if (tempval < 1.0)
			{
				tempval = 1.0;
			}
			gausnoise->SetParameter(0,tempval);
			gausnoise->SetParameter(1,sensornoisecenter);
			gausnoise->SetParameter(2,sensornoise*1.5);

			// limits on the parameters:
			// if we have five channels with only noise contributions, their width will be:
			float fivechannoisetest = sqrt(5*sensornoise*sensornoise);
			gausnoise->SetParLimits(0,0,2000);
			gausnoise->FixParameter(1,sensornoisecenter);
			gausnoise->SetParLimits(2,sensornoise,fivechannoisetest);

			// fit to 5adcs:
			trackSignal->Fit(gausnoise,"RBQ","",-50.0,5.0);
			gausnoise->Draw();

			delete gausnoise;

			// fit this guy with a landau-gaus
			//TF1 *fit_goodevents = lanGausFit(signalGoodEvents,3.0,7.0);
			TF1 *fit_goodevents = gausLanGausFitFixGausNoise(signalGoodEvents,1.0,7.0,fivechanfit->GetParameter(1),fivechanfit->GetParameter(2));
			//TF1 *fit_goodevents = gausLanGausFitFixGausNoise(signalGoodEvents,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
			TF1 *fit_goodeventselec = gausLanGausFitFixGausNoise(signalGoodEventselec,1.0,7.0,fivechannoiseelec->GetMean()*_polarity.at(irun),fivechannoiseelec->GetRMS());

			// set the graph points for signal vs voltage
			double fitparx1 = 0.0;
			double fitparx2 = 0.0;
			fitparx1 = fit_goodevents->GetParameter(4);
			fitparx2 = fit_goodevents->GetParError(4);

			double fitparx1elec = 0.0;
			double fitparx2elec = 0.0;
			fitparx1elec = fit_goodeventselec->GetParameter(4);
			fitparx2elec = fit_goodeventselec->GetParError(4);

			// cut for the histogram: only plot sensible values, drop others
			if (fitparx1 > _cut_minsignal && fitparx1 < _cut_maxsignal && ((_uselevel.at(irun) & (1<<1)) == 0))
			{

				tempcount = 0;
				tempcount = signalvoltage->GetN();

				signalvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1);
				signalvoltage->SetPointError(tempcount,_voltageerror,fitparx2);

				tempcount = 0;
				tempcount = signalvoltageelec->GetN();

				signalvoltageelec->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1elec);
				signalvoltageelec->SetPointError(tempcount,_voltageerror,fitparx2elec);

				tempcount = 0;
				tempcount = signaldistancevoltage->GetN();

				signaldistancevoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)));
				signaldistancevoltage->SetPointError(tempcount,_voltageerror,fitparx2/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)));

				tempcount = 0;
				tempcount = signaldistancevoltageelec->GetN();

				signaldistancevoltageelec->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),fitparx1elec/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)));
				signaldistancevoltageelec->SetPointError(tempcount,_voltageerror,fitparx2elec/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)));
		
				//cout << "point is " << fitparx1/(_thickness.at(irun)/cos(_dutrotation.at(irun)*PI/180.0)) << endl;
				//cout << "distance is " << cos(_dutrotation.at(irun)*PI/180.0) << endl;

				// signal to noise:
				if (sensornoise != 0)
				{
					signaltonoise = fitparx1 / sensornoise;
					signaltonoiseerror = sqrt(fitparx2*fitparx2 + sensornoiseerror*sensornoiseerror);

					tempcount = 0;
					tempcount = snvoltage->GetN();

					snvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),signaltonoise);
					snvoltage->SetPointError(tempcount,_voltageerror,signaltonoiseerror);
				} else {
					cout << "Dropping S/N point: Noise is 0!" << endl;
					cout << " " << endl;
				}
			} else {
				cout << "Dropping signal, as it is out of range! ( " << fitparx1 << " ) or wrong uselevel!" << endl;
				cout << " " << endl;
			}

			// estimate how many channels (of the 5) are giving signal and how many are noise:
			// get the fit's noise sigma and compare it to the expected width of 5 noisy channels
			float effectivechan = 0.0;
			float effectivechanerror = 0.0;
			if (sensornoise != 0)
			{
				//effectivechan = 5.0 - (fit_goodevents->GetParameter(3))*(fit_goodevents->GetParameter(3))/sensornoise/sensornoise;
				effectivechan = (fit_goodevents->GetParameter(6))*(fit_goodevents->GetParameter(6)) - 5*sensornoise*sensornoise;
				if (_debug < 1)
				{
					cout << "Effective channel count is " << effectivechan << endl;
					cout << " " << endl;
				}
			} else {
				cout << "Dropping channel point! Noise is 0!" << endl;
				cout << " " << endl;
			}
			if (sensornoiseerror != 0)
			{
				effectivechanerror = sqrt(fit_goodevents->GetParError(6)*fit_goodevents->GetParError(6) +sensornoiseerror*sensornoiseerror);
			} else {
				cout << "Dropping channel point error! Noiseerror is 0!" << endl;
				cout << " " << endl;
			}

			if ((_uselevel.at(irun) & (1<<1)) == 0)
			{

				tempcount = 0;
				tempcount = channelcountvoltage->GetN();

				channelcountvoltage->SetPoint(tempcount,_biasvoltage.at(irun)*_polarity.at(irun),effectivechan);
				channelcountvoltage->SetPointError(tempcount,_voltageerror,effectivechanerror);
			}

		} else {
			if (_debug <= 2)
			{
			    cout << "Dropping signal, not enough entries!" << endl;
			    cout << " " << endl;
			    skipsignals = true;
			}
		}

		if (skipsignals == false)
		{

			//TF1 *fit_signalmapA = gausLanGausFitFixGausNoise(signalmapA,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
			//TF1 *fit_signalmapB = gausLanGausFitFixGausNoise(signalmapB,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
			//TF1 *fit_signalmapC = gausLanGausFitFixGausNoise(signalmapC,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
			//TF1 *fit_signalmapD = gausLanGausFitFixGausNoise(signalmapD,1.0,7.0,fivechannoise->GetMean()*_polarity.at(irun),fivechannoise->GetRMS());
			TF1 *fit_signalmapA = gausLanGausFitFixGausNoise(signalmapA,1.0,7.0,fivechanfit->GetParameter(1),fivechanfit->GetParameter(2));
			TF1 *fit_signalmapB = gausLanGausFitFixGausNoise(signalmapB,1.0,7.0,fivechanfit->GetParameter(1),fivechanfit->GetParameter(2));
			TF1 *fit_signalmapC = gausLanGausFitFixGausNoise(signalmapC,1.0,7.0,fivechanfit->GetParameter(1),fivechanfit->GetParameter(2));
			TF1 *fit_signalmapD = gausLanGausFitFixGausNoise(signalmapD,1.0,7.0,fivechanfit->GetParameter(1),fivechanfit->GetParameter(2));

			//TF1 *fit_signalmapAelec = gausLanGausFitFixGausNoise(signalmapAelec,1.0,7.0,fivechannoiseelec->GetMean()*_polarity.at(irun),fivechannoiseelec->GetRMS());
			//TF1 *fit_signalmapBelec = gausLanGausFitFixGausNoise(signalmapBelec,1.0,7.0,fivechannoiseelec->GetMean()*_polarity.at(irun),fivechannoiseelec->GetRMS());
			//TF1 *fit_signalmapCelec = gausLanGausFitFixGausNoise(signalmapCelec,1.0,7.0,fivechannoiseelec->GetMean()*_polarity.at(irun),fivechannoiseelec->GetRMS());
			//TF1 *fit_signalmapDelec = gausLanGausFitFixGausNoise(signalmapDelec,1.0,7.0,fivechannoiseelec->GetMean()*_polarity.at(irun),fivechannoiseelec->GetRMS());
			TF1 *fit_signalmapAelec = gausLanGausFitFixGausNoise(signalmapAelec,1.0,7.0,fivechanfitelec->GetParameter(1)*_polarity.at(irun),fivechanfitelec->GetParameter(2));
			TF1 *fit_signalmapBelec = gausLanGausFitFixGausNoise(signalmapBelec,1.0,7.0,fivechanfitelec->GetParameter(1)*_polarity.at(irun),fivechanfitelec->GetParameter(2));
			TF1 *fit_signalmapCelec = gausLanGausFitFixGausNoise(signalmapCelec,1.0,7.0,fivechanfitelec->GetParameter(1)*_polarity.at(irun),fivechanfitelec->GetParameter(2));
			TF1 *fit_signalmapDelec = gausLanGausFitFixGausNoise(signalmapDelec,1.0,7.0,fivechanfitelec->GetParameter(1)*_polarity.at(irun),fivechanfitelec->GetParameter(2));

			double fitpara = 0.0;
			double fitparb = 0.0;
			double fitparc = 0.0;
			double fitpard = 0.0;
			double fitparera = 0.0;
			double fitparerb = 0.0;
			double fitparerc = 0.0;
			double fitparerd = 0.0;

			double fitparaelec = 0.0;
			double fitparbelec = 0.0;
			double fitparcelec = 0.0;
			double fitpardelec = 0.0;
			double fitpareraelec = 0.0;
			double fitparerbelec = 0.0;
			double fitparercelec = 0.0;
			double fitparerdelec = 0.0;

			fitpara = fit_signalmapA->GetParameter(4);
			fitparb = fit_signalmapB->GetParameter(4);
			fitparc = fit_signalmapC->GetParameter(4);
			fitpard = fit_signalmapD->GetParameter(4);
			fitparera = fit_signalmapA->GetParError(4);
			fitparerb = fit_signalmapB->GetParError(4);
			fitparerc = fit_signalmapC->GetParError(4);
			fitparerd = fit_signalmapD->GetParError(4);

			fitparaelec = fit_signalmapAelec->GetParameter(4);
			fitparbelec = fit_signalmapBelec->GetParameter(4);
			fitparcelec = fit_signalmapCelec->GetParameter(4);
			fitpardelec = fit_signalmapDelec->GetParameter(4);
			fitpareraelec = fit_signalmapAelec->GetParError(4);
			fitparerbelec = fit_signalmapBelec->GetParError(4);
			fitparercelec = fit_signalmapCelec->GetParError(4);
			fitparerdelec = fit_signalmapDelec->GetParError(4);

			int thisvoltage = _biasvoltage.at(irun)*_polarity.at(irun) / 100;

			// signalareamap is the correct histo, so we can fill it now:
			signalareamap->SetBinContent(thisvoltage,1,fitpara);
			signalareamap->SetBinContent(thisvoltage,2,fitparb);
			signalareamap->SetBinContent(thisvoltage,3,fitparc);
			signalareamap->SetBinContent(thisvoltage,4,fitpard);

			int tempcount = signalareaplot->GetN();
			signalareaplot->SetPoint(tempcount,5,fitpara);
			signalareaplot->SetPointError(tempcount,4.5,fitparera);
			tempcount = signalareaplot->GetN();
			signalareaplot->SetPoint(tempcount,15,fitparb);
			signalareaplot->SetPointError(tempcount,4.5,fitparerb);
			tempcount = signalareaplot->GetN();
			signalareaplot->SetPoint(tempcount,25,fitparc);
			signalareaplot->SetPointError(tempcount,4.5,fitparerc);
			tempcount = signalareaplot->GetN();
			signalareaplot->SetPoint(tempcount,35,fitpard);
			signalareaplot->SetPointError(tempcount,4.5,fitparerd);

			tempcount = signalareaplotelec->GetN();
			signalareaplotelec->SetPoint(tempcount,5,fitparaelec);
			signalareaplotelec->SetPointError(tempcount,4.5,fitpareraelec);
			tempcount = signalareaplotelec->GetN();
			signalareaplotelec->SetPoint(tempcount,15,fitparbelec);
			signalareaplotelec->SetPointError(tempcount,4.5,fitparerbelec);
			tempcount = signalareaplotelec->GetN();
			signalareaplotelec->SetPoint(tempcount,25,fitparcelec);
			signalareaplotelec->SetPointError(tempcount,4.5,fitparercelec);
			tempcount = signalareaplotelec->GetN();
			signalareaplotelec->SetPoint(tempcount,35,fitpardelec);
			signalareaplotelec->SetPointError(tempcount,4.5,fitparerdelec);

		}
		if (_debug <= 2)
		{
			cout << "Done signal processing, writing this run to file!" << endl;
			cout << " " << endl;
		}

		// done with the signal processing

		// starting with the output of this run

		// done with the loaded ntuple file
		f1->Close();

		// doesn't work? FIXME
		// may need a root object/application loaded at the start... 
		gStyle->SetOptFit(1111); // valgrind complains

		// write all the individual histos to file
		_outputFile->cd();

		// each run gets its own directory
		int therun = _runnumber.at(irun);
		sprintf(name,"run%i",therun);
		TDirectory* rundirectory = _outputFile->mkdir(name);

		rundirectory->cd();
		TDirectory* noisedirectory = rundirectory->mkdir("noise");
		noisedirectory->cd();
		TDirectory* noisechandirectory = noisedirectory->mkdir("channels");
		noisechandirectory->cd();
		for(int ii=0;ii<128;ii++)
		{
			noise[ii]->Write();
			noiseelec[ii]->Write();
		}
		noisedirectory->cd();
		allnoise->Write();
		allnoiseelec->Write();
		adcnoise->Write();
		adcnoiseelec->Write();
		fivechannoise->Write();
		fivechannoiseelec->Write();

		rundirectory->cd();
		TDirectory* resdirectory = rundirectory->mkdir("residuals");
		resdirectory->cd();

		residualsX->Write();
		residualsY->Write();
		residualsZ->Write();
		residualsXY->Write(); // valgrind complains
		residualsXZ->Write();
		residualsYZ->Write();
		residualsYQ->Write();
		residualsYR->Write();
		residualsYT->Write();
		residualsY_Q->Write();

		TDirectory* resdirectorysub = resdirectory->mkdir("clustersizes");
		resdirectorysub->cd();

		residualsY_clu1->Write();
		residualsY_clu2->Write();
		residualsY_clu3->Write();
		residualsY_clu4->Write();

		resdirectory->cd();

		residualsDXvsX->Write();
		residualsDXvsY->Write();
		residualsDXvsZ->Write();
		residualsDYvsX->Write();
		residualsDYvsY->Write();
		residualsDYvsZ->Write();
		residualsDZvsX->Write();
		residualsDZvsY->Write();
		residualsDZvsZ->Write();

		residualsDXvsSXU->Write();
		residualsDXvsSYU->Write();
		residualsDYvsSXU->Write();
		residualsDYvsSYU->Write();
		residualsDZvsSXU->Write();
		residualsDZvsSYU->Write();

		residualsDXvsSXD->Write();
		residualsDXvsSYD->Write();
		residualsDYvsSXD->Write();
		residualsDYvsSYD->Write();
		residualsDZvsSXD->Write();
		residualsDZvsSYD->Write();

		residualsXvsEvt->Write();
		residualsYvsEvt->Write();
		residualsZvsEvt->Write();

                for (int ipar=2;ipar<npar;ipar++)
		{
			hlist[ipar]->Write();
		}

		for (int ipar=0;ipar<3;ipar++)
		{
			delete hlist[ipar];
		}
		delete[] hlist;

		residualProfileDXvsX->Write(); // valgrind complains
		residualProfileDXvsY->Write();
		residualProfileDYvsX->Write();
		residualProfileDYvsY->Write();

		TDirectory* resdirectorycuts = resdirectory->mkdir("cuts");
		resdirectorycuts->cd();

		trackselection_discard->Write();
		trackselection_allow->Write();
		fiducial_discard->Write();
		fiducial_allow->Write();
		goodchannel_discard->Write();
		goodchannel_allow->Write();
		nohittrack->Write();

		rundirectory->cd();
		TDirectory* hitmapdirectory = rundirectory->mkdir("hitmaps");
		hitmapdirectory->cd();

		hitmapDUT->Write();
		hitmapDUTmodpitch->Write();
		DUThitmap3D->Write();

		hitmapTracks->Write();
		hitmapTracksmodpitch->Write();

		hitmapMatch->Write();
		hitmapMatchmodpitch->Write();
		hitmapMatchTracks->Write();

		scatterX->Write();
		TH2D* testing = new TH2D;
		testing = scatterX->ProjectionXY("aname","C=E");
		//testing->SetMaximum(0.5);
		testing->Write();
		scatterY->Write();
		TH2D* testing2 = new TH2D;
		testing2 = scatterY->ProjectionXY("aname2","C=E");
                //testing2->SetMaximum(0.5);
		testing2->Write();
		scatterXY->Write();
		TH2D* testing3 = new TH2D;
		testing3 = scatterXY->ProjectionXY("aname3","C=E");
                //testing3->SetMaximum(0.5);
		testing3->Write();
		kinkX->Write();
		kinkY->Write();

		hitmapChargeShared0->Write();
		hitmapChargeShared1->Write();
		hitmapChargeShared2->Write();
		hitmapChargeShared3->Write();
		hitmapChargeShared4->Write();

		hitmapClusterSizeA->Write();
		hitmapClusterSize1->Write();
		hitmapClusterSize2->Write();
		hitmapClusterSize3->Write();
		hitmapClusterSize4->Write();

		signalClusterSize1->Write();
		signalClusterSize1elec->Write();
		signalClusterSize2->Write();
		signalClusterSize2elec->Write();
		signalClusterSize3->Write();
		signalClusterSize3elec->Write();
		signalClusterSize4->Write();
		signalClusterSize4elec->Write();

		rundirectory->cd();
		TDirectory* etadirectory = rundirectory->mkdir("eta");
		etadirectory->cd();

		// a canvas for the eta integrals
		if (_cut_drawetaintegralcomparison == true)
		{

			TCanvas *canv_integral = new TCanvas("Eta Integrals","transparent pad",600,400);
			gStyle->SetOptStat(kFALSE);
			unmatchedEtaIntegral->SetStats(kFALSE);
			unmatchedEtaIntegral->Draw();
			canv_integral->Update();

			Float_t rightmax = 1.05*matchedEtaIntegral->GetMaximum();
			Float_t scale = gPad->GetUymax()/rightmax;
			matchedEtaIntegral->SetLineColor(kRed);
			matchedEtaIntegral->Scale(scale);
			matchedEtaIntegral->SetStats(kFALSE);
			matchedEtaIntegral->Draw("same");

			//draw an axis on the right side
			TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
			axis->SetLineColor(kRed);
			axis->SetLabelColor(kRed);
			axis->SetTitle("Events");
			axis->Draw();

			canv_integral->Write();

			delete canv_integral;
			delete axis;

		}

		matchedEta->Write();
		matchedEta2D->Write();
		matchedEtaIntegral->Write();
		unmatchedEta->Write();
		unmatchedEta2D->Write();
		unmatchedEtaIntegral->Write();
		striphitEta->Write();
		striphitEta2D->Write();
		striphitEtaIntegral->Write();
		noiseEta->Write();
		noisesubtractedEta->Write();
		TDirectory* etashiftdirectory = etadirectory->mkdir("shifts");
		etashiftdirectory->cd();
		EtaL50->Write();
		EtaL40->Write();
		EtaL30->Write();
		EtaL20->Write();
		EtaL10->Write();
		EtaR10->Write();
		EtaR20->Write();
		EtaR30->Write();
		EtaR40->Write();
		EtaR50->Write();

		rundirectory->cd();
		TDirectory* trackdirectory = rundirectory->mkdir("tracks");
		trackdirectory->cd();

		trackSignalMap->Write();
		trackSignalMapelec->Write();
		trackSignalTDC->Write();
		trackSignalTDCelec->Write();
		tracksperevent->Write();
		tracksvsevents->Write();
		highchanneldistribution->Write();

		trackdirectory->cd();
		
		if (_doxposcheck == true)
		{

			TCanvas* canvas_xpos = new TCanvas("XposEval","X Position Evaluation",200,10,700,500);
			canvas_xpos->cd();
			TH2F *hxpos = new TH2F("XposEval","X Position Evaluation",40,-20,20,40,0,100);
			hxpos->SetXTitle("X Position [mm]");
			hxpos->SetYTitle("Track Signal [ADCs]");
			hxpos->SetStats(0000);
			hxpos->Draw();
			overallxpos->Draw("LP");
			canvas_xpos->Update();
			canvas_xpos->Write();
			canvas_xpos->Close();
			delete canvas_xpos;
			delete hxpos;
			delete overallxpos;

			TDirectory* posxevaldir = trackdirectory->mkdir("posxeval");
			posxevaldir->cd();
			for(int ii=0;ii<100;ii++)
			{
				posxeval[ii]->Write();
			}
		}
		trackdirectory->cd();

		signalLeft2->Write();
		signalLeft2elec->Write();
		signalLeft1->Write();
		signalLeft1elec->Write();
		signalCenter->Write();
		signalCenterelec->Write();
		signalRight1->Write();
		signalRight1elec->Write();
		signalRight2->Write();
		signalRight2elec->Write();
		signalGoodEvents->Write();
		signalGoodEventselec->Write();
		signalmapA->Write();
		signalmapAelec->Write();
		signalmapB->Write();
		signalmapBelec->Write();
		signalmapC->Write();
		signalmapCelec->Write();
		signalmapD->Write();
		signalmapDelec->Write();

		if (_dotdccheck == true)
		{
			TCanvas* canvas_tdc = new TCanvas("TDCEval","TDC Evaluation",200,10,700,500);
			canvas_tdc->cd();
			TH2F *htdc = new TH2F("TDCEval","TDC Evaluation",20,0,100,40,0,200);
			htdc->SetXTitle("TDC Range [ns]");
			htdc->SetYTitle("Track Signal [ADCs]");
			htdc->SetStats(0000);
			htdc->Draw();
			overalltdc->Draw("LP");
			canvas_tdc->Update();
			canvas_tdc->Write();
			canvas_tdc->Close();
			delete canvas_tdc;
			delete htdc;
			delete overalltdc;

			TDirectory* tdcevaldir = trackdirectory->mkdir("tdceval");
			tdcevaldir->cd();
			for(int ii=0;ii<91;ii++)
			{
				tdceval[ii]->Write();
			}
		}

		trackdirectory->cd();
		TDirectory* trackcutdirectory = trackdirectory->mkdir("cuts");
		trackcutdirectory->cd();

		goodevent_discard->Write();
		goodevent_allow->Write();
		timecut_discard->Write();
		timecut_allow->Write();
		highchannel_discard->Write();
		highchannel_allow->Write();

		rundirectory->cd();
		TDirectory* clusterdirectory = rundirectory->mkdir("clusters");
		clusterdirectory->cd();

		clustersizehisto->Write();
		clusteretadistribution->Write();
		clustersignalhisto->Write();

		rundirectory->cd();
		tdcdistri->Write();
		tempdistri->Write();

		// done with cluster file
		f2->Close();

		// inc the counter for runs to be plotted vs voltage
		// note: residual, noise and signal have cuts, so this number may be too high!
		voltagepoint[vsvoltage_histonr]++;

	} // done run loop

	if (_debug <= 3)
	{
		cout << " " << endl;
		cout << "**********************************************" << endl;
		cout << " " << endl;
		cout << "Done run loop!" << endl;
	}

	// The multi-run things

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Signal Area Mapping..." << endl;
	}

	// signal area map
	_outputFile->cd();
	TDirectory* signalareamapdir = _outputFile->mkdir("Signal vs. Area Mapping");
	signalareamapdir->cd();
	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalareamap_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TH2D* signalareamap = dynamic_cast<TH2D*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);

		signalareamap->SetXTitle("Bias Voltage [|V|]");
		signalareamap->SetYTitle("Track mod Pitch [Pitch]");
		signalareamap->SetTitle(title);

		
		signalareamap->Write();

	}

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Signal Area Plot Comparison..." << endl;
	}

	// draw all the signal area plots
	// do this for each sensor and each fluence
	// get the inner loop count (aka voltage points) from clustersize, as this is uncut and comes directly from the ntuple

	// the position counter in the runlist, since we will loop over all from runcount
	int pointcounter = 0;

	_outputFile->cd();
	TDirectory* signalareaplotdir = _outputFile->mkdir("Signal vs. Area Plots");
	signalareaplotdir->cd();

	// the count on sensors, fluences and rotations, so the number of histograms
	for (unsigned int jj=0;jj<voltagegraphs.size();jj++)
	{

		sstream.str(string());
		histoname = "clustersizevoltage_";
		sstream << histoname << voltagegraphs.at(jj);
		TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		int looplimit = 0;
		looplimit = clustersizevoltage->GetN();

		sprintf(name, "Signal Area Plots_%i", jj);
		TCanvas* canv_signalareaplot = new TCanvas(name,"transparent pad",600,400);
		char thename [100];
		sprintf(thename, "hsignalareaplot_%i", jj);
		TH2F* hsignalareaplot = new TH2F(thename,"Signal vs. Track Distance from Strip Centre",20,0,40,20,0,100);
		hsignalareaplot->SetXTitle("Distance from Strip Centre [#mum]");
		hsignalareaplot->SetYTitle("Landau MPV [ADCs]");
		hsignalareaplot->SetStats(0000);
		hsignalareaplot->Draw();

		TLegend* lsignalareaplot = new TLegend(0.59,0.55,0.90,0.85);
		lsignalareaplot->SetBorderSize(1);
		lsignalareaplot->SetFillColor(0);
		lsignalareaplot->SetFillStyle(0);

		// for counting subruns
		int tempint = 0;
		for (int i=pointcounter;i<pointcounter+looplimit;i++)
		{

			sstream.str(string());
			histoname = "signalareaplot_";
			sstream << histoname << _runnumber.at(i);
			TGraphErrors* signalareaplot = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);

			gStyle->SetOptStat(kFALSE); // valgrind complains

			signalareaplot->SetLineColor(tempint+1);
			signalareaplot->SetLineWidth(3);

			// jump yellow
			if ( tempint >= 4 )
			{
				signalareaplot->SetLineColor(tempint+2);
			}

			char o [100];
			char title[300];
			char title2[300];
			sprintf(o,"Epi");
			sstream.str(string());
			string k = _sensorname.at(i);
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			double fl = _irradfluence.at(i);
			int m = _biasvoltage.at(i);
			int n = _dutrotation.at(i);

			sprintf(title2, "%s, F=%.1e, %irot", sensname.c_str(),fl,n);
			if (tempint==0)
			{
				signalareaplot->Draw("P");
				signalareaplot->SetTitle(title2);

			} else {
				signalareaplot->Draw("P");
			}

			canv_signalareaplot->Update();

			sprintf(title, "%s, F=%.1e, %iV, %irot", sensname.c_str(),fl,m,n);
			lsignalareaplot->AddEntry(signalareaplot,title,"lep");

			tempint++;

		} // done subrun loop

		lsignalareaplot->Draw();
		canv_signalareaplot->Write(); // valgrind complains
		canv_signalareaplot->Close();
		delete lsignalareaplot;
		delete hsignalareaplot;
		delete canv_signalareaplot;

		pointcounter=pointcounter+looplimit;

	} // done voltagegraph loop
	pointcounter = 0;

	for (unsigned int jj=0;jj<voltagegraphs.size();jj++)
	{

		sstream.str(string());
		histoname = "clustersizevoltage_";
		sstream << histoname << voltagegraphs.at(jj);
		TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		int looplimit = 0;
		looplimit = clustersizevoltage->GetN();

		sprintf(name, "Signal Area Plots Electrons_%i", jj);
		TCanvas* canv_signalareaplotelec = new TCanvas(name,"transparent pad",600,400);
		char thename [100];
		sprintf(thename, "hsignalareaplotelec_%i", jj);
		TH2F* hsignalareaplotelec = new TH2F(thename,"Signal vs. Track Distance from Strip Centre, Electrons",20,0,40,20,0,15000);
		hsignalareaplotelec->SetXTitle("Distance from Strip Centre [#mum]");
		hsignalareaplotelec->SetYTitle("Landau MPV [Electrons]");
		hsignalareaplotelec->SetStats(0000);
		hsignalareaplotelec->Draw();

		TLegend* lsignalareaplotelec = new TLegend(0.59,0.55,0.90,0.85);
		lsignalareaplotelec->SetBorderSize(1);
		lsignalareaplotelec->SetFillColor(0);
		lsignalareaplotelec->SetFillStyle(0);

		// for counting subruns
		int tempint = 0;
		for (int i=pointcounter;i<pointcounter+looplimit;i++)
		{

			sstream.str(string());
			histoname = "signalareaplotelec_";
			sstream << histoname << _runnumber.at(i);
			TGraphErrors* signalareaplotelec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);

			gStyle->SetOptStat(kFALSE); // valgrind complains

			signalareaplotelec->SetLineColor(tempint+1);
			signalareaplotelec->SetLineWidth(3);

			// jump yellow
			if ( tempint >= 4 )
			{
				signalareaplotelec->SetLineColor(tempint+2);
			}

			char o [100];
			char title[300];
			char title2[300];
			sprintf(o,"Epi");
			sstream.str(string());
			string k = _sensorname.at(i);
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			double fl = _irradfluence.at(i);
			int m = _biasvoltage.at(i);
			int n = _dutrotation.at(i);

			sprintf(title2, "%s, F=%.1e, %irot", sensname.c_str(),fl,n);
			if (tempint==0)
			{
				signalareaplotelec->Draw("P");
				signalareaplotelec->SetTitle(title2);

			} else {
				signalareaplotelec->Draw("P");
			}

			canv_signalareaplotelec->Update();

			sprintf(title, "%s, F=%.1e, %iV, %irot", sensname.c_str(),fl,m,n);
			lsignalareaplotelec->AddEntry(signalareaplotelec,title,"lep");

			tempint++;

		} // done subrun loop

		lsignalareaplotelec->Draw();
		canv_signalareaplotelec->Write(); // valgrind complains
		canv_signalareaplotelec->Close();
		delete lsignalareaplotelec;
		delete hsignalareaplotelec;
		delete canv_signalareaplotelec;

		pointcounter=pointcounter+looplimit;

	} // done voltagegraph loop

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Eta Comparison..." << endl;
	}

	// draw all the eta distributions
	// do this for each sensor and each fluence
	// get the inner loop count (aka voltage points) from clustersize, as this is uncut and comes directly from the ntuple

	// the position counter in the runlist, since we will loop over all from runcount
	pointcounter = 0;

	std::vector<TLegend*> leta;

	_outputFile->cd();
	TDirectory* etacomparisondir = _outputFile->mkdir("Eta Comparisons");
	etacomparisondir->cd();

	// the count on sensors, fluences and rotations, so the number of histograms
	for (unsigned int jj=0;jj<voltagegraphs.size();jj++)
	{

		sstream.str(string());
		histoname = "clustersizevoltage_";
		sstream << histoname << voltagegraphs.at(jj);
		TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());

		int looplimit = 0;
		looplimit = clustersizevoltage->GetN();

		sprintf(name, "Eta Distributions_%i", jj);
		TCanvas* canv_eta = new TCanvas(name,"transparent pad",600,400);

		// get the highest integral of this canvas for scale
		// get the ymax of the normal distri for axis range
		Double_t highestentrycount = 0;
		Double_t yscale = 0;

		for (int i=pointcounter;i<pointcounter+looplimit;i++)
		{
			sstream.str(string());
			histoname = "unmatchedEtaIntegral_";
			sstream << histoname << _runnumber.at(i);
			TH1D* thiseta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
			Double_t temp = thiseta->GetMaximum();
			if (temp>highestentrycount)
			{
				highestentrycount = temp;
			}

			sstream.str(string());
			histoname = "unmatchedEta_";
			sstream << histoname << _runnumber.at(i);
			TH1D* thiseta2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);
			Double_t temp2 = thiseta2->GetMaximum();
			if (temp2>yscale)
			{
				yscale = temp2;
			}

		}

		TLegend* leta = new TLegend(0.59,0.55,0.90,0.85);
		leta->SetBorderSize(1);
		leta->SetFillColor(0);
		leta->SetFillStyle(0);

		// for counting subruns
		int tempint = 0;

		for (int i=pointcounter;i<pointcounter+looplimit;i++)
		{

			sstream.str(string());
			histoname = "unmatchedEta_";
			sstream << histoname << _runnumber.at(i);
			TH1D* thiseta = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);

			sstream.str(string());
			histoname = "unmatchedEtaIntegral_";
			sstream << histoname << _runnumber.at(i);
			TH1D* thiseta2 = dynamic_cast<TH1D*> ( _rootObjectMap[sstream.str()]);

			Float_t scale = highestentrycount / thiseta2->GetMaximum();
			//cout << "the scale here is : " << scale << endl;

			thiseta->Scale(scale);

			gStyle->SetOptStat(kFALSE);
			thiseta->SetStats(kFALSE);
			thiseta->SetLineColor(tempint+1);

			// jump yellow
			if ( tempint >= 4 )
			{
				thiseta->SetLineColor(tempint+2);
			}

			char o [100];
			char title[300];
			char title2[300];
			sprintf(o,"Epi");
			sstream.str(string());
			string k = _sensorname.at(i);
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			double fl = _irradfluence.at(i);
			int m = _biasvoltage.at(i);
			int n = _dutrotation.at(i);

			sprintf(title2, "%s, F=%.1e, %irot", sensname.c_str(),fl,n);
			if (tempint==0)
			{
				thiseta->Draw("E");
				thiseta->SetTitle(title2);
			} else {
				thiseta->Draw("E same");
			}
			thiseta->SetMaximum(yscale);

			canv_eta->Update();

			sprintf(title, "%s, F=%.1e, %iV, %irot", sensname.c_str(),fl,m,n);
			leta->AddEntry(thiseta,title,"lep");

			tempint++;

			double counts = 0.0;
			double integral = 0.0;
			for (int kk = 1; kk<thiseta->GetNbinsX(); kk++)
			{
				counts = thiseta->GetBinContent(kk);
				integral += counts ;
			}

		} // done subrun loop

		leta->Draw();
		canv_eta->Write();
		canv_eta->Close();
		delete leta;
		delete canv_eta;

		pointcounter=pointcounter+looplimit;

	} // done voltagegraph loop
	_outputFile->cd();

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Residual X vs. Voltage Plot..." << endl;
	}

	// residual vs voltage
	_outputFile->cd();
	TDirectory* resixdirectory = _outputFile->mkdir("Residuals X");
	resixdirectory->cd();

	// first x
	TCanvas* canvas_residualsvsvoltage_X = new TCanvas("canvas_residualsvsvoltage_X","Residuals in X vs. Voltage",200,10,700,500);
	canvas_residualsvsvoltage_X->cd();
	TH2F *hresixvolt = new TH2F("hresixvolt","Residuals in X vs. Voltage",20,0,1000,20,_cut_minXresidual,_cut_maxXresidual);
	hresixvolt->SetXTitle("Bias Voltage [|V|]");
	hresixvolt->SetYTitle("Residual in X [mm]");
	hresixvolt->SetStats(0000);
	hresixvolt->Draw();

	TLegend *lresixvolt = new TLegend(0.59,0.55,0.90,0.85);
	lresixvolt->SetBorderSize(1);
	lresixvolt->SetFillColor(0);
	lresixvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsXvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* residualsXvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		residualsXvoltage->SetMarkerStyle(temp3+20);
		residualsXvoltage->SetMarkerColor(temp2+1);
		residualsXvoltage->SetMarkerSize(2);
		residualsXvoltage->SetLineColor(temp2+1);
		residualsXvoltage->SetLineWidth(2);
		residualsXvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		residualsXvoltage->SetTitle(title);
		residualsXvoltage->Draw("LP");
		lresixvolt->AddEntry(residualsXvoltage,title,"lep");

	}

	lresixvolt->Draw();
	canvas_residualsvsvoltage_X->Update();
	resixdirectory->cd();
	canvas_residualsvsvoltage_X->Write();
	canvas_residualsvsvoltage_X->Close();

	delete canvas_residualsvsvoltage_X;
	delete lresixvolt;
	delete hresixvolt;

	TCanvas* canvas_residualsvsvoltage_X_0 = new TCanvas("canvas_residualsvsvoltage_X_0","Residuals in X vs. Voltage, 0 rot",200,10,700,500);
	canvas_residualsvsvoltage_X_0->cd();
	TH2F *hresixvolt_0 = new TH2F("hresixvolt_0","Residuals in X vs. Voltage, 0 rot",20,0,1000,20,_cut_minXresidual,_cut_maxXresidual);
	hresixvolt_0->SetXTitle("Bias Voltage [|V|]");
	hresixvolt_0->SetYTitle("Residual in X [mm]");
	hresixvolt_0->SetStats(0000);
	hresixvolt_0->Draw();

	TLegend *lresixvolt_0 = new TLegend(0.59,0.55,0.90,0.85);
	lresixvolt_0->SetBorderSize(1);
	lresixvolt_0->SetFillColor(0);
	lresixvolt_0->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsXvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		// all angles smaller 20deg
		if (n < 20)
		{

			TGraphErrors* residualsXvoltage_0 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			residualsXvoltage_0->SetMarkerStyle(temp3+20);
			residualsXvoltage_0->SetMarkerColor(temp2+1);
			residualsXvoltage_0->SetMarkerSize(2);
			residualsXvoltage_0->SetLineColor(temp2+1);
			residualsXvoltage_0->SetLineWidth(2);
			residualsXvoltage_0->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			residualsXvoltage_0->SetTitle(title);
			residualsXvoltage_0->Draw("LP");
			lresixvolt_0->AddEntry(residualsXvoltage_0,title,"lep");
		}

	}

	lresixvolt_0->Draw();
	canvas_residualsvsvoltage_X_0->Update();
	resixdirectory->cd();
	canvas_residualsvsvoltage_X_0->Write();
	canvas_residualsvsvoltage_X_0->Close();

	delete canvas_residualsvsvoltage_X_0;
	delete lresixvolt_0;
	delete hresixvolt_0;

	TCanvas* canvas_residualsvsvoltage_X_25 = new TCanvas("canvas_residualsvsvoltage_X_25","Residuals in X vs. Voltage, 25 rot",200,10,700,500);
	canvas_residualsvsvoltage_X_25->cd();
	TH2F *hresixvolt_25 = new TH2F("hresixvolt_25","Residuals in X vs. Voltage, 25 rot",20,0,1000,20,_cut_minXresidual,_cut_maxXresidual);
	hresixvolt_25->SetXTitle("Bias Voltage [|V|]");
	hresixvolt_25->SetYTitle("Residual in X [mm]");
	hresixvolt_25->SetStats(0000);
	hresixvolt_25->Draw();

	TLegend *lresixvolt_25 = new TLegend(0.59,0.55,0.90,0.85);
	lresixvolt_25->SetBorderSize(1);
	lresixvolt_25->SetFillColor(0);
	lresixvolt_25->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsXvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		// all angles between 20deg and 30deg
		if (n >= 20 && n < 30)
		{

			TGraphErrors* residualsXvoltage_25 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			residualsXvoltage_25->SetMarkerStyle(temp3+20);
			residualsXvoltage_25->SetMarkerColor(temp2+1);
			residualsXvoltage_25->SetMarkerSize(2);
			residualsXvoltage_25->SetLineColor(temp2+1);
			residualsXvoltage_25->SetLineWidth(2);
			residualsXvoltage_25->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			residualsXvoltage_25->SetTitle(title);
			residualsXvoltage_25->Draw("LP");
			lresixvolt_25->AddEntry(residualsXvoltage_25,title,"lep");
		}

	}

	lresixvolt_25->Draw();
	canvas_residualsvsvoltage_X_25->Update();
	resixdirectory->cd();
	canvas_residualsvsvoltage_X_25->Write();
	canvas_residualsvsvoltage_X_25->Close();

	delete canvas_residualsvsvoltage_X_25;
	delete lresixvolt_25;
	delete hresixvolt_25;

	TCanvas* canvas_residualsvsvoltage_X_51 = new TCanvas("canvas_residualsvsvoltage_X_51","Residuals in X vs. Voltage, other rot",200,10,700,500);
	canvas_residualsvsvoltage_X_51->cd();
	TH2F *hresixvolt_51 = new TH2F("hresixvolt_51","Residuals in X vs. Voltage, other rot",20,0,1000,20,_cut_minXresidual,_cut_maxXresidual);
	hresixvolt_51->SetXTitle("Bias Voltage [|V|]");
	hresixvolt_51->SetYTitle("Residual in X [mm]");
	hresixvolt_51->SetStats(0000);
	hresixvolt_51->Draw();

	TLegend *lresixvolt_51 = new TLegend(0.59,0.55,0.90,0.85);
	lresixvolt_51->SetBorderSize(1);
	lresixvolt_51->SetFillColor(0);
	lresixvolt_51->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsXvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		// all angles over 30deg
		if (n >= 30)
		{

			TGraphErrors* residualsXvoltage_51 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			residualsXvoltage_51->SetMarkerStyle(temp3+20);
			residualsXvoltage_51->SetMarkerColor(temp2+1);
			residualsXvoltage_51->SetMarkerSize(2);
			residualsXvoltage_51->SetLineColor(temp2+1);
			residualsXvoltage_51->SetLineWidth(2);
			residualsXvoltage_51->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			residualsXvoltage_51->SetTitle(title);
			residualsXvoltage_51->Draw("LP");
			lresixvolt_51->AddEntry(residualsXvoltage_51,title,"lep");
		}

	}

	lresixvolt_51->Draw();
	canvas_residualsvsvoltage_X_51->Update();
	resixdirectory->cd();
	canvas_residualsvsvoltage_X_51->Write();
	canvas_residualsvsvoltage_X_51->Close();

	delete canvas_residualsvsvoltage_X_51;
	delete lresixvolt_51;
	delete hresixvolt_51;

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Residual Y vs. Voltage Plot..." << endl;
	}

	// then y
	_outputFile->cd();
	TDirectory* resiydirectory = _outputFile->mkdir("Residuals Y");
	resiydirectory->cd();

	TCanvas* canvas_residualsvsvoltage_Y = new TCanvas("canvas_residualsvsvoltage_Y","Residuals in Y vs. Voltage",200,10,700,500);
	canvas_residualsvsvoltage_Y->cd();
	TH2F *hresiyvolt = new TH2F("hresiyvolt","Residuals in Y vs. Voltage",20,0,1000,20,_cut_minYresidual,_cut_maxYresidual);
	hresiyvolt->SetXTitle("Bias Voltage [|V|]");
	hresiyvolt->SetYTitle("Residual in Y [mm]");
	hresiyvolt->SetStats(0000);
	hresiyvolt->Draw();

	// for the y plot we also add a simple line to guide the eye where the binary resolution should be.
	// binary is : pitch / sqrt(12) plus telescope resolution
	float binaryresolution = (sqrt(80.0*80.0/12 + _cut_telescoperesolution*_cut_telescoperesolution))/1000.0;
	TLine *binaryline = new TLine(1,binaryresolution,999,binaryresolution);
	binaryline->SetLineColor(kRed);
	binaryline->SetLineStyle(3);
	binaryline->SetLineWidth(3);
	binaryline->Draw();
	delete binaryline;

	TLegend *lresiyvolt = new TLegend(0.59,0.55,0.90,0.85);
	lresiyvolt->SetBorderSize(1);
	lresiyvolt->SetFillColor(0);
	lresiyvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsYvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* residualsYvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		residualsYvoltage->SetMarkerStyle(temp3+20);
		residualsYvoltage->SetMarkerColor(temp2+1);
		residualsYvoltage->SetMarkerSize(2);
		residualsYvoltage->SetLineColor(temp2+1);
		residualsYvoltage->SetLineWidth(2);
		residualsYvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		residualsYvoltage->SetTitle(title);
		residualsYvoltage->Draw("LP");
		lresiyvolt->AddEntry(residualsYvoltage,title,"lep");

	}

	lresiyvolt->Draw();
	canvas_residualsvsvoltage_Y->Update();
	resiydirectory->cd();
	canvas_residualsvsvoltage_Y->Write();
	canvas_residualsvsvoltage_Y->Close();

	delete canvas_residualsvsvoltage_Y;
	delete lresiyvolt;
	delete hresiyvolt;

	TCanvas* canvas_residualsvsvoltage_Y_0 = new TCanvas("canvas_residualsvsvoltage_Y_0","Residuals in Y vs. Voltage, 0 deg",200,10,700,500);
	canvas_residualsvsvoltage_Y_0->cd();
	TH2F *hresiyvolt_0 = new TH2F("hresiyvolt_0","Residuals in Y vs. Voltage, 0 deg",20,0,1000,20,_cut_minYresidual,_cut_maxYresidual);
	hresiyvolt_0->SetXTitle("Bias Voltage [|V|]");
	hresiyvolt_0->SetYTitle("Residual in Y [mm]");
	hresiyvolt_0->SetStats(0000);
	hresiyvolt_0->Draw();

	// for the y plot we also add a simple line to guide the eye where the binary resolution should be.
	TLine *binaryline_0 = new TLine(1,binaryresolution,999,binaryresolution);
	binaryline_0->SetLineColor(kRed);
	binaryline_0->SetLineStyle(3);
	binaryline_0->SetLineWidth(3);
	binaryline_0->Draw();
	delete binaryline_0;

	TLegend *lresiyvolt_0 = new TLegend(0.59,0.55,0.90,0.85);
	lresiyvolt_0->SetBorderSize(1);
	lresiyvolt_0->SetFillColor(0);
	lresiyvolt_0->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsYvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* residualsYvoltage_0 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			residualsYvoltage_0->SetMarkerStyle(temp3+20);
			residualsYvoltage_0->SetMarkerColor(temp2+1);
			residualsYvoltage_0->SetMarkerSize(2);
			residualsYvoltage_0->SetLineColor(temp2+1);
			residualsYvoltage_0->SetLineWidth(2);
			residualsYvoltage_0->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			residualsYvoltage_0->SetTitle(title);
			residualsYvoltage_0->Draw("LP");
			lresiyvolt_0->AddEntry(residualsYvoltage_0,title,"lep");
		}

	}

	lresiyvolt_0->Draw();
	canvas_residualsvsvoltage_Y_0->Update();
	resiydirectory->cd();
	canvas_residualsvsvoltage_Y_0->Write();
	canvas_residualsvsvoltage_Y_0->Close();

	delete canvas_residualsvsvoltage_Y_0;
	delete lresiyvolt_0;
	delete hresiyvolt_0;

	TCanvas* canvas_residualsvsvoltage_Y_25 = new TCanvas("canvas_residualsvsvoltage_Y_25","Residuals in Y vs. Voltage, 25 deg",200,10,700,500);
	canvas_residualsvsvoltage_Y_25->cd();
	TH2F *hresiyvolt_25 = new TH2F("hresiyvolt_25","Residuals in Y vs. Voltage, 25 deg",20,0,1000,20,_cut_minYresidual,_cut_maxYresidual);
	hresiyvolt_25->SetXTitle("Bias Voltage [|V|]");
	hresiyvolt_25->SetYTitle("Residual in Y [mm]");
	hresiyvolt_25->SetStats(0000);
	hresiyvolt_25->Draw();

	// for the y plot we also add a simple line to guide the eye where the binary resolution should be.
	TLine *binaryline_25 = new TLine(1,binaryresolution,999,binaryresolution);
	binaryline_25->SetLineColor(kRed);
	binaryline_25->SetLineStyle(3);
	binaryline_25->SetLineWidth(3);
	binaryline_25->Draw();
	delete binaryline_25;

	TLegend *lresiyvolt_25 = new TLegend(0.59,0.55,0.90,0.85);
	lresiyvolt_25->SetBorderSize(1);
	lresiyvolt_25->SetFillColor(0);
	lresiyvolt_25->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsYvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* residualsYvoltage_25 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			residualsYvoltage_25->SetMarkerStyle(temp3+20);
			residualsYvoltage_25->SetMarkerColor(temp2+1);
			residualsYvoltage_25->SetMarkerSize(2);
			residualsYvoltage_25->SetLineColor(temp2+1);
			residualsYvoltage_25->SetLineWidth(2);
			residualsYvoltage_25->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			residualsYvoltage_25->SetTitle(title);
			residualsYvoltage_25->Draw("LP");
			lresiyvolt_25->AddEntry(residualsYvoltage_25,title,"lep");
		}

	}

	lresiyvolt_25->Draw();
	canvas_residualsvsvoltage_Y_25->Update();
	resiydirectory->cd();
	canvas_residualsvsvoltage_Y_25->Write();
	canvas_residualsvsvoltage_Y_25->Close();

	delete canvas_residualsvsvoltage_Y_25;
	delete lresiyvolt_25;
	delete hresiyvolt_25;

	TCanvas* canvas_residualsvsvoltage_Y_51 = new TCanvas("canvas_residualsvsvoltage_Y_51","Residuals in Y vs. Voltage, other rot",200,10,700,500);
	canvas_residualsvsvoltage_Y_51->cd();
	TH2F *hresiyvolt_51 = new TH2F("hresiyvolt_51","Residuals in Y vs. Voltage, other rot",20,0,1000,20,_cut_minYresidual,_cut_maxYresidual);
	hresiyvolt_51->SetXTitle("Bias Voltage [|V|]");
	hresiyvolt_51->SetYTitle("Residual in Y [mm]");
	hresiyvolt_51->SetStats(0000);
	hresiyvolt_51->Draw();

	// for the y plot we also add a simple line to guide the eye where the binary resolution should be.
	TLine *binaryline_51 = new TLine(1,binaryresolution,999,binaryresolution);
	binaryline_51->SetLineColor(kRed);
	binaryline_51->SetLineStyle(3);
	binaryline_51->SetLineWidth(3);
	binaryline_51->Draw();
	delete binaryline_51;

	TLegend *lresiyvolt_51 = new TLegend(0.59,0.55,0.90,0.85);
	lresiyvolt_51->SetBorderSize(1);
	lresiyvolt_51->SetFillColor(0);
	lresiyvolt_51->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsYvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n > 30)
		{

			TGraphErrors* residualsYvoltage_51 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			residualsYvoltage_51->SetMarkerStyle(temp3+20);
			residualsYvoltage_51->SetMarkerColor(temp2+1);
			residualsYvoltage_51->SetMarkerSize(2);
			residualsYvoltage_51->SetLineColor(temp2+1);
			residualsYvoltage_51->SetLineWidth(2);
			residualsYvoltage_51->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			residualsYvoltage_51->SetTitle(title);
			residualsYvoltage_51->Draw("LP");
			lresiyvolt_51->AddEntry(residualsYvoltage_51,title,"lep");
		}

	}

	lresiyvolt_51->Draw();
	canvas_residualsvsvoltage_Y_51->Update();
	resiydirectory->cd();
	canvas_residualsvsvoltage_Y_51->Write();
	canvas_residualsvsvoltage_Y_51->Close();

	delete canvas_residualsvsvoltage_Y_51;
	delete lresiyvolt_51;
	delete hresiyvolt_51;

	

	_outputFile->cd();
	TDirectory* resiydirectoryAngle = _outputFile->mkdir("Angular Residuals");
	resiydirectoryAngle->cd();

	TCanvas* canvas_residualsvsvoltage_YAngle = new TCanvas("canvas_residualsvsvoltage_YAngle","Residuals in Y vs. DUT Rotation Angle",200,10,700,500);
	canvas_residualsvsvoltage_YAngle->cd();
	TH2F *hresiyvoltAngle = new TH2F("hresiyvoltAngle","Residuals in Y vs. DUT Rotation Angle",20,-5,95,20,_cut_minYresidual,_cut_maxYresidual);
	hresiyvoltAngle->SetXTitle("DUT Rotation Angle [#circ]");
	hresiyvoltAngle->SetYTitle("Residual in Y [mm]");
	hresiyvoltAngle->SetStats(0000);
	hresiyvoltAngle->Draw();

	TLegend *lresiyvoltAngle = new TLegend(0.59,0.55,0.90,0.85);
	lresiyvoltAngle->SetBorderSize(1);
	lresiyvoltAngle->SetFillColor(0);
	lresiyvoltAngle->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "residualsYvoltageAngle_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* residualsYvoltageAngle = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		residualsYvoltageAngle->SetMarkerStyle(temp3+20);
		residualsYvoltageAngle->SetMarkerColor(temp2+1);
		residualsYvoltageAngle->SetMarkerSize(2);
		residualsYvoltageAngle->SetLineColor(temp2+1);
		residualsYvoltageAngle->SetLineWidth(2);
		residualsYvoltageAngle->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		residualsYvoltageAngle->SetTitle(title);
		residualsYvoltageAngle->Draw("LP");
		lresiyvoltAngle->AddEntry(residualsYvoltageAngle,title,"lep");

	}

	lresiyvoltAngle->Draw();
	canvas_residualsvsvoltage_YAngle->Update();
	resiydirectoryAngle->cd();
	canvas_residualsvsvoltage_YAngle->Write();
	canvas_residualsvsvoltage_YAngle->Close();

	delete canvas_residualsvsvoltage_YAngle;
	delete lresiyvoltAngle;
	delete hresiyvoltAngle;


	// done residual vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Signal vs. Voltage Plot..." << endl;
	}

	// signal vs voltage
	_outputFile->cd();
	TDirectory* signalsdirectory = _outputFile->mkdir("Signals");
	signalsdirectory->cd();

	TCanvas* canvas_signalvoltage = new TCanvas("canvas_signalvoltage","Signal vs. Voltage",200,10,700,500);
	canvas_signalvoltage->cd();
	TH2F *hsigvolt = new TH2F("hsigvolt","Signal vs. Voltage",20,0,1000,20,0,100);
	hsigvolt->SetXTitle("Bias Voltage [|V|]");
	hsigvolt->SetYTitle("Landau MPV [ADCs]");
	hsigvolt->SetStats(0000);
	hsigvolt->Draw();

	TLegend *lsigvolt = new TLegend(0.59,0.55,0.90,0.85);
	lsigvolt->SetBorderSize(1);
	lsigvolt->SetFillColor(0);
	lsigvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* signalvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		signalvoltage->SetMarkerStyle(temp3+20);
		signalvoltage->SetMarkerColor(temp2+1);
		signalvoltage->SetMarkerSize(2);
		signalvoltage->SetLineColor(temp2+1);
		signalvoltage->SetLineWidth(2);
		signalvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		signalvoltage->SetTitle(title);
		signalvoltage->Draw("LP");
		lsigvolt->AddEntry(signalvoltage,title,"lep");

	}

	lsigvolt->Draw();
	canvas_signalvoltage->Update();
	signalsdirectory->cd();
	canvas_signalvoltage->Write();
	canvas_signalvoltage->Close();

	delete canvas_signalvoltage;
	delete hsigvolt;
	delete lsigvolt;

	TCanvas* canvas_signalvoltage_0 = new TCanvas("canvas_signalvoltage_0","Signal vs. Voltage, 0 deg",200,10,700,500);
	canvas_signalvoltage_0->cd();
	TH2F *hsigvolt_0 = new TH2F("hsigvolt_0","Signal vs. Voltage, 0 deg",20,0,1000,20,0,100);
	hsigvolt_0->SetXTitle("Bias Voltage [|V|]");
	hsigvolt_0->SetYTitle("Landau MPV [ADCs]");
	hsigvolt_0->SetStats(0000);
	hsigvolt_0->Draw();

	TLegend *lsigvolt_0 = new TLegend(0.59,0.55,0.90,0.85);
	lsigvolt_0->SetBorderSize(1);
	lsigvolt_0->SetFillColor(0);
	lsigvolt_0->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* signalvoltage_0 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			signalvoltage_0->SetMarkerStyle(temp3+20);
			signalvoltage_0->SetMarkerColor(temp2+1);
			signalvoltage_0->SetMarkerSize(2);
			signalvoltage_0->SetLineColor(temp2+1);
			signalvoltage_0->SetLineWidth(2);
			signalvoltage_0->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			signalvoltage_0->SetTitle(title);
			signalvoltage_0->Draw("LP");
			lsigvolt_0->AddEntry(signalvoltage_0,title,"lep");
		}

	}

	lsigvolt_0->Draw();
	canvas_signalvoltage_0->Update();
	signalsdirectory->cd();
	canvas_signalvoltage_0->Write();
	canvas_signalvoltage_0->Close();

	delete canvas_signalvoltage_0;
	delete hsigvolt_0;
	delete lsigvolt_0;

	TCanvas* canvas_signalvoltage_25 = new TCanvas("canvas_signalvoltage_25","Signal vs. Voltage, 25 deg",200,10,700,500);
	canvas_signalvoltage_25->cd();
	TH2F *hsigvolt_25 = new TH2F("hsigvolt_25","Signal vs. Voltage, 25 deg",20,0,1000,20,0,100);
	hsigvolt_25->SetXTitle("Bias Voltage [|V|]");
	hsigvolt_25->SetYTitle("Landau MPV [ADCs]");
	hsigvolt_25->SetStats(0000);
	hsigvolt_25->Draw();

	TLegend *lsigvolt_25 = new TLegend(0.59,0.55,0.90,0.85);
	lsigvolt_25->SetBorderSize(1);
	lsigvolt_25->SetFillColor(0);
	lsigvolt_25->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* signalvoltage_25 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			signalvoltage_25->SetMarkerStyle(temp3+20);
			signalvoltage_25->SetMarkerColor(temp2+1);
			signalvoltage_25->SetMarkerSize(2);
			signalvoltage_25->SetLineColor(temp2+1);
			signalvoltage_25->SetLineWidth(2);
			signalvoltage_25->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			signalvoltage_25->SetTitle(title);
			signalvoltage_25->Draw("LP");
			lsigvolt_25->AddEntry(signalvoltage_25,title,"lep");
		}

	}

	lsigvolt_25->Draw();
	canvas_signalvoltage_25->Update();
	signalsdirectory->cd();
	canvas_signalvoltage_25->Write();
	canvas_signalvoltage_25->Close();

	delete canvas_signalvoltage_25;
	delete hsigvolt_25;
	delete lsigvolt_25;

	TCanvas* canvas_signalvoltage_51 = new TCanvas("canvas_signalvoltage_51","Signal vs. Voltage, other rot",200,10,700,500);
	canvas_signalvoltage_51->cd();
	TH2F *hsigvolt_51 = new TH2F("hsigvolt_51","Signal vs. Voltage, other rot",20,0,1000,20,0,100);
	hsigvolt_51->SetXTitle("Bias Voltage [|V|]");
	hsigvolt_51->SetYTitle("Landau MPV [ADCs]");
	hsigvolt_51->SetStats(0000);
	hsigvolt_51->Draw();

	TLegend *lsigvolt_51 = new TLegend(0.59,0.55,0.90,0.85);
	lsigvolt_51->SetBorderSize(1);
	lsigvolt_51->SetFillColor(0);
	lsigvolt_51->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 30)
		{

			TGraphErrors* signalvoltage_51 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			signalvoltage_51->SetMarkerStyle(temp3+20);
			signalvoltage_51->SetMarkerColor(temp2+1);
			signalvoltage_51->SetMarkerSize(2);
			signalvoltage_51->SetLineColor(temp2+1);
			signalvoltage_51->SetLineWidth(2);
			signalvoltage_51->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			signalvoltage_51->SetTitle(title);
			signalvoltage_51->Draw("LP");
			lsigvolt_51->AddEntry(signalvoltage_51,title,"lep");
		}

	}

	lsigvolt_51->Draw();
	canvas_signalvoltage_51->Update();
	signalsdirectory->cd();
	canvas_signalvoltage_51->Write();
	canvas_signalvoltage_51->Close();

	delete canvas_signalvoltage_51;
	delete hsigvolt_51;
	delete lsigvolt_51;

	TCanvas* canvas_signalvoltageelec = new TCanvas("canvas_signalvoltageelec","Signal vs. Voltage, Electrons",200,10,700,500);
	canvas_signalvoltageelec->cd();
	TH2F *hsigvoltelec = new TH2F("hsigvoltelec","Signal vs. Voltage, Electrons",20,0,1000,20,0,15000);
	hsigvoltelec->SetXTitle("Bias Voltage [|V|]");
	hsigvoltelec->SetYTitle("Landau MPV [Electrons]");
	hsigvoltelec->SetStats(0000);
	hsigvoltelec->Draw();

	TLegend *lsigvoltelec = new TLegend(0.59,0.55,0.90,0.85);
	lsigvoltelec->SetBorderSize(1);
	lsigvoltelec->SetFillColor(0);
	lsigvoltelec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* signalvoltageelec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		signalvoltageelec->SetMarkerStyle(temp3+20);
		signalvoltageelec->SetMarkerColor(temp2+1);
		signalvoltageelec->SetMarkerSize(2);
		signalvoltageelec->SetLineColor(temp2+1);
		signalvoltageelec->SetLineWidth(2);
		signalvoltageelec->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		signalvoltageelec->SetTitle(title);
		signalvoltageelec->Draw("LP");
		lsigvoltelec->AddEntry(signalvoltageelec,title,"lep");

	}

	lsigvoltelec->Draw();
	canvas_signalvoltageelec->Update();
	signalsdirectory->cd();
	canvas_signalvoltageelec->Write();
	canvas_signalvoltageelec->Close();

	delete canvas_signalvoltageelec;
	delete hsigvoltelec;
	delete lsigvoltelec;

	TCanvas* canvas_signalvoltage_0elec = new TCanvas("canvas_signalvoltage_0elec","Signal vs. Voltage, 0 deg, Electrons",200,10,700,500);
	canvas_signalvoltage_0elec->cd();
	TH2F *hsigvolt_0elec = new TH2F("hsigvolt_0elec","Signal vs. Voltage, 0 deg, Electrons",20,0,1000,20,0,15000);
	hsigvolt_0elec->SetXTitle("Bias Voltage [|V|]");
	hsigvolt_0elec->SetYTitle("Landau MPV [Electrons]");
	hsigvolt_0elec->SetStats(0000);
	hsigvolt_0elec->Draw();

	TLegend *lsigvolt_0elec = new TLegend(0.59,0.55,0.90,0.85);
	lsigvolt_0elec->SetBorderSize(1);
	lsigvolt_0elec->SetFillColor(0);
	lsigvolt_0elec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* signalvoltage_0elec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			signalvoltage_0elec->SetMarkerStyle(temp3+20);
			signalvoltage_0elec->SetMarkerColor(temp2+1);
			signalvoltage_0elec->SetMarkerSize(2);
			signalvoltage_0elec->SetLineColor(temp2+1);
			signalvoltage_0elec->SetLineWidth(2);
			signalvoltage_0elec->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			signalvoltage_0elec->SetTitle(title);
			signalvoltage_0elec->Draw("LP");
			lsigvolt_0elec->AddEntry(signalvoltage_0elec,title,"lep");
		}

	}

	lsigvolt_0elec->Draw();
	canvas_signalvoltage_0elec->Update();
	signalsdirectory->cd();
	canvas_signalvoltage_0elec->Write();
	canvas_signalvoltage_0elec->Close();

	delete canvas_signalvoltage_0elec;
	delete hsigvolt_0elec;
	delete lsigvolt_0elec;

	TCanvas* canvas_signalvoltage_25elec = new TCanvas("canvas_signalvoltage_25elec","Signal vs. Voltage, 25 deg, Electrons",200,10,700,500);
	canvas_signalvoltage_25elec->cd();
	TH2F *hsigvolt_25elec = new TH2F("hsigvolt_25elec","Signal vs. Voltage, 25 deg, Electrons",20,0,1000,20,0,15000);
	hsigvolt_25elec->SetXTitle("Bias Voltage [|V|]");
	hsigvolt_25elec->SetYTitle("Landau MPV [Electrons]");
	hsigvolt_25elec->SetStats(0000);
	hsigvolt_25elec->Draw();

	TLegend *lsigvolt_25elec = new TLegend(0.59,0.55,0.90,0.85);
	lsigvolt_25elec->SetBorderSize(1);
	lsigvolt_25elec->SetFillColor(0);
	lsigvolt_25elec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* signalvoltage_25elec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			signalvoltage_25elec->SetMarkerStyle(temp3+20);
			signalvoltage_25elec->SetMarkerColor(temp2+1);
			signalvoltage_25elec->SetMarkerSize(2);
			signalvoltage_25elec->SetLineColor(temp2+1);
			signalvoltage_25elec->SetLineWidth(2);
			signalvoltage_25elec->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			signalvoltage_25elec->SetTitle(title);
			signalvoltage_25elec->Draw("LP");
			lsigvolt_25elec->AddEntry(signalvoltage_25elec,title,"lep");
		}

	}

	lsigvolt_25elec->Draw();
	canvas_signalvoltage_25elec->Update();
	signalsdirectory->cd();
	canvas_signalvoltage_25elec->Write();
	canvas_signalvoltage_25elec->Close();

	delete canvas_signalvoltage_25elec;
	delete hsigvolt_25elec;
	delete lsigvolt_25elec;

	TCanvas* canvas_signalvoltage_51elec = new TCanvas("canvas_signalvoltage_51elec","Signal vs. Voltage, other rot, Electrons",200,10,700,500);
	canvas_signalvoltage_51elec->cd();
	TH2F *hsigvolt_51elec = new TH2F("hsigvolt_51elec","Signal vs. Voltage, other rot, Electrons",20,0,1000,20,0,15000);
	hsigvolt_51elec->SetXTitle("Bias Voltage [|V|]");
	hsigvolt_51elec->SetYTitle("Landau MPV [Electrons]");
	hsigvolt_51elec->SetStats(0000);
	hsigvolt_51elec->Draw();

	TLegend *lsigvolt_51elec = new TLegend(0.59,0.55,0.90,0.85);
	lsigvolt_51elec->SetBorderSize(1);
	lsigvolt_51elec->SetFillColor(0);
	lsigvolt_51elec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signalvoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 30)
		{

			TGraphErrors* signalvoltage_51elec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			signalvoltage_51elec->SetMarkerStyle(temp3+20);
			signalvoltage_51elec->SetMarkerColor(temp2+1);
			signalvoltage_51elec->SetMarkerSize(2);
			signalvoltage_51elec->SetLineColor(temp2+1);
			signalvoltage_51elec->SetLineWidth(2);
			signalvoltage_51elec->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			signalvoltage_51elec->SetTitle(title);
			signalvoltage_51elec->Draw("LP");
			lsigvolt_51elec->AddEntry(signalvoltage_51elec,title,"lep");
		}

	}

	lsigvolt_51elec->Draw();
	canvas_signalvoltage_51elec->Update();
	signalsdirectory->cd();
	canvas_signalvoltage_51elec->Write();
	canvas_signalvoltage_51elec->Close();

	delete canvas_signalvoltage_51elec;
	delete hsigvolt_51elec;
	delete lsigvolt_51elec;

	// done signal vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Signal/Noise vs. Voltage Plot..." << endl;
	}

	// signal/noise vs voltage
	_outputFile->cd();
	TDirectory* signalsnoisedirectory = _outputFile->mkdir("Signal to Noise");
	signalsnoisedirectory->cd();

	TCanvas* canvas_snvoltage = new TCanvas("canvas_snvoltage","Signal/Noise vs. Voltage",200,10,700,500);
	canvas_snvoltage->cd();
	TH2F *hsnvolt = new TH2F("hsnvolt","Signal/noise vs. Voltage",20,0,1000,20,0,25);
	hsnvolt->SetXTitle("Bias Voltage [|V|]");
	hsnvolt->SetYTitle("Signal to Noise [1]");
	hsnvolt->SetStats(0000);
	hsnvolt->Draw();

	TLegend *lsnvolt = new TLegend(0.59,0.55,0.90,0.85);
	lsnvolt->SetBorderSize(1);
	lsnvolt->SetFillColor(0);
	lsnvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "snvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* snvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		snvoltage->SetMarkerStyle(temp3+20);
		snvoltage->SetMarkerColor(temp2+1);
		snvoltage->SetMarkerSize(2);
		snvoltage->SetLineColor(temp2+1);
		snvoltage->SetLineWidth(2);
		snvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		snvoltage->SetTitle(title);
		snvoltage->Draw("LP");
		lsnvolt->AddEntry(snvoltage,title,"lep");

	}

	lsnvolt->Draw();
	canvas_snvoltage->Update();
	signalsnoisedirectory->cd();
	canvas_snvoltage->Write();
	canvas_snvoltage->Close();

	delete canvas_snvoltage;
	delete hsnvolt;
	delete lsnvolt;

	TCanvas* canvas_snvoltage_0 = new TCanvas("canvas_snvoltage_0","Signal/Noise vs. Voltage, 0 deg",200,10,700,500);
	canvas_snvoltage_0->cd();
	TH2F *hsnvolt_0 = new TH2F("hsnvolt_0","Signal/noise vs. Voltage, 0 deg",20,0,1000,20,0,25);
	hsnvolt_0->SetXTitle("Bias Voltage [|V|]");
	hsnvolt_0->SetYTitle("Signal to Noise [1]");
	hsnvolt_0->SetStats(0000);
	hsnvolt_0->Draw();

	TLegend *lsnvolt_0 = new TLegend(0.59,0.55,0.90,0.85);
	lsnvolt_0->SetBorderSize(1);
	lsnvolt_0->SetFillColor(0);
	lsnvolt_0->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "snvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* snvoltage_0 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			snvoltage_0->SetMarkerStyle(temp3+20);
			snvoltage_0->SetMarkerColor(temp2+1);
			snvoltage_0->SetMarkerSize(2);
			snvoltage_0->SetLineColor(temp2+1);
			snvoltage_0->SetLineWidth(2);
			snvoltage_0->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			snvoltage_0->SetTitle(title);
			snvoltage_0->Draw("LP");
			lsnvolt_0->AddEntry(snvoltage_0,title,"lep");
		}

	}

	lsnvolt_0->Draw();
	canvas_snvoltage_0->Update();
	signalsnoisedirectory->cd();
	canvas_snvoltage_0->Write();
	canvas_snvoltage_0->Close();

	delete canvas_snvoltage_0;
	delete hsnvolt_0;
	delete lsnvolt_0;

	TCanvas* canvas_snvoltage_25 = new TCanvas("canvas_snvoltage_25","Signal/Noise vs. Voltage, 25 deg",200,10,700,500);
	canvas_snvoltage_25->cd();
	TH2F *hsnvolt_25 = new TH2F("hsnvolt_25","Signal/noise vs. Voltage, 25 deg",20,0,1000,20,0,25);
	hsnvolt_25->SetXTitle("Bias Voltage [|V|]");
	hsnvolt_25->SetYTitle("Signal to Noise [1]");
	hsnvolt_25->SetStats(0000);
	hsnvolt_25->Draw();

	TLegend *lsnvolt_25 = new TLegend(0.59,0.55,0.90,0.85);
	lsnvolt_25->SetBorderSize(1);
	lsnvolt_25->SetFillColor(0);
	lsnvolt_25->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "snvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* snvoltage_25 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			snvoltage_25->SetMarkerStyle(temp3+20);
			snvoltage_25->SetMarkerColor(temp2+1);
			snvoltage_25->SetMarkerSize(2);
			snvoltage_25->SetLineColor(temp2+1);
			snvoltage_25->SetLineWidth(2);
			snvoltage_25->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			snvoltage_25->SetTitle(title);
			snvoltage_25->Draw("LP");
			lsnvolt_25->AddEntry(snvoltage_25,title,"lep");
		}

	}

	lsnvolt_25->Draw();
	canvas_snvoltage_25->Update();
	signalsnoisedirectory->cd();
	canvas_snvoltage_25->Write();
	canvas_snvoltage_25->Close();

	delete canvas_snvoltage_25;
	delete hsnvolt_25;
	delete lsnvolt_25;

	TCanvas* canvas_snvoltage_51 = new TCanvas("canvas_snvoltage_51","Signal/Noise vs. Voltage, other rot",200,10,700,500);
	canvas_snvoltage_51->cd();
	TH2F *hsnvolt_51 = new TH2F("hsnvolt_51","Signal/noise vs. Voltage, other rot",20,0,1000,20,0,25);
	hsnvolt_51->SetXTitle("Bias Voltage [|V|]");
	hsnvolt_51->SetYTitle("Signal to Noise [1]");
	hsnvolt_51->SetStats(0000);
	hsnvolt_51->Draw();

	TLegend *lsnvolt_51 = new TLegend(0.59,0.55,0.90,0.85);
	lsnvolt_51->SetBorderSize(1);
	lsnvolt_51->SetFillColor(0);
	lsnvolt_51->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "snvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 30)
		{

			TGraphErrors* snvoltage_51 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			snvoltage_51->SetMarkerStyle(temp3+20);
			snvoltage_51->SetMarkerColor(temp2+1);
			snvoltage_51->SetMarkerSize(2);
			snvoltage_51->SetLineColor(temp2+1);
			snvoltage_51->SetLineWidth(2);
			snvoltage_51->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			snvoltage_51->SetTitle(title);
			snvoltage_51->Draw("LP");
			lsnvolt_51->AddEntry(snvoltage_51,title,"lep");
		}

	}

	lsnvolt_51->Draw();
	canvas_snvoltage_51->Update();
	signalsnoisedirectory->cd();
	canvas_snvoltage_51->Write();
	canvas_snvoltage_51->Close();

	delete canvas_snvoltage_51;
	delete hsnvolt_51;
	delete lsnvolt_51;

	// done signal vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Signal/Effective Thickness vs. Voltage Plot..." << endl;
	}

	// signal/noise vs voltage

	TCanvas* canvas_signaleffthickness = new TCanvas("canvas_signaleffthickness","Signal/Effective Thickness vs. Voltage",200,10,700,500);
	canvas_signaleffthickness->cd();
	TH2F *hseffthickvolt = new TH2F("hseffthickvolt","Signal/Effective Thickness vs. Voltage",20,0,1000,20,0,1);
	hseffthickvolt->SetXTitle("Bias Voltage [|V|]");
	hseffthickvolt->SetYTitle("Signal / Effective Thickness [ADC / #mum]");
	hseffthickvolt->SetStats(0000);
	hseffthickvolt->Draw();

	TLegend *leffthickvolt = new TLegend(0.59,0.55,0.90,0.85);
	leffthickvolt->SetBorderSize(1);
	leffthickvolt->SetFillColor(0);
	leffthickvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signaldistancevoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* signaldistancevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		signaldistancevoltage->SetMarkerStyle(temp3+20);
		signaldistancevoltage->SetMarkerColor(temp2+1);
		signaldistancevoltage->SetMarkerSize(2);
		signaldistancevoltage->SetLineColor(temp2+1);
		signaldistancevoltage->SetLineWidth(2);
		signaldistancevoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		signaldistancevoltage->SetTitle(title);
		signaldistancevoltage->Draw("LP");
		leffthickvolt->AddEntry(signaldistancevoltage,title,"lep");

	}

	leffthickvolt->Draw();
	canvas_signaleffthickness->Update();
	_outputFile->cd();
	canvas_signaleffthickness->Write();
	canvas_signaleffthickness->Close();

	delete canvas_signaleffthickness;
	delete hseffthickvolt;
	delete leffthickvolt;

	TCanvas* canvas_signaleffthicknesselec = new TCanvas("canvas_signaleffthicknesselec","Signal/Effective Thickness vs. Voltage, Electrons",200,10,700,500);
	canvas_signaleffthicknesselec->cd();
	TH2F *hseffthickvoltelec = new TH2F("hseffthickvoltelec","Signal/Effective Thickness vs. Voltage, Electrons",20,0,1000,20,0,200);
	hseffthickvoltelec->SetXTitle("Bias Voltage [|V|]");
	hseffthickvoltelec->SetYTitle("Signal / Effective Thickness [Electrons / #mum]");
	hseffthickvoltelec->SetStats(0000);
	hseffthickvoltelec->Draw();

	TLegend *leffthickvoltelec = new TLegend(0.59,0.55,0.90,0.85);
	leffthickvoltelec->SetBorderSize(1);
	leffthickvoltelec->SetFillColor(0);
	leffthickvoltelec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "signaldistancevoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* signaldistancevoltageelec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		signaldistancevoltageelec->SetMarkerStyle(temp3+20);
		signaldistancevoltageelec->SetMarkerColor(temp2+1);
		signaldistancevoltageelec->SetMarkerSize(2);
		signaldistancevoltageelec->SetLineColor(temp2+1);
		signaldistancevoltageelec->SetLineWidth(2);
		signaldistancevoltageelec->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		signaldistancevoltageelec->SetTitle(title);
		signaldistancevoltageelec->Draw("LP");
		leffthickvoltelec->AddEntry(signaldistancevoltageelec,title,"lep");

	}

	leffthickvoltelec->Draw();
	canvas_signaleffthicknesselec->Update();
	_outputFile->cd();
	canvas_signaleffthicknesselec->Write();
	canvas_signaleffthicknesselec->Close();

	delete canvas_signaleffthicknesselec;
	delete hseffthickvoltelec;
	delete leffthickvoltelec;

	// done signal vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Noise vs. Voltage Plot..." << endl;
	}

	// noise vs voltage
	_outputFile->cd();
	TDirectory* overallnoisedirectory = _outputFile->mkdir("Noise");
	overallnoisedirectory->cd();

	TCanvas* canvas_noisevoltage = new TCanvas("canvas_noisevoltage","Noise vs. Voltage",200,10,700,500);
	canvas_noisevoltage->cd();
	TH2F *hnoisevolt = new TH2F("hnoisevolt","Noise vs. Voltage",20,0,1000,50,0,50);
	hnoisevolt->SetXTitle("Bias Voltage [|V|]");
	hnoisevolt->SetYTitle("Average Good Channel Noise [ADCs]");
	hnoisevolt->SetStats(0000);
	hnoisevolt->Draw();

	TLegend *lnoisevolt = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevolt->SetBorderSize(1);
	lnoisevolt->SetFillColor(0);
	lnoisevolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* noisevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		noisevoltage->SetMarkerStyle(temp3+20);
		noisevoltage->SetMarkerColor(temp2+1);
		noisevoltage->SetMarkerSize(2);
		noisevoltage->SetLineColor(temp2+1);
		noisevoltage->SetLineWidth(2);
		noisevoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		noisevoltage->SetTitle(title);
		noisevoltage->Draw("LP");
		lnoisevolt->AddEntry(noisevoltage,title,"lep");

	}

	lnoisevolt->Draw();
	canvas_noisevoltage->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltage->Write();
	canvas_noisevoltage->Close();

	delete canvas_noisevoltage;
	delete lnoisevolt;
	delete hnoisevolt;

	TCanvas* canvas_noisevoltage_0 = new TCanvas("canvas_noisevoltage_0","Noise vs. Voltage, 0 deg",200,10,700,500);
	canvas_noisevoltage_0->cd();
	TH2F *hnoisevolt_0 = new TH2F("hnoisevolt_0","Noise vs. Voltage, 0 deg",20,0,1000,20,0,50);
	hnoisevolt_0->SetXTitle("Bias Voltage [|V|]");
	hnoisevolt_0->SetYTitle("Average Good Channel Noise [ADCs]");
	hnoisevolt_0->SetStats(0000);
	hnoisevolt_0->Draw();

	TLegend *lnoisevolt_0 = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevolt_0->SetBorderSize(1);
	lnoisevolt_0->SetFillColor(0);
	lnoisevolt_0->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* noisevoltage_0 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			noisevoltage_0->SetMarkerStyle(temp3+20);
			noisevoltage_0->SetMarkerColor(temp2+1);
			noisevoltage_0->SetMarkerSize(2);
			noisevoltage_0->SetLineColor(temp2+1);
			noisevoltage_0->SetLineWidth(2);
			noisevoltage_0->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			noisevoltage_0->SetTitle(title);
			noisevoltage_0->Draw("LP");
			lnoisevolt_0->AddEntry(noisevoltage_0,title,"lep");
		}

	}

	lnoisevolt_0->Draw();
	canvas_noisevoltage_0->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltage_0->Write();
	canvas_noisevoltage_0->Close();

	delete canvas_noisevoltage_0;
	delete lnoisevolt_0;
	delete hnoisevolt_0;

	TCanvas* canvas_noisevoltage_25 = new TCanvas("canvas_noisevoltage_25","Noise vs. Voltage, 25 deg",200,10,700,500);
	canvas_noisevoltage_25->cd();
	TH2F *hnoisevolt_25 = new TH2F("hnoisevolt_25","Noise vs. Voltage, 25 deg",20,0,1000,20,0,50);
	hnoisevolt_25->SetXTitle("Bias Voltage [|V|]");
	hnoisevolt_25->SetYTitle("Average Good Channel Noise [ADCs]");
	hnoisevolt_25->SetStats(0000);
	hnoisevolt_25->Draw();

	TLegend *lnoisevolt_25 = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevolt_25->SetBorderSize(1);
	lnoisevolt_25->SetFillColor(0);
	lnoisevolt_25->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* noisevoltage_25 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			noisevoltage_25->SetMarkerStyle(temp3+20);
			noisevoltage_25->SetMarkerColor(temp2+1);
			noisevoltage_25->SetMarkerSize(2);
			noisevoltage_25->SetLineColor(temp2+1);
			noisevoltage_25->SetLineWidth(2);
			noisevoltage_25->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			noisevoltage_25->SetTitle(title);
			noisevoltage_25->Draw("LP");
			lnoisevolt_25->AddEntry(noisevoltage_25,title,"lep");
		}

	}

	lnoisevolt_25->Draw();
	canvas_noisevoltage_25->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltage_25->Write();
	canvas_noisevoltage_25->Close();

	delete canvas_noisevoltage_25;
	delete lnoisevolt_25;
	delete hnoisevolt_25;

	TCanvas* canvas_noisevoltage_51 = new TCanvas("canvas_noisevoltage_51","Noise vs. Voltage, other rot",200,10,700,500);
	canvas_noisevoltage_51->cd();
	TH2F *hnoisevolt_51 = new TH2F("hnoisevolt_51","Noise vs. Voltage, other rot",20,0,1000,20,0,50);
	hnoisevolt_51->SetXTitle("Bias Voltage [|V|]");
	hnoisevolt_51->SetYTitle("Average Good Channel Noise [ADCs]");
	hnoisevolt_51->SetStats(0000);
	hnoisevolt_51->Draw();

	TLegend *lnoisevolt_51 = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevolt_51->SetBorderSize(1);
	lnoisevolt_51->SetFillColor(0);
	lnoisevolt_51->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 30)
		{

			TGraphErrors* noisevoltage_51 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			noisevoltage_51->SetMarkerStyle(temp3+20);
			noisevoltage_51->SetMarkerColor(temp2+1);
			noisevoltage_51->SetMarkerSize(2);
			noisevoltage_51->SetLineColor(temp2+1);
			noisevoltage_51->SetLineWidth(2);
			noisevoltage_51->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			noisevoltage_51->SetTitle(title);
			noisevoltage_51->Draw("LP");
			lnoisevolt_51->AddEntry(noisevoltage_51,title,"lep");
		}

	}

	lnoisevolt_51->Draw();
	canvas_noisevoltage_51->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltage_51->Write();
	canvas_noisevoltage_51->Close();

	delete canvas_noisevoltage_51;
	delete lnoisevolt_51;
	delete hnoisevolt_51;

	TCanvas* canvas_noisevoltageelec = new TCanvas("canvas_noisevoltageelec","Noise vs. Voltage, Electrons",200,10,700,500);
	canvas_noisevoltageelec->cd();
	TH2F *hnoisevoltelec = new TH2F("hnoisevoltelec","Noise vs. Voltage, Electrons",20,0,1000,20,0,7500);
	hnoisevoltelec->SetXTitle("Bias Voltage [|V|]");
	hnoisevoltelec->SetYTitle("Average Good Channel Noise [Electrons]");
	hnoisevoltelec->SetStats(0000);
	hnoisevoltelec->Draw();

	TLegend *lnoisevoltelec = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevoltelec->SetBorderSize(1);
	lnoisevoltelec->SetFillColor(0);
	lnoisevoltelec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* noisevoltageelec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		noisevoltageelec->SetMarkerStyle(temp3+20);
		noisevoltageelec->SetMarkerColor(temp2+1);
		noisevoltageelec->SetMarkerSize(2);
		noisevoltageelec->SetLineColor(temp2+1);
		noisevoltageelec->SetLineWidth(2);
		noisevoltageelec->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		noisevoltageelec->SetTitle(title);
		noisevoltageelec->Draw("LP");
		lnoisevoltelec->AddEntry(noisevoltageelec,title,"lep");

	}

	lnoisevoltelec->Draw();
	canvas_noisevoltageelec->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltageelec->Write();
	canvas_noisevoltageelec->Close();

	delete canvas_noisevoltageelec;
	delete lnoisevoltelec;
	delete hnoisevoltelec;

	TCanvas* canvas_noisevoltage_0elec = new TCanvas("canvas_noisevoltage_0elec","Noise vs. Voltage, 0 deg, Electrons",200,10,700,500);
	canvas_noisevoltage_0elec->cd();
	TH2F *hnoisevolt_0elec = new TH2F("hnoisevolt_0elec","Noise vs. Voltage, 0 deg, Electrons",20,0,1000,20,0,7500);
	hnoisevolt_0elec->SetXTitle("Bias Voltage [|V|]");
	hnoisevolt_0elec->SetYTitle("Average Good Channel Noise [Electrons]");
	hnoisevolt_0elec->SetStats(0000);
	hnoisevolt_0elec->Draw();

	TLegend *lnoisevolt_0elec = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevolt_0elec->SetBorderSize(1);
	lnoisevolt_0elec->SetFillColor(0);
	lnoisevolt_0elec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* noisevoltage_0elec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			noisevoltage_0elec->SetMarkerStyle(temp3+20);
			noisevoltage_0elec->SetMarkerColor(temp2+1);
			noisevoltage_0elec->SetMarkerSize(2);
			noisevoltage_0elec->SetLineColor(temp2+1);
			noisevoltage_0elec->SetLineWidth(2);
			noisevoltage_0elec->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			noisevoltage_0elec->SetTitle(title);
			noisevoltage_0elec->Draw("LP");
			lnoisevolt_0elec->AddEntry(noisevoltage_0elec,title,"lep");
		}

	}

	lnoisevolt_0elec->Draw();
	canvas_noisevoltage_0elec->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltage_0elec->Write();
	canvas_noisevoltage_0elec->Close();

	delete canvas_noisevoltage_0elec;
	delete lnoisevolt_0elec;
	delete hnoisevolt_0elec;

	TCanvas* canvas_noisevoltage_25elec = new TCanvas("canvas_noisevoltage_25elec","Noise vs. Voltage, 25 deg, Electrons",200,10,700,500);
	canvas_noisevoltage_25elec->cd();
	TH2F *hnoisevolt_25elec = new TH2F("hnoisevolt_25elec","Noise vs. Voltage, 25 deg, Electrons",20,0,1000,20,0,7500);
	hnoisevolt_25elec->SetXTitle("Bias Voltage [|V|]");
	hnoisevolt_25elec->SetYTitle("Average Good Channel Noise [Electrons]");
	hnoisevolt_25elec->SetStats(0000);
	hnoisevolt_25elec->Draw();

	TLegend *lnoisevolt_25elec = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevolt_25elec->SetBorderSize(1);
	lnoisevolt_25elec->SetFillColor(0);
	lnoisevolt_25elec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* noisevoltage_25elec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			noisevoltage_25elec->SetMarkerStyle(temp3+20);
			noisevoltage_25elec->SetMarkerColor(temp2+1);
			noisevoltage_25elec->SetMarkerSize(2);
			noisevoltage_25elec->SetLineColor(temp2+1);
			noisevoltage_25elec->SetLineWidth(2);
			noisevoltage_25elec->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			noisevoltage_25elec->SetTitle(title);
			noisevoltage_25elec->Draw("LP");
			lnoisevolt_25elec->AddEntry(noisevoltage_25elec,title,"lep");
		}

	}

	lnoisevolt_25elec->Draw();
	canvas_noisevoltage_25elec->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltage_25elec->Write();
	canvas_noisevoltage_25elec->Close();

	delete canvas_noisevoltage_25elec;
	delete lnoisevolt_25elec;
	delete hnoisevolt_25elec;

	TCanvas* canvas_noisevoltage_51elec = new TCanvas("canvas_noisevoltage_51elec","Noise vs. Voltage, other rot, Electrons",200,10,700,500);
	canvas_noisevoltage_51elec->cd();
	TH2F *hnoisevolt_51elec = new TH2F("hnoisevolt_51elec","Noise vs. Voltage, other rot, Electrons",20,0,1000,20,0,7500);
	hnoisevolt_51elec->SetXTitle("Bias Voltage [|V|]");
	hnoisevolt_51elec->SetYTitle("Average Good Channel Noise [Electrons]");
	hnoisevolt_51elec->SetStats(0000);
	hnoisevolt_51elec->Draw();

	TLegend *lnoisevolt_51elec = new TLegend(0.59,0.55,0.90,0.85);
	lnoisevolt_51elec->SetBorderSize(1);
	lnoisevolt_51elec->SetFillColor(0);
	lnoisevolt_51elec->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "noisevoltageelec_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 30)
		{

			TGraphErrors* noisevoltage_51elec = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			noisevoltage_51elec->SetMarkerStyle(temp3+20);
			noisevoltage_51elec->SetMarkerColor(temp2+1);
			noisevoltage_51elec->SetMarkerSize(2);
			noisevoltage_51elec->SetLineColor(temp2+1);
			noisevoltage_51elec->SetLineWidth(2);
			noisevoltage_51elec->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			noisevoltage_51elec->SetTitle(title);
			noisevoltage_51elec->Draw("LP");
			lnoisevolt_51elec->AddEntry(noisevoltage_51elec,title,"lep");
		}

	}

	lnoisevolt_51elec->Draw();
	canvas_noisevoltage_51elec->Update();
	overallnoisedirectory->cd();
	canvas_noisevoltage_51elec->Write();
	canvas_noisevoltage_51elec->Close();

	delete canvas_noisevoltage_51elec;
	delete lnoisevolt_51elec;
	delete hnoisevolt_51elec;

	// done noise vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting RGH vs. Voltage Plot..." << endl;
	}

	// rgh vs voltage
	TCanvas* canvas_rghvoltage = new TCanvas("canvas_rghvoltage","Noise vs. Voltage",200,10,700,500);
	canvas_rghvoltage->cd();
	TH2F *hrghvolt = new TH2F("hrghvolt","RGH Percentage vs. Voltage",20,0,1000,20,0,_cut_maxrghplot);
	hrghvolt->SetXTitle("Bias Voltage [|V|]");
	hrghvolt->SetYTitle("Percentage of Random Ghost Hits [%]");
	hrghvolt->SetStats(0000);
	hrghvolt->Draw();

	TLegend *lrghvolt = new TLegend(0.59,0.55,0.90,0.85);
	lrghvolt->SetBorderSize(1);
	lrghvolt->SetFillColor(0);
	lrghvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "rghvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* rghvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		rghvoltage->SetMarkerStyle(temp3+20);
		rghvoltage->SetMarkerColor(temp2+1);
		rghvoltage->SetMarkerSize(2);
		rghvoltage->SetLineColor(temp2+1);
		rghvoltage->SetLineWidth(2);
		rghvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		rghvoltage->SetTitle(title);
		rghvoltage->Draw("LP");
		lrghvolt->AddEntry(rghvoltage,title,"lep");

	}

	lrghvolt->Draw();
	canvas_rghvoltage->Update();
	_outputFile->cd();
	canvas_rghvoltage->Write();
	canvas_rghvoltage->Close();

	delete canvas_rghvoltage;
	delete hrghvolt;
	delete lrghvolt;

	// done rgh vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Clustersize vs. Voltage Plot..." << endl;
	}

	// clustersize vs voltage
	TCanvas* canvas_clustersizevoltage = new TCanvas("canvas_clustersizevoltage","Cluster Size vs. Voltage",200,10,700,500);
	canvas_clustersizevoltage->cd();
	TH2F *hclustervolt = new TH2F("hclustervolt","Cluster Size vs. Voltage",20,0,1000,20,0,5);
	hclustervolt->SetXTitle("Bias Voltage [|V|]");
	hclustervolt->SetYTitle("Cluster Size [1]");
	hclustervolt->SetStats(0000);
	hclustervolt->Draw();

	TLegend *lclustervolt = new TLegend(0.59,0.55,0.90,0.85);
	lclustervolt->SetBorderSize(1);
	lclustervolt->SetFillColor(0);
	lclustervolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "clustersizevoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* clustersizevoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		clustersizevoltage->SetMarkerStyle(temp3+20);
		clustersizevoltage->SetMarkerColor(temp2+1);
		clustersizevoltage->SetMarkerSize(2);
		clustersizevoltage->SetLineColor(temp2+1);
		clustersizevoltage->SetLineWidth(2);
		clustersizevoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		clustersizevoltage->SetTitle(title);
		clustersizevoltage->Draw("LP");
		lclustervolt->AddEntry(clustersizevoltage,title,"lep");

	}

	lclustervolt->Draw();
	canvas_clustersizevoltage->Update();
	_outputFile->cd();
	canvas_clustersizevoltage->Write();
	canvas_clustersizevoltage->Close();

	delete canvas_clustersizevoltage;
	delete hclustervolt;
	delete lclustervolt;

	// done clustersize vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Matched Clusters vs. Voltage Plot..." << endl;
	}

	// clustercount vs voltage
	TCanvas* canvas_clustercountvoltage = new TCanvas("canvas_clustercountvoltage","Cluster Count vs. Voltage",200,10,700,500);
	canvas_clustercountvoltage->cd();
	TH2F *hclustercountvolt = new TH2F("hclustercountvolt","Matched Cluster Count vs. Voltage",20,0,1000,20,0,_cut_maxclustercount);
	hclustercountvolt->SetXTitle("Bias Voltage [|V|]");
	hclustercountvolt->SetYTitle("Matched Clusters per good Track [1]");
	hclustercountvolt->SetStats(0000);
	hclustercountvolt->Draw();

	TLegend *lclustercountvolt = new TLegend(0.59,0.55,0.90,0.85);
	lclustercountvolt->SetBorderSize(1);
	lclustercountvolt->SetFillColor(0);
	lclustercountvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "clustercountvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* clustercountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		clustercountvoltage->SetMarkerStyle(temp3+20);
		clustercountvoltage->SetMarkerColor(temp2+1);
		clustercountvoltage->SetMarkerSize(2);
		clustercountvoltage->SetLineColor(temp2+1);
		clustercountvoltage->SetLineWidth(2);
		clustercountvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		clustercountvoltage->SetTitle(title);
		clustercountvoltage->Draw("LP");
		lclustercountvolt->AddEntry(clustercountvoltage,title,"lep");

	}

	lclustercountvolt->Draw();
	canvas_clustercountvoltage->Update();
	_outputFile->cd();
	canvas_clustercountvoltage->Write();
	canvas_clustercountvoltage->Close();

	delete canvas_clustercountvoltage;
	delete lclustercountvolt;
	delete hclustercountvolt;

	// done clustercount vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Effective Channels vs. Voltage Plot..." << endl;
	}

	// clustercount vs voltage
	TCanvas* canvas_channelcountvoltage = new TCanvas("canvas_channelcountvoltage","Channel Count vs. Voltage",200,10,700,500);
	canvas_channelcountvoltage->cd();
	TH2F *hchannelcountvolt = new TH2F("hchannelcountvolt","Effective Channels in Signal vs. Voltage",20,0,1000,20,0,500);
	hchannelcountvolt->SetXTitle("Bias Voltage [|V|]");
	hchannelcountvolt->SetYTitle("Effective Channels in Signal [1]");
	hchannelcountvolt->SetStats(0000);
	hchannelcountvolt->Draw();

	TLegend *lchannelcountvolt = new TLegend(0.59,0.55,0.90,0.85);
	lchannelcountvolt->SetBorderSize(1);
	lchannelcountvolt->SetFillColor(0);
	lchannelcountvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "channelcountvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* channelcountvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		channelcountvoltage->SetMarkerStyle(temp3+20);
		channelcountvoltage->SetMarkerColor(temp2+1);
		channelcountvoltage->SetMarkerSize(2);
		channelcountvoltage->SetLineColor(temp2+1);
		channelcountvoltage->SetLineWidth(2);
		channelcountvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		channelcountvoltage->SetTitle(title);
		channelcountvoltage->Draw("LP");
		lchannelcountvolt->AddEntry(channelcountvoltage,title,"lep");

	}

	lchannelcountvolt->Draw();
	canvas_channelcountvoltage->Update();
	_outputFile->cd();
	canvas_channelcountvoltage->Write();
	canvas_channelcountvoltage->Close();

	delete canvas_channelcountvoltage;
	delete hchannelcountvolt;
	delete lchannelcountvolt;

	// done channelcount vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Current vs. Voltage Plot..." << endl;
	}

	// current vs voltage
	TCanvas* canvas_currentvoltage = new TCanvas("canvas_currentvoltage","Sensor Current vs. Voltage",200,10,700,500);
	canvas_currentvoltage->cd();
	TH2F *hcurrentvolt = new TH2F("hcurrentvolt","Sensor Current vs. Voltage",20,0,1000,20,0,1000);
	hcurrentvolt->SetXTitle("Bias Voltage [|V|]");
	hcurrentvolt->SetYTitle("Sensor Current (at 253K) [| #muA|]");
	hcurrentvolt->SetStats(0000);
	hcurrentvolt->Draw();

	TLegend *lcurrentcountvolt = new TLegend(0.59,0.55,0.90,0.85);
	lcurrentcountvolt->SetBorderSize(1);
	lcurrentcountvolt->SetFillColor(0);
	lcurrentcountvolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "currentvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* currentvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		currentvoltage->SetMarkerStyle(temp3+20);
		currentvoltage->SetMarkerColor(temp2+1);
		currentvoltage->SetMarkerSize(2);
		currentvoltage->SetLineColor(temp2+1);
		currentvoltage->SetLineWidth(2);
		currentvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);

		// skip rotated sensors, as they have the same current
		if (n < 10)
		{
			currentvoltage->SetTitle(title);
			currentvoltage->Draw("LP");
			lcurrentcountvolt->AddEntry(currentvoltage,title,"lep");
		}

	}

	lcurrentcountvolt->Draw();
	canvas_currentvoltage->Update();
	_outputFile->cd();
	canvas_currentvoltage->Write();
	canvas_currentvoltage->Close();

	delete canvas_currentvoltage;
	delete hcurrentvolt;
	delete lcurrentcountvolt;

	// done current vs voltage

	// current vs voltage
	TCanvas* canvas_volumecurrent = new TCanvas("canvas_volumecurrent","Sensor Current / Thickness vs. Voltage",200,10,700,500);
	canvas_volumecurrent->cd();
	TH2F *hvolumecurrent = new TH2F("hvolumecurrent","Sensor Current vs. Voltage",20,0,1000,20,0,100000);
	hvolumecurrent->SetXTitle("Bias Voltage [|V|]");
	hvolumecurrent->SetYTitle("Sensor Current (at 253K) [| #muA| / cm^{-3}]");
	hvolumecurrent->SetStats(0000);
	hvolumecurrent->Draw();

	TLegend *lvolumecurrent = new TLegend(0.59,0.55,0.90,0.85);
	lvolumecurrent->SetBorderSize(1);
	lvolumecurrent->SetFillColor(0);
	lvolumecurrent->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "volumecurrent_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* volumecurrent = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		volumecurrent->SetMarkerStyle(temp3+20);
		volumecurrent->SetMarkerColor(temp2+1);
		volumecurrent->SetMarkerSize(2);
		volumecurrent->SetLineColor(temp2+1);
		volumecurrent->SetLineWidth(2);
		volumecurrent->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);

		// skip rotated sensors, as they have the same current
		if (n < 10)
		{
			volumecurrent->SetTitle(title);
			volumecurrent->Draw("LP");
			lvolumecurrent->AddEntry(volumecurrent,title,"lep");
		}

	}

	lvolumecurrent->Draw();
	canvas_volumecurrent->Update();
	_outputFile->cd();
	canvas_volumecurrent->Write();
	canvas_volumecurrent->Close();

	delete canvas_volumecurrent;
	delete hvolumecurrent;
	delete lvolumecurrent;

	// done current vs voltage

	// current vs voltage
	TCanvas* canvas_surfacecurrent = new TCanvas("canvas_surfacecurrent","Sensor Current / Surface vs. Voltage",200,10,700,500);
	canvas_surfacecurrent->cd();
	TH2F *hsurfacecurrent = new TH2F("hsurfacecurrent","Sensor Current vs. Voltage",20,0,1000,20,0,1000);
	hsurfacecurrent->SetXTitle("Bias Voltage [|V|]");
	hsurfacecurrent->SetYTitle("Sensor Current (at 253K) [| #muA| / cm^{-2}]");
	hsurfacecurrent->SetStats(0000);
	hsurfacecurrent->Draw();

	TLegend *lsurfacecurrent = new TLegend(0.59,0.55,0.90,0.85);
	lsurfacecurrent->SetBorderSize(1);
	lsurfacecurrent->SetFillColor(0);
	lsurfacecurrent->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "surfacecurrent_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* surfacecurrent = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		surfacecurrent->SetMarkerStyle(temp3+20);
		surfacecurrent->SetMarkerColor(temp2+1);
		surfacecurrent->SetMarkerSize(2);
		surfacecurrent->SetLineColor(temp2+1);
		surfacecurrent->SetLineWidth(2);
		surfacecurrent->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);

		// skip rotated sensors, as they have the same current
		if (n < 10)
		{
			surfacecurrent->SetTitle(title);
			surfacecurrent->Draw("LP");
			lsurfacecurrent->AddEntry(surfacecurrent,title,"lep");
		}

	}

	lsurfacecurrent->Draw();
	canvas_surfacecurrent->Update();
	_outputFile->cd();
	canvas_surfacecurrent->Write();
	canvas_surfacecurrent->Close();

	delete canvas_surfacecurrent;
	delete hsurfacecurrent;
	delete lsurfacecurrent;

	// done current vs voltage

	if (_debug <= 1)
	{
		cout << " " << endl;
		cout << "Starting Charge Sharing vs. Voltage Plot..." << endl;
	}

	// chargesharing vs voltage
	_outputFile->cd();
	TDirectory* csharedirectory = _outputFile->mkdir("Charge Sharing Threshold");
	csharedirectory->cd();
	
	TCanvas* canvas_chargesharingvoltage = new TCanvas("canvas_chargesharingvoltage","Shared Charge vs. Voltage",200,10,700,500);
	canvas_chargesharingvoltage->cd();
	TH2F *hchargesharevolt = new TH2F("hchargesharevolt","Shared Charge vs. Voltage",20,0,1000,20,0,100);
	hchargesharevolt->SetXTitle("Bias Voltage [|V|]");
	hchargesharevolt->SetYTitle("Shared Charge [%]");
	hchargesharevolt->SetStats(0000);
	hchargesharevolt->Draw();

	TLegend *lchargesharevolt = new TLegend(0.59,0.55,0.90,0.85);
	lchargesharevolt->SetBorderSize(1);
	lchargesharevolt->SetFillColor(0);
	lchargesharevolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "chargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* chargesharingvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		chargesharingvoltage->SetMarkerStyle(temp3+20);
		chargesharingvoltage->SetMarkerColor(temp2+1);
		chargesharingvoltage->SetMarkerSize(2);
		chargesharingvoltage->SetLineColor(temp2+1);
		chargesharingvoltage->SetLineWidth(2);
		chargesharingvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		chargesharingvoltage->SetTitle(title);
		chargesharingvoltage->Draw("LP");
		lchargesharevolt->AddEntry(chargesharingvoltage,title,"lep");

	}

	lchargesharevolt->Draw();
	canvas_chargesharingvoltage->Update();
	csharedirectory->cd();
	canvas_chargesharingvoltage->Write();
	canvas_chargesharingvoltage->Close();

	delete canvas_chargesharingvoltage;
	delete hchargesharevolt;
	delete lchargesharevolt;

	TCanvas* canvas_chargesharingvoltage_0 = new TCanvas("canvas_chargesharingvoltage_0","Shared Charge vs. Voltage, 0 deg",200,10,700,500);
	canvas_chargesharingvoltage_0->cd();
	TH2F *hchargesharevolt_0 = new TH2F("hchargesharevolt_0","Shared Charge vs. Voltage, 0 deg",20,0,1000,20,0,100);
	hchargesharevolt_0->SetXTitle("Bias Voltage [|V|]");
	hchargesharevolt_0->SetYTitle("Shared Charge [%]");
	hchargesharevolt_0->SetStats(0000);
	hchargesharevolt_0->Draw();

	TLegend *lchargesharevolt_0 = new TLegend(0.59,0.55,0.90,0.85);
	lchargesharevolt_0->SetBorderSize(1);
	lchargesharevolt_0->SetFillColor(0);
	lchargesharevolt_0->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "chargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* chargesharingvoltage_0 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			chargesharingvoltage_0->SetMarkerStyle(temp3+20);
			chargesharingvoltage_0->SetMarkerColor(temp2+1);
			chargesharingvoltage_0->SetMarkerSize(2);
			chargesharingvoltage_0->SetLineColor(temp2+1);
			chargesharingvoltage_0->SetLineWidth(2);
			chargesharingvoltage_0->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			chargesharingvoltage_0->SetTitle(title);
			chargesharingvoltage_0->Draw("LP");
			lchargesharevolt_0->AddEntry(chargesharingvoltage_0,title,"lep");
		}

	}

	lchargesharevolt_0->Draw();
	canvas_chargesharingvoltage_0->Update();
	csharedirectory->cd();
	canvas_chargesharingvoltage_0->Write();
	canvas_chargesharingvoltage_0->Close();

	delete canvas_chargesharingvoltage_0;
	delete hchargesharevolt_0;
	delete lchargesharevolt_0;

	TCanvas* canvas_chargesharingvoltage_25 = new TCanvas("canvas_chargesharingvoltage_25","Shared Charge vs. Voltage, 25 deg",200,10,700,500);
	canvas_chargesharingvoltage_25->cd();
	TH2F *hchargesharevolt_25 = new TH2F("hchargesharevolt_25","Shared Charge vs. Voltage, 25 deg",20,0,1000,20,0,100);
	hchargesharevolt_25->SetXTitle("Bias Voltage [|V|]");
	hchargesharevolt_25->SetYTitle("Shared Charge [%]");
	hchargesharevolt_25->SetStats(0000);
	hchargesharevolt_25->Draw();

	TLegend *lchargesharevolt_25 = new TLegend(0.59,0.55,0.90,0.85);
	lchargesharevolt_25->SetBorderSize(1);
	lchargesharevolt_25->SetFillColor(0);
	lchargesharevolt_25->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "chargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* chargesharingvoltage_25 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			chargesharingvoltage_25->SetMarkerStyle(temp3+20);
			chargesharingvoltage_25->SetMarkerColor(temp2+1);
			chargesharingvoltage_25->SetMarkerSize(2);
			chargesharingvoltage_25->SetLineColor(temp2+1);
			chargesharingvoltage_25->SetLineWidth(2);
			chargesharingvoltage_25->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			chargesharingvoltage_25->SetTitle(title);
			chargesharingvoltage_25->Draw("LP");
			lchargesharevolt_25->AddEntry(chargesharingvoltage_25,title,"lep");
		}

	}

	lchargesharevolt_25->Draw();
	canvas_chargesharingvoltage_25->Update();
	csharedirectory->cd();
	canvas_chargesharingvoltage_25->Write();
	canvas_chargesharingvoltage_25->Close();

	delete canvas_chargesharingvoltage_25;
	delete hchargesharevolt_25;
	delete lchargesharevolt_25;

	TCanvas* canvas_chargesharingvoltage_51 = new TCanvas("canvas_chargesharingvoltage_51","Shared Charge vs. Voltage, other rot",200,10,700,500);
	canvas_chargesharingvoltage_51->cd();
	TH2F *hchargesharevolt_51 = new TH2F("hchargesharevolt_51","Shared Charge vs. Voltage, other rot",20,0,1000,20,0,100);
	hchargesharevolt_51->SetXTitle("Bias Voltage [|V|]");
	hchargesharevolt_51->SetYTitle("Shared Charge [%]");
	hchargesharevolt_51->SetStats(0000);
	hchargesharevolt_51->Draw();

	TLegend *lchargesharevolt_51 = new TLegend(0.59,0.55,0.90,0.85);
	lchargesharevolt_51->SetBorderSize(1);
	lchargesharevolt_51->SetFillColor(0);
	lchargesharevolt_51->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "chargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 30)
		{

			TGraphErrors* chargesharingvoltage_51 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			chargesharingvoltage_51->SetMarkerStyle(temp3+20);
			chargesharingvoltage_51->SetMarkerColor(temp2+1);
			chargesharingvoltage_51->SetMarkerSize(2);
			chargesharingvoltage_51->SetLineColor(temp2+1);
			chargesharingvoltage_51->SetLineWidth(2);
			chargesharingvoltage_51->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			chargesharingvoltage_51->SetTitle(title);
			chargesharingvoltage_51->Draw("LP");
			lchargesharevolt_51->AddEntry(chargesharingvoltage_51,title,"lep");
		}

	}

	lchargesharevolt_51->Draw();
	canvas_chargesharingvoltage_51->Update();
	csharedirectory->cd();
	canvas_chargesharingvoltage_51->Write();
	canvas_chargesharingvoltage_51->Close();

	delete canvas_chargesharingvoltage_51;
	delete hchargesharevolt_51;
	delete lchargesharevolt_51;

	// done chargesharing vs voltage

	// eta chargesharing vs voltage
	_outputFile->cd();
	TDirectory* etacsharedirectory = _outputFile->mkdir("Charge Sharing Eta");
	etacsharedirectory->cd();

	TCanvas* canvas_etachargesharingvoltage = new TCanvas("canvas_etachargesharingvoltage","Shared Charge vs. Voltage",200,10,700,500);
	canvas_etachargesharingvoltage->cd();
	TH2F *hetachargesharevolt = new TH2F("hetachargesharevolt","Shared Charge vs. Voltage",20,0,1000,20,0,100);
	hetachargesharevolt->SetXTitle("Bias Voltage [|V|]");
	hetachargesharevolt->SetYTitle("Shared Charge [%]");
	hetachargesharevolt->SetStats(0000);
	hetachargesharevolt->Draw();

	TLegend *letachargesharevolt = new TLegend(0.59,0.55,0.90,0.85);
	letachargesharevolt->SetBorderSize(1);
	letachargesharevolt->SetFillColor(0);
	letachargesharevolt->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "etachargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		TGraphErrors* etachargesharingvoltage = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
		sstream.str(string());
		etachargesharingvoltage->SetMarkerStyle(temp3+20);
		etachargesharingvoltage->SetMarkerColor(temp2+1);
		etachargesharingvoltage->SetMarkerSize(2);
		etachargesharingvoltage->SetLineColor(temp2+1);
		etachargesharingvoltage->SetLineWidth(2);
		etachargesharingvoltage->SetLineStyle(temp+1);

		char o [100];
		char title[300];
		sprintf(o,"Epi");
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else if (k.at(3) == 'A') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 40 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'B') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 80 min @ 60#circC");
				sstream << o;
			} else if (k.at(3) == 'C') {
				sstream.str(string());
				sprintf(o,"Mcz200P, 156 min @ 60#circC");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}
		string sensname = sstream.str();
		sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
		etachargesharingvoltage->SetTitle(title);
		etachargesharingvoltage->Draw("LP");
		letachargesharevolt->AddEntry(etachargesharingvoltage,title,"lep");

	}

	letachargesharevolt->Draw();
	canvas_etachargesharingvoltage->Update();
	etacsharedirectory->cd();
	canvas_etachargesharingvoltage->Write();
	canvas_etachargesharingvoltage->Close();

	delete canvas_etachargesharingvoltage;
	delete hetachargesharevolt;
	delete letachargesharevolt;

	TCanvas* canvas_etachargesharingvoltage_0 = new TCanvas("canvas_etachargesharingvoltage_0","Shared Charge vs. Voltage, 0 deg",200,10,700,500);
	canvas_etachargesharingvoltage_0->cd();
	TH2F *hetachargesharevolt_0 = new TH2F("hetachargesharevolt_0","Shared Charge vs. Voltage, 0 deg",20,0,1000,20,0,100);
	hetachargesharevolt_0->SetXTitle("Bias Voltage [|V|]");
	hetachargesharevolt_0->SetYTitle("Shared Charge [%]");
	hetachargesharevolt_0->SetStats(0000);
	hetachargesharevolt_0->Draw();

	TLegend *letachargesharevolt_0 = new TLegend(0.59,0.55,0.90,0.85);
	letachargesharevolt_0->SetBorderSize(1);
	letachargesharevolt_0->SetFillColor(0);
	letachargesharevolt_0->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "etachargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n < 20)
		{

			TGraphErrors* etachargesharingvoltage_0 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			etachargesharingvoltage_0->SetMarkerStyle(temp3+20);
			etachargesharingvoltage_0->SetMarkerColor(temp2+1);
			etachargesharingvoltage_0->SetMarkerSize(2);
			etachargesharingvoltage_0->SetLineColor(temp2+1);
			etachargesharingvoltage_0->SetLineWidth(2);
			etachargesharingvoltage_0->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			etachargesharingvoltage_0->SetTitle(title);
			etachargesharingvoltage_0->Draw("LP");
			letachargesharevolt_0->AddEntry(etachargesharingvoltage_0,title,"lep");
		}

	}

	letachargesharevolt_0->Draw();
	canvas_etachargesharingvoltage_0->Update();
	etacsharedirectory->cd();
	canvas_etachargesharingvoltage_0->Write();
	canvas_etachargesharingvoltage_0->Close();

	delete canvas_etachargesharingvoltage_0;
	delete hetachargesharevolt_0;
	delete letachargesharevolt_0;

	TCanvas* canvas_etachargesharingvoltage_25 = new TCanvas("canvas_etachargesharingvoltage_25","Shared Charge vs. Voltage, 25 deg",200,10,700,500);
	canvas_etachargesharingvoltage_25->cd();
	TH2F *hetachargesharevolt_25 = new TH2F("hetachargesharevolt_25","Shared Charge vs. Voltage, 25 deg",20,0,1000,20,0,100);
	hetachargesharevolt_25->SetXTitle("Bias Voltage [|V|]");
	hetachargesharevolt_25->SetYTitle("Shared Charge [%]");
	hetachargesharevolt_25->SetStats(0000);
	hetachargesharevolt_25->Draw();

	TLegend *letachargesharevolt_25 = new TLegend(0.59,0.55,0.90,0.85);
	letachargesharevolt_25->SetBorderSize(1);
	letachargesharevolt_25->SetFillColor(0);
	letachargesharevolt_25->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "etachargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 20 && n < 30)
		{

			TGraphErrors* etachargesharingvoltage_25 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			etachargesharingvoltage_25->SetMarkerStyle(temp3+20);
			etachargesharingvoltage_25->SetMarkerColor(temp2+1);
			etachargesharingvoltage_25->SetMarkerSize(2);
			etachargesharingvoltage_25->SetLineColor(temp2+1);
			etachargesharingvoltage_25->SetLineWidth(2);
			etachargesharingvoltage_25->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			etachargesharingvoltage_25->SetTitle(title);
			etachargesharingvoltage_25->Draw("LP");
			letachargesharevolt_25->AddEntry(etachargesharingvoltage_25,title,"lep");
		}

	}

	letachargesharevolt_25->Draw();
	canvas_etachargesharingvoltage_25->Update();
	etacsharedirectory->cd();
	canvas_etachargesharingvoltage_25->Write();
	canvas_etachargesharingvoltage_25->Close();

	delete canvas_etachargesharingvoltage_25;
	delete hetachargesharevolt_25;
	delete letachargesharevolt_25;

	TCanvas* canvas_etachargesharingvoltage_51 = new TCanvas("canvas_etachargesharingvoltage_51","Shared Charge vs. Voltage, other rot",200,10,700,500);
	canvas_etachargesharingvoltage_51->cd();
	TH2F *hetachargesharevolt_51 = new TH2F("hetachargesharevolt_51","Shared Charge vs. Voltage, other rot",20,0,1000,20,0,100);
	hetachargesharevolt_51->SetXTitle("Bias Voltage [|V|]");
	hetachargesharevolt_51->SetYTitle("Shared Charge [%]");
	hetachargesharevolt_51->SetStats(0000);
	hetachargesharevolt_51->Draw();

	TLegend *letachargesharevolt_51 = new TLegend(0.59,0.55,0.90,0.85);
	letachargesharevolt_51->SetBorderSize(1);
	letachargesharevolt_51->SetFillColor(0);
	letachargesharevolt_51->SetFillStyle(0);

	for (unsigned int i=0;i<voltagegraphs.size();i++)
	{

		sstream.str(string());
		histoname = "etachargesharingvoltage_";
		sstream << histoname << voltagegraphs.at(i);
		int temp = voltagegraphs.at(i) % _dutrotation_total;
		int n = _dutrotation_list.at(temp);
		int temp2 = (voltagegraphs.at(i) / _dutrotation_total) % _irradfluence_total;

		double l = _irradfluence_list.at(temp2);
		int temp3 = (voltagegraphs.at(i) / (_irradfluence_total*_dutrotation_total)) % _sensorname_total;
		string k = _sensorname_list.at(temp3).c_str();

		// jump yellow:
		if (temp2==4)
		{
			temp2++; 
		}

		if (n >= 30)
		{

			TGraphErrors* etachargesharingvoltage_51 = dynamic_cast<TGraphErrors*> ( _rootObjectMap[sstream.str()]);
			sstream.str(string());
			etachargesharingvoltage_51->SetMarkerStyle(temp3+20);
			etachargesharingvoltage_51->SetMarkerColor(temp2+1);
			etachargesharingvoltage_51->SetMarkerSize(2);
			etachargesharingvoltage_51->SetLineColor(temp2+1);
			etachargesharingvoltage_51->SetLineWidth(2);
			etachargesharingvoltage_51->SetLineStyle(temp+1);

			char o [100];
			char title[300];
			sprintf(o,"Epi");
			sstream << o << k;

			if (k.at(0) == '2')
			{
				if (k.at(3) == 'F')
				{
					sstream.str(string());
					sprintf(o,"Fth200Y");
					sstream << o;
				} else if (k.at(3) == 'A') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 40 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'B') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 80 min @ 60#circC");
					sstream << o;
				} else if (k.at(3) == 'C') {
					sstream.str(string());
					sprintf(o,"Mcz200P, 156 min @ 60#circC");
					sstream << o;
				} else {
					sstream.str(string());
					sprintf(o,"Mcz");
					sstream << o << k;
				}
			}
			string sensname = sstream.str();
			sprintf(title, "%s, F=%.1e, %irot", sensname.c_str(),l,n);
			etachargesharingvoltage_51->SetTitle(title);
			etachargesharingvoltage_51->Draw("LP");
			letachargesharevolt_51->AddEntry(etachargesharingvoltage_51,title,"lep");
		}

	}

	letachargesharevolt_51->Draw();
	canvas_etachargesharingvoltage_51->Update();
	etacsharedirectory->cd();
	canvas_etachargesharingvoltage_51->Write();
	canvas_etachargesharingvoltage_51->Close();

	delete canvas_etachargesharingvoltage_51;
	delete hetachargesharevolt_51;
	delete letachargesharevolt_51;


	// write out the alignment histos
	_outputFile->cd();
	TDirectory* positiondirectory = _outputFile->mkdir("Sensor Positions");
	positiondirectory->cd();
	xshift->Write();
	yshift->Write();
	zshift->Write();
	ashift->Write();
	bshift->Write();
	cshift->Write();

	// done all multi-run things

} // done read and fill function


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// this books all the histograms needed
void bookhistos()
{

	cout << "##############################################" << endl;
	cout << " " << endl;
	cout << "Booking all the histograms!" << endl;
	cout << " " << endl;

	// the names and titles for the histograms
	char name[100];
	char title[300];
	stringstream sstream;

	// we have _runcount many runs
	// details further down on which histo does what
	std::vector<std::vector<TH1D*>> noise;
	std::vector<std::vector<TH1D*>> noiseelec;
	std::vector<TH1D*> allnoise;
	std::vector<TH1D*> allnoiseelec;
	std::vector<TH1D*> adcnoise;
	std::vector<TH1D*> adcnoiseelec;
	std::vector<TH1D*> fivechannoise;
	std::vector<TH1D*> fivechannoiseelec;

	std::vector<TH1D*> tempdistri;
	std::vector<TH1D*> tdcdistri;

	std::vector<std::vector<TH1D*>> tdceval;
	std::vector<std::vector<TH1D*>> tdcevalelec;
	std::vector<TGraphErrors*> overalltdc;

	std::vector<TGraphErrors*> overallxpos;
	std::vector<std::vector<TH1D*>> posxeval;

	std::vector<TH1D*> residualsX;
	std::vector<TH1D*> residualsY;
	std::vector<TH1D*> residualsZ;
	std::vector<TH2D*> residualsXY;
	std::vector<TH2D*> residualsXZ;
	std::vector<TH2D*> residualsYZ;
	std::vector<TH2D*> residualsYQ;
	std::vector<TH2D*> residualsYR;
	std::vector<TH2D*> residualsYT;
	std::vector<TH2D*> residualsY_Q;

	std::vector<TH1D*> residualsY_clu1;
	std::vector<TH1D*> residualsY_clu2;
	std::vector<TH1D*> residualsY_clu3;
	std::vector<TH1D*> residualsY_clu4;

	std::vector<TH2D*> residualsDXvsX;
	std::vector<TH2D*> residualsDXvsY;
	std::vector<TH2D*> residualsDXvsZ;
	std::vector<TH2D*> residualsDYvsX;
	std::vector<TH2D*> residualsDYvsY;
	std::vector<TH2D*> residualsDYvsZ;
	std::vector<TH2D*> residualsDZvsX;
	std::vector<TH2D*> residualsDZvsY;
	std::vector<TH2D*> residualsDZvsZ;

	std::vector<TH2D*> residualsDXvsSXU;
	std::vector<TH2D*> residualsDXvsSYU;
	std::vector<TH2D*> residualsDYvsSXU;
	std::vector<TH2D*> residualsDYvsSYU;
	std::vector<TH2D*> residualsDZvsSXU;
	std::vector<TH2D*> residualsDZvsSYU;

	std::vector<TH2D*> residualsDXvsSXD;
	std::vector<TH2D*> residualsDXvsSYD;
	std::vector<TH2D*> residualsDYvsSXD;
	std::vector<TH2D*> residualsDYvsSYD;
	std::vector<TH2D*> residualsDZvsSXD;
	std::vector<TH2D*> residualsDZvsSYD;

	std::vector<TH2D*> residualsXvsEvt;
	std::vector<TH2D*> residualsYvsEvt;
	std::vector<TH2D*> residualsZvsEvt;

	std::vector<TH3D*> residualMapTemp;

	std::vector<TProfile*> residualProfileDXvsX;
	std::vector<TProfile*> residualProfileDXvsY;
	std::vector<TProfile*> residualProfileDYvsX;
	std::vector<TProfile*> residualProfileDYvsY;

	std::vector<TH2D*> hitmapDUT;
	std::vector<TH2D*> hitmapDUTmodpitch;

	std::vector<TH2D*> hitmapTracks;
	std::vector<TH2D*> hitmapTracksmodpitch;

	std::vector<TH2D*> hitmapMatch;
	std::vector<TH2D*> hitmapMatchmodpitch;
	std::vector<TH2D*> hitmapMatchTracks;

	std::vector<TProfile2D*> scatterX;
	std::vector<TProfile2D*> scatterY;
	std::vector<TProfile2D*> scatterXY;

	std::vector<TH1D*> kinkX;
	std::vector<TH1D*> kinkY;

	std::vector<TH2D*> hitmapChargeShared0;
	std::vector<TH2D*> hitmapChargeShared1;
	std::vector<TH2D*> hitmapChargeShared2;
	std::vector<TH2D*> hitmapChargeShared3;
	std::vector<TH2D*> hitmapChargeShared4;

	std::vector<TH1D*> hitmapClusterSizeA;
	std::vector<TH1D*> hitmapClusterSize1;
	std::vector<TH1D*> hitmapClusterSize2;
	std::vector<TH1D*> hitmapClusterSize3;
	std::vector<TH1D*> hitmapClusterSize4;

	std::vector<TH1D*> signalClusterSize1;
	std::vector<TH1D*> signalClusterSize2;
	std::vector<TH1D*> signalClusterSize3;
	std::vector<TH1D*> signalClusterSize4;
	std::vector<TH1D*> signalClusterSize1elec;
	std::vector<TH1D*> signalClusterSize2elec;
	std::vector<TH1D*> signalClusterSize3elec;
	std::vector<TH1D*> signalClusterSize4elec;

	std::vector<TH1D*> matchedEta;
	std::vector<TH2D*> matchedEta2D;
	std::vector<TH1D*> unmatchedEta;
	std::vector<TH2D*> unmatchedEta2D;
	std::vector<TH1D*> matchedEtaIntegral;
	std::vector<TH1D*> unmatchedEtaIntegral;
	std::vector<TH1D*> striphitEta;
	std::vector<TH2D*> striphitEta2D;
	std::vector<TH1D*> striphitEtaIntegral;
	std::vector<TH1D*> noiseEta;
	std::vector<TH1D*> noisesubtractedEta;
	std::vector<TH1D*> EtaL50;
	std::vector<TH1D*> EtaL40;
	std::vector<TH1D*> EtaL30;
	std::vector<TH1D*> EtaL20;
	std::vector<TH1D*> EtaL10;
	std::vector<TH1D*> EtaR10;
	std::vector<TH1D*> EtaR20;
	std::vector<TH1D*> EtaR30;
	std::vector<TH1D*> EtaR40;
	std::vector<TH1D*> EtaR50;

	std::vector<TH1D*> trackSignal;
	std::vector<TH2D*> trackSignalMap;
	std::vector<TH2D*> trackSignalTDC;
	std::vector<TH1D*> trackSignalelec;
	std::vector<TH2D*> trackSignalMapelec;
	std::vector<TH2D*> trackSignalTDCelec;
	std::vector<TH1D*> tracksperevent;
	std::vector<TH2D*> tracksvsevents;
	std::vector<TH1D*> highchanneldistribution;

	std::vector<TH1D*> signalLeft2;
	std::vector<TH1D*> signalLeft1;
	std::vector<TH1D*> signalCenter;
	std::vector<TH1D*> signalRight1;
	std::vector<TH1D*> signalRight2;
	std::vector<TH1D*> signalGoodEvents;
	std::vector<TH1D*> signalmapA;
	std::vector<TH1D*> signalmapB;
	std::vector<TH1D*> signalmapC;
	std::vector<TH1D*> signalmapD;

	std::vector<TH1D*> signalLeft2elec;
	std::vector<TH1D*> signalLeft1elec;
	std::vector<TH1D*> signalCenterelec;
	std::vector<TH1D*> signalRight1elec;
	std::vector<TH1D*> signalRight2elec;
	std::vector<TH1D*> signalGoodEventselec;
	std::vector<TH1D*> signalmapAelec;
	std::vector<TH1D*> signalmapBelec;
	std::vector<TH1D*> signalmapCelec;
	std::vector<TH1D*> signalmapDelec;

	std::vector<TGraphErrors*> signalareaplot;
	std::vector<TGraphErrors*> signalareaplotelec;

	std::vector<TH1D*> fiducial_discard;
	std::vector<TH1D*> fiducial_allow;
	std::vector<TH1D*> goodchannel_discard;
	std::vector<TH1D*> goodchannel_allow;
	std::vector<TH1D*> nohittrack;
	std::vector<TH1D*> goodevent_discard;
	std::vector<TH1D*> goodevent_allow;
	std::vector<TH1D*> trackselection_discard;
	std::vector<TH1D*> trackselection_allow;
	std::vector<TH1D*> timecut_discard;
	std::vector<TH1D*> timecut_allow;
	std::vector<TH1D*> highchannel_discard;
	std::vector<TH1D*> highchannel_allow;

	// sort by voltage:
	std::vector<TGraphErrors*> residualsXvoltage;
	std::vector<TGraphErrors*> residualsYvoltage;
	std::vector<TGraphErrors*> residualsYvoltageAngle;

	std::vector<TGraphErrors*> signalvoltage;
	std::vector<TGraphErrors*> signalvoltageelec;
	std::vector<TGraphErrors*> snvoltage;

	std::vector<TGraphErrors*> noisevoltage;
	std::vector<TGraphErrors*> noisevoltageelec;

	std::vector<TGraphErrors*> rghvoltage;

	std::vector<TGraphErrors*> clustersizevoltage;

	std::vector<TGraphErrors*> clustercountvoltage;

	std::vector<TGraphErrors*> channelcountvoltage;

	std::vector<TGraphErrors*> currentvoltage;

	std::vector<TGraphErrors*> volumecurrent;

	std::vector<TGraphErrors*> surfacecurrent;

	std::vector<TGraphErrors*> chargesharingvoltage;

	std::vector<TGraphErrors*> etachargesharingvoltage;

	std::vector<TGraphErrors*> signaldistancevoltage;
	std::vector<TGraphErrors*> signaldistancevoltageelec;

	std::vector<TH2D*> signalareamap;

	std::vector<TH3D*> signalMapTemp;

	// only once:
	TH1D* xshift = new TH1D("xshift","Sensor Offset in X;Position [mm];Entries", 100, -5.0, 5.0);
	_rootObjectMap["xshift"] = xshift;
	TH1D* yshift = new TH1D("yshift","Sensor Offset in Y;Position [mm];Entries", 100, -5.0, 5.0);
	_rootObjectMap["yshift"] = yshift;
	TH1D* zshift = new TH1D("zshift","Sensor Offset in Z;Position [mm];Entries", 100, -10.0, 10.0);
	_rootObjectMap["zshift"] = zshift;
	TH1D* ashift = new TH1D("ashift","Sensor Offset in A;Rotation [#circ];Entries", 100, -5.0, 5.0);
	_rootObjectMap["ashift"] = ashift;
	TH1D* bshift = new TH1D("bshift","Sensor Offset in B;Rotation [#circ];Entries", 100, -5.0, 5.0);
	_rootObjectMap["bshift"] = bshift;
	TH1D* cshift = new TH1D("cshift","Sensor Offset in C;Rotation [#circ];Entries", 100, -5.0, 5.0);
	_rootObjectMap["cshift"] = cshift;



	// j is the loop var, i is the actual run number
	for ( int j = 0 ; j < _actual_runcount; j++ )
	{

		// the sensor information goes into the histogram title
		int i = _runnumber.at(j);
		double l = _irradfluence.at(j);
		int m = _biasvoltage.at(j);
		int n = _dutrotation.at(j);
		char o[100] = "Fail";
		int z = _thickness.at(j);
		string k = _sensorname.at(j).c_str();

		// slight hack: material is deducted from thickness... FIXME
		if (z <=150)
		{
			sprintf(o,"Epi");
		}
		/*
		if (z > 150)
		{
			// fth sensors from runnumber
			if (i>=638 && i<=644)
			{
				sprintf(o,"Fth");
			} else {
				sprintf(o,"MCz");
			}
		}
		*/
		stringstream sstream;
		sstream << o << k;

		if (k.at(0) == '2')
		{
			if (k.at(3) == 'F')
			{
				sstream.str(string());
				sprintf(o,"Fth200Y");
				sstream << o;
			} else {
				sstream.str(string());
				sprintf(o,"Mcz");
				sstream << o << k;
			}
		}

		// the sensor name for the title

		string sensname = sstream.str();

		std::vector<TH1D*> tempvec;
		// noise * 128 channels
		for (int ii=0;ii<128;ii++)
		{
			sprintf(name, "noise_%i_chan_%i", i,ii);
			sprintf(title, "DUT Off-Beam Noise, Channel %i, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", ii, sensname.c_str(),l,m,n,i);
			tempvec.push_back(new TH1D(name, title, 1000, -500, 500));
		}
		noise.push_back(tempvec);
		for (int ii=0;ii<128;ii++)
		{
			sprintf(name, "noise_%i_chan_%i", i,ii);
			_rootObjectMap[name] = noise[j][ii] ;
		}

		std::vector<TH1D*> tempvecelec;
		// noise * 128 channels
		for (int ii=0;ii<128;ii++)
		{
			sprintf(name, "noiseelec_%i_chan_%i", i,ii);
			sprintf(title, "DUT Off-Beam Electron Noise, Channel %i, %s, F=%.1e, %iV, %irot, run %i;Signal [Electrons];Events", ii, sensname.c_str(),l,m,n,i);
			tempvecelec.push_back(new TH1D(name, title, 1000, -75000, 75000));
		}
		noiseelec.push_back(tempvecelec);
		for (int ii=0;ii<128;ii++)
		{
			sprintf(name, "noiseelec_%i_chan_%i", i,ii);
			_rootObjectMap[name] = noiseelec[j][ii] ;
		}

		std::vector<TH1D*> tempvec2;
		// tdceval * 91 bins
		for (int ii=0;ii<91;ii++)
		{
			sprintf(name, "tdceval_%i_tdc_%i", i,ii);
			sprintf(title, "Signal with TDC %i - %i, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", ii, (ii+10), sensname.c_str(),l,m,n,i);
			tempvec2.push_back(new TH1D(name, title, 1000, -500, 500));
		}
		tdceval.push_back(tempvec2);
		for (int ii=0;ii<91;ii++)
		{
			sprintf(name, "tdceval_%i_tdc_%i", i,ii);
			_rootObjectMap[name] = tdceval[j][ii] ;
		}

		std::vector<TH1D*> tempvec3;
		// posxeval * 100 bins
		for (int ii=0;ii<100;ii++)
		{
			sprintf(name, "posxeval_%i_pos_%i", i,ii);
			sprintf(title, "Signal with position %.2f - %.2f, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", (-10.0+(ii*0.25)), (-10.0+(ii*0.25)+_xevalstep), sensname.c_str(),l,m,n,i);
			tempvec3.push_back(new TH1D(name, title, 1000, -500, 500));
		}
		posxeval.push_back(tempvec3);
		for (int ii=0;ii<100;ii++)
		{
			sprintf(name, "posxeval_%i_pos_%i", i,ii);
			_rootObjectMap[name] = posxeval[j][ii] ;
		}

		// tdc of a sensor
		sprintf(name, "overalltdc_%i", i);
		overalltdc.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = overalltdc[j] ;

		sprintf(name, "overallxpos_%i", i);
		overallxpos.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = overallxpos[j] ;

		sprintf(name, "tempdistri_%i", i);
		sprintf(title, "DUT Chip Temperature, %s, F=%.1e, %iV, %irot, run %i;Temperature [#circC];Entries", sensname.c_str(),l,m,n,i);
		tempdistri.push_back(new TH1D(name, title, 100, -10, 30));
		_rootObjectMap[name] = tempdistri[j] ;

		sprintf(name, "tdcdistri_%i", i);
		sprintf(title, "DUT TDC Distribution, %s, F=%.1e, %iV, %irot, run %i;TDC [ns];Entries", sensname.c_str(),l,m,n,i);
		tdcdistri.push_back(new TH1D(name, title, 100, 0, 100));
		_rootObjectMap[name] = tdcdistri[j] ;

		// noise of a sensor
		sprintf(name, "allnoise_%i", i);
		sprintf(title, "DUT Off-Beam Noise, %s, F=%.1e, %iV, %irot, run %i;Channel;Noise [ADCs]", sensname.c_str(),l,m,n,i);
		allnoise.push_back(new TH1D(name, title, 128, 0, 127));
		_rootObjectMap[name] = allnoise[j] ;

		sprintf(name, "allnoiseelec_%i", i);
		sprintf(title, "DUT Off-Beam Noise, Electrons, %s, F=%.1e, %iV, %irot, run %i;Channel;Noise [Electrons]", sensname.c_str(),l,m,n,i);
		allnoiseelec.push_back(new TH1D(name, title, 128, 0, 127));
		_rootObjectMap[name] = allnoiseelec[j] ;

		// noise of a sensor, all channels
		sprintf(name, "adcnoise_%i", i);
		sprintf(title, "DUT Off-Beam Noise, %s, F=%.1e, %iV, %irot, run %i;Noise [ADCs];Events", sensname.c_str(),l,m,n,i);
		adcnoise.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = adcnoise[j] ;

		sprintf(name, "adcnoiseelec_%i", i);
		sprintf(title, "DUT Off-Beam Noise, Electrons, %s, F=%.1e, %iV, %irot, run %i;Noise [Electrons];Events", sensname.c_str(),l,m,n,i);
		adcnoiseelec.push_back(new TH1D(name, title, 1000, -75000, 75000));
		_rootObjectMap[name] = adcnoiseelec[j] ;

		// noise of a sensor, 5 channels summed
		sprintf(name, "fivechannoise_%i", i);
		sprintf(title, "DUT Off-Beam Noise in 5 Channels, %s, F=%.1e, %iV, %irot, run %i;Noise [ADCs];Events", sensname.c_str(),l,m,n,i);
		fivechannoise.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = fivechannoise[j];

		sprintf(name, "fivechannoiseelec_%i", i);
		sprintf(title, "DUT Off-Beam Noise in 5 Channels, Electrons, %s, F=%.1e, %iV, %irot, run %i;Noise [Electrons];Events", sensname.c_str(),l,m,n,i);
		fivechannoiseelec.push_back(new TH1D(name, title, 1000, -75000, 75000));
		_rootObjectMap[name] = fivechannoiseelec[j];

		// residuals
		sprintf(name, "residualsX_%i", i);
		sprintf(title, "DUT Residual in X, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsX.push_back(new TH1D(name, title, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsX[j] ;

		sprintf(name, "residualsY_%i", i);
		sprintf(title, "DUT Residual in Y, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsY.push_back(new TH1D(name, title, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsY[j] ;

		sprintf(name, "residualsZ_%i", i);
		sprintf(title, "DUT Residual in Z, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsZ.push_back(new TH1D(name, title, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsZ[j] ;

		sprintf(name, "residualsXY_%i", i);
		sprintf(title, "DUT Residual in XY, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsXY.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsXY[j] ;

		sprintf(name, "residualsXZ_%i", i);
		sprintf(title, "DUT Residual in XZ, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Z_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsXZ.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsXZ[j] ;

		sprintf(name, "residualsYZ_%i", i);
		sprintf(title, "DUT Residual in YZ, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Z_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsYZ.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsYZ[j] ;

		sprintf(name, "residualsYQ_%i", i);
		sprintf(title, "DUT Residual in Y vs. Cluster Charge, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Q_{hit} [ADCs];Events", sensname.c_str(),l,m,n,i);
		residualsYQ.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, 0, 100));
		_rootObjectMap[name] = residualsYQ[j] ;

		sprintf(name, "residualsYR_%i", i);
		sprintf(title, "DUT Residual in Y vs. Cluster Size, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Cluster Size [1];Events", sensname.c_str(),l,m,n,i);
		residualsYR.push_back(new TH2D(name, title, 500, -0.5, 0.5, 10, 0, 10));
		_rootObjectMap[name] = residualsYR[j] ;

		sprintf(name, "residualsYT_%i", i);
		sprintf(title, "DUT Residual in Y vs. TDC Time, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];TDC Time [ns];Events", sensname.c_str(),l,m,n,i);
		residualsYT.push_back(new TH2D(name, title, 500, -0.5, 0.5, 100, 0, 100));
		_rootObjectMap[name] = residualsYT[j] ;

		sprintf(name, "residualsY_Q_%i", i);
		sprintf(title, "DUT Residual in Y vs. Hit Charge, %s, F=%.1e, %iV, %irot, run %i;Q_{hit} [ADCs];Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsY_Q.push_back(new TH2D(name, title, 200, -100, 100, 50, -0.5, 0.5));
		_rootObjectMap[name] = residualsY_Q[j] ;

		sprintf(name, "residualsY_clu1_%i", i);
		sprintf(title, "DUT Residual in Y from Clustersize 1, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsY_clu1.push_back(new TH1D(name, title, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsY_clu1[j] ;

		sprintf(name, "residualsY_clu2_%i", i);
		sprintf(title, "DUT Residual in Y from Clustersize 2, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsY_clu2.push_back(new TH1D(name, title, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsY_clu2[j] ;

		sprintf(name, "residualsY_clu3_%i", i);
		sprintf(title, "DUT Residual in Y from Clustersize 3, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsY_clu3.push_back(new TH1D(name, title, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsY_clu3[j] ;

		sprintf(name, "residualsY_clu4_%i", i);
		sprintf(title, "DUT Residual in Y from Clustersize 4, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsY_clu4.push_back(new TH1D(name, title, 500, -0.5, 0.5));
		_rootObjectMap[name] = residualsY_clu4[j] ;

		// residuals vs xy
		sprintf(name, "residualsDXvsX_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track X, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];X_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDXvsX.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDXvsX[j] ;

		sprintf(name, "residualsDXvsY_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track Y, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDXvsY.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDXvsY[j] ;

		sprintf(name, "residualsDXvsZ_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track Z, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Z_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDXvsZ.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDXvsZ[j] ;

		sprintf(name, "residualsDYvsX_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track X, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];X_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDYvsX.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDYvsX[j] ;

		sprintf(name, "residualsDYvsY_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track Y, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDYvsY.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDYvsY[j] ;

		sprintf(name, "residualsDYvsZ_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track Z, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Z_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDYvsZ.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDYvsZ[j] ;

		sprintf(name, "residualsDZvsX_%i" ,i);
		sprintf(title, "DUT Residual in Z vs. Track X, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];X_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDZvsX.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDZvsX[j] ;

		sprintf(name, "residualsDZvsY_%i" ,i);
		sprintf(title, "DUT Residual in Z vs. Track Y, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDZvsY.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDZvsY[j] ;

		sprintf(name, "residualsDZvsZ_%i" ,i);
		sprintf(title, "DUT Residual in Z vs. Track Z, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];Z_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsDZvsZ.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -10.0, 10.0));
		_rootObjectMap[name] = residualsDZvsZ[j] ;

		sprintf(name, "residualsDXvsSXU_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track Slope X Upstream, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];X_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDXvsSXU.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDXvsSXU[j] ;

		sprintf(name, "residualsDXvsSYU_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track Slope Y Upstream, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Y_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDXvsSYU.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDXvsSYU[j] ;

		sprintf(name, "residualsDYvsSXU_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track Slope X Upstream, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];X_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDYvsSXU.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDYvsSXU[j] ;

		sprintf(name, "residualsDYvsSYU_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track Slope Y Upstream, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Y_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDYvsSYU.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDYvsSYU[j] ;

		sprintf(name, "residualsDZvsSXU_%i" ,i);
		sprintf(title, "DUT Residual in Z vs. Track Slope X Upstream, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];X_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDZvsSXU.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDZvsSXU[j] ;

		sprintf(name, "residualsDZvsSYU_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track Slope Y Upstream, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];Y_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDZvsSYU.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDZvsSYU[j] ;

		sprintf(name, "residualsDXvsSXD_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track Slope X Downstream, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];X_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDXvsSXD.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDXvsSXD[j] ;

		sprintf(name, "residualsDXvsSYD_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track Slope Y Downstream, %s, F=%.1e, %iV, %irot, run %i;X_{hit - track} [mm];Y_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDXvsSYD.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDXvsSYD[j] ;

		sprintf(name, "residualsDYvsSXD_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track Slope X Downstream, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];X_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDYvsSXD.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDYvsSXD[j] ;

		sprintf(name, "residualsDYvsSYD_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Track Slope Y Downstream, %s, F=%.1e, %iV, %irot, run %i;Y_{hit - track} [mm];Y_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDYvsSYD.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDYvsSYD[j] ;

		sprintf(name, "residualsDZvsSXD_%i" ,i);
		sprintf(title, "DUT Residual in Z vs. Track Slope X Downstream, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];X_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDZvsSXD.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDZvsSXD[j] ;

		sprintf(name, "residualsDZvsSYD_%i" ,i);
		sprintf(title, "DUT Residual in Z vs. Track Slope Y Downstream, %s, F=%.1e, %iV, %irot, run %i;Z_{hit - track} [mm];Y_{slope} [mrad];Events", sensname.c_str(),l,m,n,i);
		residualsDZvsSYD.push_back(new TH2D(name, title, 500, -0.5, 0.5, 500, -5.0, 5.0));
		_rootObjectMap[name] = residualsDZvsSYD[j] ;

		// residual profiles
		sprintf(name, "residualProfileDXvsX_%i", i);
		sprintf(title, "DUT Residual in X vs. Track X - Profile, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];X_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
		residualProfileDXvsX.push_back(new TProfile (name, title,150,-15,15,-0.2,0.2,"s"));
		_rootObjectMap[name] = residualProfileDXvsX[j] ;

		sprintf(name, "residualProfileDXvsY_%i", i);
		sprintf(title, "DUT Residual in X vs. Track Y - Profile, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];X_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
		residualProfileDXvsY.push_back(new TProfile (name, title,150,-15,15,-0.2,0.2,"s"));
		_rootObjectMap[name] = residualProfileDXvsY[j] ;

		sprintf(name, "residualProfileDYvsX_%i", i);
		sprintf(title, "DUT Residual in Y vs. Track X - Profile, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
		residualProfileDYvsX.push_back(new TProfile (name, title,150,-15,15,-0.2,0.2,"s"));
		_rootObjectMap[name] = residualProfileDYvsX[j] ;

		sprintf(name, "residualProfileDYvsY_%i", i);
		sprintf(title, "DUT Residual in Y vs. Track Y - Profile, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];Y_{hit - track} [mm]", sensname.c_str(),l,m,n,i);
		residualProfileDYvsY.push_back(new TProfile (name, title,150,-15,15,-0.2,0.2,"s"));
		_rootObjectMap[name] = residualProfileDYvsY[j] ;

		// residuals vs event
		sprintf(name, "residualsXvsEvt_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Event Nr., %s, F=%.1e, %iV, %irot, run %i;Event Nr.;X_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsXvsEvt.push_back(new TH2D(name, title, 500, 0, 5e5, 100, -0.5, 0.5));
		_rootObjectMap[name] = residualsXvsEvt[j] ;

		sprintf(name, "residualsYvsEvt_%i" ,i);
		sprintf(title, "DUT Residual in Y vs. Event Nr., %s, F=%.1e, %iV, %irot, run %i;Event Nr.;Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsYvsEvt.push_back(new TH2D(name, title, 500, 0, 5e5, 100, -0.5, 0.5));
		_rootObjectMap[name] = residualsYvsEvt[j] ;

		sprintf(name, "residualsZvsEvt_%i" ,i);
		sprintf(title, "DUT Residual in Z vs. Event Nr., %s, F=%.1e, %iV, %irot, run %i;Event Nr.;Z_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualsZvsEvt.push_back(new TH2D(name, title, 500, 0, 5e5, 100, -0.5, 0.5));
		_rootObjectMap[name] = residualsZvsEvt[j] ;

		// residualmaps
		sprintf(name, "residualMapTemp_%i" ,i);
		sprintf(title, "DUT Residual in X vs. Track Position, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		residualMapTemp.push_back(new TH3D(name, title, 50, -10, 10, 50, -10, 10, 100, -0.1, 0.1));
		_rootObjectMap[name] = residualMapTemp[j] ;

		// DUT hitmap
		sprintf(name, "hitmapDUT_%i" ,i);
		sprintf(title, "Measured Hits on the DUT, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [mm];Y_{hit} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapDUT.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapDUT[j] ;

		sprintf(name, "hitmapDUTmodpitch_%i" ,i);
		sprintf(title, "Measured Hits on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [Pitch];Y_{hit} [Pitch];Events", sensname.c_str(),l,m,n,i);
		hitmapDUTmodpitch.push_back(new TH2D(name, title, 100, 0, 1, 100, 0, 1));
		_rootObjectMap[name] = hitmapDUTmodpitch[j] ;

		// track hitmap
		sprintf(name, "hitmapTracks_%i" ,i);
		sprintf(title, "Tracks on the DUT, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapTracks.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapTracks[j] ;

		sprintf(name, "hitmapTracksmodpitch_%i" ,i);
		sprintf(title, "Tracks on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapTracksmodpitch.push_back(new TH2D(name, title, 100, 0, 1, 100, 0, 1));
		_rootObjectMap[name] = hitmapTracksmodpitch[j] ;

		// matched hitmap
		sprintf(name, "hitmapMatch_%i" ,i);
		sprintf(title, "Matched Hits on the DUT, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [mm];Y_{hit} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapMatch.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapMatch[j] ;

		sprintf(name, "hitmapMatchmodpitch_%i" ,i);
		sprintf(title, "Matched Hits on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;X_{hit} [Pitch];Y_{hit} [Pitch];Events", sensname.c_str(),l,m,n,i);
		hitmapMatchmodpitch.push_back(new TH2D(name, title, 100, 0, 1, 100, 0, 1));
		_rootObjectMap[name] = hitmapMatchmodpitch[j];

		sprintf(name, "hitmapMatchTracks_%i" ,i);
		sprintf(title, "Matched Hits on the DUT, Track Positions, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapMatchTracks.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapMatchTracks[j] ;

		// scattering of the beam
		sprintf(name, "scatterX_%i" ,i);
		sprintf(title, "Beam Scattering at DUT in X, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Kink_{X} [mrad]", sensname.c_str(),l,m,n,i);
		scatterX.push_back(new TProfile2D(name, title, 150, -15, 15, 150, -15, 15, -1.5, 1.5, "s"));
		_rootObjectMap[name] = scatterX[j] ;

		sprintf(name, "scatterY_%i" ,i);
		sprintf(title, "Beam Scattering at DUT in Y, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Kink_{Y} [mrad]", sensname.c_str(),l,m,n,i);
		scatterY.push_back(new TProfile2D(name, title, 150, -15, 15, 150, -15, 15, -1.5, 1.5, "s"));
		_rootObjectMap[name] = scatterY[j] ;

		sprintf(name, "scatterXY_%i" ,i);
		sprintf(title, "Beam Scattering at DUT in XY, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];<Kink^{2}> [mrad^{2}]", sensname.c_str(),l,m,n,i);
		scatterXY.push_back(new TProfile2D(name, title, 150, -15, 15, 150, -15, 15, 0, 5, "s"));
		_rootObjectMap[name] = scatterXY[j] ;

		sprintf(name, "kinkX_%i" ,i);
		sprintf(title, "Track Kink at the DUT in X, %s, F=%.1e, %iV, %irot, run %i;Kink_{X} [mrad];Events", sensname.c_str(),l,m,n,i);
		kinkX.push_back(new TH1D(name, title, 100, -2.5, 2.5));
		_rootObjectMap[name] = kinkX[j] ;

		sprintf(name, "kinkY_%i" ,i);
		sprintf(title, "Track Kink at the DUT in Y, %s, F=%.1e, %iV, %irot, run %i;Kink_{Y} [mrad];Events", sensname.c_str(),l,m,n,i);
		kinkY.push_back(new TH1D(name, title, 100, -2.5, 2.5));
		_rootObjectMap[name] = kinkY[j] ;

		// charge sharing to 0, 1, 2, 3, 4 neighbours
		sprintf(name, "hitmapChargeShared0_%i" ,i);
		sprintf(title, "Tracks with Charge Sharing to No Neighbour, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapChargeShared0.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapChargeShared0[j] ;

		sprintf(name, "hitmapChargeShared1_%i" ,i);
		sprintf(title, "Tracks with Charge Sharing to 1 Neighbour, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapChargeShared1.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapChargeShared1[j] ;

		sprintf(name, "hitmapChargeShared2_%i" ,i);
		sprintf(title, "Tracks with Charge Sharing to 2 Neighbours, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapChargeShared2.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapChargeShared2[j] ;

		sprintf(name, "hitmapChargeShared3_%i" ,i);
		sprintf(title, "Tracks with Charge Sharing to 3 Neighbours, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapChargeShared3.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapChargeShared3[j] ;

		sprintf(name, "hitmapChargeShared4_%i" ,i);
		sprintf(title, "Tracks with Charge Sharing to 4 Neighbours, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Events", sensname.c_str(),l,m,n,i);
		hitmapChargeShared4.push_back(new TH2D(name, title, 1000, -15, 15, 1000, -15, 15));
		_rootObjectMap[name] = hitmapChargeShared4[j] ;

		sprintf(name, "hitmapClusterSizeA_%i" ,i);
		sprintf(title, "Any cluster size on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		hitmapClusterSizeA.push_back(new TH1D(name, title, 100, 0, 1));
		_rootObjectMap[name] = hitmapClusterSizeA[j] ;

		sprintf(name, "hitmapClusterSize1_%i" ,i);
		sprintf(title, "CS1 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		hitmapClusterSize1.push_back(new TH1D(name, title, 100, 0, 1));
		_rootObjectMap[name] = hitmapClusterSize1[j] ;

		sprintf(name, "hitmapClusterSize2_%i" ,i);
		sprintf(title, "CS2 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		hitmapClusterSize2.push_back(new TH1D(name, title, 100, 0, 1));
		_rootObjectMap[name] = hitmapClusterSize2[j] ;

		sprintf(name, "hitmapClusterSize3_%i" ,i);
		sprintf(title, "CS3 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		hitmapClusterSize3.push_back(new TH1D(name, title, 100, 0, 1));
		_rootObjectMap[name] = hitmapClusterSize3[j] ;

		sprintf(name, "hitmapClusterSize4_%i" ,i);
		sprintf(title, "CS4 on the DUT mod Pitch, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		hitmapClusterSize4.push_back(new TH1D(name, title, 100, 0, 1));
		_rootObjectMap[name] = hitmapClusterSize4[j] ;

		sprintf(name, "signalClusterSize1_%i" ,i);
		sprintf(title, "CS1 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize1.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = signalClusterSize1[j] ;

		sprintf(name, "signalClusterSize1elec_%i" ,i);
		sprintf(title, "CS1 Signal, Electrons, %s, F=%.1e, %iV, %irot, run %i;Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize1elec.push_back(new TH1D(name, title, 1000, -75000, 75000));
		_rootObjectMap[name] = signalClusterSize1elec[j] ;

		sprintf(name, "signalClusterSize2_%i" ,i);
		sprintf(title, "CS2 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize2.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = signalClusterSize2[j] ;

		sprintf(name, "signalClusterSize2elec_%i" ,i);
		sprintf(title, "CS2 Signal, Electrons, %s, F=%.1e, %iV, %irot, run %i;Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize2elec.push_back(new TH1D(name, title, 1000, -75000, 75000));
		_rootObjectMap[name] = signalClusterSize2elec[j] ;

		sprintf(name, "signalClusterSize3_%i" ,i);
		sprintf(title, "CS3 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize3.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = signalClusterSize3[j] ;

		sprintf(name, "signalClusterSize3elec_%i" ,i);
		sprintf(title, "CS3 Signal, Electrons, %s, F=%.1e, %iV, %irot, run %i;Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize3elec.push_back(new TH1D(name, title, 1000, -75000, 75000));
		_rootObjectMap[name] = signalClusterSize3elec[j] ;

		sprintf(name, "signalClusterSize4_%i" ,i);
		sprintf(title, "CS4 Signal, %s, F=%.1e, %iV, %irot, run %i;Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize4.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = signalClusterSize4[j] ;

		sprintf(name, "signalClusterSize4elec_%i" ,i);
		sprintf(title, "CS4 Signal, Electrons, %s, F=%.1e, %iV, %irot, run %i;Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalClusterSize4elec.push_back(new TH1D(name, title, 1000, -75000, 75000));
		_rootObjectMap[name] = signalClusterSize4elec[j] ;

		// eta distribution if a hit is matched to a track
		sprintf(name, "matchedEta_%i", i);
		sprintf(title, "Eta Distribution of Tracks, if DUT Hit Matched, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		matchedEta.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = matchedEta[j] ;

		sprintf(name, "matchedEta2D_%i" ,i);
		sprintf(title, "Eta Distribution of Tracks, if DUT Hit Matched vs. Track mod Pitch , %s, F=%.1e, %iV, %irot, run %i;Eta [1];Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		matchedEta2D.push_back(new TH2D(name, title, 120, -1, 2, 100, 0, 1));
		_rootObjectMap[name] = matchedEta2D[j];

		sprintf(name, "matchedEtaIntegral_%i", i);
		sprintf(title, "Eta Distribution of Tracks, if DUT Hit Matched, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		matchedEtaIntegral.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = matchedEtaIntegral[j] ;

		// eta distribution if a hit is NOT matched to a track
		sprintf(name, "unmatchedEta_%i", i);
		sprintf(title, "Eta Distribution of Tracks, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		unmatchedEta.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = unmatchedEta[j] ;

		sprintf(name, "unmatchedEta2D_%i" ,i);
		sprintf(title, "Eta Distribution of Tracks vs. Track mod Pitch , %s, F=%.1e, %iV, %irot, run %i;Eta [1];Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		unmatchedEta2D.push_back(new TH2D(name, title, 120, -1, 2, 100, 0, 1));
		_rootObjectMap[name] = unmatchedEta2D[j];

		sprintf(name, "unmatchedEtaIntegral_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		unmatchedEtaIntegral.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = unmatchedEtaIntegral[j] ;

		// eta distribution if a track hits a strip
		sprintf(name, "striphitEta_%i", i);
		sprintf(title, "Eta Distribution of Tracks Hitting a Strip, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		striphitEta.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = striphitEta[j] ;

		sprintf(name, "striphitEta2D_%i" ,i);
		sprintf(title, "Eta Distribution of Tracks Hitting a Strip vs. Track mod Pitch , %s, F=%.1e, %iV, %irot, run %i;Eta [1];Y_{track} [Pitch];Events", sensname.c_str(),l,m,n,i);
		striphitEta2D.push_back(new TH2D(name, title, 120, -1, 2, 100, 0, 1));
		_rootObjectMap[name] = striphitEta2D[j];

		sprintf(name, "striphitEtaIntegral_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		striphitEtaIntegral.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = striphitEtaIntegral[j] ;

		// eta based on noise
		sprintf(name, "noiseEta_%i", i);
		sprintf(title, "Eta Distribution of Noise, Integral, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		noiseEta.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = noiseEta[j] ;

		sprintf(name, "noisesubtractedEta_%i", i);
		sprintf(title, "Eta Distribution of Tracks with noise subtracted, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		noisesubtractedEta.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = noisesubtractedEta[j] ;

		sprintf(name, "EtaL50_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted left 0.5*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaL50.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaL50[j] ;

		sprintf(name, "EtaL40_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted left 0.4*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaL40.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaL40[j] ;

		sprintf(name, "EtaL30_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted left 0.3*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaL30.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaL30[j] ;

		sprintf(name, "EtaL20_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted left 0.2*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaL20.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaL20[j] ;

		sprintf(name, "EtaL10_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted left 0.1*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaL10.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaL10[j] ;

		sprintf(name, "EtaR10_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.1*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaR10.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaR10[j] ;

		sprintf(name, "EtaR20_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.2*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaR20.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaR20[j] ;

		sprintf(name, "EtaR30_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.3*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaR30.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaR30[j] ;

		sprintf(name, "EtaR40_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.4*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaR40.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaR40[j] ;

		sprintf(name, "EtaR50_%i", i);
		sprintf(title, "Eta Distribution of Tracks, Track impact shifted right 0.5*pitch, %s, F=%.1e, %iV, %irot, run %i;Eta [1];Events", sensname.c_str(),l,m,n,i);
		EtaR50.push_back(new TH1D(name, title, 120, -1, 2));
		_rootObjectMap[name] = EtaR50[j] ;

		// the adc count under a track
		sprintf(name, "trackSignal_%i", i);
		sprintf(title, "Signal Under a Track, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		trackSignal.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = trackSignal[j] ;

		// the adc count under a track
		sprintf(name, "trackSignalelec_%i", i);
		sprintf(title, "Signal Under a Track, Electrons, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		trackSignalelec.push_back(new TH1D(name, title, 1000, -500, 500));
		_rootObjectMap[name] = trackSignalelec[j] ;

		// the adc count under a track vs position
		sprintf(name, "trackSignalMap_%i" ,i);
		sprintf(title, "ADCs Under a Track vs. Track Position, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		trackSignalMap.push_back(new TH2D(name, title, 1000, -20, 20, 1000, -500, 500));
		_rootObjectMap[name] = trackSignalMap[j];

		sprintf(name, "trackSignalMapelec_%i" ,i);
		sprintf(title, "ADCs Under a Track vs. Track Position, Electrons, %s, F=%.1e, %iV, %irot, run %i;Y_{track} [mm];(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		trackSignalMapelec.push_back(new TH2D(name, title, 1000, -20, 20, 1000, -75000, 75000));
		_rootObjectMap[name] = trackSignalMapelec[j];

		// the adc count vs tdc of the event
		sprintf(name, "trackSignalTDC_%i" ,i);
		sprintf(title, "ADCs Under a Track vs. Event TDC, %s, F=%.1e, %iV, %irot, run %i;TDC [ns];(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		trackSignalTDC.push_back(new TH2D(name, title, 100, 0, 100, 1000, -500, 500));
		_rootObjectMap[name] = trackSignalTDC[j];

		sprintf(name, "trackSignalTDCelec_%i" ,i);
		sprintf(title, "ADCs Under a Track vs. Event TDC, Electrons, %s, F=%.1e, %iV, %irot, run %i;TDC [ns];(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		trackSignalTDCelec.push_back(new TH2D(name, title, 100, 0, 100, 1000, -75000, 75000));
		_rootObjectMap[name] = trackSignalTDCelec[j];

		// the track count per event
		sprintf(name, "tracksperevent_%i", i);
		sprintf(title, "Tracks per Event, %s, F=%.1e, %iV, %irot, run %i;Tracks;Events", sensname.c_str(),l,m,n,i);
		tracksperevent.push_back(new TH1D(name, title, 20, 0, 20));
		_rootObjectMap[name] = tracksperevent[j] ;

		sprintf(name, "tracksvsevents_%i" ,i);
		sprintf(title, "Tracks per Event vs. Event Nr., %s, F=%.1e, %iV, %irot, run %i;Event Nr.;Tracks per Event;Events", sensname.c_str(),l,m,n,i);
		tracksvsevents.push_back(new TH2D(name, title, 500, 0, 5e5, 15, 0, 14));
		_rootObjectMap[name] = tracksvsevents[j];

		// the distribution of the highest channel in a track
		sprintf(name, "highchanneldistribution_%i", i);
		sprintf(title, "Channel in a Track with Highest Signal, %s, F=%.1e, %iV, %irot, run %i;Channel Distance from Seed;Events", sensname.c_str(),l,m,n,i);
		highchanneldistribution.push_back(new TH1D(name, title, 21, -10, 10));
		_rootObjectMap[name] = highchanneldistribution[j] ;

		// the adc count under a track, interstrip sections
		sprintf(name, "signalLeft2_%i", i);
		sprintf(title, "Signal Under Selected Track, Left2 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalLeft2.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalLeft2[j] ;

		sprintf(name, "signalLeft1_%i", i);
		sprintf(title, "Signal Under Selected Track, Left1 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalLeft1.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalLeft1[j] ;

		sprintf(name, "signalCenter_%i", i);
		sprintf(title, "Signal Under Selected Track, Center Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalCenter.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalCenter[j] ;

		sprintf(name, "signalRight1_%i", i);
		sprintf(title, "Signal Under Selected Track, Right1 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalRight1.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalRight1[j] ;

		sprintf(name, "signalRight2_%i", i);
		sprintf(title, "Signal Under Selected Track, Right2 Strip, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalRight2.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalRight2[j] ;

		sprintf(name, "signalLeft2elec_%i", i);
		sprintf(title, "Signal Under Selected Track, Left2 Strip, Electrons, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalLeft2elec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalLeft2elec[j] ;

		sprintf(name, "signalLeft1elec_%i", i);
		sprintf(title, "Signal Under Selected Track, Left1 Strip, Electrons, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalLeft1elec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalLeft1elec[j] ;

		sprintf(name, "signalCenterelec_%i", i);
		sprintf(title, "Signal Under Selected Track, Center Strip, Electrons, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalCenterelec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalCenterelec[j] ;

		sprintf(name, "signalRight1elec_%i", i);
		sprintf(title, "Signal Under Selected Track, Right1 Strip, Electrons, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalRight1elec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalRight1elec[j] ;

		sprintf(name, "signalRight2elec_%i", i);
		sprintf(title, "Signal Under Selected Track, Right2 Strip, Electrons, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalRight2elec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalRight2elec[j] ;

		sprintf(name, "signalGoodEvents_%i", i);
		sprintf(title, "Signal Under Selected Track, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalGoodEvents.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalGoodEvents[j] ;

		sprintf(name, "signalGoodEventselec_%i", i);
		sprintf(title, "Signal Under Selected Track, Electrons, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalGoodEventselec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalGoodEventselec[j] ;

		sprintf(name, "signalmapA_%i", i);
		sprintf(title, "Signal Under Selected Track, Area A, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalmapA.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalmapA[j] ;

		sprintf(name, "signalmapB_%i", i);
		sprintf(title, "Signal Under Selected Track, Area B, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalmapB.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalmapB[j] ;

		sprintf(name, "signalmapC_%i", i);
		sprintf(title, "Signal Under Selected Track, Area C, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalmapC.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalmapC[j] ;

		sprintf(name, "signalmapD_%i", i);
		sprintf(title, "Signal Under Selected Track, Area D, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [ADCs];Events", sensname.c_str(),l,m,n,i);
		signalmapD.push_back(new TH1D(name, title, 500, -500, 500));
		_rootObjectMap[name] = signalmapD[j] ;

		sprintf(name, "signalmapAelec_%i", i);
		sprintf(title, "Signal Under Selected Track, Electrons, Area A, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalmapAelec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalmapAelec[j] ;

		sprintf(name, "signalmapBelec_%i", i);
		sprintf(title, "Signal Under Selected Track, Electrons, Area B, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalmapBelec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalmapBelec[j] ;

		sprintf(name, "signalmapCelec_%i", i);
		sprintf(title, "Signal Under Selected Track, Electrons, Area C, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalmapCelec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalmapCelec[j] ;

		sprintf(name, "signalmapDelec_%i", i);
		sprintf(title, "Signal Under Selected Track, Electrons, Area D, %s, F=%.1e, %iV, %irot, run %i;(-1) * Signal [Electrons];Events", sensname.c_str(),l,m,n,i);
		signalmapDelec.push_back(new TH1D(name, title, 500, -75000, 75000));
		_rootObjectMap[name] = signalmapDelec[j] ;

		sprintf(name, "signalareaplot_%i",i);
		signalareaplot.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = signalareaplot[j];

		sprintf(name, "signalareaplotelec_%i",i);
		signalareaplotelec.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = signalareaplotelec[j];

		// histos to see what we have cut:
		sprintf(name, "fiducial_discard_%i", i);
		sprintf(title, "Tracks discarded by fiducial cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		fiducial_discard.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = fiducial_discard[j] ;

		sprintf(name, "goodchannel_discard_%i", i);
		sprintf(title, "Tracks discarded by good channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		goodchannel_discard.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = goodchannel_discard[j] ;

		sprintf(name, "goodevent_discard_%i", i);
		sprintf(title, "Tracks discarded by good event cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		goodevent_discard.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = goodevent_discard[j] ;

		sprintf(name, "trackselection_discard_%i", i);
		sprintf(title, "Tracks discarded by track selection cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		trackselection_discard.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = trackselection_discard[j] ;

		sprintf(name, "timecut_discard_%i", i);
		sprintf(title, "Tracks discarded by time cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		timecut_discard.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = timecut_discard[j] ;

		sprintf(name, "highchannel_discard_%i", i);
		sprintf(title, "Tracks discarded by central high channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		highchannel_discard.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = highchannel_discard[j] ;

		sprintf(name, "nohittrack_%i", i);
		sprintf(title, "Tracks with no DUT hit, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		nohittrack.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = nohittrack[j] ;

		// and not cut:
		sprintf(name, "fiducial_allow_%i", i);
		sprintf(title, "Tracks allowed by fiducial cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		fiducial_allow.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = fiducial_allow[j] ;

		sprintf(name, "goodchannel_allow_%i", i);
		sprintf(title, "Tracks allowed by good channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		goodchannel_allow.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = goodchannel_allow[j] ;

		sprintf(name, "goodevent_allow_%i", i);
		sprintf(title, "Tracks allowed by good event cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		goodevent_allow.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = goodevent_allow[j] ;

		sprintf(name, "trackselection_allow_%i", i);
		sprintf(title, "Tracks allowed by track selection cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		trackselection_allow.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = trackselection_allow[j] ;

		sprintf(name, "timecut_allow_%i", i);
		sprintf(title, "Tracks allowed by time cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		timecut_allow.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = timecut_allow[j] ;

		sprintf(name, "highchannel_allow_%i", i);
		sprintf(title, "Tracks allowed by central high channel cut, %s, F=%.1e, %iV, %irot, run %i;1;Events", sensname.c_str(),l,m,n,i);
		highchannel_allow.push_back(new TH1D(name, title, 10, 0, 10));
		_rootObjectMap[name] = highchannel_allow[j] ;

		sprintf(name, "signalMapTemp_%i" ,i);
		sprintf(title, "Signal vs. Track Position, %s, F=%.1e, %iV, %irot, run %i;X_{track} [mm];Y_{track} [mm];Y_{hit - track} [mm];Events", sensname.c_str(),l,m,n,i);
		signalMapTemp.push_back(new TH3D(name, title, 50, -10, 10, 50, -10, 10, 350, -50, 300));
		_rootObjectMap[name] = signalMapTemp[j] ;

	} // done loop over all runs

	// only once -> write in main dir
	_outputFile->cd();

	for (int i = 0;i<_bookinglimit;i++)
	{
		// resolution vs voltage
		sprintf(name, "residualsXvoltage_%i",i);
		residualsXvoltage.push_back(new TGraphErrors(0)); // valgrind complains
		_rootObjectMap[name] = residualsXvoltage[i];

		sprintf(name, "residualsYvoltage_%i",i);
		residualsYvoltage.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = residualsYvoltage[i];

		sprintf(name, "residualsYvoltageAngle_%i",i);
		residualsYvoltageAngle.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = residualsYvoltageAngle[i];

		// signal vs. voltage
		sprintf(name, "signalvoltage_%i",i);
		signalvoltage.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = signalvoltage[i];

		sprintf(name, "signalvoltageelec_%i",i);
		signalvoltageelec.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = signalvoltageelec[i];

		// signal/noise vs. voltage
		sprintf(name, "snvoltage_%i",i);
		snvoltage.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = snvoltage[i];

		// signal/distance vs. voltage
		sprintf(name, "signaldistancevoltage_%i",i);
		signaldistancevoltage.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = signaldistancevoltage[i];

		sprintf(name, "signaldistancevoltageelec_%i",i);
		signaldistancevoltageelec.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = signaldistancevoltageelec[i];

		// noise vs. voltage
		sprintf(name, "noisevoltage_%i",i);
		noisevoltage.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = noisevoltage[i];

		sprintf(name, "noisevoltageelec_%i",i);
		noisevoltageelec.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = noisevoltageelec[i];

		// rgh vs. voltage
		sprintf(name, "rghvoltage_%i",i);
		rghvoltage.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = rghvoltage[i];

		// clustersize vs. voltage
		sprintf(name, "clustersizevoltage_%i",i);
		clustersizevoltage.push_back(new TGraphErrors()); // valgrind complains
		_rootObjectMap[name] = clustersizevoltage[i];

		// clustercount vs. voltage
		sprintf(name, "clustercountvoltage_%i",i);
		clustercountvoltage.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = clustercountvoltage[i];

		// channelcount vs. voltage
		sprintf(name, "channelcountvoltage_%i",i);
		channelcountvoltage.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = channelcountvoltage[i];

		// IV
		sprintf(name, "currentvoltage_%i",i);
		currentvoltage.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = currentvoltage[i];

		sprintf(name, "volumecurrent_%i",i);
		volumecurrent.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = volumecurrent[i];

		sprintf(name, "surfacecurrent_%i",i);
		surfacecurrent.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = surfacecurrent[i];

		// charge sharing
		sprintf(name, "chargesharingvoltage_%i",i);
		chargesharingvoltage.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = chargesharingvoltage[i];

		// eta charge sharing
		sprintf(name, "etachargesharingvoltage_%i",i);
		etachargesharingvoltage.push_back(new TGraphErrors); // valgrind complains
		_rootObjectMap[name] = etachargesharingvoltage[i];

		// signal area map
		sprintf(name, "signalareamap_%i",i);
		signalareamap.push_back(new TH2D(name, "a title", 10,0,1000,4,0,1));
		_rootObjectMap[name] = signalareamap[i];

	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// the landau gaussian convolution
Double_t langaufun(Double_t* x, Double_t* par) 
{

	//Fit parameters:
	//par[0]=Width (scale) parameter of Landau density
	//par[1]=Most Probable (MP, location) parameter of Landau density
	//par[2]=Total area (integral -inf to inf, normalization constant)
	//par[3]=Width (sigma) of convoluted Gaussian function
	//
	//In the Landau distribution (represented by the CERNLIB approximation), 
	//the maximum is located at x=-0.22278298 with the location parameter=0.
	//This shift is corrected within this function, so that the actual
	//maximum is identical to the MP parameter.

	// Numeric constants
	Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	Double_t mpshift  = -0.22278298;       // Landau maximum location

	// Control constants
	Double_t np = 100.0;      // number of convolution steps
	Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

	// Variables
	Double_t xx;
	Double_t mpc;
	Double_t fland;
	Double_t sum = 0.0;
	Double_t xlow,xupp;
	Double_t step;
	Double_t i;

	// MP shift correction
	mpc = par[1] - mpshift * par[0]; 

	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];

	step = (xupp-xlow) / np;

	// Convolution integral of Landau and Gaussian by sum
	for(i=1.0; i<=np/2; i++)
	{
		xx = xlow + (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);

		xx = xupp - (i-.5) * step;
		fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
	}
	return (par[2] * step * sum * invsq2pi / par[3]);
}

// the function performs a gaussian fit around the highest bin content of the histo and uses the gaus sigma and center to set the 
// fit range and starting parameters of the landau gaus fit
// the parameters are: histo to be fitted, number of sigma in the negative direction for the fit range and the same in positive direction
// function to be used
TF1* lanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit) 
{
	int histMaxBin = 0;
	histMaxBin = inHist->GetMaximumBin();
	double histMax = 0.0;
	histMax = inHist->GetXaxis()->GetBinCenter(histMaxBin);

	const double minMaxPosition = 0.5; // minimum position accepted for the max
	if(histMax < minMaxPosition) // if the max is too low it probably it is in the noise
	{
		double tempdouble = 0.0;
		tempdouble = 1 + (minMaxPosition - inHist->GetXaxis()->GetXmin()) / inHist->GetXaxis()->GetBinWidth(5);
		tempdouble = tempdouble + 0.5;
		int startBin = (int) tempdouble;
		int endBin = 0;
		endBin = inHist->GetXaxis()->GetNbins() - 2;
		double yMax = 0;
		for(int iBin = startBin; iBin < endBin; ++iBin)
		{
			if(inHist->GetBinContent(iBin) >= yMax)
			{
				yMax = inHist->GetBinContent(iBin);
				histMax = inHist->GetXaxis()->GetBinCenter(iBin);
			}
		}
	}

	double halfRange = 20; // guess for the range of the gaus fit
	if(inHist->GetXaxis()->GetXmax() > 10000)
	{
		halfRange = 20 * 200; // charge in electrons
	}

	TF1* gausFit = new TF1("gausFit", "gaus", histMax - halfRange, histMax + halfRange);
	gausFit->SetLineColor(kBlue);

	inHist->Fit(gausFit, "RQL");

	double* Gpar = gausFit->GetParameters();

	double fitR1 = 0.0;
	double fitR2 = 0.0;
	fitR1 = Gpar[1] - Gpar[2] * negSigmaFit; // fit range from the gaus mean and sigma
	fitR2 = Gpar[1] + Gpar[2] * posSigmaFit; //...
	double gausSig = 0.0;
	gausSig = Gpar[2]; // initial guess of the gaus sigma

	// Setting fit range and start values
	Double_t fr[2] = {0.0};
	Double_t sv[4] = {0.0};
	Double_t pllo[4] = {0.0};
	Double_t plhi[4] = {0.0};
	fr[0]=fitR1;
	fr[1]=fitR2;

	// find mpv and integral start value
	int binMin = 0; // max and min bin number (corresponding to the range)
	int binMax = 0;
	double intStart = 0.0; // start value of the integral
	double mpvStart = 0.0;
	mpvStart = Gpar[1]; // start value of the mpv

	double binW = 0.0;
	binW = inHist->GetXaxis()->GetBinWidth(5); // bin width from a random bin
	double xMin = 0;
	xMin = inHist->GetXaxis()->GetXmin();
	double tempdouble =0.0;
	if (binW != 0.0)
	{
		tempdouble = 1 + (fitR1 - xMin) / binW;
	}
	tempdouble = tempdouble + 0.5;
	binMin = (int) tempdouble;
	double tempdouble2 = 0.0;
	if (binW != 0.0)
	{
		tempdouble2 = 1 + (fitR2 - xMin) / binW;
	}
	tempdouble2 = tempdouble2 + 0.5;
	binMax = (int) tempdouble2;

	double binCont = 0.0;
	// double yMax = 0; // variable used to look for the maximum (mpv start)
	for(int iBn = binMin; iBn < binMax; ++iBn) // integral in the fit range to get the start value
	{
		binCont = inHist->GetBinContent(iBn);
		intStart += binCont;
		// if(binCont > yMax) 
		// 	 {
		// 	   yMax = binCont;
		// 	   mpvStart = inHist->GetXaxis()->GetBinCenter(iBn);
	    // 	 }
	}

	// starting parameters
	sv[0] = 5;//landau width
	if(inHist->GetXaxis()->GetXmax() > 10000)
	{
		sv[0] = 5 * 200; // charge in electrons
	}
	sv[1] = mpvStart; // mpv landau
	sv[2] = intStart; // integral
	if(gausSig > sv[0])
	{
		sv[3] = sqrt(gausSig * gausSig - sv[0] * sv[0]); // gaussian width
	} else {
		sv[3] = gausSig;
	}

	// parameter limits
	pllo[0]=0.01; pllo[1]=-15.0; pllo[2]=1.0; pllo[3]=gausSig * 0.1;
	plhi[0]=20.0; plhi[1]=200.0; plhi[2]=10000000.0; plhi[3]=gausSig;

	if(inHist->GetXaxis()->GetXmax() > 10000) // charge in electrons
	{
		pllo[0]=0.01; pllo[1]=-15.0; pllo[2]=1.0; pllo[3]=gausSig * 0.1;
		plhi[0]=20.0 * 200; plhi[1]=200.0 * 200; plhi[2]=10000000.0; plhi[3]=gausSig;
	}

	/*
	cout << "Landau width " << sv[0] << " " << pllo[0] << " " << plhi[0] << endl;
	cout << "MPV          " << sv[1] << " " << pllo[1] << " " << plhi[1] << endl;
	cout << "Area         " << sv[2] << " " << pllo[2] << " " << plhi[2] << endl;
	cout << "Gaus sigma   " << sv[3] << " " << pllo[3] << " " << plhi[3] << endl;
	*/

	TF1* ffit = new TF1("lanGausFit", langaufun, fr[0], fr[1], 4);
	//  ffit->SetNpx(1e4);
	ffit->SetParameters(sv);
	ffit->SetParNames("Width","MPV","Area","GSigma");
	ffit->SetLineColor(kRed);

	for (int i = 0; i < 4; i++)
	{
		ffit->SetParLimits(i, pllo[i], plhi[i]);
	}

	//ffit->SetRange(15,150);
	inHist->Fit(ffit,"RQL"); // fit within specified range // valgrind: Conditional jump or move depends on uninitialised value(s) here...

	return ffit;

	delete gausFit;
	delete ffit;
}


Double_t gausLangaufun(Double_t* x, Double_t* par) // a peak at 0 and a landau gaussian convolution
{
	Double_t gausPart = par[0] * TMath::Gaus(*x, par[1], par[5]);
	Double_t langauPart = langaufun(x, &par[2]); // par 0 to 2 belong to the gauss part

	return langauPart + gausPart;
}


// the sigma of the convolution and the one of the noise are constrained to be the same
TF1* gausLanGausFit(TH1* inHist, double negSigmaFit, double posSigmaFit)
{
	TF1* gausFunc = new TF1("gausFunc", "gaus", -30, 5); // referred as g0 in the next comments
	inHist->Fit(gausFunc, "RQL");

	TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

	const int nPars = 6;
	double par[nPars] = {0};

	for(int i = 0; i < 2; ++i) par[i] = gausFunc->GetParameter(i); // get the gaus fit parameters
	for(int i = 2; i < nPars; ++i) par[i] = langauFunc->GetParameter(i - 2); // get parameters form langaus fit

	double parLimHi[nPars] = {0};
	double parLimLo[nPars] = {0};

	parLimLo[0] = -5; // g0 const
	parLimHi[0] = 1000000;
	parLimLo[1] = par[1] - 0.5 * gausFunc->GetParameter(2); // g0 mean
	parLimHi[1] = par[1] + 0.5 * gausFunc->GetParameter(2);

	for(int i = 2; i < nPars; ++i) // allow a 50% variation on the already fitted parameters
	{
		parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
		parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
	}

	const char* parNames[nPars] = {"ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma"};
	std::cout << "Start parameters and limits\n";
	for(int i = 0; i < nPars; ++i)
	    std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
	std::cout << std::endl;

	double fitR1 = inHist->GetXaxis()->GetXmin();
	double fitR2 = inHist->GetXaxis()->GetXmax();

	TF1* gausLang = new TF1("gausLang", gausLangaufun, fitR1, fitR2, nPars);
	gausLang->SetNpx(1000);
	gausLang->SetParameters(par);
	gausLang->SetParNames("ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma");
	for(int i = 0; i < nPars; ++i)
	{
		gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);
	}

	inHist->Fit(gausLang, "RL");

	return gausLang;

	delete gausFunc;
	delete gausLang;
}


// the sigma of the convolution and the one of the noise are constrained to be the same, mean and sigma of the noise are fixed
TF1* gausLanGausFitFixGaus(TH1* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma) // gauss parameters (mean and sigma) from another histo
{
	TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

	const int nPars = 6;
	double par[nPars] = {0};

	// set the gaus fit parameters
	par[0] = inHist->GetBinContent(inHist->FindBin(0)); // constant gets the value of the bin at 0
	par[1] = mean; // these 2 remain fixed
	par[5] = sigma;
	for(int i = 2; i < nPars - 1; ++i)
	{
		par[i] = langauFunc->GetParameter(i - 2); // get parameters form langaus fit, except gaus sigma
	}

	double parLimHi[nPars] = {0};
	double parLimLo[nPars] = {0};

	parLimLo[0] = -5; // g0 const
	parLimHi[0] = 1000000;
	parLimLo[1] = mean; // g0 mean
	parLimHi[1] = mean;
	parLimLo[5] = sigma; // g0 sigma
	parLimHi[5] = sigma;

	for(int i = 2; i < nPars - 1; ++i) // allow a 50% variation on the already fitted parameters, fix gaus sigma
	{
		parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
		parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
	}
	/*
	const char* parNames[nPars] = {"ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma"};
	std::cout << "Start parameters and limits\n";
	for(int i = 0; i < nPars; ++i)
	{
		std::cout << parNames[i] << "\t\t" << par[i] << "    " << parLimLo[i] << "   " << parLimHi[i] << " \n";
	}
	std::cout << std::endl;
	*/
	double fitR1 = inHist->GetXaxis()->GetXmin();
	double fitR2 = inHist->GetXaxis()->GetXmax();

	TF1* gausLang = new TF1("gausLang", gausLangaufun, fitR1, fitR2, nPars);
	gausLang->SetNpx(10000);
	gausLang->SetParameters(par);
	gausLang->SetParNames("ConstG0", "MeanG0", "Width", "MPV", "Area", "GSigma");
	for(int i = 0; i < nPars; ++i)
	{
		gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);
	}

	inHist->Fit(gausLang, "RQL");

	return gausLang;

	delete gausLang;
}


Double_t gausNoiseLangaufun(Double_t* x, Double_t* par) // a peak at 0 and a landau gaussian convolution
{
	Double_t gausPart = par[0] * TMath::Gaus(*x, par[1], par[2]);
	Double_t langauPart = langaufun(x, &par[3]); // par 0 to 2 belong to the gauss part

	return langauPart + gausPart;
}


// just fix the gaus parameters of the gaussian close to 0, the rest is free
TF1* gausLanGausFitFixGausNoise(TH1D* inHist, double negSigmaFit, double posSigmaFit, double mean, double sigma) // gauss parameters (mean and sigma) from another histo
{
	// subtract the noise from the starting histo to determine the landau gauss parameters
	TF1* noise = new TF1("noise", "gaus", inHist->GetXaxis()->GetXmin(), inHist->GetXaxis()->GetXmax());
	noise->SetParameter(0, inHist->GetBinContent(inHist->FindBin(0)));
	noise->SetParameter(1, mean);
	noise->SetParameter(2, sigma);
	TH1D* subHist = new TH1D(*inHist);
	subHist->Add(noise , -1);

	TF1* langauFunc = lanGausFit(subHist, negSigmaFit, posSigmaFit);
	//TF1* langauFunc = lanGausFit(inHist, negSigmaFit, posSigmaFit);

	const int nPars = 7;
	double par[nPars] = {0};

	// set the gaus fit parameters
	// constant gets the value of the bin at 0
	par[0] = inHist->GetBinContent(inHist->FindBin(0));
	// if the bin is empty, parameter should be at least 1 to prevent error from the lower limit, which is 1.
	/*
	if (par[0] < 1.0)
	{
		par[0] = 1.0;
	}
	*/

	par[1] = mean; // these 2 remain fixed
	par[2] = sigma;
	for(int i = 3; i < nPars; ++i)
	{
		par[i] = langauFunc->GetParameter(i - 3); // get parameters form langaus fit, except gaus sigma
	}

	double parLimHi[nPars] = {0};
	double parLimLo[nPars] = {0};

	parLimLo[0] = -5; // g0 const
	parLimHi[0] = 1000000;
	parLimLo[1] = mean; // g0 mean
	parLimHi[1] = mean;
	parLimLo[2] = sigma; // g0 sigma
	parLimHi[2] = sigma;

	for(int i = 3; i < nPars; ++i) // allow a 50% variation on the already fitted parameters, fix gaus sigma
	{
		parLimLo[i] = par[i] - 0.5 * fabs(par[i]);
		parLimHi[i] = par[i] + 0.5 * fabs(par[i]);
	}

	/*
	const char* parNames[nPars] = {"ConstG0", "MeanG0", "SigmaG0", "Width", "MPV", "Area", "GSigma"};
	cout << "Start parameters and limits:" << endl;
	for(int i = 0; i < nPars; ++i)
	{
		cout << " " << parNames[i] << ": Start value: " << par[i] << " , low limit: " << parLimLo[i] << " , high limit: " << parLimHi[i] << endl;
	}
	*/

	double fitR1 = inHist->GetXaxis()->GetXmin();
	double fitR2 = inHist->GetXaxis()->GetXmax();

	TF1* gausLang = new TF1("gausLang", gausNoiseLangaufun, fitR1, fitR2, nPars);
	gausLang->SetNpx(10000);
	gausLang->SetParameters(par);
	gausLang->SetParNames("ConstG0", "MeanG0", "SigmaG0", "Width", "MPV", "Area", "GSigma");
	for(int i = 0; i < nPars; ++i)
	{
		gausLang->SetParLimits(i, parLimLo[i], parLimHi[i]);
	}

	inHist->Fit(gausLang, "RQL");

	return gausLang;

	delete gausLang;
}
