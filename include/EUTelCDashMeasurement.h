/* $Id$ */
/*
  Class for outputing variables in the format expected by CDash, the
  server hosting CTest results. See http://public.kitware.com/Dart/HTML/Tests.shtml

  The variable can either be of type integer, double or string.
  If the string ends in ".png" it will be considered to be the path to a png image
  which would be uploaded by CTest. Other formats are known to CDash but not yet implemented here.

  All measurements will be accessible on the CDash website for the respective test.

  The output is only generated if the DO_TESTING precompiler flag is set.

  Example output:
  <DartMeasurementFile name="TestImage" type="image/png">/home/test/Testing/Temporary/TestTexturedSphere.png</DartMeasurementFile>


  Example usage:

  CDashMeasurement meas1("integer_test",3);
  cout << meas1;

  CDashMeasurement meas2("double_test",3.14);
  cout << meas2;

  CDashMeasurement meas3("string_test","all went good");
  cout << meas3;

  CDashMeasurement meas4("pngfile_test","all_went_good.png");
  cout << meas4;


  Contact: Hanno Perrey <hanno.perrey@desy.de>

*/

#ifndef CDASHM_H_SEEN
#define CDASHM_H_SEEN

#include <iostream>
#include <sstream>

using namespace std;

class CDashMeasurement
{
  stringstream mymeasurement;
  stringstream mtype;
  stringstream mname;
  bool isImage;
  
  void init(){
    mtype.clear();
    mname.clear();
    mymeasurement.clear();
    isImage = false;    
  }

  void setInteger( int m ){
    init();
    mymeasurement << m;
    mtype << "numeric/integer";
  }

  void setDouble( double m ){
    init();
    mymeasurement << m;
    mtype << "numeric/double";
  }

  void setString( string m ){
    init();
    mymeasurement << m;
    string ending=".png";
    // check if string is path to png image
    if (m.length() > ending.length())
      if (m.compare (m.length() - ending.length(), ending.length(), ending))
	isImage = true;
    if (isImage)
      mtype << "image/png";
    else
      mtype << "text/string";
  }

public:
  CDashMeasurement( string name, int value )   { mname << name; setInteger(value); }
  CDashMeasurement( string name, double value ){ mname << name; setDouble(value); }
  CDashMeasurement( string name, string value ){ mname << name; setString(value); }

  // only show output when precompiler flag is set
#ifdef DO_TESTING
  friend ostream& operator<<( ostream& os, const CDashMeasurement& cdm )
  {

    // example output:
    // <DartMeasurementFile name="TestImage" type="image/png">/home/test/Testing/Temporary/TestTexturedSphere.png</DartMeasurementFile>
    os << "<DartMeasurement";
    if (cdm.isImage) os << "File";
    os << " name=\"" << cdm.mname.str() << "\" type=\"" << cdm.mtype.str() << "\">" << cdm.mymeasurement.str()<< "</DartMeasurement";
    if (cdm.isImage) os << "File";
    os << ">" << endl;

    return os;
  }
#else
  // no output case (if testing precompiler flag is not set:)
  friend ostream& operator<<( ostream& os, const CDashMeasurement& )
  {
    return os;
  }

#endif


};


#endif  // CDASHM_H_SEEN
