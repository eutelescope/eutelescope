// Version: $Id$
//! Precluster for Tracker Detailed Simulation
/*! 
  Simple rectangular precluster from the charge deposits
  stored in map.
  A few useful methods for charge distribution analysis are present.

  @author Piotr Niezurawski

  Date: 2008-12-11
*/

/* 
   Description: Precluster for Tracker Detailed Simulation.

   Author: Piotr Niezurawski

   Date: 2009-01-14
*/

#include <TDSPrecluster.h>


using namespace TDS;
using namespace std;



//! Vector of pixels' charges sorted in charge in descending order
std::vector<double> TDSPrecluster::getVecCharges_DescendingInCharge()
{
  std::vector<double> theVecCharges_DescendingInCharge;
  theVecCharges_DescendingInCharge.reserve(vectorOfPixels.size());

  for (std::vector<TDSPixel>::iterator i = vectorOfPixels.begin(); i != vectorOfPixels.end(); ++i)
    {
      theVecCharges_DescendingInCharge.push_back( i->getCharge() );
    }
  return theVecCharges_DescendingInCharge;
}


//! Vector of pixels' charges sorted in |charge| in descending order
std::vector<double> TDSPrecluster::getVecCharges_DescendingInAbsCharge()
{
  // Temporary multimap for sorting
  typedef std::multimap<double, std::vector<TDSPixel>::iterator, std::greater<double> /* descending order */> type_MultimapForSorting;
  type_MultimapForSorting tempMultimap;

  // Multimap sorts pixels iterators according to |charge| !!!
  for (std::vector<TDSPixel>::iterator i = vectorOfPixels.begin(); i != vectorOfPixels.end(); ++i)
    {
      tempMultimap.insert( make_pair( std::abs(i->getCharge()), i) );
    }

  // Final vector
  std::vector<double> theVecCharges_DescendingInAbsCharge;
  theVecCharges_DescendingInAbsCharge.reserve(vectorOfPixels.size());

  for (type_MultimapForSorting::iterator j = tempMultimap.begin(); j != tempMultimap.end(); ++j)
    {
      theVecCharges_DescendingInAbsCharge.push_back( j->second->getCharge() );
    }
      
  return theVecCharges_DescendingInAbsCharge;
}



//! Vector of pixels' charges sorted in charge/distance_from_seed in descending order
std::vector<double> TDSPrecluster::getVecCharges_DescendingInChargeByDistance()
{
  // Temporary multimap for sorting
  typedef std::multimap<double, std::vector<TDSPixel>::iterator, std::greater<double> /* descending order */> type_MultimapForSorting;
  type_MultimapForSorting tempMultimap;

  // Multimap sorts pixels iterators according to charge/distance !!!
  double distance_from_seed, ratio;
  unsigned long int deltaL, deltaW;
  for (std::vector<TDSPixel>::iterator i = vectorOfPixels.begin(); i != vectorOfPixels.end(); ++i)
    {
      deltaL = pixelL - i->getIndexAlongL();
      deltaW = pixelW - i->getIndexAlongW();
      distance_from_seed = std::sqrt( deltaL*deltaL + deltaW*deltaW );
      if ( distance_from_seed > 0. )
	{
	  ratio = i->getCharge() / distance_from_seed;  // distance in number of pixels between different pixels is always >= 1.
	}
      else
	{
	  ratio = i->getCharge(); // seed
	}
      tempMultimap.insert( make_pair( ratio, i) );
    }

  // Final vector
  std::vector<double> theVecCharges_DescendingInChargeByDistance;
  theVecCharges_DescendingInChargeByDistance.reserve(vectorOfPixels.size());

  for (type_MultimapForSorting::iterator j = tempMultimap.begin(); j != tempMultimap.end(); ++j)
    {
      theVecCharges_DescendingInChargeByDistance.push_back( j->second->getCharge() );
    }
      
  return theVecCharges_DescendingInChargeByDistance;
}



//! Vector of pixels' charges sorted in |charge|/distance_from_seed in descending order
std::vector<double> TDSPrecluster::getVecCharges_DescendingInAbsChargeByDistance()
{
  // Temporary multimap for sorting
  typedef std::multimap<double, std::vector<TDSPixel>::iterator, std::greater<double> /* descending order */> type_MultimapForSorting;
  type_MultimapForSorting tempMultimap;

  // Multimap sorts pixels iterators according to |charge|/distance !!!
  double distance_from_seed, ratio;
  unsigned long int deltaL, deltaW;
  for (std::vector<TDSPixel>::iterator i = vectorOfPixels.begin(); i != vectorOfPixels.end(); ++i)
    {
      deltaL = pixelL - i->getIndexAlongL();
      deltaW = pixelW - i->getIndexAlongW();
      distance_from_seed = std::sqrt( deltaL*deltaL + deltaW*deltaW );
      if ( distance_from_seed > 0 )
	{
	  ratio = std::abs( i->getCharge() ) / distance_from_seed; // distance in number of pixels between different pixels is always >= 1
	}
      else
	{
	  ratio = std::abs( i->getCharge() ); // seed
	}
      tempMultimap.insert( make_pair( ratio, i) );
    }

  // Final vector
  std::vector<double> theVecCharges_DescendingInAbsChargeByDistance;
  theVecCharges_DescendingInAbsChargeByDistance.reserve(vectorOfPixels.size());

  for (type_MultimapForSorting::iterator j = tempMultimap.begin(); j != tempMultimap.end(); ++j)
    {
      theVecCharges_DescendingInAbsChargeByDistance.push_back( j->second->getCharge() );
    }
      
  return theVecCharges_DescendingInAbsChargeByDistance;
}


//! Print to stdout info about precluster
void TDSPrecluster::print()
{
  std::cout << "pixelL=" << pixelL << " " << "pixelW=" << pixelW << std::endl;
  std::cout << "coordL=" << coordL << " " << "coordW=" << coordW << std::endl;
  std::cout << "coordL_chargeCenter=" << coordL_chargeCenter << " " << "coordW_chargeCenter=" << coordW_chargeCenter << std::endl;
  std::cout << "charge=" << charge << std::endl;
  std::cout << "vectorOfPixels.size()=" << vectorOfPixels.size() << std::endl;
  std::vector<TDSPixel>::iterator i;
  for (i=vectorOfPixels.begin(); i!=vectorOfPixels.end(); ++i)
    {
      std::cout << "   iL, iW, cL, cW, charge = " << 
	i->getIndexAlongL() << "\t" << i->getIndexAlongW() << "\t" <<  
	i->getCoordL() << "\t" <<  i->getCoordW() << "\t" <<  i->getCharge() << std::endl;
    }
}



