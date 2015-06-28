//STL
#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

//Eigen
#include <Eigen/Core>
#include <Eigen/Dense>

//GTest
#include "gtest/gtest.h"

//EUTelescope
#include "eutelgeotest.h"

#define PI 3.14159265

namespace eugeo = eutelescope::geo;

// The fixture for testing class eutelgeotest. From google test primer.
class eutelgeotestTest : public ::testing::Test {
protected:
	// You can remove any or all of the following functions if its body is empty.

	eutelgeotestTest() {
		// You can do set-up work for each test here.
		seed = std::chrono::system_clock::now().time_since_epoch().count();
	}

	virtual ~eutelgeotestTest() {
		// You can do clean-up work that doesn't throw exceptions here.
	}

	// If the constructor and destructor are not enough for setting up
	// and cleaning up each test, you can define the following methods:
	virtual void SetUp() {
	
		// Code here will be called immediately after the constructor (right before each test).
		generator.seed( seed );
		std::uniform_real_distribution<double> distribution(0.1,8.0);
	
		xVec[0] = distribution(generator); 
		xVec[1] = 0; 
		xVec[2] = 0;

		yVec[0] = 0; 
		yVec[1] = distribution(generator); 
		yVec[2] = 0;

		zVec[0] = 0; 
		zVec[1] = 0; 
		zVec[2] = distribution(generator);

		rVec[0] = 0; 
		rVec[1] = 0; 
		rVec[2] = 0;

		xVecN[0] = 1; 
		xVecN[1] = 0; 
		xVecN[2] = 0;

		yVecN[0] = 0; 
		yVecN[1] = 1;
		yVecN[2] = 0;

		zVecN[0] = 0; 
		zVecN[1] = 0; 
		zVecN[2] = 1; 
		}

	virtual void TearDown() {
		// Code here will be called immediately after each test (right before the destructor).
	}

	// Objects declared here can be used by all tests in the test case for eutelgeotest.
	//creating this will create the gear mgr
	eutelgeotest p;
	unsigned seed;
	std::default_random_engine generator;
	double xVec [3];
	double yVec [3];
	double zVec [3];
	
	double xVecN [3];
	double yVecN [3];
	double zVecN [3];
	
	double rVec [3];
};

// Test case must be called the class above
// Also note: use TEST_F instead of TEST to access the test fixture (from google test primer)

/** This test will take randomly generated vectors, transform them into the global frame
 *  and back into the local one again. This must return the original vector again. It is
 *  done for every plance 100 times.
 */ 
TEST_F(eutelgeotestTest, RandomBackAndForthVectorTrans) {

	auto sensorIdVec = eugeo::gGeometry().sensorIDsVec();

	std::uniform_real_distribution<double> distribution(0.0,3.0);
	
	double const abs_err = 1e-13;

	double vecInitial [3] = {0, 0, 0};
	double vecTrans [3] = {0, 0, 0};
	double vecFinal [3] = {0, 0, 0};

	for(auto sensorID: sensorIdVec) {
		for(size_t i = 0; i < 100; i++) {
			vecInitial[0] = distribution(generator);
			vecInitial[1] = distribution(generator);
			vecInitial[2] = distribution(generator);
		
			eugeo::gGeometry().local2MasterVec(sensorID, vecInitial, vecTrans);
			eugeo::gGeometry().master2LocalVec(sensorID, vecTrans, vecFinal);
		
			ASSERT_NEAR(vecInitial[0], vecFinal[0], abs_err);
			ASSERT_NEAR(vecInitial[1], vecFinal[1], abs_err);
			ASSERT_NEAR(vecInitial[2], vecFinal[2], abs_err);
		}	
	}
}
/** The same as the RandomBackAndForthVectorTrans test. but in this case we transform points
 *  and not vectors.
 */
TEST_F(eutelgeotestTest, RandomBackAndForthPointTrans) {

	auto sensorIdVec = eugeo::gGeometry().sensorIDsVec();

	std::uniform_real_distribution<double> distribution(0.0,3.0);
	
	double const abs_err = 1e-13;

	double pointInitial [3] = {0, 0, 0};
	double pointTrans [3] = {0, 0, 0};
	double pointFinal [3] = {0, 0, 0};

	for(auto sensorID: sensorIdVec) {
		for(size_t i = 0; i < 100; i++) {
			pointInitial[0] = distribution(generator);
			pointInitial[1] = distribution(generator);
			pointInitial[2] = distribution(generator);
		
			eugeo::gGeometry().local2Master(sensorID, pointInitial, pointTrans);
			eugeo::gGeometry().master2Local(sensorID, pointTrans, pointFinal);
		
			ASSERT_NEAR(pointInitial[0], pointFinal[0], abs_err);
			ASSERT_NEAR(pointInitial[1], pointFinal[1], abs_err);
			ASSERT_NEAR(pointInitial[2], pointFinal[2], abs_err);
		}	
	}
}

TEST_F(eutelgeotestTest, SpecificVectorTransYAxis) {

	double const abs_err = 1e-8;
	//SensorID 1 corresponds to only Y-axis rotation
	int const sensorID = 1;
	eugeo::gGeometry().local2MasterVec(sensorID, xVec, rVec);
	ASSERT_NEAR(rVec[0], xVec[0]*cos(45*PI/180), abs_err);
	ASSERT_NEAR(rVec[1], 0, abs_err);
	ASSERT_NEAR(rVec[2], xVec[0]*sin(-45*PI/180), abs_err);
	
	eugeo::gGeometry().local2MasterVec(sensorID, yVec, rVec);
	ASSERT_NEAR(rVec[0], yVec[0], abs_err);
	ASSERT_NEAR(rVec[1], yVec[1], abs_err);
	ASSERT_NEAR(rVec[2], yVec[2], abs_err);

	eugeo::gGeometry().local2MasterVec(sensorID, zVec, rVec);
	ASSERT_NEAR(rVec[0], zVec[2]*cos(45*PI/180), abs_err);
	ASSERT_NEAR(rVec[1], 0, abs_err);
	ASSERT_NEAR(rVec[2], zVec[2]*sin(45*PI/180), abs_err);	
}

TEST_F(eutelgeotestTest, SpecificVectorTransXAxis) {

	double const abs_err = 1e-8;
	//SensorID 2 corresponds to only X-axis rotation
	int const sensorID = 2;

	eugeo::gGeometry().local2MasterVec(sensorID, xVec, rVec);
	ASSERT_NEAR(rVec[0], xVec[0], abs_err);
	ASSERT_NEAR(rVec[1], xVec[1], abs_err);
	ASSERT_NEAR(rVec[2], xVec[2], abs_err);
	
	eugeo::gGeometry().local2MasterVec(sensorID, yVec, rVec);
	ASSERT_NEAR(rVec[0], 0, abs_err);
	ASSERT_NEAR(rVec[1], yVec[1]*cos(45*PI/180), abs_err);
	ASSERT_NEAR(rVec[2], yVec[1]*sin(45*PI/180), abs_err);

	eugeo::gGeometry().local2MasterVec(sensorID, zVec, rVec);
	ASSERT_NEAR(rVec[0], 0, abs_err);
	ASSERT_NEAR(rVec[1], zVec[2]*sin(-45*PI/180), abs_err);
	ASSERT_NEAR(rVec[2], zVec[2]*cos(45*PI/180), abs_err);	
}

TEST_F(eutelgeotestTest, SpecificVectorTransZAxis) {

	double const abs_err = 1e-8;
	//SensorID 3 corresponds to only Z-axis rotation
	int const sensorID = 3;

	eugeo::gGeometry().local2MasterVec(sensorID, xVec, rVec);
	ASSERT_NEAR(rVec[0], xVec[0]*cos(45*PI/180), abs_err);
	ASSERT_NEAR(rVec[1], xVec[0]*sin(45*PI/180), abs_err);
	ASSERT_NEAR(rVec[2], 0, abs_err);
	
	eugeo::gGeometry().local2MasterVec(sensorID, yVec, rVec);
	ASSERT_NEAR(rVec[0], yVec[1]*sin(-45*PI/180), abs_err);
	ASSERT_NEAR(rVec[1], yVec[1]*cos(45*PI/180), abs_err);
	ASSERT_NEAR(rVec[2], 0, abs_err);

	eugeo::gGeometry().local2MasterVec(sensorID, zVec, rVec);
	ASSERT_NEAR(rVec[0], zVec[0], abs_err);
	ASSERT_NEAR(rVec[1], zVec[1], abs_err);
	ASSERT_NEAR(rVec[2], zVec[2], abs_err);	
}

TEST_F(eutelgeotestTest, NormalVectorTest) {

	double const abs_err = 1e-8;
	
	auto sensorIDVec = eugeo::gGeometry().sensorIDsVec();
	for( auto sensorID: sensorIDVec ) {
		eugeo::gGeometry().local2MasterVec(sensorID, xVecN, rVec);
		Eigen::Vector3d xV( rVec[0], rVec[1], rVec[2] );
		
		eugeo::gGeometry().local2MasterVec(sensorID, yVecN, rVec);
		Eigen::Vector3d yV( rVec[0], rVec[1], rVec[2] );

		eugeo::gGeometry().local2MasterVec(sensorID, zVecN, rVec);
		Eigen::Vector3d zV( rVec[0], rVec[1], rVec[2] );
		
		Eigen::Vector3d zVComp = xV.cross(yV);

		TVector3 normGeoFW = eugeo::gGeometry().siPlaneNormal( sensorID );

		Eigen::Vector3d zVGeoFW ( normGeoFW[0], normGeoFW[1], normGeoFW[2] );
		std::cout << "X: " << eugeo::gGeometry().siPlaneXAxis( sensorID )[0] << " Y: " << eugeo::gGeometry().siPlaneYAxis( sensorID )[1] << std::endl;
		std::cout << sensorID << std::endl;
		for(size_t i = 0; i < 3; i++) {
			ASSERT_NEAR(zVGeoFW[i], zVComp[i], abs_err);	
			ASSERT_NEAR(zVGeoFW[i], zV[i], abs_err);
			ASSERT_NEAR(zVGeoFW.norm(), 1, abs_err);
		}
	}
}	
// }  // namespace - could surround eutelgeotestTest in a namespace
