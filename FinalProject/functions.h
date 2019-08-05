/*
 * -----------
 * functions.h
 * -----------	
 * */
#include <iostream>
//#include <cmath> // Need to include more from main?
//#include <string>

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

using namespace std;
using namespace Eigen;

/*
 *	--------------
 *	DEFINE STRUCTS
 *	--------------
 */


// Define methods to read/write Eigen matrix to binary file
// Source: https://stackoverflow.com/questions/25389480/how-to-write-read-an-eigen-matrix-from-binary-file
//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix_MxN;

// Define struct to recieve getUserInput() output
struct UserSetParameters{
	string dataSource;
	string calculationType;
	int NUM_LEVELS_M;
	int NUM_PARTITIONS_J;
	int NUM_KNOTS_r;
	double offsetPercentage;
	int NUM_LEVELS_SERIAL_S;
};

// Define struct to recieve loadData() out
struct LoadDataOut{
	MatrixXd outputDataMat;
	Vector4d domainBoundaries; 
	int nObs;
};

// Define struct to recieve buildStructure out
struct BuildStructureOut{
	vector<MatrixXd> knots;
	vector<MatrixXd> partitions;
	VectorXi nRegions;
	vector<vector<double>> outputData;	
};

// Define structure to recieve createKnots out
struct CreateKnotsOut{
	MatrixXd knots;
	//MatrixXd Y;
};

// Define struct to recieve indexMatrix
struct CreateIndexMatrixOut{
	VectorXi indexMatrix;
	MatrixXi matrixOfRegionsAtFirstParallelLevel;
};


// generic print vector function
void printVector(vector<int> thisVec){
    for(int i=0; i < thisVec.size(); i++){
	cout << " ";
        cout << thisVec.at(i) <<  " ";
	}
}




/*
 *	----------------------
 *	MAJOR MODEL PROCEDURES
 *	---------------------- 
 */

// getUserInput prototype, its declaration
UserSetParameters getUserInput();

// loadData prototype
LoadDataOut loadData(const std::string& dataSource, double offsetPercentage);

// buildStructure prototype
BuildStructureOut buildStructure(int NUM_LEVELS_M, int NUM_PARTITIONS_J, int NUM_KNOTS_r, int NUM_LEVELS_SERIAL_S, int numCores, int myRank, int nObsData, double offsetPercentage, Vector4d domainBoundaries, MatrixXd *data); 



/*
 *	-----------------
 *	SMALLER FUNCTIONS
 *	-----------------
 */

// createIndexMatrix prototype
CreateIndexMatrixOut createIndexMatrix(int NUM_LEVELS_M, int NUM_PARTITIONS_J, int numCores, int myRank, int NUM_LEVELS_SERIAL_S, int nTotalRegionsAssignedToEachWorker, VectorXi nRegions, VectorXi cumulativeRegions);
// findIndex() prototype
int findIndex(int level, int tileNum, VectorXi nRegions);
// findAncestry() prototype
VectorXi findAncestry(int index, VectorXi nRegions, VectorXi cumulativeRegions, int NUM_PARTITIONS_J);
// findParent() prototype
Array3i findParent(int index, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions);
// findDescendants() prototype
VectorXi findDescendants(int index, int NUM_LEVELS_M, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions);
// findLevelAndTile() prototype
Array2i findLevelAndTile(int index, VectorXi nRegions);
// createKnots() prototype
MatrixXd createKnots(double xMin, double xMax, double yMin, double yMax, double offsetPercentage, int nX, int nY);
// createPartition() prototype
MatrixXd createPartition(double xMin0, double xMax0, double yMin0, double yMax0, int NUM_PARTITIONS_J);
// findChildren() prototype
VectorXi findChildren(int index, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions);
// findChild() prototype
VectorXi findChild(int index, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions);
// validateUserInput prototype
void validateUserInput(string calculationType, int NUM_LEVELS_M, int NUM_PARTITIONS_J, int NUM_KNOTS_r, double offsetPercentage, int numCores, int NUM_LEVELS_SERIAL_S);




#endif
