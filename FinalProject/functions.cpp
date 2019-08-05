/* -----------
 * functions.cpp
 * ----------
 * */


#include <iostream>
#include <mpi.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <typeinfo>
//#include <random>
//#include <cstdlib>
//#include <ctime>
//#include <cmath>
//#include <numeric>
//#include <algorithm>
//#include <string>

// Include own files
#include "functions.h"

using namespace Eigen;
using namespace std;


/*
 *	----------------------
 *	MAJOR MODEL PROCEDURES
 *	----------------------
 * */

// getUserInput() implementation
UserSetParameters  getUserInput(){
	// Open user input file
	ifstream userInputFile;
	userInputFile.open("userInput.txt");
	if(userInputFile.fail()){
		std::cout << "Error in getUserInput(): Loading userInputFile failed. MPI Aborting... " << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	// Extract variables from userInput.txt
	string temp; // declare temporary string; 
	getline(userInputFile, temp); string dataSource = temp;
	getline(userInputFile, temp); string calculationType = temp;
	getline(userInputFile, temp); double NUM_LEVELS_M = stod(temp); // convert string to double
	getline(userInputFile, temp); double NUM_PARTITIONS_J = stod(temp);
	getline(userInputFile, temp); double NUM_KNOTS_r = stod(temp);
	getline(userInputFile, temp); double offsetPercentage = stod(temp);
	getline(userInputFile, temp); double NUM_LEVELS_SERIAL_S = stod(temp);
	temp.clear(); // Clear the temp variable
	// Close userInput.txt
	userInputFile.close();
	// Assign output to struct
	UserSetParameters collectedUserInput = {dataSource, calculationType, NUM_LEVELS_M, NUM_PARTITIONS_J, NUM_KNOTS_r, offsetPercentage, NUM_LEVELS_SERIAL_S};
	// Return the output struct
	return collectedUserInput;
}

// loadData() function implementation
LoadDataOut loadData(const std::string& dataSource, double offsetPercentage){
	const char *ds = dataSource.c_str();
	//cout << "dataSource: " << dataSource << endl;	
	// Open data as binary file	
	ifstream dataInput(dataSource, ios::in | ios::binary);
	if(dataInput.fail()){
		cout << "Error in loadData(): loading data failed. MPI Aborting..." << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	// get length of file
	dataInput.seekg(0, dataInput.end);
	int length = dataInput.tellg(); // in bytes
	dataInput.seekg(0, dataInput.beg);
	// calulate the number of rows (assumed data comes in (x,y,z) coordinates)
	int nRows = length/(3*sizeof(double));
	int nElements = length/(sizeof(double));
	MatrixXd outputDataMat(nRows,3);
	//cout << "The length of the file is: " << length << endl;
	//cout << "There are " << nRows << " rows of data" << endl;
	double * dataArray = new double[length];
	// read the data as a block - stores in column major order
	dataInput.read((char *)dataArray, length);
	// Check to make sure data was read correctly
	if(dataInput){
		// If file contents are read correctly, assign them to the dataArray
		for(int i = 0; i < nElements; ++i ){
			outputDataMat(i) = (double)dataArray[i]; // outputDataMat is store in column major format, so loop through entire contents
		}	
		// Check to make sure all entires are valid
		int numNanElements = isnan(outputDataMat.col(2).array()).count(); // First count how many nan or inf entires there are
		int numInfElements = isinf(outputDataMat.col(2).array()).count();
		if(numNanElements > 0 || numInfElements > 0){
			cout << "Error (loadData): there are NaN or Inf values in the data set. \n Please remove these values from the data before loading it. ";
			dataInput.close();
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	else{ // Else MPI_Abort
		cout << "Error (loadData): only " << dataInput.gcount() << "  could be read" << endl;
		dataInput.close();
		MPI_Abort(MPI_COMM_WORLD,1);	
	}
	// Delete the buffer. Contains entire file
	delete [] dataArray;
	// Determine the boundaries of the domain spaned by the data
	double xmin0  = outputDataMat.col(0).minCoeff(); double xmax0 = outputDataMat.col(0).maxCoeff();
	double ymin0 = outputDataMat.col(1).minCoeff(); double ymax0 = outputDataMat.col(1).maxCoeff();
	Vector4d domainBoundaries(xmin0, xmax0, ymin0, ymax0);
	// Assign output to a LoadDataOut struct instance
	LoadDataOut loadedDataOutput = {outputDataMat, domainBoundaries, nRows};
	// Return outuput struct
	return loadedDataOutput;
}

/*
 *	-----------------
 *	SMALLER FUNCTIONS
 * 	-----------------
 * */

// createParititon() implementation
MatrixXd createPartition(double xMin0, double xMax0, double yMin0, double yMax0, int NUM_PARTITIONS_J){
	MatrixXd partition(NUM_PARTITIONS_J,4); // (NUM_PARTITIONS_J, 4)?
	switch(NUM_PARTITIONS_J){
		case 2:
			if((xMax0-xMin0) >= (yMax0-yMin0)){
				VectorXd tempx = VectorXd::LinSpaced(NUM_PARTITIONS_J + 1, xMin0, xMax0);
				partition.col(0) = tempx.head(NUM_PARTITIONS_J).replicate(NUM_PARTITIONS_J/2,1);
				partition.col(1) = tempx.tail(NUM_PARTITIONS_J).replicate(NUM_PARTITIONS_J/2,1);
				VectorXd tempy = VectorXd::LinSpaced(NUM_PARTITIONS_J/2 + 1, yMin0, yMax0);
				partition.col(2) = tempy.head(NUM_PARTITIONS_J-1).replicate(NUM_PARTITIONS_J/2+1, 1);
				partition.col(3) = tempy.tail(NUM_PARTITIONS_J-1).replicate(NUM_PARTITIONS_J/2+1, 1);					
				
			}else{
				VectorXd tempx = VectorXd::LinSpaced(NUM_PARTITIONS_J/2 + 1, xMin0, xMax0);
				partition.col(0) = tempx.head(NUM_PARTITIONS_J - 1).replicate(NUM_PARTITIONS_J/2 + 1, 1);
				partition.col(1) = tempx.tail(NUM_PARTITIONS_J - 1).replicate(NUM_PARTITIONS_J/2 + 1, 1);
				VectorXd tempy = VectorXd::LinSpaced(NUM_PARTITIONS_J + 1, yMin0, yMax0);
				partition.col(2) = tempy.head(NUM_PARTITIONS_J).replicate(NUM_PARTITIONS_J/2, 1);
				partition.col(3) = tempy.tail(NUM_PARTITIONS_J).replicate(NUM_PARTITIONS_J/2, 1); 			
			}
			break;
		case 4:{
			VectorXd tempx = VectorXd::LinSpaced(NUM_PARTITIONS_J/2 + 1, xMin0, xMax0);
			partition.col(0) = tempx.head(NUM_PARTITIONS_J/2).replicate(NUM_PARTITIONS_J/2, 1);
			partition.col(1) = tempx.tail(NUM_PARTITIONS_J/2).replicate(NUM_PARTITIONS_J/2, 1);
			VectorXd tempy = VectorXd::LinSpaced(NUM_PARTITIONS_J/2 + 1, yMin0, yMax0);	
			partition.col(2) = tempy.head(NUM_PARTITIONS_J/2).replicate(NUM_PARTITIONS_J/2, 1);
			partition.col(3) = tempy.tail(NUM_PARTITIONS_J/2).replicate(NUM_PARTITIONS_J/2, 1);
			break;
			}
		default:
			cout << "Error: (createPartitions): J must be either 2 or 4." << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
	}
	// Return the partition
	return partition;
}

// createKnots() implementation
MatrixXd createKnots(double xMin, double xMax, double yMin, double yMax, double offsetPercentage, int nX, int nY){
	// Calculate and declare intial quantities and Eigen objects
	double offsetX = (xMax-xMin)*offsetPercentage/100.; // can mult by exp(1). keep for testing
	double offsetY = (yMax-yMin)*offsetPercentage/100.;
	VectorXd xlinspace =  VectorXd::LinSpaced(nX, xMin + offsetX, xMax - offsetX);
	VectorXd ylinspace = VectorXd::LinSpaced(nY, yMin + offsetY, yMax - offsetY);	
	MatrixXd knots(nX*nY, 2); //MatrixXd Y(nX, nY);
	// Fill matrices
	for(int i = 0; i < nX; ++i){
		for(int j = 0; j < nY; ++j){
			knots.row(i*nX + j) << xlinspace(i), ylinspace(j); // LB: investigate peformance versus directly indexing in. present seems notationally elegeant.
			//knots(i*nX + j, 1) = ylinspace(j);   
		}
	}
	// Return knots
	return knots;
}

// createIndexMatrix() implementation
CreateIndexMatrixOut createIndexMatrix(int NUM_LEVELS_M, int NUM_PARTITIONS_J, int numCores, int myRank, int NUM_LEVELS_SERIAL_S, int nTotalRegionsAssignedToEachWorker, VectorXi nRegions, VectorXi cumulativeRegions){
	// Calculate inital quantities
	int NUM_LEVEL_BEGIN_PARALLEL = NUM_LEVELS_SERIAL_S + 1;	
	//VectorXi cumulativeRegions(NUM_LEVELS_M); // Calculate cummulativeRegions vector with for-loop
	//for(int i = 0; i < nRegions.array().size(); ++i){cumulativeRegions(i) = nRegions.head(i+1).sum();}
	VectorXi vectorOfRegionsAtFirstParallelLevel = VectorXi::LinSpaced(cumulativeRegions(NUM_LEVELS_SERIAL_S) - nRegions(NUM_LEVELS_SERIAL_S) + 1, nRegions(NUM_LEVELS_SERIAL_S) , cumulativeRegions(NUM_LEVELS_SERIAL_S));
	//cout << " vectorOfRegionsAtFirstParallelLevel: " << endl <<  vectorOfRegionsAtFirstParallelLevel << endl;
	int matRows = nRegions(NUM_LEVEL_BEGIN_PARALLEL - 1)/numCores; // subtract 1 from index b/c arrays start @ 0
	// LB: double check to see if worth finding most efficient way of doing this
	Map<MatrixXi> matrixOfRegionsAtFirstParallelLevel(vectorOfRegionsAtFirstParallelLevel.data(), matRows, numCores);
	//if(myRank == 0 ) {cout << "matrixOfRegionsAtFirstParallelLevel: " << endl << matrixOfRegionsAtFirstParallelLevel << endl;}
	VectorXi indexMatrix(nTotalRegionsAssignedToEachWorker);
	
	if(numCores == nRegions(NUM_LEVELS_SERIAL_S)){
		indexMatrix(NUM_LEVEL_BEGIN_PARALLEL) = vectorOfRegionsAtFirstParallelLevel(myRank);
		indexMatrix.head(NUM_LEVELS_SERIAL_S) = findAncestry(vectorOfRegionsAtFirstParallelLevel(myRank), nRegions, cumulativeRegions, NUM_PARTITIONS_J);
		indexMatrix.tail(nTotalRegionsAssignedToEachWorker - NUM_LEVELS_SERIAL_S) = findDescendants(vectorOfRegionsAtFirstParallelLevel(myRank) , NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
		
	}else if(numCores < nRegions(NUM_LEVELS_SERIAL_S)){
		//cout << "numCores < nRegion(NUM_LEVELS_SERIAL_S) is true" << endl;
		//  Assign the regions allocated to this worker in a vector
		VectorXi thisWorkerParallelAssignment = matrixOfRegionsAtFirstParallelLevel.col(myRank);	
		//cout << "thisWorkerParallelAssignment: " << endl << thisWorkerParallelAssignment << endl;
		// Find all parents of the first parallel level
		int thisRegion, thisRegionParent;
		// For every J regions, there is one parent
		vector<int> tempVecOfParentsOfFirstParallelLevel(thisWorkerParallelAssignment.array().size()); // use std::vector to use unique(). automatically dynamically allocated
		for(int i = 0; i < thisWorkerParallelAssignment.array().size(); ++i){
			thisRegion = thisWorkerParallelAssignment(i);
			Array3i foundParent = findParent(thisRegion, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
			thisRegionParent = foundParent(2);
			tempVecOfParentsOfFirstParallelLevel[i] = thisRegionParent;
		}
		// Select the unique entries of the vector
		tempVecOfParentsOfFirstParallelLevel.erase(unique(tempVecOfParentsOfFirstParallelLevel.begin(), tempVecOfParentsOfFirstParallelLevel.end()), tempVecOfParentsOfFirstParallelLevel.end());
		// Loop through the parent regions of the first level at which to compute in parallel.
		// Find all decendants of the parallel regions parents.
		// These are unique to each column. Moreover, find the ancestors of these regions and store them both in vectors with repeats.	
		int nSingleParentDescendants = 0; // Calculate Parents Descendants
		for(int m = NUM_LEVELS_SERIAL_S; m <= NUM_LEVELS_M; ++m){nSingleParentDescendants += pow(NUM_PARTITIONS_J, (m - NUM_LEVELS_SERIAL_S));}
		int nTilesAssignedToWorker = nRegions(NUM_LEVELS_SERIAL_S)/numCores;
		int nParentsOfTilesAssignedToWorker = tempVecOfParentsOfFirstParallelLevel.size();	
		// The total number of parents descendants is based off how many descendants a single parent has times how many parents of tiles are assigned to a worker.
		// We must also account for how many regions there are at this level 
		// LB: DIRECTLY BELOW NEEDS THOROUGH TESTING
		vector<int> tempVecOfParentsDescendants; //(nSingleParentDescendants*nParentsOfTilesAssignedToWorker - nRegions(NUM_LEVELS_SERIAL_S - 1));
		vector<int> tempVecOfParentsAncestors; //((NUM_LEVELS_SERIAL_S - 1)*nParentsOfTilesAssignedToWorker);
		for(int j = 0; j < nParentsOfTilesAssignedToWorker; ++j){
			// Parent descendants
			VectorXi theseDescendants = findDescendants(tempVecOfParentsOfFirstParallelLevel[j], NUM_LEVELS_M, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
			int nTheseDescendants = theseDescendants.tail(theseDescendants.array().size() - 1).array().size();
			vector<int> theseDescendants2(theseDescendants.tail(theseDescendants.size() - 1).data(), theseDescendants.tail(theseDescendants.size() - 1).data() + nTheseDescendants);
			tempVecOfParentsDescendants.insert(tempVecOfParentsDescendants.end(), theseDescendants2.begin(), theseDescendants2.end());	
			// Parent ancestors
			VectorXi theseAncestors = findAncestry(tempVecOfParentsOfFirstParallelLevel[j], nRegions, cumulativeRegions, NUM_PARTITIONS_J);
			vector<int> theseAncestors2(theseAncestors.data(), theseAncestors.data() + theseAncestors.size());
			tempVecOfParentsAncestors.insert(tempVecOfParentsAncestors.end(), theseAncestors2.begin(), theseAncestors2.end());
		}
		// Sort vectors and select unique entries - potential area for optimization for differing n 
		// (https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector)
		sort(tempVecOfParentsDescendants.begin(), tempVecOfParentsDescendants.end());
		tempVecOfParentsDescendants.erase(unique(tempVecOfParentsDescendants.begin(), tempVecOfParentsDescendants.end()), tempVecOfParentsDescendants.end());
		sort(tempVecOfParentsAncestors.begin(), tempVecOfParentsAncestors.end());
		tempVecOfParentsAncestors.erase(unique(tempVecOfParentsAncestors.begin(), tempVecOfParentsAncestors.end()), tempVecOfParentsAncestors.end());
		// Need to map std::vector<int> into Eigen::VectorXi to get everything into indexMatrix
		Map<VectorXi> parentsAncestors(tempVecOfParentsAncestors.data(), tempVecOfParentsAncestors.size());
		Map<VectorXi> parentsOfFirstParallelLevel(tempVecOfParentsOfFirstParallelLevel.data(), tempVecOfParentsOfFirstParallelLevel.size());
		Map<VectorXi> parentsDescendants(tempVecOfParentsDescendants.data(), tempVecOfParentsDescendants.size());
		// Put everything into the indexMatrix to return
		indexMatrix << parentsAncestors, parentsOfFirstParallelLevel, parentsDescendants; 
	}else{
		cout << "Error (createIndexMatrix): numCores > nRegions(NUM_LEVELS_SERIAL_S - 1)" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	// Return output
	CreateIndexMatrixOut createdIndexMatrixOutput;
	createdIndexMatrixOutput = {indexMatrix, matrixOfRegionsAtFirstParallelLevel};
	return createdIndexMatrixOutput;
}

// findIndex() implementation - tested
int findIndex(int level, int tileNum, VectorXi nRegions){return nRegions.head(level-1).sum() + tileNum;}

// findAncestry() implementation - tested
VectorXi findAncestry(int index, VectorXi nRegions, VectorXi cumulativeRegions, int NUM_PARTITIONS_J){
	VectorXi indexAncestry;
	if(index == 1){return indexAncestry; /*empty object*/ }
	else{
		int indexSmaller = (cumulativeRegions.array() < index).count();
		indexAncestry.resize(indexSmaller);
		// Fill ancestry by looping through parent's parents
		Array3i foundParent; int indexParent;
		for(int k = 0; k < indexSmaller; ++k){
			foundParent = findParent(index, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
			indexParent = foundParent(2);
			indexAncestry((indexSmaller - 1) - k) = indexParent; // Fill from last entry until the zeroth to order from coarsest to finest
			index = indexParent; // update index
		}
	}
	return indexAncestry;
}

// findParent() implementation - tested
Array3i findParent(int index, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions){
	Array3i out(3); // out: <levelParent = indexSmaller, tileParent, indexParent>
	if(index == 1){out << -1, -1, -1; /*set to nonsenical values*/ }
	else{
		int indexSmaller = (cumulativeRegions.array() < index).count();
		int level = indexSmaller + 1; // present level
		int tile = index - cumulativeRegions(indexSmaller - 1);
		int tileParent = ceil((double)tile/(double)NUM_PARTITIONS_J);
		int indexParent = nRegions.head(level - 2).sum() + tileParent;
		out << indexSmaller, tileParent, indexParent;				
	}
	return out;
}

// findDescendants() implementation
VectorXi findDescendants(int index, int NUM_LEVELS_M, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions){
	// Find the level of this index
	Array2i foundLevelAndTile = findLevelAndTile(index, cumulativeRegions);
	int thisLevel = foundLevelAndTile(0);
	// Calculate how many descendants this index has
	int nDescendants = 0;
	for(int i = 0; i <= (NUM_LEVELS_M - thisLevel); ++i){nDescendants += pow(NUM_PARTITIONS_J, i);}
	// Allocate vector for this index's descendants
	VectorXi indexDescendants(nDescendants);
	indexDescendants(0) = index; // All regions at finer resolutions have this index as a common ancestor
	// Loop to calculate descendants
	int thisRegion = 0, counter = 0;
	for(int iThisLevel = 1; iThisLevel <= (NUM_LEVELS_M - thisLevel + 1); ++iThisLevel){
		for(int jRegion = 0; jRegion < pow(NUM_PARTITIONS_J, (iThisLevel - 1)); ++jRegion){
			thisRegion = index*pow(NUM_PARTITIONS_J, (iThisLevel - 1)) + jRegion;
			indexDescendants(counter) = thisRegion;	
			counter++;
		}
	}
	return indexDescendants;
}

// findLevelAndTile() implementation - tested
Array2i findLevelAndTile(int index, VectorXi cumulativeRegions){
	Array2i out; //<level, tile>
	if(index == 1){out << 1, 1;}
	else{
		int indexSmaller = (cumulativeRegions.array() < index).count();
		out << (indexSmaller + 1), (index - cumulativeRegions(indexSmaller - 1));	
	}
	return out;
}

// findChildren() implementation - tested
VectorXi findChildren(int index, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions){
	if(cumulativeRegions.array().size() == 1){
		cout << "Error (findChildren): Only one region. No children if only one region. " << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	VectorXi indexChildren;
	if(index > cumulativeRegions(cumulativeRegions.array().size() - 2)){
	}
	else{
		if(index == 1){
			// Special case for the coarsest resolution
			indexChildren = VectorXi::LinSpaced(NUM_PARTITIONS_J, 2, cumulativeRegions(1));
		}
		else{
			indexChildren = findChild(index, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
		}
	}
	return indexChildren;
}

// findChild() - tested
VectorXi findChild(int index, int NUM_PARTITIONS_J, VectorXi nRegions, VectorXi cumulativeRegions){
	Array2i foundLevelTile = findLevelAndTile(index, cumulativeRegions);
	int level = foundLevelTile(0); int tile = foundLevelTile(1);
	int levelChildren = level + 1; int nChildren = (NUM_PARTITIONS_J * tile - NUM_PARTITIONS_J * (tile-1));
	VectorXi tileChildren = VectorXi::LinSpaced(nChildren, NUM_PARTITIONS_J * (tile-1) + 1, NUM_PARTITIONS_J * tile);
	VectorXi indexChildren(nChildren);
	for(int k = 0; k < nChildren; ++k){
		indexChildren(k) = findIndex(levelChildren, tileChildren(k), nRegions);
	}
	return indexChildren;
}

// validateUserInput implementation
void validateUserInput(string calculationType, int NUM_LEVELS_M, int NUM_PARTITIONS_J, int NUM_KNOTS_r, double offsetPercentage, int numCores, int NUM_LEVELS_SERIAL_S){

	if(!strcmp(calculationType.c_str(), "build_structure")){
		cout << "Error (validateUserInput): build_structure only computational mode available" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if(!(floor(NUM_LEVELS_M) == NUM_LEVELS_M)){
		cout << "Error (validateUserInput): NUM_LEVELS_M must be a positive integer" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if(NUM_PARTITIONS_J != 2 & NUM_PARTITIONS_J != 4){
		cout << "Error (validateUserInput): NUM_PARTITIONS_J must be either 2 or 4" << endl;	
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if(!(floor(NUM_KNOTS_r) == NUM_KNOTS_r)){
		cout << "Error (validateUserInput): NUM_KNOTS_r must be a positive integer" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);	
	}
	if((offsetPercentage < 0 || offsetPercentage >= 1)){
		cout << "Error (validateUserInput): offsetPercentage must be between 0 and 1" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if(abs(round(log(numCores)/log(NUM_PARTITIONS_J)) - log(numCores)/log(NUM_PARTITIONS_J))){
		cout << "Error (validateUserInput:) numCores (given in .slurm file) must be a power of J" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}	
	if(!(floor(NUM_LEVELS_SERIAL_S) == NUM_LEVELS_SERIAL_S) || NUM_LEVELS_SERIAL_S < 1 || NUM_LEVELS_SERIAL_S >= NUM_LEVELS_M){
		cout << "Error (validateUserInput): NUM_LEVELS_SERIAL_S must be a positive integer <= NUM_LEVELS_M" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}


}



