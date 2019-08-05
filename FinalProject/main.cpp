// Include other files
#include <iostream>
#include <mpi.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
//#include <random>
//#include <cstdlib>
//#include <ctime>
#include <cmath>
//#include <numeric>
//#include <algorithm>
//#include <string>

// Includemy own files
#include "functions.h"
#include "functions.cpp"
#include "buildStructure.cpp"

// Declare namespaces
using namespace Eigen;
using namespace std;

// main function
int main(int argc, char * argv[]){
	/*
 *	-------
 *	PREAMBLE
 *	--------
 * */
	// Initialize MPI & get status
	MPI_Init(&argc, &argv); // LB : Can I do anything cool with argc & argv here?
	MPI_Status status;
	// Start the timer
	double startTime = MPI_Wtime();
	// set precision and scientific notation
	std::cout.precision(4);
	std::cout << scientific;
	// Get cluster info
	int myRank, numCores;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numCores);
	// Define "useful" variables for main
	int ierror = 1; // returning 1 implies something went wrong	

	/*
 *	-----------
 *	GET USER INPUT
 *	-----------
 * */
	UserSetParameters userInput = getUserInput();
	string dataSource = userInput.dataSource;
	string calculationsType = userInput.calculationType;
	int NUM_LEVELS_M = userInput.NUM_LEVELS_M;
	int NUM_PARTITIONS_J = userInput.NUM_PARTITIONS_J;
	int NUM_KNOTS_r = userInput.NUM_KNOTS_r;
	int offsetPercentage = userInput.offsetPercentage;
	int NUM_LEVELS_SERIAL_S = userInput.NUM_LEVELS_SERIAL_S;
	// Progress indicators
	if(myRank == 0){cout << "User input loaded. \n Loading data... \n";}

	/*
 *	-------------------
 *	VALIDATE USER INPUT
 *	-------------------
 * */
	void validateUserInput(string calculationType, int NUM_LEVELS_M, int NUM_PARTITIONS_J, int NUM_KNOTS_r, double offsetPercentage, int numCores, int NUM_LEVELS_SERIAL_S);




	/*
 * 	-----------
 *	LOAD THE DATA
 *	-----------
 * */
	LoadDataOut loadedDataOutput = loadData(dataSource, offsetPercentage);
	MatrixXd data = loadedDataOutput.outputDataMat;
	Vector4d domainBoundaries = loadedDataOutput.domainBoundaries;
	int nObsData = loadedDataOutput.nObs;
	if(myRank == 0){cout << "Data loaded. \n Building structure...\n";}

	/*
 *	-------------------
 *	BUILD THE STRUCTURE
 *	-------------------
 * */
	
	BuildStructureOut builtStructureOutput = buildStructure(NUM_LEVELS_M, NUM_PARTITIONS_J, NUM_KNOTS_r, NUM_LEVELS_SERIAL_S, numCores, myRank, nObsData, offsetPercentage, domainBoundaries, &data);
	if(myRank == 0){cout << "Structure built. \n";}
	
	/*
 *	--------
 *	CLEAN UP & MPI FINALIZE
 *	--------
 * */

	double endTime = MPI_Wtime() - startTime;
	double totalTime = 0;
	//MPI_Reduce(&endTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	// Have master report total execution time
	if(myRank == 0){
		cout << "Total MPI_Wtime: " << endTime << endl;
	//	cout << "T_ave: " << totalTime/(double)numCores << endl;
	}
	// Finalize MPI
	MPI_Finalize();


	
}
