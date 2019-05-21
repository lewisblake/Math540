#include <iostream>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <fstream>
#include <Eigen/Dense>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <numeric>
#include <algorithm>
using namespace std;
using namespace Eigen;


int main(int argc, char * argv[] ){
	/*
 *	--------
 *	PREAMBLE
 *	--------
 * */
	
	// Initalize MPI & get status
	MPI_Init(&argc, &argv);
	MPI_Status status;
	// Start the timer
	double startTime = MPI_Wtime();
	// Set precision and scientific notation
	cout.precision(4);
	cout << scientific;
	// Get cluster info
	int myRank, numCoresP;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numCoresP);
	// Initalize random number generator s.t. each core has a different seed
	srand((unsigned)time(0) + myRank);		
	// Define model parameters
	int J = 2;
	int L = J * numCoresP;

	/*
 *	------
 *	STEP 1
 *	------
 * 	*/

	// Each core declares and fills its random weight vector
	//vector<double> randomWeightVec(J);
	double randomWeightVec[J];
	double upper, lower; // declare variables used in for-loop
	double k = myRank + 1; // define k in terms of rank (offset by 1)
	//double localSum;
	for(int ii = 0; ii < J; ++ii){
		upper = ((double)ii + 1) * k + 1;
		lower = k / ((double)ii + 1);
		randomWeightVec[ii] = lower + static_cast <double> (rand()) / static_cast <double> (RAND_MAX/(upper-lower));
	}

	// Master allocates radomWeightVecs and normalizedVecW
	double * randomWeightVecs; // null pointer by default
	double * normalizedVecW = new double[L]; // all cores must allocate normalizedVecW to use with MPI_Bcast()
	//double * globalSumS;
	if(myRank == 0){
		//globalSumS = 0; // For MPI_Reduce()
		randomWeightVecs = new double[L]; // L = J * numCoresP;
	}
	
	// MPI_Gather()
	// Note: the size of recv buff is the size of a *single* recieve, in this case J.
	MPI_Gather(randomWeightVec, J, MPI_DOUBLE, randomWeightVecs, J, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(myRank == 0){
		double S = 0;	
		for(int i = 0; i < L; ++i){
			S += randomWeightVecs[i];
		}
		normalizedVecW = new double[L];
		// Weight normalizedVecW	
		double sumWi = 0;
		for(int i = 0; i < L; ++i){
			normalizedVecW[i] = randomWeightVecs[i]/S;
			sumWi += normalizedVecW[i];
		}

		// Error handling: Check to make sure w_i's sum to 1, within tolerance.
		if (abs(sumWi - 1) > pow(10, -10)){
			cout << "Error: (Step 1), sumWi not equal to 1. MPI Aborting...." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}	
	
	// MPI_Barrier to ensure all cores execute MPI_Bcat() at same time.
	MPI_Barrier(MPI_COMM_WORLD);
	// Broadcast normalizedVecW to all processing cores in MPI_COMM_WORLD
	MPI_Bcast(normalizedVecW, L, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

	/* 
 *	------
 *	STEP 2
 *	------
 * 	*/


	// Each core allocates two J * L matrices for creating J rows of X and C
	MatrixXd mySubMatX(J,L);
	MatrixXd mySubMatC(J,L);
	
	// Each cores build its subset of the X matrix 
	for(int i = 0; i < J; ++i){
		for(k = 0; k < L; ++k){
			upper = ((double)(i+1)*(myRank+1)) / (k+1);
			lower = -((double)(i+1)*(myRank+1)) * (k+1);  
			mySubMatX(i,k) = lower + static_cast <double> (rand()) / static_cast <double> (RAND_MAX/(upper-lower));
		}
	}

	// Each core computes J sample means
	double mySampleMeans[J];
	for(int i = 0; i < J; ++i){
		mySampleMeans[i] = mySubMatX.row(i).mean();
	}

	// Each core computes their J rows of C
	// First fill rows of C and then weights it
	
	// Block Diagonal
	double temp = 0;
	for(int i = 0; i < J; ++i){
		for(int k = 0; k <=i; ++k){// int k = 0; k <= i; ++k
			for(int l = 0; l < L; ++l){
				temp += (mySubMatX(i,l) - mySampleMeans[i])*(mySubMatX(k,l) - mySampleMeans[k])*normalizedVecW[l];
				//mySubMatC(i, myRank * J + k) += (mySubMatX(i,l) - mySampleMeans[i])*(mySubMatX(k,l) - mySampleMeans[k])*normalizedVecW[l];
			
			}	
			// Assign temp to mySubMatC
			mySubMatC(i, myRank * J + k) = temp;
			if(i != k){
				mySubMatC(k, myRank * J + i) = temp;
			}
			// Clear temp
			temp = 0;	
		}
	}
	
	
	// Off Block Diagonal
	// Deal with MPI_Send() first
	if(myRank != 0){ // Master core (0) does not need to send
		for(int iCore = myRank -1; iCore >= 0; --iCore){
			MPI_Send(mySubMatX.data(), J*L, MPI_DOUBLE, iCore, myRank, MPI_COMM_WORLD);
		}	
	}

	// Then deal with MPI_Recv()	
	// Declare recieve buffer
	MatrixXd recvSubMatX(J,L);

	for(int jCore = myRank + 1; jCore < numCoresP; ++jCore){ // Each core recives from only cores with greater rank	
		// Each core recieves rows of X from every core with greater rank
		MPI_Recv(recvSubMatX.data(), J*L, MPI_DOUBLE, jCore, jCore, MPI_COMM_WORLD, &status);	
		// Create the sub block of mySubMatC
		mySubMatC.block(0, jCore * J, J, J) = mySubMatX * recvSubMatX.transpose();	
	}
	

	// Now each core weights it's subset of the C matrix
	double normSquaredNormalizedVecW;
	for(int i = 0 ; i < L; ++i){
		normSquaredNormalizedVecW += pow(normalizedVecW[i], 2.0);
	}
	double weight = 1 / (1 - normSquaredNormalizedVecW); 
	mySubMatC = weight * mySubMatC;


	/* 
 *
 *	Find Min, Max, and locations
 * */

	// Each processing core computes the min & max of its J * L block of C and identifies corresponding locations in the block
	// Note: locations are given by coordinates assuming indices start at 0
	// Get location of maximum
	MatrixXd::Index maxRow, maxCol;	
	double max = mySubMatC.block(0, myRank * J, J, L - (myRank * J)).maxCoeff(&maxRow, &maxCol);
	// Get location of minimum
	MatrixXd::Index minRow, minCol;
	double min =  mySubMatC.block(0, myRank * J, J, L - (myRank * J)).minCoeff(&minRow, &minCol);
	// Shift column locations of min and max to account for checking submatrix
	minCol = minCol + myRank * J;
	maxCol = maxCol + myRank * J;
	minRow = minRow + myRank * J;
	maxRow = maxRow + myRank * J;


	double * minMaxData;
	double * minMaxRecvBuffer;
	if(myRank == 0 ){
		minMaxData = new double[numCoresP * 7]; // 7 cols corresponding to myRank, min, minRow, minCol, max, maxRow, maxCol
		minMaxRecvBuffer[7];	
	}
	// Each core fills it's send buffer
	double minMaxSendBuffer[7];
	minMaxSendBuffer[0] = myRank; 
	minMaxSendBuffer[1] = min; minMaxSendBuffer[2] = minRow; minMaxSendBuffer[3] = minCol;
	minMaxSendBuffer[4] = max; minMaxSendBuffer[5] = maxRow; minMaxSendBuffer[6] = maxCol;
	
	// Gather min's, max's, and associated data onto the master
	MPI_Gather(minMaxSendBuffer, 7, MPI_DOUBLE, minMaxData, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// Master takes gathered data, sorts it into an Eigen::MatrixXd 
	if(myRank ==0){
		// Transform minMaxData into Eigen::MatrixXd
		MatrixXd gatheredData = Map<Matrix<double, Dynamic, Dynamic, RowMajor>> (minMaxData, numCoresP, 7);
		// Master now finds (C_max, core_max), (C_min, core_min), (i_max, j_max), (i_min, j_min)
		MatrixXd::Index globalMaxRow, globalMinRow;
		double C_max = gatheredData.col(4).maxCoeff(&globalMaxRow);// col w/ index 4 has max's
		double C_min = gatheredData.col(1).minCoeff(&globalMinRow);
		int core_max = gatheredData(globalMaxRow,0);
		int core_min  = gatheredData(globalMinRow,0);
		int i_max = gatheredData(globalMaxRow, 5);
		int j_max = gatheredData(globalMaxRow, 6);
		int i_min = gatheredData(globalMinRow, 2);
		int j_min = gatheredData(globalMinRow, 3);	

		// Outout (C_max, core_max), (C_min, core_min), (i_max, j_max), (i_min, j_min) to console.
		cout << "(C_max, core_max) =  (" << C_max << ", " << core_max << ")" << endl;
		cout << "(C_min, core_min) = (" << C_min << ", " << core_min << ")" << endl;
		cout << "NOTE: (i_max, j_max) and (i_min, j_min) are given assuming array indices start at zero. \n Additionally, only one set of coordinates are given. \n Since C is a symmetric covariance matrix, the same maximums and minimums are obtained in the locations of C given by swapping the row and column coordinates of the extrema presented." << endl;
		cout << "(i_max, j_max) = (" << i_max << ", " << j_max << ")" << endl;
		cout << "(i_min, j_min) = (" << i_min << ", " << j_min << ")" << endl;

	}

	// SPECIAL CASE: In case where J == 2 && numCoresP == 4, have each core print their portion of the C matrix
	if(J == 2 && numCoresP == 4){
		MatrixXd gatheredC(L,L); // Initalized as matrix of zeros
		// Eigen matrices are stored in column major format, so take transpose for sending.
		mySubMatC.transposeInPlace();

		// Gather all the subsets of C onto master and store it in a L * L matrix
		MPI_Gather(mySubMatC.data(), J*L, MPI_DOUBLE, gatheredC.data(), J*L, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if(myRank == 0){
			gatheredC.transposeInPlace();
			// Ensure symmetric since C is not typically ever stored in one location
			for(int i = 0; i < L; ++i){
				for(int j = 0; j < L; ++j){
					if(gatheredC(j,i) != gatheredC(i,j)){
						gatheredC(j,i) = gatheredC(i,j);
					}
				}
			}
			// Print out C matrix for special case
			cout << "The C matrix is: " << endl <<  gatheredC << endl;
		}
	}	

	/* 
 *	--------
 *	CLEAN UP
 * 	--------
 * 	*/

	delete [] normalizedVecW;
	double endTime = MPI_Wtime() - startTime;
	// Delete dynamically created arrays
	if(myRank == 0){ // Master created randomWeightVecs so it must also delete it. On all other cores it's a null pointer. 
		delete [] randomWeightVecs;	
		if(J != 2 || numCoresP != 4){				
			cout << "MPI walltime: " << endTime << endl;
		}

	}

	/*
* 	-----------
*	MPI FINALIZE
*	-----------
*	*/

	// Finalize MPI
	MPI_Finalize();

}
