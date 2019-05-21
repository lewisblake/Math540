
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <mpi.h>
using namespace std;

// Cosine integrand implementation
double cosineIntegrand(double waveNumber, double x ){
    return cos(100 * x - waveNumber * sin(x));
}

// Function to create linear space like MATLAB
vector<double > linspace(double start, double end, int nPoints){
    // Initialize vector for linear space
    vector<double > linearSpace;
    // Error handle special cases
    if(nPoints == 0) {return linearSpace;}
    if(nPoints == 1){
        linearSpace.push_back(start);
	return linearSpace;
    }
    // If passed special cases, make delta
    double delta = (end - start)/((double)nPoints - 1);
    // Loop through each point less one
    for(int i=0; i < nPoints - 1; i++){ // linearSpace.size()
        linearSpace.push_back(start + delta*(double)i);
    }
    // Force end to be last entry in linspace
    linearSpace.push_back(end);
    // Return the linearSpace vector
    return linearSpace;
}

// Problem 1
vector<double > evaluateA(double (fun)(double, double) , double waveNum, vector<double > points, int mPoints){
   	// Intialize output vector
    vector<double > out(3, 0.0); // 3 outputs are M, T, S	
	// Intialize M, T, S (default initalization is 0, but include for clarity)
    double M = 0; double T = 0; double S = 0;
	double len, mpts;
	
	for (int j = 1; j < points.size(); j++){ 
        	// Define len
       		len = points[j] - points[j-1];
        	mpts = (points[j] + points[j-1])/2;
        	// Calculate M and T
        	M += fun(waveNum, mpts) * len;
        	T += (0.5) * (fun(waveNum, points[j-1]) + fun(waveNum, points[j])) * len;
   	 }	
    	// Calculate S with M and T
    	S = (2.0/3.0) * M + (1.0/3.0) * T;
    	// Assign M, T, and S to output vector
    	out[0] = M; out[1] = T; out[2] = S;
   	M = 0; T = 0; S = 0; // For safety reassign to 0 
	// Return out vector
	return out;
}

// evaluateFunction needed for Problem 2 of assignment 2
vector<double > evaluateFunction(vector<double > points, int ptsize, double waveNum){
   vector<double > fun_val_points(ptsize);
   for(int i=0; i < points.size(); i++){
       fun_val_points[i] = cos(100 * points[i] - waveNum * sin(points[i]) );
   }
   return fun_val_points;
}

// getSubIntervalsFunction
vector<double > getSubIntervals(vector<double > kIntervals, int myRank, int numCores, int K){
	vector<double > theseSubIntervals;

	if(myRank == numCores - 1){
		for (int i= myRank *(K/numCores)-1; i<kIntervals.size(); i++ ){
			theseSubIntervals.push_back(kIntervals.at(i));
		}
		
        }
	else if(myRank == 0){
		for (int i = 0; i < K/numCores; i++){
			theseSubIntervals.push_back(kIntervals.at(i));
		}
	}	
	else{	
		for (int i= myRank *(K/numCores)-1; i < (myRank + 1)*(K/numCores); i++ ){
			theseSubIntervals.push_back(kIntervals.at(i));
		}	
    	}
	return theseSubIntervals;
}


// MAIN PROGRAM
int main(int argc, char ** argv){
	// Set precision and scientific notation
	cout.precision(16);
	std::cout << std::scientific;	
	// Declare output file for CSV
	std::ofstream outputCSV;
        outputCSV.open("savedData.csv", ios::out); //output savedData to csv file
        double savedData[3][9901];	

	// Initialize C MPI
	MPI_Init(&argc, &argv);
	MPI_Status status;
	// Get cluster info
	int myRank, numCores;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numCores);

	// Loop through K in [100, 10000] in serial.
	for(int K = 100; K <= 10000; K++){
		// Parition [0, pi] into K subintervals
		vector<double > kIntervals = linspace(0, M_PI, K);
		// Intialize vector of subintervals for this worker
		// Give each core it's appropriate subintervals
		vector<double > theseSubIntervals = getSubIntervals(kIntervals, myRank, numCores, K);	

		double thisWorkersM = 0, thisWorkersT = 0, thisWorkersS = 0;
		vector<double > quadraturePts(101);
		vector<double > thisSubIntEvaluatedA(3);
		double thisSubIntStart = 0, thisSubIntEnd = 0;
		// Each core loops through its subintervals, placing 101 quadratutre points
		for(int iSubInt = 1; iSubInt < theseSubIntervals.size(); iSubInt++){
			thisSubIntStart = theseSubIntervals[iSubInt-1];
			thisSubIntEnd = theseSubIntervals[iSubInt];
	
			quadraturePts = linspace(thisSubIntStart, thisSubIntEnd, 101);
			
			thisSubIntEvaluatedA = evaluateA(cosineIntegrand, K, quadraturePts, 101);
			
			thisWorkersM += thisSubIntEvaluatedA[0];
			thisWorkersT += thisSubIntEvaluatedA[1];
			thisWorkersS += thisSubIntEvaluatedA[2];

		}
			
		// Send / Recieve local M, T, S
		if(myRank == numCores -1){ // Last core (master) recieves data from all other cores (workers)
			// Declare array to recieve the data from workers
			double recvWorkerDataArray[3] = {0, 0, 0};

			// Loop through each core collecting data
			for(int iCore = 0; iCore < numCores - 1; iCore++){
				//MPI::COMM_WORLD.Recv(&recvWorkerDataArray, 3, MPI::DOUBLE, iCore, iCore);
				MPI_Recv(recvWorkerDataArray, 3, MPI_DOUBLE, iCore, iCore, MPI_COMM_WORLD, &status);
				thisWorkersM += recvWorkerDataArray[0];
				thisWorkersT += recvWorkerDataArray[1];
				thisWorkersS += recvWorkerDataArray[2];
			}

			// Save data
			savedData[0][K-100] = thisWorkersM;
			savedData[1][K-100] = thisWorkersT;
			savedData[2][K-100] = thisWorkersS;
			
		}else{
			double thisWorkersDataArray[3] = {thisWorkersM, thisWorkersT, thisWorkersS};
			MPI_Send(thisWorkersDataArray, 3, MPI_DOUBLE, (numCores - 1), myRank, MPI_COMM_WORLD);			
		}
	}
	
	// Have master write files
	if(myRank == numCores - 1){	
		// Print out table
		cout << "The following approximations to oscillatory integrals H(K) were computed with P = " << numCores << " cores." << endl;
		cout << endl << "---TABLE---" << endl;
		cout << "  K  " << "        H_{M,100}(K)        " << "        H_{T,100}(K)        " << "        H_{S,100}(K)" << endl;
		cout << "---------------------------------------------------------------------------------" << endl;  
		for(int i = 0; i < 10; i++){
			int index = i+1;
			int waveNumIWant = index*1000;
			cout << waveNumIWant << " , " << savedData[0][waveNumIWant - 100] << " , " << savedData[1][waveNumIWant - 100] << " , " << savedData[2][waveNumIWant - 100] << endl;	

		}	
		// output savedData to csv file
		for(int iRow = 0; iRow < 3; iRow++){
			for(int jCol = 0; jCol < 9901; jCol++){
				outputCSV << savedData[iRow][jCol] << ",";
			}
			outputCSV << "\n";
		}

	}
	// MPI Finalize	
	int MPI_Finalize();
}



    
    

