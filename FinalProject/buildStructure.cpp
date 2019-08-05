#include <iostream>
#include <mpi.h>

// include own files
//#include "functions.cpp"

// buildStructure() implementation
BuildStructureOut buildStructure(int NUM_LEVELS_M, int NUM_PARTITIONS_J, int NUM_KNOTS_r, int NUM_LEVELS_SERIAL_S, int numCores, int myRank, int nObsData, double offsetPercentage, Vector4d domainBoundaries, MatrixXd *data){
	// Set finest knot level - generlaize
	int finestKnotLevel = NUM_LEVELS_M - 1;
	int indexEndFinestKnotLevel = pow(NUM_PARTITIONS_J, finestKnotLevel) - 1;
	// Calculate quanities of interest
	VectorXi mLevels = VectorXi::LinSpaced(NUM_LEVELS_M, 0, NUM_LEVELS_M - 1);
	VectorXi nRegions = pow(NUM_PARTITIONS_J, mLevels.array());
	VectorXi cumulativeRegions(NUM_LEVELS_M); // Calculate cummulativeRegions vector with for-loop
        for(int i = 0; i < nRegions.array().size(); ++i){cumulativeRegions(i) = nRegions.head(i+1).sum();}
	// Calculate the number of knots in each direction
	long double sr = sqrt((double)NUM_KNOTS_r);
	int nKnotsX0, nKnotsX, nKnotsY0, nKnotsY;
	if((sr - floor(sr)) == 0){
		nKnotsX0 = sr; nKnotsX = sr; nKnotsY0 = sr; nKnotsY = sr;
	}else{
		nKnotsX0 = ceil(sr); nKnotsX = ceil(sr); nKnotsY0 = NUM_KNOTS_r/nKnotsX0; nKnotsY = NUM_KNOTS_r/nKnotsX;
	}	
	
	// More quantities of interest
	int nRegionsAtFinestLevelForEachWorker = nRegions(NUM_LEVELS_M - 1)/numCores; // Assume numCores is a power of J. Substract 1 b/c array indices begin @ 0
	// Calculate quantities needed for nTotalRegionsAssginedToEachWorker
	int  maxLevelOnASingleRow = (nRegions.array() <= numCores).count();
	VectorXi counter = VectorXi::LinSpaced(NUM_LEVELS_M - maxLevelOnASingleRow, 1, NUM_LEVELS_M - maxLevelOnASingleRow);
	int nTotalRegionsAssignedToEachWorker = maxLevelOnASingleRow + pow(NUM_PARTITIONS_J, counter.array()).sum();

	// Create IndexMatrix to store continuous index for all regions
	CreateIndexMatrixOut createdIndexMatrixOutput  = createIndexMatrix(NUM_LEVELS_M, NUM_PARTITIONS_J, numCores, myRank, NUM_LEVELS_SERIAL_S, nTotalRegionsAssignedToEachWorker, nRegions, cumulativeRegions); // LB: consider adding cumulativeRegions as an argument so only need to calculate once
	VectorXi indexMatrix = createdIndexMatrixOutput.indexMatrix;
	MatrixXi matrixOfRegionsAtFirstParallelLevel = createdIndexMatrixOutput.matrixOfRegionsAtFirstParallelLevel;

	//Find the last index within the indexMatrix corresponding to the finest level at which the knots are not set to the data.
	int indexOfFinestKnotLevelWithinIndexMatrix = (indexMatrix.array() < indexEndFinestKnotLevel).count(); // LB: 4/27 removed - 1 
	int nRowsWithRepeatedEntriesInIndexMatrix = (nRegions.array() < numCores).count();

	// Allocate space for objects created when building multi-resolution structure
	vector<MatrixXd> knots(nTotalRegionsAssignedToEachWorker);
	vector<MatrixXd> partitions(indexOfFinestKnotLevelWithinIndexMatrix + 1);
	//vector<vector<double>> predictionLocations(nRegionsAtFinestLevelForEachWorker);

	/* Construct the zeroth level */
	double xMin0 = domainBoundaries(0); double xMax0 = domainBoundaries(1);
	double yMin0 = domainBoundaries(2); double yMax0 = domainBoundaries(3);
	// Edge buffer added to xMax0 and yMax0 to include all observations on the zeroth level
	xMax0 = xMax0 + (offsetPercentage/2.)*(xMax0 - xMin0);
	yMax0 = yMax0 + (offsetPercentage/2.)*(yMax0 - yMin0);
	// Create the knots at the coarsest resolution	
	MatrixXd knots0  = createKnots(xMin0, xMax0, yMin0, yMax0, offsetPercentage, nKnotsX0, nKnotsY0);
	knots[0] = knots0;
	// Create partitions at the coarsest level
	MatrixXd partition0 = createPartition(xMin0, xMax0, yMin0, yMax0, NUM_PARTITIONS_J);
	partitions[0] = partition0;
	
	/*Loop through indexMatrix building the multi-resolution hierarchical structure*/	
	//cout << "nRowsWithRepeatedEntriesInIndexMatrix: " << nRowsWithRepeatedEntriesInIndexMatrix <<  endl;
	int indexCurrent, indexParent, thisIndexParentInIndexMatrix, thesePartitionBoundaries;
	Array3i foundParent; MatrixXd thisParentsPartition;
	VectorXd xMin, xMax, yMin, yMax; VectorXi foundChildren;
	double thisXMin, thisXMax, thisYMin, thisYMax;
	// Loop up to level M-1 creating partitions and placing knots
	for(int iRow = 1; iRow <= indexOfFinestKnotLevelWithinIndexMatrix; ++iRow){ // 4/28 LB: changed from nTotalRegionsAssignedToEachWorker(?)
		// Find this region's index
		indexCurrent = indexMatrix(iRow);
		// Find this region's parent and assign it to an int
		foundParent = findParent(indexCurrent, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
		indexParent = foundParent(2);
		// Find this region's parent's location in indexMatrix
		thisIndexParentInIndexMatrix = (indexMatrix.array() < indexParent).count();	
		// Get partition coordinates of parent
		thisParentsPartition = partitions[thisIndexParentInIndexMatrix];
		xMin = thisParentsPartition.col(0); xMax = thisParentsPartition.col(1); yMin = thisParentsPartition.col(2); yMax = thisParentsPartition.col(3);
		// Assign correct xMin,...,yMax for thisRegion
		foundChildren = findChildren(indexParent, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
		thesePartitionBoundaries = (foundChildren.array() < indexCurrent ).maxCoeff();
		thisXMin = xMin(thesePartitionBoundaries); thisXMax = xMax(thesePartitionBoundaries);
            	thisYMin = yMin(thesePartitionBoundaries); thisYMax = yMax(thesePartitionBoundaries);
		// Since computing over at levels up to M-1, create knots and partition the domain for next level
		// Create knots
		knots[iRow] = createKnots(thisXMin, thisXMax, thisYMin, thisYMax, offsetPercentage, nKnotsX, nKnotsY);			
		// Create partition
		partitions[iRow] = createPartition(thisXMin, thisXMax, thisYMin, thisYMax, NUM_PARTITIONS_J);		
	}	
		
	/*Prep for assigning data at the finest resolution*/
	// LB 4/28: Make special case for serial where entire domain is assigned to one core

	VectorXi thisWorkerParallelAssignment = matrixOfRegionsAtFirstParallelLevel.col(myRank);
	int beginParallelAssignment = (indexMatrix.array() < thisWorkerParallelAssignment(0)).count();
	int endParallelAssignment = (indexMatrix.array() < thisWorkerParallelAssignment(thisWorkerParallelAssignment.size() - 1)).count();
	// Check to make sure that they're not exactly the same so that we actually extract the correct subset of the vector
	if(beginParallelAssignment == endParallelAssignment){endParallelAssignment+=1;}	
	vector<MatrixXd> thisWorkersAssignmentPartitions(partitions.begin() + beginParallelAssignment, partitions.begin() + endParallelAssignment);
	int thisIndexInIndexMatrix;
	double workerXMin, workerXMax, workerYMin, workerYMax;
	workerXMin = thisWorkersAssignmentPartitions[0].col(0).minCoeff(); workerXMax = thisWorkersAssignmentPartitions[0].col(1).maxCoeff();
	workerYMin = thisWorkersAssignmentPartitions[0].col(2).minCoeff(); workerYMax = thisWorkersAssignmentPartitions[0].col(3).minCoeff();
	// Loop through all of this workers parallel assignmet regions to find the entire subset of the domain this worker is reponsible for
	for(int i = 1; i < endParallelAssignment - beginParallelAssignment; ++i){
		if(workerXMin > thisWorkersAssignmentPartitions[i].col(0).minCoeff()){workerXMin = thisWorkersAssignmentPartitions[i].col(0).minCoeff();}
		if(workerXMax < thisWorkersAssignmentPartitions[i].col(1).maxCoeff()){workerXMax = thisWorkersAssignmentPartitions[i].col(1).maxCoeff();}
		if(workerYMin > thisWorkersAssignmentPartitions[i].col(2).minCoeff()){workerYMin = thisWorkersAssignmentPartitions[i].col(2).minCoeff();}
		if(workerYMax < thisWorkersAssignmentPartitions[i].col(3).maxCoeff()){workerYMax = thisWorkersAssignmentPartitions[i].col(3).maxCoeff();}
	}

	// Loop through the data, selecting only the subset needed by this worker
	MatrixXd processedData(nObsData,3); // Allocate space for potentially up to all the data, then delete rows not needed for thie worker after loop
	int countNumObsForThisWorker = 0;
	for(int iDataRow = 0; iDataRow < nObsData; ++iDataRow){
		// If this data is within the XY-coords assigned to this worker
		if((data->col(0)(iDataRow) >= workerXMin) & (data->col(0)(iDataRow) <= workerXMax) & (data->col(1)(iDataRow) >= workerYMin) & (data->col(1)(iDataRow) <= workerYMax)){
				// Then assign this loaded data to the processedData for this worker and increment the number of obs assigned to this worker
				processedData.row(countNumObsForThisWorker) = data->row(iDataRow);
				countNumObsForThisWorker++;
		}
	}
	//cout << "countNumObsForThisWorker: " << countNumObsForThisWorker << endl; 
	// Resize the processedData Matrix to only contain rows with data assigned to them
	processedData.conservativeResize(countNumObsForThisWorker,3);

	// Loop through the finest resolution (i.e., level m = M), assigning the data to the knots	
	vector<vector<double>> outputData(nRegionsAtFinestLevelForEachWorker);
	vector<double> thisRegionsKnotsX, thisRegionsKnotsY, thisRegionsObs;
	for(int iRow = indexOfFinestKnotLevelWithinIndexMatrix + 1; iRow < nTotalRegionsAssignedToEachWorker; ++iRow){	
		// Find this region's index
		indexCurrent = indexMatrix(iRow);
		// Find this region's parent and assign it to an int
		foundParent = findParent(indexCurrent, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
		indexParent = foundParent(2);
		// Find this region's parent's location in indexMatrix
		thisIndexParentInIndexMatrix = (indexMatrix.array() < indexParent).count();	
		// Get partition coordinates of parent
		thisParentsPartition = partitions[thisIndexParentInIndexMatrix];
		xMin = thisParentsPartition.col(0); xMax = thisParentsPartition.col(1); yMin = thisParentsPartition.col(2); yMax = thisParentsPartition.col(3);
		// Assign correct xMin,...,yMax for thisRegion
		foundChildren = findChildren(indexParent, NUM_PARTITIONS_J, nRegions, cumulativeRegions);
		thesePartitionBoundaries = (foundChildren.array() < indexCurrent ).maxCoeff();
		thisXMin = xMin(thesePartitionBoundaries); thisXMax = xMax(thesePartitionBoundaries);
            	thisYMin = yMin(thesePartitionBoundaries); thisYMax = yMax(thesePartitionBoundaries);
		// Collect the data needed for this region // 4/29 LB: This is potentially a huge area for optimization
		
		for(int j = 0; j < countNumObsForThisWorker; ++j){
			// If this data observation has coordinates contained within this region, assign it locations to thisRegionsKnotsMatrix
			
			if((processedData(j, 0) >= thisXMin) & (processedData(j, 0) <= thisXMax) & (processedData(j, 1) >= thisYMin) & (processedData(j, 1) <= thisYMin)){
				thisRegionsKnotsX.push_back(processedData(j,0));
				thisRegionsKnotsY.push_back(processedData(j,1)); 			
				thisRegionsObs.push_back(processedData(j,2));
			}	
			
			
		}
		// Assign knots and observations to their respective data structures
		Map<VectorXd> thisRegionsKnotsXVec(thisRegionsKnotsX.data(), thisRegionsKnotsX.size());	
		Map<VectorXd> thisRegionsKnotsYVec(thisRegionsKnotsY.data(), thisRegionsKnotsY.size());
		knots[iRow].resize(thisRegionsKnotsY.size(), 2);
		knots[iRow] << thisRegionsKnotsXVec, thisRegionsKnotsYVec;
		outputData[iRow - (indexOfFinestKnotLevelWithinIndexMatrix + 1)] = thisRegionsObs;
		// Clear vectors for next time through the loop
		thisRegionsKnotsX.clear(); thisRegionsKnotsY.clear(); thisRegionsObs.clear();
	}
	






	// Allocate space for objects used in building the structure Define final output structure
	BuildStructureOut buildStructureOutput  = {knots, partitions, nRegions, outputData};
	// Return the final output
	return buildStructureOutput;
}


