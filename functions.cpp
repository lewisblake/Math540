// functions.cpp
#include <iostream>
//#include <vector>
#include <cmath>
#include "functions.h"
using namespace std;

// getN() function implementation
int getN(){
    // Declare unsgined int N
    int N;
    // Get user input with do-while loop
    do
    {   // Prompt user to give N in command window
        cout << "Enter an even positive integer, N: " << endl;
        cin >> N; // Get user input
        // Check to make sure inputted N is not odd
        if(N %2 == 1) {cout << "Supplied integer is odd, please try again." << endl;}
    } while (N % 2 == 1); // Repeat prompt while input is odd
    // Return N from user input
    return N;
}


// createA() function implementation
std::vector<std::vector<double> > createA(int N){
    // Initalize space for vector
    vector<vector<double> > Amat(N, vector<double>(N));
    // For-loop to fill Amat
    // faster to loop to N rather than using .size() function?
    for (double i = 0; i < Amat.size(); i++) // declare i, j as doubles for arithmetic
    {
        for (double j = 0; j < Amat.size(); j++)
        {
            Amat[i][j] = (1 / ((i + 1) + (j + 1) - 1));
        }
    }
    // Return the A matrix
    return Amat;
}
/*
// Create a struct to wrap output from createB() into
struct Foo {
    vector<vector<double> > Bmatrix;
    int countN;
}; */

// createB() function implementation
//std::vector<std::vector<std::vector<double> > >
Foo createB(int N){
    // Initalize space for output vector
    //std::vector<vector<vector<double> > > createBOutput(2);
    // Initalize space for B matrix
    vector<vector<double> > Bmat(N, vector<double>(N));
    // Initalize countN
    unsigned int countN = 0;
    // Initialize upperBound
    double upperBound = pow(10.0,100.0);
    // For-loop to fill Bmat
    for (double i = 0; i < Bmat.size(); i++){
        for(double j = 0; j < Bmat.size(); j++){
            // Define alpha's and beta's
            // Arrays start @ 0 in C++
            double alpha1 = (double)N + (i+1) -1; double alpha2 = N + (j+1) - 1; double alpha3 = (i+1) + (j+1) - 2;
            double beta1 = (double)N - (j+1); double beta2 = (double)N - (i+1); double beta3 = (i+1) - 1;
            // Calculate b_ij
            double firstTerm = exp(lgamma(alpha1 + 1) - lgamma(beta1 + 1) - lgamma(alpha1 - beta1 +1));
            double secondTerm = exp(lgamma(alpha2 + 1) - lgamma(beta2 + 1) - lgamma(alpha2 - beta2 + 1));
            double thirdTerm = exp(lgamma(alpha3 + 1) - lgamma(beta3 + 1) - lgamma(alpha3 - beta3 + 1));
            double b = pow(-1.0, (i+1) + (j+1)) * ((i+1) + (j+1) -1) * firstTerm * secondTerm * pow(thirdTerm, 2.0);

           if(abs(b) <= upperBound){
               Bmat[i][j] = b; // Fill Bmat
           }
           else{
               double firstTermTilde = pow((sin(alpha1)/(cos(beta1) * tan(alpha1 + beta1))), 3.0);
               double secondTermTilde = pow((sin(alpha2)/(cos(beta2) * tan(alpha2 + beta2))), 4.0);
               double thirdTerm = pow((sin(alpha3)/(cos(beta3) * tan(alpha3 + beta3))), 5.0);
               double btilde = pow(-1.0, (i+1) + (j+1)) * pow((i+1) + (j+1) - 1, 2.0)  * firstTermTilde * secondTermTilde * thirdTerm;
               Bmat[i][j] = btilde;
               // Increment countN
               countN++;
           }
        }
    };
    Foo createBOutput = {Bmat, countN};
    // Return the vector containing Bmat and countN
    return createBOutput;
}

// createAB() function implementation
vector<vector<double> > createAB(int N, vector<vector<double> > Amat, vector<vector<double> > Bmat){
    // Initalize space for AB matrix
    vector<vector<double> > ABmat(N, vector<double>(N));

    // Perform matrix multiplications the old-fashioned way...
    for (unsigned int iRow = 0; iRow < N; iRow++){
        for (unsigned int jColumn = 0; jColumn < N; jColumn ++){
            double theSum = 0;
            for (unsigned int k = 0; k < N; k++){
                theSum = theSum + Amat[iRow][k] * Bmat[k][jColumn];
            }
            // Fill ABmat
            ABmat[iRow][jColumn] = theSum;
        }
    }

    return ABmat;
}

// createC() function implementation
vector<vector<double> > createC(int N, vector<vector<double> > ABmat){
    double upperBound = pow(10.0,100.0);
    // Initialize C matrix
    vector<vector<double> > Cmat(N, vector<double>(N));

    // For-loop to fill Cmat
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(ABmat[i][j] <= upperBound){
                Cmat[i][j] = ABmat[i][j];

            }
            else if(abs(ABmat[i][j]) > upperBound || isnan(ABmat[i][j])){
                Cmat[i][j] = (double)(i^j);
            }
        }
    }
    return Cmat;
}

// printMatrix() function implementation
void printMatrix(vector<std::vector<double> > myMatrix){
    cout << endl;
    for (int i = 0; i < myMatrix.size(); i++)
    {
        for (int j = 0; j < myMatrix.size(); j++)
        {
            cout << myMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

}