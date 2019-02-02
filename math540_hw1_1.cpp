#include <iostream>
#include <vector>
#include <iomanip> // std::setprecision
#include <chrono>
#include <ratio>
#include "functions.cpp"
using namespace std;
using namespace std::chrono;

int main()
{
    // Initialize timer
    using clock = chrono::steady_clock;
    // Set precision for 14 decimal places - better way?
    cout.precision(16);
    // Get N from userInput();
    unsigned int N = getN();
    // Start Timer
    auto t1 = std::chrono::high_resolution_clock::now();
    // Create A matrix with createA() function
    vector<vector<double> > Amat = createA(N);
    vector<vector<double> > Amat2 = createA(N/2);
    // Check outputs from createA()
    //cout << "Amat:" << endl;
    //printMatrix(Amat);

    // Create B matrix with createB() function
    Foo createBOutput = createB(N);
    vector<vector<double> > Bmat = createBOutput.Bmatrix;
    int countN = createBOutput.countN;

    Foo createBOutput2 = createB(N/2);
    vector<vector<double> > Bmat2 = createBOutput2.Bmatrix;
    int countN2 = createBOutput2.countN;

    
    // Check outputs from createB()
    //cout << "count N is: " << countN << endl;
    //cout << "Bmat:" << endl;
    //printMatrix(Bmat);

    // Create AB matrix with createAB()
    vector<vector<double> > ABmat =  createAB(N, Amat, Bmat);
    vector<vector<double> > ABmat2 =  createAB(N/2, Amat2, Bmat2);

    // Check outpouts from createAB()
    //cout << "ABmat:" << endl;
    //printMatrix(ABmat);

    // Create C matrix with createC()
    vector<vector<double> > Cmat = createC(N, ABmat);
    vector<vector<double> > Cmat2 = createC(N/2, ABmat2);
    // Check outputs from createC() 
    //cout << "Cmat:" << endl;
    //printMatrix(Cmat);
    // cout << Cmat[99][99] << endl; // Checking C_N(N/2, N/2)
    // cout << Cmat[199][199] << endl; // Checking C_N(N,N)

    //End Timer
    auto t2 = std::chrono::high_resolution_clock::now();
    // floating-point duration: no duration_cast needed
    std::chrono::duration<double, std::milli> fp_s = t2 - t1;


    // Print out information for problem 1
    cout << "The value of N is: " << N << endl;
    cout << "The value of C_N(N/2, N/2) is: " <<  Cmat[N/2 - 1][N/2 -1] << endl; // Check this !!!
    cout << "The value of C_N(N, N) is: " << Cmat[N-1][N-1] << endl; // Check this as well !!!
    cout << "The value of count(N) is: " << countN << endl;
    cout << "The value of T_{1, Q_1}(N) is: "  << fp_s.count() << " miliseconds" << endl;

}