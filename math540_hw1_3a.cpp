#include <iostream>
#include <vector>
#include <iomanip> // std::setprecision
#include <chrono>
#include <ratio>
#include <fstream>
#include "functions.cpp"
using namespace std;
using namespace std::chrono;

int main()
{
    // Initialize timer
    using clock = chrono::steady_clock;
    // Set precision for 14 decimal places - better way?
    cout.precision(16);
    // Delcare outputfiles
    ofstream AmatOutput; // for Amat
    ofstream BmatOutput; // for Bmat
    // Open the output files for Am
    AmatOutput.open("Amat.bin", ios::out | ios::binary);
    BmatOutput.open("Bmat.bin", ios::out | ios::binary);
    // Get N from userInput();
    int N = getN();
    
    //cout << "Enter a even positive integer: " << endl;
    //int N;
    //cin >> N;
    
    // Start Timer
    auto t1 = std::chrono::high_resolution_clock::now();

    // Check to make sure output files are open before creating A, B matrices
    if (!(AmatOutput.is_open() && BmatOutput.is_open()))
    {
        cout << "Error creating output files." << endl;
        return -1;
    }
    // Create A matrix with createA() function
    vector<vector<double> > Amat = createA(N);
    // Create B matrix with createB() function
    Foo createBOutput = createB(N);
    vector<vector<double> > Bmat = createBOutput.Bmatrix;
    int countN = createBOutput.countN;

    // Declare array to put contents of Amat, Bmat
    //double AmatArray[N][N];
    //double BmatArray[N][N];

    double **AmatArray = new double*[N];
    double **BmatArray = new double*[N];

    for (int i = 0; i< N; i++){
        AmatArray[i] = new double[N];
        BmatArray[i] = new double[N];
    }

    // Assign contents of Amat, Bmat to AmatArray, BmatArray resp.
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            AmatArray[i][j] = Amat[i][j];
            BmatArray[i][j] = Bmat[i][j];
            //AmatArray[i][j] = Amat[i][j];
            //BmatArray[i][j] = Bmat[i][j];
        }
    }

    // Write contents of AmartArray, BmatArray to binary files Amat.bin, Bmat.bin
    for (int i = 0; i < N; i++)
    {
        AmatOutput.write((char *)AmatArray[i], N * sizeof(double));
        BmatOutput.write((char *)BmatArray[i], N * sizeof(double));
    }

    // Deallocate AmatArray2, BmatArray2 memory
    // after writing to binary files
    for (int i = 0; i< N; i++)
    {
        delete [] AmatArray[i];
        delete []BmatArray[i];
    }
    delete [] AmatArray;
    delete [] BmatArray;

    // Close the output files
    AmatOutput.close();
    BmatOutput.close();

    //End Timer
    auto t2 = std::chrono::high_resolution_clock::now();
    // floating-point duration: no duration_cast needed
    std::chrono::duration<double, std::milli> fp_s = t2 - t1;

    // Output to console
    cout << "The value of N is: " << N << endl;
    cout << "The value of count(N) is: " << countN << endl;
    cout << "The value of T_{1, Q_{3a}}(N) is: " << fp_s.count() << " milliseconds" << endl;


    return 0;
}