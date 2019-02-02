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
    // Get N from userInput();
    int N = getN();

    // Initiailze timer
    using clock = chrono::steady_clock;
    // Start Timer
    auto t1 = std::chrono::high_resolution_clock::now();

    // Set precision for 14 decimal places - better way?
    cout.precision(16);

    // Open Amat.bin and Bmat.bin in one line
    ifstream AmatInput("Amat.bin", ios::in | ios::binary);
    ifstream BmatInput("Bmat.bin", ios::in | ios::binary);
    // Initialize space for Amat and Bmat
    vector<vector<double> > Amat(N, vector<double>(N));
    vector<vector<double> > Bmat(N, vector<double>(N));
    
    // Initialize arrays to put contents of AmatInput, BmatInput
    //double AmatArray[N][N];
    //double BmatArray[N][N];

    double **AmatArray = new double*[N];
    double **BmatArray = new double*[N];
    for (int i = 0; i< N; i++){
        AmatArray[i] = new double[N];
        BmatArray[i] = new double[N];
    }

    // Read the binary file the same way you wrote it
    // Assign read contents into AmatArray, BmatArray
    for (int i = 0; i < N; i++)
    {
        AmatInput.read((char *)AmatArray[i], N * sizeof(double));
        BmatInput.read((char *)BmatArray[i], N * sizeof(double));
    }

    // Assign contents of AmatArray, BmatArray into Amat, Bmat
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Amat[i][j] = AmatArray[i][j];
            Bmat[i][j] = BmatArray[i][j];
        }
    }
    
    // Deallocate AmatArray, BmatArray memory
    // after writing to binary files
    for (int i = 0; i< N; i++)
    {
        delete [] AmatArray[i];
        delete []BmatArray[i];
    }
    delete [] AmatArray;
    delete [] BmatArray;

    // Create AB matrix with createAB()
    vector<vector<double> > ABmat = createAB(N, Amat, Bmat);

    // Create C matrix with createC()
    vector<vector<double> > Cmat = createC(N, ABmat);

    //End Timer
    auto t2 = std::chrono::high_resolution_clock::now();
    // floating-point duration: no duration_cast needed
    std::chrono::duration<double, std::milli> fp_s = t2 - t1;

    
    // Output to console
    cout << "The value of N is: " << N << endl;
    cout << "The value of C_N(N/2, N/2) is: " <<  Cmat[N/2 - 1][N/2 -1] << endl; // Check this !!!
    cout << "The value of C_N(N, N) is: " << Cmat[N-1][N-1] << endl; // Check this as well !!!
    cout << "The value of T_{1, Q_{3b}}(N) is: "  << fp_s.count() << " milliseconds" << endl;
    

    return 0;
}