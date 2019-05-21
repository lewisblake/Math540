// functions.h
#include<iostream>
#include <cmath> // need to incluce <cmath> here?
using namespace std;
#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_


// Foo struct declaration -  move to functions.cpp??
struct Foo {
    vector<vector<double> > Bmatrix;
    int countN;
};
// getN() function
int getN(); // Prompts user input to get N

// createA() function
std::vector<std::vector<double> > createA(int N); // createA() prototype, it's declaration

// createB() function
//std::vector<std::vector<std::vector<double> > >
Foo createB(int N); // createB() prototype, it's declaration

// createAB() function
vector<vector<double> > createAB(int N, vector<vector<double> > Amat, vector<vector<double> > Bmat);

// createC() function prototype
vector<vector<double> > createC(int N, vector<vector<double> > ABmat);

// printMatrix() function prototype
void printMatrix(vector<std::vector<double> > myMatrix);

#endif
