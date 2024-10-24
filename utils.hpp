#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <iostream>

// Function to print a matrix
void printMatrix(const std::vector<std::vector<double>>& matrix);

// Function to ask user for a choice (integer input)
int askChoice(const std::string& prompt, int min, int max);

// Function to get input for a matrix
std::vector<std::vector<double>> inputMatrix(int rows, int cols);

// Function to get an input vector
std::vector<double> inputVector(int size);

// Utility to print a vector
void printVector(const std::vector<double>& vec);

#endif
