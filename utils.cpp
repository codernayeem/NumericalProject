#include "utils.hpp"
#include <iostream>
#include <vector>

// Function to print a matrix
void printMatrix(const std::vector<std::vector<double>>& matrix) {
    std::cout << "Matrix:\n";
    for (const auto& row : matrix) {
        for (double elem : row) {
            std::cout << elem << "\t";
        }
        std::cout << "\n";
    }
}

// Function to ask user for a choice (integer input)
int askChoice(const std::string& prompt, int min, int max) {
    int choice;
    do {
        std::cout << prompt << " (" << min << " - " << max << "): ";
        std::cin >> choice;
        if (choice < min || choice > max) {
            std::cout << "Invalid choice. Try again.\n";
        }
    } while (choice < min || choice > max);
    return choice;
}

// Function to get input for a matrix
std::vector<std::vector<double>> inputMatrix(int rows, int cols) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    std::cout << "Enter the matrix elements row by row:\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << "Element [" << i + 1 << "][" << j + 1 << "]: ";
            std::cin >> matrix[i][j];
        }
    }
    return matrix;
}

// Function to get input vector
std::vector<double> inputVector(int size) {
    std::vector<double> vec(size);
    std::cout << "Enter the elements of the vector:\n";
    for (int i = 0; i < size; ++i) {
        std::cout << "Element [" << i + 1 << "]: ";
        std::cin >> vec[i];
    }
    return vec;
}

// Utility to print a vector
void printVector(const std::vector<double>& vec) {
    std::cout << "Vector: [ ";
    for (double elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << "]\n";
}
