#include "matrix_inversion.hpp"
#include "utils.hpp"
#include <iostream>
#include <vector>
#include <iomanip>

void displayMatrixInversionMenu() {
    clearScreen();
    cout << "\n";
    cout << "\t ==============================" << endl;
    cout << "\t|                              |" << endl;
    cout << "\t       ";
    printText("Matrix Inversion", 0, 2, true);
    cout << endl;
    cout << "\t      ";
    printText("> (n x n) Matrix <", 0, 2, true);
    cout << endl;
    cout << "\t|                              |" << endl;
    cout << "\t ============================== " << endl << endl;
}

vector<vector<double>> inputMatrix(int n, int m) {
    vector<vector<double>> mat(n, vector<double>(m));
    cout << "[?] - Enter the matrix (" << n << "x" << n  << "): " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            cin >> mat[i][j];
        }
    }
    return mat;
}

vector<vector<double>> invert(vector<vector<double>> A) {
    int n = A.size();
    vector<vector<double>> inv(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i) {
        inv[i][i] = 1;
    }

    for (int i = 0; i < n; ++i) {
        double pivot = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = 0; k < n; ++k) {
                    A[j][k] -= factor * A[i][k];
                    inv[j][k] -= factor * inv[i][k];
                }
            }
        }
    }

    return inv;
}

bool isInvertible(const vector<vector<double>>& mat) {
    int n = mat.size();
    vector<vector<double>> temp = mat;
    double det = 1;

    for (int i = 0; i < n; ++i) {
        if (temp[i][i] == 0) {
            bool swapped = false;
            for (int j = i + 1; j < n; ++j) {
                if (temp[j][i] != 0) {
                    swap(temp[i], temp[j]);
                    det *= -1;
                    swapped = true;
                    break;
                }
            }
            if (!swapped) return false;
        }

        det *= temp[i][i];
        for (int j = i + 1; j < n; ++j) {
            double factor = temp[j][i] / temp[i][i];
            for (int k = i; k < n; ++k) {
                temp[j][k] -= factor * temp[i][k];
            }
        }
    }

    return det != 0;
}

void printMatrix(const vector<vector<double>>& mat) {
    int n = mat.size();
    int m = mat[0].size();
    // Print top border
    cout << " ";
    for (int j = 0; j < m; ++j) {
        cout << "---------";
    }
    cout << "-" << endl;

    for (int i = 0; i < n; ++i) {
        cout << "| ";
        for (int j = 0; j < m; ++j) {
            cout << setw(8) << mat[i][j] << " ";
        }
        cout << "|" << endl;
    }
    
    // Print bottom border
    cout << " ";
    for (int j = 0; j < m; ++j) {
        cout << "---------";
    }
    cout << "-" << endl;

}

void invertMatrix() {
    displayMatrixInversionMenu();

    int n;

    cout << "[?] - Enter n: ";
    cin >> n;
    
    vector<vector<double>> mat = inputMatrix(n, n);

    cout << endl << endl;
    cout << "\t -- Input Matrix --- " << endl << endl;
    printMatrix(mat);

    cout << endl << endl;

    if(!isInvertible(mat)) {
        cout << "\t --- Matrix is not invertible ---" << endl;
    }else{
        vector<vector<double>> inv = invert(mat);
        cout << "\t --- Inverted Matrix ---" << endl << endl;
        printMatrix(inv);
    }

    cout << endl << endl << "Press any key to continue..." << endl;
    getChar();
}
