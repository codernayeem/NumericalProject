#include "linear_eqs.hpp"
#include "utils.hpp"
#include<iostream>
#include<cmath>
#include<vector>
#include<iomanip>
using namespace std;

void displayLinearEquationsMenu() {
    clearScreen();
    cout << "\n";
    cout << "\t ==============================" << endl;
    cout << "\t|                              |" << endl;
    cout << "\t      ";
    printText("Solve Linear Equations", 0, 2, true);
    cout << endl;
    cout << "\t|                              |" << endl;
    cout << "\t ============================== " << endl << endl;
}

void choice()
{
    cout << "\t>> Choose a option:\n\n";
    cout << "\t1. Jacobi Iterative Method\n";
    cout << "\t2. Gauss-Siedel Iterative Method\n";
    cout << "\t3. Gauss Elimination\n";
    cout << "\t4. Gauss-Jordan Elimination\n";
    cout << "\t5. LU Factorization\n";
    cout << "\t6. Go Back to Main Menu\n";
    cout << "\n";
    printText("\t<< Enter your choice (1-6) >>", 0, 2);
    cout << endl;
}

void inputMatrix(vector<vector<double>> &A, vector<double> &b, int n, int iterations = 1000, double tolerance = 1e-6)
{
    cout << "Enter the coefficient matrix A" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the right-hand side vector b:" << endl;
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }
}

bool Jacobi(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int n, int iterations = 1000, double tolerance = 1e-6){
    vector<double> temp(n);
    vector<double> error(n);
    int itr = iterations;

    while(itr--){
        temp = x;

        for (int i = 0; i < n; i++) {
            double sum = 0;

            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i][j] * temp[j];
                }
            }

            x[i] = (b[i] - sum) / A[i][i];
            error[i] = abs(x[i] - temp[i]);
        }
    }
    double total_error = 0;
    for(auto it: error){
        total_error += it;
    }
    if(total_error < tolerance) return true;

    return false;
}

bool Gauss_Seidel(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int n, int iterations = 1000, double tolerance = 1e-6){
    vector<double> temp(n);
    vector<double> error(n);
    int itr = iterations;
    
    while(itr--){
        temp = x;

        for (int i = 0; i < n; i++) {
            double sum = 0;

            for (int j = 0; j < n; j++) {
                if(i == j) continue;
                sum += A[i][j] * x[j];
            }

            x[i] = (b[i] - sum) / A[i][i];
            error[i] = abs(x[i] - temp[i]);
        }

    }
    double total_error = 0;
    for(auto it: error){
        total_error += it;
    }
    if(total_error < tolerance) return true;

    return false;
}

vector<double> Gauss_Elimination(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    printText("Gauss Eliminated Form of the Matrix\n", 0, 2, true);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            std::cout << std::fixed << std::setprecision(4) << A[i][j] << " ";
        }
        cout << endl;
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return x;
}

void Gauss_Jordan_Elimination(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        double pivot = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= pivot;
        }
        b[i] /= pivot;

        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }
    printText("Gauss-Jordan Eliminated Form of the Matrix\n", 0, 2, true);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << fixed << setprecision(4);
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

void display_answers(vector<double> &x)
{
    int n = x.size();

    printText("Answers\n", 0, 2, true);
    for(int i=0; i<n; i++){
        cout << "Root " << i+1 << " = " << x[i] << endl; 
    }
}

void solveLinearEquations() {
    displayLinearEquationsMenu();
    choice();
    
    int option;
    while(true) {
        option = getChar();
        if (option >= '1' && option <= '8') {
            option -= '0';
            break;
        }
    }

    if(option == 6) {
        return;
    }

    cout << endl;
    switch (option) {
        case 1:
            cout << "[+] - Jacobi Iterative Method selected." << endl;
            break;
        case 2:
            cout << "[+] - Gauss-Siedel Iterative Method." << endl;
            break;
        case 3:
            cout << "[+] - Gauss Elimination." << endl;
            break;
        case 4:
            cout << "[+] - Gauss-Jordan Elimination." << endl;
            break;
        case 5:
            cout << "[+] - LU Factorization." << endl;
            break;
    }

    int n;
    cout << "[?] - Enter the number of variables: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    vector<double> x(n, 0.0);
    bool solvable;

    switch (option){
        case 1:
            inputMatrix(A, b, n);
            solvable = Jacobi(A, b, x, n);
            if(solvable){
                display_answers(x);
            }
            else{
                cout << "[!] - Equations did not converge " << endl;
            }
            break;
        case 2:
            inputMatrix(A, b, n);
            solvable = Gauss_Seidel(A, b, x, n);
            if(solvable){
                display_answers(x);
            }
            else{
                cout << "[!] - Equations did not converge " << endl;
            }
            break;
        case 3:
            inputMatrix(A, b, n);
            x = Gauss_Elimination(A, b);
            display_answers(x);
            break;
        case 4:
            inputMatrix(A, b, n);
            Gauss_Jordan_Elimination(A, b);
            display_answers(b);
            break;
        case 5:
            cout << "[+] - Coming SOON" << endl;
            break;
    }

    cout << endl << "Press any key to continue..." << endl;
    getChar();
}