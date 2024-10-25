#include "linear_eqs.hpp"
#include<iostream>
#include<cmath>
#include<vector>
using namespace std;
const double tolerance = 1e-6;
const int iterations = 1000;

void choice()
{
    cout << "Numerical Methods Console Application\n";
    cout << "1. Jacobi Iterative Method\n";
    cout << "2. Gauss-Siedel Iterative Method\n";
    cout << "3. Gauss Elimination\n";
    cout << "4. Gauss-Jordan Elimination\n";
    cout << "5. Exit\n";
}

void input(vector<vector<double>> &A, vector<double> &b, int n)
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

bool Jacobi(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int n){
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
bool Gauss_Seidel(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int n){
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
void Gauss_Elimination(){};
void Gauss_Jordan_Elimination(){};

void display_answers(vector<double> &x)
{
    int n = x.size();

    cout << "Answers ->" << endl; 
    for(int i=0; i<n; i++){
        cout << "Root " << i+1 << " = " << x[i] << endl; 
    }
}


void solveLinearEquations() {
    std::cout << "Solving Linear Equations (Jacobi, Gauss-Seidel, etc.)...\n";
    std::cout << "Enter Your Choice\n";
    choice();

    int option, n;
    cin >> option;

    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);
    vector<double> x(n, 0.0);

    input(A, b, n);

    if(option == 1){
        bool solvable = Jacobi(A, b, x, n);
        if(solvable){
            display_answers(x);
        }
        else{
            cout << "Equations did not converge " << endl;
        }
    }
    else if(option == 2){
        bool solvable = Gauss_Seidel(A, b, x, n);
        if(solvable){
            display_answers(x);
        }
        else{
            cout << "Equations did not converge " << endl;
        }
    }

    return;
}