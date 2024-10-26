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

void inputMatrix(vector<vector<double>> &A, vector<double> &b, int n, int iterations = 1000, double tolerance = 1e-4)
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


bool Jacobi(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int n, int iterations = 1000, double tolerance = 1e-4){
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
        double total_error = 0;
        for(auto it: error){
            total_error += it;
        }
        if(total_error / n < tolerance) {
            cout << "[+] - Equations converged in " << iterations - itr << " iterations." << endl;
            return true;
        }
    }

    return false;
}

bool Gauss_Seidel(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int n, int iterations = 1000, double tolerance = 1e-4){
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
        double total_error = 0;
        for(auto it: error){
            total_error += it;
        }
        if(total_error / n < tolerance){
            cout << "[+] - Equations converged in " << iterations - itr << " iterations." << endl;
            return true;
        }
    }
    return false;
}

int Gauss_Elimination(vector<vector<double>>& A, vector<double>& B, vector<double>& x) {
    int n = A.size();
    double EPSILON = 1e-9;

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int i = col + 1; i < n; ++i) {
            if (fabs(A[i][col]) > fabs(A[pivot][col])) {
                pivot = i;
            }
        }

        if (fabs(A[pivot][col]) < EPSILON) {
            continue;
        }

        swap(A[col], A[pivot]);
        swap(B[col], B[pivot]);

        for (int i = col + 1; i < n; ++i) {
            double factor = A[i][col] / A[col][col];
            for (int j = col; j < n; ++j) {
                A[i][j] -= factor * A[col][j];
            }
            B[i] -= factor * B[col];
        }
    }

    for (int i = 0; i < n; ++i) {
        bool allZero = true;
        for (int j = 0; j < n; ++j) {
            if (fabs(A[i][j]) > EPSILON) {
                allZero = false;
                break;
            }
        }
        if (allZero && fabs(B[i]) > EPSILON) {
            return -1;
        }
        if (allZero && fabs(B[i]) < EPSILON) {
            return 1;
        }
    }
    
    printText("Gauss Eliminated Form of the Matrix\n", 0, 2, true);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << fixed << setprecision(4) << A[i][j] << " ";
        }
        cout << endl;
    }

    for (int i = n - 1; i >= 0; --i) {
        x[i] = B[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return 0;
}

int Gauss_Jordan_Elimination(vector<vector<double>>& A, vector<double>& B) {
    int n = A.size();
    int m = A[0].size();
    double EPSILON = 1e-9;

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int i = col + 1; i < n; ++i) {
            if (fabs(A[i][col]) > fabs(A[pivot][col])) {
                pivot = i;
            }
        }

        if (fabs(A[pivot][col]) < EPSILON) {
            continue;
        }

        swap(A[col], A[pivot]);
        swap(B[col], B[pivot]);

        double pivotValue = A[col][col];
        for (int j = 0; j < m; ++j) {
            A[col][j] /= pivotValue;
        }
        B[col] /= pivotValue;

        for (int i = 0; i < n; ++i) {
            if (i != col) {
                double factor = A[i][col];
                for (int j = 0; j < m; ++j) {
                    A[i][j] -= factor * A[col][j];
                }
                B[i] -= factor * B[col];
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        bool allZero = true;
        for (int j = 0; j < m; ++j) {
            if (fabs(A[i][j]) > EPSILON) {
                allZero = false;
                break;
            }
        }
        if (allZero && fabs(B[i]) > EPSILON) {
            return -1;
        }
        if (allZero && fabs(B[i]) < EPSILON) {
            return 1;
        }
    }

    printText("Gauss-Jordan Eliminated Form of the Matrix\n", 0, 2, true);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << fixed << setprecision(4);
            cout << abs(A[i][j]) << " ";
        }
        cout << endl;
    }
    return 0;
}



bool LU_Factorization(vector<vector<double>> &A, vector<double> &b, vector<double> &x){
    int n = A.size();

    if(n > 4){
        cout << "LU Factorization is not available for n > 4" << endl;
        return false;
    }

    if(n == 1){
        if(A[0][0] == 0){
            cout << "Determinant is zero, cannot solve the equations" << endl;
            return false;
        }
        x[0] = b[0] / A[0][0];
    }else if(n == 2){
        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if(det == 0){
            cout << "Determinant is zero, cannot solve the equations" << endl;
            return false;
        }

        x[0] = (A[1][1] * b[0] - A[0][1] * b[1]) / det;
        x[1] = (A[0][0] * b[1] - A[1][0] * b[0]) / det;
    }else if(n == 3){
        double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
        if(det == 0){
            cout << "Determinant is zero, cannot solve the equations" << endl;
            return false;
        }

        x[0] = (b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2]) + A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2])) / det;
        x[1] = (A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2]) - b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0])) / det;
        x[2] = (A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1]) - A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0]) + b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])) / det;
    }else if(n == 4){
        double det = A[0][0] * (A[1][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) + A[1][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])) - A[0][1] * (A[1][0] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])) + A[0][2] * (A[1][0] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) - A[1][1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])) - A[0][3] * (A[1][0] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]) - A[1][1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0]) + A[1][2] * (A[2][0] * A[3][1] - A[2][1] * A[3][0]));
        if(det == 0){
            cout << "Determinant is zero, cannot solve the equations" << endl;
            return false;
        }

        x[0] = (b[0] * (A[1][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) + A[1][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])) - A[0][1] * (b[1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (b[2] * A[3][3] - A[2][3] * b[3]) + A[1][3] * (b[2] * A[3][2] - A[2][2] * b[3])) + A[0][2] * (b[1] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) - A[1][1] * (b[2] * A[3][3] - A[2][3] * b[3]) + A[1][3] * (b[2] * A[3][1] - A[2][1] * b[3])) - A[0][3] * (b[1] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]) - A[1][1] * (b[2] * A[3][2] - A[2][2] * b[3]) + A[1][2] * (b[2] * A[3][1] - A[2][1] * b[3]))) / det;
        x[1] = (A[0][0] * (b[1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (b[2] * A[3][3] - A[2][3] * b[3]) + A[1][3] * (b[2] * A[3][2] - A[2][2] * b[3])) - b[0] * (A[1][0] * (A[2][2] * A[3][3] - A[2][3] * A[3][2]) - A[1][2] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])) + A[0][2] * (A[1][0] * (b[2] * A[3][3] - A[2][3] * b[3]) - b[1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * b[3] - b[2] * A[3][0])) - A[0][3] * (A[1][0] * (b[2] * A[3][2] - A[2][2] * b[3]) - A[1][2] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0]))) / det;
        x[2] = (A[0][0] * (A[1][1] * (b[2] * A[3][3] - A[2][3] * b[3]) - b[1] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) + A[1][3] * (A[2][1] * b[3] - b[2] * A[3][1])) - A[0][1] * (A[1][0] * (b[2] * A[3][3] - A[2][3] * b[3]) - b[1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * b[3] - b[2] * A[3][0])) + b[0] * (A[1][0] * (A[2][1] * A[3][3] - A[2][3] * A[3][1]) - A[1][1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0]) + A[1][3] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])) - A[0][3] * (A[1][0] * (A[2][1] * b[3] - b[2] * A[3][1]) - A[1][1] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][1] - A[2][1] * A[3][0]))) / det;
        x[3] = (A[0][0] * (A[1][1] * (A[2][2] * b[3] - b[2] * A[3][2]) - A[1][2] * (A[2][1] * b[3] - b[2] * A[3][1]) + b[1] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])) - A[0][1] * (A[1][0] * (A[2][2] * b[3] - b[2] * A[3][2]) - A[1][2] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])) + A[0][2] * (A[1][0] * (A[2][1] * b[3] - b[2] * A[3][1]) - A[1][1] * (A[2][0] * b[3] - b[2] * A[3][0]) + b[1] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])) - b[0] * (A[1][0] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]) - A[1][1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0]) + A[1][2] * (A[2][0] * A[3][1] - A[2][1] * A[3][0]))) / det;
    }
    return true;
};

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
    int ans;

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
            ans = Gauss_Elimination(A, b, x);
            if(ans == -1){
                printText("No solution exists.\n", 0, 2, true);
            }
            else if(ans == 1){
                printText("Infinitely many solutions exist.\n", 0, 2, true);
            }
            else{
                display_answers(x);
            }
            break;
        case 4:
            inputMatrix(A, b, n);
            ans = Gauss_Jordan_Elimination(A, b);
            if(ans == -1){
                printText("No solution exists.\n", 0, 2, true);
            }
            else if(ans == 1){
                printText("Infinitely many solutions exist.\n", 0, 2, true);
            }
            else{
                display_answers(b);
            }
            break;
        case 5:
            inputMatrix(A, b, n);
            if(LU_Factorization(A, b, x)){
                display_answers(x);
            }
            break;
    }

    cout << endl << "Press any key to continue..." << endl;
    getChar();
}