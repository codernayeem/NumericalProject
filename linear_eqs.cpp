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

void makeDiagonallyDominant(vector<vector<double>>& A, vector<double>& B) {
    int n = A.size();
    
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        double maxElement = fabs(A[i][i]);

        for (int j = i + 1; j < n; ++j) {
            if (fabs(A[j][i]) > maxElement) {
                maxElement = fabs(A[j][i]);
                maxRow = j;
            }
        }

        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(B[i], B[maxRow]);
        }

        double sum = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                sum += fabs(A[i][j]);
            }
        }

        if (fabs(A[i][i]) < sum) {
            return ;
        }
    }
    return;
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

    // for max 4x4 matrix
    double L21, L31, L41, L32, L42, L43;
    double U11, U12, U13, U14, U22, U23, U24, U33, U34, U44;
    double y1, y2, y3, y4;
    double x1, x2, x3, x4;
    if (n == 1){
        x[0] = b[0] / A[0][0];
    }if(n == 2){
        L21 = A[1][0] / A[0][0];
        U11 = A[0][0];
        U12 = A[0][1];
        U22 = A[1][1] - L21 * U12;

        y1 = b[0];
        y2 = b[1] - L21 * y1;

        x2 = y2 / U22;
        x1 = (y1 - U12 * x2) / U11;

        x[0] = x1, x[1] = x2;
    }else if(n == 3){
        L21 = A[1][0] / A[0][0];
        L31 = A[2][0] / A[0][0];
        L32 = A[2][1] / A[1][1];

        U11 = A[0][0];
        U12 = A[0][1];
        U13 = A[0][2];
        U22 = A[1][1] - L21 * U12;
        U23 = A[1][2] - L21 * U13;
        U33 = A[2][2] - L31 * U13 - L32 * U23;

        y1 = b[0];
        y2 = b[1] - L21 * y1;
        y3 = b[2] - L31 * y1 - L32 * y2;

        x3 = y3 / U33;
        x2 = (y2 - U23 * x3) / U22;
        x1 = (y1 - U12 * x2 - U13 * x3) / U11;

        x[0] = x1, x[1] = x2, x[2] = x3;
    }else if(n == 4){
        L21 = A[1][0] / A[0][0];
        L31 = A[2][0] / A[0][0];
        L41 = A[3][0] / A[0][0];
        L32 = A[2][1] / A[1][1];
        L42 = A[3][1] / A[1][1];
        L43 = A[3][2] / A[2][2];

        U11 = A[0][0];
        U12 = A[0][1];
        U13 = A[0][2];
        U14 = A[0][3];
        U22 = A[1][1] - L21 * U12;
        U23 = A[1][2] - L21 * U13;
        U24 = A[1][3] - L21 * U14;
        U33 = A[2][2] - L31 * U13;
        U34 = A[2][3] - L31 * U14;
        U44 = A[3][3] - L41 * U14 - L42 * U24 - L43 * U34;

        y1 = b[0];
        y2 = b[1] - L21 * y1;
        y3 = b[2] - L31 * y1 - L32 * y2;
        y4 = b[3] - L41 * y1 - L42 * y2 - L43 * y3;

        x4 = y4 / U44;
        x3 = (y3 - U34 * x4) / U33;
        x2 = (y2 - U23 * x3 - U24 * x4) / U22;
        x1 = (y1 - U12 * x2 - U13 * x3 - U14 * x4) / U11;

        x[0] = x1, x[1] = x2, x[2] = x3, x[3] = x4;
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
            makeDiagonallyDominant(A, b);
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
            makeDiagonallyDominant(A, b);
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