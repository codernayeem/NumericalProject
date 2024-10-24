#include <iostream>
#include "linear_eqs.hpp"
#include "non_linear_eqs.hpp"
#include "diff_eqs.hpp"
#include "matrix_inversion.hpp"

using namespace std;

void displayMenu() {
    cout << "Numerical Methods Console Application\n";
    cout << "1. Solve Linear Equations\n";
    cout << "2. Solve Non-Linear Equations\n";
    cout << "3. Solve Differential Equations\n";
    cout << "4. Matrix Inversion\n";
    cout << "5. Exit\n";
}

int main() {
    int choice;
    do {
        displayMenu();
        cout << "Enter your choice: ";
        cin >> choice;
        
        switch(choice) {
            case 1:
                solveLinearEquations();
                break;
            case 2:
                solveNonLinearEquations();
                break;
            case 3:
                solveDifferentialEquations();
                break;
            case 4:
                invertMatrix();
                break;
            case 5:
                cout << "Exiting...\n";
                break;
            default:
                cout << "Invalid choice. Try again.\n";
        }
    } while (choice != 5);

    return 0;
}
