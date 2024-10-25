#include <iostream>
#include "linear_eqs.hpp"
#include "non_linear_eqs.hpp"
#include "diff_eqs.hpp"
#include "matrix_inversion.hpp"
#include "utils.hpp"

using namespace std;

void displayMenu() {
    clearScreen();
    cout << "\n";
    cout << "\t ==============================" << endl;
    cout << "\t|                              |" << endl;
    cout << "\t    ";
    printText("Numerical Methods Project", 0, 2, true);
    cout << endl;
    cout << "\t|                              |" << endl;
    cout << "\t ============================== " << endl;

    cout << "\n";
    cout << "\t  1. Solve Linear Equations" << endl;
    cout << "\t  2. Solve Non-Linear Equations" << endl;
    cout << "\t  3. Solve Differential Equations" << endl;
    cout << "\t  4. Matrix Inversion" << endl;
    cout << "\t  5. Quit" << endl;
    cout << "\n";

    printText("Press a key (1-5) to select an Option", 0, 2);

    cout << "\n";
}

void exitScreen() {
    cout << "\n\n";
    cout << "\t ==============================" << endl;
    cout << "\t|                              |" << endl;
    cout << "\t      ";
    printText("Thank you ... (-_-)", 0, 2);
    cout << endl;
    cout << "\t|                              |" << endl;
    cout << "\t ============================== " << endl;
    cout << "\n";
    exit(0);
}


int main() {
    char choice;

    while (true) {
        displayMenu();
        choice = getChar();
        clearScreen();
        
        switch(choice) {
            case '1':
                solveLinearEquations();
                break;
            case '2':
                solveNonLinearEquations();
                break;
            case '3':
                solveDifferentialEquations();
                break;
            case '4':
                invertMatrix();
                break;
            case '5':
                exitScreen();
                break;
        }
    }

    return 0;
}
