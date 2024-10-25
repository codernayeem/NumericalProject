#include "non_linear_eqs.hpp"
#include "utils.hpp"
#include <bits/stdc++.h>
#include <math.h>
const double tolerance = 0.00001;
using namespace std;

void bisection(vector<double> coefficient_matrix)
{
    double xMax = sqrt((coefficient_matrix[1] / coefficient_matrix[0]) * (coefficient_matrix[1] / coefficient_matrix[0]) - 2 * (coefficient_matrix[2] / coefficient_matrix[0]));
    double a = -xMax;
    double b = xMax;
    cout << "[+] - Search interval for root = [" << a << ", " << b << "]" << endl;
    if (isnan(a) || isnan(b))
    {
        cerr << "[!] - No Solution" << endl;
        return;
    }
    if (f(coefficient_matrix, a) * f(coefficient_matrix, b) >= 0)
    {
        cerr << "[!] - The function does not change sign on the interval" << endl;
        return;
    }
    int itr = 0;
    while (true)
    {
        itr++;
        double x = (a + b) / 2;
        double f_x = f(coefficient_matrix, x);

        if (fabs(f_x) <= tolerance)
        {
            cout << "[+] - Root : " << x << endl;
            cout << "[+] - Iteration Required : " << itr << endl;
            return;
        }

        if (f(coefficient_matrix, a) * f_x < 0)
        {
            b = x;
        }
        else
        {
            a = x;
        }
    }
}

void false_position(vector<double> coefficient_matrix, double tolerance = 0.0001)
{
    double xMax = sqrt((coefficient_matrix[1] / coefficient_matrix[0]) * (coefficient_matrix[1] / coefficient_matrix[0]) - 2 * (coefficient_matrix[2] / coefficient_matrix[0]));
    double a = -xMax;
    double b = xMax;
    cout << "[+] - Search interval for root = [" << a << ", " << b << "]" << endl;
    if (f(coefficient_matrix, a) * f(coefficient_matrix, b) >= 0)
    {
        cerr << "[!] - The function does not change sign on the interval" << endl;
        return;
    }
    int itr = 0;
    while (true)
    {
        itr++;
        double x = (a * f(coefficient_matrix, b) - b * f(coefficient_matrix, a)) / (f(coefficient_matrix, b) - f(coefficient_matrix, a));
        double f_x = f(coefficient_matrix, x);

        if (fabs(f_x) <= tolerance)
        {
            cout << "[+] - Root : " << x << endl;
            cout << "[+] - Iteration Required : " << itr << endl;
            return;
        }

        if (f(coefficient_matrix, a) * f_x < 0)
        {
            b = x;
        }
        else
        {
            a = x;
        }
    }
}

void newton_raphson(vector<double> coefficient_matrix, double tolerance = 0.0001)
{
    int itr = 0;
    double x = 0.0;
    while (true)
    {
        itr++;
        double f_x = f(coefficient_matrix, x);
        double f_prime_x = fprime(coefficient_matrix, x);
        // Newton-Raphson update formula
        double x_new = x - f_x / f_prime_x;

        if (fabs(x_new - x) <= tolerance || fabs(f_x) <= tolerance)
        {
            cout << "[+] - Root : " << x_new << endl;
            cout << "[+] - Iteration Required : " << itr << endl;
            return;
        }
        x = x_new;
    }
}

void secant(vector<double> coefficient_matrix, double tolerance = 0.0001)
{
    double x0 = 0.0;
    double x1 = 1.0;
    int itr = 0;
    while (true)
    {
        itr++;
        double f_x0 = f(coefficient_matrix, x0);
        double f_x1 = f(coefficient_matrix, x1);
        double x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
        if (fabs(x2 - x1) <= tolerance || fabs(f_x1) <= tolerance)
        {
            cout << "[+] - Root : " << x2 << endl;
            cout << "[+] - Iteration Required : " << itr << endl;
            return;
        }
        x0 = x1;
        x1 = x2;
    }
}
void displayNonLinearEquationsMenu()
{
    clearScreen();
    cout << "\n";
    cout << "\t ==============================" << endl;
    cout << "\t|                              |" << endl;
    cout << "\t  ";
    printText("\033[35m Solve Non-Linear Equations \033[0m", 0, 2, true);
    cout << endl;
    cout << "\t|                              |" << endl;
    cout << "\t ============================== " << endl
         << endl;
}

void solveNonLinearEquations()
{
    displayNonLinearEquationsMenu();
    cout << "\t>> Choose a option:\n\n";
    cout << "\t1. Bi-section Method\n";
    cout << "\t2. False Position Method\n";
    cout << "\t3. Newton-Raphson Method\n";
    cout << "\t4. Secant Method\n";
    cout << "\t5. Go Back to Main Menu\n";
    cout << "\n";
    printText("\t\033[35m<< Enter your choice (1-5) \033>>>", 0, 2);
    cout << endl;
    int choice;
    while (true)
    {
        choice = getChar();
        if (choice >= '1' && choice <= '8')
        {
            choice -= '0';
            break;
        }
    }
    if (choice == 5)
    {
        return;
    }
    cout << endl;
    if (choice == 1)
        cout << "[+] - Bi-section Method selected." << endl;
    else if (choice == 2)
        cout << "[+] - False Position selected." << endl;
    else if (choice == 3)
        cout << "[+] - Newton-Raphson Method selected." << endl;
    else if (choice == 4)
        cout << "[+] - Secant Method selected." << endl;
    vector<double> coefficient_matrix;
    coefficient_matrix = inputPolynomial();
    cout << "[+] - Equation: ";
    printPolynomial(coefficient_matrix);
    cout << endl;
    if (choice == 1)
        bisection(coefficient_matrix);
    else if (choice == 2)
        false_position(coefficient_matrix);
    else if (choice == 3)
        newton_raphson(coefficient_matrix);
    else if (choice == 4)
        secant(coefficient_matrix);
    cout << endl
         << "Press any key to continue..." << endl;
    getChar();
}