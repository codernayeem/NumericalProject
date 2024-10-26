#include "non_linear_eqs.hpp"
#include "utils.hpp"
#include <iostream>
#include <math.h>
#include <iomanip>
#include <limits>
const double tolerance = 0.0000001;
using namespace std;

void bisection(vector<double> coef,double a,double b)
{
    cout << "[+] - Search interval for root = [" << a << ", " << b << "]" << endl;
    if (isnan(a) || isnan(b))
    {
        cout << "[!] - No Solution" << endl;
        return;
    }
    if (f(coef, a) * f(coef, b) >= 0)
    {
        cout << "[!] - The function does not change sign on the interval" << endl;
        return;
    }
    double old_x=a;
    int itr = 0;
    while (true)
    {
        itr++;
        double x = (a + b) / 2;
        double f_x = f(coef, x);
        if (fabs(x-old_x) <= tolerance|| f_x==0.0)
        {
            setColor(2, 0);
            cout <<fixed<<setprecision(4)<< "[+] - Root : " << x << endl;
            setColor(7, 0);
            cout << "[+] - itr Required : " << itr << endl;
            return;
        }

        if (f(coef, a) * f_x < 0)
        {
            b = x;
        }
        else if (f(coef, b) * f_x < 0)
        {
            a = x;
        }
        old_x=x;
    }
}

void false_position(vector<double> coef,double a,double b)
{
    cout << "[+] - Search interval for root = [" << a << ", " << b << "]" << endl;
    if (f(coef, a) * f(coef, b) >= 0)
    {
        cout << "[!] - The function does not change sign on the interval" << endl;
        return;
    }
    double old_x=a;
    int itr = 0;
    while (true)
    {
        itr++;
        double x = (a * f(coef, b) - b * f(coef, a)) / (f(coef, b) - f(coef, a));
        double f_x = f(coef, x);
        if (fabs(x-old_x) <= tolerance|| f_x==0.0)
        {
            setColor(2, 0);
            cout <<fixed<<setprecision(4)<< "[+] - Root : " << x << endl;
            setColor(7, 0);
            cout << "[+] - itr Required : " << itr << endl;
            return;
        }

        if (f(coef, a) * f_x < 0)
        {
            b = x;
        }
        else if (f(coef, b) * f_x < 0)
        {
            a = x;
        }
        old_x=x;
    }

}

// Synthetic division to deflate the polynomial
void syntheticDivision(vector<double>& coef, double root) {
    vector<double> new_coef(coef.size() - 1);
    new_coef[0] = coef[0];
    for (int i = 1; i < new_coef.size(); ++i) {
        new_coef[i] = coef[i] + root * new_coef[i - 1];
    }
    coef = new_coef;
}

// Newton-Raphson method to find a root
double newton_raphson(const vector<double>& coef, double initial_guess, int maxIter = 500, double tolerance = 1e-10) {
    double x = initial_guess, i;

    for (i = 0; i < maxIter; ++i) {
        double fx = f(coef, x);
        double fpx = fprime(coef, x);

        if (fabs(fpx) < tolerance){
            cout << "[!] - Division by zero. Cannot Find Solution" << endl;
            return numeric_limits<double>::infinity();
        };

        double x_new = x - fx / fpx;

        if (fabs(f(coef, x_new)) < tolerance){
            setColor(2, 0);
            cout << fixed << setprecision(4) << "[+] - Root : " << x_new << endl;
            setColor(7, 0);
            cout << "[+] - itr Required : " << i+1 << endl;
            return x_new;
        }
        x = x_new;
    }
    
    cout << "[!] - Maximum iterations reached. Cannot Find Solution" << endl;
    return numeric_limits<double>::infinity();
}

void findRoots_newton_raphson(vector<double> coef) {
    double x;
    cout<<"Enter initial guess x0: ";
    cin>>x;
    while (coef.size() > 1) {
        double r = newton_raphson(coef, x);
        if (r == numeric_limits<double>::infinity()) {
            return;
        }
        syntheticDivision(coef, r);
    }
}

double secant(vector<double> coef, double x0, double x1, double maxIter = 2000, double tolerance = 1e-7)
{
    int itr = 0;
    while (itr < maxIter)
    {
        itr++;
        double x = x1 - ((f(coef,x1) * (x1 - x0)) / (f(coef,x1) - f(coef,x0)));
        
        if (fabs(f(coef, x)) <= tolerance)
        {
            setColor(2, 0);
            cout  << "[+] - Root : " << x << endl;
            setColor(7, 0);
            cout << "[+] - itr Required : " << itr << endl;
            return x;
        }
        x0 = x1;
        x1 = x;
    }
    cout << "[!] - Maximum iterations reached. Cannot Find Solution" << endl;
    return numeric_limits<double>::infinity();
}

void findRoots_secant(vector<double> coef) {
    double x0,x1;
    cout<<"Enter initial guess x0 & x1: ";
    cin>>x0>>x1;
    while (coef.size() > 1) {
        double r = secant(coef, x0, x1);
        if (r == numeric_limits<double>::infinity()) {
            return;
        }
        syntheticDivision(coef, r);
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
    vector<double> coef;
    coef = inputPolynomial();
    cout << "[+] - Equation: ";
    printPolynomial(coef);
    cout << endl;
    if (choice == 1){
        double xMax = sqrt((coef[1] / coef[0]) * (coef[1] / coef[0]) - 2 * (coef[2] / coef[0]));
        double a,b;
        int root=1;
        for(double j=-xMax; j<=xMax; j+=0.1)
        {
            if(f(coef,j)*f(coef,j+0.1)<0)
            {
                a=j;
                b=j+0.1;
                cout<<"For root "<<root++<<endl;
                bisection(coef,a,b);
            }
        }   
    }
    else if (choice == 2){
        double xMax = sqrt((coef[1] / coef[0]) * (coef[1] / coef[0]) - 2 * (coef[2] / coef[0]));
        double a,b;
        int root =1;
        for(double j=-xMax; j<=xMax; j+=0.1)
        { 
            if(f(coef,j)*f(coef,j+0.1)<0)
            {
                a=j;
                b=j+0.1;
                cout<<"For root "<<root++<<endl;
                false_position(coef,a,b);
            }
        } 
    }
    else if (choice == 3){
        findRoots_newton_raphson(coef);
    }
    else if (choice == 4){
        findRoots_secant(coef);
    }

    cout << endl << "Press any key to continue..." << endl;
    getChar();
}