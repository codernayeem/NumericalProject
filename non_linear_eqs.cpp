#include "non_linear_eqs.hpp"
#include "utils.hpp"
#include <bits/stdc++.h>
#include <math.h>
const double tolerance = 0.0000001;
using namespace std;

void bisection(vector<double> coef,double a,double b)
{
    cout << "[+] - Search interval for root = [" << a << ", " << b << "]" << endl;
    if (isnan(a) || isnan(b))
    {
        cerr << "[!] - No Solution" << endl;
        return;
    }
    if (f(coef, a) * f(coef, b) >= 0)
    {
        cerr << "[!] - The function does not change sign on the interval" << endl;
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
            cout <<fixed<<setprecision(4)<< "[+] - Root : " << x << endl;
            cout << "[+] - Iteration Required : " << itr << endl;
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
        cerr << "[!] - The function does not change sign on the interval" << endl;
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
            cout <<fixed<<setprecision(4)<< "[+] - Root : " << x << endl;
            cout << "[+] - Iteration Required : " << itr << endl;
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

void newton_raphson(vector<double> coef)
{
     int x;
    cout<<"Enter initial guess x0: ";
    cin>>x;
    int itr = 0;
    while (true)
    {
        itr++;
        double f_x = f(coef, x);
        double f_prime_x = fprime(coef, x);
        // Newton-Raphson update formula
        double x_new = x - f_x / f_prime_x;
         if (fabs(x_new) < tolerance)
            x = 0.0;
        if (fabs(x_new - x) <= tolerance ||f(coef,x_new) == 0.0)
        {
            cout  <<fixed<<setprecision(4)<< "[+] - Root : " << x_new << endl;
            cout << "[+] - Iteration Required : " << itr << endl;
            return;
        }
        x = x_new;
    }
}

void secant(vector<double> coef)
{
    int x0,x1;
    cout<<"Enter initial guess x0 & x1: ";
    cin>>x0>>x1;
    int itr = 0;
    while (true)
    {
        itr++;
        double f_x0 = f(coef, x0);
        double f_x1 = f(coef, x1);
        double x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
        if (fabs(x2) < tolerance)
            x2 = 0.0;
        if (fabs(x2 - x1) <= tolerance || fabs(f(coef, x2)) <= tolerance)
        {
            cout <<fixed<<setprecision(4)<< "[+] - Root : " << x2 << endl;
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
       
        newton_raphson(coef);
    }
    else if (choice == 4){
      
        secant(coef);
    }
    cout << endl
         << "Press any key to continue..." << endl;
    getChar();
}