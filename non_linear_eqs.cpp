#include "non_linear_eqs.hpp"
#include <bits/stdc++.h>
using namespace std;
const double max_tolerence = 0.00001;

double a=1,b=-3,c=0,d=2;
double f(double x)
{
    return x * x - 3 * x + 2; // f(x) = x^2 - 3x + 2
}
double f_prime(double x)
{
    return 2 * x - 3; // f'(x) = 2x - 3
}
void bisection(double a, double b)
{
    cout << "Using Bisection Method:" << endl;
    if (f(a) * f(b) >= 0)
    {
        cerr << "The function does not change sign on the interval." << endl;
        return;
    }
    int iteration = 0;
    int flag = 0;
    double old_x = a;
    while (true)
    {
        iteration++;
        double x = (a + b) / 2;
        // cout << "a = " << a << " " << "b = " << b << " ";
        // cout << "Approximate Root = " << x << " " << f(x);
        // if (flag)
        // {
        //     cout << " error = " << x - old_x;
        // }
        cout << endl;
        if (f(x) == 0.0)
        {
            cout << " Root = " << x << endl;
           // cout << "Iteration = " << iteration << endl;
            return;
        }
        else if (f(x) * f(a) < 0)
            b = x;
        else if (f(x) * f(b) < 0)
            a = x;
        if (fabs(x - old_x) <= max_tolerence)
        {
            cout << " Root = " << x << endl;
            //cout << "Iteration = " << iteration << endl;
            return;
        }
        old_x = x;
        flag = 1;
    }
}

void false_position(double a, double b)
{
    cout << "Using False Position Method:" << endl;
    if (f(a) * f(b) >= 0)
    {
        cerr << "The function does not change sign on the interval." << endl;
        return;
    }
    int iteration = 0, flag = 0;
    double old_x = a;
    while (true)
    {
        iteration++;
        double x = (a * f(b) - b * f(a)) / (f(b) - f(a));
        // cout << "a = " << a << " " << "b = " << b << " ";
        // cout << "Approximate Root = " << x<<" "<<f(x);
        // if (flag)
        // {
        //     cout << " error = " << x - old_x;
        // }
        cout << endl;
        if (f(x) == 0.0)
        {
            cout << " Root = " << x << endl;
           // cout << "Iteration = " << iteration << endl;
            return;
        }
        else if (f(x) * f(a) < 0)
            b = x;
        else if (f(x) * f(b) < 0)
            a = x;
        if (fabs(x - old_x) <= max_tolerence)
        {
            cout << " Root = " << x << endl;
            //cout << "Iteration = " << iteration << endl;
            return;
        }
        old_x = x;
        flag = 1;
    }
}

void newton_raphson()
{
    cout << "Using Newton Raphson Method:" << endl;
    int iteration = 0;
    double x = 1.5; // Initial guess
    double old_x = 0;

    while (true)
    {
        iteration++;
        old_x = x;

        // Newton-Raphson update formula
        x = old_x - (f(old_x) / f_prime(old_x));

        // If x is close enough to zero, set it to zero for precision
        if (fabs(x) < max_tolerence)
        {
            x = 0.0;
        }

        // Print the approximation and f(x)
        // cout << "Approximate Root = " << x << " " << f(x) << std::endl;

        // Check for convergence to root
        if (fabs(f(x)) < max_tolerence || fabs(x - old_x) <= max_tolerence)
        {
            cout << "Root = " << x << std::endl;
          //  cout << "Iterations = " << iteration << std::endl;
            return;
        }
    }
}

void secant()
{
    cout << "Using Secant Method:" << endl;
    int iteration = 0;
    double x1 = 0.5, x2 = 1.5, old_x = 0;
    while (true)
    {
        iteration++;
        double x = x2 - ((f(x2) * (x2 - x1)) / (f(x2) - f(x1)));
        if (fabs(x) < max_tolerence)
            x = 0.0;
        // cout << "Approximate Root = " << x << " " << f(x) << endl;
        if (f(x) == 0.0)
        {
            cout << " Root = " << x << endl;
         //   cout << "Iteration = " << iteration << endl;
            return;
        }
        if (fabs(x - x2) <= max_tolerence)
        {
            cout << "Root = " << x << endl;
        //  cout << "Iteration = " << iteration << endl;
            return;
        }
        x1 = x2;
        old_x = x2 = x;
    }
}
void solveNonLinearEquations()
{
    cout<<"Solving Non-Linear Equations"<<endl;
    cout<<"1. Bi-section Method"<<endl;
    cout<<"2. False Position Method"<<endl;
    cout<<"3. Newton-Raphson Method"<<endl;
    cout<<"4. Secant Method"<<endl;
    cout<<"Which method do you want to use?"<<endl;
    int choice;
    cin>>choice;

    if(choice==1) {
        double xMax=abs(sqrt(((b/a)*(b/a))-(2*(c/a))));
        for(int i=-xMax; i<xMax; i++)
        {
            if (f(i) * f(i + 1) < 0)
            {
                a = i;
                b = i + 1;
                break;
            }
        }
    bisection(a,b);
    }
    else if(choice==2){
        double xMax=abs(sqrt(((b/a)*(b/a))-(2*(c/a))));
        for(int i=-xMax; i<xMax; i++)
        {
            if (f(i) * f(i + 1) < 0)
            {
                a = i;
                b = i + 1;
                break;
            }
        }
        false_position(a,b);
    }
    else if(choice==3)newton_raphson();
    else if(choice==4)secant();
    return;
}
