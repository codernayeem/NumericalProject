#include "diff_eqs.hpp"
#include "utils.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <functional>
#include <sstream>

using namespace std;

double aSinbx(double a, double b, double x) {
    return a * sin(b * x);
}

double aCosbx(double a, double b, double x) {
    return a * cos(b * x);
}

double aTanbx(double a, double b, double x) {
    return a * tan(b * x);
}

double aLogbx(double a, double b, double x) {
    return a * log(b * x);
}

double aX2_bXY_cY2(double a, double b, double c, double x, double y) {
    return a * x * x + b * x * y + c * y * y;
}

double aX_bY_c(double a, double b, double c, double x, double y) {
    return a * x + b * y + c;
}

string formatEquation(double a, double b, const string& trigFunc) {
    ostringstream oss;
    if (a == 1) {
        oss << trigFunc << "(" << b << "x)";
    } else if (a == -1) {
        oss << "-" << trigFunc << "(" << b << "x)";
    } else if (a != 0) {
        oss << a << trigFunc << "(" << b << "x)";
    }
    return oss.str();
}

string formatEquation(double a, double b, double c, const string& term1, const string& term2, const string& term3) {
    ostringstream oss;
    bool firstTerm = true;

    if (a != 0) {
        if (a == 1) {
            oss << term1;
        } else if (a == -1) {
            oss << "-" << term1;
        } else {
            oss << a << term1;
        }
        firstTerm = false;
    }

    if (b != 0) {
        if (!firstTerm) oss << " + ";
        if (b == 1) {
            oss << term2;
        } else if (b == -1) {
            oss << "-" << term2;
        } else {
            oss << b << term2;
        }
        firstTerm = false;
    }

    if (c != 0) {
        if (!firstTerm) oss << " + ";
        if (c == 1) {
            oss << term3;
        } else if (c == -1) {
            oss << "-" << term3;
        } else {
            oss << c << term3;
        }
    }

    return oss.str();
}

void rk_method(function<double(double, double)> func, double x0, double y0, double h, double range) {
    double x = x0, y = y0;

    cout << "+---------------+---------------+" << endl;
    cout << "|       x       |    y = f(x)   |" << endl;
    cout << "+---------------+---------------+" << endl;
    cout << "| " << setw(13) << x << " | " << setw(13) << y << " |" << endl;

    
    const double epsilon = 1e-9;
    while ((h > 0 && x < range - epsilon) || (h < 0 && x > range + epsilon)) {
        double k1 = h * func(x, y);
        double k2 = h * func(x + h / 2, y + k1 / 2);
        double k3 = h * func(x + h / 2, y + k2 / 2);
        double k4 = h * func(x + h, y + k3);

        y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        x = x + h;

        cout << "| " << setw(13) << x << " | " << setw(13) << y << " |" << endl;
    }
    cout << "+---------------+---------------+" << endl;
}

void displayDifferentialEquationsMenu() {
    clearScreen();
    cout << "\n";
    cout << "\t ==============================" << endl;
    cout << "\t|                              |" << endl;
    cout << "\t  ";
    printText("Solve Differential Equations", 0, 2, true);
    cout << endl;
    cout << "\t      ";
    printText("> Runge-Kutta Method <", 0, 2, true);
    cout << endl;
    cout << "\t|                              |" << endl;
    cout << "\t ============================== " << endl << endl;
}

void solveDifferentialEquations() {
    displayDifferentialEquationsMenu();
    double x0, y0, h, range;

    int choice;
    cout << "\t>> Choose the function to be used:\n\n";
    cout << "\t1. Polynomial\n";
    cout << "\t2. aSin(bx)\n";
    cout << "\t3. aCos(bx)\n";
    cout << "\t4. aTan(bx)\n";
    cout << "\t5. aLog(bx)\n";
    cout << "\t6. aX^2 + bXY + cY^2\n";
    cout << "\t7. aX + bY + c\n";
    cout << "\t8. Go Back to Main Menu\n";
    cout << "\n";
    printText("\t<< Enter your choice (1-8) >>", 0, 2);
    cout << endl;

    while(true) {
        choice = getChar();
        if (choice >= '1' && choice <= '8') {
            choice -= '0';
            break;
        }
    }

    if(choice == 8) {
        return;
    }

    switch (choice) {
        case 1:
            cout << "[+] - Polynomial selected." << endl;
            break;
        case 2:
            cout << "[+] - aSin(bx) selected." << endl;
            break;
        case 3:
            cout << "[+] - aCos(bx) selected." << endl;
            break;
        case 4:
            cout << "[+] - aTan(bx) selected." << endl;
            break;
        case 5:
            cout << "[+] - aLog(bx) selected." << endl;
            break;
        case 6:
            cout << "[+] - aX^2 + bXY + cY^2 selected." << endl;
            break;
        case 7:
            cout << "[+] - aX + bY + c selected." << endl;
            break;
    }

    double a, b, c;
    function<double(double, double)> func;
    vector<double> coef;

    switch (choice) {
        case 1:
            coef = inputPolynomial();
            cout << "[+] - Equation: ";
            printPolynomial(coef);
            func = [coef](double x, double y) { return f(coef, x); };
            break;
        case 2:
            cout << "[?] - Enter a: ";
            cin >> a;
            cout << "[?] - Enter b: ";
            cin >> b;
            cout << "[+] - Equation: " << formatEquation(a, b, "sin") << endl;
            func = [a, b](double x, double y) { return aSinbx(a, b, x); };
            break;
        case 3:
            cout << "[?] - Enter a: ";
            cin >> a;
            cout << "[?] - Enter b: ";
            cin >> b;
            cout << "[+] - Equation: " << formatEquation(a, b, "cos") << endl;
            func = [a, b](double x, double y) { return aCosbx(a, b, x); };
            break;
        case 4:
            cout << "[?] - Enter a: ";
            cin >> a;
            cout << "[?] - Enter b: ";
            cin >> b;
            cout << "[+] - Equation: " << formatEquation(a, b, "tan") << endl;
            func = [a, b](double x, double y) { return aTanbx(a, b, x); };
            break;
        case 5:
            cout << "[?] - Enter a: ";
            cin >> a;
            cout << "[?] - Enter b: ";
            cin >> b;
            cout << "[+] - Equation: " << formatEquation(a, b, "log") << endl;
            func = [a, b](double x, double y) { return aLogbx(a, b, x); };
            break;
        case 6:
            cout << "[?] - Enter a: ";
            cin >> a;
            cout << "[?] - Enter b: ";
            cin >> b;
            cout << "[?] - Enter c: ";
            cin >> c;
            cout << "[+] - Equation: " << formatEquation(a, b, c, "x^2", "xy", "y^2") << endl;
            func = [a, b, c](double x, double y) { return aX2_bXY_cY2(a, b, c, x, y); };
            break;
        case 7:
            cout << "[?] - Enter a: ";
            cin >> a;
            cout << "[?] - Enter b: ";
            cin >> b;
            cout << "[?] - Enter c: ";
            cin >> c;
            cout << "[+] - Equation: " << formatEquation(a, b, c, "x", "y", "") << endl;
            func = [a, b, c](double x, double y) { return aX_bY_c(a, b, c, x, y); };
            break;
    }

    cout << endl;
    cout << "[?] - Enter initial value of x (x0): ";
    cin >> x0;
    cout << "[?] - Enter initial value of y (y0): ";
    cin >> y0;
    do {
        cout << "[?] - Enter step size (h) (h != 0) : ";
        cin >> h;
        if (h == 0) {
            cout << "[!] - Step size must not be 0. Please try again." << endl;
        }
    } while (h == 0);
    cout << "[?] - Enter the range ( " << (h > 0 ? ">" : "<") << " " << x0 << ")       : ";
    cin >> range;

    cout << endl;
    
    rk_method(func, x0, y0, h, range);

    cout << endl << "Press any key to continue..." << endl;
    getChar();
}
