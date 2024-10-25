#include "utils.hpp"
#include <iostream>
#include <vector>
#include <conio.h>
#include <windows.h>

using namespace std;

void setColor(int textColor, int bgColor) {
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), (bgColor << 4) | textColor);
}

void printText(string text, int bgColor = 0, int textColor = 7, bool bold = false) {
    setColor(textColor, bgColor);
    if(bold) {
        cout << "\e[32;1m";
    }
    cout << text;
    if(bold) {
        cout << "\e[0m";
    }
    setColor(7, 0);
}

void printText(string text, int bgColor = 0, int textColor = 7) {
    printText(text, bgColor, textColor, false);
}

void clearScreen() {
    system("cls");
}

char getChar() {
    return _getch();
}

double f(const vector<double>& coef, double x) {
    double result = 0.0;
    double power = 1.0;
    for (int i = coef.size() - 1; i >= 0; --i) {
        result += coef[i] * power;
        power *= x;
    }
    return result;
}

double fprime(const vector<double>& coef, double x) {
    double result = 0.0;
    double power = 1.0;
    for (int i = coef.size() - 1; i > 0; --i) {
        result += (coef.size() - i) * coef[i] * power;
        power *= x;
    }
    return result;
}

void printPolynomial(const vector<double>& coef) {
    bool first = true;
    for (int i = 0; i < coef.size(); ++i) {
        if (coef[i] != 0) {
            if (!first && coef[i] > 0) {
                cout << " + ";
            } else if (coef[i] < 0) {
                cout << " - ";
            }
            if (abs(coef[i]) != 1 || i == coef.size() - 1) {
                cout << abs(coef[i]);
            }
            if (i < coef.size() - 1) {
                cout << "x";
                if (i < coef.size() - 2) {
                    cout << "^" << (coef.size() - 1 - i);
                }
            }
            first = false;
        }
    }
    cout << endl;
}


vector<double> inputPolynomial() {
    int degree;
    cout << "[+] - Enter the degree of the polynomial: ";
    cin >> degree;
    cout << "[+] - Enter the " << degree+1 << " coefficients (from highest to lowest degree)\n  >> ";
    vector<double> coef(degree + 1);
    for (int i = 0; i <= degree; ++i) {
        cin >> coef[i];
    }
    return coef;
}
