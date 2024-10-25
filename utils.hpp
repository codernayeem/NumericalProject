#ifndef UTILS_HPP
#define UTILS_HPP
#include <string>
#include <vector>
using namespace std;

void setColor(int textColor, int bgColor);
void printText(string text, int bgColor, int textColor, bool bold);
void printText(string text, int bgColor, int textColor);
void clearScreen();
char getChar();
double f(const vector<double>& coef, double x);
double fprime(const vector<double>& coef, double x);
void printPolynomial(const vector<double>& coef);
vector<double> inputPolynomial();
void displayMatrixRain();
  

#endif
