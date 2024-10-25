#include "utils.hpp"
#include <iostream>
#include <vector>
#include <conio.h>
#include <windows.h>
#include <chrono>
#include <thread>
#include <vector>
#include <cstdlib>
const int WIDTH = 60;
const int HEIGHT = 20;
const int DURATION = 15; // Adjust total frames for the animation duration
const int TEXT_DELAY = 10; // Delay before showing the welcome message

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

void displayMatrixRain() {
    vector<int> lengths(WIDTH, 0);
    vector<int> positions(WIDTH, 0);
    vector<int> speeds(WIDTH, 0);

    // Initialize columns with random properties for the rain effect
    for (int i = 0; i < WIDTH; ++i) {
        lengths[i] = rand() % HEIGHT;
        positions[i] = -rand() % HEIGHT;
        speeds[i] = 1 + rand() % 3;
    }

    for (int frame = 0; frame < DURATION; ++frame) {
        cout << "\033[2J\033[1;1H"; // Clear the console

        for (int y = 0; y < HEIGHT; ++y) {
            for (int x = 0; x < WIDTH; ++x) {
                if (positions[x] == y) {
                    // Display binary characters ('0' or '1') in green
                    cout << "\033[32m" << (rand() % 2 == 0 ? '0' : '1');
                } else if (positions[x] - lengths[x] <= y && y < positions[x]) {
                    // Trailing effect with light green dots
                    cout << "\033[92m.";
                } else {
                    cout << " ";
                }
            }
            cout << "\n";
        }

        // Update column positions to create the falling effect
        for (int x = 0; x < WIDTH; ++x) {
            positions[x] += speeds[x];
            if (positions[x] - lengths[x] >= HEIGHT) {
                positions[x] = -rand() % HEIGHT; // Reset column
                lengths[x] = rand() % HEIGHT;    // Reset length
            }
        }

        // Display the welcome message after a brief delay
        if (frame >= TEXT_DELAY) {
            string message = "Welcome to Our Application";
            int startX = (WIDTH - message.length()) / 2;
            int startY = HEIGHT / 2;

            // Move cursor to center and display the welcome message in pink
            cout << "\033[" << startY << ";" << startX << "H\033[35m" << message;
        }

        cout << "\033[0m"; // Reset color
        this_thread::sleep_for(chrono::milliseconds(100)); // Control frame speed
    }
    cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\";
}


   
