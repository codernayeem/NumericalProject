@echo off

:: Step 1: Compile each .cpp file into object files
g++ -c main.cpp linear_eqs.cpp non_linear_eqs.cpp diff_eqs.cpp matrix_inversion.cpp utils.cpp -std=c++11

:: Step 2: Link all the object files into a single executable
g++ -o numerical_methods_app main.o linear_eqs.o non_linear_eqs.o diff_eqs.o matrix_inversion.o utils.o

:: Step 3: Run the application
numerical_methods_app.exe
