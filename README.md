# Numerical Methods Console Application

This console application implements various numerical methods for solving linear equations, non-linear equations, and differential equations, along with matrix inversion. Each algorithm is designed to perform specific calculations and is encapsulated in separate files for modularity and clarity.

## Table of Contents
- [Usage](#usage)
- __Algorithms Explanation__
  - Linear Equations
    - Jacobi Method
    - Gauss-Seidel Method
    - Gauss Elimination
    - Gauss-Jordan Elimination
    - LU Factorization
  - Non-linear Equations
    - Bisection Method
    - False Position Method
    - Secant Method
    - Newton-Raphson Method
  - Differential Equations
    - Runge-Kutta Method
  - Matrix Inversion
- [Contributors](#contributors)

## Usage
1. Clone the repository:
   - `git clone https://github.com/codernayeem/NumericalProject.git`
2. Navigate to the project directory:
   - `cd NumericalProject`
3. Compile & Run the application:
   - `./run.bat`
4. (Optional) Run the application _(After building once)_:
   - `./numerical_methods`

# Working Principle of Algorithms and Features

### 1. Jacobi Iterative Method
The Jacobi method iteratively solves the system ( `Ax = b` ) by separating ( `A` ) into a diagonal component. 

**Steps**:
1. Initialize `x^(0)`.
2. For each iteration `k`:
   - Update each element `x_i^(k+1)` using:
     ```
     x_i^(k+1) = (b_i - sum(a_ij * x_j^(k) for j != i)) / a_ii
     ```
3. Check for convergence by evaluating if the average error falls below a tolerance level.

**Features**:
- Convergence check based on tolerance.
- Diagonal dominance check for matrix stability.

---

### 2. Gauss-Seidel Iterative Method
The Gauss-Seidel method is similar to Jacobi but uses the most recent values of ( x ) during each update.

**Steps**:
1. Initialize \( x^{(0)} \).
2. For each iteration \( k \):
   - Update \( x_i^{(k+1)} \) using previously updated values:
   ```
   x_i^{(k+1)} = (b_i - sum(a_{ij} * x_j^{(k+1)} for j < i) - sum(a_{ij} * x_j^{(k)} for j > i)) / a_{ii}
   ```
3. Convergence is achieved if the average error is less than the tolerance level.

**Features**:
- Faster convergence by using updated values immediately.

---

### 3. Gauss Elimination
Gauss Elimination transforms the matrix ( A ) into an upper triangular form to solve the system.

**Steps**:
1. Perform partial pivoting to select the largest element as the pivot.
2. For each row ( i ):
   - Eliminate variables below the pivot:
   `A[j] = A[j] - (A[j][i] / A[i][i]) * A[i]`
3. Use back substitution to solve for ( x ).

**Features**:
- Allows handling of singular matrices with a result indicating infinite or no solutions.

---

### 4. Gauss-Jordan Elimination
An extension of Gauss Elimination, Gauss-Jordan transforms \( A \) into a reduced row-echelon form.

**Steps**:
1. For each column:
   - Select a pivot, divide the row by the pivot to make it 1.
2. Eliminate other rows:
   - For each row \( i \), subtract multiples of the pivot row to zero out the column.
3. Solution vector \( x \) is obtained directly.

**Features**:
- Provides a unique solution if it exists.
- Checks for singular matrices.

---

### 5. LU Factorization
LU Factorization decomposes \( A \) as \( A = LU \), where \( L \) is a lower triangular matrix and \( U \) is an upper triangular matrix.

**Steps**:
1. Decompose \( A \) such that:
   \[
   A = LU
   \]
2. Solve \( Ly = b \) for \( y \) using forward substitution.
3. Solve \( Ux = y \) for \( x \) using back substitution.

**Features**:
- Efficient for systems with multiple right-hand side vectors.
- Limited to \( n \leq 4 \) in this implementation.

---

### Key Mathematical Terms

- **Diagonal Dominance**: A matrix \( A \) is diagonally dominant if \( |a_{ii}| > \sum_{j \neq i} |a_{ij}| \).
- **Convergence Condition**: For iterative methods, convergence is reached if the average error, calculated as:
  \[
  \text{Error} = \frac{\sum_{i=1}^n |x_i^{(k+1)} - x_i^{(k)}|}{n}
  \]
  is below the tolerance.
- **Partial Pivoting**: Selecting the largest absolute value in the column to avoid numerical instability.

The code includes tolerance settings, iterative limits, and back substitution routines for enhanced functionality.


# Working Principle of Differential Equation Solver (Runge-Kutta Method)

The application solves differential equations using the Runge-Kutta method, a widely used technique for its balance between accuracy and computational efficiency.

## Step-by-Step Process

### 1. **Function Selection**
   - The user is presented with a menu to choose a mathematical function type to define the differential equation, such as:
     - Polynomial
     - Trigonometric functions (sine, cosine, tangent)
     - Logarithmic function
     - Custom functions of `x` and `y`, including terms like `x^2`, `xy`, and `y^2`
   - Based on the user’s choice, the program prompts for coefficients or parameters required to define the chosen function.
   - Once the parameters are input, the application formulates the equation and displays it to the user.

### 2. **Initial Conditions and Step Size Setup**
   - The user provides:
     - The initial values of `x` and `y` (denoted as `x0` and `y0`).
     - A non-zero step size `h` that controls the increments for the `x` values.
     - The `range` for which the solution is computed, marking the endpoint of the calculation.

### 3. **Runge-Kutta Calculation**
   - The application uses the Runge-Kutta method (fourth-order) to solve the differential equation iteratively from `x0` to the specified range.
   - **Intermediate Calculations**:
     - At each step, the method computes four intermediate slopes:
       - `k1 = h * f(x, y)`
       - `k2 = h * f(x + h/2, y + k1/2)`
       - `k3 = h * f(x + h/2, y + k2/2)`
       - `k4 = h * f(x + h, y + k3)`
     - These slopes represent estimates of the derivative at different points within the current interval.
   - **Updating the Solution**:
     - The next value of `y` is calculated by averaging these slopes:
       - `y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6`
     - The `x` value is incremented by `h` to proceed to the next step.

### 4. **Result Display**
   - At each iteration, the new values of `x` and `y` are printed in a table format for easy tracking of the solution progression.

### 5. **Completion**
   - The process repeats until the solution reaches the specified range, providing a step-by-step solution for the differential equation using the selected function and initial values.

This structured approach ensures that the solution is computed accurately over the desired interval.
# Working Principle of Matrix Inversion (Gauss-Jordan Method)

This application calculates the inverse of a square matrix using the Gauss-Jordan method. The process includes checking if the matrix is invertible before attempting inversion.

## Step-by-Step Process

### 1. **Invertibility Check**
   - Before inverting, the function checks if the matrix is invertible by attempting to calculate its determinant.
   - Steps:
     - For each pivot position, the function ensures a non-zero value by swapping rows if needed. If no non-zero pivot is found in a column, the matrix is deemed singular (non-invertible), and the function returns `false`.
     - If a pivot exists, the function reduces the matrix to an upper triangular form, calculating the determinant in the process.
     - If the determinant is non-zero, the matrix is invertible.

### 2. **Initial Setup for Inversion**
   - If the matrix is confirmed to be invertible, an identity matrix of the same size is created. This identity matrix will be transformed into the inverse of `A` as the matrix is reduced.
   - The matrix `A` is augmented with this identity matrix, forming an augmented matrix `[A | I]`.

### 3. **Row Operations to Create Pivots**
   - For each row `i`, the algorithm identifies the pivot element and normalizes it:
     - **Pivot Normalization**:
       - Divide the entire row by the pivot element so that the pivot becomes 1. This operation is applied to both `A` and `inv` (initially the identity matrix).
     - **Zeroing Elements Above and Below the Pivot**:
       - For all rows except the pivot row `i`, adjust each element in the column to zero by subtracting a multiple of the pivot row.
       - This step ensures that all elements above and below each pivot are zero, moving the matrix closer to reduced row echelon form.

### 4. **Iterate Over Each Row**
   - Repeat the row operations for each row in `A`, moving down the diagonal until the entire left side of the augmented matrix becomes the identity matrix.
   - As `A` transforms into the identity matrix, `inv` transforms into the inverse of `A`.

### 5. **Return the Inverse Matrix**
   - Once the left half of the augmented matrix is reduced to the identity matrix, the right half contains the inverse of the original matrix `A`.
   - The `invert` function returns this inverse matrix.

### 6. **Failure Case Handling**
   - If the matrix is not invertible (as detected by `isInvertible`), the application can return an error message indicating the matrix is singular.

This structured approach ensures that only invertible matrices are processed, and the result is computed accurately using Gauss-Jordan elimination.

## Contributors (Team No 18)
- Ayesha Mehereen - 2107039  [@Mehereen-1](https://github.com/Mehereen-1)
- Shormi Ghosh - 21070109 [@ShormiGhosh](https://github.com/ShormiGhosh)
- Md. Nayeem - 2107050 [@codernayeem](https://github.com/codernayeem)
