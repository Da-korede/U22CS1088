# Simplex Algorithm in C++
#Matric No= U22CS1088

## Description
This project implements the **Simplex Algorithm** in C++ to solve linear programming problems. The Simplex method is used to maximize or minimize a linear objective function subject to linear constraints.

## Features
- Supports standard form linear programming problems (maximization).
- Uses the **tableau method** for pivoting and iterating towards optimality.
- Detects unbounded solutions.
- Extracts the optimal values for the original variables.

## Problem Formulation
The algorithm solves linear programming problems of the form:

**Maximize:**
\[ Z = c_1x_1 + c_2x_2 + ... + c_nx_n \]

## How It Works
1. **Initialize the tableau** by adding slack variables to convert inequalities to equalities.
2. **Select the entering variable** (most negative coefficient in the objective row).
3. **Select the leaving variable** using the minimum ratio test.
4. **Pivot** to update the tableau.
5. **Repeat** until all coefficients in the objective row are non-negative.
6. **Extract the optimal solution** and return the results.

## Compilation & Execution
To compile and run the program, use the following commands:

```sh
g++ simplex.cpp -o simplex
./simplex
```

