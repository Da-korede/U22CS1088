#include <iostream>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
// Matric No = U22CS1088
//Owner Name= Muhammad Abass Mudasir

using namespace std;

// Returns a pair: {optimal_value, solution vector}
pair<double, vector<double>> simplex(const vector<vector<double>> &A,
                                      const vector<double> &b,
                                      const vector<double> &c) {
    int m = A.size();       // number of constraints
    int n = c.size();       // number of original variables

    // Build the simplex tableau.
    // Dimensions: (m + 1) rows and (n + m + 1) columns.
    // The first m rows correspond to constraints; the last row is for the objective.
    vector<vector<double>> tableau(m + 1, vector<double>(n + m + 1, 0.0));

    // Set up the constraint rows.
    for (int i = 0; i < m; i++) {
        // Coefficients for original variables.
        for (int j = 0; j < n; j++) {
            tableau[i][j] = A[i][j];
        }
        // Slack variable coefficient (identity matrix).
        tableau[i][n + i] = 1.0;
        // Right-hand side value.
        tableau[i].back() = b[i];
    }

    // Set up the objective row.
    // For maximization, we put -c in the objective row.
    for (int j = 0; j < n; j++) {
        tableau[m][j] = -c[j];
    }

    // Main simplex iterations.
    while (true) {
        // 1. Choose the entering variable (pivot column)
        int pivot_col = -1;
        for (int j = 0; j < n + m; j++) {
            if (tableau[m][j] < -1e-9) { // use a small threshold for numerical stability
                pivot_col = j;
                break; // choose the first column with a negative coefficient
            }
        }
        // If no negative coefficients remain, we have reached optimality.
        if (pivot_col == -1) break;

        // 2. Choose the leaving variable (pivot row) using the minimum ratio test.
        double min_ratio = numeric_limits<double>::infinity();
        int pivot_row = -1;
        for (int i = 0; i < m; i++) {
            if (tableau[i][pivot_col] > 1e-9) { // only consider positive entries
                double ratio = tableau[i].back() / tableau[i][pivot_col];
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    pivot_row = i;
                }
            }
        }

        // If no pivot row is found, the LP is unbounded.
        if (pivot_row == -1) {
            throw runtime_error("Linear program is unbounded.");
        }

        // 3. Pivot: make the pivot element 1 and adjust the tableau.
        double pivot_value = tableau[pivot_row][pivot_col];
        // Normalize the pivot row.
        for (int j = 0; j < n + m + 1; j++) {
            tableau[pivot_row][j] /= pivot_value;
        }
        // Zero out the pivot column in the other rows.
        for (int i = 0; i < m + 1; i++) {
            if (i == pivot_row) continue;
            double factor = tableau[i][pivot_col];
            for (int j = 0; j < n + m + 1; j++) {
                tableau[i][j] -= factor * tableau[pivot_row][j];
            }
        }
    }

    // Extract the solution for the original variables.
    vector<double> solution(n, 0.0);
    for (int j = 0; j < n; j++) {
        int basic_row = -1;
        bool is_basic = true;
        // Check if column j is a unit column.
        for (int i = 0; i < m; i++) {
            if (fabs(tableau[i][j] - 1.0) < 1e-9) {
                if (basic_row == -1)
                    basic_row = i;
                else {
                    is_basic = false;
                    break;
                }
            } else if (fabs(tableau[i][j]) > 1e-9) {
                is_basic = false;
                break;
            }
        }
        if (is_basic && basic_row != -1) {
            solution[j] = tableau[basic_row].back();
        } else {
            solution[j] = 0.0;
        }
    }

    double optimal_value = tableau[m].back();
    return {optimal_value, solution};
}

int main() {
    // Example linear program:
    // Maximize: 3x + 2y
    // Subject to:
    //    x + y <= 4
    //    x     <= 2
    //         y <= 3
    //    x, y >= 0

    // Coefficients for constraints (A * x <= b)
    vector<vector<double>> A = {
        {1, 1},
        {1, 0},
        {0, 1}
    };

    // Right-hand side values.
    vector<double> b = {4, 2, 3};

    // Coefficients for the objective function.
    vector<double> c = {3, 2};

    try {
        auto result = simplex(A, b, c);
        double optimal_value = result.first;
        vector<double> solution = result.second;

        cout << "Optimal value: " << optimal_value << endl;
        for (size_t i = 0; i < solution.size(); i++) {
            cout << "x[" << i << "] = " << solution[i] << endl;
        }
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
#include <iostream>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>

using namespace std;

// Returns a pair: {optimal_value, solution vector}
pair<double, vector<double>> simplex(const vector<vector<double>> &A,
                                      const vector<double> &b,
                                      const vector<double> &c) {
    int m = A.size();       // number of constraints
    int n = c.size();       // number of original variables

    // Build the simplex tableau.
    // Dimensions: (m + 1) rows and (n + m + 1) columns.
    // The first m rows correspond to constraints; the last row is for the objective.
    vector<vector<double>> tableau(m + 1, vector<double>(n + m + 1, 0.0));

    // Set up the constraint rows.
    for (int i = 0; i < m; i++) {
        // Coefficients for original variables.
        for (int j = 0; j < n; j++) {
            tableau[i][j] = A[i][j];
        }
        // Slack variable coefficient (identity matrix).
        tableau[i][n + i] = 1.0;
        // Right-hand side value.
        tableau[i].back() = b[i];
    }

    // Set up the objective row.
    // For maximization, we put -c in the objective row.
    for (int j = 0; j < n; j++) {
        tableau[m][j] = -c[j];
    }

    // Main simplex iterations.
    while (true) {
        // 1. Choose the entering variable (pivot column)
        int pivot_col = -1;
        for (int j = 0; j < n + m; j++) {
            if (tableau[m][j] < -1e-9) { // use a small threshold for numerical stability
                pivot_col = j;
                break; // choose the first column with a negative coefficient
            }
        }
        // If no negative coefficients remain, we have reached optimality.
        if (pivot_col == -1) break;

        // 2. Choose the leaving variable (pivot row) using the minimum ratio test.
        double min_ratio = numeric_limits<double>::infinity();
        int pivot_row = -1;
        for (int i = 0; i < m; i++) {
            if (tableau[i][pivot_col] > 1e-9) { // only consider positive entries
                double ratio = tableau[i].back() / tableau[i][pivot_col];
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    pivot_row = i;
                }
            }
        }

        // If no pivot row is found, the LP is unbounded.
        if (pivot_row == -1) {
            throw runtime_error("Linear program is unbounded.");
        }

        // 3. Pivot: make the pivot element 1 and adjust the tableau.
        double pivot_value = tableau[pivot_row][pivot_col];
        // Normalize the pivot row.
        for (int j = 0; j < n + m + 1; j++) {
            tableau[pivot_row][j] /= pivot_value;
        }
        // Zero out the pivot column in the other rows.
        for (int i = 0; i < m + 1; i++) {
            if (i == pivot_row) continue;
            double factor = tableau[i][pivot_col];
            for (int j = 0; j < n + m + 1; j++) {
                tableau[i][j] -= factor * tableau[pivot_row][j];
            }
        }
    }

    // Extract the solution for the original variables.
    vector<double> solution(n, 0.0);
    for (int j = 0; j < n; j++) {
        int basic_row = -1;
        bool is_basic = true;
        // Check if column j is a unit column.
        for (int i = 0; i < m; i++) {
            if (fabs(tableau[i][j] - 1.0) < 1e-9) {
                if (basic_row == -1)
                    basic_row = i;
                else {
                    is_basic = false;
                    break;
                }
            } else if (fabs(tableau[i][j]) > 1e-9) {
                is_basic = false;
                break;
            }
        }
        if (is_basic && basic_row != -1) {
            solution[j] = tableau[basic_row].back();
        } else {
            solution[j] = 0.0;
        }
    }

    double optimal_value = tableau[m].back();
    return {optimal_value, solution};
}

int main() {
    // Example linear program:
    // Maximize: 3x + 2y
    // Subject to:
    //    x + y <= 4
    //    x     <= 2
    //         y <= 3
    //    x, y >= 0

    // Coefficients for constraints (A * x <= b)
    vector<vector<double>> A = {
        {1, 1},
        {1, 0},
        {0, 1}
    };

    // Right-hand side values.
    vector<double> b = {4, 2, 3};

    // Coefficients for the objective function.
    vector<double> c = {3, 2};

    try {
        auto result = simplex(A, b, c);
        double optimal_value = result.first;
        vector<double> solution = result.second;

        cout << "Optimal value: " << optimal_value << endl;
        for (size_t i = 0; i < solution.size(); i++) {
            cout << "x[" << i << "] = " << solution[i] << endl;
        }
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
