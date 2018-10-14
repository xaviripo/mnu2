
#include <stdio.h>
#include <math.h>

// Calculates the infinity norm-induced distance of n-sized vectors x and y
double norm(double *x, double *y, size_t n) {

  double diff = 0;
  double entry_diff = 0;

  for (int j = 0; j < n; j++) {
    entry_diff = fabs(y[j]-x[j]);
    diff = diff>entry_diff ? diff : entry_diff;
  }

  return diff;

}

int main() {

  // Declaration of variables
  // Size of the vectors, established by the problem statement.
  const int n = 1e+6;
  // Error, established by the problem statement.
  // By a proposition, we know that the error to the actual solution is
  // smaller or equal to the error between two iterations multiplied by
  // \beta/(1+\beta), where \beta is the norm of the B_J matrix.
  // This norm has been calculated to be 2/3 in the previous point.
  // Then, the error to the actual solution is lower or equal than twice
  // the error between iterations, so we take 1/2 of the given tolerance,
  // or 0.5e-12 = 5e-13
  const double tol = 5e-13;
  // To store the solution in.
  static double x[n] = {0.};
  // Here we'll store the next iteration of x and then we'll copy it over x again.
  static double y[n] = {0.};

  // To ease notation and calculation.
  double t = 1/3.;
  double f = 1/4.;
  double i = 1/((double) n);

  double diff = 0;

  // Main algorithm loop. Iterate until error is small enough.
  do {

    // Copy y to x for the next iteration.
    for (int j = 0; j < n; j++) {
      x[j] = y[j];
    }

    // We'll calculate each entry of x separately.
    // There are 3 main parts:
    // 1. We calculate x_0 and x_1.
    // 2. We calculate x_2 through x_{n-3}.
    // 3. We calculate x_{n-2} and x_{n-1}.
    // This is because B_J has negative entries in the two top and two bottom rows.

    // 1. We calculate x_0 and x_1.
    y[0] = t*(x[2]-x[n-2]+2*i);
    y[1] = f*(x[3]-x[n-1]+2*i);

    // 2. We calculate x_2 through x_{n-3}.
    // Even subindices are multiplied by 1/3, odd ones by 1/4.
    // Even subindices:
    for (int j = 2; j < n-3; j += 2) {
      y[j] = t*(x[j-2] + x[j+2] + i*(j+2));
    }

    // Odd subindices:
    for (int j = 3; j < n-2; j += 2) {
      y[j] = f*(x[j-2] + x[j+2] + i*(j+1));
    }

    // 3. We calculate x_{n-2} and x_{n-1}.
    y[n-2] = t*(x[n-4]-x[0]+1);
    y[n-1] = f*(x[n-3]-x[1]+1);

    // We've finished calculating the next iteration of x.
    // Now we calculate the norm of the difference between this iteration (x)
    // and the previous one (y) to see whether we're under the given error.
  } while (norm(x, y, n) > tol);

  for (int j=0; j<n; j++){
    printf("x[%d] = %f\n", j, x[j]);
  }

}
