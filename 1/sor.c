

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
  // TODO establish actual tol needed
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

  // SOR constant
  double w = 1.0;

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
    y[0] = x[0] + w*t*(x[2]-3*x[0]-x[n-2]+2*i);
    y[1] = x[1] + w*f*(x[3]-4*x[1]-x[n-1]+2*i);

    // 2. We calculate x_2 through x_{n-3}.
    // Even subindices are multiplied by 1/3, odd ones by 1/4.
    for (int j = 2; j < n-2; j++) {
      y[j] = x[j] + w*(j%2 ? f : t)*(y[j-2] - (j%2 ? 4 : 3)*x[j] + x[j+2] + i*(j+2-(j%2)));
    }

    // 3. We calculate x_{n-2} and x_{n-1}.
    y[n-2] = x[n-2] + w*t*(x[n-4]-3*x[n-2]-y[0]+1);
    y[n-1] = x[n-1] + w*f*(x[n-3]-4*x[n-1]-y[1]+1);

    // We've finished calculating the next iteration of x.
    // Now we calculate the norm of the difference between this iteration (x)
    // and the previous one (y) to see whether we're under the given error.
  } while (norm(x, y, n) > tol);

  for (int j=0; j<n; j++){
    printf("x[%d] = %.10f\n", j, x[j]);
  }

}
