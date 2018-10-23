#include <stdio.h>

#include "common.c"

int sor(double *x, double *y, int n, double w) {

  // TODO establish actual tol needed
  const double tol = 5e-13;

  // To ease notation and calculation.
  double t = 1/3.;
  double f = 1/4.;
  double i = 1/((double) n);

  // Number of iterations needed
  int it = 0;

  // Main algorithm loop. Iterate until error is small enough.
  do {

    // One iteration more
    it++;

    // Copy y to x for the next iteration.
    for (int j = 0; j < n; j++) {
      x[j] = y[j];
    }

    // We'll calculate each entry of x separately.
    // There are 3 main parts:
    // 1. We calculate x_0 and x_1.
    // 2. We calculate x_2 through x_{n-3}.
    // 3. We calculate x_{n-2} and x_{n-1}.
    // This is because B has negative entries in the two top and two bottom rows.

    // 1. We calculate x_0 and x_1.
    y[0] = t*(x[2]-x[n-2]+2*i);
    y[1] = f*(x[3]-x[n-1]+2*i);

    // 2. We calculate x_2 through x_{n-3}.
    for (int j = 2; j < n-2; j++) {
      y[j] = x[j] + w*(j%2 ? f : t)*(y[j-2] - (j%2 ? 4 : 3)*x[j] + x[j+2] + i*(j+2-(j%2)));
    }

    // 3. We calculate x_{n-2} and x_{n-1}.
    y[n-2] = t*(x[n-4]-y[0]+1);
    y[n-1] = f*(x[n-3]-y[1]+1);

    // We've finished calculating the next iteration of x.
    // Now we calculate the norm of the difference between this iteration (x)
    // and the previous one (y) to see whether we're under the given error.
  } while (distance(x, y, n) > tol);

  return it;

}

// Try with different values for w and return the best one
double find_w(double *x, double *y, int n) {

  double w = 0.;
  double w_best = 0.;

  unsigned int it = 0;
  unsigned int it_best = -1;

  const double w_min = 0.1;
  const double w_max = 2.;
  const double w_step = 0.1;

  for (w = w_min; w < w_max; w += w_step) {

    for (int j = 0; j < n; j++) {
      x[j] = 0.;
      y[j] = 0.;
    }

    it = sor(x, y, n, w);
    if (it < it_best) {
      it_best = it;
      w_best = w;
    }

  }

  // Clean the vectors one last time
  for (int j = 0; j < n; j++) {
    x[j] = 0.;
    y[j] = 0.;
  }

  return w_best;

}

int main() {

  // Size of the vectors, established by the problem statement.
  const int n = 1e+6;

  // To store the solution in.
  static double x[n] = {0.};

  // Here we'll store the next iteration of x and then we'll copy it over x again.
  static double y[n] = {0.};

  // SOR constant; best value found experimentally to be 1.2
  double w = 1.2; // = find_w(x, y, n);

  sor(x, y, n, w);

  for (int j=0; j<n; j++){
    printf("%.12f\n", x[j]);
  }

}
