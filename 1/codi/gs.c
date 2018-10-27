#include <stdio.h>

#include "common.c"

// Size of the vectors, established by the problem statement.
#define N 1000000

int gs(double *x, double *y, int n) {

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
      y[j] = (j%2 ? f : t)*(y[j-2] + x[j+2] + i*(j+2-(j%2)));
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

int main() {

  // To store the solution in.
  static double x[N] = {0.};

  // Here we'll store the next iteration of x and then we'll copy it over x again.
  static double y[N] = {0.};

  int it = gs(x, y, N);

  for (int j = 0; j < N; j++){
    printf("%.12f\n", x[j]);
  }

  // To see the number of iterations
  //printf("%d", it);

}
