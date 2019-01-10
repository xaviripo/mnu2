#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Change following definition to 1 for debugging trace:
#define DEBUG 0
#define debug if(DEBUG) printf


typedef struct v2 {
  double x;
  double y;
} v2;

double norm(v2 a) {
  return pow(pow(a.x, 2.) + pow(a.y, 2.), .5);
}

double distance(v2 a, v2 b) {
  return pow(pow(a.x - b.x, 2.) + pow(a.y - b.y, 2.), .5);
}

double f(double x, double y) {
  return +1*x*x*x*x +
    +2*x*x*x*y +
    +3*x*x*y*y +
    +2*x*y*y*y +
    +2*y*y*y*y +
    -2.6*x*x*x +
    -3.2*x*x*y +
    -2.6*x*y*y +
    -3.2*y*y*y +
    -7.22*x*x +
    -16*x*y +
    -15.22*y*y +
    +20.8*x +
    +25.6*y +
    -5.94;
}

v2 grad_f(double x, double y) {
  return (v2) {

    +4*x*x*x +
    +6*x*x*y +
    +6*x*y*y +
    +2*y*y*y +
    -7.8*x*x +
    -6.4*x*y +
    -2.6*y*y +
    -14.44*x +
    -16*y +
    +20.8,

    +2*x*x*x +
    +6*x*x*y +
    +6*x*y*y +
    +8*y*y*y +
    -3.2*x*x +
    -5.2*x*y +
    -9.6*y*y +
    -16*x +
    -30.44*y +
    +25.6

  };
}

v2 get_initial_point(double tol) {

  // We use Newton's method to find an initial solution on y=0
  double x = 0.; // x_0
  double x_next = 0.;

  do {
    x = x_next;
    x_next = x - f(x, 0.)/grad_f(x, 0.).x;
  } while (fabs(f(x_next, 0.)) > tol);

  return (v2) {x_next, 0.};

}

v2 get_next_point(v2 point, v2 vect, double dis, double tol) {

  v2 point_next;
  point_next.x = point.x + vect.x;
  point_next.y = point.y + vect.y;

  v2 point_init;
  point_init.x = point.x;
  point_init.y = point.y;

  double dissqr = dis*dis;
  double
    f_xy,   // f(x,y)
    g_xy,   // g(x,y) = (x-x0)^2 + (y-y0)^2 - dis^2
    df_x,   // (df/dx)(x,y)
    df_y,   // (df/dy)(x,y)
    dg_x,   // (dg/dx)(x,y)
    dg_y,   // (dg/dy)(x,y)
    scalar; // -1/det(DF(x,y))
  v2 df_xy;  // (grad f)(x,y)

  do {

    debug("\tNext point: (%g, %g)\n", point_next.x, point_next.y);

    point = point_next;

    // F = (f,g)
    // f(x,y) as defined by the problem statement
    f_xy = f(point.x, point.y);
    // g(x,y) = (x-x0)^2 + (y-y0)^2 - dis^2
    g_xy = pow(point.x - point_init.x, 2.) + pow(point.y - point_init.y, 2.) - dissqr;

    // (grad F)(x,y)
    df_xy = grad_f(point.x, point.y);
    df_x = df_xy.x;
    df_y = df_xy.y;
    dg_x = 2*(point.x - point_init.x);
    dg_y = 2*(point.y - point_init.y);
    scalar = 1./(df_x*dg_y - df_y*dg_x);
    debug("\tdf_x: %g\n", df_x);
    debug("\tdg_y: %g\n", dg_y);
    debug("\tdf_y: %g\n", df_y);
    debug("\tdg_x: %g\n", dg_x);
    debug("\tScalar: %g\n", scalar);

    // Apply Newton's method with F=(f,g)
    point_next.x = point.x - scalar*(dg_y*f_xy - df_y*g_xy);
    point_next.y = point.y - scalar*(-dg_x*f_xy + df_x*g_xy);
    debug("\tf(%g, %g) = %g\n", point_next.x, point_next.y, f(point_next.x, point_next.y));

  } while (fabs(f(point_next.x, point_next.y)) > tol);

  return point_next;

}

void get_list_points() {

  FILE *f = fopen("output.txt", "w");
  if (f == NULL) {
    printf("Error opening file\n");
    exit(1);
  }

  const double tol = 1e-10;
  const double dis = .01;

  const v2 point_init = get_initial_point(tol);

  v2 point = point_init;
  v2 point_next = point;

  v2 grad, vect, vect_next;

  fprintf(f, "%.10f %.10f\n", point_next.x, point_next.y);

  double vect_norm;

  int i = 0;
  do {
    debug("-------- ITERATION %d --------\n", ++i);

    // Update point and vector from previous iteration
    debug("Update point and vector from previous iteration\n");
    point = point_next;
    vect = vect_next;
    debug("\tPoint: (%g, %g)\n", point.x, point.y);
    debug("\tVector: (%g, %g)\n", vect.x, vect.y);

    // Calculate next initial approximation
    debug("Calculate next initial approximation\n");
    grad = grad_f(point.x, point.y);
    debug("\tGradient: (%g, %g)\n", grad.x, grad.y);
    vect_next = (v2) {-grad.y, grad.x};

    // Change it to length dis
    debug("Change it to length dis\n");
    debug("\tPre-scaled vector: (%g, %g)\n", vect_next.x, vect_next.y);
    vect_norm = norm(vect_next);
    vect_next.x = dis*vect_next.x/vect_norm;
    vect_next.y = dis*vect_next.y/vect_norm;
    debug("\tScaled vector: (%g, %g)\n", vect_next.x, vect_next.y);

    // Make sure it's got the correct direction
    debug("Make sure it's got the correct direction\n");
    debug("\tOld direction: (%g, %g)\n", vect_next.x, vect_next.y);
    if (vect.x*vect_next.x + vect.y*vect_next.y < 0) {
      vect_next.x = -vect_next.x;
      vect_next.y = -vect_next.y;
    }
    debug("\tNew direction: (%g, %g)\n", vect_next.x, vect_next.y);

    // Get the next point
    debug("Get the next point\n");
    point_next = get_next_point(point, vect_next, dis, tol);
    debug("\tPoint next: (%g, %g)\n", point_next.x, point_next.y);

    fprintf(f, "%.10f %.10f\n", point_next.x, point_next.y);

  } while (distance(point_init, point_next) > dis);

  fclose(f);

}

int main() {
  debug("\x1B[36m");
  get_list_points();
  debug("\x1B[0m");
}
