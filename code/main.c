#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.67408e-11  // gravitational constant
#define DT 0.001  // time step
#define MAX_ITERATIONS 10000  // maximum number of iterations

// structure to represent a celestial body
typedef struct {
  double mass;  // mass of the body
  double x;  // x-coordinate of the body's position
  double y;  // y-coordinate of the body's position
  double vx;  // x-component of the body's velocity
  double vy;  // y-component of the body's velocity
} Body;

// function to calculate the acceleration of a body due to the gravitational forces of all other bodies
void calc_acceleration(Body *bodies, int n, int i, double *ax, double *ay) {
  *ax = 0;
  *ay = 0;
  for (int j = 0; j < n; j++) {
    if (i == j) {
      continue;
    }
    double dx = bodies[j].x - bodies[i].x;
    double dy = bodies[j].y - bodies[i].y;
    double r2 = dx * dx + dy * dy;
    double r = sqrt(r2);
    double F = G * bodies[i].mass * bodies[j].mass / r2;
    *ax += F * dx / r;
    *ay += F * dy / r;
  }
}

int main() {
  // create an array of celestial bodies
  int n = 3;
  Body bodies[3] = {
    {1.0, 0.0, 0.0, 0.0, 0.0},  // sun
    {1.0e-3, 1.0, 0.0, 0.0, 2 * M_PI},  // earth
    {2.0e-3, 0.0, 1.0, -2 * M_PI, 0.0}  // moon
  };

  // simulate the motion of the bodies
  for (int t = 0; t < MAX_ITERATIONS; t++) {
    for (int i = 0; i < n; i++) {
      double ax, ay;
      calc_acceleration(bodies, n, i, &ax, &ay);
      bodies[i].vx += ax * DT;
      bodies[i].vy += ay * DT;
      bodies[i].x += bodies[i].vx * DT;
      bodies[i].y += bodies[i].vy * DT;
    }
  }

  // print the final positions of the bodies
  for (int i = 0; i < n; i++) {
    printf("Body %d: (%.2f, %.2f)\n", i, bodies[i].x, bodies[i].y);
  }

  return 0;
}
