#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.67408e-11  // gravitational constant
#define DT 0.001  // time step
#define MAX_ITERATIONS 10000  // maximum number of iterations

typedef struct {
    double mass;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
} Body;


void cleanup(Body **bodies, int num_bodies) {
    for(int i = 0; i < num_bodies; i++) {
        free(bodies[i]);
    }
}


Body **generate_initial_population(int n, double max_mass, double max_pos, double max_vel) {
    Body **bodies = malloc(sizeof(Body *) * n);

    for(int i = 0; i < n; i++) {
        bodies[i] = malloc(sizeof(Body));
        bodies[i]->mass = (double)rand() / RAND_MAX * max_mass;
        bodies[i]->x = (double)rand() / RAND_MAX * max_pos;
        bodies[i]->y = (double)rand() / RAND_MAX * max_pos;
        bodies[i]->z = (double)rand() / RAND_MAX * max_pos;
        bodies[i]->vx = (double)rand() / RAND_MAX * max_vel;
        bodies[i]->vy = (double)rand() / RAND_MAX * max_vel;
        bodies[i]->vz = (double)rand() / RAND_MAX * max_vel;
    }
    return bodies;
}

void print_boddies(Body **bodies, int n, int iteration) {
    if (iteration == 0)
        printf("iteration,mass,x,y,z,vx,vy,vz\n");
    for(int i = 0; i < n; i++) 
        printf("%d,%f,%f,%f,%f,%f,%f,%f\n", iteration, bodies[i]->mass, bodies[i]->x, bodies[i]->y, bodies[i]->z, bodies[i]->vx, bodies[i]->vy, bodies[i]->vz);
}


Body **calculate_iteration(Body** bodies, int num_bodies) {
    Body** new_bodies = (Body**) malloc(num_bodies * sizeof(Body*));
    for (int i = 0; i < num_bodies; i++) {
        new_bodies[i] = malloc(sizeof(Body));
        double forceX = 0;
        double forceY = 0;
        double forceZ = 0;
        for (int j = 0; j < num_bodies; j++) {
            double rx = bodies[i]->x - bodies[j]->x;
            double ry = bodies[i]->y - bodies[j]->y;
            double rz = bodies[i]->z - bodies[j]->z;
            double distance = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

            forceX += bodies[j]->mass * rx / (pow(distance, 3) + G);// G je epsilon xd
            forceY += bodies[j]->mass * ry / (pow(distance, 3) + G); // G je epsilon xd
            forceZ += bodies[j]->mass * rz / (pow(distance, 3) + G); // G je epsilon xd
        }      

        forceX *= G * bodies[i]->mass;  
        forceY *= G * bodies[i]->mass;  
        forceZ *= G * bodies[i]->mass;  

        double ax = forceX / bodies[i]->mass;
        double ay = forceY / bodies[i]->mass;
        double az = forceZ / bodies[i]->mass;

        new_bodies[i]->mass = bodies[i]->mass;
        new_bodies[i]->x = bodies[i]->mass + bodies[i]->vx * DT + 0.5 * ax * pow(DT, 2);
        new_bodies[i]->y = bodies[i]->mass + bodies[i]->vy * DT + 0.5 * ay * pow(DT, 2);
        new_bodies[i]->z = bodies[i]->mass + bodies[i]->vz * DT + 0.5 * az * pow(DT, 2);

        new_bodies[i]->vx = bodies[i]->vx + ax * DT;
        new_bodies[i]->vy = bodies[i]->vy + ay * DT;
        new_bodies[i]->vz = bodies[i]->vz + az * DT;

    }

    return new_bodies;
}


int main() {
    Body **bodies;
    bodies = generate_initial_population(100, 100, 100, 100);
    for(int i = 0; i < 100; i ++) {
        bodies = calculate_iteration(bodies, 100);
        print_boddies(bodies,100, i);
    }
    return 0;
}
