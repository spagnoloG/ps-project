#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define G 1 // gravitational constant
#define DT 0.001  // time step
#define N_ITER 100000
#define N_BODIES 10
#define EPSILON 1

struct timeval  tv1, tv2;

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
    srand(42);
    for(int i = 0; i < n; i++) {
        bodies[i] = malloc(sizeof(Body));
        bodies[i]->mass = 20;
        bodies[i]->x = (double)rand() / RAND_MAX * max_pos;
        bodies[i]->y = (double)rand() / RAND_MAX * max_pos;
        bodies[i]->z = (double)rand() / RAND_MAX * max_pos;
        if(i == 0) {
            bodies[i]->vx = 0;
            bodies[i]->vy = 0;
            bodies[i]->vz = 0;
        } else {
            bodies[i]->vx = (double)rand() / RAND_MAX * max_vel;
            bodies[i]->vy = (double)rand() / RAND_MAX * max_vel;
            bodies[i]->vz = (double)rand() / RAND_MAX * max_vel;
        }
    }
    return bodies;
}

void print_boddies(Body **bodies, int n, int iteration) {
    //if (iteration == 0)
    //    printf("iteration,mass,x,y,z,vx,vy,vz\n");
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
            double rx = bodies[j]->x - bodies[i]->x;
            double ry = bodies[j]->y - bodies[i]->y;
            double rz = bodies[j]->z - bodies[i]->z;
            double distance = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

            forceX += bodies[j]->mass * rx / (pow(distance, 3) + EPSILON);// G je epsilon xd
            forceY += bodies[j]->mass * ry / (pow(distance, 3) + EPSILON); // G je epsilon xd
            forceZ += bodies[j]->mass * rz / (pow(distance, 3) + EPSILON); // G je epsilon xd
        }      

        forceX *= G * bodies[i]->mass;  
        forceY *= G * bodies[i]->mass;  
        forceZ *= G * bodies[i]->mass;  

        double ax = forceX / bodies[i]->mass;
        double ay = forceY / bodies[i]->mass;
        double az = forceZ / bodies[i]->mass;

        new_bodies[i]->mass = bodies[i]->mass;
        new_bodies[i]->x = bodies[i]->x + bodies[i]->vx * DT + 0.5 * ax * pow(DT, 2);
        new_bodies[i]->y = bodies[i]->y + bodies[i]->vy * DT + 0.5 * ay * pow(DT, 2);
        new_bodies[i]->z = bodies[i]->z + bodies[i]->vz * DT + 0.5 * az * pow(DT, 2);

        new_bodies[i]->vx = bodies[i]->vx + ax * DT;
        new_bodies[i]->vy = bodies[i]->vy + ay * DT;
        new_bodies[i]->vz = bodies[i]->vz + az * DT;
    }

    cleanup(bodies, num_bodies);
    return new_bodies;
}

int main() {
    Body **bodies;
    bodies = generate_initial_population(N_BODIES, 100000000, 3, 3);
    gettimeofday(&tv1, NULL);
    for(int i = 0; i < N_ITER; i ++) {
        print_boddies(bodies,N_BODIES, i);
        bodies = calculate_iteration(bodies, N_BODIES);
    }
    gettimeofday(&tv2, NULL);
    //printf ("CPU: %f seconds\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
    //     (double) (tv2.tv_sec - tv1.tv_sec)); // Do not time for now

    print_boddies(bodies,N_BODIES, N_ITER);
    return 0;
}
