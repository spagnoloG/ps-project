/*
 * PS final projcet. Problem: n-bodies
 * Authors: Marcel Rucigoj, Aleks Stepancic, Gasper Spagnolo
 * FRI
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <errno.h>
#include <execinfo.h>
#include <signal.h>

#define G 1                         // gravitational contant
#define DT 0.001                    // time derivative
#define EPSILON 1                   // epsilon to avoid division by 0
#define LOG_FILE "output_mp.txt"// logfile location

struct timeval  tv1, tv2;
struct winsize w;

typedef struct {
    float mass;
    float x;
    float y;
    float z;
    float vx;
    float vy;
    float vz;
} Body;

/*
 * error handling 
 */
void segfault_handler(int sig) {
    void *array[10];
    size_t size;

    size = backtrace(array, 10);

    fprintf(stderr, "[-] Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

void throw_err(int line, int err_code) {
    char err_buf[256];

    sprintf(err_buf, "[-] Error: %s -> %d", __FILE__, line);
    perror(err_buf);
    exit(err_code);
}

/*
 * misc functions
 */
void cleanup(Body **bodies, int num_bodies) {
    for(int i = 0; i < num_bodies; i++) {
        free(bodies[i]);
    }
}

void print_boddies(Body **bodies, int n, int iteration, FILE *logfile) {
    for(int i = 0; i < n; i++) {
        fprintf(logfile,"%d,%f,%f,%f,%f,%f,%f,%f\n", iteration, bodies[i]->mass, bodies[i]->x, bodies[i]->y, bodies[i]->z, bodies[i]->vx, bodies[i]->vy, bodies[i]->vz);
    }
}

void progress_bar(int current, int max) {
    int i;
    char pb_buffer[256]; 
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    int width = w.ws_col;
    float progress = (float)current / (float)max; // Calculate progress as a fraction
    int chars_printed = (int)(progress * width); // Calculate number of characters to print
                                                 //
    sprintf(pb_buffer, "[ %d/%d ] ", current, max);
    long pb_buf_len = strlen(pb_buffer);
    for(int i = 0; i < pb_buf_len; i++) {
        printf("%c", pb_buffer[i]);
    }
    
    for (i = 0; i < chars_printed - pb_buf_len; i++) {
        printf("=");
    }

    printf("\r");
    fflush(stdout);
}

/*
 * Algorithm
 */
Body **generate_initial_population(int n, float max_mass, float max_pos, float max_vel) {
    Body **bodies = malloc(sizeof(Body *) * n);
    if(bodies == NULL)
        throw_err(__LINE__, 1);
    
    srand(42);
    for(int i = 0; i < n; i++) {
        bodies[i] = malloc(sizeof(Body));
        if(bodies[i] == NULL)
            throw_err(__LINE__, 1);

        bodies[i]->mass = (float)rand() / RAND_MAX * max_mass;
        bodies[i]->x = (float)rand() / RAND_MAX * max_pos;
        bodies[i]->y = (float)rand() / RAND_MAX * max_pos;
        bodies[i]->z = (float)rand() / RAND_MAX * max_pos;
        if(i == 0) {
            bodies[i]->vx = 0;
            bodies[i]->vy = 0;
            bodies[i]->vz = 0;
        } else {
            bodies[i]->vx = (float)rand() / RAND_MAX * max_vel;
            bodies[i]->vy = (float)rand() / RAND_MAX * max_vel;
            bodies[i]->vz = (float)rand() / RAND_MAX * max_vel;
        }
    }
    return bodies;
}

Body **calculate_iteration(Body** bodies, int num_bodies) {
    Body** new_bodies = (Body**) malloc(num_bodies * sizeof(Body*));
    if(new_bodies == NULL)
        throw_err(__LINE__, 1);
    #pragma omp parallel for schedule(guided, 50)
    for (int i = 0; i < num_bodies; i++) {
        new_bodies[i] = malloc(sizeof(Body));
        if(new_bodies[i] == NULL)
            throw_err(__LINE__, 1);
        float forceX = 0;
        float forceY = 0;
        float forceZ = 0;
        for (int j = 0; j < num_bodies; j++) {
            float rx = bodies[j]->x - bodies[i]->x;
            float ry = bodies[j]->y - bodies[i]->y;
            float rz = bodies[j]->z - bodies[i]->z;
            float distance = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

            forceX += bodies[j]->mass * rx / (pow(distance, 3) + EPSILON);// G je epsilon xd
            forceY += bodies[j]->mass * ry / (pow(distance, 3) + EPSILON); // G je epsilon xd
            forceZ += bodies[j]->mass * rz / (pow(distance, 3) + EPSILON); // G je epsilon xd
        }      

        forceX *= G * bodies[i]->mass;  
        forceY *= G * bodies[i]->mass;  
        forceZ *= G * bodies[i]->mass;  

        float ax = forceX / bodies[i]->mass;
        float ay = forceY / bodies[i]->mass;
        float az = forceZ / bodies[i]->mass;

        new_bodies[i]->mass = bodies[i]->mass;
        new_bodies[i]->vx = bodies[i]->vx + ax * DT;
        new_bodies[i]->vy = bodies[i]->vy + ay * DT;
        new_bodies[i]->vz = bodies[i]->vz + az * DT;
        new_bodies[i]->x = bodies[i]->x + bodies[i]->vx * DT + 0.5 * ax * pow(DT, 2);
        new_bodies[i]->y = bodies[i]->y + bodies[i]->vy * DT + 0.5 * ay * pow(DT, 2);
        new_bodies[i]->z = bodies[i]->z + bodies[i]->vz * DT + 0.5 * az * pow(DT, 2);

    }

    cleanup(bodies, num_bodies);
    return new_bodies;
}


int main(int argc, char *argv[]) {
    Body **bodies;
    FILE *fp;
    int N_BODIES, N_ITER;

    if(argc != 3) {
        printf("Usage: %s <num_bodies> <num_iterations>\n", argv[0]);
        exit(1);
    }

    #pragma omp parallel
    #pragma omp master

    N_BODIES = atoi(argv[1]);
    N_ITER = atoi(argv[2]);

    fp = fopen(LOG_FILE, "wa");

    if (fp == NULL)
        throw_err(__LINE__, errno);

    bodies = generate_initial_population(N_BODIES, 30, 3, 3);

    gettimeofday(&tv1, NULL);
    for(int i = 0; i < N_ITER; i ++) { 
        bodies = calculate_iteration(bodies, N_BODIES);
    }
    gettimeofday(&tv2, NULL);

    printf ("CPU: %f seconds\n", (float) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (float) (tv2.tv_sec - tv1.tv_sec));

    print_boddies(bodies,N_BODIES, N_ITER, fp);

    cleanup(bodies, N_BODIES);
    fclose(fp);
    return 0;
}
