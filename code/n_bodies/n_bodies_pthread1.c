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
#include <pthread.h>

#define G 1                                 // gravitational contant
#define DT 0.001                            // time derivative
#define EPSILON 1                           // epsilon to avoid division by 0
#define LOG_FILE "output_pthread1.txt"       // logfile location

struct timeval  tv1, tv2;
struct winsize w;

typedef struct {
    double mass;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
} Body;

typedef struct {
    Body *bodies;
    Body *new_bodies;
    int n_bodies;
    int n_iter;
    int n_threads;
    int start;
    int end;
    int tid;
} ThreadData;

ThreadData *thread_data_g;
pthread_barrier_t barrier;
FILE *fp;

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

void print_boddies(Body *bodies, int n, int iteration, FILE *logfile) {
    for(int i = 0; i < n; i++) {
        fprintf(logfile,"%d,%f,%f,%f,%f,%f,%f,%f\n", iteration, bodies[i].mass, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].vx, bodies[i].vy, bodies[i].vz);
        //printf("%d,%f,%f,%f,%f,%f,%f,%f\n", iteration, bodies[i].mass, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].vx, bodies[i].vy, bodies[i].vz);
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
Body *generate_initial_population(int n, double max_mass, double max_pos, double max_vel) {
    Body *bodies = malloc(sizeof(Body) * n);
    if(bodies == NULL)
        throw_err(__LINE__, 1);
    
    srand(42);

    for(int i = 0; i < n; i++) {

        bodies[i].mass = (double)rand() / RAND_MAX * max_mass;
        bodies[i].x = (double)rand() / RAND_MAX * max_pos;
        bodies[i].y = (double)rand() / RAND_MAX * max_pos;
        bodies[i].z = (double)rand() / RAND_MAX * max_pos;
        if(i == 0) {
            bodies[i].vx = 0;
            bodies[i].vy = 0;
            bodies[i].vz = 0;
        } else {
            bodies[i].vx = (double)rand() / RAND_MAX * max_vel;
            bodies[i].vy = (double)rand() / RAND_MAX * max_vel;
            bodies[i].vz = (double)rand() / RAND_MAX * max_vel;
        }
    }
    return bodies;
}

/*
 * @param bodies: array of all bodies
 * @param num_boides: number of all bodies
 * @param start: start index of bodies array
 * @param end: end index of bodies array
 */
Body *calculate_iteration(Body* bodies, int num_bodies, int start, int end) {
    int n_bodies_to_process = end - start;
    Body* new_bodies = (Body*) malloc((n_bodies_to_process) * sizeof(Body));
    if(new_bodies == NULL)
        throw_err(__LINE__, 1);
    for (int i = start; i < end; i++) {
        double forceX = 0;
        double forceY = 0;
        double forceZ = 0;
        for (int j = 0; j < num_bodies; j++) {
            double rx = bodies[j].x - bodies[i].x;
            double ry = bodies[j].y - bodies[i].y;
            double rz = bodies[j].z - bodies[i].z;
            double distance = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

            forceX += bodies[j].mass * rx / (pow(distance, 3) + EPSILON);// G je epsilon xd
            forceY += bodies[j].mass * ry / (pow(distance, 3) + EPSILON); // G je epsilon xd
            forceZ += bodies[j].mass * rz / (pow(distance, 3) + EPSILON); // G je epsilon xd
        }      

        forceX *= G * bodies[i].mass;  
        forceY *= G * bodies[i].mass;  
        forceZ *= G * bodies[i].mass;  

        double ax = forceX / bodies[i].mass;
        double ay = forceY / bodies[i].mass;
        double az = forceZ / bodies[i].mass;

        new_bodies[i-start].mass = bodies[i].mass;
        new_bodies[i-start].x = bodies[i].x + bodies[i].vx * DT + 0.5 * ax * pow(DT, 2);
        new_bodies[i-start].y = bodies[i].y + bodies[i].vy * DT + 0.5 * ay * pow(DT, 2);
        new_bodies[i-start].z = bodies[i].z + bodies[i].vz * DT + 0.5 * az * pow(DT, 2);

        new_bodies[i-start].vx = bodies[i].vx + ax * DT;
        new_bodies[i-start].vy = bodies[i].vy + ay * DT;
        new_bodies[i-start].vz = bodies[i].vz + az * DT;
    }
    return new_bodies;
}

void *parallel_iteration(void *arg) {
    ThreadData *thread_descriptor = (ThreadData *) arg;
    ThreadData master_descriptor = thread_data_g[0];

    fflush(stdout);

    for(int i = 0; i < thread_descriptor->n_iter; i++) { 
        thread_descriptor->new_bodies = calculate_iteration(thread_descriptor->bodies, thread_descriptor->n_bodies, thread_descriptor->start, thread_descriptor->end); 
        
        // copy all bodies from new_bodies to bodies of the master thread
        memcpy(&master_descriptor.bodies[thread_descriptor->start], thread_descriptor->new_bodies, (thread_descriptor->end - thread_descriptor->start) * sizeof(Body));

        free(thread_descriptor->new_bodies);
        
        thread_descriptor->bodies = master_descriptor.bodies;

        pthread_barrier_wait(&barrier);

    }


    return NULL;
}


int main(int argc, char *argv[]) {
    Body *bodies;
    int N_BODIES, N_ITER, N_THREADS;

    if(argc != 4) {
        printf("Usage: %s <num_bodies> <num_iterations> <n_threads>\n", argv[0]);
        exit(1);
    }

    N_BODIES = atoi(argv[1]);
    N_ITER = atoi(argv[2]);
    N_THREADS = atoi(argv[3]);
    pthread_t threads[N_THREADS];
        
    // malloc insteod on stac
    thread_data_g = (ThreadData *)malloc(sizeof(ThreadData) * N_THREADS);

    // init barrier
    pthread_barrier_init(&barrier, NULL, N_THREADS);

    fp = fopen(LOG_FILE, "wa");

    if (fp == NULL)
        throw_err(__LINE__, errno);

    bodies = generate_initial_population(N_BODIES, 30, 3, 3);
    
    // initialize thread data 
    for(int i = 0; i < N_THREADS; i++) {
        thread_data_g[i].tid = i;
        thread_data_g[i].bodies = bodies;
        thread_data_g[i].n_bodies = N_BODIES;
        thread_data_g[i].start = i * N_BODIES / N_THREADS;
        thread_data_g[i].end = (i + 1) * N_BODIES / N_THREADS;
        thread_data_g[i].n_iter = N_ITER;
        thread_data_g[i].n_threads = N_THREADS;
    }

    fflush(stdout);

    gettimeofday(&tv1, NULL);

    // Compute iteration
    for(int j = 0; j < N_THREADS; j++)
        if(pthread_create(&threads[j], NULL, parallel_iteration, (void *) &thread_data_g[j]) != 0)
            throw_err(__LINE__, errno);

    fflush(stdout);

    // Join threads
    for(int j = 0; j < N_THREADS; j++)
        if(pthread_join(threads[j], NULL) != 0)
            throw_err(__LINE__, errno);

    fflush(stdout);

    gettimeofday(&tv2, NULL);

    printf ("CPU: %f seconds\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));

    print_boddies(bodies,N_BODIES, N_ITER, fp);

    fclose(fp);
    return 0;
}
