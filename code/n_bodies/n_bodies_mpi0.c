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
#include <omp.h>
#include <mpi.h>
#include <stddef.h>

#define G 1                   // gravitational contant
#define DT 0.001              // time derivative
#define EPSILON 1             // epsilon to avoid division by 0
#define LOG_FILE "n_bodies_mpi" // logfile location

struct timeval tv1, tv2;
struct winsize w;

typedef struct
{
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
void segfault_handler(int sig)
{
    void *array[10];
    size_t size;

    size = backtrace(array, 10);

    fprintf(stderr, "[-] Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

void throw_err(int line, int err_code)
{
    char err_buf[256];

    sprintf(err_buf, "[-] Error: %s -> %d", __FILE__, line);
    perror(err_buf);
    exit(err_code);
}

/*
 * misc functions
 */
void cleanup(Body *bodies, int num_bodies)
{

    free(bodies);
}

void print_boddies(Body *bodies, int n, int iteration, FILE *logfile)
{
    for (int i = 0; i < n; i++)
    {
        fprintf(logfile, "%d,%f,%f,%f,%f,%f,%f,%f\n", iteration, bodies[i].mass, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].vx, bodies[i].vy, bodies[i].vz);
        printf("%d,%f,%f,%f,%f,%f,%f,%f\n", iteration, bodies[i].mass, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].vx, bodies[i].vy, bodies[i].vz);
    }
    fflush(stdout);
}

void progress_bar(int current, int max)
{
    int i;
    char pb_buffer[256];
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    int width = w.ws_col;
    float progress = (float)current / (float)max; // Calculate progress as a fraction
    int chars_printed = (int)(progress * width);  // Calculate number of characters to print
                                                  //
    sprintf(pb_buffer, "[ %d/%d ] ", current, max);
    long pb_buf_len = strlen(pb_buffer);
    for (int i = 0; i < pb_buf_len; i++)
    {
        printf("%c", pb_buffer[i]);
    }

    for (i = 0; i < chars_printed - pb_buf_len; i++)
    {
        printf("=");
    }

    printf("\r");
    fflush(stdout);
}

/*
 * Algorithm
 */

void count_displacement(int *count, int *countByte, int *displacement, int *displacementByte, int size, int NBodies)
{
    int rem = NBodies % size;
    int sum = 0;
    for (int i = 0; i < size; i++)
    {
        count[i] = NBodies / size;
        if (rem > 0)
        {
            count[i]++;
            rem--;
        }
        countByte[i] = count[i] * sizeof(Body);
        displacement[i] = sum;
        displacementByte[i] = displacement[i] * sizeof(Body);
        sum += count[i];
    }
}

Body *generate_initial_population(Body *bodies, int n, float max_mass, float max_pos, float max_vel)
{
    if (bodies == NULL)
        throw_err(__LINE__, 1);

    srand(42);
    for (int i = 0; i < n; i++)
    {

        bodies[i].mass = (float)rand() / (float) RAND_MAX * max_mass;
        bodies[i].x = (float)rand() / (float) RAND_MAX * max_pos;
        bodies[i].y = (float)rand() / (float) RAND_MAX * max_pos;
        bodies[i].z = (float)rand() / (float) RAND_MAX * max_pos;
        if (i == 0)
        {
            bodies[i].vx = 0;
            bodies[i].vy = 0;
            bodies[i].vz = 0;
        }
        else
        {
            bodies[i].vx = (float)rand() / (float) RAND_MAX * max_vel;
            bodies[i].vy = (float)rand() / (float) RAND_MAX * max_vel;
            bodies[i].vz = (float)rand() / (float) RAND_MAX * max_vel;
        }
    }
    return bodies;
}

Body *calculate_iteration(Body *bodies, int num_bodies, int start, int end)
{

    Body *new_bodies = (Body *)malloc((end - start) * sizeof(Body));
    if (new_bodies == NULL)
        throw_err(__LINE__, 1);
    int i;
    #pragma omp parallel for schedule(dynamic, 10)
    for (i = start; i < end; i++)
    {

        float forceX = 0;
        float forceY = 0;
        float forceZ = 0;
        for (int j = 0; j < num_bodies; j++)
        {
            float rx = bodies[j].x - bodies[i].x;
            float ry = bodies[j].y - bodies[i].y;
            float rz = bodies[j].z - bodies[i].z;
            float distance = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

            forceX += bodies[j].mass * rx / (pow(distance, 3) + EPSILON); // G je epsilon xd
            forceY += bodies[j].mass * ry / (pow(distance, 3) + EPSILON); // G je epsilon xd
            forceZ += bodies[j].mass * rz / (pow(distance, 3) + EPSILON); // G je epsilon xd
        }

        forceX *= G * bodies[i].mass;
        forceY *= G * bodies[i].mass;
        forceZ *= G * bodies[i].mass;

        float ax = forceX / bodies[i].mass;
        float ay = forceY / bodies[i].mass;
        float az = forceZ / bodies[i].mass;

        new_bodies[i - start].mass = bodies[i].mass;
        new_bodies[i - start].x = bodies[i].x + bodies[i].vx * DT + 0.5 * ax * pow(DT, 2);
        new_bodies[i - start].y = bodies[i].y + bodies[i].vy * DT + 0.5 * ay * pow(DT, 2);
        new_bodies[i - start].z = bodies[i].z + bodies[i].vz * DT + 0.5 * az * pow(DT, 2);

        new_bodies[i - start].vx = bodies[i].vx + ax * DT;
        new_bodies[i - start].vy = bodies[i].vy + ay * DT;
        new_bodies[i - start].vz = bodies[i].vz + az * DT;
    }
    return new_bodies;
}


int main(int argc, char *argv[])
{
    FILE *fp;
    int N_BODIES, N_ITER;

    if (argc != 3)
    {
        printf("Usage: %s <num_bodies> <num_iterations>\n", argv[0]);
        exit(1);
    }

    #pragma omp parallel
    #pragma omp master

    N_BODIES = atoi(argv[1]);
    N_ITER = atoi(argv[2]);
    int myid, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int *displacement = (int *)malloc(sizeof(int) * size);
    int *displacementByte = (int *)malloc(sizeof(int) * size);
    int *count = (int *)malloc(sizeof(int) * size);
    int *countByte = (int *)malloc(sizeof(int) * size);
    Body *bodies = (Body *)malloc(sizeof(Body) * N_BODIES);
    // float *recv = (float *)malloc(count[myid] * sizeof(float));

    count_displacement(count, countByte, displacement, displacementByte, size, N_BODIES);

    fp = fopen(LOG_FILE, "wa");

    if (fp == NULL)
        throw_err(__LINE__, errno);
    if (myid == 0)
        bodies = generate_initial_population(bodies, N_BODIES, 30, 3, 3);
    
    // MPI DATATYPE
    const int n_items=7;
    int blocklengths[7] = {1,1,1,1,1,1,1};
    MPI_Datatype types[7] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_body_type;
    MPI_Aint offsets[7];
    
    offsets[0] = offsetof(Body, mass);
    offsets[1] = offsetof(Body, x);
    offsets[2] = offsetof(Body, y);
    offsets[3] = offsetof(Body, z);
    offsets[4] = offsetof(Body, vx);
    offsets[5] = offsetof(Body, vy);
    offsets[6] = offsetof(Body, vz);

    MPI_Type_create_struct(n_items, blocklengths, offsets, types, &mpi_body_type);
    MPI_Type_commit(&mpi_body_type);

    
    MPI_Bcast(bodies, N_BODIES, mpi_body_type, 0, MPI_COMM_WORLD);
    
    fflush(fp);
    gettimeofday(&tv1, NULL);
    for (int i = 0; i < N_ITER; i++)
    {
        
        Body *new_bodies = calculate_iteration(bodies, N_BODIES, displacement[myid], displacement[myid] + count[myid]);
        
        MPI_Allgatherv(new_bodies, count[myid], mpi_body_type, bodies, count, displacement, mpi_body_type, MPI_COMM_WORLD);
        free(new_bodies);
    }
    gettimeofday(&tv2, NULL);

    printf("CPU: %f seconds\n", (float)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
                                    (float)(tv2.tv_sec - tv1.tv_sec));

    // cleanup(bodies, N_BODIES);
    fclose(fp);

    MPI_Finalize();
    return 0;
}
