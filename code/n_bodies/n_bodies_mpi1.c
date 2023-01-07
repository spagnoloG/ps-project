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

#define G 1                         // gravitational contant
#define DT 0.001                    // time derivative
#define EPSILON 1                   // epsilon to avoid division by 0
#define LOG_FILE "n_bodies_mpi1.txt" // logfile location

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

typedef struct  {
    float Fx;
    float Fy;
    float Fz;
} Force;

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
        bodies[i].z = (float)rand() /(float) RAND_MAX * max_pos;
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


void compute_forces(Force *forces, Body *new_bodies, int start, int end){
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = start; i < end; i++) {
        forces[i - start].Fx = 0;
        forces[i - start].Fy = 0;
        forces[i - start].Fz = 0;
        for(int j = start; j <  end; j++) {
            float rx = new_bodies[i - start].x - new_bodies[j - start].x;
            float ry = new_bodies[i - start].y - new_bodies[j - start].y;
            float rz = new_bodies[i - start].z - new_bodies[j - start].z; 
            float distance = sqrt(pow(rx,2) + pow(ry,2) + pow(rz,2));
            
            forces[i - start].Fx += new_bodies[i - start].mass * rx / (pow(distance , 3)+ EPSILON);
            forces[i - start].Fy += new_bodies[i - start].mass * ry / (pow(distance , 3)+ EPSILON);
            forces[i - start].Fz += new_bodies[i - start].mass * rz / (pow(distance , 3)+ EPSILON);
        }
    }
}

void compute_partial_iteration(Body *bodies, Force *forces, Body *new_bodies,  int num_bodies, int start, int end){
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = start; i < end; i++) {
        for(int j = 0; j < num_bodies && (j < start || j > end); j++) {
            float rx = bodies[j].x - bodies[i].x;
            float ry = bodies[j].y - bodies[i].y;
            float rz = bodies[j].z - bodies[i].z;
            float distance = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

            forces[i - start].Fx += bodies[j].mass * rx / (pow(distance, 3) + EPSILON);
            forces[i - start].Fy += bodies[j].mass * ry / (pow(distance, 3) + EPSILON);
            forces[i - start].Fz += bodies[j].mass * rz / (pow(distance, 3) + EPSILON);
        }

        forces[i - start].Fx *= bodies[i].mass * G;
        forces[i - start].Fy *= bodies[i].mass * G;
        forces[i - start].Fz *= bodies[i].mass * G;

        float ax = forces[i - start].Fx / bodies[i].mass;
        float ay = forces[i - start].Fy / bodies[i].mass;
        float az = forces[i - start].Fz / bodies[i].mass;

        new_bodies[i - start].mass = bodies[i].mass;
        new_bodies[i - start].x = bodies[i].x + bodies[i].vx * DT + 0.5 * ax * pow(DT, 2);
        new_bodies[i - start].y = bodies[i].y + bodies[i].vy * DT + 0.5 * ay * pow(DT, 2);
        new_bodies[i - start].z = bodies[i].z + bodies[i].vz * DT + 0.5 * az * pow(DT, 2);

        new_bodies[i - start].vx = bodies[i].vx + ax * DT;
        new_bodies[i - start].vy = bodies[i].vy + ay * DT;
        new_bodies[i - start].vz = bodies[i].vz + az * DT;
    }
}

Body *calculate_iteration(Body *bodies, int num_bodies, int start, int end)
{

    Body *new_bodies = (Body *)malloc((end - start) * sizeof(Body));
    if (new_bodies == NULL)
        throw_err(__LINE__, 1);
    int i;
    #pragma omp parallel for schedule(dynamic, 1)
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
    int myid, size;
    MPI_Request request;

    printf("start init\n");
    fflush(stdout);

    MPI_Init(&argc, &argv);
    printf("init done\n");
    fflush(stdout);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("init MPI ok\n");
    fflush(stdout);

    if (argc != 3)
    {
        printf("Usage: %s <num_bodies> <num_iterations>\n", argv[0]);
        exit(1);
    }

    N_BODIES = atoi(argv[1]);
    N_ITER = atoi(argv[2]);

    printf("args parsed\n");
    fflush(stdout);

    int *displacement = (int *)malloc(sizeof(int) * size);
    int *displacementByte = (int *)malloc(sizeof(int) * size);
    int *count = (int *)malloc(sizeof(int) * size);
    int *countByte = (int *)malloc(sizeof(int) * size);
    Body *bodies = (Body *)malloc(sizeof(Body) * N_BODIES);

    printf("end init\n");
    fflush(stdout);
    count_displacement(count, countByte, displacement, displacementByte, size, N_BODIES);
    printf("counted\n");
    fflush(stdout);

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
    // MPI DATATYPE

    MPI_Bcast(bodies, N_BODIES, mpi_body_type, 0, MPI_COMM_WORLD);
 
    gettimeofday(&tv1, NULL);
    Body *new_bodies;
    Force *my_forces;
    printf("h1\n");
    fflush(stdout);
    #pragma omp parallel
    #pragma omp master
    for (int i = 0; i < N_ITER; i++)
    {
        //if (i % 1000 == 0 && myid == 0)
        //{
        //        print_boddies(bodies, N_BODIES, i, fp);
        //}

        if( i == 0){
            printf("h2\n");
            fflush(stdout);
            new_bodies = calculate_iteration(bodies, N_BODIES, displacement[myid], displacement[myid] + count[myid]);
            my_forces = malloc(sizeof(Force) * count[myid]);
            print_boddies(new_bodies, count[myid], i, fp);
            printf("h3\n");
            fflush(stdout);
        }
        else { // compute partial results out of new_bodies
            printf("h4\n");
            fflush(stdout);
            compute_forces(my_forces, new_bodies, count[myid], N_BODIES);
            printf("h5\n");
            fflush(stdout);
            // test if  allgather request is done
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            printf("h6\n");
            fflush(stdout);
            // compute new bodies
            compute_partial_iteration(bodies, my_forces, new_bodies, N_BODIES, displacement[myid], displacement[myid] + count[myid]);
            printf("h7\n");
            fflush(stdout);
        }
        printf("beforegatherv\n");
        fflush(stdout);
        // print count
        for(int i = 0; i < size; i++){
            printf("count[%d] = %d\n", i, count[i]);
        }
        // print displacement 
        for(int i = 0; i < size; i++){
            printf("displacement[%d] = %d\n", i, displacement[i]);
        }

        MPI_Allgatherv(new_bodies, count[myid], mpi_body_type, bodies, count, displacement, mpi_body_type, MPI_COMM_WORLD);
        printf("haftergatherv\n");
        fflush(stdout);
    }
    gettimeofday(&tv2, NULL);

    printf("CPU: %f seconds\n", (float)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
                                    (float)(tv2.tv_sec - tv1.tv_sec));

    fclose(fp);

    MPI_Finalize();
    return 0;
}
