#include "mpi.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>


#define GRID_X_SIZE 1000
#define GRID_Y_SIZE 1000
#define MASTER 0           /* task ID of master task */
#define NUM_TICKS 6000

static float* heat_grid1 = NULL;


float* allocate_grid(uint32_t, uint32_t);
void print_grid(float*, uint32_t, uint32_t);
void init_grid(float*);

void print_interval_step(float grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t);
float get_heat_value(float*, int, int);
float calc_heat_value(float*, int, int);

int main (int argc, char *argv[])
{
    struct timespec StartTime, EndTime;
    int ret;
    uint16_t x_per_thread;
    double	localsum, sum; 
    float* past_grid;
    float* next_grid;
    int	taskid;
    int numtasks;
    int return_code;
    long start_x;
    long end_x;
    bool is_last_thread;
    uint16_t ticks;
    uint16_t x;
    uint16_t y;
    uint32_t grid_pos;

    MPI_Status status;

    heat_grid1 = allocate_grid(GRID_X_SIZE, GRID_Y_SIZE);
    init_grid(heat_grid1);

    /* Obtain number of tasks and task ID */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    printf ("MPI task %d has started...\n", taskid);

    if (taskid == MASTER) {
        ret = clock_gettime(CLOCK_REALTIME, &StartTime);
        assert(ret == 0);
    }

    x_per_thread = GRID_X_SIZE / numtasks;
    start_x = x_per_thread * taskid;

    is_last_thread = taskid == (numtasks-1);
    if(is_last_thread) end_x = GRID_X_SIZE;
    else end_x = x_per_thread * (taskid + 1);

    past_grid = allocate_grid(x_per_thread, GRID_Y_SIZE);
    next_grid = allocate_grid(x_per_thread, GRID_Y_SIZE);

    //Magic happens
    MPI_Scatter(heat_grid1, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                past_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                MASTER, MPI_COMM_WORLD);

    printf("Thread id: %u, start: %u, end: %u\n", taskid, start_x, end_x);

    for(x = start_x; x < end_x; x++) {
        printf("%u, X: %u\n", taskid, x);
        for(y = 0; y < GRID_Y_SIZE; y++) {
            grid_pos = (x * GRID_Y_SIZE) + y;
            next_grid[grid_pos] = calc_heat_value(past_grid, (int)x, (int)y);
        }
    }
    printf("GOT OUT\n");

    MPI_Gather(next_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
               heat_grid1, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
               MASTER, MPI_COMM_WORLD);

    print_grid(heat_grid1, GRID_X_SIZE, GRID_Y_SIZE);


    /* Master prints sum */
    if (taskid == MASTER) { 
        ret = clock_gettime(CLOCK_REALTIME, &EndTime);
        assert(ret == 0);
        unsigned long long int runtime = 1000000000 * (EndTime.tv_sec - StartTime.tv_sec) + EndTime.tv_nsec - StartTime.tv_nsec; 
        printf("\nTime = %lld nanoseconds\t(%d.%09lld sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
    }     

    MPI_Finalize();

    free(heat_grid1);
    return 0;
}


//void print_interval_step(float grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t cycle) {
//    printf("Cycle: %u. (1,1): %f, (150,150): %f, (400,400): %f, (500,500): %f, (750,750): %f, (900,900): %f\n",
//           cycle,
//           heat_grid1[1][1], heat_grid1[150][150], heat_grid1[400][400], heat_grid1[500][500],
//           heat_grid1[750][750], heat_grid1[900][900]);
//}


void init_grid(float* grid) {
    const uint16_t min_cord = 200;
    const uint16_t max_cord = 800;
    const float center_temp = 500.0;
    const float border_temp = 0.0;
    bool x_in_center = false;
    bool y_in_center = false;
    uint16_t x;
    uint16_t y;
    uint32_t grid_pos;

    for(x = 0; x < GRID_X_SIZE; x++) {
        for(y = 0; y < GRID_Y_SIZE; y++) {
            grid_pos = (x * GRID_X_SIZE) + y;
            x_in_center = (x >= min_cord) && (x <= max_cord);
            y_in_center = (y >= min_cord) && (y <= max_cord);
            if( x_in_center && y_in_center ) grid[grid_pos] = center_temp; 
            else grid[grid_pos] = border_temp;
        }
    }
}


float get_heat_value(float* grid, int x, int y) {
    const float out_of_bounds_heat_value = 0.0;
    float heat_value;
    bool x_out_of_bounds = (x < 0) || (x >= GRID_X_SIZE);
    bool y_out_of_bounds = (y < 0) || (y >= GRID_Y_SIZE);

    if( x_out_of_bounds || y_out_of_bounds) heat_value = out_of_bounds_heat_value;
    else heat_value = grid[ (x * GRID_Y_SIZE) + y];

    return heat_value;
}


float calc_heat_value(float* past_grid, int x, int y) {
    const float Cx = 0.12;
    const float Cy = 0.1;
    float Uxy = get_heat_value(past_grid, x, y);
    float Uxm1y = get_heat_value(past_grid, x - 1, y);
    float Uxp1y = get_heat_value(past_grid, x + 1, y);
    float Uxym1 = get_heat_value(past_grid, x, y - 1);
    float Uxyp1 = get_heat_value(past_grid, x, y + 1);

    float result = Uxy + \
                    Cx * (Uxp1y + Uxm1y - (2.0 * Uxy)) + \
                    Cy * (Uxyp1 + Uxym1 - (2.0 * Uxy));

    return result;
}


void print_grid(float* grid, uint32_t x_size, uint32_t  y_size) {
    uint32_t x, y;
    uint32_t grid_pos;

    for(x = 0; x < x_size; x++) {
        for(y = 0; y < y_size; y++) {
            grid_pos = (x * y_size) + y;
            printf("x: %u, y: %u, %f\n", x, y, grid[grid_pos]);
        }
    }
}


float* allocate_grid(uint32_t x_size, uint32_t y_size) {
    uint32_t i;
    float* grid;

    // Make a large flat array so MPI can send the message
    grid = malloc( sizeof(float) * x_size * y_size);
    if(grid == NULL) {
        printf("Malloc of grid size X: %u and size Y: %u failed.\n", x_size, y_size);
        exit(1);
    }

    for(i = 0; i < (x_size * y_size); i++) {
        grid[i] = 0.0;
    }

    return grid;
}

