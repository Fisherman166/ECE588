//*****************************************************************************
//
// Homework 4 - ECE588
// Written by: Sean Koppenhafer (Koppen2)
// 02/03/2017
//
//*****************************************************************************

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
#define NUM_CYCLES 6000


static float* heat_grid = NULL;
static int numtasks;
static uint16_t x_per_thread;


float* allocate_grid(uint32_t, uint32_t);
void init_grid(float*);
void calc_per_thread_values(uint16_t*, uint16_t*, uint16_t*, int, int);


static bool is_first_task(int);
static bool is_last_task(int, int);
uint32_t calc_grid_position(uint16_t, uint16_t);
uint16_t start_and_end_differece(uint16_t, uint16_t);


void scatter_and_send(float*, uint16_t, int, int);
float get_heat_value(float*, int, int, int);
float calc_heat_value(float*, int, int, int);


void print_grid(float*, uint32_t, uint32_t, int);
void print_interval_step(float*, uint16_t);


//*****************************************************************************
// Functions
//*****************************************************************************
int main (int argc, char *argv[])
{
    struct timespec StartTime, EndTime;
    int	taskid;
    int ret;
    float* past_grid;
    float* next_grid;
    uint16_t start_x;
    uint16_t end_x;
    uint16_t x;
    uint16_t y;
    uint32_t grid_pos;
    uint16_t cycle;

    MPI_Status status;

    heat_grid = allocate_grid(GRID_X_SIZE, GRID_Y_SIZE);
    init_grid(heat_grid);

    /* Obtain number of tasks and task ID */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

    if (taskid == MASTER) {
        ret = clock_gettime(CLOCK_REALTIME, &StartTime);
        assert(ret == 0);
    }

    calc_per_thread_values(&x_per_thread, &start_x, &end_x, taskid, numtasks);

    // Allocate 2 extra columns for the lower and higher columns that need to be passed from
    // the lower and higher task
    past_grid = allocate_grid(x_per_thread + 2, GRID_Y_SIZE);
    next_grid = allocate_grid(x_per_thread + 2, GRID_Y_SIZE);

    for(cycle = 0; cycle <= NUM_CYCLES; cycle++) {
        if(taskid == MASTER) {
            if( (cycle % 200) == 0) print_interval_step(heat_grid, cycle);
        }

        scatter_and_send(past_grid, x_per_thread, taskid, numtasks);

        for(x = 0; x < x_per_thread; x++) {
            for(y = 0; y < GRID_Y_SIZE; y++) {
                grid_pos = calc_grid_position(x, y);
                next_grid[grid_pos] = calc_heat_value(past_grid, (int)x, (int)y, taskid);
            }
        }

        MPI_Gather(next_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                   heat_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                   MASTER, MPI_COMM_WORLD);
    }

    /* Master prints sum */
    if (taskid == MASTER) {
        ret = clock_gettime(CLOCK_REALTIME, &EndTime);
        assert(ret == 0);
        unsigned long long int runtime = 1000000000 * (EndTime.tv_sec - StartTime.tv_sec) + EndTime.tv_nsec - StartTime.tv_nsec;
        printf("\nTime = %lld nanoseconds\t(%d.%09lld sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
    }

    MPI_Finalize();

    free(heat_grid);
    free(past_grid);
    free(next_grid);
    return 0;
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
            grid_pos = calc_grid_position(x, y);
            x_in_center = (x >= min_cord) && (x <= max_cord);
            y_in_center = (y >= min_cord) && (y <= max_cord);
            if( x_in_center && y_in_center ) grid[grid_pos] = center_temp;
            else grid[grid_pos] = border_temp;
        }
    }
}


void calc_per_thread_values(uint16_t* x_per_thread, uint16_t* start_x, uint16_t* end_x,
                            int taskid, int numtasks) {
    *x_per_thread = GRID_X_SIZE / numtasks;
    *start_x = *x_per_thread * taskid;

    if( is_last_task(taskid, numtasks) ) *end_x = GRID_X_SIZE;
    else *end_x = *x_per_thread * (taskid + 1);
}


static bool is_first_task(int taskid) {
    if(taskid == 0) return 1;
    return 0;
}


static bool is_last_task(int taskid, int numtasks) {
    int last_task = numtasks - 1;
    if(taskid == last_task) return 1;
    return 0;
}


uint32_t calc_grid_position(uint16_t x, uint16_t y) {
    uint32_t grid_position = (x * GRID_Y_SIZE) + y;
    return grid_position;
}


uint16_t start_and_end_differece(uint16_t start_x, uint16_t end_x) {
    uint16_t difference = end_x - start_x;
    return difference;
}


void scatter_and_send(float* past_grid, uint16_t x_per_thread, int taskid, int num_tasks) {
    const int send_low_column = 0;
    const int send_high_column = 1;
    float* lower_column;
    float* upper_column;
    int next_taskid, last_taskid;
    uint32_t grid_pos = 0;
    uint16_t lower_column_x = x_per_thread + 1; // Lower column of current task is upper column of previous task
    uint16_t upper_column_x = x_per_thread; // Upper column of current task is lower column of previous task


    MPI_Scatter(heat_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                past_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                MASTER, MPI_COMM_WORLD);

    // For a single thread this is not required
    if(num_tasks > 1) {
        // First column in past grid is always the lower column to pass.
        grid_pos = calc_grid_position(0, 0);
        lower_column = &(past_grid[grid_pos]);

        // Last column in past grid is always the upper column to pass.
        grid_pos = calc_grid_position(x_per_thread - 1, 0);
        upper_column = &(past_grid[grid_pos]);

        if( is_first_task(taskid) ) last_taskid = numtasks - 1;
        else last_taskid = taskid - 1;
        next_taskid = ((taskid + num_tasks) + 1) % num_tasks;

        MPI_Send(lower_column, GRID_Y_SIZE, MPI_FLOAT, last_taskid, send_low_column, MPI_COMM_WORLD);
        grid_pos = calc_grid_position(lower_column_x, 0);
        MPI_Recv(&(past_grid[grid_pos]), GRID_Y_SIZE, MPI_FLOAT, next_taskid, send_low_column, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(upper_column, GRID_Y_SIZE, MPI_FLOAT, next_taskid, send_high_column, MPI_COMM_WORLD);
        grid_pos = calc_grid_position(upper_column_x, 0);
        MPI_Recv(&(past_grid[grid_pos]), GRID_Y_SIZE, MPI_FLOAT, last_taskid, send_high_column, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}


float get_heat_value(float* grid, int x, int y, int taskid) {
    const float out_of_bounds_heat_value = 0.0;
    const uint16_t low_x_position = x_per_thread;
    const uint16_t high_x_position = x_per_thread + 1;
    float heat_value;
    uint32_t grid_pos = 0;

    bool x_out_of_bounds = ((x < 1) && (taskid == 0)) ||
                           ((x >= x_per_thread) && is_last_task(taskid, numtasks));
    bool y_out_of_bounds = (y < 0) || (y >= GRID_Y_SIZE);
    bool is_low_x_border = (x == -1) && (taskid != 0);
    bool is_high_x_border = (x == x_per_thread) && !(is_last_task(taskid, numtasks));

    if( x_out_of_bounds || y_out_of_bounds) heat_value = out_of_bounds_heat_value;
    else {
        if( is_low_x_border ) grid_pos = calc_grid_position(low_x_position, y);
        else if( is_high_x_border ) grid_pos = calc_grid_position(high_x_position, y);
        else grid_pos = calc_grid_position(x, y);
        heat_value = grid[grid_pos];
    }

    return heat_value;
}


float calc_heat_value(float* past_grid, int x, int y, int taskid) {
    const float Cx = 0.12;
    const float Cy = 0.1;
    float Uxy = get_heat_value(past_grid, x, y, taskid);
    float Uxm1y = get_heat_value(past_grid, x - 1, y, taskid);
    float Uxp1y = get_heat_value(past_grid, x + 1, y, taskid);
    float Uxym1 = get_heat_value(past_grid, x, y - 1, taskid);
    float Uxyp1 = get_heat_value(past_grid, x, y + 1, taskid);

    float result = Uxy + \
                    Cx * (Uxp1y + Uxm1y - (2.0 * Uxy)) + \
                    Cy * (Uxyp1 + Uxym1 - (2.0 * Uxy));

    return result;
}


void print_grid(float* grid, uint32_t x_size, uint32_t  y_size, int taskid) {
    uint32_t x, y;
    uint32_t grid_pos;

    for(x = 0; x < x_size; x++) {
        for(y = 0; y < y_size; y++) {
            grid_pos = (x * y_size) + y;
            printf("taskid: %d, x: %u, y: %u, %f\n", taskid, x, y, grid[grid_pos]);
        }
    }
}


void print_interval_step(float* grid, uint16_t cycle) {
    const uint16_t pos1 = 1;
    const uint16_t pos150 = 150;
    const uint16_t pos400 = 400;
    const uint16_t pos500 = 500;
    const uint16_t pos750 = 750;
    const uint16_t pos900 = 900;
    uint32_t grid_pos1 = calc_grid_position(pos1, pos1);
    uint32_t grid_pos2 = calc_grid_position(pos150, pos150);
    uint32_t grid_pos3 = calc_grid_position(pos400, pos400);
    uint32_t grid_pos4 = calc_grid_position(pos500, pos500);
    uint32_t grid_pos5 = calc_grid_position(pos750, pos750);
    uint32_t grid_pos6 = calc_grid_position(pos900, pos900);

    printf("Cycle: %u. (1,1): %f, (150,150): %f, (400,400): %f, (500,500): %f, (750,750): %f, (900,900): %f\n",
           cycle, grid[grid_pos1], grid[grid_pos2], grid[grid_pos3], grid[grid_pos4],
           grid[grid_pos5], grid[grid_pos6]);
}

