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
static int numtasks;


float* allocate_grid(uint32_t, uint32_t);
void print_grid(float*, uint32_t, uint32_t);
void init_grid(float*);
static bool is_first_task(int);
static bool is_last_task(int, int);
void print_interval_step(float grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t);
uint32_t calc_grid_position(uint16_t, uint16_t);

void scatter_and_send(float*, uint16_t, int, int);
float get_heat_value(float*, int, int, int);
float calc_heat_value(float*, int, int, int);

int main (int argc, char *argv[])
{
    struct timespec StartTime, EndTime;
    int	taskid;
    int ret;
    uint16_t x_per_thread;
    float* past_grid;
    float* next_grid;
    long start_x;
    long end_x;
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

    if( is_last_task(taskid, numtasks) ) end_x = GRID_X_SIZE;
    else end_x = x_per_thread * (taskid + 1);

    past_grid = allocate_grid(x_per_thread + 2, GRID_Y_SIZE);
    next_grid = allocate_grid(x_per_thread + 2, GRID_Y_SIZE);

    scatter_and_send(past_grid, x_per_thread, taskid, numtasks);

    if(taskid == 0) {
        print_grid(past_grid, 506, GRID_Y_SIZE);
        exit(0);
    }
    else {
        while(1);
    }
    
    for(x = 0; x < x_per_thread; x++) {
        for(y = 0; y < GRID_Y_SIZE; y++) {
            grid_pos = calc_grid_position(x, y);
            next_grid[grid_pos] = calc_heat_value(past_grid, (int)x, (int)y, taskid);
        }
    }

    MPI_Gather(next_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
               heat_grid1, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
               MASTER, MPI_COMM_WORLD);

    if(taskid == 0) {
        print_grid(heat_grid1, GRID_X_SIZE, GRID_Y_SIZE);
        exit(0);
    }
    else {
        while(1);
    }


    /* Master prints sum */
    if (taskid == MASTER) { 
        ret = clock_gettime(CLOCK_REALTIME, &EndTime);
        assert(ret == 0);
        unsigned long long int runtime = 1000000000 * (EndTime.tv_sec - StartTime.tv_sec) + EndTime.tv_nsec - StartTime.tv_nsec; 
        printf("\nTime = %lld nanoseconds\t(%d.%09lld sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
    }     

    MPI_Finalize();

    free(heat_grid1);
    free(past_grid);
    free(next_grid);
    return 0;
}


//void print_interval_step(float grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t cycle) {
//    printf("Cycle: %u. (1,1): %f, (150,150): %f, (400,400): %f, (500,500): %f, (750,750): %f, (900,900): %f\n",
//           cycle,
//           heat_grid1[1][1], heat_grid1[150][150], heat_grid1[400][400], heat_grid1[500][500],
//           heat_grid1[750][750], heat_grid1[900][900]);
//}



static bool is_last_task(int taskid, int numtasks) {
    int last_task = numtasks - 1;
    if(taskid == last_task) return 1;
    return 0;
}


static bool is_first_task(int taskid) {
    if(taskid == 0) return 1;
    return 0;
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


float get_heat_value(float* grid, int x, int y, int taskid) {
    const float out_of_bounds_heat_value = 0.0;
    const uint16_t low_x_position = 500;
    const uint16_t high_x_position = 501;
    float heat_value;
    uint32_t grid_pos = 0;

    bool x_out_of_bounds = ((x < 1) && (taskid == 0)) || 
                           ((x > 499) && is_last_task(taskid, numtasks));
    bool y_out_of_bounds = (y < 0) || (y >= GRID_Y_SIZE);
    bool is_low_x_border = (x == -1) && (taskid != 0);
    bool is_high_x_border = (x == 500) && (taskid != 1);

    if( x_out_of_bounds || y_out_of_bounds) heat_value = out_of_bounds_heat_value;
    else {
        if( is_low_x_border ) grid_pos = calc_grid_position(low_x_position, y);
        else if( is_high_x_border ) grid_pos = calc_grid_position(high_x_position, y);
        else grid_pos = calc_grid_position(x, y);
        heat_value = grid[grid_pos];
    }

    if( is_low_x_border) printf("taskid %d, %d, %d, %f\n", taskid, x, y, heat_value);

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


void scatter_and_send(float* past_grid, uint16_t x_per_thread, int taskid, int num_tasks) {
    const int send_low_column = 0;
    const int send_high_column = 1;
    float* lower_column;
    float* upper_column;
    int next_taskid;
    int out;
    uint32_t grid_pos = 0;

    MPI_Scatter(heat_grid1, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                past_grid, x_per_thread * GRID_Y_SIZE, MPI_FLOAT,
                MASTER, MPI_COMM_WORLD);

    next_taskid = ((taskid + num_tasks) - 1) % num_tasks; 
    if( is_first_task(taskid) ) {
        lower_column = &(heat_grid1[0]);
        grid_pos = calc_grid_position(x_per_thread - 1, 0);
        upper_column = &(heat_grid1[grid_pos]);
    }
    else if( is_last_task(taskid, num_tasks) ) {
        grid_pos = calc_grid_position(x_per_thread, 0);
        lower_column = &(heat_grid1[grid_pos]);

        grid_pos = calc_grid_position( (x_per_thread * taskid) - 1, 0);
        upper_column = &(heat_grid1[grid_pos]);
    }

    if(taskid != 0) {
        MPI_Send(lower_column, GRID_Y_SIZE, MPI_FLOAT, next_taskid, send_low_column, MPI_COMM_WORLD);
    }
    else {
        grid_pos = calc_grid_position(x_per_thread, 0);
        MPI_Recv(&(past_grid[grid_pos]), GRID_Y_SIZE, MPI_FLOAT, next_taskid, send_low_column, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(taskid == 0) {
        MPI_Send(upper_column, GRID_Y_SIZE, MPI_FLOAT, next_taskid, send_high_column, MPI_COMM_WORLD);
    }
    else {
        grid_pos = calc_grid_position(x_per_thread + 1, 0);
        MPI_Recv(&(past_grid[grid_pos]), GRID_Y_SIZE, MPI_FLOAT, next_taskid, send_high_column, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

uint32_t calc_grid_position(uint16_t x, uint16_t y) {
    uint32_t grid_position = (x * GRID_Y_SIZE) + y;
    return grid_position;
}

