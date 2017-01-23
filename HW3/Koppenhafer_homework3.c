//*****************************************************************************
//
// Homework 2 - Problem 2 - ECE588
// Written by: Sean Koppenhafer (Koppen2)
// 01/20/2017
//
//*****************************************************************************

#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define GRID_X_SIZE 1000
#define GRID_Y_SIZE 1000

static double heat_grid1[GRID_X_SIZE][GRID_Y_SIZE];
static double heat_grid2[GRID_X_SIZE][GRID_Y_SIZE];

//*****************************************************************************
// Helper Functions
//*****************************************************************************
uint32_t parse_cmdline(int, char**, int, int);
void get_system_time(struct timespec*);
unsigned long long int calc_runtime(struct timespec, struct timespec);
void print_results(double, unsigned long long int);
void print_grid(double grid[GRID_X_SIZE][GRID_Y_SIZE]);


void init_grid(double grid[GRID_X_SIZE][GRID_Y_SIZE]);


bool run_heat_calculations(double grid1[GRID_X_SIZE][GRID_Y_SIZE],
                           double grid2[GRID_X_SIZE][GRID_Y_SIZE],
                           uint16_t);
double get_heat_value(double grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t, uint16_t);
double calc_heat_value(double past_grid[GRID_X_SIZE][GRID_Y_SIZE],
                       uint16_t x, uint16_t y);


//*****************************************************************************
// Main
//*****************************************************************************
int main(int argc, char* argv[]) {
    int args_expected = 2;
    struct timespec start_time, end_time;

    get_system_time(&start_time);
    uint32_t number_of_threads = parse_cmdline(argc, &argv[0], 1, args_expected);

    init_grid(heat_grid1);
    init_grid(heat_grid2);

    bool grid1_last_used = run_heat_calculations(heat_grid1, heat_grid2, 6000);
    print_grid(heat_grid2);

    get_system_time(&end_time);
    unsigned long long int runtime = calc_runtime(start_time, end_time);
    print_results(0.0, runtime);

    return 0;
}


//*****************************************************************************
// Helper Functions
//*****************************************************************************
uint32_t parse_cmdline(int argc, char** argv, int arg_position, int args_expected) {
    uint32_t number_on_cmdline;
    if(argc < args_expected) {
        printf("ERROR: Argc is less than %d for number of threads. Exiting.", args_expected);
        exit(-1);
    }
    if(arg_position >= argc) {
        printf("ERROR: Asked for argument number %d, which is equal to or larger than argc: %d. Exiting.", args_expected, argc);
        exit(-2);
    }
    sscanf(argv[arg_position], "%"SCNu32, &number_on_cmdline);
    return number_on_cmdline;
}


void get_system_time(struct timespec* container) {
    int retval = clock_gettime(CLOCK_REALTIME, container);
    if(retval != 0) {
        printf("ERROR: Failed to get system time.\n");
        exit(-1);
    }
}


unsigned long long int calc_runtime(struct timespec start_time, struct timespec end_time) {
    unsigned long long int runtime = 1000000000 * (end_time.tv_sec - start_time.tv_sec) + end_time.tv_nsec - start_time.tv_nsec;
    return runtime;
}


void print_results(double pi, unsigned long long int runtime) {
    printf("The value of pi is %.22f\n", pi);
    printf("Time = %llu nanoseconds\t(%llu.%09llu sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
}


void print_grid(double grid[GRID_X_SIZE][GRID_Y_SIZE]) {
    for(uint16_t x = 0; x < GRID_X_SIZE; x++) {
        for(uint16_t y = 0; y < GRID_Y_SIZE; y++) {
            printf("X: %u, Y: %u, Value: %f\n", x, y, grid[x][y]);
        }
    }
}


//*****************************************************************************
// Init Functions
//*****************************************************************************
void init_grid(double grid[GRID_X_SIZE][GRID_Y_SIZE]) {
    bool x_in_center = false;
    bool y_in_center = false;

    for(uint16_t x = 0; x < GRID_X_SIZE; x++) {
        for(uint16_t y = 0; y < GRID_Y_SIZE; y++) {
            x_in_center = (x >= 200) && (x <= 800);
            y_in_center = (y >= 200) && (y <= 800);
            if( x_in_center && y_in_center ) grid[x][y] = 500.0;
            else grid[x][y] = 0.0;
        }
    }
}


bool run_heat_calculations(double grid1[GRID_X_SIZE][GRID_Y_SIZE],
                           double grid2[GRID_X_SIZE][GRID_Y_SIZE],
                           uint16_t time_ticks) {
    static bool grid1_older = true;

    for(uint16_t tick = 0; tick < time_ticks; tick++) {
        for(uint16_t x = 1; x < GRID_X_SIZE - 1; x++) {
            for(uint16_t y = 1; y < GRID_Y_SIZE - 1; y++) {
                if(grid1_older) grid2[x][y] = calc_heat_value(grid1, x, y);
                else grid1[x][y] = calc_heat_value(grid2, x, y);
            }
        }

        if(grid1_older) grid1_older = false;
        else grid1_older = true;
    }
    return grid1_older;
}


double get_heat_value(double grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t x, uint16_t y) {
    const double out_of_bounds_heat_value = 0;
    double heat_value;
    bool x_out_of_bounds = (x < 0.0) || (x > GRID_X_SIZE);
    bool y_out_of_bounds = (y < 0.0) || (y > GRID_Y_SIZE);

    if( x_out_of_bounds || y_out_of_bounds) heat_value = out_of_bounds_heat_value;
    else heat_value = grid[x][y];

    return heat_value;
}


double calc_heat_value(double past_grid[GRID_X_SIZE][GRID_Y_SIZE],
                       uint16_t x, uint16_t y) {
    const double Cx = 0.12;
    const double Cy = 0.1;
    double Uxy = get_heat_value(past_grid, x, y);
    double Uxm1y = get_heat_value(past_grid, x - 1, y);
    double Uxp1y = get_heat_value(past_grid, x + 1, y);
    double Uxym1 = get_heat_value(past_grid, x, y - 1);
    double Uxyp1 = get_heat_value(past_grid, x, y + 1);

    double result = Uxy + \
                    Cx * (Uxp1y + Uxm1y - (2.0 * Uxy)) + \
                    Cy * (Uxyp1 + Uxym1 - (2.0 * Uxy));

    return result;
}

