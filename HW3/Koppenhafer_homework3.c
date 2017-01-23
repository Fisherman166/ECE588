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



int main(int argc, char* argv[]) {
    int args_expected = 2;
    struct timespec start_time, end_time;

    get_system_time(&start_time);
    uint32_t number_of_threads = parse_cmdline(argc, &argv[0], 1, args_expected);

    init_grid(heat_grid1);
    init_grid(heat_grid2);
    print_grid(heat_grid1);

    get_system_time(&end_time);
    unsigned long long int runtime = calc_runtime(start_time, end_time);
    //print_results(pi, runtime);

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


