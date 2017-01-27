//*****************************************************************************
//
// Homework 3 - ECE588
// Written by: Sean Koppenhafer (Koppen2)
// 01/25/2017
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

typedef struct {
    uint32_t thread_id;
    uint16_t start_x;
    uint16_t end_x;
    uint16_t time_ticks;
} thread_args;

static float heat_grid1[GRID_X_SIZE][GRID_Y_SIZE];
static float heat_grid2[GRID_X_SIZE][GRID_Y_SIZE];
static pthread_t* threads = NULL;
static thread_args* args = NULL;

static pthread_mutex_t cycle_mutex;
static pthread_cond_t enable_all_threads;
static pthread_cond_t monitor_check;

static uint32_t number_of_threads;
static uint32_t threads_done;
static uint32_t all_threads_done_mask;

//*****************************************************************************
// Helper Functions
//*****************************************************************************
uint32_t parse_cmdline(int, char**, int, int);
void get_system_time(struct timespec*);
unsigned long long int calc_runtime(struct timespec, struct timespec);
void print_interval_step(float grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t);


//init
void init_grid(float grid[GRID_X_SIZE][GRID_Y_SIZE]);
void init_mutex_and_cond();
uint32_t gen_thread_mask(uint32_t);
pthread_t* allocate_threads(uint32_t);
thread_args* create_threads(pthread_t*, uint32_t);


void run_threads(pthread_t*, uint32_t);
void* run_heat_calculations(void*);
void* monitor(void*);
float get_heat_value(float grid[GRID_X_SIZE][GRID_Y_SIZE], int, int);
float calc_heat_value(float past_grid[GRID_X_SIZE][GRID_Y_SIZE],
                       int, int);

void cleanup();


//*****************************************************************************
// Main
//*****************************************************************************
int main(int argc, char* argv[]) {
    int args_expected = 2;
    struct timespec start_time, end_time;

    get_system_time(&start_time);
    number_of_threads = parse_cmdline(argc, &argv[0], 1, args_expected);

    init_grid(heat_grid1);
    init_grid(heat_grid2);
    threads = allocate_threads(number_of_threads);
    init_mutex_and_cond();
    all_threads_done_mask = gen_thread_mask(number_of_threads);

    args = create_threads(threads, number_of_threads);
    run_threads(threads, number_of_threads);

    get_system_time(&end_time);
    unsigned long long int runtime = calc_runtime(start_time, end_time);
    printf("Time = %llu nanoseconds\t(%llu.%09llu sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
    cleanup();

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


void print_interval_step(float grid[GRID_X_SIZE][GRID_Y_SIZE], uint16_t cycle) {
    printf("Cycle: %u. (1,1): %f, (150,150): %f, (400,400): %f, (500,500): %f, (750,750): %f, (900,900): %f\n",
           cycle,
           heat_grid2[1][1], heat_grid2[150][150], heat_grid2[400][400], heat_grid2[500][500],
           heat_grid2[750][750], heat_grid2[900][900]);
}


void init_mutex_and_cond() {
    int retval = pthread_mutex_init(&cycle_mutex, NULL);
    if(retval) {
        printf("ERROR: Could not init mutex\n");
        exit(-2);
    }
    retval = pthread_cond_init(&enable_all_threads, NULL);
    if(retval) {
        printf("ERROR: Could not init cond var enable_all_threads\n");
        exit(-2);
    }
    retval = pthread_cond_init(&monitor_check, NULL);
    if(retval) {
        printf("ERROR: Could not init cond var monitor_check\n");
        exit(-2);
    }
}


uint32_t gen_thread_mask(uint32_t number_of_threads) {
    uint32_t thread_mask = 0;

    for(uint32_t i = 0; i < number_of_threads; i++) {
        thread_mask <<= 1;
        thread_mask |= 0x1;
    }
    return thread_mask;
}


//*****************************************************************************
// Init Functions
//*****************************************************************************
void init_grid(float grid[GRID_X_SIZE][GRID_Y_SIZE]) {
    const uint16_t min_cord = 200;
    const uint16_t max_cord = 800;
    const float center_temp = 500.0;
    const float border_temp = 0.0;
    bool x_in_center = false;
    bool y_in_center = false;

    for(uint16_t x = 0; x < GRID_X_SIZE; x++) {
        for(uint16_t y = 0; y < GRID_Y_SIZE; y++) {
            x_in_center = (x >= min_cord) && (x <= max_cord);
            y_in_center = (y >= min_cord) && (y <= max_cord);
            if( x_in_center && y_in_center ) grid[x][y] = center_temp;
            else grid[x][y] = border_temp;
        }
    }
}


pthread_t* allocate_threads(uint32_t number_of_threads) {
    uint32_t threads_and_monitor = number_of_threads + 1;

    pthread_t* threads = malloc( sizeof(pthread_t) * threads_and_monitor );
    if(threads == NULL) {
        printf("ERROR: Threads malloc returned NULL. Attempted to malloc %u threads\n", number_of_threads);
        exit(-1);
    }
    return threads;
}


thread_args* create_threads(pthread_t* threads, uint32_t number_of_threads) {
    thread_args* args = malloc( sizeof(thread_args) * number_of_threads);
    if(args == NULL) {
        printf("ERROR: Thread args malloc failed. Tried to malloc %u\n", number_of_threads);
        exit(-1);
    }

    const uint16_t time_ticks = 6000;
    uint16_t x_per_thread = GRID_X_SIZE / number_of_threads;
    uint32_t last_thread = number_of_threads - 1;
    for(uint32_t i = 0; i < number_of_threads; i++) {
        args[i].thread_id = i;
        args[i].start_x = x_per_thread * i;
        if(i == last_thread) args[i].end_x = GRID_X_SIZE;
        else args[i].end_x = x_per_thread * (i + 1);
        args[i].time_ticks = time_ticks;
        pthread_create(&threads[i], NULL, &run_heat_calculations, (void*)&args[i]);
    }

    // Last thread is always monitor
    pthread_create(&threads[number_of_threads], NULL, &monitor, (void*)time_ticks);

    return args;
}


void run_threads(pthread_t* threads, uint32_t number_of_threads) {
    uint32_t threads_and_monitor = number_of_threads + 1;

    for(uint32_t i = 0; i < threads_and_monitor; i++) {
        pthread_join(threads[i], NULL);
    }
}


void* run_heat_calculations(void* void_args) {
    thread_args* args = (thread_args*)void_args;
    bool grid1_is_older = 1;
    uint32_t thread_flag_bit = 1 << args->thread_id;

    for(uint16_t tick = 0; tick <= args->time_ticks; tick++) {
        for(uint16_t x = args->start_x; x < args->end_x; x++) {
            for(uint16_t y = 0; y < GRID_Y_SIZE; y++) {
                if(grid1_is_older) heat_grid2[x][y] = calc_heat_value(heat_grid1, (int)x, (int)y);
                else heat_grid1[x][y] = calc_heat_value(heat_grid2, (int)x, (int)y);
            }
        }
        grid1_is_older ^= 1;

        pthread_mutex_lock(&cycle_mutex);
            threads_done |= thread_flag_bit;
            pthread_cond_signal(&monitor_check);

            while(threads_done & thread_flag_bit) {
                pthread_cond_wait(&enable_all_threads, &cycle_mutex);
            }
        pthread_mutex_unlock(&cycle_mutex);
    }
    return 0;
}


void* monitor(void* args) {
    const uint16_t cycle_interval = 200;
    uint16_t iterations = (uint16_t)args;
    bool grid1_is_older = 1;

    for(uint16_t cycle = 0; cycle <= iterations; cycle++) {
        if( (cycle % cycle_interval) == 0 ) {
            if(grid1_is_older) print_interval_step(heat_grid2, cycle);
            else print_interval_step(heat_grid1, cycle);
        }

        pthread_mutex_lock(&cycle_mutex);
            while(threads_done != all_threads_done_mask) {
                pthread_cond_wait(&monitor_check, &cycle_mutex);
            }
            threads_done = 0;
        pthread_mutex_unlock(&cycle_mutex);
        pthread_cond_broadcast(&enable_all_threads);

        grid1_is_older ^= 1;
    }

    return 0;
}


float get_heat_value(float grid[GRID_X_SIZE][GRID_Y_SIZE], int x, int y) {
    const float out_of_bounds_heat_value = 0.0;
    float heat_value;
    bool x_out_of_bounds = (x < 0) || (x >= GRID_X_SIZE);
    bool y_out_of_bounds = (y < 0) || (y >= GRID_Y_SIZE);

    if( x_out_of_bounds || y_out_of_bounds) heat_value = out_of_bounds_heat_value;
    else heat_value = grid[x][y];

    return heat_value;
}


float calc_heat_value(float past_grid[GRID_X_SIZE][GRID_Y_SIZE],
                       int x, int y) {
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

void cleanup() {
    int retval;
    retval = pthread_mutex_destroy(&cycle_mutex);
    if(retval) {
        printf("ERROR: Could not destroy mutex\n");
        exit(-2);
    }
    retval = pthread_cond_destroy(&enable_all_threads);
    if(retval) {
        printf("ERROR: Could not destroy cond var enable_all_threads\n");
        exit(-2);
    }
    retval = pthread_cond_destroy(&monitor_check);
    if(retval) {
        printf("ERROR: Could not destroy cond var monitor_check\n");
        exit(-2);
    }

    free(threads);
    free(args);
}

