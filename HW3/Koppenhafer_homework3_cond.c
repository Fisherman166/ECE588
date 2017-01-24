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

typedef struct {
    uint32_t thread_id;
    uint16_t start_x;
    uint16_t end_x;
    uint16_t time_ticks;
} thread_args;

static double heat_grid1[GRID_X_SIZE][GRID_Y_SIZE];
static double heat_grid2[GRID_X_SIZE][GRID_Y_SIZE];
static pthread_t* threads = NULL;
static pthread_mutex_t* mutexs = NULL;
static pthread_cond_t* condition_vars = NULL;
static uint16_t* thread_cycles = NULL;
static thread_args* args = NULL;
static uint32_t number_of_threads;

//*****************************************************************************
// Helper Functions
//*****************************************************************************
uint32_t parse_cmdline(int, char**, int, int);
void get_system_time(struct timespec*);
unsigned long long int calc_runtime(struct timespec, struct timespec);
void print_interval_step(double grid[GRID_X_SIZE][GRID_Y_SIZE]);
void print_grid(double grid[GRID_X_SIZE][GRID_Y_SIZE]);


//init
void init_grid(double grid[GRID_X_SIZE][GRID_Y_SIZE]);
pthread_t* allocate_threads(uint32_t);
pthread_mutex_t* allocate_mutexs(uint32_t);
pthread_cond_t* allocate_cond_vars(uint32_t);
uint16_t* allocate_thread_cycles(uint32_t);
void init_mutexs(pthread_mutex_t*, uint32_t);
void init_cond_vars(pthread_cond_t*, uint32_t);
thread_args* create_threads(pthread_t*, uint32_t);


void run_threads(pthread_t*, uint32_t);
void* run_heat_calculations(void*);
double get_heat_value(double grid[GRID_X_SIZE][GRID_Y_SIZE], int, int);
double calc_heat_value(double past_grid[GRID_X_SIZE][GRID_Y_SIZE],
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
    mutexs = allocate_mutexs(number_of_threads);
    condition_vars = allocate_cond_vars(number_of_threads);
    init_mutexs(mutexs, number_of_threads);
    init_cond_vars(condition_vars, number_of_threads);
    thread_cycles = allocate_thread_cycles(number_of_threads);
    args = create_threads(threads, number_of_threads);

    run_threads(threads, number_of_threads);
    //print_grid(heat_grid2);

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


void print_interval_step(double grid[GRID_X_SIZE][GRID_Y_SIZE]) {
    printf("(1,1): %f, (150,150): %f, (400,400): %f, (500,500): %f, (750,750): %f, (900,900): %f\n",
           heat_grid2[1][1], heat_grid2[150][150], heat_grid2[400][400], heat_grid2[500][500],
           heat_grid2[750][750], heat_grid2[900][900]);
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


pthread_t* allocate_threads(uint32_t number_of_threads) {
    pthread_t* threads = malloc( sizeof(pthread_t) * number_of_threads );
    if(threads == NULL) {
        printf("ERROR: Threads malloc returned NULL. Attempted to malloc %u threads\n", number_of_threads);
        exit(-1);
    }
    return threads;
}


pthread_mutex_t* allocate_mutexs(uint32_t number_of_threads) {
    pthread_mutex_t* mutexs = malloc( sizeof(pthread_mutex_t) * number_of_threads );
    if(mutexs == NULL) {
        printf("ERROR: Mutexs malloc returned NULL. Attempted to malloc %u mutexs\n", number_of_threads);
        exit(-1);
    }
    return mutexs;
}


pthread_cond_t* allocate_cond_vars(uint32_t number_of_threads) {
    pthread_cond_t* cond_vars = malloc( sizeof(pthread_cond_t) * number_of_threads );
    if( cond_vars == NULL) {
        printf("ERROR: Cond vars malloc returned NULL. Attempted to malloc %u mutexs\n", number_of_threads);
        exit(-1);
    }
    return cond_vars;
}


uint16_t* allocate_thread_cycles(uint32_t number_of_threads) {
    uint16_t* thread_cycles = malloc( sizeof(uint16_t) * number_of_threads );
    if(thread_cycles == NULL) {
        printf("ERROR: thread_cycles malloc returned NULL. Attempted to malloc %u thread_cycles\n", number_of_threads);
        exit(-1);
    }
    
    // Initialize to 0
    for(uint32_t i = 0; i < number_of_threads; i++) {
        thread_cycles[i] = 0;
    }

    return thread_cycles;
}


void init_mutexs(pthread_mutex_t* mutexs, uint32_t number_of_threads) {
    int retval;
    for(uint32_t i = 0; i < number_of_threads; i++) {
        retval = pthread_mutex_init(&mutexs[i], NULL);
        if(retval) {
            printf("ERROR: Could not init mutex number %u\n", i);
            exit(-2);
        }
    }
}


void init_cond_vars(pthread_cond_t* condition_vars, uint32_t number_of_threads) {
    int retval;
    for(uint32_t i = 0; i < number_of_threads; i++) {
        retval = pthread_cond_init(&condition_vars[i], NULL);
        if(retval) {
            printf("ERROR: Could not init cond var number %u\n", i);
            exit(-2);
        }
    }
}


thread_args* create_threads(pthread_t* threads, uint32_t number_of_threads) {
    uint16_t x_per_thread = GRID_X_SIZE / number_of_threads;
    thread_args* args = malloc( sizeof(thread_args) * number_of_threads);
    if(args == NULL) {
        printf("ERROR: Thread args malloc failed. Tried to malloc %u\n", number_of_threads);
        exit(-1);
    }

    for(uint32_t i = 0; i < number_of_threads; i++) {
        args[i].thread_id = i;
        args[i].start_x = x_per_thread * i;
        args[i].end_x = x_per_thread * (i + 1);
        args[i].time_ticks = 6000;
        pthread_create(&threads[i], NULL, &run_heat_calculations, (void*)&args[i]);
    }
    return args;
}


void run_threads(pthread_t* threads, uint32_t number_of_threads) {
    for(uint32_t i = 0; i < number_of_threads; i++) {
        pthread_join(threads[i], NULL);
    }
}


void* run_heat_calculations(void* void_args) {
    thread_args* args = (thread_args*)void_args;
    bool local_grid1_older = true;
    uint16_t local_cycle_count = 0;

    for(uint16_t tick = 0; tick < args->time_ticks; tick++) {
        for(uint16_t x = args->start_x; x < args->end_x; x++) {
            for(uint16_t y = 0; y < GRID_Y_SIZE; y++) {
                if(local_grid1_older) heat_grid2[x][y] = calc_heat_value(heat_grid1, (int)x, (int)y);
                else heat_grid1[x][y] = calc_heat_value(heat_grid2, (int)x, (int)y);
            }
        }

        // Thread 0 prints the required info at 200 cycle intervals
        // Do this first so that cycle 200 is complete but we haven't flipped
        // to the new array yet - actually showing the data from cycle 200
        if(args->thread_id == 0) {
            if( (local_cycle_count % 200) == 0 ) {
                if(local_grid1_older) print_interval_step(heat_grid2);
                else print_interval_step(heat_grid1);
            }
        }

        // Flip to the other array to use as the reference data
        local_cycle_count++;
        if(local_grid1_older) local_grid1_older = false;
        else local_grid1_older = true;

        // This is the only place a thread can write to this value
        pthread_mutex_lock(&mutexs[args->thread_id]);
            thread_cycles[args->thread_id] = local_cycle_count;
            pthread_cond_broadcast(&condition_vars[args->thread_id]);
        pthread_mutex_unlock(&mutexs[args->thread_id]);


        if(args->thread_id != 0) {
            pthread_mutex_lock(&mutexs[args->thread_id - 1]);
                while( local_cycle_count != thread_cycles[args->thread_id - 1] )
                    pthread_cond_wait(&condition_vars[args->thread_id - 1], &mutexs[args->thread_id - 1]);
            pthread_mutex_unlock(&mutexs[args->thread_id - 1]);
        }

        pthread_mutex_lock(&mutexs[args->thread_id]);
            pthread_cond_broadcast(&condition_vars[args->thread_id]);
        pthread_mutex_unlock(&mutexs[args->thread_id]);

        if(args->thread_id != (number_of_threads - 1)) {
            pthread_mutex_lock(&mutexs[args->thread_id + 1]);
                while( local_cycle_count != thread_cycles[args->thread_id + 1] )
                    pthread_cond_wait(&condition_vars[args->thread_id + 1], &mutexs[args->thread_id + 1]);
            pthread_mutex_unlock(&mutexs[args->thread_id + 1]);

        }

        pthread_mutex_lock(&mutexs[args->thread_id]);
            pthread_cond_broadcast(&condition_vars[args->thread_id]);
        pthread_mutex_unlock(&mutexs[args->thread_id]);
    }

    pthread_exit(0);
}


double get_heat_value(double grid[GRID_X_SIZE][GRID_Y_SIZE], int x, int y) {
    const double out_of_bounds_heat_value = 0;
    double heat_value;
    bool x_out_of_bounds = (x < 0) || (x >= GRID_X_SIZE);
    bool y_out_of_bounds = (y < 0) || (y >= GRID_Y_SIZE);

    if( x_out_of_bounds || y_out_of_bounds) heat_value = out_of_bounds_heat_value;
    else heat_value = grid[x][y];

    return heat_value;
}


double calc_heat_value(double past_grid[GRID_X_SIZE][GRID_Y_SIZE],
                       int x, int y) {
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

void cleanup() {
    int retval;
    for(uint32_t i = 0; i < number_of_threads; i++) {
        retval = pthread_mutex_destroy(&mutexs[i]);
        if(retval) {
            printf("ERROR: Could not destroy mutex number %u\n", i);
            exit(-2);
        }
        retval = pthread_cond_destroy(&condition_vars[i]);
        if(retval) {
            printf("ERROR: Could not destroy cond var number %u\n", i);
            exit(-2);
        }
    }

    free(threads);
    free(mutexs);
    free(condition_vars);
    free(thread_cycles);
    free(args);
}

