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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>


static uint32_t num_points_per_thread;
static uint32_t global_circle_hit_count = 0;
static pthread_mutex_t circle_count_lock;


uint32_t parse_cmdline(int, char**, int, int);
void get_system_time(struct timespec*);
unsigned long long int calc_runtime(struct timespec, struct timespec);
void print_results(double, unsigned long long int);


pthread_t* allocate_threads(uint32_t);
void create_threads(pthread_t*, uint32_t, uint32_t);
void run_threads(pthread_t*, uint32_t);
void* process_points(void*);
double calculate_pi(uint32_t, uint32_t);


//*****************************************************************************
// Functions
//*****************************************************************************
int main(int argc, char* argv[]) {
    int args_expected = 3;
    double pi;
    struct timespec start_time, end_time;

    get_system_time(&start_time);
    uint32_t number_of_points = parse_cmdline(argc, &argv[0], 1, args_expected);
    uint32_t number_of_threads = parse_cmdline(argc, &argv[0], 2, args_expected);

    int retval = pthread_mutex_init(&circle_count_lock, NULL);
    if(retval) {
        printf("ERROR: Failed to init circle_count_lock mutex.\n");
        exit(-1);
    }
    pthread_t* threads = allocate_threads(number_of_threads);
    create_threads(threads, number_of_points, number_of_threads);
    run_threads(threads, number_of_threads);

    get_system_time(&end_time);
    unsigned long long int runtime = calc_runtime(start_time, end_time);
    pi = calculate_pi(global_circle_hit_count, number_of_points);
    print_results(pi, runtime);

    free(threads);
    pthread_mutex_destroy(&circle_count_lock);

    return 0;
}


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


pthread_t* allocate_threads(uint32_t number_of_threads) {
    pthread_t* threads = malloc( sizeof(pthread_t) * number_of_threads );
    if(threads == NULL) {
        printf("ERROR: Threads malloc returned NULL. Attempted to malloc %u threads\n", number_of_threads);
        exit(-1);
    }
    return threads;
}


void create_threads(pthread_t* threads, uint32_t num_points, uint32_t number_of_threads) {
    num_points_per_thread = num_points / number_of_threads;
    time_t current_time = time(NULL);
    time_t seed;

    for(uint32_t i = 0; i < number_of_threads; i++) {
        // Add in some variance so that each thread has a different seed. Just calling
        // time makes the seed the same for all threads
        seed = current_time + (time_t)i;
        pthread_create(&threads[i], NULL, &process_points, (void*)seed);
    }
}


void run_threads(pthread_t* threads, uint32_t number_of_threads) {
    for(uint32_t i = 0; i < number_of_threads; i++) {
        pthread_join(threads[i], NULL);
    }
}


// Generates a threading safe random number between 0 and 1
static double gen_random_num(unsigned int* seed) {
    double random_num = (double)( (((double)rand_r(seed)) + 1.0) / (RAND_MAX + 2.0) );
    return random_num;
}


static double square_cord(double cord) {
    double squared_cord = cord * cord;
    return squared_cord;
}


// Function that is called in many different threads count the number of points in the circle
void* process_points(void* void_args) {
    double random_num_x, random_num_y;
    double x_cord, y_cord;
    unsigned int seed = (unsigned int)void_args;
    unsigned int seed1 = seed;
    unsigned int seed2 = seed + 1;
    uint32_t in_circle_count = 0;

    for(uint32_t i = 0; i < num_points_per_thread; i++) {
        random_num_x = gen_random_num(&seed1);
        random_num_y = gen_random_num(&seed2);
        x_cord = (2.0 * random_num_x) - 1.0;
        y_cord = (2.0 * random_num_y) - 1.0;

        // This somehow figures out if the cords are in the circle
        if( (square_cord(x_cord) + square_cord(y_cord)) <= 1.0 ) {
            in_circle_count++;
        }
    }
    pthread_mutex_lock(&circle_count_lock);
        global_circle_hit_count += in_circle_count;
    pthread_mutex_unlock(&circle_count_lock);

    pthread_exit(0);
}


double calculate_pi(uint32_t total_hits, uint32_t num_points) {
    double pi; 
    pi = (4.0 * (double)total_hits) / (double)num_points;
    return pi;
}

