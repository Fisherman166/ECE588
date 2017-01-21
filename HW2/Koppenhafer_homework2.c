#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>

typedef struct {
    uint32_t num_points;
    time_t seed;
} args_for_threads;


void print_error(char*);
uint32_t parse_cmdline(int, char**, int, int);
void get_system_time(struct timespec*);
unsigned long long int calc_runtime(struct timespec, struct timespec);
void print_results(double, unsigned long long int);

pthread_t* allocate_threads(uint32_t);
args_for_threads* create_threads(pthread_t*, uint32_t, uint32_t);
uint32_t** run_threads(pthread_t*, uint32_t);
void* process_points(void*);
double calculate_pi(uint32_t**, uint32_t, uint32_t);

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

    pthread_t* threads = allocate_threads(number_of_threads);
    args_for_threads* args = create_threads(threads, number_of_points, number_of_threads);
    uint32_t** results = run_threads(threads, number_of_threads);

    pi = calculate_pi(results, number_of_points, number_of_threads);
    get_system_time(&end_time);
    unsigned long long int runtime = calc_runtime(start_time, end_time);
    print_results(pi, runtime);

    free(threads);
    free(args);
    for(uint32_t i = 0; i < number_of_threads; i++) {
        free(results[i]);
    }
    free(results);

    return 0;
}


void print_error(char* string) {
    printf("ERROR: %s\n", string);
}


uint32_t parse_cmdline(int argc, char** argv, int arg_position, int args_expected) {
    char error_str[100];
    uint32_t number_on_cmdline;
    if(argc < args_expected) {
        sprintf(&error_str[0], "Argc is less than %d for number of threads. Exiting.", args_expected);
        print_error(&error_str[0]);
        exit(-1);
    }
    if(arg_position >= argc) {
        sprintf(&error_str[0], "Asked for argument number %d, which is equal to or larger than argc: %d. Exiting.", args_expected, argc);
        print_error(&error_str[0]);
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


args_for_threads* create_threads(pthread_t* threads, uint32_t num_points, uint32_t number_of_threads) {
    uint32_t points_per_thread = num_points / number_of_threads;
    time_t current_time = time(NULL);
    args_for_threads* args = malloc( sizeof(args_for_threads) * number_of_threads );
    if(args == NULL) {
        printf("ERROR: args malloc returned NULL. Attempted to malloc %u args\n", number_of_threads);
        exit(-1);
    }

    for(uint32_t i = 0; i < number_of_threads; i++) {
        args[i].num_points = points_per_thread;
        // Add in some variance so that each thread has a different seed. Just calling
        // time makes the seed the same for all threads
        args[i].seed = current_time + (time_t)i;
        pthread_create(&(threads[i]), NULL, &process_points, (void*)&(args[i]));
    }
    return args;
}


uint32_t** run_threads(pthread_t* threads, uint32_t number_of_threads) {
    uint32_t** results = malloc( sizeof(uint32_t*) * number_of_threads );
    if(results == NULL) {
        printf("ERROR: results malloc returned NULL. Attempted to malloc %u args\n", number_of_threads);
        exit(-1);
    }

    for(uint32_t i = 0; i < number_of_threads; i++) {
        pthread_join(threads[i], (void**)&results[i]);
    }
    
    return results;
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


// Function that is called in many different threads to process the results
void* process_points(void* void_args) {
    double random_num_x, random_num_y;
    double x_cord, y_cord;
    args_for_threads* args = (args_for_threads*)void_args;
    unsigned int seed1 = (unsigned int)(args->seed);
    unsigned int seed2 = (unsigned int)((args->seed) + 1);
    uint32_t* in_circle_count = malloc( sizeof(uint32_t) );
    if(in_circle_count == NULL) {
        printf("ERROR: in_circle_count malloc returned NULL.\n");
        exit(-1);
    }

    for(uint32_t i = 0; i < args->num_points; i++) {
        random_num_x = gen_random_num(&seed1);
        random_num_y = gen_random_num(&seed2);
        x_cord = (2.0 * random_num_x) - 1.0;
        y_cord = (2.0 * random_num_y) - 1.0;

        // This somehow figures out if the cords are in the circle
        if( (square_cord(x_cord) + square_cord(y_cord)) <= 1.0 ) {
            *in_circle_count += 1;
        }
    }
    return in_circle_count;
}


double calculate_pi(uint32_t** circle_counts, uint32_t num_points, uint32_t num_of_threads) {
    double pi; 
    uint32_t total_score = 0;

    for(uint32_t i = 0; i < num_of_threads; i++) {
        total_score += *circle_counts[i];
    }
    pi = (4.0 * (double)(total_score)) / (double)(num_points);

    return pi;
}

