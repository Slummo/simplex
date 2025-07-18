#include "problem.h"
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_errno.h>

#include <sys/time.h>
#include <sys/resource.h>

void print_performance_report(struct timeval* restrict t_start, struct timeval* restrict t_end,
                              struct rusage* restrict usage) {
    puts("\n\nPerformance report:");
    if (getrusage(RUSAGE_SELF, usage) != 0) {
        perror("getrusage");
    } else {
        printf("  -Peak used RAM: %ld KB\n", usage->ru_maxrss);
    }

    if (gettimeofday(t_end, NULL) != 0) {
        perror("gettimeofday");
    } else {
        double elapsed = (t_end->tv_sec - t_start->tv_sec) + (t_end->tv_usec - t_start->tv_usec) * 1e-6;
        printf("  -Elapsed time: %lf\n", elapsed);
    }
}

int main(int argc, char** args) {
    struct timeval t_start, t_end;
    struct rusage usage;

    if (gettimeofday(&t_start, NULL) != 0) {
        perror("gettimeofday");
        return EXIT_FAILURE;
    }

    if (argc > 2) {
        fprintf(stderr, "Wrong number of arguments! Usage: 'zmax <filename>' or just 'zmax' to read from stdin\n");
        return EXIT_FAILURE;
    }

    gsl_set_error_handler_off();

    FILE* stream = argc == 2 ? fopen(args[1], "r") : stdin;
    if (!stream) {
        perror("Failed to determine stream");
        return EXIT_FAILURE;
    }

    problem_t problem = {0};
    if (!problem_from_stream(&problem, stream)) {
        fprintf(stderr, "Failed to create problem\n");
        return EXIT_FAILURE;
    }

    problem_print(&problem, "Problem");

    solution_t solution = {0};
    if (!problem_solve(&problem, &solution)) {
        fprintf(stderr, "Failed to solve problem\n");
        problem_free(&problem);
        solution_free(&solution);
        return EXIT_FAILURE;
    }

    solution_print(&solution, "Solution");

    problem_free(&problem);
    solution_free(&solution);
    print_performance_report(&t_start, &t_end, &usage);

    return EXIT_SUCCESS;
}
