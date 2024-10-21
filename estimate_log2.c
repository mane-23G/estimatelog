/*******************************************************************************
  Title          : estimate_log.c
  Author         : Mane Galstyan
  Created on     : March 6, 2024
  Description    : Interactively estimate log
  Purpose        : To use MPI broadcast and MPI reduce
  Usage          : mpirun --use-hwthread-cpus estimate_log <log_number> <num_intervals>
  Build with     : mpicc -Wall -o estimat_log estimate_log.c -lm
  Modifications  :
 
*******************************************************************************/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <mpi.h>

#define ROOT 0
#define USG_MSG "usage: %s <logarithm number> <number of intervals>"

double approximate_log(double log_number, int num_segments, int id, int p);
/**  approximate_log()
 *  This returns an approximation of ln(
 *  To compute pi, we can approximate arctan(1.0). We use the fact that
 *      d/dx(ln(x) = 1/(x)
 *  and compute the area under the curve of 1/(x) in the interval
 *  from 0 to log_number.
 *  We use the rectangle rule to approximate the area under the curve.
 *  We divide [0,log_number] by num_segments, and compute the
 *  value at the midoint of each segment with dx * (i-0.5) + 1 By summing these
 *  values and multiplying by dx we have the approximation.
*/

int main( int argc, char *argv[] )
{
    int id;             /* rank of executing process   */
    int p;              /* number of processes         */
    double log_estimate; /* estimated value of pi       */
    double local_log;    /* each process's contribution */
    double elapsed_time; /* used to measure the time */
    double log_number;
    int num_intervals;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    
    /* processor 0 is the only one to read and wrtie*/
    
        /*check to see if there is the correct number of argument*/
        if(argc != 3) {
            if(id == ROOT)
                printf("usage: %s <logarithm number> <number of intervals>\n",argv[0]);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    
        /*parse through command line first argument*/
        errno = 0;
        log_number = strtod(argv[1], NULL);
        if(errno == ERANGE || errno == EINVAL || log_number < 1) {
            if(id == ROOT)
                printf("Found invalid argument %s\n", argv[1]);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        /*parse through command line second argument*/
        errno = 0;
        num_intervals = strtol(argv[2], NULL, 10);
        if(errno == ERANGE || num_intervals < 1) {
            if(id == ROOT)
                printf("Found invalid argument %s\n", argv[2]);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

    /*start timing the program*/
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = MPI_Wtime();
    
    /*no need to broadcast since all the processors process the command line*/
//    MPI_Bcast(&num_intervals, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
//    MPI_Bcast(&log_number, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    /* each process calcultes the approximation for their intervals */
    local_log = approximate_log(log_number, num_intervals, id, p);

        /* MPI_Reduce collects the local estimate from each process into
           a global value, log_estimate. The ROOT is the process that receives
           all values in the reduction. The reduce operator is MPI_SUM. */
    MPI_Reduce(&local_log, &log_estimate, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        
    /*end timing the program */
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    /* ROOT does the printing. The error is calculated by comparing the actual log(log_number) */
    if(ROOT == id) {
        printf("%.16g   %.16f   %.16f   %.6f seconds\n", log_number, log_estimate, fabs(log_estimate - log(log_number)),elapsed_time);
        fflush(stdout);
    }

    /*call MPI_Finalize() to close MPI */
    MPI_Finalize();
    return 0;
}

double approximate_log (double log_number, int num_segments, int id, int p)
{
    double dx, sum, x;
    int i;
    
    /* Set dx to the width of each segments */
    dx = (log_number - 1.0) / num_segments;
    
    /* Initialize sum */
    sum = 0.0;
    
    /* Each process will compute its share of the segments. If the
     segments are numbered 1, 2, 3, ...,n, from left to right, then
     process id computes segment k if id = (k-1) % p, or equivalently
     it computes segments id+1, id+p+1, id+2p+1, ... up to id+mp+1,
     where m is the largest number such that id+mp+1 <= num_segments. */
    for (i = id + 1; i <= num_segments; i += p) {
        
        x = dx * ((double)i - 0.5); /* x is midpoint of segment i */
        sum += 1.0 / (x + 1.0);   /* add new area to sum */
    }
    /* we multiply sum by dx because we are computing an integral and dx is the differential */
    return dx * sum;
}
