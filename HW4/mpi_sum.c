#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>

#define MAX 1000000000     /* number to compute sum */
#define MASTER 0           /* task ID of master task */

int main (int argc, char *argv[])
{
struct timespec StartTime, EndTime;
int ret;
double	localsum, sum; /* local sum and global sum */
long long localsum_int; /* temporary integer sum used in this task */

int	taskid,	        /* task ID */
	numtasks,       /* number of tasks */
	rc;             /* return code */
long    start, end, localcount, i;

  MPI_Status status;

  /* Obtain number of tasks and task ID */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  printf ("MPI task %d has started...\n", taskid);
  
  if (taskid == MASTER) {
    ret = clock_gettime(CLOCK_REALTIME, &StartTime);
    assert(ret == 0);
  }

  localcount = MAX/numtasks;
  start = localcount*taskid + 1;
  end = start + localcount - 1;
  if (taskid == (numtasks-1))  /* last task */
    end = MAX;

  localsum_int = 0;

  for (i = start; i <= end; i++) {
     localsum_int += i;
  }

  localsum = (double)localsum_int;

   /* Use MPI_Reduce to sum values of localsum across all tasks 
    * Master will store the accumulated value in sum 
    * - localsum is the send buffer
    * - sum is the receive buffer (used by the receiving task only)
    * - the size of the message is sizeof(double)
    * - MASTER is the task that will receive the result of the reduction
    *   operation
    * - MPI_SUM is a pre-defined reduction function (double-precision
    *   floating-point vector addition).  Must be declared extern.
    * - MPI_COMM_WORLD is the group of tasks that will participate.
    */

   rc = MPI_Reduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM,
                   MASTER, MPI_COMM_WORLD);
   if (rc != MPI_SUCCESS)
      printf("%d: failure on mpi_reduce\n", taskid);

   /* Master prints sum */
   if (taskid == MASTER) { 
      ret = clock_gettime(CLOCK_REALTIME, &EndTime);
      assert(ret == 0);
      printf("   sum = %40.1f\n", sum);
      unsigned long long int runtime = 1000000000 * (EndTime.tv_sec - StartTime.tv_sec) + EndTime.tv_nsec - StartTime.tv_nsec; 
      printf("\nTime = %lld nanoseconds\t(%d.%09lld sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
   }     

  MPI_Finalize();
  return 0;
}
