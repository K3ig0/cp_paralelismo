
/*Ejecucion:

mpicc mpi-pi.c -o mpi-pi
mpirun -np 5 ./mpi-pi

*/

#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int i, n;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x;

    int numprocs, my_id;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    printf("my_id %d numprocs %d\n", my_id, numprocs);

    if (numprocs > 1) {

        if (my_id == 0) {
            printf("Enter the number of intervals: (0 quits) \n");
            scanf("%d",&n);
        }
        while (1)
        {
            if (n == 0) break;
            

            else {



                //codigo



            }
        }
    }

    // Solo un proceso
    else {
        while (1) {
            printf("No paralelism. Only 1 proccess running..");
            printf("Enter the number of intervals: (0 quits) \n");
            scanf("%d",&n);
        
            if (n == 0) break;
      
            h   = 1.0 / (double) n;
            sum = 0.0;
            for (i = 1; i <= n; i++) {
                x = h * ((double)i - 0.5);
                sum += 4.0 / (1.0 + x*x);
            }
            pi = h * sum;

            printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize();

}
