
/*Ejecucion:

mpicc mpi-pi.c -o mpi-pi
mpirun -np 5 ./mpi-pi

*/

#include <stdio.h>
#include <math.h>
#include <mpi.h>


int MPI_FattreeColectiva(void * buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm){
    int my_id, numprocs, i, ret;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if (my_id == 0) {
        for (i=1; i < numprocs; i++){
            MPI_Send(buffer, count, datatype, i, root, comm);
        }
    }
    else {
        ret = MPI_Recv(buffer, count, datatype, root, 0, comm, &status);
    }
    return ret;
}

int generate_steps(int numprocs){
    return (int)ceil(log10(numprocs) / log10(2));
}

int MPI_BinomialColectiva(void * buffer, int count, MPI_Datatype datatype,
        int root, MPI_Comm comm){
    int my_id, numprocs, i, ret, steps;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    steps = generate_steps(numprocs);
    for (i=1; i <= steps; i++){
        if ((my_id < (pow(2, (i-1))))) {
            if (numprocs<=(my_id+(pow(2, (i-1))))) break;
            MPI_Send(buffer, count, datatype, (my_id+(pow(2, (i-1)))), 0, comm);
         }
        else if (my_id < (pow(2, i))){
            ret = MPI_Recv(buffer, count, datatype, (my_id-(pow(2, (i-1)))), 0, comm, &status);
         }
    }
    return ret; 
}

int main(int argc, char *argv[])
{
    int i, n;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x, sum1;

    int numprocs, my_id;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    printf("my_id %d numprocs %d\n", my_id, numprocs);

    if (numprocs > 1) { //Empieza código paralelizado
        while (1) {
            if (my_id == 0) {
                printf("Enter the number of intervals: (0 quits) \n");
                scanf("%d",&n);
            }
            MPI_BinomialColectiva(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); 
            if (n == 0) break;
            /*Start alg*/ 
            h   = 1.0 / (double) n;
            sum = 0.0;
            for (i = my_id+1; i <= n; i+=numprocs) {
                x = h * ((double)i - 0.5);
                sum += 4.0 / (1.0 + x*x);
            }
            /*Finish alg*/
            MPI_Reduce(&sum, &sum1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (my_id == 0) {
                pi = h * sum1;
                printf("(%d) pi is approximately %.16f, Error is %.16f\n", my_id,
                        pi, fabs(pi - PI25DT));
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
