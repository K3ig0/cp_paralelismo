/*The Mandelbrot set is a fractal that is defined as the set of points c
  in the complex plane for which the sequence z_{n+1} = z_n^2 + c
  with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

/*Ejecucion:

  mpicc mandel.c -o mandel
  mpirun -np 4 ./mandel > mandel.txt
  python2 view.py mandel.txt

*/

#include <stdio.h>
#include <mpi.h>

#define DEBUG 1

#define         X_RESN  1024     /* x resolution */
#define         Y_RESN  1024       /* y resolution */
#define         X_MIN   -2.0	/*Boundaries of the mandelbrot set*/
#define         X_MAX    2.0
#define         Y_MIN   -2.0
#define         Y_MAX    2.0
#define		maxIterations	1000 /*Cuanto mayor, más detalle en la imagen y mayor coste computacional*/
#define WORK_TAG 1
#define DIE_TAG 2

typedef struct complextype
{
    float real, imag;
} Compl;


int main (int argc, char *argv[] )
{

    /* Mandelbrot variables */
    int i, j, k, p, t, s, tag;
    Compl   z, c;
    float   lengthsq, temp;
    int res[X_RESN][Y_RESN]; 
    int temp_row[Y_RESN];

    // Spliting the job
    int numprocs, my_id, n;
    unsigned long bal = 0;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    n = X_RESN;
    if (my_id == 0) {
        k=0;
        t=0;
        for (p=1; p < numprocs; p++){
            if (t>n-1) break;
            MPI_Send(&t, 1, MPI_INT, p, WORK_TAG, MPI_COMM_WORLD);
            t++;
        }
        while (k<n) {
            MPI_Recv(&temp_row, Y_RESN, MPI_INT, MPI_ANY_SOURCE, 
                    MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            for (s=0; s<Y_RESN; s++) res[status.MPI_TAG][s] = temp_row[s];
            k++;
            if (t >= n) tag = DIE_TAG;
            else tag = WORK_TAG;
            MPI_Send(&t, 1, MPI_INT, status.MPI_SOURCE, tag, MPI_COMM_WORLD);
            t++;
        }
    }
    else{
        while (1){
            MPI_Recv(&i, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == DIE_TAG) break;
            for (j=0; j< Y_RESN; j++) {
                z.real = z.imag = 0.0;
                c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
                c.imag = Y_MAX - i * (Y_MAX - Y_MIN)/Y_RESN;
                k = 0;

                do  {    /* iterate for pixel color */

                    temp = z.real*z.real - z.imag*z.imag + c.real;
                    z.imag = 2.0*z.real*z.imag + c.imag;
                    z.real = temp;
                    lengthsq = z.real*z.real+z.imag*z.imag;
                    k++;

                } while (lengthsq < 4.0 && k < maxIterations);
                bal = bal + (k*10); //10 operations

                if (k >= maxIterations) res[i][j] = 0;
                else res[i][j] = k;
            }
            MPI_Send(&res[i], Y_RESN, MPI_INT, 0, i, MPI_COMM_WORLD);
        }
    }
    //fprintf(stderr, "The proccess %d has done %lu float point operations \n", my_id, bal);
    if( DEBUG && (my_id == 0)) {
        for(i=0;i<X_RESN;i++) {
            for(j=0;j<Y_RESN-1;j++) {
                printf("%d\t", res[i][j]);
            }
            printf("%d\n", res[i][Y_RESN-1]);}
    }

    MPI_Finalize();
}
