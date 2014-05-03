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



typedef struct complextype
{
    float real, imag;
} Compl;


int main (int argc, char *argv[] )
{

    /* Mandelbrot variables */
    int i, j, k;
    Compl   z, c;
    float   lengthsq, temp;
    int res[X_RESN][Y_RESN]; 

    // Spliting the job
    int numprocs, my_id, area, restantes, h;
    unsigned long bal = 0;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    int recvcounts[numprocs];
    int displs[numprocs];

    area = (int) (X_RESN / numprocs);
    restantes = X_RESN % numprocs;

    // Generating recvcounts and displs arrays
    // All process do it because we use it in the for loop
    recvcounts[0] = (area + restantes)*Y_RESN;
    displs[0] = 0;
    for (h=1; h<numprocs; h++){
        recvcounts[h] = area * Y_RESN;
        displs[h] = displs[h-1] + recvcounts[h-1];
    }
    
    for (i=(displs[my_id]/Y_RESN); i < ((displs[my_id]/Y_RESN)+(recvcounts[my_id]/Y_RESN)); i++)
        for(j=0; j < Y_RESN; j++) {
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
    fprintf(stderr, "The proccess %d has done %lu float point operations \n", my_id, bal);

    if (my_id == 0)
        MPI_Gatherv(MPI_IN_PLACE, (area+restantes)*Y_RESN, MPI_INT, 
                res, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);   
    else 
        MPI_Gatherv(res[(my_id*area)+restantes], area*Y_RESN, MPI_INT, 
                res, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);   	
    if( DEBUG && (my_id == 0)) {
        for(i=0;i<X_RESN;i++) {
            for(j=0;j<Y_RESN-1;j++) {
                printf("%d\t", res[i][j]);
            }
            printf("%d\n", res[i][Y_RESN-1]);}
    }
    MPI_Finalize();
}
