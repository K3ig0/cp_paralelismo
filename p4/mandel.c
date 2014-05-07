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
    int i, j, k, p;
    Compl   z, c;
    float   lengthsq, temp;
    int res[X_RESN][Y_RESN]; 

    // Spliting the job
    int numprocs, my_id, n;
    unsigned long bal = 0;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    n = X_RESN / (numprocs-1);

    if (my_id == 0) {
    k=0;
        while (k < n) {
            for (p=0; p < numprocs-1; p++){
                TODO: Enviar su parte
                k++;
            }
        }
    }
    else{
        while (0){
            if (noReciboNada) MPI_Finalize();
TODO: Asignar i para que ejecute la parte que le manden
          for (j=0; i< Y_RESN; j++) {
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
TODO: Enviar los datos calculados al master of the puppets
        }
    }
    //fprintf(stderr, "The proccess %d has done %lu float point operations \n", my_id, bal);

    /*if (my_id == 0)
      MPI_Gather(MPI_IN_PLACE, area*Y_RESN, MPI_INT, 
      res, area*Y_RESN, MPI_INT, 0, MPI_COMM_WORLD);   
      else 
      MPI_Gather(res[my_id*area], area*Y_RESN, MPI_INT, 
      res, area*Y_RESN, MPI_INT, 0, MPI_COMM_WORLD);
      */
    if( DEBUG && (my_id == 0)) {
        for(i=0;i<X_RESN;i++) {
            for(j=0;j<Y_RESN-1;j++) {
                printf("%d\t", res[i][j]);
            }
            printf("%d\n", res[i][Y_RESN-1]);}
    }
    MPI_Finalize();
}
