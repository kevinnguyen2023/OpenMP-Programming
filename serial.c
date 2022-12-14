/* File:      serial.c
 * Purpose:  Estimates the value of the natural logarithm of 2: ln(2)
 *
 *              ln(2) = 1 -1/2 + 1/3 - 1/4 ....
 *              nth term is (-1)^(n+1) 1/n
 *
 * Run:      serial  <n>
 *           n is the number of terms of the series to use
 *
 * Input:    none
 * Output:   The estimate of ln(2) and the value of ln(2) computed by the
 *           in the math library
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> 

void Usage(char* prog_name);

int main(int argc, char* argv[]) {
   long n, i;
   int sign = 1.0;
   double sum = 0.0;
   
   double start, end; 

   if (argc != 2) Usage(argv[0]);
   n = strtoll(argv[1], NULL, 10);

   /* START OF CODE TO BE PARALLELIZED */
   
   start = omp_get_wtime();
   
   //double x = 1.0;
   // num_threads = number of threads with specific number 
   #   pragma omp parallel num_threads(2) firstprivate(sign) reduction(+:sum) 
   {

#  pragma omp for
   for (i = 1; i < n; i++) {
      sum += (double)sign*1/i;
      if (i % 2 == 0)
         sign = 1.0; 
      else 
        sign = -1.0; 
       //x += sign/(2*i + 1); 
      //sign = -sign 
   }   
   
}  /* END OF CODE TO BE PARALLELIZED */
   end = omp_get_wtime(); 
   printf("Elapsed time %f \n", end - start);

   printf("   Program estimate of ln(2) = %.14f\n", sum);
   printf("     From Math Library ln(2) = %.14f\n", log(2));
   return 0;
}  /* main */

/*------------------------------------------------------------------
 * Function:  Usage
 * Purpose:   Print a message explaining how to run the program
 * In arg:    prog_name
 */
void Usage(char* prog_name) {
   fprintf(stderr, "usage: %s  <n>\n", prog_name);  /* Change */
   fprintf(stderr, "   n is the number of terms and should be >= 1\n");
   exit(0);
}  /* Usage */
