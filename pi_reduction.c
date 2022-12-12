/*
	CS-770 Big Data Analytics
	HW3	OpenMP reduction   
	Niko Grisel Todorov
*/
#include <omp.h>
#include <stdio.h>
static int num_cores = 32;
static long num_steps = 1000000000;
double step;

void main()
{
	printf("requested cores = %d\n", num_cores);
	printf("requested steps = %d\n", num_steps);

	int i;
	double x, sum, pi, extime, start;
	step = 1.0 / (double)num_steps;

	for (i = 1; i <= num_cores; i++) {
		sum = 0.0;
		omp_set_num_threads(i);
		start = omp_get_wtime();

#pragma omp parallel for reduction(+:sum)
		for (i = 1; i < num_steps; i++) {
			x = (i + 0.5) * step; 
			sum += 4.0 / (1.0 + x * x);
		}
		// back to serial
		pi = sum * step;
		extime = omp_get_wtime() - start;
		printf("%7.5fs:%2d pi = %30.28f\n", extime, i, pi);
	}
}
