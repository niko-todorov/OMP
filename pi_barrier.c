/*
	CS-770 Big Data Analytics
	HW3	OpenMP barrier   
	Niko Grisel Todorov
*/
#include <omp.h>
#include <stdio.h>
static int num_cores = 4;
static long num_steps = 1000000000;
//double step;

void main()
{
	printf("requested cores = %d\n", num_cores);
	printf("requested steps = %d\n\n", num_steps);

	double pi = 0.0;
	omp_set_num_threads(num_cores);
	double extime, start = omp_get_wtime();
#pragma omp parallel
	{
		int i, tid, num; 
		double x, sum, step;
		step = 1.0 / (double)num_steps;
		tid = omp_get_thread_num();
		num = omp_get_num_threads();
		for (i = tid, sum = 0.0; i < num_steps; i = i + num) {
			x = (i + 0.5) * step; 
			sum += 4.0 / (1.0 + x * x);
		}
		// printf("thread #%2d sum = %30.28f\n", tid, sum*step); // bad
#pragma omp barrier
		printf("thread #%2d sum = %30.28f\n", tid, sum * step); // good
		pi += sum * step;
		// printf("pi = %30.28f from thread #%d\n", pi, tid);
	}
	extime = omp_get_wtime() - start;
	printf("\n%7.5fs:%2d pi = %30.28f\n", extime, num_cores, pi);
}
