/*
	CS-770 Big Data Analytics
	HW3	OpenMP reduction   
	Niko Grisel Todorov
*/
#include <omp.h>
#include <stdio.h>
static int _Num_threads = 32;
static long _Num_steps = 1000000000;
double step;

void main()
{
	int _Num_CPUs = omp_get_num_procs();
	printf("available cores = %d\n", _Num_CPUs);
	printf("set use threads = %d\n", _Num_threads);
	printf("requested steps = %d\n", _Num_steps);

	int i;
	double x, sum, pi, extime, start;
	step = 1.0 / (double)_Num_steps;

	for (i = 1; i <= _Num_threads; i++) {
		sum = 0.0;
		omp_set_num_threads(i);
		start = omp_get_wtime();

#pragma omp parallel for reduction(+:sum)
		for (i = 1; i < _Num_steps; i++) {
			x = (i + 0.5) * step; 
			sum += 4.0 / (1.0 + x * x);
		}
		// back to serial
		pi = sum * step;
		extime = omp_get_wtime() - start;
		printf("%7.5fs:%2d pi = %30.28f\n", extime, i, pi);
	}
}
