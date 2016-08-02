#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <complex.h>
#include "Maxfiles.h"
#include "MaxSLiCInterface.h"


const int N = 64; //Poisson_N;

// Helper macro
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void four1(double data[], unsigned long nn, int isign);

void four1Call(double complex data[], unsigned long nn, int isign){
	four1((double*)(data)-1,nn,isign);
}

// Get time in seconds
static
double timesec()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return ((double) (t.tv_sec * 1e6L + t.tv_usec)) / 1.0e6;
}

// Allocate the specified number of bytes in memory and initialise to 0
static
void *allocate_mem(size_t bytes)
{
	void *ptr = calloc(bytes, 1);
	if (ptr == NULL) {
		fprintf(stderr, "Error: Unable to allocate %lu bytes\n", bytes);
		abort();
	}
	return ptr;
}

// Check if a and b are equal according to the specified relative tolerance.
static
bool crelequal(float complex a, float complex b, float relative_tolerance)
{
    float diff = cabsf(a - b);
    float largest = MAX(cabsf(a), cabsf(b));

    if (diff <= largest * relative_tolerance)
        return true;
    return false;
}

// CPU FFT.
// Perform an FFT on the input data and return the output in a newly allocated buffer.
// The buffer will need to be deallocated using free().
static
float complex *fft_cpu(float complex *input)
{
	double complex f[N];

	// FFT dim1 of rho
	for(int i=0; i<N;i++){
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++)
				f[k] = input[(i*N+j)*N+k];
			four1Call(f,N,1);
			for (int k = 0; k < N; k++)
				input[(i*N+j)*N+k] = f[k];

		}
	}

	// FFT dim2 of rho
	for (int k = 0; k < N; k++) {
		for(int i=0;i<N;i++){
			for (int j = 0; j < N; j++)
				f[j] = input[(i*N+j)*N+k];
			four1Call(f,N,1);
			for (int j = 0; j < N; j++)
				input[(i*N+j)*N+k] = f[j];
		}
	}

	// FFT dim3 of rho
	for (int j = 0; j < N; j++) {
		for(int k=0;k<N;k++){
			for (int i = 0; i < N; i++)
				f[i] = input[(i*N+j)*N+k];
			four1Call(f,N,1);
			for (int i = 0; i < N; i++)
				input[(i*N+j)*N+k] = f[i];
		}
	}

	// solve equation in Fourier space
	float h=1/(float)N;
	double pi = 4 * atan(1.0);
	float complex W = cos(2.0 * pi / N) +I*sin(2.0 * pi / N);
	float complex Wl=1.0, Wm = 1.0, Wn = 1.0;
	for(int l=0; l<N; l++){
		for (int m = 0; m < N; m++) {
			for (int n = 0; n < N; n++) {
				float complex denom = Wm + 1.0 / Wm + Wn + 1.0 / Wn + Wl + 1.0 / Wl - 6.0;
				if (denom != 0.0)
					input[(l*N+m)*N+n] *= h * h / denom;
				Wn *= W;
			}
			Wm *= W;
		}
		Wl *= W;
	}

	// inverse FFT dim1 of rho
	for(int i=0; i<N; i++){
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++)
				f[k] = input[(i*N+j)*N+k];
			four1Call(f,N,-1);
			for (int k = 0; k < N; k++)
				input[(i*N+j)*N+k] = f[k];
		}
	}

	// inverse FFT dim2 of rho
	for (int k = 0; k < N; k++) {
		for(int i=0;i<N;i++){
			for (int j = 0; j < N; j++)
				f[j] = input[(i*N+j)*N+k];
			four1Call(f,N,-1);
			for (int j = 0; j < N; j++)
				input[(i*N+j)*N+k] = f[j];
		}
	}

	// inverse FFT dim3 of rho
	for (int j = 0; j < N; j++) {
		for(int k=0;k<N;k++){
			for (int i = 0; i < N; i++)
				f[i] = input[(i*N+j)*N+k];
			four1Call(f,N,-1);
			for (int i = 0; i < N; i++)
				input[(i*N+j)*N+k] = f[i];
		}
	}

	//Normalizacija potrebna za ovakvu diskretnu furijeovu transformaciju
	for(int i=0; i<N;i++){
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++){
				float complex factor = 1/sqrt(N*N*N);
				input[(i*N+j)*N+k] *= factor;
			}
		}
	}

	return input;
}

// DFE FFT.
// Perform an FFT on the input data and return the output in a newly allocated buffer.
// The buffer will need to be deallocated using free().
static
float complex *fft_dfe(float complex *input)
{

	const size_t transfer_size = sizeof(float complex) * N * N * N;
	float complex *output = (float complex *)allocate_mem(transfer_size);
	max_file_t *maxfile = Poisson_init();
	max_engine_t *engine = max_load(maxfile, "*");

	max_actions_t* act = max_actions_init(maxfile, "MyInterface");

	max_set_param_uint64t(act,"N",N);

	max_queue_input(act,  "fft_in",  input,  transfer_size);
	max_queue_output(act, "fft_out", output, transfer_size);

	max_run(engine, act);
	max_unload(engine);

	return output;
}

int main(void)
{
	FILE *fileInput, *fileOutput, *dataFile;
			fileInput=fopen("input","r");
			fileOutput=fopen("output","w");
			dataFile=fopen("potential.data","w");

	float h=1/(float)N;

	printf("Poisson equation example - %d x %d x %d points\n", N, N, N);

	unsigned short *rng = seed48((unsigned short[]){1,1,1});
	float complex *input = (float complex *)allocate_mem(sizeof(float complex) * N * N * N);
	float complex *cpu_output;

	float q=-10.0;
	for (int j = 0; j < N * N * N; j++) {
		//float real, imag;
		//fscanf(fileInput,"%f",&real);
		//fscanf(fileInput,"%f",&imag);
		//input[j] = real + I*imag;

		input[j] = erand48(rng) + I*erand48(rng);

		fprintf(fileInput,"%.3f %.3f\n",creal(input[j]),cimag(input[j]));
	}

	double start_time, end_time;

	printf("Executing FFT on DFE...\n");
	start_time = timesec();
		float complex *dfe_output = fft_dfe(input);
	end_time = timesec();
	printf("DFE took %.6fs\n", end_time - start_time);


	printf("Executing FFT on CPU...\n");
	start_time = timesec();
		cpu_output = fft_cpu(input);
	end_time = timesec();
	printf("CPU took %.6fs\n", end_time - start_time);


	printf("Checking that output data matches...\n"); fflush(stdout);
	bool pass = true;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				int index = (i*N + j)*N+k;

				if (!crelequal(cpu_output[index], dfe_output[index],1e-2)) {
					printf("Mismatch @ [%5d][%5d][%5d]: cpu = [% .8f % .8f] dfe = [% .8f, % .8f]\n",
							i, j, k,
							crealf(cpu_output[index]),
							cimagf(cpu_output[index]),
							crealf(dfe_output[index]),
							cimagf(dfe_output[index])
					);
					pass = false;
				}
			}
		}
	}

	printf("Writing all output data to file: output\n");
	for(int j = 0; j < N * N * N; j++){
		fprintf(fileOutput,"[%d][%d][%d]\t(%.3f,%.3f)\n",j/(N*N),(j/N)%N,j%N,creal(dfe_output[j]),cimag(dfe_output[j]));
	}

	printf("Writing plottable data to file: potential.data\n");
	int i=N/2;
	for (int j = 0; j < N; j++) {
		double x = j * h;
		for (int k = 0; k < N; k++) {
			double y = k * h;
			fprintf(dataFile,"%.3f\t%.3f\t%.3f\n",x,y,creal(dfe_output[(i*N+j)*N+k]));
		}
		fprintf(dataFile,"\n");
	}

	printf("Cleaning up...\n");
	free(input);
	free(dfe_output);

	printf("Data check %s.\n", pass ? "passed" : "failed");
	return !pass;
}
