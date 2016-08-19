#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MaxSLiCInterface.h"
#include "Maxfiles.h"

#define ND Poisson_N
#define MD Poisson_M
#define LD Poisson_L


/**
 * Wrapper for malloc in order to check if malloc returns NULL and exit if this happens.
 *
 * @param size Size of the memory to allocate.
 * @return pointer to the allocated memory.
 */
void * mallocWrapper(size_t size) {
	void * pointer = malloc(size);
	if (pointer == NULL) {
		fprintf(stderr, "Was not able to allocate enough memory");
		exit(-1);
	}
	return pointer;
}

/**
 * Generates twiddle factors needed for Poisson solver
 *
 * @return pointer to allocated and filled memory
 */
double complex * create_twiddles() {
	double pi = 4 * atan(1.0);
	double complex* contents = mallocWrapper(ND * sizeof(double complex));;
	for (int b = 0; b <2*ND; b+=2){
		contents[b/2]   = cos(pi/ND*b) + I * sin(pi/ND*b);
	}

	return contents;
}

/**
 * Compares the result calculated by the dfe with the result of the cpu version.
 * The signal to noise ratio gets calculated and if it is to low the actual results as well
 * as the expected results get printed.
 *
 * @param size Number of samples.
 * @param expected Expected result.
 * @param result Actual result.
 * @return 0 if SNR is ok. 1 if not.
 */
int check(const int size, const float complex *expected, const float complex *result)
{
	// calculate SNR
	double S = 0.0;
	double N = 0.0;

	for (int i = 0; i < size; i++) {
		float complex res = result[i];
		float complex exp = expected[i];

		float complex err = exp - res;

		S += (creal(exp) * creal(exp)) + (cimag(exp) * cimag(exp));
		N += (creal(err) * creal(err)) + (cimag(err) * cimag(err));
	}

	double SNR = (10 * log10(S / N));
	printf("SNR: %f\n", SNR);
	// if result is wrong print everything
	if (SNR < 69.0) {
		for (int i = 0; i < size; i++) {
			float complex res = result[i];
			float complex exp = expected[i];

			printf("Index %d: Is: %f + %f * i\tExpected: %f + %f * i\n", i, creal(res), cimag(res),
					creal(exp), cimag(exp));
		}
	}

	return SNR < 69.0;
}

/**
 * randomly generate data.
 *
 * @param size Number of samples.
 * @param data Pointer to the array used to store the data.
 */
void generateTestData(const int size, float complex *data) {
	srand(time(NULL));
	for (int i = 0; i < size; i++) {
		float real = (float)rand()/(float)RAND_MAX * 10;
		float imag = (float)rand()/(float)RAND_MAX * 10;
		int signReal = rand() % 2;
		int signImag = rand() % 2;

		real = signReal == 0 ? real : -real;
		imag = signImag == 0 ? imag : -imag;

		data[i] = real + I * imag;
		data[i] = i;
	}
}

/**
 * Get time in seconds
 */
double timesec() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return ((double) (t.tv_sec * 1e6L + t.tv_usec)) / 1.0e6;
}

/**
 * CPU code to solve Poisson equation.
 *
 * @param values Input samples.
 */
void poissonCPU(float complex* input, double complex* contents, float complex* expected) {
	float h=1/(float)ND;
	float complex Wl=1.0, Wm = 1.0, Wn = 1.0;
	for(int l=0; l<ND; l++){
		Wl = contents[l];
		for (int m = 0; m < ND; m++) {
			Wm = contents[m];
			for (int n = 0; n < ND; n++) {
				Wn = contents[n];
				float complex denom = Wm + 1.0 / Wm + Wn + 1.0 / Wn + Wl + 1.0 / Wl - 6.0;
				if (denom != 0.0)
					expected[(l*ND+m)*ND+n] = input[(l*ND+m)*ND+n] * h * h / denom;
			}
		}
	}
}

/**
 *
 */
void poissonCPUWrapper(const int size, float complex* inputData, double complex* contents, float complex* expectedData) {
	//fftCPUWrapper(size, inputData, expectedData);
	poissonCPU(inputData,contents,expectedData);
	//ifftCPUWrapper(size, expectedData, expectedData);
}

/**
 * Code to solve Poisson equation on DFE.
 *
 * @param size Number of samples.
 * @param input Array of samples.
 * @param result Array for the coefficients.
 */
void poissonDFE(const int size, float complex* input, double complex* contents, float complex* result) {
	printf("Running on DFE.\n");
	float h=1/(float)ND;
	float dh = h * h;
	Poisson(size / 4, dh, input, size * sizeof(float complex), result, size * sizeof(float complex), (double*)contents);
}

int main(void)
{
	const int size = ND * MD * LD;
	float complex* inputData = mallocWrapper(size * sizeof(float complex));
	float complex* expectedData = mallocWrapper(size * sizeof(float complex));
	float complex* resultData = mallocWrapper(size * sizeof(float complex));
	double complex* contents = create_twiddles();

	generateTestData(size, inputData);

	poissonCPUWrapper(size, inputData, contents, expectedData);

	poissonDFE(size, inputData, contents, resultData);

	int status = check(size, expectedData, resultData);
	if(status != 0)
		printf("Test failed!\n");
	else
		printf("Test passed!\n");

	printf("Done.\n");

	free(inputData);
	free(expectedData);
	free(resultData);
	return status;
}
