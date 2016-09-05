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
 * Randomly generate data or reads input from file.
 *
 * @param size 		Number of samples.
 * @param data 		Pointer to the array used to store the data.
 * @param type 		True if input should be read from file
 * @param fileName 	Name of the input file
 */
void generateTestData(const int size, float complex *data, int type, char *fileName) {
	srand(time(NULL));

	FILE *input = NULL;
	if(type!=0){
		 input = fopen(fileName,"r");
	}

	for (int i = 0; i < size; i++) {
		float real, imag;

		if(type != 0){ //input from file
			fscanf(input,"%f",&real);
			fscanf(input,"%f",&imag);
		}
		else{ //random input
			real = (float)rand()/(float)RAND_MAX * 10;
			imag = (float)rand()/(float)RAND_MAX * 10;
			int signReal = rand() % 2;
			int signImag = rand() % 2;

			real = signReal == 0 ? real : -real;
			imag = signImag == 0 ? imag : -imag;
		}

		data[i] = real + I * imag;
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
 * CPU code to calculate a 1D fft.
 * This function is not written to be fast. Only to be simple.
 *
 * @param n Size of the fft (has to be power of 2).
 * @param values Input samples. Also used to store the coefficients.
 */
void fftCPU(const int n, float complex* values) {
	if (n > 1) {
		float complex *g = (float complex*) mallocWrapper((n / 2) * sizeof(float complex));
		float complex *u = (float complex*) mallocWrapper((n / 2) * sizeof(float complex));


		// reorder data
		for (int i = 0; i < n / 2; i++) {
			g[i] = values[i * 2];
			u[i] = values[i * 2 + 1];
		}

		// calculate fft recursively
		fftCPU(n / 2, g);
		fftCPU(n / 2, u);

		// combine results
		for (int k = 0; k < n / 2; k++) {
			float complex expFactor = cexpf(-2.0 * M_PI * I * (float complex) k / (float complex) n);
			values[k]         = g[k] + u[k] * expFactor;
			values[k + n / 2] = g[k] - u[k] * expFactor;
		}

		// free allocated memory
		free(g);
		free(u);
	}
}

/**
 * CPU code to calculate a 1D ifft.
 * This function is not written to be fast. Only to be simple.
 *
 * @param n Size of the ifft (has to be power of 2).
 * @param values Input samples. Also used to store the coefficients.
 */
void ifftCPU(const int n, float complex* values) {
	if (n > 1) {
		float complex *g = (float complex*) mallocWrapper((n / 2) * sizeof(float complex));
		float complex *u = (float complex*) mallocWrapper((n / 2) * sizeof(float complex));

		// reorder data
		for (int i = 0; i < n / 2; i++) {
			g[i] = values[i * 2];
			u[i] = values[i * 2 + 1];
		}

		// calculate fft recursively
		ifftCPU(n / 2, g);
		ifftCPU(n / 2, u);

		// combine results
		for (int k = 0; k < n / 2; k++) {
			float complex expFactor = cexpf(+2.0 * M_PI * I * (float complex) k / (float complex) n);
			values[k] = g[k] + u[k] * expFactor;
			values[k + n / 2] = g[k] - u[k] * expFactor;
		}

		// free allocated memory
		free(g);
		free(u);
	}
}

/**
 * Function to transpose a 2D Array.
 *
 * @param firstDimension Size of the first dimension.
 * @param secondDimension Size of the second dimension.
 * @param data The array to transpose. It will be replaced by the transposed array.
 */
void transposeData(const int firstDimension, const int secondDimension, float complex*** data) {
	float complex** newData = mallocWrapper(firstDimension * sizeof(float complex*));
	for (int i = 0; i < firstDimension; i++) {
		newData[i] = mallocWrapper(secondDimension * sizeof(float complex));
		for (int j = 0; j < secondDimension; j++) {
			newData[i][j] = (*data)[j][i];
		}
	}

	// free space again
	for (int i = 0; i < secondDimension; i++) {
		free((*data)[i]);
	}
	free(*data);
	*data = newData;
}

/**
 * Wrapper function for the CPU implementation of the 1D fft.
 * In the 1D case the input data gets copied into a new array in order to do not change the input.
 * In the 2D case the data has to be reordered, transposed and multiple 1D ffts have to be executed.
 * The same principle applies for the 3D case
 *
 * @param size number of samples.
 * @param inputData Sample array.
 * @param expectedData Array for the coefficients.
 */
void fftCPUWrapper(const int size, float complex* inputData, float complex* expectedData) {
	if (MD == 1) { // 1D FFT
		memcpy((void*) expectedData, (void*) inputData, size * sizeof(float complex));
		fftCPU(size, expectedData);
	} else if (LD == 1) { // 2D FFT
		// To make things clearer again we copy our data into a two dimensional array
		float complex** inputData2D = mallocWrapper(MD * sizeof(float complex*));
		for (int i = 0; i < MD; i++) {
			inputData2D[i] = mallocWrapper(ND * sizeof(float complex));
			for (int j = 0; j < ND; j++) {
				inputData2D[i][j] = inputData[i * ND + j];
			}
		}

		// Calculate fft on each row
		for (int i = 0; i < MD; i++) {
			fftCPU(ND, inputData2D[i]);
		}

		// transpose
		transposeData(ND, MD, &inputData2D);

		// calculate fft on each column
		for (int i = 0; i < ND; i++) {
			fftCPU(MD, inputData2D[i]);
		}

		// transpose back
		transposeData(MD, ND, &inputData2D);

		for (int i = 0; i < MD; i++) {
			for (int j = 0; j < ND; j++) {
				expectedData[i * ND + j] = inputData2D[i][j];
			}
		}

		for (int i = 0; i < MD; i++) {
			free(inputData2D[i]);
		}
		free(inputData2D);
	} else { // 3D FFT
		// To make things clearer again we copy our data into a three dimensional Array
		float complex*** inputData3D = mallocWrapper(LD * sizeof(float complex**));
		for (int i = 0; i < LD; i++) {
			inputData3D[i] = mallocWrapper(MD * sizeof(float complex*));
			for (int j = 0; j < MD; j++) {
				inputData3D[i][j] = mallocWrapper(ND * sizeof(float complex));
				for (int k = 0; k < ND; k++) {
					inputData3D[i][j][k] = inputData[i * ND * MD + j * ND + k];
				}
			}
		}

		// Now we can first calculate the fft on each row
		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < MD; j++) {
				fftCPU(ND, inputData3D[i][j]);
			}
		}

		// Transpose N and M Dimension
		for (int i = 0; i < LD; i++) {
			transposeData(ND, MD, &(inputData3D[i]));
		}

		// Calculate fft on each column
		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < ND; j++) {
				fftCPU(MD, inputData3D[i][j]);
			}
		}

		// Transpose N and M Dimension back
		for (int i = 0; i < LD; i++) {
			transposeData(MD, ND, &(inputData3D[i]));
		}

		// Now do fft in the L Dimension
		float complex* buffer = mallocWrapper(LD * sizeof(float complex));
		for (int i = 0; i < MD; i++) {
			for (int j = 0; j < ND; j++) {
				for (int k = 0; k < LD; k++) {
					buffer[k] = inputData3D[k][i][j];
				}
				fftCPU(LD, buffer);
				for (int k = 0; k < LD; k++) {
					inputData3D[k][i][j] = buffer[k];
				}
			}
		}
		free(buffer);

		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < MD; j++) {
				for (int k = 0; k < ND; k++) {
					expectedData[i * ND * MD + j * ND + k] = inputData3D[i][j][k];
				}
			}
		}

		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < MD; j++) {
				free(inputData3D[i][j]);
			}
			free(inputData3D[i]);
		}
		free(inputData3D);
	}
}

/**
 * Wrapper function for the CPU implementation of the 1D ifft.
 * In the 1D case the input data gets copied into a new array in order to do not change the input.
 * In the 2D case the data has to be reordered, transposed and multiple 1D ffts have to be executed.
 * The same principle applies for the 3D case
 *
 * @param size number of samples.
 * @param inputData Sample array.
 * @param expectedData Array for the coefficients.
 */
void ifftCPUWrapper(const int size, float complex* inputData, float complex* expectedData) {
	if (MD == 1) { // 1D FFT
		memcpy((void*) expectedData, (void*) inputData, size * sizeof(float complex));
		ifftCPU(size, expectedData);
	} else if (LD == 1) { // 2D FFT
		// To make things clearer again we copy our data into a two dimensional array
		float complex** inputData2D = mallocWrapper(MD * sizeof(float complex*));
		for (int i = 0; i < MD; i++) {
			inputData2D[i] = mallocWrapper(ND * sizeof(float complex));
			for (int j = 0; j < ND; j++) {
				inputData2D[i][j] = inputData[i * ND + j];
			}
		}

		// Calculate ifft on each row
		for (int i = 0; i < MD; i++) {
			ifftCPU(ND, inputData2D[i]);
		}

		// transpose
		transposeData(ND, MD, &inputData2D);

		// calculate ifft on each column
		for (int i = 0; i < ND; i++) {
			ifftCPU(MD, inputData2D[i]);
		}

		// transpose back
		transposeData(MD, ND, &inputData2D);

		for (int i = 0; i < MD; i++) {
			for (int j = 0; j < ND; j++) {
				expectedData[i * ND + j] = inputData2D[i][j];
			}
		}

		for (int i = 0; i < MD; i++) {
			free(inputData2D[i]);
		}
		free(inputData2D);
	} else { // 3D IFFT
		// To make things clearer again we copy our data into a three dimensional Array
		float complex*** inputData3D = mallocWrapper(LD * sizeof(float complex**));
		for (int i = 0; i < LD; i++) {
			inputData3D[i] = mallocWrapper(MD * sizeof(float complex*));
			for (int j = 0; j < MD; j++) {
				inputData3D[i][j] = mallocWrapper(ND * sizeof(float complex));
				for (int k = 0; k < ND; k++) {
					inputData3D[i][j][k] = inputData[i * ND * MD + j * ND + k];
				}
			}
		}

		// Now we can first calculate the ifft on each row
		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < MD; j++) {
				ifftCPU(ND, inputData3D[i][j]);
			}
		}

		// Transpose N and M Dimension
		for (int i = 0; i < LD; i++) {
			transposeData(ND, MD, &(inputData3D[i]));
		}

		// Calculate ifft on each column
		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < ND; j++) {
				ifftCPU(MD, inputData3D[i][j]);
			}
		}

		// Transpose N and M Dimension back
		for (int i = 0; i < LD; i++) {
			transposeData(MD, ND, &(inputData3D[i]));
		}

		// Now do ifft in the L Dimension
		float complex* buffer = mallocWrapper(LD * sizeof(float complex));
		for (int i = 0; i < MD; i++) {
			for (int j = 0; j < ND; j++) {
				for (int k = 0; k < LD; k++) {
					buffer[k] = inputData3D[k][i][j];
				}
				ifftCPU(LD, buffer);
				for (int k = 0; k < LD; k++) {
					inputData3D[k][i][j] = buffer[k];
				}
			}
		}
		free(buffer);

		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < MD; j++) {
				for (int k = 0; k < ND; k++) {
					expectedData[i * ND * MD + j * ND + k] = inputData3D[i][j][k] / (LD * MD * ND);
				}
			}
		}

		for (int i = 0; i < LD; i++) {
			for (int j = 0; j < MD; j++) {
				free(inputData3D[i][j]);
			}
			free(inputData3D[i]);
		}
		free(inputData3D);
	}
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
	fftCPUWrapper(size, inputData, expectedData);
	poissonCPU(expectedData,contents, expectedData);
	ifftCPUWrapper(size, expectedData, expectedData);
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
	Poisson(size / 4, size / 4, size / 4, dh, input, size * sizeof(float complex), result, size * sizeof(float complex), (double*)contents);
}

/**
 * Writes all data to output file.
 *
 * @param size Number of samples
 * @param data Array of samples
 * @param fileName Name of the output file
 */
void writeAllToFile(const int size, const float complex *data, char *fileName){
	FILE *output = fopen(fileName, "w");

	printf("Writing all output data to file: %s.\n", fileName);
	for(int i = 0; i < size; i++){
		fprintf(output,"[%d][%d][%d]\t(%.3f,%.3f)\n",i/(size*size),(i/size)%size,i%size,creal(data[i]),cimag(data[i]));
	}
}

/**
 * Writes plottable data to output file.
 *
 * @param size Number of samples
 * @param data Array of samples
 * @param fileName Name of the output file
 */
void writePlottableToFile(const int size, const float complex *data, char *fileName){
	FILE *output = fopen(fileName, "w");
	float h=1/(float)ND;

	printf("Writing plottable data to file: %s.\n", fileName);
	int i=ND/2;
	for (int j = 0; j < ND; j++) {
		double x = j * h;
		for (int k = 0; k < ND; k++) {
			double y = k * h;
			fprintf(output,"%.3f\t%.3f\t%.3f\n",x,y,creal(data[(i*ND+j)*ND+k]));
		}
		fprintf(output,"\n");
	}
}

int main(void)
{
	const int size = ND * MD * LD;
	float complex* inputData = mallocWrapper(size * sizeof(float complex));
	float complex* expectedData = mallocWrapper(size * sizeof(float complex));
	float complex* resultData = mallocWrapper(size * sizeof(float complex));
	double complex* contents = create_twiddles();

	generateTestData(size, inputData, 1, "samples/input_pointCharge");

	poissonCPUWrapper(size, inputData, contents, expectedData);

	poissonDFE(size, inputData, contents, resultData);

	int status = check(size, expectedData, resultData);
	if(status != 0)
		printf("Test failed!\n");
	else
		printf("Test passed!\n");

	writeAllToFile(size, resultData, "output/all_data");

	writePlottableToFile(size, resultData, "output/potential.data");

	printf("Done.\n");

	free(inputData);
	free(expectedData);
	free(resultData);

	return status;
}
