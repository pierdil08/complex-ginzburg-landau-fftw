#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * laplacian function computes the laplacian using the pseudospectral method. When dubugging, it'll output A.out once. This is used with the cgldebug.m script to analyze results. This is computed after the after the forward FFT, before applying the laplacian.
 */

void laplacian(fftw_complex *A_raw, fftw_complex *lapA_raw,fftw_plan *p_forward, fftw_plan *p_inverse double *k2, int N) {

    int NN = N * N;

    // Allocate fresh memory for temporary step
    fftw_complex *Atempij = fftw_malloc(sizeof(fftw_complex) * NN);
    fftw_complex *A_fft = fftw_malloc(sizeof(fftw_complex) * NN);

    // Cast to double* for real/imag indexing - easier for math
    double *A_r = (double *)A_raw;
    double *Atemp_r = (double *)Atempij;
    double *Afft_r = (double *)A_fft;
    double *lapA_r = (double *)lapA_raw;

    // Explicit copy: A_raw -> Atempij
    for (int i = 0; i < NN; i++) {
        Atemp_r[2*i]     = A_r[2*i];     // real part
        Atemp_r[2*i + 1] = A_r[2*i + 1]; // imag part
    }

    // Forward FFT
    fftw_execute(p_forward);

    // Apply Laplacian operator: multiply each Fourier component by -k^2
    for (int i = 0; i < NN; i++) {
        Afft_r[2*i]     *= k2[i]; // real part scaled
        Afft_r[2*i + 1] *= k2[i]; // imag part scaled
    }

    // Inverse FFT
    fftw_execute(p_inverse);

    // Normalize and explicitly copy result into lapA_raw
    for (int i = 0; i < NN; i++) {
        lapA_r[2*i]     = Atemp_r[2*i]     / NN; // real part
        lapA_r[2*i + 1] = Atemp_r[2*i + 1] / NN; // imag part
    }

    // Free temporary FFT plans and arrays
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_inverse);
    fftw_free(Atempij);
    fftw_free(A_fft);
}

/*
int main(int argc, char* argv[])

Inputs:
argv[1] = N - grid points per dimension
argv[2] = c1
argv[3] = c3
argv[4] = M - # time steps

Output:
cgl.out - contains 11 data points for A. Those being at 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 steps.
*/

int main(int argc, char* argv[]) {
    clock_t start = clock();

    // input arguments
    int argi = 0;
    int N = atoi(argv[++argi]);
    double c1 = atof(argv[++argi]);
    double c3 = atof(argv[++argi]);
    double M = atof(argv[++argi]);
    double T = 10000.0;
    double dt = T / M;
    double L = 128 * M_PI;
//print input arguments
    printf("N = %d, c1 = %f, c3 = %f, T = %f\n", N, c1, c3, T);

    //if no seed specified at argv[5] then generate seed using seconds since Jan 1st 1970. 
    long int seed;
    if (argi < argc - 1) {
        seed = atol(argv[++argi]);
    } else {
        seed = (long int) time(NULL);
    }
    srand48(seed);
    printf("seed = %ld\n", seed);

    // Allocate memory for fftw_complex type. this type is need to do the forward and back fourier transforms. 
    fftw_complex *A_raw = fftw_malloc(N * N * sizeof(fftw_complex));
    fftw_complex *lapA_raw = fftw_malloc(N * N * sizeof(fftw_complex));
    fftw_complex *a1_raw = fftw_malloc(N * N * sizeof(fftw_complex));
    fftw_complex *a2_raw = fftw_malloc(N * N * sizeof(fftw_complex));
    fftw_complex *a3_raw = fftw_malloc(N * N * sizeof(fftw_complex));
    //Allocate memory for complex type. This type is better for calculations. For example, cabs couldn't be used on fftw_complex type.
    double complex *A = (double complex *)A_raw;
    double complex *lapA = (double complex *)lapA_raw;
    double complex *a1 = (double complex *)a1_raw;
    double complex *a2 = (double complex *)a2_raw;
    double complex *a3 = (double complex *)a3_raw;

// Create FFT plans
    fftw_plan p_forward = fftw_plan_dft_2d(N, N, Atempij, A_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan p_inverse = fftw_plan_dft_2d(N, N, A_fft, Atempij, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Initialize A with random values
    for (int i = 0; i < N * N; i++) {
        A[i] = (drand48() * 3.0 - 1.5) + I * (drand48() * 3.0 - 1.5);
    }

    // Save initial A to file
    FILE* fid_init = fopen("CGL.out", "w");
    fwrite(A, sizeof(fftw_complex), N * N, fid_init);
    fclose(fid_init);

    // Precompute k^2
    double *k2 = malloc(sizeof(double) * N * N);
    for (int i = 0; i < N; i++) {
        int ky = (i <= N/2) ? i : i - N;
        for (int j = 0; j < N; j++) {
            int kx = (j <= N/2) ? j : j - N;
            k2[i * N + j] = -(double)(kx*kx + ky*ky);
        }
    }
    //Here the coefficient is calculated. We multiply instead of using the pow function as it is more efficient.
    double coeff = (2.0 * M_PI / L) * (2.0 * M_PI / L);
    
    //establish file 
    FILE* fileid = fopen("cgl.out", "w");

    // Time stepping loop
    for (int step = 0; step < M; ++step) {
        if (step % ((int)(M/10)) == 0) {
            fwrite(A, sizeof(fftw_complex), N*N, fileid);
        }

        // Stage 1 of RK4
        laplacian(A_raw, lapA_raw,p_forward, p_inverse k2, N);
        for (int i = 0; i < N * N; i++) {
            double complex rhs = A[i] + coeff * (1.0 + I * c1) * lapA[i] - (1.0 - I * c3) * cabs(A[i]) * cabs(A[i]) * A[i];
            a1[i] = A[i] + dt / 4.0 * rhs;
        }
        //FILE* fid1 = fopen("stage1.out", "w");
        //fwrite(a1_raw, sizeof(fftw_complex), N*N, fid1);
        //fclose(fid1);

        // Stage 2
        laplacian((fftw_complex *)a1, lapA_raw,p_forward, p_inverse k2, N);
        for (int i = 0; i < N * N; i++) {
            double complex rhs = a1[i] + coeff * (1.0 + I * c1) * lapA[i] - (1.0 - I * c3) * cabs(a1[i]) * cabs(a1[i]) * a1[i];
            a2[i] = A[i] + dt / 3.0 * rhs;
        }
        //FILE* fid2 = fopen("stage2.out", "w");
        //fwrite(a2, sizeof(fftw_complex), N*N, fid2);
        //fclose(fid2);

        // Stage 3
        laplacian((fftw_complex *)a2, lapA_raw,p_forward, p_inverse, k2, N);
        for (int i = 0; i < N * N; i++) {
            double complex rhs = a2[i] + coeff * (1.0 + I * c1) * lapA[i] - (1.0 - I * c3) * cabs(a2[i]) * cabs(a2[i]) * a2[i];
            a3[i] = A[i] + dt / 2.0 * rhs;
        }
       // FILE* fid3 = fopen("stage3.out", "w");
       // fwrite(a3, sizeof(fftw_complex), N*N, fid3);
       // fclose(fid3);

        // Stage 4
        laplacian((fftw_complex *)a3, lapA_raw,p_forward, p_inverse, k2, N);
        for (int i = 0; i < N * N; i++) {
            double complex rhs = a3[i] + coeff * (1.0 + I * c1) * lapA[i] - (1.0 - I * c3) * cabs(a3[i]) * cabs(a3[i]) * a3[i];
            A[i] = A[i] + dt * rhs;
        }
       // FILE* fid4 = fopen("stage4.out", "w");
       // fwrite(A, sizeof(fftw_complex), N*N, fid4);
       // fclose(fid4);

        //exit(0); // Stop after one step for debugging
    }
    //Write the final time step for the data
    fwrite(A, sizeof(fftw_complex), N*N, fileid);
    fclose(fileid);

    // Free memory
    fftw_free(A_raw);
    fftw_free(lapA_raw);
    fftw_free(a1_raw);
    fftw_free(a2_raw);
    fftw_free(a3_raw);
    free(k2);
    
    //Print the total time that the simulation took. 
    printf("Elapsed time: %g seconds\n", (float)(clock() - start) / CLOCKS_PER_SEC);
    return 0;
}

