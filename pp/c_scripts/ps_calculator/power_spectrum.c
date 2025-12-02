#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <omp.h>

int get_field_size(const char *filename) {
    /* Get the size of the field from the file */
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
      perror("Failed to open file");
      return -1;
    }

    // Get the file size
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    fclose(file);

    // Get the size of the mesh
    float size_in_float = (float) size;
    float Nfloats  = size_in_float / 4;
    float nmesh  = cbrt(Nfloats); // Get mesh size as nmesh^3 = Nfloats
    int int_nmesh = (int) nmesh;
    if (floor(nmesh) == nmesh) {
      printf("mesh size is %d\n", int_nmesh);
      return int_nmesh;
    } else {
      int int_Nfloats = (int) Nfloats;
      printf("Incorrectly formatted binary file: expect a size of 4*n^3 bytes (n^3 floats), but received 4*%d\n", int_Nfloats);
      exit(1);
    }
}

int extract_bin_data(const char *inputFilePath, double **field, int ntiles, int nmesh, int nbuff) {
    FILE *inputFile; // Output file handling removed since not writing to a file
    //char *inputFilePath = NULL;
    long fileSize;
    float *data, *trimmedData;
    int i, j, k;

    int boxSize = ntiles * (nmesh - 2 * nbuff) + 2 * nbuff;

    printf("Input File Path: %s\n", inputFilePath);
    printf("ntiles: %d, nmesh: %d, nbuff: %d\n", ntiles, nmesh, nbuff);

    inputFile = fopen(inputFilePath, "rb");
    if (!inputFile) {
        perror("Failed to open input file");
        return EXIT_FAILURE;
    }

    fseek(inputFile, 0L, SEEK_END);
    fileSize = ftell(inputFile);
    rewind(inputFile);

    int totalElements = fileSize / sizeof(float);
    data = (float *)malloc(totalElements * sizeof(float));
    if (!data) {
        perror("Memory allocation failed");
        fclose(inputFile);
        return EXIT_FAILURE;
    }

    fread(data, sizeof(float), totalElements, inputFile);
    fclose(inputFile);

    int trimmedSize = boxSize * boxSize * boxSize;
    trimmedData = (float *)malloc(trimmedSize * sizeof(float));
    if (!trimmedData) {
        perror("Memory allocation for trimmed data failed");
        free(data);
        return EXIT_FAILURE;
    }

    int index = 0;
    for (k = 0; k < boxSize; ++k) {
        for (j = 0; j < boxSize; ++j) {
            for (i = 0; i < boxSize; ++i) {
                if ((k * boxSize * boxSize) + (j * boxSize) + i < totalElements) {
                    trimmedData[index++] = data[(k * boxSize * boxSize) + (j * boxSize) + i];
                }
            }
        }
    }
    free(data); // Original data no longer needed

    // Allocate and convert to double
    //*field = (double *)malloc(trimmedSize * sizeof(double));
    //if (!(*field)) {
    //    perror("Memory allocation for double array failed");
    //    free(trimmedData);
    //    return EXIT_FAILURE;
    //}

    for (int idx = 0; idx < trimmedSize; ++idx) {
        (*field)[idx] = (double) trimmedData[idx];
    }
    free(trimmedData); // Free trimmedData as it's no longer needed

    printf("Data conversion complete.\n");
    return EXIT_SUCCESS;
}

void read_density_field(const char *filename, double *field, int n) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n * n * n; i++) {
        if (fscanf(file, "%lf", &field[i]) != 1) {
            perror("Error reading file");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

void compute_power_spectrum_and_write_csv(double *field, int n, double box_size, const char *output_filename) {
    fftw_complex *out;
    fftw_plan p;
    int N = n * n * n;
    double doubleN = (double) N;
    double power;
    int ii, jj, kk;
    int num_threads = 40;

    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * doubleN);
    p = fftw_plan_dft_r2c_3d(n, n, n, field, out, FFTW_ESTIMATE);
    fftw_execute(p);

    int nyquist = n / 2 + 1;
    int klen = sqrt(3.0) * nyquist;
    double *P_k = calloc(klen + 1, sizeof(double));
    double *P_0 = calloc(klen + 1, sizeof(double));

    double dk = 2.0 * M_PI / box_size;
    double weights[2];
    int l;

    omp_set_num_threads(num_threads); 
    #pragma omp parallel for private(ii, jj, kk, l, weights, power) reduction(+:P_k[:klen+1], P_0[:klen+1])
    for (int k = 0; k < n; k++) {
        int kk = (k <= nyquist) ? k : n - k;
        for (int i = 0; i < n; i++) {
            int ii = (i <= nyquist) ? i : n - i;
            for (int j = 0; j < n; j++) {
                int jj = (j <= nyquist) ? j : n - j;
                double mode = sqrt(ii*ii + jj*jj + kk*kk);
                int l = floor(mode);
                double delta = mode - l;
                weights[0] = (1.0 - delta) * (1.0 - delta);
                weights[1] = delta * delta;
    
                if (l < klen) {
                    int index = i*n*(n/2+1) + j*(n/2+1) + (k <= nyquist ? k : nyquist - (k - nyquist));
                    power = out[index][0] * out[index][0] + out[index][1] * out[index][1];
    
                    P_k[l] += weights[0] * power;
                    P_0[l] += weights[0];
                    P_k[l+1] += weights[1] * power;
                    P_0[l+1] += weights[1];
                }
            }
        }
    }

    /*
    for (int i = 0; i < n; i++) {
        int ii = (i <= nyquist) ? i : n - i;
        for (int j = 0; j < n; j++) {
            int jj = (j <= nyquist) ? j : n - j;
            for (int k = 0; k <= nyquist; k++) {
                int kk = k;  // only positive k-values needed due to symmetry in FFT
                double mode = sqrt(ii*ii + jj*jj + kk*kk);
                l = floor(mode);
                double delta = mode - l;
                weights[0] = (1.0 - delta) * (1.0 - delta);
                weights[1] = delta * delta;

                if (l < klen) {
                    power = out[i*n*(n/2+1) + j*(n/2+1) + k][0] * out[i*n*(n/2+1) + j*(n/2+1) + k][0] +
                                   out[i*n*(n/2+1) + j*(n/2+1) + k][1] * out[i*n*(n/2+1) + j*(n/2+1) + k][1];

                    P_k[l] += weights[0] * power;
                    P_0[l] += weights[0];
                    P_k[l+1] += weights[1] * power;
                    P_0[l+1] += weights[1];
                }
            }
        }
    }
    */

    FILE *fp = fopen(output_filename, "w");
    if (fp == NULL) {
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "k,Power\n");

    for (int i = 0; i <= klen; i++) {
        if (P_0[i] > 0) {
            double bin_k = i * dk;
            fprintf(fp, "%g,%g\n", bin_k, P_k[i] / P_0[i]);
        }
    }

    free(P_k);
    free(P_0);
    fftw_destroy_plan(p);
    fftw_free(out);
    fclose(fp);
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <filename> <ntiles> <nbuff> <box_size> <output_csv>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char *filename = argv[1];
    int ntiles = atoi(argv[2]);
    int nbuff = atoi(argv[3]);
    double box_size = atof(argv[4]);
    const char *output_filename = argv[5];

    // Get the size of the mesh
    int nmesh = get_field_size(filename);

    // Allocate memory for density field
    double *field = (double *)malloc(nmesh * nmesh * nmesh * sizeof(double));
    if (field == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Can either read the data from txt or from the bin file
    //read_density_field(filename, field, n);
    extract_bin_data(filename, &field, ntiles, nmesh, nbuff);

    // Compute power spectrum and write to CSV
    compute_power_spectrum_and_write_csv(field, nmesh, box_size, output_filename);

    free(field);
    return 0;
}

