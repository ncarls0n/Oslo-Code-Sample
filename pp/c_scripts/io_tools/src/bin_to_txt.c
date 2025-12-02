#include <stdio.h>
#include <stdlib.h>
#include "bin_to_txt.h"
#include "read_config.h"

int main() {
    FILE *inputFile, *outputFile; // File pointers
    char *inputFilePath = NULL, *outputFilePath = NULL;
    int ntiles, nmesh, nbuff;
    long fileSize;
    float *data, *trimmedData;
    int i, j, k;

    // Load configurations from the file BEFORE using these values
    loadConfig("../configs/config.txt", &inputFilePath, &outputFilePath, &ntiles, &nmesh, &nbuff);

    // Now that ntiles, nmesh, and nbuff have been loaded, calculate boxSize
    int boxSize = ntiles * (nmesh - 2 * nbuff) + 2 * nbuff; // Size of the box in cells Number

    // Use the variables for further processing...
    printf("Input File Path: %s\n", inputFilePath);
    printf("Output File Path: %s\n", outputFilePath);
    printf("ntiles: %d, nmesh: %d, nbuff: %d\n", ntiles, nmesh, nbuff);

    // Open input binary file
    inputFile = fopen(inputFilePath, "rb");
    if (!inputFile) {
        perror("Failed to open input file");
        return EXIT_FAILURE;
    }

    // Get file size
    fseek(inputFile, 0L, SEEK_END);
    fileSize = ftell(inputFile);
    rewind(inputFile);

    // Allocate memory for data
    int totalElements = fileSize / sizeof(float);
    data = (float *)malloc(totalElements * sizeof(float));
    if (!data) {
        perror("Memory allocation failed");
        fclose(inputFile);
        return EXIT_FAILURE;
    }

    // Read data from file
    fread(data, sizeof(float), totalElements, inputFile);
    fclose(inputFile);

    // Assuming boxSize*boxSize*boxSize equals totalElements
    // Allocate memory for trimmedData if trimming is needed
    //int trimmedSize = (boxSize - 2 * nbuff) * (boxSize - 2 * nbuff) * (boxSize - 2 * nbuff);
    int trimmedSize = boxSize * boxSize * boxSize;
    trimmedData = (float *)malloc(trimmedSize * sizeof(float));
    if (!trimmedData) {
        perror("Memory allocation for trimmed data failed");
        free(data);
        return EXIT_FAILURE;
    }

    // Copy data to trimmedData with trimming
    int index = 0;
    for (k = 0; k < boxSize; ++k) {
        for (j = 0; j < boxSize; ++j) {
            for (i = 0; i < boxSize; ++i) {
                // Ensure this indexing does not go out of bounds of the original data array
                if ((k * boxSize * boxSize) + (j * boxSize) + i < totalElements) {
                    trimmedData[index++] = data[(k * boxSize * boxSize) + (j * boxSize) + i];
                }
            }
        }
    }

    // Write trimmed data to output file
    outputFile = fopen(outputFilePath, "w");
    if (!outputFile) {
        perror("Failed to open output file");
        free(data);
        free(trimmedData);
        return EXIT_FAILURE;
    }

    for (i = 0; i < trimmedSize; ++i) {
        fprintf(outputFile, "%f\n", trimmedData[i]);
    }

    fclose(outputFile);
    free(data);
    free(trimmedData);

    printf("Data conversion complete.\n");
    return EXIT_SUCCESS;
}
