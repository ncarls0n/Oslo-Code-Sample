#include <stdio.h>
#include <stdlib.h>

// Function to count the number of lines in the file, to determine the size of the data array
long countLines(const char *filePath) {
    FILE *file = fopen(filePath, "r");
    if (!file) {
        perror("Failed to open file for line counting");
        exit(EXIT_FAILURE);
    }
    long lines = 0;
    int ch;
    while (EOF != (ch = getc(file))) {
        if (ch == '\n')
            lines++;
    }
    fclose(file);
    return lines;
}

int main() {
    const char *inputFilePath = "rhog_out.txt"; // Input text file path
    const char *outputFilePath = "output_binary.bin"; // Output binary file path
    FILE *inputFile, *outputFile;
    long totalElements;
    
    // Count the number of lines/data elements in the input file
    totalElements = countLines(inputFilePath);
    
    // Allocate memory for the data array
    float *data = (float *)malloc(totalElements * sizeof(float));
    if (!data) {
        perror("Memory allocation failed");
        return EXIT_FAILURE;
    }

    // Open the input text file for reading
    inputFile = fopen(inputFilePath, "r");
    if (!inputFile) {
        perror("Failed to open input file");
        free(data);
        return EXIT_FAILURE;
    }

    // Read data from text file into the array
    for (long i = 0; i < totalElements; i++) {
        if (fscanf(inputFile, "%f", &data[i]) != 1) {
            fprintf(stderr, "Failed to read data at line %ld\n", i + 1);
            fclose(inputFile);
            free(data);
            return EXIT_FAILURE;
        }
    }
    fclose(inputFile);

    // Open the output binary file for writing
    outputFile = fopen(outputFilePath, "wb");
    if (!outputFile) {
        perror("Failed to open output file");
        free(data);
        return EXIT_FAILURE;
    }

    // Write the array to the binary file
    fwrite(data, sizeof(float), totalElements, outputFile);
    fclose(outputFile);
    free(data);

    printf("Conversion to binary file complete.\n");
    return EXIT_SUCCESS;
}
