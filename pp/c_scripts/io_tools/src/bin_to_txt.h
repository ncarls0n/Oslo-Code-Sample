#ifndef BINARY_READER_H
#define BINARY_READER_H

/**
 * Reads binary data from a specified file and writes it in a human-readable format
 * to another file.
 * 
 * @param inputFilePath Path to the input binary file.
 * @param outputFilePath Path to the output text file.
 * @param n Dimension of the square matrix in one plane.
 * @param local_nz Number of planes in the 3D matrix.
 * @return 0 on success, -1 on failure.
 */
int readBinaryData(const char* inputFilePath, const char* outputFilePath, int n, int local_nz);

#endif 
