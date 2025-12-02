#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to trim whitespace from the start and end of a string
char* trimWhitespace(char* str);

// Function to load configurations from a file
void loadConfig(const char *configFilePath, char **inputFilePath, char **outputFilePath, int *ntiles, int *nmesh, int *nbuff); 
