#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_config.h"

// Function to trim whitespace from the start and end of a string
char* trimWhitespace(char* str) {
    char* end;

    // Trim leading space
    while(isspace((unsigned char)*str)) str++;

    // All spaces?
    if(*str == 0)
        return str;

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator
    *(end+1) = 0;

    return str;
}

// Function to load configurations from a file
void loadConfig(const char *configFilePath, char **inputFilePath, char **outputFilePath, int *ntiles, int *nmesh, int *nbuff) {
    FILE *file = fopen(configFilePath, "r");
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    if (file == NULL) {
        perror("Could not open config file");
        exit(EXIT_FAILURE);
    }

    while ((read = getline(&line, &len, file)) != -1) {
        // Splitting the line into key and value
        char *token = strtok(line, "=");
        char *key = trimWhitespace(token);
        token = strtok(NULL, "=");
        char *value = trimWhitespace(token);

        // Assigning values to variables based on the key
        if (strcmp(key, "inputFilePath") == 0) {
            *inputFilePath = realloc(*inputFilePath, strlen(value) + 1);
            strcpy(*inputFilePath, value);
        } else if (strcmp(key, "outputFilePath") == 0) {
            *outputFilePath = realloc(*outputFilePath, strlen(value) + 1);
            strcpy(*outputFilePath, value);
        } else if (strcmp(key, "ntiles") == 0) {
            *ntiles = atoi(value);
        } else if (strcmp(key, "nmesh") == 0) {
            *nmesh = atoi(value);
        } else if (strcmp(key, "nbuff") == 0) {
            *nbuff = atoi(value);
        }
    }

    fclose(file);
    if (line)
        free(line);
}
