#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
    //Switches
    int hpkvd_params, compile_hpkvd, create_filterbank, merge_params, compile_merge;
    int map_params, compile_maps, batch, submit;

    //ClusterInfo
    char machine[50], submit_command[50];

    //PeakPatchParams
    int seed;
    char run_name[100], short_name[50], runtype[10];
    double boxsize;
    int nmesh, nbuff, ntile, largerun;
    char tlimit[10];
    int nnodes, tpnode, ntasks, ncpus, nompth;
    int ievol, num_redshifts;
    double maximum_redshift, global_redshift;
    int ilpt, ioutfield, ireadfield, iwant_field_part;
    char fielddir[100], densfilein[100], densfileout[100];
    int NonGauss;
    double fNL, A_nG, B_nG, R_nG, m_phi, m_chi, phi_w, phi_p, vev, m_tach, a_e;

    //WebskyParams
    char maps[100];
    int nside_map, npix_map;
    double fov_map, zmin_map, zmax_map;
    char tabfile_map[100], tabfile_sfr[100];
    int model_map, scramble_map, center_map, chihview_map, PSZcut_map, ellmax;

    //Cosmology
    double Omx, OmB, Omvac, h, ns, As, sigma8, tau, mnu;
    char pkfile[100];

} Config;

void initializeConfig(Config *config) {
    memset(config, 0, sizeof(Config));
}



void set_config_value(Config *config, const char *key, const char *value) {
    if (strcmp(key, "hpkvd_params") == 0) {
        config->hpkvd_params = atoi(value);
    } else if (strcmp(key, "compile_hpkvd") == 0) {
        config->compile_hpkvd = atoi(value);
    } else if (strcmp(key, "create_filterbank") == 0) {
        config->create_filterbank = atoi(value);
    } else if (strcmp(key, "merge_params") == 0) {
        config->merge_params = atoi(value);
    } else if (strcmp(key, "compile_merge") == 0) {
        config->compile_merge = atoi(value);
    } else if (strcmp(key, "map_params") == 0) {
        config->map_params = atoi(value);
    } else if (strcmp(key, "compile_maps") == 0) {
        config->compile_maps = atoi(value);
    } else if (strcmp(key, "batch") == 0) {
        config->batch = atoi(value);
    } else if (strcmp(key, "submit") == 0) {
        config->submit = atoi(value);
    } else if (strcmp(key, "machine") == 0) {
        strncpy(config->machine, value, sizeof(config->machine) - 1);
    } else if (strcmp(key, "submit_command") == 0) {
        strncpy(config->submit_command, value, sizeof(config->submit_command) - 1);
    } else if (strcmp(key, "seed") == 0) {
        config->seed = atoi(value);
    } else if (strcmp(key, "run_name") == 0) {
        strncpy(config->run_name, value, sizeof(config->run_name) - 1);
    } else if (strcmp(key, "short_name") == 0) {
        strncpy(config->short_name, value, sizeof(config->short_name) - 1);
    } else if (strcmp(key, "runtype") == 0) {
        strncpy(config->runtype, value, sizeof(config->runtype) - 1);
    } else if (strcmp(key, "boxsize") == 0) {
        config->boxsize = atof(value);
    } else if (strcmp(key, "nmesh") == 0) {
        config->nmesh = atoi(value);
    } else if (strcmp(key, "nbuff") == 0) {
        config->nbuff = atoi(value);
    } else if (strcmp(key, "ntile") == 0) {
        config->ntile = atoi(value);
    } else if (strcmp(key, "largerun") == 0) {
        config->largerun = atoi(value);
    } else if (strcmp(key, "tlimit") == 0) {
        strncpy(config->tlimit, value, sizeof(config->tlimit) - 1);
    } else if (strcmp(key, "nnodes") == 0) {
        config->nnodes = atoi(value);
    } else if (strcmp(key, "tpnode") == 0) {
        config->tpnode = atoi(value);
    } else if (strcmp(key, "ntasks") == 0) {
        config->ntasks = atoi(value);
    } else if (strcmp(key, "ncpus") == 0) {
        config->ncpus = atoi(value);
    } else if (strcmp(key, "nompth") == 0) {
        config->nompth = atoi(value);
    } else if (strcmp(key, "ievol") == 0) {
        config->ievol = atoi(value);
    } else if (strcmp(key, "num_redshifts") == 0) {
        config->num_redshifts = atoi(value);
    } else if (strcmp(key, "maximum_redshift") == 0) {
        config->maximum_redshift = atof(value);
    } else if (strcmp(key, "global_redshift") == 0) {
        config->global_redshift = atof(value);
    } else if (strcmp(key, "ilpt") == 0) {
        config->ilpt = atoi(value);
    } else if (strcmp(key, "ioutfield") == 0) {
        config->ioutfield = atoi(value);
    } else if (strcmp(key, "ireadfield") == 0) {
        config->ireadfield = atoi(value);
    } else if (strcmp(key, "iwant_field_part") == 0) {
        config->iwant_field_part = atoi(value);
    } else if (strcmp(key, "fielddir") == 0) {
        strncpy(config->fielddir, value, sizeof(config->fielddir) - 1);
    } else if (strcmp(key, "densfilein") == 0) {
        strncpy(config->densfilein, value, sizeof(config->densfilein) - 1);
    } else if (strcmp(key, "densfileout") == 0) {
        strncpy(config->densfileout, value, sizeof(config->densfileout) - 1);
    } else if (strcmp(key, "NonGauss") == 0) {
        config->NonGauss = atoi(value);
    } else if (strcmp(key, "fNL") == 0) {
        config->fNL = atof(value);
    } else if (strcmp(key, "A_nG") == 0) {
        config->A_nG = atof(value);
    } else if (strcmp(key, "B_nG") == 0) {
        config->B_nG = atof(value);
    } else if (strcmp(key, "R_nG") == 0) {
        config->R_nG = atof(value);
    } else if (strcmp(key, "m_phi") == 0) {
        config->m_phi = atof(value);
    } else if (strcmp(key, "m_chi") == 0) {
        config->m_chi = atof(value);
    } else if (strcmp(key, "phi_w") == 0) {
        config->phi_w = atof(value);
    } else if (strcmp(key, "phi_p") == 0) {
        config->phi_p = atof(value);
    } else if (strcmp(key, "vev") == 0) {
        config->vev = atof(value);
    } else if (strcmp(key, "m_tach") == 0) {
        config->m_tach = atof(value);
    } else if (strcmp(key, "a_e") == 0) {
        config->a_e = atof(value);
    } else if (strcmp(key, "maps") == 0) {
        strncpy(config->maps, value, sizeof(config->maps) - 1);
    } else if (strcmp(key, "nside_map") == 0) {
        config->nside_map = atoi(value);
    } else if (strcmp(key, "npix_map") == 0) {
        config->npix_map = atoi(value);
    } else if (strcmp(key, "fov_map") == 0) {
        config->fov_map = atof(value);
    } else if (strcmp(key, "zmin_map") == 0) {
        config->zmin_map = atof(value);
    } else if (strcmp(key, "zmax_map") == 0) {
        config->zmax_map = atof(value);
    } else if (strcmp(key, "tabfile_map") == 0) {
        strncpy(config->tabfile_map, value, sizeof(config->tabfile_map) - 1);
    } else if (strcmp(key, "tabfile_sfr") == 0) {
        strncpy(config->tabfile_sfr, value, sizeof(config->tabfile_sfr) - 1);
    } else if (strcmp(key, "model_map") == 0) {
        config->model_map = atoi(value);
    } else if (strcmp(key, "scramble_map") == 0) {
        config->scramble_map = atoi(value);
    } else if (strcmp(key, "center_map") == 0) {
        config->center_map = atoi(value);
    } else if (strcmp(key, "chihview_map") == 0) {
        config->chihview_map = atoi(value);
    } else if (strcmp(key, "PSZcut_map") == 0) {
        config->PSZcut_map = atoi(value);
    } else if (strcmp(key, "ellmax") == 0) {
        config->ellmax = atoi(value);
    } else if (strcmp(key, "Omx") == 0) {
        config->Omx = atof(value);
    } else if (strcmp(key, "OmB") == 0) {
        config->OmB = atof(value);
    } else if (strcmp(key, "Omvac") == 0) {
        config->Omvac = atof(value);
    } else if (strcmp(key, "h") == 0) {
        config->h = atof(value);
    } else if (strcmp(key, "ns") == 0) {
        config->ns = atof(value);
    } else if (strcmp(key, "As") == 0) {
        config->As = atof(value);
    } else if (strcmp(key, "sigma8") == 0) {
        config->sigma8 = atof(value);
    } else if (strcmp(key, "tau") == 0) {
        config->tau = atof(value);
    } else if (strcmp(key, "mnu") == 0) {
        config->mnu = atof(value);
    } else if (strcmp(key, "pkfile") == 0) {
        strncpy(config->pkfile, value, sizeof(config->pkfile) - 1);
    } else {
        printf("Unknown key: %s\n", key);
    }
}

void readParameterFile(const char* filename, Config *config) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#' || line[0] == '\n' || isspace(line[0])) continue;
        //parseParameter(config, line);
		printf("1");
    }

    fclose(file);
}

void printExtractedParameters(const Config *config); // Assume defined elsewhere
													 
// Prototype of the function that processes the key and value
void process_parameter_entry(const char* key, const char* value);

// Function to trim leading and trailing whitespace
char* trim_whitespace(char *str) {
    char *end;

    // Trim leading space
    while(isspace((unsigned char)*str)) str++;

    // All spaces?
    if(*str == 0)  
        return str;

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator character
    *(end + 1) = '\0';

    return str;
}

// Function to parse the line and process it
void parse_and_process_line(char *line, const Config *config) {
    char *key_part = strtok(line, "=");
    char *value_part = strtok(NULL, "");

    if (key_part && value_part) {
        char *trimmed_key = trim_whitespace(key_part);
        char *trimmed_value = trim_whitespace(value_part);
        //process_parameter_entry(trimmed_key, trimmed_value);
        set_config_value(&config, key_part, trimmed_value);
    } else {
        printf("Error: Line does not contain a valid 'key = value' format.\n");
    }
}

// Example implementation of the function that processes the key and value
void process_parameter_entry(const char* key, const char* value) {
    printf("Processed Key: '%s', Value: '%s'\n", key, value);
}

void read_param_lines(const char* filename, Config *config) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return;
    }

    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        char *ptr = line;
        int is_empty = 1;  // Flag to check if line is empty
						   
		// Trim everything after the '#' character
        char *hash_pos = strchr(line, '#');
        if (hash_pos) {
        // Check if there is a newline character following the '#'
            if (*(hash_pos + 1) == '\n') {
                *(hash_pos + 1) = '\0'; // Terminate the string after the newline character
            } else {
                *hash_pos = '\0'; // Terminate the string at the position of '#'
            }
        }

        // Check each character in the line
        while (*ptr != '\0') {
            if (!isspace((unsigned char)*ptr)) {
                is_empty = 0;
                break;
            }
            ptr++;
        }

		// Add newline at the end of a line if it doesn't have one
		if (line[strlen(line) - 1] != '\n') strcat(line, "\n");

        // If line is not empty, process it
        if (!is_empty) {
            printf("%s", line);
		    parse_and_process_line(line, config);
			printf("\n");
        }
    }

    fclose(file);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: %s <parameter_file>\n", argv[0]);
        return 1;
    }

    Config config;
    initializeConfig(&config);
    //readParameterFile(argv[1], &config);
	read_param_lines(argv[1], &config);
    //printExtractedParameters(&config);

    return 0;
}
void printExtractedParameters(const Config *config) {
    printf("Machine: %s\n", config->machine);
    printf("Submit Command: %s\n", config->submit_command);
    printf("HPKVD Params: %d\n", config->hpkvd_params);
    printf("Compile HPKVD: %d\n", config->compile_hpkvd);
    printf("Create Filterbank: %d\n", config->create_filterbank);
    printf("Merge Params: %d\n", config->merge_params);
    printf("Compile Merge: %d\n", config->compile_merge);
    printf("Map Params: %d\n", config->map_params);
    printf("Compile Maps: %d\n", config->compile_maps);
    printf("Batch: %d\n", config->batch);
    printf("Submit: %d\n", config->submit);
    printf("Seed: %d\n", config->seed);
    printf("Run Name: %s\n", config->run_name);
    printf("Short Name: %s\n", config->short_name);
    printf("Runtype: %s\n", config->runtype);
    printf("Boxsize: %f\n", config->boxsize);
    printf("Nmesh: %d\n", config->nmesh);
    printf("Nbuff: %d\n", config->nbuff);
    printf("Ntile: %d\n", config->ntile);
    printf("Largerun: %d\n", config->largerun);
    printf("Tlimit: %s\n", config->tlimit);
    printf("Nnodes: %d\n", config->nnodes);
    printf("Tpnode: %d\n", config->tpnode);
    printf("Ntasks: %d\n", config->ntasks);
    printf("Ncpus: %d\n", config->ncpus);
    printf("Nompth: %d\n", config->nompth);
    printf("Ievol: %d\n", config->ievol);
    printf("Num Redshifts: %d\n", config->num_redshifts);
    printf("Maximum Redshift: %f\n", config->maximum_redshift);
    printf("Global Redshift: %f\n", config->global_redshift);
    printf("Ilpt: %d\n", config->ilpt);
    printf("Ioutfield: %d\n", config->ioutfield);
    printf("Ireadfield: %d\n", config->ireadfield);
    printf("Iwant Field Part: %d\n", config->iwant_field_part);
    printf("Field Directory: %s\n", config->fielddir);
    printf("Density File In: %s\n", config->densfilein);
    printf("Density File Out: %s\n", config->densfileout);
    printf("NonGauss: %d\n", config->NonGauss);
    printf("FNL: %f\n", config->fNL);
    printf("A_nG: %f\n", config->A_nG);
    printf("B_nG: %f\n", config->B_nG);
    printf("R_nG: %f\n", config->R_nG);
    printf("M_phi: %f\n", config->m_phi);
    printf("M_chi: %f\n", config->m_chi);
    printf("Phi_w: %f\n", config->phi_w);
    printf("Phi_p: %f\n", config->phi_p);
    printf("Vev: %f\n", config->vev);
    printf("M_tach: %f\n", config->m_tach);
    printf("A_e: %f\n", config->a_e);
    printf("Maps: %s\n", config->maps);
    printf("Nside Map: %d\n", config->nside_map);
    printf("Npix Map: %d\n", config->npix_map);
    printf("FOV Map: %f\n", config->fov_map);
    printf("Zmin Map: %f\n", config->zmin_map);
    printf("Zmax Map: %f\n", config->zmax_map);
    printf("Tabfile Map: %s\n", config->tabfile_map);
    printf("Tabfile SFR: %s\n", config->tabfile_sfr);
    printf("Model Map: %d\n", config->model_map);
    printf("Scramble Map: %d\n", config->scramble_map);
    printf("Center Map: %d\n", config->center_map);
    printf("Chihview Map: %d\n", config->chihview_map);
    printf("PSZcut Map: %d\n", config->PSZcut_map);
    printf("Ellmax: %d\n", config->ellmax);
    printf("Omx: %f\n", config->Omx);
    printf("OmB: %f\n", config->OmB);
    printf("Omvac: %f\n", config->Omvac);
    printf("H: %f\n", config->h);
    printf("Ns: %f\n", config->ns);
    printf("As: %f\n", config->As);
    printf("Sigma8: %f\n", config->sigma8);
    printf("Tau: %f\n", config->tau);
    printf("Mnu: %f\n", config->mnu);
    printf("Pkfile: %s\n", config->pkfile);
    
}
