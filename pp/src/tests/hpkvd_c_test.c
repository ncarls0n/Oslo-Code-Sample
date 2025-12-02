#include <stdio.h>

extern void c_run_hpkvd(float *input_field, int *maketable_flag, int *seedFFT);

int main() {
    int maketable_flag = 1; // Example value
    int seedFFT = 13579;    // Example value
    float input_field[10][10][10] = {{{0}}}; // Example 3D input field

    // Call Fortran subroutine, you can pass NULL if no input field is required
    c_run_hpkvd(&input_field[0][0][0], &maketable_flag, &seedFFT);

    return 0;
}
