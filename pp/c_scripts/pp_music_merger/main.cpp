#include <iostream>

extern "C" {
    void initialize_array();
    void get_array_pointer(void** ptr);
}

int main() {
    initialize_array();
    void* ptr;
    get_array_pointer(&ptr);

    double* array = static_cast<double*>(ptr);

    int n = 3;

    // Print array as seen by C++ (row-major order)
    std::cout << "Array as seen by C++ (row-major order):" << std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                std::cout << "array[" << i << "][" << j << "][" << k << "] = " 
                          << array[i * n * n + j * n + k] << std::endl;
            }
        }
    }

    // Print array respecting Fortran's column-major order
    std::cout << "Array respecting Fortran's column-major order:" << std::endl;
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                std::cout << "array[" << i << "][" << j << "][" << k << "] = " 
                          << array[i + j * n + k * n * n] << std::endl;
            }
        }
    }

    return 0;
}

