// File used to generate reference outputs used in tests
// Download HIIR from https://ldesoras.fr/src/hiir-1.50.zip
//
// Run with:
// clang++ -std=c++14 -I. -o test_hiir.exe test_hiir.cpp hiir/*.cpp -O2 && test_hiir

#include <iostream>
#include <iomanip>
#include "hiir/PolyphaseIir2Designer.h"
#include "hiir/Downsampler2xFpu.h"
#include "hiir/Upsampler2xFpu.h"

void print_arr(float arr [], int len) {
    for (int i = 0; i < len; ++i) {
        std::cout << "    " << arr[i] << "," << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    const int precision = std::numeric_limits<float>::max_digits10;
    std::cout << std::setprecision(precision);

    // Coefficient calculation
    const int n_coef = 8;
    double coefs[n_coef];

    std::cout << "=== Computing coefficients ===" << std::endl;
    std::cout << "Number of coefficients: " << n_coef << std::endl;
    std::cout << "Transition bandwidth: 0.01" << std::endl << std::endl;

    hiir::PolyphaseIir2Designer::coefs_spec_spec_order_tbw(
        coefs, n_coef, 0.01
    );

    // Convert to float for printing
    float coefs_f[n_coef];
    for (int i = 0; i < n_coef; ++i) {
        coefs_f[i] = static_cast<float>(coefs[i]);
    }

    std::cout << "Coefficients:" << std::endl;
    print_arr(coefs_f, n_coef);

    // Test downsampler with simple input
    std::cout << "=== Downsampler ===" << std::endl;
    hiir::Downsampler2xFpu<n_coef> downsampler;
    downsampler.set_coefs(coefs);

    // Create test signal
    float src[16] = {1.0f, 0.5f, 0.0f, -0.5f, -1.0f, -0.5f, 0.0f, 0.5f, 1.0f, 0.5f, 0.0f, -0.5f, -1.0f, -0.5f, 0.0f, 0.5f};
    float dst[8];

    std::cout << "Input:" << std::endl;
    print_arr(src, 16);

    downsampler.process_block(dst, src, 8);

    std::cout << "Output (downsampled):" << std::endl;
    print_arr(dst, 8);

    // Test upsampler
    std::cout << "=== Upsampler ===" << std::endl;
    hiir::Upsampler2xFpu<n_coef> upsampler;
    upsampler.set_coefs(coefs);

    float src_up[8] = {1.0f, 0.5f, 0.0f, -0.5f, -1.0f, -0.5f, 0.0f, 0.5f};
    float dst_up[16];

    std::cout << "Input:" << std::endl;
    print_arr(src_up, 8);

    upsampler.process_block(dst_up, src_up, 8);

    std::cout << "Output (upsampled):" << std::endl;
    print_arr(dst_up, 16);

    return 0;
}
