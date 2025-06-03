#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846

void fft(float complex *x, int n, int inverse) {
    if (n <= 1) return;

    float complex *even = malloc(n/2 * sizeof *even);
    float complex *odd  = malloc(n/2 * sizeof *odd);
    for (int i = 0; i < n/2; i++) {
        even[i] = x[2*i];
        odd[i]  = x[2*i + 1];
    }

    fft(even, n/2, inverse);
    fft(odd,  n/2, inverse);

    // Heart of the FFT algorithm
    // Combine results
    for (int k = 0; k < n/2; k++) {
        float sign = inverse ? +1.0f : -1.0f;
        float complex t = cexpf(sign * 2.0f * PI * I * k / n) * odd[k]; // cexpf(...) computes e^(iθ) = cos(θ) + isin(θ). This is the "twiddle factor" multiplied by the odd-indexed value
        x[k] = even[k] + t;
        x[k + n/2] = even[k] - t;
    }

    free(even);
    free(odd);

    // if (!inverse) {
    //     // normalize forward FFT so that a subsequent inverse
    //     // (with the same function) will restore original amplitude
    //     for (int i = 0; i < n; i++)
    //         x[i] /= 2;  // or divide by n at the top level instead
    // }
}

int main() {
    int n = 16;
    float complex *signal = malloc(n * sizeof *signal);

    printf("Input signal:\n");
    for (int i = 0; i < n; i++)
        signal[i] = i + 0.0f * I,
        printf("%2.0f + %2.0fi\n", crealf(signal[i]), cimagf(signal[i]));

    fft(signal, n, 0);

    printf("\nFFT results:\n");
    for (int i = 0; i < n; i++)
        printf("% .5f + % .5fi\n", crealf(signal[i]), cimagf(signal[i]));


    // fft(signal, n, 1);  // Perform inverse FFT (second parameter is 1)

    // printf("\nInverse FFT results (should match input):\n");
    // for (int i = 0; i < n; i++)
    //     printf("% .5f + % .5fi\n", crealf(signal[i]), cimagf(signal[i]));

    free(signal);
     // Add this line to pause the program
    // printf("\nPress Enter to continue...");
    // getchar();
    return 0;
}
