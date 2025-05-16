#include "fft.h"
#include "gtest/gtest.h"
#include <complex.h>
#include <fftw3.h>
#include <iostream>
#include <tuple>
#include <typeinfo>

namespace {
  TEST(FFTTest, Small) {
    int num_samples = 16;
    int window_size = 8;

    // Calculate actual results.
    rarray<std::complex<double>, 2> stft;
    rarray<std::complex<double>, 1> signal(num_samples);
    for (size_t i = 0; i < num_samples; i++)
      signal[i] = std::complex<double>(double(i), 0.0);
    stft = fft::stft_fft(signal, window_size, 1);

    // Check dimension of results.
    EXPECT_EQ(stft.shape()[0], num_samples - window_size + 1) << "Incorrect row size.";
    EXPECT_EQ(stft.shape()[1], window_size) << "Incorrect column size.";

    // Calculate expected results using FFTW.
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    p = fftw_plan_dft_1d(window_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    for (int j = 0; j < num_samples - window_size + 1; j++)
    {
      for (int i = 0; i < window_size; i++)
      {
        in[i][0] = double(i + j);
        in[i][1] = 0.0;
      }
      fftw_execute(p);
      for (int i = 0; i < window_size; i++)
      {
        EXPECT_NEAR(out[i][0], stft[j][i].real(), 1e-7) << "Incorrect real part at (" << j << "," << i << ").";
        EXPECT_NEAR(out[i][1], stft[j][i].imag(), 1e-7) << "Incorrect imaginary part at (" << j << "," << i << ").";
      }
    }
    fftw_destroy_plan(p);
    fftw_free(in); 
    fftw_free(out);
  }
}
