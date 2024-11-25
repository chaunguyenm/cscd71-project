#ifndef FFT_H
#define FFT_H

#include <complex>
#include <vector>
#include <string>
#include "rarray"

namespace fft
{
  std::vector<std::vector<std::complex<double>>> stft_dft(
      std::vector<std::complex<double>> &vec,
      size_t window_size, size_t window_step);
  std::vector<std::vector<std::complex<double>>> stft_fft(
      std::vector<std::complex<double>> &vec, size_t window_size, size_t window_step);
  std::vector<std::vector<std::complex<double>>> stft2_fft(
      std::vector<std::complex<double>> &vec, size_t window_size, size_t window_step);
  std::vector<std::complex<double>> dft(
      const std::vector<std::complex<double>> &signal,
      size_t vec_begin, size_t vec_size);
  size_t reverse_bits(size_t val, int width);
  bool transformRadix2(std::vector<std::complex<double>> &vec);
  bool transformRadix2(std::vector<std::complex<double>> &vec,
                       rarray<std::complex<double>, 2> &fft);
  std::vector<std::complex<double>> fft(
      const std::vector<std::complex<double>> &signal,
      size_t vec_begin, size_t vec_size);
}

#endif
