#ifndef FFT_H
#define FFT_H

#include <complex>
#include <vector>
#include <string>
#include "rarray"

namespace fft
{
  rarray<std::complex<double>, 2> stft_dft(
      rarray<std::complex<double>, 1> &vec,
      size_t window_size, size_t window_step);
  rarray<std::complex<double>, 2> stft_fft(
      rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step);
  rarray<std::complex<double>, 2> stft_ff(
      rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step);

  rarray<std::complex<double>, 1> dft(
      const rarray<std::complex<double>, 1> &signal,
      size_t vec_begin, size_t vec_size);
  rarray<std::complex<double>, 1> fft(
      const rarray<std::complex<double>, 1> &signal,
      size_t vec_begin, size_t vec_size);

  size_t reverse_bits(size_t val, int width);
  bool transform(rarray<std::complex<double>, 1> &vec);
  bool transform(rarray<std::complex<double>, 1> &vec,
                 rarray<std::complex<double>, 2> &fft,
                 rarray<std::complex<double>, 2> &tw);
}

#endif
