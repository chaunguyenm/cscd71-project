#include "fft.h"
#include <fstream>
#include <iostream>

#define PI 3.14159265

std::vector<std::vector<std::complex<double>>> fft::stft_dft(
    std::vector<std::complex<double>> &vec, size_t window_size, size_t window_step)
{
  std::vector<std::vector<std::complex<double>>> spectrogram = std::vector<std::vector<std::complex<double>>>();
  for (size_t begin = 0; begin < vec.size() - window_size; begin += window_step)
  {
    spectrogram.push_back(fft::dft(vec, begin, window_size));
  }
  return spectrogram;
}

std::vector<std::vector<std::complex<double>>> fft::stft_fft(
    std::vector<std::complex<double>> &vec, size_t window_size, size_t window_step)
{
  std::vector<std::vector<std::complex<double>>> spectrogram = std::vector<std::vector<std::complex<double>>>();
  for (size_t begin = 0; begin < vec.size() - window_size; begin += window_step)
  {
    spectrogram.push_back(fft::fft(vec, begin, window_size));
  }
  return spectrogram;
}

std::vector<std::vector<std::complex<double>>> fft::stft2_fft(
    std::vector<std::complex<double>> &vec, size_t window_size, size_t window_step)
{
  size_t num_stages = log2(window_size) + 1;
  rarray<std::complex<double>, 2> M(window_size / 2, window_size / 2 - 1),
      M_prime(window_size / 2, window_size / 2 - 1), fft0(window_size, num_stages);
  
  // Initialize M
  transformRadix2(vec, fft0);
  for (size_t s = 0; s < num_stages; s++)
    for (size_t k = 0; k < window_size / pow(2, s + 1); k++)
      M[s][k] = fft0[s][2 * k + 1];
}

std::vector<std::complex<double>> fft::dft(
    const std::vector<std::complex<double>> &signal,
    size_t vec_begin, size_t vec_size)
{
  std::vector<std::complex<double>> output;
  for (size_t k = 0; k < vec_size / 2.; k++)
  {
    std::complex<double> sum = 0.;
    for (size_t t = 0; t < vec_size; t++)
    {
      float angle = 2 * PI * t * k / vec_size;
      sum += signal[t + vec_begin] * std::exp(std::complex<double>(0, -angle));
    }
    output.push_back(sum);
  }
  return output;
}

std::vector<std::complex<double>> fft::fft(const std::vector<std::complex<double>> &signal, size_t vec_begin, size_t vec_size)
{
  std::vector<std::complex<double>>::const_iterator first = signal.begin() + vec_begin;
  std::vector<std::complex<double>>::const_iterator last = signal.begin() + vec_begin + vec_size;
  std::vector<std::complex<double>> subsignal(first, last);

  rarray<std::complex<double>, 2> fft(8, 4);
fft:
  transformRadix2(subsignal, fft);
  return subsignal;
}

bool fft::transformRadix2(std::vector<std::complex<double>> &vec)
{
  // Length variables
  size_t n = vec.size();
  int levels = 0; // Compute levels = floor(log2(n))
  for (size_t temp = n; temp > 1U; temp >>= 1)
    levels++;
  if ((size_t)1U << levels != n)
    return false; // n is not a power of 2

  // Trigonometric tables
  std::vector<std::complex<double>> exptable(n / 2);
  for (size_t i = 0; i < n / 2; i++)
    exptable[i] = std::polar(1.0, -2. * PI * i / n);

  // Bit-reversed addressing permutation
  for (size_t i = 0; i < n; i++)
  {
    size_t j = fft::reverse_bits(i, levels);
    if (j > i)
    {
      std::complex<double> temp = vec[i];
      vec[i] = vec[j];
      vec[j] = temp;
    }
  }

  // Cooley-Tukey decimation-in-time radix-2 FFT
  for (size_t size = 2; size <= n; size *= 2)
  {
    size_t halfsize = size / 2;
    size_t tablestep = n / size;
    for (size_t i = 0; i < n; i += size)
    {
      for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep)
      {
        size_t l = j + halfsize;
        std::complex<double> temp = vec[l] * exptable[k];
        vec[l] = vec[j] - temp;
        vec[j] += temp;
      }
    }
    if (size == n) // Prevent overflow in 'size *= 2'
      break;
  }

  return true;
}

bool fft::transformRadix2(std::vector<std::complex<double>> &vec,
                          rarray<std::complex<double>, 2> &fft)
{
  // Length variables
  size_t n = vec.size();
  int levels = 0; // Compute levels = floor(log2(n))
  for (size_t temp = n; temp > 1U; temp >>= 1)
    levels++;
  if ((size_t)1U << levels != n)
    return false; // n is not a power of 2

  // Trigonometric tables
  std::vector<std::complex<double>> exptable(n / 2);
  for (size_t i = 0; i < n / 2; i++)
    exptable[i] = std::polar(1.0, -2. * PI * i / n);

  // Bit-reversed addressing permutation
  for (size_t i = 0; i < n; i++)
  {
    size_t j = fft::reverse_bits(i, levels);
    if (j > i)
    {
      std::complex<double> temp = vec[i];
      vec[i] = vec[j];
      vec[j] = temp;
    }
  }
  for (size_t i = 0; i < n; i++)
    fft[i][0] = vec[i];

  // Cooley-Tukey decimation-in-time radix-2 FFT
  for (size_t size = 2; size <= n; size *= 2)
  {
    size_t halfsize = size / 2;
    size_t tablestep = n / size;
    for (size_t i = 0; i < n; i += size)
    {
      for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep)
      {
        size_t l = j + halfsize;
        std::complex<double> temp = vec[l] * exptable[k];
        vec[l] = vec[j] - temp;
        vec[j] += temp;
        fft[j][log2(size)] = vec[j];
        fft[l][log2(size)] = vec[l];
      }
    }
    if (size == n) // Prevent overflow in 'size *= 2'
      break;
  }

  return true;
}

size_t fft::reverse_bits(size_t val, int width)
{
  size_t result = 0;
  for (int i = 0; i < width; i++, val >>= 1)
    result = (result << 1) | (val & 1U);
  return result;
}
