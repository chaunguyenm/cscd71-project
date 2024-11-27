#include "fft.h"
#include <fstream>
#include <iostream>

#define PI 3.14159265

rarray<std::complex<double>, 2> fft::stft_dft(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t spectrogram_size = ceil((vec.size() - window_size) / float(window_step));
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);

#pragma omp parallel for default(none) shared(vec, window_size, window_step, spectrogram, spectrogram_size)
  for (size_t i = 0; i < spectrogram_size; i++)
  {
    rarray<std::complex<double>, 1> dft = fft::dft(vec, i * window_step, window_size);

#pragma omp parallel for default(none) shared(window_size, dft, i, spectrogram)
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = dft[j];
  }
  return spectrogram;
}

rarray<std::complex<double>, 1> fft::dft(
    const rarray<std::complex<double>, 1> &signal,
    size_t vec_begin, size_t vec_size)
{
  if (vec_size == 0)
    vec_size = signal.size();
  size_t vec_end = vec_begin + vec_size;

  rarray<std::complex<double>, 1> output(vec_size);
  for (size_t k = 0; k < vec_size; k++)
  {
    std::complex<double> sum = 0.;
    for (size_t t = 0; t < vec_size; t++)
    {
      float angle = 2 * PI * t * k / vec_size;
      sum += signal[t + vec_begin] * std::exp(std::complex<double>(0, -angle));
    }
    output[k] = sum;
  }
  return output;
}

rarray<std::complex<double>, 2> fft::stft_fft(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t spectrogram_size = ceil((vec.size() - window_size) / float(window_step));
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);

#pragma omp parallel for default(none) shared(vec, window_size, window_step, spectrogram, spectrogram_size)
  for (size_t i = 0; i < spectrogram_size; i++)
  {
    rarray<std::complex<double>, 1> fft = fft::fft(vec, i * window_step, window_size);

#pragma omp parallel for default(none) shared(window_size, fft, i, spectrogram)
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = fft[j];
  }
  return spectrogram;
}

rarray<std::complex<double>, 2> fft::stft(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  /*
   * STEP 1: Create buffers M and M'
   */
  size_t num_stages = log2(window_size) + 1;
  size_t spectrogram_size = ceil((vec.size() - window_size) / float(window_step));
  size_t buffer_rows = window_size / 2;
  size_t buffer_cols = window_size / 2 - 1;

  rarray<std::vector<std::complex<double>>, 2> M(buffer_rows, buffer_cols),
      M_prime(buffer_rows, buffer_cols);
  rarray<std::complex<double>, 2> fft0(window_size, num_stages),
      spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 1> exptable(window_size / 2);

  /*
   * STEP 2: Initialize M
   */
  // Perform FFT on 1st N signals
  rarray<std::complex<double>, 1> subvec(window_size);
  for (size_t i = 0; i < window_size; i++)
    subvec[i] = vec[i];
  transform(subvec, fft0, exptable);
  std::cout << fft0 << "\n";
  // Set M to initial FFT
  for (size_t i = 0; i < window_size / 2; i++)
    M[i][0].push_back(vec[i + window_size / 2]);
  for (size_t s = 1; s < num_stages - 1; s++)
    for (size_t k = 0; k < window_size / pow(2, s + 1); k++)
      for (size_t m = 0; m < pow(2, s); m++)
        M[k][s].push_back(fft0[(2 * k + 1) * pow(2, s) + m][s]);

  for (size_t r = 0; r < buffer_rows; r++)
  {
    for (size_t c = 0; c < buffer_cols; c++)
    {
      std::cout << r << "," << c << "=[";
      for (const auto &elem : M[r][c])
        std::cout << elem << ",";
      std::cout << "] ";
    }
    std::cout << "\n";
  }

  /*
   * STEP 3: Compute batches of N/2 FFTs
   */
  for (size_t n = 1; n < vec.size() - window_size + 1; n++)
    for (size_t s = 1; s < num_stages; s++)
    {
      // Initialize temporary input buffer
      rarray<std::vector<std::complex<double>>, 1> B0(window_size / 2);
      rarray<std::vector<std::complex<double>>, 1> B1(window_size / 2);
      rarray<size_t, 1> B2(window_size / 2);

      for (size_t i = 0; i < window_size / 2; i++)
      {
        B0[i].push_back(vec[i + n + window_size - 1]);
        B1[i] = M[i][s - 1];
        B2[i] = s;

        // std::cout << "B0[" << i << "]=" << vec[i + n + window_size - 1] << "\n";
        // std::cout << "B1[" << i << "]=";
        // for (const auto &element : M[i][s - 1])
        //   std::cout << element << ", ";
        // std::cout << "\n";
        // std::cout << "B2[" << i << "]=" << s << "\n";
      }

      // std::cout << B0 << "\n";
      // std::cout << B1 << "\n";
      // std::cout << B2 << "\n";

      // Parallel
      for (size_t k = 0; k < pow(2, s) - 2; k++)
      {
        std::complex<double> xr = B0[0][k] * std::polar(1.0, -2. * PI * 0 / n);
        std::complex<double> x0 = B1[0][k] + xr;
        std::complex<double> x1 = B1[0][k] - xr;

        // std::cout << x0 << " " << x1 << "\n";
      }
    }

  return spectrogram;
}

rarray<std::complex<double>, 1> fft::fft(
    const rarray<std::complex<double>, 1> &signal,
    size_t vec_begin, size_t vec_size)
{
  rarray<std::complex<double>, 1> subsignal(vec_size);
  for (size_t i = 0; i < vec_size; i++)
    subsignal[i] = signal[vec_begin + i];

  fft::transform(subsignal);
  return subsignal;
}

bool fft::transform(rarray<std::complex<double>, 1> &vec)
{
  // Length variables
  size_t n = vec.size();
  int levels = 0; // Compute levels = floor(log2(n))
  for (size_t temp = n; temp > 1U; temp >>= 1)
    levels++;
  if ((size_t)1U << levels != n)
    return false; // n is not a power of 2

  // Trigonometric tables
  rarray<std::complex<double>, 1> exptable(n / 2);
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

bool fft::transform(rarray<std::complex<double>, 1> &vec,
                    rarray<std::complex<double>, 2> &fft,
                    rarray<std::complex<double>, 1> &exptable)
{
  // Length variables
  size_t n = vec.size();
  int levels = 0; // Compute levels = floor(log2(n))
  for (size_t temp = n; temp > 1U; temp >>= 1)
    levels++;
  if ((size_t)1U << levels != n)
    return false; // n is not a power of 2

  // Trigonometric tables
  // rarray<std::complex<double>, 1> exptable(n / 2);
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
        // std::cout << j << "," << log2(size) << "=" << exptable[k] << " ";
        // std::cout << l << "," << log2(size) << "=" << exptable[k] << " ";
      }
      // std::cout << "\n";
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
