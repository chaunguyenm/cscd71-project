#include "fft.h"
#include <fstream>
#include <iostream>
#include <mpi.h>

#define PI 3.14159265

rarray<std::complex<double>, 2> fft::stft_dft(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 1> dft;

#pragma omp parallel for default(none) shared(vec, window_size, window_step, spectrogram, spectrogram_size)
  for (size_t i = 0; i < spectrogram_size; i++)
  {
    dft = fft::dft(vec, i * window_step, window_size);

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
  rarray<std::complex<double>, 1> output(vec_size);
#pragma omp parallel for default(none) shared(signal, vec_begin, vec_size, output)
  for (size_t k = 0; k < vec_size; k++)
  {
    std::complex<double> sum = 0.;
    for (size_t t = 0; t < vec_size; t++)
    {
      double angle = 2 * PI * t * k / vec_size;
      sum += signal[t + vec_begin] * std::exp(std::complex<double>(0, -angle));
    }
    output[k] = sum;
  }
  return output;
}

rarray<std::complex<double>, 2> fft::stft_fft(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 1> fft;

#pragma omp parallel for default(none) shared(vec, window_size, window_step, spectrogram, spectrogram_size)
  for (size_t i = 0; i < spectrogram_size; i++)
  {
    fft = fft::fft(vec, i * window_step, window_size);

#pragma omp parallel for default(none) shared(window_size, fft, i, spectrogram)
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = fft[j];
  }
  return spectrogram;
}

rarray<std::complex<double>, 2> fft::stft_ff(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t num_stages = log2(window_size) + 1;
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> fft(window_size, num_stages),
      spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 2> tw(window_size, num_stages - 1);
  rarray<std::complex<double>, 2> buffer(window_size / 2, num_stages);

  // Perform FFT on 1st N signals
  rarray<std::complex<double>, 1> subvec(window_size);
#pragma omp parallel for default(none) shared(window_size, subvec, vec)
  for (size_t i = 0; i < window_size; i++)
    subvec[i] = vec[i];
  transform(subvec, fft, tw);

#pragma omp parallel for default(none) shared(window_size, spectrogram, fft, num_stages)
  for (size_t i = 0; i < window_size; i++)
    spectrogram[0][i] = fft[i][num_stages - 1];

    // Prepare buffer
#pragma omp paralel for default(none) shared(num_stages, window_size, buffer, fft)
  for (size_t stage = 0; stage < num_stages; stage++)
  {
    size_t count = 1 << stage;
    for (size_t r = 0; (2 * r + 1) * count < window_size; r++)
    {
      size_t buffer_idx = reverse_bits(r, log2(window_size / 2));
      for (size_t c = 0; c < count; c++)
        buffer[buffer_idx + c][stage] = fft[(2 * r + 1) * count + c][stage];
    }
  }

  for (size_t i = 1; i + window_size - 1 < vec.size(); i++)
  {
    for (size_t s = 1; s < num_stages; s++)
    {
      size_t count = 1 << (s - 1);
      size_t m = ((i - 1) % (window_size / (1 << s))) * count;
      for (size_t k = 0; k < count; k++)
      {
        std::complex<double> xr;
        if (s == 1)
        {
          xr = vec[i + window_size - 1];
          fft[window_size - (1 << (s - 1)) + k][s - 1] = xr;
        }
        else
          xr = fft[window_size - (1 << (s - 1)) + k][s - 1];
        fft[window_size - (1 << s) + k][s] =
            buffer[m + k][s - 1] + tw[window_size - (1 << s) + k][s - 1] * xr;
        fft[window_size - (1 << s) + k + count][s] =
            buffer[m + k][s - 1] - tw[window_size - (1 << s) + k + count][s - 1] * xr;
      }
      for (size_t k = 0; k < count; k++)
        buffer[m + k][s - 1] = fft[window_size - (1 << (s - 1)) + k][s - 1];
    }
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = fft[j][num_stages - 1];
  }

  return spectrogram;
}

rarray<std::complex<double>, 2> fft::stft_qpff(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t num_stages = log2(window_size) + 1;
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> fft(window_size, num_stages),
      spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 2> tw(window_size, num_stages - 1);
  rarray<std::complex<double>, 2> buffer(window_size / 2, num_stages);

  rarray<std::complex<double>, 3>
      M(window_size / 2, num_stages - 1, window_size / 2),
      M_prime(window_size / 2, num_stages - 1, window_size / 2);

  rarray<std::complex<double>, 3> B(window_size / 2, 2, window_size / 2);
  rarray<std::complex<double>, 2> f(window_size / 2, window_size);

  // Perform FFT on 1st N signals
  rarray<std::complex<double>, 1> subvec(window_size);
#pragma omp parallel for default(none) shared(window_size, subvec, vec)
  for (size_t i = 0; i < window_size; i++)
    subvec[i] = vec[i];
  transform(subvec, fft, tw);

#pragma omp parallel for default(none) shared(window_size, spectrogram, fft, num_stages)
  for (size_t i = 0; i < window_size; i++)
    spectrogram[0][i] = fft[i][num_stages - 1];

    // Initialize M
#pragma omp parallel for default(none) shared(num_stages, window_size, M, fft)
  for (size_t s = 0; s < num_stages - 1; s++)
  {
    size_t num_rows = window_size / (1 << (s + 1));
    for (size_t k = 0; k < num_rows; k++)
    {
      size_t row_idx = reverse_bits(k, log2(num_rows));
      for (size_t m = 0; m < (1 << s); m++)
      {
        M[row_idx][s][m] = fft[(2 * k + 1) * (1 << s) + m][s];
      }
    }
  }

  for (size_t n = 1; window_size / 2 - 1 + n <= vec.size() - window_size; n += window_size / 2)
  {
#pragma omp parallel for default(none) shared(window_size, spectrogram, M, M_prime, vec, n, num_stages, M_prime)
    for (size_t s = 0; s < num_stages - 1; s++)
    {
      for (size_t i = 0; i < window_size / 2; i++)
      {
        for (size_t idx = 0; idx < window_size / 2; idx++)
          B[i][1][idx] = M[i][s][idx];

        if (s == 0)
        {
          B[i][0][0] = vec[n + window_size - 1 + i];
          M_prime[i][s][0] = vec[n + window_size - 1 + i];
        }
        else
        {
          size_t start_idx = i + window_size / (1 << (s + 1));
          if (start_idx >= window_size / 2)
          {
            start_idx = start_idx % (window_size / 2);
            for (size_t idx = 0; idx < window_size / 2; idx++)
              B[i][0][idx] = M_prime[start_idx][s][idx];
          }
          else
            for (size_t idx = 0; idx < window_size / 2; idx++)
              B[i][0][idx] = M[start_idx][s][idx];
        }
      }

      // Parallel B
      f = fft::stft_qpff_batch(B, tw, window_size, s);

#pragma omp parallel for default(none) shared(spectrogram, M_prime, M, f, num_stages, window_size)
      for (size_t i = 0; i < window_size / 2; i++)
      {
        size_t col_idx = s + 1;

        if (col_idx < num_stages - 1)
        {
          size_t start_idx = i + (window_size / (1 << (col_idx + 1)));

          if (start_idx >= window_size / 2)
          {
            start_idx = start_idx % (window_size / 2);
            for (size_t idx = 0; idx < window_size / 2; idx++)
              M_prime[start_idx][col_idx][idx] = f[i][idx];
          }
          else
          {
            for (size_t idx = 0; idx < window_size / 2; idx++)
              M[start_idx][col_idx][idx] = f[i][idx];
          }
        }
        else
          for (size_t idx = 0; idx < window_size; idx++)
            spectrogram[n + i][idx] = f[i][idx];
      }
    }

    M = M_prime.copy();
  }
  return spectrogram;
}

rarray<std::complex<double>, 2> fft::stft_qpff_batch(
    rarray<std::complex<double>, 3> B, rarray<std::complex<double>, 2> tw,
    size_t window_size, size_t s)
{
  size_t count = (1 << s);
  rarray<std::complex<double>, 2> f(window_size / 2, window_size);

  for (size_t i = 0; i < window_size / 2; i++)
  {
    for (size_t k = 0; k < count; k++)
    {
      std::complex<double> xr, x0, x1;
      xr = B[i][0][k] * tw[window_size - (1 << (s + 1)) + k][s];
      x0 = B[i][1][k] + xr;
      x1 = B[i][1][k] - xr;

      f[i][k] = x0;
      f[i][k + count] = x1;
    }
  }

  return f;
}

rarray<std::complex<double>, 1> fft::fft(
    const rarray<std::complex<double>, 1> &signal,
    size_t vec_begin, size_t vec_size)
{
  rarray<std::complex<double>, 1> subsignal(vec_size);
#pragma omp parallel for default(none) shared(signal, vec_begin, vec_size, subsignal)
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
#pragma omp parallel for default(none) shared(exptable, n)
  for (size_t i = 0; i < n / 2; i++)
    exptable[i] = std::polar(1.0, -2. * PI * i / n);

    // Bit-reversed addressing permutation
#pragma omp parallel for default(none) shared(levels, vec, n)
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
                    rarray<std::complex<double>, 2> &tw)
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
#pragma omp parallel for default(none) shared(exptable, n)
  for (size_t i = 0; i < n / 2; i++)
    exptable[i] = std::polar(1.0, -2. * PI * i / n);

    // Bit-reversed addressing permutation
#pragma omp parallel for default(none) shared(vec, levels, n)
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
#pragma omp parallel for default(none) shared(fft, vec, n)
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
        tw[j][log2(size) - 1] = exptable[k];
        tw[l][log2(size) - 1] = exptable[k];
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

rarray<std::complex<double>, 2> fft_mpi::stft_dft(
    rarray<std::complex<double>, 1> &vec,
    size_t window_size, size_t window_step)
{
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 1> dft;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t chunk_size = spectrogram_size / size;
  size_t start_idx = rank * chunk_size;
  size_t end_idx = (rank == size - 1) ? spectrogram_size : (rank + 1) * chunk_size;

  for (size_t i = start_idx; i < end_idx; i++)
  {
    dft = fft::dft(vec, i * window_step, window_size);
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = dft[j];
  }

  if (rank == 0)
    for (int p = 1; p < size; p++)
    {
      size_t received_chunk_size = (p == size - 1) ? spectrogram_size - p * chunk_size : chunk_size;
      MPI_Recv(&spectrogram[p * chunk_size][0],
               received_chunk_size * window_size, MPI_DOUBLE_COMPLEX, p, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  else
    MPI_Send(&spectrogram[start_idx][0],
             (end_idx - start_idx) * window_size, MPI_DOUBLE_COMPLEX, 0, 0,
             MPI_COMM_WORLD);

  return spectrogram;
}

rarray<std::complex<double>, 2> fft_mpi::stft_fft(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 1> fft;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t chunk_size = spectrogram_size / size;
  size_t start_idx = rank * chunk_size;
  size_t end_idx = (rank == size - 1) ? spectrogram_size : (rank + 1) * chunk_size;

  for (size_t i = start_idx; i < end_idx; i++)
  {
    fft = fft::fft(vec, i * window_step, window_size);
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = fft[j];
  }

  if (rank == 0)
    for (int p = 1; p < size; p++)
    {
      size_t received_chunk_size = (p == size - 1) ? spectrogram_size - p * chunk_size : chunk_size;
      MPI_Recv(&spectrogram[p * chunk_size][0],
               received_chunk_size * window_size, MPI_DOUBLE_COMPLEX, p, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  else
    MPI_Send(&spectrogram[start_idx][0],
             (end_idx - start_idx) * window_size, MPI_DOUBLE_COMPLEX, 0, 0,
             MPI_COMM_WORLD);

  return spectrogram;
}

rarray<std::complex<double>, 2> fft_mpi::stft_qpff(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t num_stages = log2(window_size) + 1;
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> fft(window_size, num_stages),
      spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 2> tw(window_size, num_stages - 1);
  rarray<std::complex<double>, 2> buffer(window_size / 2, num_stages);

  rarray<std::complex<double>, 3>
      M(window_size / 2, num_stages - 1, window_size / 2),
      M_prime(window_size / 2, num_stages - 1, window_size / 2);

  rarray<std::complex<double>, 3> B(window_size / 2, 2, window_size / 2);
  rarray<std::complex<double>, 2> f(window_size / 2, window_size);

  // Perform FFT on 1st N signals
  rarray<std::complex<double>, 1> subvec(window_size);
  for (size_t i = 0; i < window_size; i++)
    subvec[i] = vec[i];
  fft::transform(subvec, fft, tw);

  for (size_t i = 0; i < window_size; i++)
    spectrogram[0][i] = fft[i][num_stages - 1];

  // Initialize M
  for (size_t s = 0; s < num_stages - 1; s++)
  {
    size_t num_rows = window_size / (1 << (s + 1));
    for (size_t k = 0; k < num_rows; k++)
    {
      size_t row_idx = fft::reverse_bits(k, log2(num_rows));
      for (size_t m = 0; m < (1 << s); m++)
      {
        M[row_idx][s][m] = fft[(2 * k + 1) * (1 << s) + m][s];
      }
    }
  }

  for (size_t n = 1; window_size / 2 - 1 + n <= vec.size() - window_size; n += window_size / 2)
  {
    for (size_t s = 0; s < num_stages - 1; s++)
    {
      for (size_t i = 0; i < window_size / 2; i++)
      {
        for (size_t idx = 0; idx < window_size / 2; idx++)
          B[i][1][idx] = M[i][s][idx];

        if (s == 0)
        {
          B[i][0][0] = vec[n + window_size - 1 + i];
          M_prime[i][s][0] = vec[n + window_size - 1 + i];
        }
        else
        {
          size_t start_idx = i + window_size / (1 << (s + 1));
          if (start_idx >= window_size / 2)
          {
            start_idx = start_idx % (window_size / 2);
            for (size_t idx = 0; idx < window_size / 2; idx++)
              B[i][0][idx] = M_prime[start_idx][s][idx];
          }
          else
            for (size_t idx = 0; idx < window_size / 2; idx++)
              B[i][0][idx] = M[start_idx][s][idx];
        }
      }

      // Parallel B
      fft_mpi::stft_qpff_batch(B, tw, window_size, s, f);

      // Collect results
      for (size_t i = 0; i < window_size / 2; i++)
      {
        size_t col_idx = s + 1;

        if (col_idx < num_stages - 1)
        {
          size_t start_idx = i + (window_size / (1 << (col_idx + 1)));

          if (start_idx >= window_size / 2)
          {
            start_idx = start_idx % (window_size / 2);
            for (size_t idx = 0; idx < window_size / 2; idx++)
              M_prime[start_idx][col_idx][idx] = f[i][idx];
          }
          else
          {
            for (size_t idx = 0; idx < window_size / 2; idx++)
              M[start_idx][col_idx][idx] = f[i][idx];
          }
        }
        else
        {
          for (size_t idx = 0; idx < window_size; idx++)
            spectrogram[n + i][idx] = f[i][idx];
        }
      }
    }

    M = M_prime.copy();
  }

  return spectrogram;
}

void fft_mpi::stft_qpff_batch(
    rarray<std::complex<double>, 3> B, rarray<std::complex<double>, 2> tw,
    size_t window_size, size_t s, rarray<std::complex<double>, 2> f)
{
  int rank, size, chunksize, remaining;
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Bcast(B.data(), B.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

  chunksize = (window_size / 2) / size;
  remaining = (window_size / 2) % size;

  size_t count = (1 << s);

  int local_chunksize = rank < remaining ? chunksize + 1 : chunksize;
  int offset = rank * chunksize + (remaining > rank ? rank : remaining);
  std::complex<double> local_f[local_chunksize * window_size];

  // Compute local f for this process
  for (size_t i = 0; i < local_chunksize; i++)
  {
    for (size_t k = 0; k < count; k++)
    {
      std::complex<double> xr, x0, x1;
      xr = B[i + offset][0][k] * tw[window_size - (1 << (s + 1)) + k][s];
      x0 = B[i + offset][1][k] + xr;
      x1 = B[i + offset][1][k] - xr;

      local_f[i * window_size + k] = x0;
      local_f[i * window_size + (k + count)] = x1;
    }
  }

  if (rank == 0)
  {
    // Copy from local_f of rank 0 to global f
    for (size_t i = 0; i < local_chunksize; i++)
    {
      for (size_t k = 0; k < count; k++)
      {
        f[i][k] = local_f[i * window_size + k];
        f[i][k + count] = local_f[i * window_size + (k + count)];
      }
    }

    // Receive from remaining processes
    for (int p = 1; p < size; p++)
    {
      int recv_chunksize = p < remaining ? chunksize + 1 : chunksize;
      int recv_offset = p * chunksize + (remaining > p ? p : remaining);
      std::complex<double> *recv_f = new std::complex<double>[recv_chunksize * window_size];
      MPI_Recv(recv_f, recv_chunksize * window_size, MPI_DOUBLE_COMPLEX, p, 0,
               MPI_COMM_WORLD, &status);

      // Copy from recv array to global f
      for (size_t i = 0; i < recv_chunksize; i++)
      {
        for (size_t k = 0; k < count; k++)
        {
          f[i + recv_offset][k] = recv_f[i * window_size + k];
          f[i + recv_offset][k + count] = recv_f[i * window_size + (k + count)];
        }
      }
    }
  }
  else
  {
    MPI_Send(local_f, local_chunksize * window_size,
             MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD);
  }
}

rarray<std::complex<double>, 2> fft_hybrid::stft_dft(
    rarray<std::complex<double>, 1> &vec,
    size_t window_size, size_t window_step)
{
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 1> dft;

  int rank, size, provided;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t chunk_size = spectrogram_size / size;
  size_t start_idx = rank * chunk_size;
  size_t end_idx = (rank == size - 1) ? spectrogram_size : (rank + 1) * chunk_size;

#pragma omp parallel for default(none) shared(start_idx, end_idx, vec, window_step, window_size, spectrogram)
  for (size_t i = start_idx; i < end_idx; i++)
  {
    dft = fft::dft(vec, i * window_step, window_size);
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = dft[j];
  }

  if (rank == 0)
    for (int p = 1; p < size; p++)
    {
      size_t received_chunk_size = (p == size - 1) ? spectrogram_size - p * chunk_size : chunk_size;
      MPI_Recv(&spectrogram[p * chunk_size][0],
               received_chunk_size * window_size, MPI_DOUBLE_COMPLEX, p, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  else
    MPI_Send(&spectrogram[start_idx][0],
             (end_idx - start_idx) * window_size, MPI_DOUBLE_COMPLEX, 0, 0,
             MPI_COMM_WORLD);

  return spectrogram;
}

rarray<std::complex<double>, 2> fft_hybrid::stft_fft(
    rarray<std::complex<double>, 1> &vec, size_t window_size, size_t window_step)
{
  size_t spectrogram_size = floor((vec.size() - window_size) / float(window_step)) + 1;
  rarray<std::complex<double>, 2> spectrogram(spectrogram_size, window_size);
  rarray<std::complex<double>, 1> fft;

  int rank, size, provided;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t chunk_size = spectrogram_size / size;
  size_t start_idx = rank * chunk_size;
  size_t end_idx = (rank == size - 1) ? spectrogram_size : (rank + 1) * chunk_size;

#pragma omp parallel for default(none) shared(start_idx, end_idx, vec, window_step, window_size, spectrogram)
  for (size_t i = start_idx; i < end_idx; i++)
  {
    fft = fft::fft(vec, i * window_step, window_size);
    for (size_t j = 0; j < window_size; j++)
      spectrogram[i][j] = fft[j];
  }

  if (rank == 0)
    for (int p = 1; p < size; p++)
    {
      size_t received_chunk_size = (p == size - 1) ? spectrogram_size - p * chunk_size : chunk_size;
      MPI_Recv(&spectrogram[p * chunk_size][0],
               received_chunk_size * window_size, MPI_DOUBLE_COMPLEX, p, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  else
    MPI_Send(&spectrogram[start_idx][0],
             (end_idx - start_idx) * window_size, MPI_DOUBLE_COMPLEX, 0, 0,
             MPI_COMM_WORLD);

  return spectrogram;
}
