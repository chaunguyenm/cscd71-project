#include "stft.h"
#include "rarray"
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <mpi.h>
#include <map>
#include <string>
#include <functional>
#include <iostream>
#include <cstring>

rarray<std::complex<double>, 2> compute_stft(
    rarray<std::complex<double>, 1> signal,
    unsigned long window_size,
    const char *algorithm_cstr,
    const char *parallel_cstr)
{
  using namespace std;
  using STFTFunc = function<rarray<complex<double>, 2>(rarray<complex<double>, 1> &, size_t, size_t)>;

  string algorithm = algorithm_cstr ? string(algorithm_cstr) : "";
  string parallel = parallel_cstr ? string(parallel_cstr) : "";

  auto normalize = [](string s)
  {
    for (auto &c : s)
      c = tolower(c);
    return s;
  };
  algorithm = normalize(algorithm);
  parallel = normalize(parallel);

  // key = (algorithm, parallel)
  const map<pair<string, string>, STFTFunc> stft_dispatch = {
      {{"dft", "mpi"}, stft_mpi::stft_dft},
      {{"dft", "omp"}, stft::stft_dft},
      {{"dft", "hybrid"}, stft_mpi::stft_dft},

      {{"fft", "mpi"}, stft_mpi::stft_fft},
      {{"fft", "omp"}, stft::stft_fft},
      {{"fft", "hybrid"}, stft_mpi::stft_fft},

      {{"ff", ""}, stft::stft_ff}, // no parallel options

      {{"qpff", "mpi"}, stft_mpi::stft_qpff},
      {{"qpff", "omp"}, stft::stft_qpff},
      {{"qpff", "hybrid"}, stft_mpi::stft_qpff},
  };

  rarray<complex<double>, 2> stft;

  auto it = stft_dispatch.find({algorithm, parallel});
  if (it == stft_dispatch.end() && parallel.empty())
  {
    it = stft_dispatch.find({algorithm, ""});
  }

  if (it != stft_dispatch.end())
  {
    rarray<complex<double>, 1> sig_copy = signal;
    return it->second(sig_copy, window_size, 1);
  }

  std::cerr << "Unrecognized algorithm or parallel scheme. See -h for usage.\n";
  return stft;
}

int main(int argc, char *argv[])
{
  int helpFlag = 0;
  int silentFlag = 0;
  char *algorithm = NULL;
  char *parallel = NULL;
  char *outputFile = NULL;
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "hsa:o:p:")) != -1)
    switch (c)
    {
    case 'h':
      helpFlag = 1;
      break;
    case 's':
      silentFlag = 1;
      break;
    case 'a':
      algorithm = optarg;
      break;
    case 'o':
      outputFile = optarg;
      break;
    case 'p':
      parallel = optarg;
      break;
    case '?':
      if (optopt == 'a' || optopt == 'o')
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint(optopt))
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
      return 1;
    default:
      abort();
    }

  if (helpFlag)
  {
    fprintf(stdout, "usage: %s [-h] [-s] [-a algorithm] [-p parallel] "
                    "[-o output] [num-samples] [window-size]\n",
            argv[0]);
    fprintf(stdout, "num-samples\tlong, size of input signal\n");
    fprintf(stdout, "window-size\tlong, less than num-samples, power of 2, size of stft window\n");
    fprintf(stdout, "-h\thelp\n");
    fprintf(stdout, "-s\tsilent (do not print any output)\n");
    fprintf(stdout, "-a\tshort time fourier transform algorithm: "
                    "dft, fft (default), ff\n");
    fprintf(stdout, "-p\tparallel scheme: omp (default), "
                    "mpi (must be run with mpirun)\n");
    fprintf(stdout, "-o\toutput file: "
                    "[/path/to/file]\n");
    return 1;
  }

  if (optind + 2 != argc)
  {
    fprintf(stderr, "Invalid arguments. See -h for usage.\n");
    return 1;
  }

  unsigned long num_samples = 0;
  unsigned long window_size = 0;
  try
  {
    num_samples = std::stoul(*(argv + optind));
  }
  catch (std::invalid_argument const &arg)
  {
    fprintf(stderr, "Invalid arguments [num-samples]. See -h for usage.\n");
    return 1;
  }
  catch (std::out_of_range const &arg)
  {
    fprintf(stderr, "Out of range arguments [num-samples]. See -h for usage.\n");
    return 1;
  }
  try
  {
    window_size = std::stoul(*(argv + optind + 1));
  }
  catch (std::invalid_argument const &arg)
  {
    fprintf(stderr, "Invalid arguments [window-size]. See -h for usage.\n");
    return 1;
  }
  catch (std::out_of_range const &arg)
  {
    fprintf(stderr, "Out of range arguments [window-size]. See -h for usage.\n");
    return 1;
  }

  if (num_samples < window_size)
  {
    fprintf(stderr, "num-samples must be greater than window-size. See -h for usage.\n");
    return 1;
  }
  if (window_size == 0 || num_samples == 0)
  {
    fprintf(stderr, "num-samples and window-size must be positive. See -h for usage.\n");
    return 1;
  }
  if ((window_size & (window_size - 1)) != 0)
  {
    fprintf(stderr, "window-size must be power of 2. See -h for usage.\n");
    return 1;
  }

  rarray<std::complex<double>, 1> signal(num_samples);
  for (size_t i = 0; i < num_samples; i++)
    signal[i] = std::complex<double>(double(i), 0.0);
  // std::cout << signal << "\n";

  int rank, size, provided;
  if (parallel != NULL && strncmp(parallel, "mpi", strlen("mpi")) == 0)
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }
  else if (parallel != NULL && strncmp(parallel, "hybrid", strlen("hybrid")) == 0)
  {
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }
  else
  {
    rank = 0;
    size = 0;
    provided = 0;
  }

  rarray<std::complex<double>, 2> stft;
  stft = compute_stft(signal, window_size, (const char *)algorithm, (const char *)parallel);

  if (rank == 0)
  {
    if (outputFile != NULL)
    {
      std::ofstream outfile;
      outfile.open(outputFile);
      if (outfile.is_open())
      {
        outfile << stft << "\n";
        outfile.close();
      }
      else
      {
        fprintf(stderr, "Error opening output file.\n");
        return 1;
      }
    }
    else if (silentFlag == 0)
      std::cout << stft << "\n";
  }

  if (parallel != NULL && strncmp(parallel, "mpi", strlen("mpi")) == 0)
    MPI_Finalize();
  else if (parallel != NULL && strncmp(parallel, "hybrid", strlen("hybrid")) == 0)
    MPI_Finalize();
  return 0;
}
