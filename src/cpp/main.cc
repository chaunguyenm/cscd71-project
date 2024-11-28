#include "fft.h"
#include <iostream>
#include <unistd.h>

int main(int argc, char *argv[])
{
  int helpFlag = 0;
  char *algorithm = NULL;
  char *parallel = NULL;
  char *outputFile = NULL;
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "ha:o:p:")) != -1)
    switch (c)
    {
    case 'h':
      helpFlag = 1;
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
    fprintf(stdout, "usage: %s [-h] [-a algorithm] [-p parallel] "
                    "[-o output] [num-samples] [window-size]\n",
            argv[0]);
    fprintf(stdout, "num-samples\tlong long, size of input signal\n");
    fprintf(stdout, "window-size\tlong long, less than num-samples, power of 2, size of stft window\n");
    fprintf(stdout, "-h\thelp\n");
    fprintf(stdout, "-a\tshort time fourier transform algorithm: "
                    "dft, fft (default), ff\n");
    fprintf(stdout, "-p\tparallel scheme: omp (default), "
		    "mpi (must be run with mpirun), "
		    "hybrid (omp + mpi, must be run with mpirun)\n");
    fprintf(stdout, "-o\toutput file: "
                    "[/path/to/file]\n");
    return 1;
  }

  if (optind + 2 != argc)
  {
    fprintf(stderr, "Invalid arguments. See -h for usage.\n");
    return 1;
  }

  long long num_samples = -1;
  long long window_size = -1;
  try
  {
    num_samples = std::stoll(*(argv + optind));
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
    window_size = std::stoll(*(argv + optind + 1));
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

  rarray<std::complex<double>, 1> signal(num_samples);
  for (size_t i = 0; i < num_samples; i++)
    signal[i] = std::complex<double>(double(i), 0.0);
  // std::cout << signal << "\n";

  if (outputFile != NULL)
  {
  }
  else
  {
    rarray<std::complex<double>, 2> stft;
    if (algorithm != NULL && strncmp((const char *)algorithm,
				     "dft", strlen("dft")) == 0)
    {
       // std::cout << "dft\n";
       if (parallel != NULL && strncmp((const char *)parallel,
				       "mpi", strlen("mpi")) == 0)
         stft = fft_mpi::stft_dft(signal, window_size, 1);
       else
         stft = fft::stft_dft(signal, window_size, 1);
    }
    else if (algorithm != NULL && strncmp((const char *)algorithm,
				          "fft", strlen("fft")) == 0)
    {
      // std::cout << "fft\n"; 
      if (parallel != NULL && strncmp((const char *)parallel,
				      "mpi", strlen("mpi")) == 0)
        stft = fft_mpi::stft_fft(signal, window_size, 1);
      else
        stft = fft::stft_fft(signal, window_size, 1);
    }
    else if (algorithm != NULL && strncmp((const char *)algorithm,
					  "ff", strlen("ff")) == 0)
    {
      // std::cout << "ff\n";
      stft = fft::stft_ff(signal, window_size, 1);
    }
    else
    {
      std::cout << "Unrecognized algorithm. See -h for usage.\n";
      return 1;
    }
    // std::cout << stft << "\n";
  }
  
  return 0;
}

