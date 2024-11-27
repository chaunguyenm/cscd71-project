#include "fft.h"
#include <iostream>
#include <chrono>

int main()
{
  size_t size = 20;
  rarray<std::complex<double>, 1> signal(size);
  for (size_t i = 0; i < size; i++)
    signal[i] = std::complex<double>(double(i), 0.0);
  std::cout << signal << "\n";

  rarray<std::complex<double>, 2> stft_dft = fft::stft_dft(signal, 4, 1);
  std::cout << stft_dft << "\n";

  rarray<std::complex<double>, 2> stft_fft = fft::stft_fft(signal, 4, 1);
  std::cout << stft_fft << "\n";

  rarray<std::complex<double>, 2> stft_ff = fft::stft_ff(signal, 4, 1);
  std::cout << stft_ff << "\n";

  return 0; 
}
