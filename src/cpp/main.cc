#include "fft.h"
#include <iostream>
#include <chrono>

int main()
{
  size_t size = 16;
  rarray<std::complex<double>, 1> signal(size);
  for (size_t i = 0; i < size; i++)
    signal[i] = std::complex<double>(double(i), 0.0);
  // std::cout << signal << "\n";

  // rarray<std::complex<double>, 2> stft_dft = fft::stft_dft(signal, 8, 2);
  // std::cout << stft_dft << "\n";

  rarray<std::complex<double>, 2> stft_fft = fft::stft_fft(signal, 8, 1);
  std::cout << stft_fft << "\n";

  fft::stft_ff(signal, 8, 1);
  // std::cout << stft_fft << "\n";
  // std::vector<std::complex<double>> signal = {
  //   std::complex<double>(0.0, 0.0),
  //   std::complex<double>(1.0, 0.0),
  //   std::complex<double>(2.0, 0.0),
  //   std::complex<double>(3.0, 0.0),
  //   std::complex<double>(4.0, 0.0),
  //   std::complex<double>(5.0, 0.0),
  //   std::complex<double>(6.0, 0.0),
  //   std::complex<double>(7.0, 0.0),
  //   std::complex<double>(8.0, 0.0),
  //   std::complex<double>(9.0, 0.0),
  //   std::complex<double>(10.0, 0.0),
  //   std::complex<double>(11.0, 0.0),
  //   std::complex<double>(12.0, 0.0),
  //   std::complex<double>(13.0, 0.0),
  //   std::complex<double>(14.0, 0.0),
  //   std::complex<double>(15.0, 0.0)
  // };
  // std::vector<std::vector<std::complex<double>>> spectrogram1 = fft::stft_dft(signal, 8, 2);
  // std::cout << "Spectrogram matrix:" << std::endl;
  // for (const auto& row : spectrogram1) {
  //   for (const auto& element : row) {
  //     std::cout << "(" << element.real() << ", " << element.imag() << ") ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::vector<std::vector<std::complex<double>>> spectrogram2 = fft::stft_fft(signal, 8, 2);
  // std::cout << "Spectrogram matrix:" << std::endl;
  // for (const auto& row : spectrogram2) {
  //   for (const auto& element : row) {
  //     std::cout << "(" << element.real() << ", " << element.imag() << ") ";
  //   }
  //   std::cout << std::endl;
  // }
  // return 0; 
}
