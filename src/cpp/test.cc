#include <complex>
#include "rarray"

int main()
{
  rarray<std::complex<double>, 1> signal(2);
  signal[0] = std::complex<double>(0., 0.);
  signal[1] = std::complex<double>(1., 0.);
  std::cout << signal << "\n";
  return 0;
}