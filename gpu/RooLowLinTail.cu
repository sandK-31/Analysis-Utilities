#include "RooFitPhotopeakKernels.hpp"

#include <cuda_runtime.h>

namespace {

__global__ void RooLowLinTailKernel(double *output, const double *x_vals,
                                    double mu, double slope,
                                    double inv_sqrt2_sigma, size_t n) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double y = x_vals[i] - mu;
    double lin = 1.0 + slope * y;
    if (lin < 0.0) {
      lin = 0.0;
    }
    double v = lin * erfc(y * inv_sqrt2_sigma);
    output[i] = v > 1e-300 ? v : 1e-300;
  }
}

} // namespace

void RooLowLinTail_launchKernel(double *output, const double *x_vals, double mu,
                                double slope, double inv_sqrt2_sigma,
                                size_t n) {
  const int threads = 256;
  const int blocks = (n + threads - 1) / threads;
  RooLowLinTailKernel<<<blocks, threads>>>(output, x_vals, mu, slope,
                                           inv_sqrt2_sigma, n);
}
