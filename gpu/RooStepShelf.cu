#include "RooFitPhotopeakKernels.hpp"

#include <cuda_runtime.h>

namespace {

__device__ inline double SigmoidNegDev(double z) {
  if (z < 0.0) {
    return 1.0 / (1.0 + exp(z));
  }
  return exp(-z) / (1.0 + exp(-z));
}

__global__ void RooStepShelfKernel(double *output, const double *x_vals,
                                   double mu, double inv_sigma, size_t n) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double z = (x_vals[i] - mu) * inv_sigma;
    double s = SigmoidNegDev(z);
    double v = s * s;
    output[i] = v > 1e-300 ? v : 1e-300;
  }
}

} // namespace

void RooStepShelf_launchKernel(double *output, const double *x_vals, double mu,
                               double inv_sigma, size_t n) {
  const int threads = 256;
  const int blocks = (n + threads - 1) / threads;
  RooStepShelfKernel<<<blocks, threads>>>(output, x_vals, mu, inv_sigma, n);
}
