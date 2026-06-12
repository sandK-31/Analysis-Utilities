#include "RooFitPhotopeakKernels.hpp"

#include <cuda_runtime.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

// exp(a)*erfc(b) without intermediate overflow; see ExpErfc in
// src/RooFitPhotopeakPdfs.cpp for the derivation. In the tail a and b share the
// sign of z and b >= a/sqrt(2), so when exp(a) would overflow erfc(b)
// underflows and the product tends to 0.
__device__ double ExpErfc(double a, double b) {
  if (a <= 700.0)
    return exp(a) * erfc(b);
  double t = 1.0 / (b * b);
  double erfcx = (1.0 / (b * sqrt(M_PI))) *
                 (1.0 - 0.5 * t + 0.75 * t * t - 1.875 * t * t * t);
  return exp(a - b * b) * erfcx;
}

__global__ void RooHighExpTailKernel(double *output, const double *x_vals,
                                     double mu, double inv_tau,
                                     double inv_sqrt2_sigma, size_t n) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double z = mu - x_vals[i];
    double v = ExpErfc(z * inv_tau, z * inv_sqrt2_sigma);
    output[i] = v > 1e-300 ? v : 1e-300;
  }
}

} // namespace

void RooHighExpTail_launchKernel(double *output, const double *x_vals,
                                 double mu, double inv_tau,
                                 double inv_sqrt2_sigma, size_t n) {
  const int threads = 256;
  const int blocks = (n + threads - 1) / threads;
  RooHighExpTailKernel<<<blocks, threads>>>(output, x_vals, mu, inv_tau,
                                            inv_sqrt2_sigma, n);
}
