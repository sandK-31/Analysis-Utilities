#ifndef ROOFITPHOTOPEAKKERNELS_H
#define ROOFITPHOTOPEAKKERNELS_H

#include <cstddef>

// Host-callable launch wrappers for the RooFit photopeak PDF CUDA kernels.
// Each launcher takes device pointers and the hoisted constants computed
// once per Minuit evaluation, dispatches the kernel, and returns once the
// stream has accepted the launch (no implicit synchronization).

void RooLowExpTail_launchKernel(double *output, const double *x_vals, double mu,
                                double inv_tau, double inv_sqrt2_sigma,
                                size_t n);

void RooStepShelf_launchKernel(double *output, const double *x_vals, double mu,
                               double inv_sigma, size_t n);

void RooLowLinTail_launchKernel(double *output, const double *x_vals, double mu,
                                double slope, double inv_sqrt2_sigma, size_t n);

void RooHighExpTail_launchKernel(double *output, const double *x_vals,
                                 double mu, double inv_tau,
                                 double inv_sqrt2_sigma, size_t n);

#endif
