#include "RooFitPhotopeakPdfs.hpp"

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <cmath>

#ifdef AU_ROOFIT_BACKEND_CUDA
#include "RooFitPhotopeakKernels.hpp"
#include <cuda_runtime.h>

namespace {
bool IsDevicePointer(const void *ptr) {
  cudaPointerAttributes attrs;
  cudaError_t err = cudaPointerGetAttributes(&attrs, ptr);
  if (err != cudaSuccess) {
    // Clear the error state — querying a host pointer is not a hard failure.
    cudaGetLastError();
    return false;
  }
  return attrs.type == cudaMemoryTypeDevice ||
         attrs.type == cudaMemoryTypeManaged;
}
} // namespace
#endif

ClassImp(RooStepShelf);
ClassImp(RooLowExpTail);
ClassImp(RooLowLinTail);
ClassImp(RooHighExpTail);

namespace {

Double_t SoftPlusNeg(Double_t z) {
  if (z > 0)
    return std::log1p(std::exp(-z));
  return -z + std::log1p(std::exp(z));
}

Double_t SigmoidNeg(Double_t z) {
  if (z < 0)
    return 1.0 / (1.0 + std::exp(z));
  return std::exp(-z) / (1.0 + std::exp(-z));
}

Double_t StepAntideriv(Double_t z) { return -SoftPlusNeg(z) + SigmoidNeg(z); }

// exp(a) * erfc(b) computed without intermediate overflow. In the exponential
// tails a = y/tau and b = y/(sqrt(2)*sigma) always share the sign of y, and
// b = a * tau/(sqrt(2)*sigma) = a * ratio/sqrt(2) with ratio >= 1, so whenever
// exp(a) would overflow (a large positive) b is also large positive and
// erfc(b) underflows -- the mathematical product tends to 0. Fold the two
// exponentials together and use the large-argument expansion of the scaled
// complementary error function erfcx(b) = exp(b*b)*erfc(b) ~
//   1/(b*sqrt(pi)) * (1 - 1/(2 b^2) + 3/(4 b^4) - 15/(8 b^6)).
// The direct branch (a <= 700) keeps exp(a) finite (exp(700) ~ 1e304).
Double_t ExpErfc(Double_t a, Double_t b) {
  if (a <= 700.0)
    return std::exp(a) * std::erfc(b);
  Double_t t = 1.0 / (b * b);
  Double_t erfcx = (1.0 / (b * std::sqrt(M_PI))) *
                   (1.0 - 0.5 * t + 0.75 * t * t - 1.875 * t * t * t);
  return std::exp(a - b * b) * erfcx;
}

// Smallest density any component pdf is allowed to return. A pdf that evaluates
// to exactly 0 across the fit range integrates to 0, and RooFit's extended
// Range() projection then forms a 0/0 ratio (RooRatio::rangeProj) -> NaN. A
// tiny positive floor keeps every component's normalized integral well-defined
// while being utterly negligible wherever the real density is non-zero.
const Double_t kDensityFloor = 1e-300;

Double_t ExpTailDensity(Double_t y, Double_t sigma, Double_t tau) {
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  Double_t v = ExpErfc(y / tau, y / sqrt2_sigma);
  return v > kDensityFloor ? v : kDensityFloor;
}

Double_t ExpTailAntideriv(Double_t y, Double_t sigma, Double_t tau) {
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  Double_t offset = sigma * sigma / tau;
  Double_t gauss_corr = std::exp(sigma * sigma / (2.0 * tau * tau));
  return tau * (ExpErfc(y / tau, y / sqrt2_sigma) +
                gauss_corr * std::erf((y - offset) / sqrt2_sigma));
}

Double_t LowLinDensity(Double_t y, Double_t sigma, Double_t slope) {
  Double_t lin = 1.0 + slope * y;
  if (lin < 0.0)
    lin = 0.0;
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  Double_t v = lin * std::erfc(y / sqrt2_sigma);
  return v > kDensityFloor ? v : kDensityFloor;
}

Double_t LowLinAntideriv(Double_t y, Double_t sigma, Double_t slope) {
  Double_t sqrt2_sigma = std::sqrt(2.0) * sigma;
  Double_t beta = 1.0 / sqrt2_sigma;
  Double_t by = beta * y;
  Double_t erfc_by = std::erfc(by);
  Double_t erf_by = std::erf(by);
  Double_t gauss = std::exp(-(by * by));
  Double_t inv_beta_sqrt_pi = 1.0 / (beta * std::sqrt(M_PI));
  Double_t flat = y * erfc_by - gauss * inv_beta_sqrt_pi;
  Double_t lin = 0.5 * y * y * erfc_by - y * gauss * inv_beta_sqrt_pi / 2.0 +
                 erf_by / (4.0 * beta * beta);
  return flat + slope * lin;
}

Double_t LowLinIntegral(Double_t y_lo, Double_t y_hi, Double_t sigma,
                        Double_t slope) {
  Double_t eff_lo = y_lo;
  Double_t eff_hi = y_hi;
  if (slope > 1e-10) {
    Double_t threshold = -1.0 / slope;
    if (threshold > eff_lo)
      eff_lo = threshold;
    if (eff_lo >= eff_hi)
      return 1e-300;
  } else if (slope < -1e-10) {
    Double_t threshold = -1.0 / slope;
    if (threshold < eff_hi)
      eff_hi = threshold;
    if (eff_lo >= eff_hi)
      return 1e-300;
  }
  Double_t val = LowLinAntideriv(eff_hi, sigma, slope) -
                 LowLinAntideriv(eff_lo, sigma, slope);
  if (!std::isfinite(val) || val < 1e-300)
    return 1e-300;
  return val;
}

} // namespace

RooStepShelf::RooStepShelf(const char *name, const char *title, RooAbsReal &x,
                           RooAbsReal &mu, RooAbsReal &sigma)
    : RooAbsPdf(name, title), x_("x", "x", this, x), mu_("mu", "mu", this, mu),
      sigma_("sigma", "sigma", this, sigma) {}

RooStepShelf::RooStepShelf(const RooStepShelf &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_) {}

Double_t RooStepShelf::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return kDensityFloor;
  Double_t z = ((Double_t)x_ - (Double_t)mu_) / sigma;
  Double_t s = SigmoidNeg(z);
  Double_t v = s * s;
  return v > kDensityFloor ? v : kDensityFloor;
}

void RooStepShelf::doEval(RooFit::EvalContext &ctx) const {
  std::span<const double> x_vals = ctx.at(x_);
  std::span<double> output = ctx.output();
  size_t n = output.size();
  Double_t sigma = ctx.at(sigma_)[0];
  Double_t mu = ctx.at(mu_)[0];
  if (sigma <= 0) {
    for (size_t i = 0; i < n; ++i)
      output[i] = kDensityFloor;
    return;
  }
  Double_t inv_sigma = 1.0 / sigma;
#ifdef AU_ROOFIT_BACKEND_CUDA
  if (IsDevicePointer(output.data())) {
    RooStepShelf_launchKernel(output.data(), x_vals.data(), mu, inv_sigma, n);
    return;
  }
#endif
#pragma omp parallel for simd
  for (size_t i = 0; i < n; ++i) {
    Double_t z = (x_vals[i] - mu) * inv_sigma;
    Double_t s = SigmoidNeg(z);
    Double_t v = s * s;
    output[i] = v > kDensityFloor ? v : kDensityFloor;
  }
}

Int_t RooStepShelf::getAnalyticalIntegral(RooArgSet &allVars,
                                          RooArgSet &analVars,
                                          const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooStepShelf::analyticalIntegral(Int_t code,
                                          const char *rangeName) const {
  if (code != 1)
    return 1e-300;
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return 1e-300;
  Double_t mu = (Double_t)mu_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t z_lo = (x_lo - mu) / sigma;
  Double_t z_hi = (x_hi - mu) / sigma;
  Double_t val = sigma * (StepAntideriv(z_hi) - StepAntideriv(z_lo));
  if (!std::isfinite(val) || val < 1e-300)
    return 1e-300;
  return val;
}

RooLowExpTail::RooLowExpTail(const char *name, const char *title, RooAbsReal &x,
                             RooAbsReal &mu, RooAbsReal &sigma, RooAbsReal &tau)
    : RooAbsPdf(name, title), x_("x", "x", this, x), mu_("mu", "mu", this, mu),
      sigma_("sigma", "sigma", this, sigma), tau_("tau", "tau", this, tau) {}

RooLowExpTail::RooLowExpTail(const RooLowExpTail &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_),
      tau_("tau", this, other.tau_) {}

Double_t RooLowExpTail::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  Double_t ratio = (Double_t)tau_;
  if (sigma <= 0 || ratio <= 0)
    return 0.0;
  Double_t tau = ratio * sigma;
  Double_t y = (Double_t)x_ - (Double_t)mu_;
  return ExpTailDensity(y, sigma, tau);
}

void RooLowExpTail::doEval(RooFit::EvalContext &ctx) const {
  std::span<const double> x_vals = ctx.at(x_);
  std::span<double> output = ctx.output();
  size_t n = output.size();
  Double_t sigma = ctx.at(sigma_)[0];
  Double_t mu = ctx.at(mu_)[0];
  Double_t ratio = ctx.at(tau_)[0];
  if (sigma <= 0 || ratio <= 0) {
    for (size_t i = 0; i < n; ++i)
      output[i] = kDensityFloor;
    return;
  }
  Double_t tau = ratio * sigma;
  Double_t inv_tau = 1.0 / tau;
  Double_t inv_sqrt2_sigma = 1.0 / (std::sqrt(2.0) * sigma);
#ifdef AU_ROOFIT_BACKEND_CUDA
  if (IsDevicePointer(output.data())) {
    RooLowExpTail_launchKernel(output.data(), x_vals.data(), mu, inv_tau,
                               inv_sqrt2_sigma, n);
    return;
  }
#endif
#pragma omp parallel for simd
  for (size_t i = 0; i < n; ++i) {
    Double_t y = x_vals[i] - mu;
    Double_t v = ExpErfc(y * inv_tau, y * inv_sqrt2_sigma);
    output[i] = v > kDensityFloor ? v : kDensityFloor;
  }
}

Int_t RooLowExpTail::getAnalyticalIntegral(RooArgSet &allVars,
                                           RooArgSet &analVars,
                                           const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooLowExpTail::analyticalIntegral(Int_t code,
                                           const char *rangeName) const {
  if (code != 1)
    return 1e-300;
  Double_t sigma = (Double_t)sigma_;
  Double_t ratio = (Double_t)tau_;
  if (sigma <= 0 || ratio <= 0)
    return 1e-300;
  Double_t tau = ratio * sigma;
  Double_t mu = (Double_t)mu_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t y_lo = x_lo - mu;
  Double_t y_hi = x_hi - mu;
  Double_t val =
      ExpTailAntideriv(y_hi, sigma, tau) - ExpTailAntideriv(y_lo, sigma, tau);
  if (!std::isfinite(val) || val < 1e-300)
    return 1e-300;
  return val;
}

RooLowLinTail::RooLowLinTail(const char *name, const char *title, RooAbsReal &x,
                             RooAbsReal &mu, RooAbsReal &sigma,
                             RooAbsReal &slope)
    : RooAbsPdf(name, title), x_("x", "x", this, x), mu_("mu", "mu", this, mu),
      sigma_("sigma", "sigma", this, sigma),
      slope_("slope", "slope", this, slope) {}

RooLowLinTail::RooLowLinTail(const RooLowLinTail &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_),
      slope_("slope", this, other.slope_) {}

Double_t RooLowLinTail::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return 0.0;
  Double_t y = (Double_t)x_ - (Double_t)mu_;
  Double_t slope = (Double_t)slope_;
  return LowLinDensity(y, sigma, slope);
}

void RooLowLinTail::doEval(RooFit::EvalContext &ctx) const {
  std::span<const double> x_vals = ctx.at(x_);
  std::span<double> output = ctx.output();
  size_t n = output.size();
  Double_t sigma = ctx.at(sigma_)[0];
  Double_t mu = ctx.at(mu_)[0];
  Double_t slope = ctx.at(slope_)[0];
  if (sigma <= 0) {
    for (size_t i = 0; i < n; ++i)
      output[i] = kDensityFloor;
    return;
  }
  Double_t inv_sqrt2_sigma = 1.0 / (std::sqrt(2.0) * sigma);
#ifdef AU_ROOFIT_BACKEND_CUDA
  if (IsDevicePointer(output.data())) {
    RooLowLinTail_launchKernel(output.data(), x_vals.data(), mu, slope,
                               inv_sqrt2_sigma, n);
    return;
  }
#endif
#pragma omp parallel for simd
  for (size_t i = 0; i < n; ++i) {
    Double_t y = x_vals[i] - mu;
    Double_t lin = 1.0 + slope * y;
    if (lin < 0.0)
      lin = 0.0;
    Double_t v = lin * std::erfc(y * inv_sqrt2_sigma);
    output[i] = v > kDensityFloor ? v : kDensityFloor;
  }
}

Int_t RooLowLinTail::getAnalyticalIntegral(RooArgSet &allVars,
                                           RooArgSet &analVars,
                                           const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooLowLinTail::analyticalIntegral(Int_t code,
                                           const char *rangeName) const {
  if (code != 1)
    return 0.0;
  Double_t sigma = (Double_t)sigma_;
  if (sigma <= 0)
    return 0.0;
  Double_t mu = (Double_t)mu_;
  Double_t slope = (Double_t)slope_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t y_lo = x_lo - mu;
  Double_t y_hi = x_hi - mu;
  return LowLinIntegral(y_lo, y_hi, sigma, slope);
}

RooHighExpTail::RooHighExpTail(const char *name, const char *title,
                               RooAbsReal &x, RooAbsReal &mu, RooAbsReal &sigma,
                               RooAbsReal &tau)
    : RooAbsPdf(name, title), x_("x", "x", this, x), mu_("mu", "mu", this, mu),
      sigma_("sigma", "sigma", this, sigma), tau_("tau", "tau", this, tau) {}

RooHighExpTail::RooHighExpTail(const RooHighExpTail &other, const char *name)
    : RooAbsPdf(other, name), x_("x", this, other.x_),
      mu_("mu", this, other.mu_), sigma_("sigma", this, other.sigma_),
      tau_("tau", this, other.tau_) {}

Double_t RooHighExpTail::evaluate() const {
  Double_t sigma = (Double_t)sigma_;
  Double_t ratio = (Double_t)tau_;
  if (sigma <= 0 || ratio <= 0)
    return 0.0;
  Double_t tau = ratio * sigma;
  Double_t z = (Double_t)mu_ - (Double_t)x_;
  return ExpTailDensity(z, sigma, tau);
}

void RooHighExpTail::doEval(RooFit::EvalContext &ctx) const {
  std::span<const double> x_vals = ctx.at(x_);
  std::span<double> output = ctx.output();
  size_t n = output.size();
  Double_t sigma = ctx.at(sigma_)[0];
  Double_t mu = ctx.at(mu_)[0];
  Double_t ratio = ctx.at(tau_)[0];
  if (sigma <= 0 || ratio <= 0) {
    for (size_t i = 0; i < n; ++i)
      output[i] = kDensityFloor;
    return;
  }
  Double_t tau = ratio * sigma;
  Double_t inv_tau = 1.0 / tau;
  Double_t inv_sqrt2_sigma = 1.0 / (std::sqrt(2.0) * sigma);
#ifdef AU_ROOFIT_BACKEND_CUDA
  if (IsDevicePointer(output.data())) {
    RooHighExpTail_launchKernel(output.data(), x_vals.data(), mu, inv_tau,
                                inv_sqrt2_sigma, n);
    return;
  }
#endif
#pragma omp parallel for simd
  for (size_t i = 0; i < n; ++i) {
    Double_t z = mu - x_vals[i];
    Double_t v = ExpErfc(z * inv_tau, z * inv_sqrt2_sigma);
    output[i] = v > kDensityFloor ? v : kDensityFloor;
  }
}

Int_t RooHighExpTail::getAnalyticalIntegral(RooArgSet &allVars,
                                            RooArgSet &analVars,
                                            const char * /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_))
    return 1;
  return 0;
}

Double_t RooHighExpTail::analyticalIntegral(Int_t code,
                                            const char *rangeName) const {
  if (code != 1)
    return 1e-300;
  Double_t sigma = (Double_t)sigma_;
  Double_t ratio = (Double_t)tau_;
  if (sigma <= 0 || ratio <= 0)
    return 1e-300;
  Double_t tau = ratio * sigma;
  Double_t mu = (Double_t)mu_;
  Double_t x_lo = x_.min(rangeName);
  Double_t x_hi = x_.max(rangeName);
  Double_t z_lo = mu - x_lo;
  Double_t z_hi = mu - x_hi;
  Double_t val =
      ExpTailAntideriv(z_lo, sigma, tau) - ExpTailAntideriv(z_hi, sigma, tau);
  if (!std::isfinite(val) || val < 1e-300)
    return 1e-300;
  return val;
}
