#ifndef FITTINGUTILS_H
#define FITTINGUTILS_H

#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iomanip>

// Forward declaration (defined in InteractiveFitEditor.hpp)
Bool_t LaunchInteractiveFitEditor(TH1 *hist, TF1 *fit_func, Double_t range_low,
                                  Double_t range_high, Int_t num_peaks = 1,
                                  const TString &info_label = "");

namespace FittingFunctions {
Double_t Gaussian(Double_t *x, Double_t *par);
Double_t LinearBackground(Double_t *x, Double_t *par);
Double_t Step(Double_t *x, Double_t *par);
Double_t LowTail(Double_t *x, Double_t *par);
Double_t HighTail(Double_t *x, Double_t *par);
Double_t PeakFunction(Double_t *x, Double_t *par);
Double_t DoublePeakFunction(Double_t *x, Double_t *par);
Double_t TriplePeakFunction(Double_t *x, Double_t *par);
} // namespace FittingFunctions

struct PeakFitResult {
  Float_t mu = -1, mu_error = -1;
  Float_t sigma = -1, sigma_error = -1;
  Float_t gaus_amplitude = -1, gaus_amplitude_error = -1;
  Float_t step_amplitude = -1, step_amplitude_error = -1;
  Float_t low_exp_tail_amplitude = -1, low_exp_tail_amplitude_error = -1;
  Float_t low_exp_tail_ratio = -1, low_exp_tail_ratio_error = -1;
  Float_t low_lin_tail_amplitude = -1, low_lin_tail_amplitude_error = -1;
  Float_t low_lin_tail_slope = -1, low_lin_tail_slope_error = -1;
  Float_t high_exp_tail_amplitude = -1, high_exp_tail_amplitude_error = -1;
  Float_t high_exp_tail_ratio = -1, high_exp_tail_ratio_error = -1;
};

struct FitResult {
  std::vector<PeakFitResult> peaks; // 1-3 entries supported
  Float_t bkg_constant = -1, bkg_constant_error = -1;
  Float_t lin_bkg_slope = -1, lin_bkg_slope_error = -1;
  Float_t reduced_chi2 = -1;
  Bool_t valid = kFALSE;
};

class FittingUtils {
private:
  TF1 *fit_function_;
  TH1 *working_hist_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;

  Bool_t use_flat_background_;
  Bool_t use_step_;
  Bool_t use_low_exp_tail_;
  Bool_t use_low_lin_tail_;
  Bool_t use_high_exp_tail_;

  Bool_t use_manual_init_;
  Bool_t interactive_;
  std::vector<Double_t> manual_params_;

  Double_t EstimateBackground();
  Double_t ClampToBounds(Int_t param_index, Double_t value);

  void SaveInteractiveParams(const TString &input_name,
                             const TString &peak_name);
  Bool_t LoadInteractiveParams(const TString &input_name,
                               const TString &peak_name);

  void SortPeaksByMu(Int_t num_peaks);
  void AppendPeakGraphs(std::vector<TGraph *> &components, Int_t param_offset,
                        Style_t line_style, TF1 *background, Int_t npts,
                        Double_t x_step);

public:
  FittingUtils(TH1 *working_hist, Float_t fit_range_low, Float_t fit_range_high,
               Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
               Bool_t use_low_exp_tail = kFALSE,
               Bool_t use_low_lin_tail = kFALSE,
               Bool_t use_high_exp_tail = kFALSE);
  ~FittingUtils();

  void SetBackgroundModel(Bool_t use_flat_background) {
    use_flat_background_ = use_flat_background;
  }
  void SetStep(Bool_t use_step = kTRUE) { use_step_ = use_step; }
  void SetLowExpTail(Bool_t use_low_exp_tail = kTRUE) {
    use_low_exp_tail_ = use_low_exp_tail;
  }
  void SetLowLinTail(Bool_t use_low_lin_tail = kTRUE) {
    use_low_lin_tail_ = use_low_lin_tail;
  }
  void SetHighExpTail(Bool_t use_high_exp_tail = kTRUE) {
    use_high_exp_tail_ = use_high_exp_tail;
  }
  void SetInteractive(Bool_t interactive = kTRUE) {
    interactive_ = interactive;
  }

  void SetManualParameters(const std::vector<Double_t> &params);
  void SetManualParameter(Int_t index, Double_t value);
  void ClearManualParameters() {
    use_manual_init_ = kFALSE;
    manual_params_.clear();
  }

  TF1 *GetFitFunction() { return fit_function_; }
  void SetFitFunction(TF1 *func) { fit_function_ = func; }

  void PlotFitSinglePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitDoublePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitTriplePeak(const TString input_name, const TString peak_name,
                         const TString label = "");

  FitResult FitSinglePeak(const TString input_name, const TString peak_name);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          const PeakFitResult &constrained_peak,
                          Double_t mu2_init);
  FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                          const FitResult &constrained_peaks,
                          Double_t mu3_init);
};

#endif
