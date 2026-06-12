#ifndef ROOFITUTILS_H
#define ROOFITUTILS_H

#include "FittingUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitPhotopeakPdfs.hpp"

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooCmdArg.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooMsgService.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>

#include <TBranch.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

// Pick the best-available RooFit evaluation backend. Defaults to CPU (batched
// via doEval). Compile with -DAU_ROOFIT_BACKEND_CUDA=1 once ROOT is built with
// CUDA support to switch to GPU evaluation.
inline RooCmdArg BestAvailableBackend() {
#if defined(AU_ROOFIT_BACKEND_CUDA) && AU_ROOFIT_BACKEND_CUDA
  return RooFit::EvalBackend::Cuda();
#else
  return RooFit::EvalBackend::Cpu();
#endif
}

namespace RooFitFunctions {
RooAbsPdf *MakeGaussian(const TString &name, RooRealVar &x, RooRealVar &mu,
                        RooRealVar &sigma);
RooAbsPdf *MakeStepShelf(const TString &name, RooRealVar &x, RooRealVar &mu,
                         RooRealVar &sigma);
RooAbsPdf *MakeLowExpTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                          RooRealVar &sigma, RooRealVar &tau);
RooAbsPdf *MakeLowLinTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                          RooRealVar &sigma, RooRealVar &slope);
RooAbsPdf *MakeHighExpTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                           RooRealVar &sigma, RooRealVar &tau);
RooAbsPdf *MakeLinearBackground(const TString &name, RooRealVar &x,
                                RooRealVar &slope);
} // namespace RooFitFunctions

struct RooFitPeakModel {
  RooRealVar *mu = nullptr;
  RooRealVar *sigma = nullptr;
  RooRealVar *gaus_yield = nullptr;
  RooRealVar *ratio_step = nullptr;
  RooRealVar *ratio_low_exp = nullptr;
  RooRealVar *tau_low_exp = nullptr;
  RooRealVar *ratio_low_lin = nullptr;
  RooRealVar *slope_low_lin = nullptr;
  RooRealVar *ratio_high_exp = nullptr;
  RooRealVar *tau_high_exp = nullptr;

  RooAbsPdf *gauss_pdf = nullptr;
  RooAbsPdf *step_pdf = nullptr;
  RooAbsPdf *low_exp_pdf = nullptr;
  RooAbsPdf *low_lin_pdf = nullptr;
  RooAbsPdf *high_exp_pdf = nullptr;

  RooFormulaVar *step_yield = nullptr;
  RooFormulaVar *low_exp_yield = nullptr;
  RooFormulaVar *low_lin_yield = nullptr;
  RooFormulaVar *high_exp_yield = nullptr;
};

struct RooFitBackgroundModel {
  RooRealVar *bkg_yield = nullptr;
  RooRealVar *bkg_slope = nullptr;
  RooAbsPdf *bkg_pdf = nullptr;
};

struct RooFitChannelConfig {
  TString name;
  TH1 *hist;
  std::vector<Double_t> events;
  Float_t fit_range_low;
  Float_t fit_range_high;
  Float_t display_bin_width_kev;
  Int_t num_peaks;
  std::vector<Double_t> mu_inits;
  std::vector<Bool_t> mu_fixed;
  std::vector<Bool_t> use_step_per_peak;
  Bool_t bkg_yield_fixed = kFALSE;
  Bool_t bkg_slope_fixed = kFALSE;
  Bool_t lock_shape_after_seed = kFALSE;
  Bool_t use_flat_background;
  Bool_t use_step;
  Bool_t use_low_exp_tail;
  Bool_t use_low_lin_tail;
  Bool_t use_high_exp_tail;
};

struct RooFitParamLink {
  TString target_channel;
  TString target_param;
  TString source_channel;
  TString source_param;
};

class RooFitUtils {
private:
  TH1 *working_hist_;
  std::vector<Double_t> events_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;
  Float_t display_bin_width_kev_;

  Bool_t use_flat_background_;
  Bool_t use_step_;
  Bool_t use_low_exp_tail_;
  Bool_t use_low_lin_tail_;
  Bool_t use_high_exp_tail_;

  Bool_t use_manual_init_;
  Bool_t interactive_;
  Bool_t fit_debug_;
  std::vector<Double_t> manual_params_;

  RooRealVar *x_;
  RooDataSet *unbinned_data_;
  RooAddPdf *total_pdf_;
  Int_t num_peaks_;

  std::vector<RooFitPeakModel> peaks_;
  RooFitBackgroundModel bkg_;
  std::vector<RooAbsArg *> owned_args_;

  std::vector<RooFitChannelConfig> sim_channels_;
  std::vector<RooFitParamLink> sim_links_;
  std::map<TString, FitResult> sim_seeds_;
  std::map<TString, std::vector<RooFitPeakModel>> sim_channel_peaks_;
  std::map<TString, RooFitBackgroundModel> sim_channel_bkg_;
  std::map<TString, RooAbsPdf *> sim_channel_pdfs_;
  std::map<TString, RooDataSet *> sim_channel_data_;
  std::map<TString, TString> sim_channel_range_names_;
  RooCategory *sim_category_;
  RooSimultaneous *sim_pdf_;
  RooDataSet *sim_combined_data_;
  Bool_t sim_mode_;

  static constexpr const char *kFitRangeName = "fitrange";

  void InitState();
  void BuildDisplayHistogram();
  void BuildUnbinnedData();
  static RooDataSet *BuildUnbinnedDataFrom(const std::vector<Double_t> &events,
                                           RooRealVar *x);
  RooRealVar *ResolveOrCreate(const TString &channel, const TString &param_name,
                              std::map<TString, RooRealVar *> &registry,
                              Double_t init_val, Double_t lo, Double_t hi);
  void BuildChannelModel(const RooFitChannelConfig &cfg,
                         std::map<TString, RooRealVar *> &registry);
  void ApplySeedToChannel(const TString &channel);
  void ApplyChannelMuLocks();
  void ApplyChannelBkgLocks();
  void ApplyChannelShapeLocks();
  void SaveSimInteractiveParams(const TString &input_name,
                                const TString &base_label);
  Bool_t LoadSimInteractiveParams(const TString &input_name,
                                  const TString &base_label);
  TString ParamFullName(const TString &channel, const TString &param);
  TString SourceForTarget(const TString &target);
  Double_t ComputeChannelChi2(const TString &channel,
                              const std::vector<RooFitPeakModel> &peaks,
                              const RooFitBackgroundModel &bkg, Int_t &ndof);
  void PlotChannel(const TString &channel, Int_t num_peaks,
                   const std::vector<RooFitPeakModel> &peaks,
                   const RooFitBackgroundModel &bkg, const TString &input_name,
                   const TString &base_label, const TString &chi2_label);
  PeakFitResult ExtractPeakResultFor(const RooFitPeakModel &p);

  Double_t EstimateBackground();

  void BuildPeak(Int_t peak_idx, Double_t mu_init, Double_t sigma_init,
                 Double_t peak_height, Double_t range_width);
  void BuildBackground(Double_t bkg_estimate, Double_t peak_height,
                       Double_t range_width);
  void BuildTotalModel();
  void ConfigureComponentFlagsForPeak(Int_t peak_idx);

  void FixComponent(Int_t peak_idx, const TString &component);
  void ReleaseComponent(Int_t peak_idx, const TString &component);
  Bool_t ComponentIsActive(Int_t peak_idx, const TString &component);

  std::vector<RooRealVar *> CollectFloatingParams();
  std::vector<RooRealVar *> CollectAllParams();
  RooFitResult *RunFit(Bool_t quiet);
  Double_t ComputeReducedChi2(RooFitResult *fit_result, Int_t &ndof);

  void SnapshotParams(std::vector<Double_t> &vals, std::vector<Double_t> &errs,
                      std::vector<Bool_t> &consts);
  void RestoreParams(const std::vector<Double_t> &vals,
                     const std::vector<Double_t> &errs,
                     const std::vector<Bool_t> &consts);
  void TestLowSideGroup(Int_t peak_idx, Double_t &best_chi2,
                        std::vector<Double_t> &best_vals,
                        std::vector<Double_t> &best_errs,
                        std::vector<Bool_t> &best_const);
  void TestHighTailIndependent(Int_t peak_idx, Double_t &best_chi2,
                               std::vector<Double_t> &best_vals,
                               std::vector<Double_t> &best_errs,
                               std::vector<Bool_t> &best_const);

  PeakFitResult ExtractPeakResult(Int_t peak_idx);

  void SaveInteractiveParams(const TString &input_name,
                             const TString &peak_name);
  Bool_t LoadInteractiveParams(const TString &input_name,
                               const TString &peak_name);
  void AdoptSavedRange(const TString &input_name, const TString &peak_name);
  void AdoptSavedSimRange(const TString &input_name, const TString &base_label);

  void SortPeaksByMu(Int_t num_peaks);
  void AppendPeakGraphs(std::vector<TGraph *> &components, Int_t peak_idx,
                        Style_t line_style, RooAbsPdf *background_pdf,
                        Double_t bkg_yield_val, Int_t npts, Double_t x_step,
                        Double_t bin_width);

  void RegisterOwned(RooAbsArg *arg);

public:
  RooFitUtils();
  RooFitUtils(const std::vector<Double_t> &events, Float_t fit_range_low,
              Float_t fit_range_high, Float_t display_bin_width_kev,
              Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
              Bool_t use_low_exp_tail = kFALSE,
              Bool_t use_low_lin_tail = kFALSE,
              Bool_t use_high_exp_tail = kFALSE);
  ~RooFitUtils();

  static std::vector<Double_t> LoadEventsFromTree(TTree *tree,
                                                  const TString &branch_name);
  static TH1F *BuildDisplayHistogramFrom(const std::vector<Double_t> &events,
                                         Float_t fit_range_low,
                                         Float_t fit_range_high,
                                         Float_t display_bin_width_kev);
  static void RefillDisplayHistogram(TH1 *hist,
                                     const std::vector<Double_t> &events,
                                     Float_t fit_range_low,
                                     Float_t fit_range_high,
                                     Float_t display_bin_width_kev);

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
  // When on, the simultaneous fit un-suppresses RooFit evaluation errors and
  // prints the seed NLL, so an invalid-NLL failure names the offending pdf.
  // Defaults to off; enable per-fit with sim.SetFitDebug() before
  // FitSimultaneous.
  void SetFitDebug(Bool_t fit_debug = kTRUE) { fit_debug_ = fit_debug; }

  void SetManualParameters(const std::vector<Double_t> &params);
  void SetManualParameter(Int_t index, Double_t value);
  void ClearManualParameters() {
    use_manual_init_ = kFALSE;
    manual_params_.clear();
  }

  void PlotFitSinglePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitDoublePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitTriplePeak(const TString input_name, const TString peak_name,
                         const TString label = "");

  FitResult FitSinglePeak(const TString input_name, const TString peak_name);
  // link_sigma: tie the two peaks to a single Gaussian width. Peak 2's PDFs are
  // rebuilt to reference peak 1's sigma and Sigma2 is frozen, so the doublet
  // shares one resolution (e.g. close, same-detector lines like the Pb K-beta
  // group) without Sigma2 floating and trading against the background.
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init,
                          Bool_t link_sigma = kFALSE);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          const PeakFitResult &constrained_peak,
                          Double_t mu2_init);
  FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                          const FitResult &constrained_peaks,
                          Double_t mu3_init);

  void AddChannel(
      const TString &name, const std::vector<Double_t> &events,
      Float_t fit_range_low, Float_t fit_range_high,
      Float_t display_bin_width_kev, Int_t num_peaks,
      const std::vector<Double_t> &mu_inits,
      Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
      Bool_t use_low_exp_tail = kFALSE, Bool_t use_low_lin_tail = kFALSE,
      Bool_t use_high_exp_tail = kFALSE,
      const std::vector<Bool_t> &mu_fixed = std::vector<Bool_t>(),
      Bool_t bkg_yield_fixed = kFALSE, Bool_t bkg_slope_fixed = kFALSE,
      Bool_t lock_shape_after_seed = kFALSE,
      const std::vector<Bool_t> &use_step_per_peak = std::vector<Bool_t>());
  void LinkParameter(const TString &target, const TString &source);
  void LinkPeakShape(const TString &target_channel, Int_t target_peak,
                     const TString &source_channel, Int_t source_peak);
  void SeedChannel(const TString &channel_name, const FitResult &result);
  std::vector<FitResult> FitSimultaneous(const TString &input_name,
                                         const TString &base_label);
};

#endif
