#include "RooFitUtils.hpp"

#include "InteractiveRooFitEditor.hpp"
#include "InteractiveSimultaneousFitEditor.hpp"

#include <cmath>

RooAbsPdf *RooFitFunctions::MakeGaussian(const TString &name, RooRealVar &x,
                                         RooRealVar &mu, RooRealVar &sigma) {
  return new RooGaussian(name.Data(), name.Data(), x, mu, sigma);
}

RooAbsPdf *RooFitFunctions::MakeStepShelf(const TString &name, RooRealVar &x,
                                          RooRealVar &mu, RooRealVar &sigma) {
  return new RooStepShelf(name.Data(), name.Data(), x, mu, sigma);
}

RooAbsPdf *RooFitFunctions::MakeLowExpTail(const TString &name, RooRealVar &x,
                                           RooRealVar &mu, RooRealVar &sigma,
                                           RooRealVar &tau) {
  return new RooLowExpTail(name.Data(), name.Data(), x, mu, sigma, tau);
}

RooAbsPdf *RooFitFunctions::MakeLowLinTail(const TString &name, RooRealVar &x,
                                           RooRealVar &mu, RooRealVar &sigma,
                                           RooRealVar &slope) {
  return new RooLowLinTail(name.Data(), name.Data(), x, mu, sigma, slope);
}

RooAbsPdf *RooFitFunctions::MakeHighExpTail(const TString &name, RooRealVar &x,
                                            RooRealVar &mu, RooRealVar &sigma,
                                            RooRealVar &tau) {
  return new RooHighExpTail(name.Data(), name.Data(), x, mu, sigma, tau);
}

RooAbsPdf *RooFitFunctions::MakeLinearBackground(const TString &name,
                                                 RooRealVar &x,
                                                 RooRealVar &slope) {
  return new RooPolynomial(name.Data(), name.Data(), x, RooArgList(slope));
}

void RooFitUtils::RegisterOwned(RooAbsArg *arg) { owned_args_.push_back(arg); }

void RooFitUtils::InitState() {
  working_hist_ = nullptr;
  events_.clear();
  fit_range_low_ = 0;
  fit_range_high_ = 0;
  display_bin_width_kev_ = 0;
  use_flat_background_ = kFALSE;
  use_step_ = kFALSE;
  use_low_exp_tail_ = kFALSE;
  use_low_lin_tail_ = kFALSE;
  use_high_exp_tail_ = kFALSE;
  use_manual_init_ = kFALSE;
  interactive_ = kFALSE;
  fit_debug_ = kFALSE;

  x_ = nullptr;
  unbinned_data_ = nullptr;
  total_pdf_ = nullptr;
  num_peaks_ = 0;

  sim_category_ = nullptr;
  sim_pdf_ = nullptr;
  sim_combined_data_ = nullptr;
  sim_mode_ = kFALSE;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  for (Int_t i = 0; i < RooMsgService::instance().numStreams(); i++) {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::InputArguments);
  }
  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::Ignore);
  RooRealVar::enableSilentClipping();
}

std::vector<Double_t>
RooFitUtils::LoadEventsFromTree(TTree *tree, const TString &branch_name) {
  std::vector<Double_t> out;
  if (!tree)
    return out;
  Float_t energy_f = 0;
  Double_t energy_d = 0;
  TBranch *br = tree->GetBranch(branch_name.Data());
  if (!br) {
    std::cerr << "ERROR: LoadEventsFromTree: branch '" << branch_name
              << "' not found" << std::endl;
    return out;
  }
  TLeaf *leaf = br->GetLeaf(branch_name.Data());
  Bool_t is_double = (leaf && TString(leaf->GetTypeName()) == "Double_t");
  if (is_double)
    tree->SetBranchAddress(branch_name.Data(), &energy_d);
  else
    tree->SetBranchAddress(branch_name.Data(), &energy_f);

  Long64_t n = tree->GetEntries();
  out.reserve(n);
  for (Long64_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    out.push_back(is_double ? energy_d : (Double_t)energy_f);
  }
  tree->ResetBranchAddresses();
  return out;
}

TH1F *RooFitUtils::BuildDisplayHistogramFrom(
    const std::vector<Double_t> &events, Float_t fit_range_low,
    Float_t fit_range_high, Float_t display_bin_width_kev) {
  Float_t hist_lo = 0.85f * fit_range_low;
  Float_t hist_hi = 1.15f * fit_range_high;
  Int_t nbins = TMath::Max(
      1, (Int_t)TMath::Nint((hist_hi - hist_lo) / display_bin_width_kev));
  hist_hi = hist_lo + nbins * display_bin_width_kev;
  TString hname = PlottingUtils::GetRandomName();
  TH1F *hist = new TH1F(hname,
                        TString::Format("; Energy [keV]; Counts / %.0f eV",
                                        display_bin_width_kev * 1000.0),
                        nbins, hist_lo, hist_hi);
  hist->SetDirectory(0);
  for (size_t i = 0; i < events.size(); i++) {
    Double_t e = events[i];
    if (e >= hist_lo && e < hist_hi)
      hist->Fill(e);
  }
  return hist;
}

void RooFitUtils::RefillDisplayHistogram(TH1 *hist,
                                         const std::vector<Double_t> &events,
                                         Float_t fit_range_low,
                                         Float_t fit_range_high,
                                         Float_t display_bin_width_kev) {
  Float_t hist_lo = 0.85f * fit_range_low;
  Float_t hist_hi = 1.15f * fit_range_high;
  Int_t nbins = TMath::Max(
      1, (Int_t)TMath::Nint((hist_hi - hist_lo) / display_bin_width_kev));
  hist_hi = hist_lo + nbins * display_bin_width_kev;
  hist->SetBins(nbins, hist_lo, hist_hi);
  hist->Reset();
  for (size_t i = 0; i < events.size(); i++) {
    Double_t e = events[i];
    if (e >= hist_lo && e < hist_hi)
      hist->Fill(e);
  }
}

RooDataSet *
RooFitUtils::BuildUnbinnedDataFrom(const std::vector<Double_t> &events,
                                   RooRealVar *x) {
  RooArgSet vars(*x);
  RooDataSet *ds = new RooDataSet("unbinned_data", "unbinned_data", vars);
  Double_t xmin = x->getMin();
  Double_t xmax = x->getMax();
  for (size_t i = 0; i < events.size(); i++) {
    Double_t e = events[i];
    if (e < xmin || e > xmax)
      continue;
    x->setVal(e);
    ds->add(vars);
  }
  return ds;
}

void RooFitUtils::BuildDisplayHistogram() {
  delete working_hist_;
  working_hist_ = BuildDisplayHistogramFrom(
      events_, fit_range_low_, fit_range_high_, display_bin_width_kev_);
}

void RooFitUtils::BuildUnbinnedData() {
  delete unbinned_data_;
  unbinned_data_ = BuildUnbinnedDataFrom(events_, x_);
}

RooFitUtils::RooFitUtils() {
  InitState();
  sim_mode_ = kTRUE;
}

RooFitUtils::RooFitUtils(const std::vector<Double_t> &events,
                         Float_t fit_range_low, Float_t fit_range_high,
                         Float_t display_bin_width_kev,
                         Bool_t use_flat_background, Bool_t use_step,
                         Bool_t use_low_exp_tail, Bool_t use_low_lin_tail,
                         Bool_t use_high_exp_tail) {
  InitState();
  events_ = events;
  fit_range_low_ = fit_range_low;
  fit_range_high_ = fit_range_high;
  display_bin_width_kev_ = display_bin_width_kev;
  use_flat_background_ = use_flat_background;
  use_step_ = use_step;
  use_low_exp_tail_ = use_low_exp_tail;
  use_low_lin_tail_ = use_low_lin_tail;
  use_high_exp_tail_ = use_high_exp_tail;

  BuildDisplayHistogram();

  std::cout << "Fit configuration:" << std::endl;
  std::cout << std::endl;
  if (use_flat_background_) {
    std::cout << "Background: FLAT" << std::endl;
  } else {
    std::cout << "Background: LINEAR" << std::endl;
  }
  std::cout << "Step function: " << (use_step_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "Low exponential tail: "
            << (use_low_exp_tail_ ? "ENABLED" : "DISABLED") << std::endl;
  std::cout << "Low linear tail: "
            << (use_low_lin_tail_ ? "ENABLED" : "DISABLED") << std::endl;
  std::cout << "High exponential tail: "
            << (use_high_exp_tail_ ? "ENABLED" : "DISABLED") << std::endl;
  std::cout << "Unbinned events: " << events_.size() << std::endl;
}

RooFitUtils::~RooFitUtils() {
  for (size_t i = 0; i < owned_args_.size(); i++) {
    delete owned_args_[i];
  }
  owned_args_.clear();
  delete unbinned_data_;
  delete sim_combined_data_;
  delete sim_pdf_;
  delete sim_category_;
  std::map<TString, RooDataSet *>::iterator dit;
  for (dit = sim_channel_data_.begin(); dit != sim_channel_data_.end(); ++dit) {
    delete dit->second;
  }
  sim_channel_data_.clear();
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    delete sim_channels_[i].hist;
  }
  delete working_hist_;
}

void RooFitUtils::SetManualParameters(const std::vector<Double_t> &params) {
  Int_t expected = num_peaks_ * 10 + 2;
  if (num_peaks_ == 0 || params.size() != (size_t)expected) {
    std::cerr << "ERROR: Manual parameters size (" << params.size()
              << ") doesn't match expected (" << expected << ")" << std::endl;
    return;
  }

  manual_params_ = params;
  use_manual_init_ = kTRUE;

  std::vector<RooRealVar *> all = CollectAllParams();
  for (size_t i = 0; i < all.size(); i++) {
    all[i]->setVal(params[i]);
    all[i]->setConstant(kTRUE);
  }

  std::cout << "Manual parameters set:" << std::endl;
  for (size_t i = 0; i < all.size(); i++) {
    std::cout << "  Par[" << i << "] " << all[i]->GetName() << " = "
              << params[i] << std::endl;
  }
}

void RooFitUtils::SetManualParameter(Int_t index, Double_t value) {
  std::vector<RooRealVar *> all = CollectAllParams();
  if (index < 0 || index >= (Int_t)all.size()) {
    std::cerr << "ERROR: Parameter index " << index << " out of range [0, "
              << (Int_t)all.size() - 1 << "]" << std::endl;
    return;
  }

  if (!use_manual_init_) {
    manual_params_.resize(all.size(), 0.0);
    use_manual_init_ = kTRUE;
  }

  manual_params_[index] = value;
  all[index]->setVal(value);
  all[index]->setConstant(kTRUE);

  std::cout << "Set Par[" << index << "] " << all[index]->GetName() << " = "
            << value << std::endl;
}

Double_t RooFitUtils::EstimateBackground() {
  Int_t left_bin = working_hist_->FindBin(fit_range_low_);
  Int_t right_bin = working_hist_->FindBin(fit_range_high_);

  Int_t n_sideband = (right_bin - left_bin) / 10;
  Double_t left_avg = 0;
  Double_t right_avg = 0;

  for (Int_t i = 0; i < n_sideband; i++) {
    left_avg += working_hist_->GetBinContent(left_bin + i);
    right_avg += working_hist_->GetBinContent(right_bin - i);
  }

  return (left_avg + right_avg) / (2.0 * n_sideband);
}

void RooFitUtils::BuildPeak(Int_t peak_idx, Double_t mu_init,
                            Double_t sigma_init, Double_t peak_height,
                            Double_t range_width) {
  RooFitPeakModel p;
  TString suffix = TString::Format("%d", peak_idx + 1);

  p.mu = new RooRealVar("Mu" + suffix, "Mu" + suffix, mu_init, fit_range_low_,
                        fit_range_high_);
  Double_t sigma_lo = working_hist_->GetBinWidth(1);
  Double_t sigma_hi = range_width * 0.5;
  if (sigma_init < sigma_lo)
    sigma_init = 2.0 * sigma_lo;
  p.sigma = new RooRealVar("Sigma" + suffix, "Sigma" + suffix, sigma_init,
                           sigma_lo, sigma_hi);

  Int_t mu_bin = working_hist_->FindBin(mu_init);
  Double_t local_height = working_hist_->GetBinContent(mu_bin);
  Double_t bkg_floor = EstimateBackground();
  Double_t net_height = local_height - bkg_floor;
  if (net_height < 0.1 * local_height)
    net_height = 0.1 * local_height;
  Double_t total_init =
      net_height * sigma_init * TMath::Sqrt(2.0 * TMath::Pi());
  p.gaus_yield =
      new RooRealVar("GausAmplitude" + suffix, "GausAmplitude" + suffix,
                     total_init, 0, peak_height * range_width * 10.0);

  p.ratio_step = new RooRealVar("StepAmplitude" + suffix,
                                "StepAmplitude" + suffix, 0.0, 0.0, 0.5);
  p.ratio_low_exp =
      new RooRealVar("LowExpTailAmplitude" + suffix,
                     "LowExpTailAmplitude" + suffix, 0.0, 0.0, 0.5);
  p.tau_low_exp = new RooRealVar("LowExpTailRatio" + suffix,
                                 "LowExpTailRatio" + suffix, 1.5, 1.0, 100.0);
  p.ratio_low_lin =
      new RooRealVar("LowLinTailAmplitude" + suffix,
                     "LowLinTailAmplitude" + suffix, 0.0, 0.0, 0.5);
  p.slope_low_lin = new RooRealVar("LowLinTailSlope" + suffix,
                                   "LowLinTailSlope" + suffix, 0.0, -0.1, 0.1);
  p.ratio_high_exp =
      new RooRealVar("HighExpTailAmplitude" + suffix,
                     "HighExpTailAmplitude" + suffix, 0.0, 0.0, 0.5);
  p.tau_high_exp = new RooRealVar("HighExpTailRatio" + suffix,
                                  "HighExpTailRatio" + suffix, 1.5, 1.0, 100.0);

  RegisterOwned(p.mu);
  RegisterOwned(p.sigma);
  RegisterOwned(p.gaus_yield);
  RegisterOwned(p.ratio_step);
  RegisterOwned(p.ratio_low_exp);
  RegisterOwned(p.tau_low_exp);
  RegisterOwned(p.ratio_low_lin);
  RegisterOwned(p.slope_low_lin);
  RegisterOwned(p.ratio_high_exp);
  RegisterOwned(p.tau_high_exp);

  p.gauss_pdf =
      RooFitFunctions::MakeGaussian("gauss_pdf" + suffix, *x_, *p.mu, *p.sigma);
  p.step_pdf =
      RooFitFunctions::MakeStepShelf("step_pdf" + suffix, *x_, *p.mu, *p.sigma);
  p.low_exp_pdf = RooFitFunctions::MakeLowExpTail(
      "low_exp_pdf" + suffix, *x_, *p.mu, *p.sigma, *p.tau_low_exp);
  p.low_lin_pdf = RooFitFunctions::MakeLowLinTail(
      "low_lin_pdf" + suffix, *x_, *p.mu, *p.sigma, *p.slope_low_lin);
  p.high_exp_pdf = RooFitFunctions::MakeHighExpTail(
      "high_exp_pdf" + suffix, *x_, *p.mu, *p.sigma, *p.tau_high_exp);

  RegisterOwned(p.gauss_pdf);
  RegisterOwned(p.step_pdf);
  RegisterOwned(p.low_exp_pdf);
  RegisterOwned(p.low_lin_pdf);
  RegisterOwned(p.high_exp_pdf);

  p.step_yield = new RooFormulaVar("step_yield" + suffix, "@0*@1",
                                   RooArgList(*p.gaus_yield, *p.ratio_step));
  p.low_exp_yield =
      new RooFormulaVar("low_exp_yield" + suffix, "@0*@1",
                        RooArgList(*p.gaus_yield, *p.ratio_low_exp));
  p.low_lin_yield =
      new RooFormulaVar("low_lin_yield" + suffix, "@0*@1",
                        RooArgList(*p.gaus_yield, *p.ratio_low_lin));
  p.high_exp_yield =
      new RooFormulaVar("high_exp_yield" + suffix, "@0*@1",
                        RooArgList(*p.gaus_yield, *p.ratio_high_exp));

  RegisterOwned(p.step_yield);
  RegisterOwned(p.low_exp_yield);
  RegisterOwned(p.low_lin_yield);
  RegisterOwned(p.high_exp_yield);

  peaks_.push_back(p);
}

void RooFitUtils::BuildBackground(Double_t bkg_estimate, Double_t peak_height,
                                  Double_t range_width) {
  bkg_.bkg_yield =
      new RooRealVar("BkgConstant", "BkgConstant", bkg_estimate * range_width,
                     0, peak_height * range_width * 10.0);
  // RooPolynomial evaluates as 1 + slope*x, so positivity over the fit range
  // forces slope >= -1/fit_range_high_. No constraint on the positive side, so
  // let it grow with the fit window instead of mirroring the tight neg bound.
  Double_t slope_lo = -0.9 / fit_range_high_;
  Double_t slope_hi = 5.0 / (fit_range_high_ - fit_range_low_);
  bkg_.bkg_slope =
      new RooRealVar("BkgSlope", "BkgSlope", 0.0, slope_lo, slope_hi);

  RegisterOwned(bkg_.bkg_yield);
  RegisterOwned(bkg_.bkg_slope);

  bkg_.bkg_pdf =
      RooFitFunctions::MakeLinearBackground("bkg_pdf", *x_, *bkg_.bkg_slope);
  if (use_flat_background_) {
    bkg_.bkg_slope->setVal(0.0);
    bkg_.bkg_slope->setConstant(kTRUE);
  }
  RegisterOwned(bkg_.bkg_pdf);
}

void RooFitUtils::BuildTotalModel() {
  RooArgList pdf_list;
  RooArgList coef_list;
  for (size_t pi = 0; pi < peaks_.size(); pi++) {
    pdf_list.add(*peaks_[pi].gauss_pdf);
    coef_list.add(*peaks_[pi].gaus_yield);
    pdf_list.add(*peaks_[pi].step_pdf);
    coef_list.add(*peaks_[pi].step_yield);
    pdf_list.add(*peaks_[pi].low_exp_pdf);
    coef_list.add(*peaks_[pi].low_exp_yield);
    pdf_list.add(*peaks_[pi].low_lin_pdf);
    coef_list.add(*peaks_[pi].low_lin_yield);
    pdf_list.add(*peaks_[pi].high_exp_pdf);
    coef_list.add(*peaks_[pi].high_exp_yield);
  }
  pdf_list.add(*bkg_.bkg_pdf);
  coef_list.add(*bkg_.bkg_yield);

  total_pdf_ = new RooAddPdf("total_pdf", "total_pdf", pdf_list, coef_list);
  RegisterOwned(total_pdf_);
}

void RooFitUtils::ConfigureComponentFlagsForPeak(Int_t peak_idx) {
  RooFitPeakModel &p = peaks_[peak_idx];

  if (use_step_) {
    p.ratio_step->setVal(0.0);
    p.ratio_step->setConstant(kFALSE);
  } else {
    p.ratio_step->setVal(0.0);
    p.ratio_step->setConstant(kTRUE);
  }

  if (use_low_exp_tail_) {
    p.ratio_low_exp->setVal(0.1);
    p.ratio_low_exp->setConstant(kFALSE);
    p.tau_low_exp->setVal(1.5);
    p.tau_low_exp->setConstant(kFALSE);
  } else {
    p.ratio_low_exp->setVal(0.0);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(1.0);
    p.tau_low_exp->setConstant(kTRUE);
  }

  if (use_low_lin_tail_) {
    p.ratio_low_lin->setVal(0.1);
    p.ratio_low_lin->setConstant(kFALSE);
    p.slope_low_lin->setVal(0.0);
    p.slope_low_lin->setConstant(kFALSE);
  } else {
    p.ratio_low_lin->setVal(0.0);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(0.0);
    p.slope_low_lin->setConstant(kTRUE);
  }

  if (use_high_exp_tail_) {
    p.ratio_high_exp->setVal(0.1);
    p.ratio_high_exp->setConstant(kFALSE);
    p.tau_high_exp->setVal(1.5);
    p.tau_high_exp->setConstant(kFALSE);
  } else {
    p.ratio_high_exp->setVal(0.0);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
}

void RooFitUtils::FixComponent(Int_t peak_idx, const TString &component) {
  RooFitPeakModel &p = peaks_[peak_idx];
  if (component == "step") {
    p.ratio_step->setVal(0.0);
    p.ratio_step->setConstant(kTRUE);
  } else if (component == "low_exp") {
    p.ratio_low_exp->setVal(0.0);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(1.0);
    p.tau_low_exp->setConstant(kTRUE);
  } else if (component == "low_lin") {
    p.ratio_low_lin->setVal(0.0);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(0.0);
    p.slope_low_lin->setConstant(kTRUE);
  } else if (component == "high_exp") {
    p.ratio_high_exp->setVal(0.0);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
}

void RooFitUtils::ReleaseComponent(Int_t peak_idx, const TString &component) {
  RooFitPeakModel &p = peaks_[peak_idx];
  if (component == "step") {
    p.ratio_step->setConstant(kFALSE);
    p.ratio_step->setVal(0.15);
  } else if (component == "low_exp") {
    p.ratio_low_exp->setConstant(kFALSE);
    p.tau_low_exp->setConstant(kFALSE);
    p.ratio_low_exp->setVal(0.15);
    p.tau_low_exp->setVal(1.5);
  } else if (component == "low_lin") {
    p.ratio_low_lin->setConstant(kFALSE);
    p.slope_low_lin->setConstant(kFALSE);
    p.ratio_low_lin->setVal(0.15);
    p.slope_low_lin->setVal(0.0);
  } else if (component == "high_exp") {
    p.ratio_high_exp->setConstant(kFALSE);
    p.tau_high_exp->setConstant(kFALSE);
    p.ratio_high_exp->setVal(0.15);
    p.tau_high_exp->setVal(1.5);
  }
}

std::vector<RooRealVar *> RooFitUtils::CollectAllParams() {
  std::vector<RooRealVar *> out;
  for (size_t pi = 0; pi < peaks_.size(); pi++) {
    out.push_back(peaks_[pi].mu);
    out.push_back(peaks_[pi].sigma);
    out.push_back(peaks_[pi].gaus_yield);
    out.push_back(peaks_[pi].ratio_step);
    out.push_back(peaks_[pi].ratio_low_exp);
    out.push_back(peaks_[pi].tau_low_exp);
    out.push_back(peaks_[pi].ratio_low_lin);
    out.push_back(peaks_[pi].slope_low_lin);
    out.push_back(peaks_[pi].ratio_high_exp);
    out.push_back(peaks_[pi].tau_high_exp);
  }
  out.push_back(bkg_.bkg_yield);
  out.push_back(bkg_.bkg_slope);
  return out;
}

std::vector<RooRealVar *> RooFitUtils::CollectFloatingParams() {
  std::vector<RooRealVar *> all = CollectAllParams();
  std::vector<RooRealVar *> out;
  for (size_t i = 0; i < all.size(); i++) {
    if (!all[i]->isConstant()) {
      out.push_back(all[i]);
    }
  }
  return out;
}

RooFitResult *RooFitUtils::RunFit(Bool_t quiet) {
  Int_t print_level = quiet ? -1 : 0;
  RooFitResult *result = total_pdf_->fitTo(
      *unbinned_data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE),
      RooFit::Range(kFitRangeName), RooFit::SumW2Error(kFALSE),
      RooFit::PrintLevel(print_level), RooFit::PrintEvalErrors(-1),
      RooFit::Strategy(2), RooFit::Minimizer("Minuit2", "migrad"),
      BestAvailableBackend());
  return result;
}

Double_t RooFitUtils::ComputeReducedChi2(RooFitResult *fit_result,
                                         Int_t &ndof) {
  RooArgSet nset(*x_);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  Double_t bin_width = working_hist_->GetBinWidth(1);
  Double_t saved = x_->getVal();

  Double_t chi2 = 0;
  Int_t nbins_in_range = 0;
  Int_t nbins_hist = working_hist_->GetNbinsX();
  for (Int_t i = 1; i <= nbins_hist; i++) {
    Double_t xv = working_hist_->GetBinCenter(i);
    if (xv < fit_range_low_ || xv > fit_range_high_)
      continue;
    Double_t data = working_hist_->GetBinContent(i);
    Double_t error = working_hist_->GetBinError(i);
    if (error <= 0 || data <= 0)
      continue;
    x_->setVal(xv);
    Double_t fit_val = total_exp * total_pdf_->getVal(&nset) * bin_width;
    Double_t residual = (data - fit_val) / error;
    chi2 += residual * residual;
    nbins_in_range++;
  }
  x_->setVal(saved);

  Int_t npars = fit_result ? fit_result->floatParsFinal().size()
                           : (Int_t)CollectFloatingParams().size();
  ndof = nbins_in_range - npars;
  if (ndof <= 0)
    return -1;
  return chi2 / ndof;
}

void RooFitUtils::SnapshotParams(std::vector<Double_t> &vals,
                                 std::vector<Double_t> &errs,
                                 std::vector<Bool_t> &consts) {
  std::vector<RooRealVar *> all = CollectAllParams();
  vals.resize(all.size());
  errs.resize(all.size());
  consts.resize(all.size());
  for (size_t i = 0; i < all.size(); i++) {
    vals[i] = all[i]->getVal();
    errs[i] = all[i]->getError();
    consts[i] = all[i]->isConstant();
  }
}

void RooFitUtils::RestoreParams(const std::vector<Double_t> &vals,
                                const std::vector<Double_t> &errs,
                                const std::vector<Bool_t> &consts) {
  std::vector<RooRealVar *> all = CollectAllParams();
  for (size_t i = 0; i < all.size() && i < vals.size(); i++) {
    all[i]->setVal(vals[i]);
    all[i]->setError(errs[i]);
    all[i]->setConstant(consts[i]);
  }
}

void RooFitUtils::TestLowSideGroup(Int_t peak_idx, Double_t &best_chi2,
                                   std::vector<Double_t> &best_vals,
                                   std::vector<Double_t> &best_errs,
                                   std::vector<Bool_t> &best_const) {
  Bool_t any_low_side = use_step_ || use_low_exp_tail_ || use_low_lin_tail_;
  if (!any_low_side)
    return;

  std::cout << "Testing low-side component group for peak " << peak_idx + 1
            << "..." << std::endl;
  if (use_step_)
    ReleaseComponent(peak_idx, "step");
  if (use_low_exp_tail_)
    ReleaseComponent(peak_idx, "low_exp");
  if (use_low_lin_tail_)
    ReleaseComponent(peak_idx, "low_lin");

  RooFitResult *group_fit = RunFit(kTRUE);
  Bool_t group_ok = group_fit && group_fit->status() == 0;
  Int_t tmp_ndof = 0;
  Double_t chi2_group = group_ok ? ComputeReducedChi2(group_fit, tmp_ndof) : -1;
  delete group_fit;

  if (group_ok && chi2_group < best_chi2) {
    std::cout << "Low-side group peak " << peak_idx + 1
              << " ACCEPTED, pruning..." << std::endl;
    best_chi2 = chi2_group;
    SnapshotParams(best_vals, best_errs, best_const);

    const TString comps[3] = {"step", "low_exp", "low_lin"};
    Bool_t enabled[3] = {use_step_, use_low_exp_tail_, use_low_lin_tail_};
    for (Int_t ci = 0; ci < 3; ci++) {
      if (!enabled[ci])
        continue;
      FixComponent(peak_idx, comps[ci]);
      RooFitResult *pf = RunFit(kTRUE);
      Bool_t ok = pf && pf->status() == 0;
      Int_t nd = 0;
      Double_t c2 = ok ? ComputeReducedChi2(pf, nd) : -1;
      delete pf;
      if (ok && c2 <= best_chi2) {
        std::cout << "  " << comps[ci] << " peak " << peak_idx + 1 << " pruned"
                  << std::endl;
        best_chi2 = c2;
        SnapshotParams(best_vals, best_errs, best_const);
      } else {
        std::cout << "  " << comps[ci] << " peak " << peak_idx + 1
                  << " retained" << std::endl;
        ReleaseComponent(peak_idx, comps[ci]);
        RestoreParams(best_vals, best_errs, best_const);
      }
    }
  } else {
    std::cout << "Low-side group peak " << peak_idx + 1 << " REJECTED"
              << std::endl;
    if (use_step_)
      FixComponent(peak_idx, "step");
    if (use_low_exp_tail_)
      FixComponent(peak_idx, "low_exp");
    if (use_low_lin_tail_)
      FixComponent(peak_idx, "low_lin");
    RestoreParams(best_vals, best_errs, best_const);
  }
}

void RooFitUtils::TestHighTailIndependent(Int_t peak_idx, Double_t &best_chi2,
                                          std::vector<Double_t> &best_vals,
                                          std::vector<Double_t> &best_errs,
                                          std::vector<Bool_t> &best_const) {
  if (!use_high_exp_tail_)
    return;

  std::cout << "Testing high exponential tail for peak " << peak_idx + 1
            << "..." << std::endl;
  ReleaseComponent(peak_idx, "high_exp");
  RooFitResult *htail_fit = RunFit(kTRUE);
  Bool_t htail_ok = htail_fit && htail_fit->status() == 0;
  Int_t tmp_ndof = 0;
  Double_t chi2_htail = htail_ok ? ComputeReducedChi2(htail_fit, tmp_ndof) : -1;
  delete htail_fit;
  if (htail_ok && chi2_htail < best_chi2) {
    std::cout << "High exp tail peak " << peak_idx + 1 << " ACCEPTED"
              << std::endl;
    best_chi2 = chi2_htail;
    SnapshotParams(best_vals, best_errs, best_const);
  } else {
    std::cout << "High exp tail peak " << peak_idx + 1 << " REJECTED"
              << std::endl;
    FixComponent(peak_idx, "high_exp");
    RestoreParams(best_vals, best_errs, best_const);
  }
}

PeakFitResult RooFitUtils::ExtractPeakResult(Int_t peak_idx) {
  RooFitPeakModel &p = peaks_[peak_idx];
  PeakFitResult result;
  result.mu = p.mu->getVal();
  result.mu_error = p.mu->getError();
  result.sigma = p.sigma->getVal();
  result.sigma_error = p.sigma->getError();
  result.gaus_amplitude = p.gaus_yield->getVal();
  result.gaus_amplitude_error = p.gaus_yield->getError();

  Double_t ga = result.gaus_amplitude;
  result.step_amplitude = p.ratio_step->getVal() * ga;
  result.step_amplitude_error = p.ratio_step->getError() * ga;
  result.low_exp_tail_amplitude = p.ratio_low_exp->getVal() * ga;
  result.low_exp_tail_amplitude_error = p.ratio_low_exp->getError() * ga;
  result.low_exp_tail_ratio = p.tau_low_exp->getVal();
  result.low_exp_tail_ratio_error = p.tau_low_exp->getError();
  result.low_lin_tail_amplitude = p.ratio_low_lin->getVal() * ga;
  result.low_lin_tail_amplitude_error = p.ratio_low_lin->getError() * ga;
  result.low_lin_tail_slope = p.slope_low_lin->getVal();
  result.low_lin_tail_slope_error = p.slope_low_lin->getError();
  result.high_exp_tail_amplitude = p.ratio_high_exp->getVal() * ga;
  result.high_exp_tail_amplitude_error = p.ratio_high_exp->getError() * ga;
  result.high_exp_tail_ratio = p.tau_high_exp->getVal();
  result.high_exp_tail_ratio_error = p.tau_high_exp->getError();
  return result;
}

void RooFitUtils::SortPeaksByMu(Int_t num_peaks) {
  for (Int_t i = 0; i < num_peaks - 1; i++) {
    for (Int_t j = 0; j < num_peaks - i - 1; j++) {
      Double_t mu_j = peaks_[j].mu->getVal();
      Double_t mu_next = peaks_[j + 1].mu->getVal();
      if (mu_j > mu_next) {
        std::cout << "Sorting peaks: swapping peak " << j + 1 << " (mu=" << mu_j
                  << ") and peak " << j + 2 << " (mu=" << mu_next << ")"
                  << std::endl;
        RooFitPeakModel tmp = peaks_[j];
        peaks_[j] = peaks_[j + 1];
        peaks_[j + 1] = tmp;
      }
    }
  }
}

void RooFitUtils::AdoptSavedSimRange(const TString &input_name,
                                     const TString &base_label) {
  if (!interactive_)
    return;
  TString filename = PlottingUtils::GetPlotsBaseDir() + "/fits/" + base_label +
                     "_" + input_name + ".simroofits";
  std::ifstream in(filename.Data());
  if (!in.is_open())
    return;
  std::string token;
  if (!(in >> token) || token != "RANGE")
    return;
  Double_t rlo, rhi;
  if (!(in >> rlo >> rhi))
    return;
  if (!(rlo < rhi))
    return;
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    sim_channels_[i].fit_range_low = rlo;
    sim_channels_[i].fit_range_high = rhi;
    delete sim_channels_[i].hist;
    sim_channels_[i].hist =
        BuildDisplayHistogramFrom(sim_channels_[i].events, rlo, rhi,
                                  sim_channels_[i].display_bin_width_kev);
  }
}

void RooFitUtils::AdoptSavedRange(const TString &input_name,
                                  const TString &peak_name) {
  if (!interactive_)
    return;
  TString filename = PlottingUtils::GetPlotsBaseDir() + "/fits/" + peak_name +
                     "_" + input_name + ".roofits";
  std::ifstream in(filename.Data());
  if (!in.is_open())
    return;
  std::string token;
  if (!(in >> token) || token != "RANGE")
    return;
  Double_t rlo, rhi;
  if (!(in >> rlo >> rhi))
    return;
  if (!(rlo < rhi))
    return;
  fit_range_low_ = rlo;
  fit_range_high_ = rhi;
  BuildDisplayHistogram();
}

void RooFitUtils::SaveInteractiveParams(const TString &input_name,
                                        const TString &peak_name) {
  TString fits_dir = PlottingUtils::GetPlotsBaseDir() + "/fits";
  gSystem->mkdir(fits_dir, kTRUE);
  TString filename = fits_dir + "/" + peak_name + "_" + input_name + ".roofits";
  std::ofstream out(filename.Data());
  if (!out.is_open()) {
    std::cerr << "WARNING: Could not save interactive params to " << filename
              << std::endl;
    return;
  }
  out << std::setprecision(17);
  out << "RANGE " << fit_range_low_ << " " << fit_range_high_ << "\n";
  std::vector<RooRealVar *> all = CollectAllParams();
  for (size_t i = 0; i < all.size(); i++) {
    out << all[i]->GetName() << " " << all[i]->getVal() << " "
        << all[i]->getError() << " " << (all[i]->isConstant() ? 1 : 0) << "\n";
  }
  out.close();
  std::cout << "Saved interactive params to " << filename << std::endl;
}

Bool_t RooFitUtils::LoadInteractiveParams(const TString &input_name,
                                          const TString &peak_name) {
  TString filename = PlottingUtils::GetPlotsBaseDir() + "/fits/" + peak_name +
                     "_" + input_name + ".roofits";
  std::ifstream in(filename.Data());
  if (!in.is_open())
    return kFALSE;

  std::vector<RooRealVar *> all = CollectAllParams();
  std::string token;
  Int_t idx = 0;

  in >> token;
  if (token == "RANGE") {
    Double_t rlo, rhi;
    in >> rlo >> rhi;
    fit_range_low_ = rlo;
    fit_range_high_ = rhi;
    x_->setRange(rlo, rhi);
    x_->setRange(kFitRangeName, rlo, rhi);
    BuildDisplayHistogram();
  }

  Double_t value, error;
  Int_t fixed;
  while (in >> token >> value >> error >> fixed && idx < (Int_t)all.size()) {
    all[idx]->setVal(value);
    all[idx]->setError(error);
    all[idx]->setConstant(fixed ? kTRUE : kFALSE);
    idx++;
  }
  in.close();

  if (idx != (Int_t)all.size()) {
    std::cerr << "WARNING: Parameter count mismatch in " << filename
              << " (expected " << all.size() << ", got " << idx << ")"
              << std::endl;
    return kFALSE;
  }

  std::cout << "Loaded interactive params from " << filename << std::endl;
  return kTRUE;
}

void RooFitUtils::SaveSimInteractiveParams(const TString &input_name,
                                           const TString &base_label) {
  TString fits_dir = PlottingUtils::GetPlotsBaseDir() + "/fits";
  gSystem->mkdir(fits_dir, kTRUE);
  TString filename =
      fits_dir + "/" + base_label + "_" + input_name + ".simroofits";
  std::ofstream out(filename.Data());
  if (!out.is_open()) {
    std::cerr << "WARNING: Could not save sim interactive params to "
              << filename << std::endl;
    return;
  }
  out << std::setprecision(17);
  out << "RANGE " << x_->getMin("fitrange") << " " << x_->getMax("fitrange")
      << "\n";

  std::set<RooRealVar *> seen;
  std::vector<RooRealVar *> ordered;
  for (size_t ci = 0; ci < sim_channels_.size(); ci++) {
    const TString &cname = sim_channels_[ci].name;
    std::vector<RooFitPeakModel> &peaks = sim_channel_peaks_[cname];
    for (size_t pi = 0; pi < peaks.size(); pi++) {
      RooFitPeakModel &p = peaks[pi];
      RooRealVar *vars[10] = {p.mu,
                              p.sigma,
                              p.gaus_yield,
                              p.ratio_step,
                              p.ratio_low_exp,
                              p.tau_low_exp,
                              p.ratio_low_lin,
                              p.slope_low_lin,
                              p.ratio_high_exp,
                              p.tau_high_exp};
      for (Int_t k = 0; k < 10; k++) {
        if (seen.insert(vars[k]).second)
          ordered.push_back(vars[k]);
      }
    }
    RooFitBackgroundModel &bkg = sim_channel_bkg_[cname];
    if (seen.insert(bkg.bkg_yield).second)
      ordered.push_back(bkg.bkg_yield);
    if (seen.insert(bkg.bkg_slope).second)
      ordered.push_back(bkg.bkg_slope);
  }

  for (size_t i = 0; i < ordered.size(); i++) {
    out << ordered[i]->GetName() << " " << ordered[i]->getVal() << " "
        << ordered[i]->getError() << " " << (ordered[i]->isConstant() ? 1 : 0)
        << "\n";
  }
  out.close();
  std::cout << "Saved sim interactive params to " << filename << std::endl;
}

Bool_t RooFitUtils::LoadSimInteractiveParams(const TString &input_name,
                                             const TString &base_label) {
  TString filename = PlottingUtils::GetPlotsBaseDir() + "/fits/" + base_label +
                     "_" + input_name + ".simroofits";
  std::ifstream in(filename.Data());
  if (!in.is_open())
    return kFALSE;

  std::map<std::string, RooRealVar *> by_name;
  for (size_t ci = 0; ci < sim_channels_.size(); ci++) {
    const TString &cname = sim_channels_[ci].name;
    std::vector<RooFitPeakModel> &peaks = sim_channel_peaks_[cname];
    for (size_t pi = 0; pi < peaks.size(); pi++) {
      RooFitPeakModel &p = peaks[pi];
      RooRealVar *vars[10] = {p.mu,
                              p.sigma,
                              p.gaus_yield,
                              p.ratio_step,
                              p.ratio_low_exp,
                              p.tau_low_exp,
                              p.ratio_low_lin,
                              p.slope_low_lin,
                              p.ratio_high_exp,
                              p.tau_high_exp};
      for (Int_t k = 0; k < 10; k++) {
        by_name[vars[k]->GetName()] = vars[k];
      }
    }
    RooFitBackgroundModel &bkg = sim_channel_bkg_[cname];
    by_name[bkg.bkg_yield->GetName()] = bkg.bkg_yield;
    by_name[bkg.bkg_slope->GetName()] = bkg.bkg_slope;
  }

  std::string token;
  in >> token;
  if (token == "RANGE") {
    Double_t rlo, rhi;
    in >> rlo >> rhi;
    x_->setRange(rlo, rhi);
    x_->setRange("fitrange", rlo, rhi);
  }

  Double_t value, error;
  Int_t fixed;
  while (in >> token >> value >> error >> fixed) {
    std::map<std::string, RooRealVar *>::iterator it = by_name.find(token);
    if (it == by_name.end()) {
      std::cerr << "WARNING: param " << token
                << " from .simroofits not found in model" << std::endl;
      continue;
    }
    it->second->setVal(value);
    it->second->setError(error);
    it->second->setConstant(fixed ? kTRUE : kFALSE);
  }
  in.close();
  std::cout << "Loaded sim interactive params from " << filename << std::endl;
  return kTRUE;
}

void RooFitUtils::AppendPeakGraphs(std::vector<TGraph *> &components,
                                   Int_t peak_idx, Style_t line_style,
                                   RooAbsPdf *background_pdf,
                                   Double_t bkg_yield_val, Int_t npts,
                                   Double_t x_step, Double_t bin_width) {
  RooFitPeakModel &p = peaks_[peak_idx];
  Width_t line_width = PlottingUtils::GetLineWidth();
  RooArgSet nset(*x_);

  TGraph *peak_graph = new TGraph(npts);
  Double_t gy = p.gaus_yield->getVal();
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = gy * p.gauss_pdf->getVal(&nset) * bin_width;
    Double_t bkg_v = bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
    peak_graph->SetPoint(i, xv, y + bkg_v);
  }
  peak_graph->SetLineColor(kBlack);
  peak_graph->SetLineStyle(line_style);
  peak_graph->SetLineWidth(line_width);
  components.push_back(peak_graph);

  if (TMath::Abs(p.ratio_step->getVal()) > 1e-6) {
    TGraph *step_graph = new TGraph(npts);
    Double_t sy = p.step_yield->getVal();
    for (Int_t i = 0; i < npts; i++) {
      Double_t xv = fit_range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t y = sy * p.step_pdf->getVal(&nset) * bin_width;
      Double_t bkg_v =
          bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
      step_graph->SetPoint(i, xv, y + bkg_v);
    }
    step_graph->SetLineColor(kGray);
    step_graph->SetLineStyle(line_style);
    step_graph->SetLineWidth(line_width);
    components.push_back(step_graph);
  }

  if (TMath::Abs(p.ratio_low_exp->getVal()) > 1e-6 ||
      TMath::Abs(p.ratio_low_lin->getVal()) > 1e-6) {
    TGraph *low_tail_graph = new TGraph(npts);
    Double_t lexp_y = p.low_exp_yield->getVal();
    Double_t llin_y = p.low_lin_yield->getVal();
    for (Int_t i = 0; i < npts; i++) {
      Double_t xv = fit_range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t y_exp = lexp_y * p.low_exp_pdf->getVal(&nset) * bin_width;
      Double_t y_lin = llin_y * p.low_lin_pdf->getVal(&nset) * bin_width;
      Double_t bkg_v =
          bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
      low_tail_graph->SetPoint(i, xv, y_exp + y_lin + bkg_v);
    }
    low_tail_graph->SetLineColor(kRed);
    low_tail_graph->SetLineStyle(line_style);
    low_tail_graph->SetLineWidth(line_width);
    components.push_back(low_tail_graph);
  }

  if (TMath::Abs(p.ratio_high_exp->getVal()) > 1e-6) {
    TGraph *high_tail_graph = new TGraph(npts);
    Double_t hexp_y = p.high_exp_yield->getVal();
    for (Int_t i = 0; i < npts; i++) {
      Double_t xv = fit_range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t y = hexp_y * p.high_exp_pdf->getVal(&nset) * bin_width;
      Double_t bkg_v =
          bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
      high_tail_graph->SetPoint(i, xv, y + bkg_v);
    }
    high_tail_graph->SetLineColor(kOrange);
    high_tail_graph->SetLineStyle(line_style);
    high_tail_graph->SetLineWidth(line_width);
    components.push_back(high_tail_graph);
  }
}

void RooFitUtils::PlotFitSinglePeak(const TString input_name,
                                    const TString peak_name,
                                    const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();
  Double_t bin_width = working_hist_->GetBinWidth(1);
  RooArgSet nset(*x_);

  TGraph *total_graph = new TGraph(npts);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = total_exp * total_pdf_->getVal(&nset) * bin_width;
    total_graph->SetPoint(i, xv, y);
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  Double_t bkg_yield_val = bkg_.bkg_yield->getVal();
  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = bkg_yield_val * bkg_.bkg_pdf->getVal(&nset) * bin_width;
    background_graph->SetPoint(i, xv, y);
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);
  AppendPeakGraphs(components, 0, 1, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                   bin_width);

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);

  delete total_graph;
  for (Int_t i = 0; i < (Int_t)components.size(); i++) {
    delete components[i];
  }
}

void RooFitUtils::PlotFitDoublePeak(const TString input_name,
                                    const TString peak_name,
                                    const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();
  Double_t bin_width = working_hist_->GetBinWidth(1);
  RooArgSet nset(*x_);

  TGraph *total_graph = new TGraph(npts);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = total_exp * total_pdf_->getVal(&nset) * bin_width;
    total_graph->SetPoint(i, xv, y);
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  Double_t bkg_yield_val = bkg_.bkg_yield->getVal();
  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = bkg_yield_val * bkg_.bkg_pdf->getVal(&nset) * bin_width;
    background_graph->SetPoint(i, xv, y);
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);
  AppendPeakGraphs(components, 0, 1, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                   bin_width);
  AppendPeakGraphs(components, 1, 3, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                   bin_width);

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);

  delete total_graph;
  for (Int_t i = 0; i < (Int_t)components.size(); i++) {
    delete components[i];
  }
}

void RooFitUtils::PlotFitTriplePeak(const TString input_name,
                                    const TString peak_name,
                                    const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();
  Double_t bin_width = working_hist_->GetBinWidth(1);
  RooArgSet nset(*x_);

  TGraph *total_graph = new TGraph(npts);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = total_exp * total_pdf_->getVal(&nset) * bin_width;
    total_graph->SetPoint(i, xv, y);
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  Double_t bkg_yield_val = bkg_.bkg_yield->getVal();
  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = bkg_yield_val * bkg_.bkg_pdf->getVal(&nset) * bin_width;
    background_graph->SetPoint(i, xv, y);
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);
  AppendPeakGraphs(components, 0, 1, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                   bin_width);
  AppendPeakGraphs(components, 1, 3, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                   bin_width);
  AppendPeakGraphs(components, 2, 4, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                   bin_width);

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);

  delete total_graph;
  for (Int_t i = 0; i < (Int_t)components.size(); i++) {
    delete components[i];
  }
}

FitResult RooFitUtils::FitSinglePeak(const TString input_name,
                                     const TString peak_name) {
  FitResult results;
  results.peaks.emplace_back();

  AdoptSavedRange(input_name, peak_name);

  num_peaks_ = 1;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t mu_init = (fit_range_low_ + fit_range_high_) / 2;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  BuildPeak(0, mu_init, sigma_init, peak_height, range_width);
  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  BuildUnbinnedData();
  x_->setRange(fit_range_low_, fit_range_high_);

  ConfigureComponentFlagsForPeak(0);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(
              working_hist_, &events_, display_bin_width_kev_, total_pdf_, x_,
              unbinned_data_, &peaks_, &bkg_, fit_range_low_, fit_range_high_,
              peak_name + " / " + input_name)) {
        fit_range_low_ = x_->getMin("fitrange");
        fit_range_high_ = x_->getMax("fitrange");
        BuildDisplayHistogram();
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(0, "step");
    FixComponent(0, "low_exp");
    FixComponent(0, "low_lin");
    FixComponent(0, "high_exp");

    if (use_manual_init_) {
      std::cout << "Using manually initialized parameters" << std::endl;
      std::vector<RooRealVar *> all = CollectAllParams();
      for (size_t i = 0; i < manual_params_.size() && i < all.size(); i++) {
        all[i]->setVal(manual_params_[i]);
      }
    }

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial fit failed" << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(0, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(0, best_chi2, best_vals, best_errs, best_const);

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitSinglePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Fit did not converge" << std::endl;
  }

  return results;
}

FitResult RooFitUtils::FitDoublePeak(const TString input_name,
                                     const TString peak_name, Double_t mu1_init,
                                     Double_t mu2_init, Bool_t link_sigma) {
  FitResult results;
  results.peaks.emplace_back();
  results.peaks.emplace_back();

  AdoptSavedRange(input_name, peak_name);

  if (mu1_init > mu2_init) {
    std::cout << "WARNING: mu1_init > mu2_init, swapping initial values"
              << std::endl;
    Double_t tmp = mu1_init;
    mu1_init = mu2_init;
    mu2_init = tmp;
  }

  num_peaks_ = 2;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  BuildPeak(0, mu1_init, sigma_init, peak_height, range_width);
  BuildPeak(1, mu2_init, sigma_init, peak_height, range_width);

  // Share one Gaussian width across the doublet: rebuild peak 2's PDFs to use
  // peak 1's sigma var, then freeze Sigma2 so it neither floats nor is saved as
  // an independent width. Sigma2 stays in the peak struct (param count, save/
  // load, and CollectAllParams contracts unchanged) but is inert.
  if (link_sigma) {
    RooRealVar *s = peaks_[0].sigma;
    RooFitPeakModel &p = peaks_[1];
    p.sigma->setVal(s->getVal());
    p.sigma->setConstant(kTRUE);
    p.gauss_pdf =
        RooFitFunctions::MakeGaussian("gauss_pdf2_linked", *x_, *p.mu, *s);
    p.step_pdf =
        RooFitFunctions::MakeStepShelf("step_pdf2_linked", *x_, *p.mu, *s);
    p.low_exp_pdf = RooFitFunctions::MakeLowExpTail("low_exp_pdf2_linked", *x_,
                                                    *p.mu, *s, *p.tau_low_exp);
    p.low_lin_pdf = RooFitFunctions::MakeLowLinTail(
        "low_lin_pdf2_linked", *x_, *p.mu, *s, *p.slope_low_lin);
    p.high_exp_pdf = RooFitFunctions::MakeHighExpTail(
        "high_exp_pdf2_linked", *x_, *p.mu, *s, *p.tau_high_exp);
    RegisterOwned(p.gauss_pdf);
    RegisterOwned(p.step_pdf);
    RegisterOwned(p.low_exp_pdf);
    RegisterOwned(p.low_lin_pdf);
    RegisterOwned(p.high_exp_pdf);
  }

  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  Double_t mu_midpoint = 0.5 * (mu1_init + mu2_init);
  peaks_[0].mu->setMax(mu_midpoint);
  peaks_[1].mu->setMin(mu_midpoint);

  BuildUnbinnedData();
  x_->setRange(fit_range_low_, fit_range_high_);

  ConfigureComponentFlagsForPeak(0);
  ConfigureComponentFlagsForPeak(1);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(
              working_hist_, &events_, display_bin_width_kev_, total_pdf_, x_,
              unbinned_data_, &peaks_, &bkg_, fit_range_low_, fit_range_high_,
              peak_name + " / " + input_name)) {
        fit_range_low_ = x_->getMin("fitrange");
        fit_range_high_ = x_->getMax("fitrange");
        BuildDisplayHistogram();
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(0, "step");
    FixComponent(0, "low_exp");
    FixComponent(0, "low_lin");
    FixComponent(0, "high_exp");
    FixComponent(1, "step");
    FixComponent(1, "low_exp");
    FixComponent(1, "low_lin");
    FixComponent(1, "high_exp");

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial double peak fit failed" << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(0, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(1, best_chi2, best_vals, best_errs, best_const);

    {
      std::cout
          << "Testing inter-peak group (peak1 high tail + peak2 low-side)..."
          << std::endl;
      if (use_high_exp_tail_)
        ReleaseComponent(0, "high_exp");
      if (use_step_)
        ReleaseComponent(1, "step");
      if (use_low_exp_tail_)
        ReleaseComponent(1, "low_exp");
      if (use_low_lin_tail_)
        ReleaseComponent(1, "low_lin");

      RooFitResult *group_fit = RunFit(kTRUE);
      Bool_t ok = group_fit && group_fit->status() == 0;
      Int_t nd = 0;
      Double_t c2 = ok ? ComputeReducedChi2(group_fit, nd) : -1;
      delete group_fit;

      if (ok && c2 < best_chi2) {
        std::cout << "Inter-peak group ACCEPTED, pruning..." << std::endl;
        best_chi2 = c2;
        SnapshotParams(best_vals, best_errs, best_const);

        struct CompRef {
          Int_t peak_idx;
          TString comp;
          Bool_t enabled;
        };
        CompRef refs[4] = {
            {0, "high_exp", use_high_exp_tail_},
            {1, "step", use_step_},
            {1, "low_exp", use_low_exp_tail_},
            {1, "low_lin", use_low_lin_tail_},
        };
        for (Int_t ri = 0; ri < 4; ri++) {
          if (!refs[ri].enabled)
            continue;
          FixComponent(refs[ri].peak_idx, refs[ri].comp);
          RooFitResult *pf = RunFit(kTRUE);
          Bool_t pok = pf && pf->status() == 0;
          Int_t pnd = 0;
          Double_t pc2 = pok ? ComputeReducedChi2(pf, pnd) : -1;
          delete pf;
          if (pok && pc2 <= best_chi2) {
            std::cout << "  " << refs[ri].comp << " peak "
                      << refs[ri].peak_idx + 1 << " pruned" << std::endl;
            best_chi2 = pc2;
            SnapshotParams(best_vals, best_errs, best_const);
          } else {
            std::cout << "  " << refs[ri].comp << " peak "
                      << refs[ri].peak_idx + 1 << " retained" << std::endl;
            ReleaseComponent(refs[ri].peak_idx, refs[ri].comp);
            RestoreParams(best_vals, best_errs, best_const);
          }
        }
      } else {
        std::cout << "Inter-peak group REJECTED" << std::endl;
        if (use_high_exp_tail_)
          FixComponent(0, "high_exp");
        if (use_step_)
          FixComponent(1, "step");
        if (use_low_exp_tail_)
          FixComponent(1, "low_exp");
        if (use_low_lin_tail_)
          FixComponent(1, "low_lin");
        RestoreParams(best_vals, best_errs, best_const);
      }
    }

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Double peak fit converged successfully" << std::endl;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    } else {
      std::cout << "ERROR: Double peak fit failed to converge" << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    SortPeaksByMu(2);
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitDoublePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.peaks[1] = ExtractPeakResult(1);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Double peak fit failed" << std::endl;
  }

  return results;
}

FitResult RooFitUtils::FitDoublePeak(const TString input_name,
                                     const TString peak_name,
                                     const PeakFitResult &constrained_peak,
                                     Double_t mu2_init) {
  FitResult results;
  results.peaks.emplace_back();
  results.peaks.emplace_back();

  AdoptSavedRange(input_name, peak_name);

  num_peaks_ = 2;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  BuildPeak(0, constrained_peak.mu, constrained_peak.sigma, peak_height,
            range_width);
  BuildPeak(1, mu2_init, sigma_init, peak_height, range_width);
  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  BuildUnbinnedData();
  x_->setRange(fit_range_low_, fit_range_high_);

  {
    RooFitPeakModel &p = peaks_[0];
    Double_t cga = constrained_peak.gaus_amplitude;
    p.mu->setVal(constrained_peak.mu);
    p.mu->setConstant(kTRUE);
    p.sigma->setVal(constrained_peak.sigma);
    p.sigma->setConstant(kTRUE);
    p.gaus_yield->setVal(cga);
    p.gaus_yield->setRange(0, peak_height * range_width * 10.0);
    p.gaus_yield->setConstant(kFALSE);
    p.ratio_step->setVal(constrained_peak.step_amplitude / cga);
    p.ratio_step->setConstant(kTRUE);
    p.ratio_low_exp->setVal(constrained_peak.low_exp_tail_amplitude / cga);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(constrained_peak.low_exp_tail_ratio > 0
                              ? constrained_peak.low_exp_tail_ratio
                              : 1.0);
    p.tau_low_exp->setConstant(kTRUE);
    p.ratio_low_lin->setVal(constrained_peak.low_lin_tail_amplitude / cga);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(constrained_peak.low_lin_tail_slope);
    p.slope_low_lin->setConstant(kTRUE);
    p.ratio_high_exp->setVal(constrained_peak.high_exp_tail_amplitude / cga);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(constrained_peak.high_exp_tail_ratio > 0
                               ? constrained_peak.high_exp_tail_ratio
                               : 1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
  ConfigureComponentFlagsForPeak(1);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(
              working_hist_, &events_, display_bin_width_kev_, total_pdf_, x_,
              unbinned_data_, &peaks_, &bkg_, fit_range_low_, fit_range_high_,
              peak_name + " / " + input_name)) {
        fit_range_low_ = x_->getMin("fitrange");
        fit_range_high_ = x_->getMax("fitrange");
        BuildDisplayHistogram();
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(1, "step");
    FixComponent(1, "low_exp");
    FixComponent(1, "low_lin");
    FixComponent(1, "high_exp");

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial constrained double peak fit failed"
                << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(1, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(1, best_chi2, best_vals, best_errs, best_const);

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    } else {
      std::cout << "ERROR: Constrained double peak fit failed" << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    SortPeaksByMu(2);
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitDoublePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.peaks[1] = ExtractPeakResult(1);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  }

  return results;
}

FitResult RooFitUtils::FitTriplePeak(const TString input_name,
                                     const TString peak_name,
                                     const FitResult &constrained_peaks,
                                     Double_t mu3_init) {
  FitResult results;
  results.peaks.emplace_back();
  results.peaks.emplace_back();
  results.peaks.emplace_back();

  AdoptSavedRange(input_name, peak_name);

  num_peaks_ = 3;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  for (Int_t pi = 0; pi < 2; pi++) {
    const PeakFitResult &cp = constrained_peaks.peaks[pi];
    BuildPeak(pi, cp.mu, cp.sigma, peak_height, range_width);
  }
  BuildPeak(2, mu3_init, sigma_init, peak_height, range_width);
  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  BuildUnbinnedData();
  x_->setRange(fit_range_low_, fit_range_high_);

  for (Int_t pi = 0; pi < 2; pi++) {
    const PeakFitResult &cp = constrained_peaks.peaks[pi];
    RooFitPeakModel &p = peaks_[pi];
    Double_t cga = cp.gaus_amplitude;
    p.mu->setVal(cp.mu);
    p.mu->setConstant(kTRUE);
    p.sigma->setVal(cp.sigma);
    p.sigma->setConstant(kTRUE);
    p.gaus_yield->setVal(cga);
    p.gaus_yield->setRange(0, peak_height * range_width * 10.0);
    p.gaus_yield->setConstant(kFALSE);
    p.ratio_step->setVal(cp.step_amplitude / cga);
    p.ratio_step->setConstant(kTRUE);
    p.ratio_low_exp->setVal(cp.low_exp_tail_amplitude / cga);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(cp.low_exp_tail_ratio > 0 ? cp.low_exp_tail_ratio
                                                    : 1.0);
    p.tau_low_exp->setConstant(kTRUE);
    p.ratio_low_lin->setVal(cp.low_lin_tail_amplitude / cga);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(cp.low_lin_tail_slope);
    p.slope_low_lin->setConstant(kTRUE);
    p.ratio_high_exp->setVal(cp.high_exp_tail_amplitude / cga);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(cp.high_exp_tail_ratio > 0 ? cp.high_exp_tail_ratio
                                                      : 1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
  ConfigureComponentFlagsForPeak(2);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(
              working_hist_, &events_, display_bin_width_kev_, total_pdf_, x_,
              unbinned_data_, &peaks_, &bkg_, fit_range_low_, fit_range_high_,
              peak_name + " / " + input_name)) {
        fit_range_low_ = x_->getMin("fitrange");
        fit_range_high_ = x_->getMax("fitrange");
        BuildDisplayHistogram();
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(2, "step");
    FixComponent(2, "low_exp");
    FixComponent(2, "low_lin");
    FixComponent(2, "high_exp");

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial triple peak fit failed" << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(2, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(2, best_chi2, best_vals, best_errs, best_const);

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Triple peak fit converged successfully" << std::endl;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    } else {
      std::cout << "ERROR: Triple peak fit failed to converge" << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    SortPeaksByMu(3);
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitTriplePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.peaks[1] = ExtractPeakResult(1);
    results.peaks[2] = ExtractPeakResult(2);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  }

  return results;
}

TString RooFitUtils::ParamFullName(const TString &channel,
                                   const TString &param) {
  return channel + ":" + param;
}

TString RooFitUtils::SourceForTarget(const TString &target) {
  for (size_t i = 0; i < sim_links_.size(); i++) {
    TString full_target =
        ParamFullName(sim_links_[i].target_channel, sim_links_[i].target_param);
    if (full_target == target) {
      return ParamFullName(sim_links_[i].source_channel,
                           sim_links_[i].source_param);
    }
  }
  return "";
}

RooRealVar *
RooFitUtils::ResolveOrCreate(const TString &channel, const TString &param_name,
                             std::map<TString, RooRealVar *> &registry,
                             Double_t init_val, Double_t lo, Double_t hi) {
  TString full = ParamFullName(channel, param_name);
  TString source = SourceForTarget(full);
  if (source.Length() > 0) {
    std::map<TString, RooRealVar *>::iterator it = registry.find(source);
    if (it == registry.end()) {
      std::cerr << "ERROR: source param " << source << " for link target "
                << full << " not yet built; check channel order." << std::endl;
      return nullptr;
    }
    registry[full] = it->second;
    return it->second;
  }
  RooRealVar *v = new RooRealVar(full.Data(), full.Data(), init_val, lo, hi);
  RegisterOwned(v);
  registry[full] = v;
  return v;
}

void RooFitUtils::AddChannel(const TString &name,
                             const std::vector<Double_t> &events,
                             Float_t fit_range_low, Float_t fit_range_high,
                             Float_t display_bin_width_kev, Int_t num_peaks,
                             const std::vector<Double_t> &mu_inits,
                             Bool_t use_flat_background, Bool_t use_step,
                             Bool_t use_low_exp_tail, Bool_t use_low_lin_tail,
                             Bool_t use_high_exp_tail,
                             const std::vector<Bool_t> &mu_fixed,
                             Bool_t bkg_yield_fixed, Bool_t bkg_slope_fixed,
                             Bool_t lock_shape_after_seed,
                             const std::vector<Bool_t> &use_step_per_peak) {
  if (!sim_mode_) {
    std::cerr << "ERROR: AddChannel called on a single-channel RooFitUtils "
                 "instance; construct with the default ctor for sim mode."
              << std::endl;
    return;
  }
  if ((Int_t)mu_inits.size() != num_peaks) {
    std::cerr << "ERROR: mu_inits size (" << mu_inits.size()
              << ") must match num_peaks (" << num_peaks << ")" << std::endl;
    return;
  }
  if (!mu_fixed.empty() && (Int_t)mu_fixed.size() != num_peaks) {
    std::cerr << "ERROR: mu_fixed size (" << mu_fixed.size()
              << ") must match num_peaks (" << num_peaks << ")" << std::endl;
    return;
  }
  if (!use_step_per_peak.empty() &&
      (Int_t)use_step_per_peak.size() != num_peaks) {
    std::cerr << "ERROR: use_step_per_peak size (" << use_step_per_peak.size()
              << ") must match num_peaks (" << num_peaks << ")" << std::endl;
    return;
  }
  RooFitChannelConfig cfg;
  cfg.name = name;
  cfg.events = events;
  cfg.hist = BuildDisplayHistogramFrom(events, fit_range_low, fit_range_high,
                                       display_bin_width_kev);
  cfg.fit_range_low = fit_range_low;
  cfg.fit_range_high = fit_range_high;
  cfg.display_bin_width_kev = display_bin_width_kev;
  cfg.num_peaks = num_peaks;
  cfg.mu_inits = mu_inits;
  cfg.mu_fixed =
      mu_fixed.empty() ? std::vector<Bool_t>(num_peaks, kFALSE) : mu_fixed;
  cfg.use_step_per_peak = use_step_per_peak;
  cfg.bkg_yield_fixed = bkg_yield_fixed;
  cfg.bkg_slope_fixed = bkg_slope_fixed;
  cfg.lock_shape_after_seed = lock_shape_after_seed;
  cfg.use_flat_background = use_flat_background;
  cfg.use_step = use_step;
  cfg.use_low_exp_tail = use_low_exp_tail;
  cfg.use_low_lin_tail = use_low_lin_tail;
  cfg.use_high_exp_tail = use_high_exp_tail;
  sim_channels_.push_back(cfg);
}

void RooFitUtils::LinkParameter(const TString &target, const TString &source) {
  Ssiz_t t_sep = target.Index(":");
  Ssiz_t s_sep = source.Index(":");
  if (t_sep < 0 || s_sep < 0) {
    std::cerr << "ERROR: LinkParameter expects 'channel:Param' format"
              << std::endl;
    return;
  }
  RooFitParamLink lk;
  lk.target_channel = TString(target(0, t_sep));
  lk.target_param = TString(target(t_sep + 1, target.Length() - t_sep - 1));
  lk.source_channel = TString(source(0, s_sep));
  lk.source_param = TString(source(s_sep + 1, source.Length() - s_sep - 1));
  sim_links_.push_back(lk);
}

void RooFitUtils::LinkPeakShape(const TString &target_channel,
                                Int_t target_peak,
                                const TString &source_channel,
                                Int_t source_peak) {
  TString t_suffix = TString::Format("%d", target_peak + 1);
  TString s_suffix = TString::Format("%d", source_peak + 1);
  const char *shape_params[9] = {"Mu",
                                 "Sigma",
                                 "StepAmplitude",
                                 "LowExpTailAmplitude",
                                 "LowExpTailRatio",
                                 "LowLinTailAmplitude",
                                 "LowLinTailSlope",
                                 "HighExpTailAmplitude",
                                 "HighExpTailRatio"};
  for (Int_t i = 0; i < 9; i++) {
    LinkParameter(target_channel + ":" + shape_params[i] + t_suffix,
                  source_channel + ":" + shape_params[i] + s_suffix);
  }
}

void RooFitUtils::SeedChannel(const TString &channel_name,
                              const FitResult &result) {
  sim_seeds_[channel_name] = result;
}

void RooFitUtils::BuildChannelModel(const RooFitChannelConfig &cfg,
                                    std::map<TString, RooRealVar *> &registry) {
  Double_t range_width = cfg.fit_range_high - cfg.fit_range_low;
  Double_t sigma_lo = cfg.hist->GetBinWidth(1);
  Double_t sigma_hi = range_width * 0.1;
  Double_t sigma_init = range_width * 0.05;
  if (sigma_init < sigma_lo)
    sigma_init = 2.0 * sigma_lo;
  if (sigma_init > sigma_hi)
    sigma_init = 0.5 * sigma_hi;
  Double_t peak_height = cfg.hist->GetBinContent(cfg.hist->GetMaximumBin());
  Double_t hist_xmin = cfg.hist->GetXaxis()->GetXmin();
  Double_t hist_xmax = cfg.hist->GetXaxis()->GetXmax();

  if (x_ == nullptr) {
    Double_t global_xmin = hist_xmin;
    Double_t global_xmax = hist_xmax;
    for (size_t i = 0; i < sim_channels_.size(); i++) {
      Double_t cmin = sim_channels_[i].hist->GetXaxis()->GetXmin();
      Double_t cmax = sim_channels_[i].hist->GetXaxis()->GetXmax();
      if (cmin < global_xmin)
        global_xmin = cmin;
      if (cmax > global_xmax)
        global_xmax = cmax;
    }
    x_ = new RooRealVar("x", "x", global_xmin, global_xmax);
    RegisterOwned(x_);
  }

  TString range_name = TString("fitrange_") + cfg.name;
  x_->setRange(range_name.Data(), cfg.fit_range_low, cfg.fit_range_high);
  sim_channel_range_names_[cfg.name] = range_name;

  std::vector<RooFitPeakModel> peaks;
  for (Int_t pi = 0; pi < cfg.num_peaks; pi++) {
    RooFitPeakModel p;
    TString suffix = TString::Format("%d", pi + 1);
    Double_t mu_init = cfg.mu_inits[pi];

    Int_t mu_bin = cfg.hist->FindBin(mu_init);
    Double_t local_height = cfg.hist->GetBinContent(mu_bin);
    Double_t bkg_left = 0;
    Double_t bkg_right = 0;
    Int_t lbin = cfg.hist->FindBin(cfg.fit_range_low);
    Int_t rbin = cfg.hist->FindBin(cfg.fit_range_high);
    Int_t nside = (rbin - lbin) / 10;
    if (nside < 1)
      nside = 1;
    for (Int_t i = 0; i < nside; i++) {
      bkg_left += cfg.hist->GetBinContent(lbin + i);
      bkg_right += cfg.hist->GetBinContent(rbin - i);
    }
    Double_t bkg_floor = (bkg_left + bkg_right) / (2.0 * nside);
    Double_t net_height = local_height - bkg_floor;
    if (net_height < 0.1 * local_height)
      net_height = 0.1 * local_height;
    Double_t total_init =
        net_height * sigma_init * TMath::Sqrt(2.0 * TMath::Pi());

    p.mu = ResolveOrCreate(cfg.name, "Mu" + suffix, registry, mu_init,
                           cfg.fit_range_low, cfg.fit_range_high);
    p.sigma = ResolveOrCreate(cfg.name, "Sigma" + suffix, registry, sigma_init,
                              sigma_lo, sigma_hi);
    p.gaus_yield =
        ResolveOrCreate(cfg.name, "GausAmplitude" + suffix, registry,
                        total_init, 0, peak_height * range_width * 10.0);
    p.ratio_step = ResolveOrCreate(cfg.name, "StepAmplitude" + suffix, registry,
                                   0.0, 0.0, 0.5);
    p.ratio_low_exp = ResolveOrCreate(cfg.name, "LowExpTailAmplitude" + suffix,
                                      registry, 0.0, 0.0, 0.5);
    p.tau_low_exp = ResolveOrCreate(cfg.name, "LowExpTailRatio" + suffix,
                                    registry, 1.5, 1.0, 100.0);
    p.ratio_low_lin = ResolveOrCreate(cfg.name, "LowLinTailAmplitude" + suffix,
                                      registry, 0.0, 0.0, 0.5);
    p.slope_low_lin = ResolveOrCreate(cfg.name, "LowLinTailSlope" + suffix,
                                      registry, 0.0, -0.1, 0.1);
    p.ratio_high_exp = ResolveOrCreate(
        cfg.name, "HighExpTailAmplitude" + suffix, registry, 0.0, 0.0, 0.5);
    p.tau_high_exp = ResolveOrCreate(cfg.name, "HighExpTailRatio" + suffix,
                                     registry, 1.5, 1.0, 100.0);

    TString pdf_suffix = "_" + cfg.name + "_" + suffix;
    p.gauss_pdf = RooFitFunctions::MakeGaussian("gauss_pdf" + pdf_suffix, *x_,
                                                *p.mu, *p.sigma);
    p.step_pdf = RooFitFunctions::MakeStepShelf("step_pdf" + pdf_suffix, *x_,
                                                *p.mu, *p.sigma);
    p.low_exp_pdf = RooFitFunctions::MakeLowExpTail(
        "low_exp_pdf" + pdf_suffix, *x_, *p.mu, *p.sigma, *p.tau_low_exp);
    p.low_lin_pdf = RooFitFunctions::MakeLowLinTail(
        "low_lin_pdf" + pdf_suffix, *x_, *p.mu, *p.sigma, *p.slope_low_lin);
    p.high_exp_pdf = RooFitFunctions::MakeHighExpTail(
        "high_exp_pdf" + pdf_suffix, *x_, *p.mu, *p.sigma, *p.tau_high_exp);
    RegisterOwned(p.gauss_pdf);
    RegisterOwned(p.step_pdf);
    RegisterOwned(p.low_exp_pdf);
    RegisterOwned(p.low_lin_pdf);
    RegisterOwned(p.high_exp_pdf);

    p.step_yield =
        new RooFormulaVar(("step_yield" + pdf_suffix).Data(), "@0*@1",
                          RooArgList(*p.gaus_yield, *p.ratio_step));
    p.low_exp_yield =
        new RooFormulaVar(("low_exp_yield" + pdf_suffix).Data(), "@0*@1",
                          RooArgList(*p.gaus_yield, *p.ratio_low_exp));
    p.low_lin_yield =
        new RooFormulaVar(("low_lin_yield" + pdf_suffix).Data(), "@0*@1",
                          RooArgList(*p.gaus_yield, *p.ratio_low_lin));
    p.high_exp_yield =
        new RooFormulaVar(("high_exp_yield" + pdf_suffix).Data(), "@0*@1",
                          RooArgList(*p.gaus_yield, *p.ratio_high_exp));
    RegisterOwned(p.step_yield);
    RegisterOwned(p.low_exp_yield);
    RegisterOwned(p.low_lin_yield);
    RegisterOwned(p.high_exp_yield);

    Bool_t peak_use_step = cfg.use_step;
    if (pi < (Int_t)cfg.use_step_per_peak.size())
      peak_use_step = cfg.use_step_per_peak[pi];
    if (!peak_use_step) {
      p.ratio_step->setVal(0.0);
      p.ratio_step->setConstant(kTRUE);
    }
    if (!cfg.use_low_exp_tail) {
      p.ratio_low_exp->setVal(0.0);
      p.ratio_low_exp->setConstant(kTRUE);
      p.tau_low_exp->setVal(1.0);
      p.tau_low_exp->setConstant(kTRUE);
    }
    if (!cfg.use_low_lin_tail) {
      p.ratio_low_lin->setVal(0.0);
      p.ratio_low_lin->setConstant(kTRUE);
      p.slope_low_lin->setVal(0.0);
      p.slope_low_lin->setConstant(kTRUE);
    }
    if (!cfg.use_high_exp_tail) {
      p.ratio_high_exp->setVal(0.0);
      p.ratio_high_exp->setConstant(kTRUE);
      p.tau_high_exp->setVal(1.0);
      p.tau_high_exp->setConstant(kTRUE);
    }

    peaks.push_back(p);
  }

  if (cfg.num_peaks > 1) {
    std::vector<Int_t> sorted_idx(cfg.num_peaks);
    for (Int_t i = 0; i < cfg.num_peaks; i++)
      sorted_idx[i] = i;
    std::sort(sorted_idx.begin(), sorted_idx.end(), [&](Int_t a, Int_t b) {
      return cfg.mu_inits[a] < cfg.mu_inits[b];
    });
    for (Int_t k = 0; k < cfg.num_peaks - 1; k++) {
      Int_t left = sorted_idx[k];
      Int_t right = sorted_idx[k + 1];
      Double_t midpoint = 0.5 * (cfg.mu_inits[left] + cfg.mu_inits[right]);
      peaks[left].mu->setMax(midpoint);
      peaks[right].mu->setMin(midpoint);
    }
  }

  RooFitBackgroundModel bkg;
  Double_t bkg_estimate = 0;
  Int_t lbin = cfg.hist->FindBin(cfg.fit_range_low);
  Int_t rbin = cfg.hist->FindBin(cfg.fit_range_high);
  Int_t nside = (rbin - lbin) / 10;
  if (nside < 1)
    nside = 1;
  Double_t bkg_left = 0;
  Double_t bkg_right = 0;
  for (Int_t i = 0; i < nside; i++) {
    bkg_left += cfg.hist->GetBinContent(lbin + i);
    bkg_right += cfg.hist->GetBinContent(rbin - i);
  }
  bkg_estimate = (bkg_left + bkg_right) / (2.0 * nside);

  bkg.bkg_yield = ResolveOrCreate(cfg.name, "BkgConstant", registry,
                                  bkg_estimate * range_width, 0,
                                  peak_height * range_width * 10.0);
  Double_t slope_lo = -0.9 / cfg.fit_range_high;
  Double_t slope_hi = 5.0 / range_width;
  bkg.bkg_slope =
      ResolveOrCreate(cfg.name, "BkgSlope", registry, 0.0, slope_lo, slope_hi);
  TString bkg_pdf_name = "bkg_pdf_" + cfg.name;
  bkg.bkg_pdf =
      RooFitFunctions::MakeLinearBackground(bkg_pdf_name, *x_, *bkg.bkg_slope);
  if (cfg.use_flat_background) {
    bkg.bkg_slope->setVal(0.0);
    bkg.bkg_slope->setConstant(kTRUE);
  }
  RegisterOwned(bkg.bkg_pdf);

  RooArgList pdf_list;
  RooArgList coef_list;
  for (size_t pi = 0; pi < peaks.size(); pi++) {
    pdf_list.add(*peaks[pi].gauss_pdf);
    coef_list.add(*peaks[pi].gaus_yield);
    pdf_list.add(*peaks[pi].step_pdf);
    coef_list.add(*peaks[pi].step_yield);
    pdf_list.add(*peaks[pi].low_exp_pdf);
    coef_list.add(*peaks[pi].low_exp_yield);
    pdf_list.add(*peaks[pi].low_lin_pdf);
    coef_list.add(*peaks[pi].low_lin_yield);
    pdf_list.add(*peaks[pi].high_exp_pdf);
    coef_list.add(*peaks[pi].high_exp_yield);
  }
  pdf_list.add(*bkg.bkg_pdf);
  coef_list.add(*bkg.bkg_yield);

  TString sum_name = "total_pdf_" + cfg.name;
  RooAddPdf *sum =
      new RooAddPdf(sum_name.Data(), sum_name.Data(), pdf_list, coef_list);
  RegisterOwned(sum);
  sim_channel_pdfs_[cfg.name] = sum;
  sim_channel_peaks_[cfg.name] = peaks;
  sim_channel_bkg_[cfg.name] = bkg;

  RooArgSet vars(*x_);
  RooDataSet *ds = new RooDataSet(("data_" + cfg.name).Data(),
                                  ("data_" + cfg.name).Data(), vars);
  Double_t xmin = x_->getMin();
  Double_t xmax = x_->getMax();
  for (size_t i = 0; i < cfg.events.size(); i++) {
    Double_t e = cfg.events[i];
    if (e < xmin || e > xmax)
      continue;
    x_->setVal(e);
    ds->add(vars);
  }
  sim_channel_data_[cfg.name] = ds;
}

void RooFitUtils::ApplySeedToChannel(const TString &channel) {
  std::map<TString, FitResult>::iterator it = sim_seeds_.find(channel);
  if (it == sim_seeds_.end())
    return;
  const FitResult &seed = it->second;
  std::vector<RooFitPeakModel> &peaks = sim_channel_peaks_[channel];
  RooFitBackgroundModel &bkg = sim_channel_bkg_[channel];
  for (size_t pi = 0; pi < peaks.size() && pi < seed.peaks.size(); pi++) {
    const PeakFitResult &cp = seed.peaks[pi];
    RooFitPeakModel &p = peaks[pi];
    Double_t cga = cp.gaus_amplitude;
    if (cp.mu >= 0)
      p.mu->setVal(cp.mu);
    if (cp.sigma >= 0)
      p.sigma->setVal(cp.sigma);
    if (cga >= 0)
      p.gaus_yield->setVal(cga);
    if (cp.step_amplitude >= 0 && cga > 0)
      p.ratio_step->setVal(cp.step_amplitude / cga);
    if (cp.low_exp_tail_amplitude >= 0 && cga > 0)
      p.ratio_low_exp->setVal(cp.low_exp_tail_amplitude / cga);
    if (cp.low_exp_tail_ratio > 0)
      p.tau_low_exp->setVal(cp.low_exp_tail_ratio);
    if (cp.low_lin_tail_amplitude >= 0 && cga > 0)
      p.ratio_low_lin->setVal(cp.low_lin_tail_amplitude / cga);
    if (cp.low_lin_tail_slope > -1)
      p.slope_low_lin->setVal(cp.low_lin_tail_slope);
    if (cp.high_exp_tail_amplitude >= 0 && cga > 0)
      p.ratio_high_exp->setVal(cp.high_exp_tail_amplitude / cga);
    if (cp.high_exp_tail_ratio > 0)
      p.tau_high_exp->setVal(cp.high_exp_tail_ratio);
  }
  if (seed.bkg_constant >= 0)
    bkg.bkg_yield->setVal(seed.bkg_constant);
  if (seed.lin_bkg_slope > -1e6)
    bkg.bkg_slope->setVal(seed.lin_bkg_slope);
}

void RooFitUtils::ApplyChannelMuLocks() {
  for (size_t ci = 0; ci < sim_channels_.size(); ci++) {
    const RooFitChannelConfig &cfg = sim_channels_[ci];
    std::vector<RooFitPeakModel> &peaks = sim_channel_peaks_[cfg.name];
    Int_t n_lock = TMath::Min((Int_t)peaks.size(), (Int_t)cfg.mu_fixed.size());
    for (Int_t pi = 0; pi < n_lock; pi++) {
      if (!cfg.mu_fixed[pi])
        continue;
      Double_t target = cfg.mu_inits[pi];
      peaks[pi].mu->setVal(target);
      peaks[pi].mu->setConstant(kTRUE);
    }
  }
}

void RooFitUtils::ApplyChannelBkgLocks() {
  for (size_t ci = 0; ci < sim_channels_.size(); ci++) {
    const RooFitChannelConfig &cfg = sim_channels_[ci];
    RooFitBackgroundModel &bkg = sim_channel_bkg_[cfg.name];
    if (cfg.bkg_yield_fixed && bkg.bkg_yield)
      bkg.bkg_yield->setConstant(kTRUE);
    if (cfg.bkg_slope_fixed && bkg.bkg_slope)
      bkg.bkg_slope->setConstant(kTRUE);
  }
}

void RooFitUtils::ApplyChannelShapeLocks() {
  for (size_t ci = 0; ci < sim_channels_.size(); ci++) {
    const RooFitChannelConfig &cfg = sim_channels_[ci];
    if (!cfg.lock_shape_after_seed)
      continue;
    std::vector<RooFitPeakModel> &peaks = sim_channel_peaks_[cfg.name];
    for (size_t pi = 0; pi < peaks.size(); pi++) {
      RooFitPeakModel &p = peaks[pi];
      if (p.sigma)
        p.sigma->setConstant(kTRUE);
      if (p.gaus_yield)
        p.gaus_yield->setConstant(kTRUE);
      if (p.ratio_step)
        p.ratio_step->setConstant(kTRUE);
      if (p.ratio_low_exp)
        p.ratio_low_exp->setConstant(kTRUE);
      if (p.tau_low_exp)
        p.tau_low_exp->setConstant(kTRUE);
      if (p.ratio_low_lin)
        p.ratio_low_lin->setConstant(kTRUE);
      if (p.slope_low_lin)
        p.slope_low_lin->setConstant(kTRUE);
      if (p.ratio_high_exp)
        p.ratio_high_exp->setConstant(kTRUE);
      if (p.tau_high_exp)
        p.tau_high_exp->setConstant(kTRUE);
    }
  }
}

Double_t RooFitUtils::ComputeChannelChi2(
    const TString &channel, const std::vector<RooFitPeakModel> & /*peaks*/,
    const RooFitBackgroundModel & /*bkg*/, Int_t &ndof) {
  RooAbsPdf *pdf = sim_channel_pdfs_[channel];
  RooFitChannelConfig const *cfg = nullptr;
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    if (sim_channels_[i].name == channel) {
      cfg = &sim_channels_[i];
      break;
    }
  }
  if (!cfg) {
    ndof = 0;
    return -1;
  }
  Double_t saved_val = x_->getVal();

  RooArgSet nset(*x_);
  Double_t total_exp = pdf->expectedEvents(&nset);
  Double_t bin_width = cfg->hist->GetBinWidth(1);

  Double_t chi2 = 0;
  Int_t nbins_in_range = 0;
  Int_t nbins_hist = cfg->hist->GetNbinsX();
  for (Int_t i = 1; i <= nbins_hist; i++) {
    Double_t xv = cfg->hist->GetBinCenter(i);
    if (xv < cfg->fit_range_low || xv > cfg->fit_range_high)
      continue;
    Double_t data = cfg->hist->GetBinContent(i);
    Double_t error = cfg->hist->GetBinError(i);
    if (error <= 0 || data <= 0)
      continue;
    x_->setVal(xv);
    Double_t fit_val = total_exp * pdf->getVal(&nset) * bin_width;
    Double_t residual = (data - fit_val) / error;
    chi2 += residual * residual;
    nbins_in_range++;
  }
  x_->setVal(saved_val);

  Int_t npars = 0;
  RooArgSet *params = pdf->getParameters(RooArgSet(*x_));
  for (Int_t i = 0; i < (Int_t)params->size(); i++) {
    RooRealVar *v = dynamic_cast<RooRealVar *>((*params)[i]);
    if (v && !v->isConstant())
      npars++;
  }
  delete params;
  ndof = nbins_in_range - npars;
  if (ndof <= 0)
    return -1;
  return chi2 / ndof;
}

PeakFitResult RooFitUtils::ExtractPeakResultFor(const RooFitPeakModel &p) {
  PeakFitResult r;
  r.mu = p.mu->getVal();
  r.mu_error = p.mu->getError();
  r.sigma = p.sigma->getVal();
  r.sigma_error = p.sigma->getError();
  r.gaus_amplitude = p.gaus_yield->getVal();
  r.gaus_amplitude_error = p.gaus_yield->getError();
  Double_t ga = r.gaus_amplitude;
  r.step_amplitude = p.ratio_step->getVal() * ga;
  r.step_amplitude_error = p.ratio_step->getError() * ga;
  r.low_exp_tail_amplitude = p.ratio_low_exp->getVal() * ga;
  r.low_exp_tail_amplitude_error = p.ratio_low_exp->getError() * ga;
  r.low_exp_tail_ratio = p.tau_low_exp->getVal();
  r.low_exp_tail_ratio_error = p.tau_low_exp->getError();
  r.low_lin_tail_amplitude = p.ratio_low_lin->getVal() * ga;
  r.low_lin_tail_amplitude_error = p.ratio_low_lin->getError() * ga;
  r.low_lin_tail_slope = p.slope_low_lin->getVal();
  r.low_lin_tail_slope_error = p.slope_low_lin->getError();
  r.high_exp_tail_amplitude = p.ratio_high_exp->getVal() * ga;
  r.high_exp_tail_amplitude_error = p.ratio_high_exp->getError() * ga;
  r.high_exp_tail_ratio = p.tau_high_exp->getVal();
  r.high_exp_tail_ratio_error = p.tau_high_exp->getError();
  return r;
}

void RooFitUtils::PlotChannel(const TString &channel, Int_t num_peaks,
                              const std::vector<RooFitPeakModel> &peaks,
                              const RooFitBackgroundModel &bkg,
                              const TString &input_name,
                              const TString &base_label,
                              const TString &chi2_label) {
  RooFitChannelConfig const *cfg = nullptr;
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    if (sim_channels_[i].name == channel) {
      cfg = &sim_channels_[i];
      break;
    }
  }
  if (!cfg)
    return;

  TH1 *saved_hist = working_hist_;
  Float_t saved_lo = fit_range_low_;
  Float_t saved_hi = fit_range_high_;
  std::vector<RooFitPeakModel> saved_peaks = peaks_;
  RooFitBackgroundModel saved_bkg = bkg_;
  Int_t saved_np = num_peaks_;
  RooAddPdf *saved_total = total_pdf_;

  working_hist_ = cfg->hist;
  fit_range_low_ = cfg->fit_range_low;
  fit_range_high_ = cfg->fit_range_high;
  peaks_ = peaks;
  bkg_ = bkg;
  num_peaks_ = num_peaks;
  total_pdf_ = static_cast<RooAddPdf *>(sim_channel_pdfs_[channel]);

  TString plot_name = base_label + "_" + channel;
  if (num_peaks == 1)
    PlotFitSinglePeak(input_name, plot_name, chi2_label);
  else if (num_peaks == 2)
    PlotFitDoublePeak(input_name, plot_name, chi2_label);
  else if (num_peaks == 3)
    PlotFitTriplePeak(input_name, plot_name, chi2_label);

  working_hist_ = saved_hist;
  fit_range_low_ = saved_lo;
  fit_range_high_ = saved_hi;
  peaks_ = saved_peaks;
  bkg_ = saved_bkg;
  num_peaks_ = saved_np;
  total_pdf_ = saved_total;
}

std::vector<FitResult> RooFitUtils::FitSimultaneous(const TString &input_name,
                                                    const TString &base_label) {
  std::vector<FitResult> results;
  if (sim_channels_.empty()) {
    std::cerr << "ERROR: FitSimultaneous called with no channels" << std::endl;
    return results;
  }

  AdoptSavedSimRange(input_name, base_label);

  std::map<TString, RooRealVar *> registry;
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    BuildChannelModel(sim_channels_[i], registry);
  }

  Float_t union_lo = sim_channels_[0].fit_range_low;
  Float_t union_hi = sim_channels_[0].fit_range_high;
  for (size_t i = 1; i < sim_channels_.size(); i++) {
    if (sim_channels_[i].fit_range_low < union_lo)
      union_lo = sim_channels_[i].fit_range_low;
    if (sim_channels_[i].fit_range_high > union_hi)
      union_hi = sim_channels_[i].fit_range_high;
  }
  x_->setRange(union_lo, union_hi);
  x_->setRange("fitrange", union_lo, union_hi);

  for (size_t i = 0; i < sim_channels_.size(); i++) {
    ApplySeedToChannel(sim_channels_[i].name);
  }
  ApplyChannelMuLocks();
  ApplyChannelBkgLocks();
  ApplyChannelShapeLocks();

  sim_category_ = new RooCategory("channel", "channel");
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    sim_category_->defineType(sim_channels_[i].name.Data());
  }

  std::map<std::string, RooDataSet *> data_map;
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    data_map[sim_channels_[i].name.Data()] =
        sim_channel_data_[sim_channels_[i].name];
  }
  sim_combined_data_ =
      new RooDataSet("combined_data", "combined_data", RooArgSet(*x_),
                     RooFit::Index(*sim_category_), RooFit::Import(data_map));

  sim_pdf_ = new RooSimultaneous("sim_pdf", "sim_pdf", *sim_category_);
  for (size_t i = 0; i < sim_channels_.size(); i++) {
    sim_pdf_->addPdf(*sim_channel_pdfs_[sim_channels_[i].name],
                     sim_channels_[i].name.Data());
  }

  std::cout << "Running simultaneous fit over " << sim_channels_.size()
            << " channels" << std::endl;

  RooFitResult *fit_result = nullptr;
  Bool_t sim_valid = kFALSE;

  if (interactive_) {
    if (LoadSimInteractiveParams(input_name, base_label)) {
      sim_valid = kTRUE;
      Float_t loaded_lo = x_->getMin("fitrange");
      Float_t loaded_hi = x_->getMax("fitrange");
      for (size_t i = 0; i < sim_channels_.size(); i++) {
        sim_channels_[i].fit_range_low = loaded_lo;
        sim_channels_[i].fit_range_high = loaded_hi;
        delete sim_channels_[i].hist;
        sim_channels_[i].hist = BuildDisplayHistogramFrom(
            sim_channels_[i].events, loaded_lo, loaded_hi,
            sim_channels_[i].display_bin_width_kev);
      }
      ApplyChannelMuLocks();
      ApplyChannelBkgLocks();
      ApplyChannelShapeLocks();
    } else {
      std::vector<SimEditorChannelView> views;
      for (size_t i = 0; i < sim_channels_.size(); i++) {
        const RooFitChannelConfig &cfg = sim_channels_[i];
        SimEditorChannelView v;
        v.name = cfg.name;
        v.hist = cfg.hist;
        v.events = &sim_channels_[i].events;
        v.display_bin_width_kev = cfg.display_bin_width_kev;
        v.pdf = sim_channel_pdfs_[cfg.name];
        v.data = sim_channel_data_[cfg.name];
        v.peaks = &sim_channel_peaks_[cfg.name];
        v.bkg = &sim_channel_bkg_[cfg.name];
        v.num_peaks = cfg.num_peaks;
        views.push_back(v);
      }
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      TString info = base_label + " / " + input_name;
      Bool_t accepted = LaunchInteractiveSimultaneousFitEditor(
          sim_pdf_, sim_combined_data_, x_, views, union_lo, union_hi, info,
          fit_debug_);
      gROOT->SetBatch(was_batch);
      sim_valid = accepted;
      if (accepted) {
        Float_t edited_lo = x_->getMin("fitrange");
        Float_t edited_hi = x_->getMax("fitrange");
        for (size_t i = 0; i < sim_channels_.size(); i++) {
          sim_channels_[i].fit_range_low = edited_lo;
          sim_channels_[i].fit_range_high = edited_hi;
          delete sim_channels_[i].hist;
          sim_channels_[i].hist = BuildDisplayHistogramFrom(
              sim_channels_[i].events, edited_lo, edited_hi,
              sim_channels_[i].display_bin_width_kev);
        }
        ApplyChannelMuLocks();
        ApplyChannelBkgLocks();
        ApplyChannelShapeLocks();
        SaveSimInteractiveParams(input_name, base_label);
      } else {
        std::cout << "Interactive sim fit cancelled" << std::endl;
      }
    }
  } else {
    if (fit_debug_) {
      // Un-suppress RooFit eval errors (the fit normally forces -1) so the
      // offending pdf is named, and report the seed NLL up front.
      RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);
      RooAbsReal *nll = sim_pdf_->createNLL(
          *sim_combined_data_, RooFit::Extended(kTRUE),
          RooFit::Range("fitrange"), RooFit::SplitRange(kTRUE));
      Double_t nll_seed = (nll != 0) ? nll->getVal() : 0.0;
      std::cout << "=== AU_ROOFIT_FIT_DEBUG: seed NLL = " << nll_seed
                << (std::isfinite(nll_seed) ? "" : "   <== NON-FINITE")
                << " ===" << std::endl;
      if (nll != 0)
        delete nll;
    }

    // SplitRange normalizes each channel over its own "fitrange_<name>" window
    // (set in BuildChannelModel) instead of the union range. Without it, a
    // channel whose peak sits below union_hi has its linear background and
    // exponential tails evaluated far past their own fit range: the background
    // polynomial 1+slope*x can go negative and the tails overflow, poisoning
    // the NLL once high statistics populate the out-of-range bins.
    Int_t print_level = fit_debug_ ? 1 : 0;
    Int_t eval_errors = fit_debug_ ? 10 : -1;
    fit_result = sim_pdf_->fitTo(
        *sim_combined_data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE),
        RooFit::Range("fitrange"), RooFit::SplitRange(kTRUE),
        RooFit::SumW2Error(kFALSE), RooFit::PrintLevel(print_level),
        RooFit::PrintEvalErrors(eval_errors), RooFit::Strategy(1),
        RooFit::Minimizer("Minuit2", "migrad"), BestAvailableBackend());
    sim_valid = (fit_result && fit_result->status() == 0);
    if (!sim_valid) {
      std::cout << "WARNING: simultaneous fit did not converge cleanly"
                << std::endl;
    }
  }

  for (size_t i = 0; i < sim_channels_.size(); i++) {
    const RooFitChannelConfig &cfg = sim_channels_[i];
    Int_t ndof = 0;
    Double_t chi2 = ComputeChannelChi2(cfg.name, sim_channel_peaks_[cfg.name],
                                       sim_channel_bkg_[cfg.name], ndof);
    std::cout << "Channel '" << cfg.name << "' chi2/ndf = " << chi2
              << " (ndof=" << ndof << ")" << std::endl;

    FitResult cr;
    for (Int_t pi = 0; pi < cfg.num_peaks; pi++) {
      cr.peaks.push_back(
          ExtractPeakResultFor(sim_channel_peaks_[cfg.name][pi]));
    }
    cr.bkg_constant = sim_channel_bkg_[cfg.name].bkg_yield->getVal();
    cr.bkg_constant_error = sim_channel_bkg_[cfg.name].bkg_yield->getError();
    cr.lin_bkg_slope = sim_channel_bkg_[cfg.name].bkg_slope->getVal();
    cr.lin_bkg_slope_error = sim_channel_bkg_[cfg.name].bkg_slope->getError();
    cr.reduced_chi2 = chi2;
    cr.valid = sim_valid;
    results.push_back(cr);

    TString chi2_label = Form("#chi^{2}/ndf = %.3f", chi2);
    PlotChannel(cfg.name, cfg.num_peaks, sim_channel_peaks_[cfg.name],
                sim_channel_bkg_[cfg.name], input_name, base_label, chi2_label);
  }

  delete fit_result;
  return results;
}
